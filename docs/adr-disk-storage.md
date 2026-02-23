# ADR: Disk-Backed Snapshot Storage for Long Simulations

**Status:** Accepted
**Date:** 2026-02-16
**Context:** TIDAL simulation pipeline stores entire time history in RAM, killing long-running simulations.

## Problem Statement

TIDAL stores **every snapshot in RAM** twice:

1. **During simulation**: py-pde's `MemoryStorage` accumulates all snapshots in a Python list.
2. **During measurement**: `SimulationData` eagerly loads the full NPZ into contiguous arrays.

Memory formula: `n_snapshots x (n_fields + n_momenta) x prod(grid_shape) x 8 bytes`

| Example                    | Grid | Snapshots | Slots | Memory     |
| -------------------------- | ---- | --------- | ----- | ---------- |
| coupled_proca (2+1D)       | 96^2 | 1,250     | 10    | **880 MB** |
| gravitational_waves (3+1D) | 64^3 | 500       | 14    | **14 GB**  |

At 3+1D grid sizes or long time series, the process is OOM-killed well before the simulation completes.

## Requirements

- **O(1) memory** in snapshot count for both writing and reading.
- **Framework-agnostic**: format must be readable from Python, Julia, C, Rust (future solver migration).
- **HPC-ready**: works on shared filesystems (NFS, Lustre, GPFS), crash-resilient (HPC wall-time kills).
- **Zero new dependencies**: only numpy.
- **Backward compatible**: existing NPZ workflows must still work.

## Alternatives Evaluated

### 1. HDF5 (PDEBench approach)

[PDEBench](https://github.com/pdebench/PDEBench) uses HDF5 for PDE simulation datasets with chunked access and per-dataset compression.

**Pros:** Mature, multi-language (C/Fortran/Julia/Python), chunked random access, widely used in HPC.
**Cons:** Requires `h5py` dependency (C library build), breaks zero-new-deps principle. File locking can be problematic on some shared filesystems.
**Decision:** Rejected. Adds non-trivial dependency; HDF5 file locking has known issues on NFS.

### 2. Zarr

[Zarr](https://zarr.dev/) offers chunked + compressed N-dimensional array storage with cloud support (S3, GCS).

**Pros:** Modern Python API, excellent multi-axis slicing, compression (Blosc, zstd), cloud-friendly, append-friendly.
**Cons:** Adds `zarr` + `numcodecs` dependencies, Python-centric (no native Julia/C reader), compression provides poor ratio on float64 physics data (typically 1.5-2x vs 8x for integer data).
**Decision:** Deferred as future extension. If cloud storage or compression becomes needed, a `ZarrWriter` can implement the same `append()`/`close()` API without changing measurement code. See [mmap vs Zarr/HDF5 analysis](https://pythonspeed.com/articles/mmap-vs-zarr-hdf5/).

### 3. py-pde FileStorage

py-pde provides built-in `FileStorage` for HDF5-based streaming during simulation.

**Pros:** Direct integration with solver, streaming.
**Cons:** Still needs `h5py`, couples storage format to py-pde (blocks Julia migration), py-pde-specific API.
**Decision:** Rejected. Framework coupling is a non-starter for planned Julia migration.

### 4. Incremental Measurement During Simulation

Compute energy/conversion at each snapshot, discard field data immediately. O(1) memory.

**Pros:** Absolute minimum memory.
**Cons:** User cannot do post-hoc analysis (new measurement types, debugging, visualization). Defeats the purpose of `tidal measure` CLI.
**Decision:** Rejected. Full field data must be available for flexible post-hoc analysis.

### 5. NPZ with Memory Mapping

numpy supports `mmap_mode` when loading NPZ files.

**Pros:** No format change needed.
**Cons:** Individual arrays within a ZIP-based NPZ are not contiguous on disk. Memory mapping through the ZIP layer is unreliable for large files. The per-snapshot keying (`phi_0_t0`, `phi_0_t1`, ...) prevents true contiguous time-axis access.
**Decision:** Rejected. NPZ's ZIP container prevents efficient mmap.

### 6. Directory of `.npy` Files with Memory Mapping (Chosen)

Pre-allocate one `.npy` file per field with shape `(n_snapshots, *grid_shape)`, write via `numpy.memmap(mode='w+')`, read via `numpy.load(mmap_mode='r')`.

**Pros:**

- Zero dependencies (pure numpy).
- C-contiguous time-first layout matches our sequential access pattern (snapshot-by-snapshot iteration).
- Julia reads `.npy` natively (`NPZ.jl`, `Mmap.mmap()`). C/Rust have `.npy` parsers.
- Pre-allocated files provide crash resilience: all flushed snapshots survive HPC wall-time kills.
- Exact snapshot count is known at start (`int(t_end / interval) + 1`), so no upper-bound estimation or truncation needed.
- Works on any POSIX filesystem, Windows, VMs, HPC clusters.

**Cons:**

- No compression (float64 data compresses poorly anyway; ~1.5-2x with Blosc).
- No random-access chunking (but our access is sequential along time axis, so this is fine).
- One file per field (manageable; typical simulations have 4-10 fields).

**Decision:** Accepted.

## Design

### Directory Layout

```
output_dir/
  metadata.json          # Grid, parameters, field list, snapshot count
  times.npy              # shape (n_snapshots,), float64
  phi_0.npy              # shape (n_snapshots, *grid_shape), float64
  chi_0.npy              # shape (n_snapshots, *grid_shape), float64
  pi_phi_0.npy           # shape (n_snapshots, *grid_shape), float64
  pi_chi_0.npy           # shape (n_snapshots, *grid_shape), float64
```

### metadata.json Schema

```json
{
  "version": 1,
  "n_snapshots": 1250,
  "grid_spacing": [0.5208, 0.5208],
  "grid_bounds": [
    [0, 50],
    [0, 50]
  ],
  "grid_shape": [96, 96],
  "periodic": [true, true],
  "parameters": { "mA": 1.0, "mB": 1.0, "g": 0.1 },
  "spec_path": "examples/data/coupled_proca.json",
  "fields": ["A_1", "A_2", "B_1", "B_2"],
  "momenta": ["A_1", "A_2", "B_1", "B_2"],
  "dtype": "float64"
}
```

### Exact Snapshot Count

Since `t_end` and `snapshot_interval` are both known at simulation start:

```python
n_snapshots = int(t_end / snapshot_interval) + 1   # +1 for initial state at t=0
```

This is exact. No upper-bound estimation, no truncation needed. The `.npy` files are pre-allocated with the correct shape.

### Writing: SnapshotWriter

`SnapshotWriter` uses `numpy.lib.format.open_memmap()` to pre-allocate proper `.npy` files (with headers), then writes each snapshot in-place:

```python
# Each append() is O(grid_size) — writes to the next row in the memmap
writer = SnapshotWriter(output_dir, field_names, momentum_names,
                         grid_shape, n_snapshots, ...)
for t, fields, momenta in simulation:
    writer.append(t, fields, momenta)  # O(1) memory
writer.close()  # writes metadata.json
```

Key properties:

- Each `append()` writes to pre-allocated memmap rows. O(1) memory, O(grid_size) I/O.
- `flush()` on each append for crash resilience.
- `metadata.json` written at `close()`. If missing (crash), snapshot count is recoverable from `times.npy` by finding the last non-zero entry.

### Reading: SimulationData.from_directory()

```python
data = SimulationData.from_directory(output_dir, spec)
# data.fields["phi_0"] is a np.memmap — O(1) RAM
# Only pages actually accessed are loaded by the OS page cache
```

`mmap_mode="r"` gives zero-copy, read-only access. The OS manages which pages are in RAM via the page cache. Sequential access (iterate snapshots) is optimal because the data is C-contiguous along the time axis.

### CLI Integration

- `tidal simulate spec.json --output /tmp/results/` creates a snapshot directory (disk-backed, O(1) memory).
- `tidal measure /tmp/results/ --spec spec.json` loads from the snapshot directory.
- `SimulationData.load(path, spec)` is the universal entry point (directory only; NPZ support was removed).

### Crash Recovery

When `metadata.json` is missing (writer wasn't closed due to crash/kill):

1. Look for `times.npy` in the directory.
2. Memory-map it read-only.
3. Find the last index where `times[i] > 0` (unwritten entries are 0.0 from memmap pre-allocation).
4. Use that count as `n_snapshots`. All prior snapshots are intact.

### Access Pattern Analysis

| Measurement Function             | Access Pattern                 | Mmap Behavior                                   |
| -------------------------------- | ------------------------------ | ----------------------------------------------- |
| `compute_energy_timeseries`      | Sequential over time           | Optimal: sequential page faults along time axis |
| `compute_conversion_probability` | Sequential over time           | Optimal: same as above                          |
| `compute_spectrum`               | All snapshots of one field     | Full scan: pages loaded sequentially            |
| `compute_dispersion`             | Full field + temporal FFT      | Full scan: all pages loaded                     |
| `compute_mixing_length`          | Derived from conversion result | No direct field access                          |

All access patterns are sequential along the first (time) axis, which is exactly what C-contiguous memmap optimizes for.

## HPC and Cross-Platform Compatibility

- **No dependencies beyond numpy**: Works on any system with Python + numpy.
- **POSIX + Windows**: `numpy.memmap` abstracts OS-level mmap calls.
- **Shared filesystems** (NFS, Lustre, GPFS): Sequential I/O pattern (write forward, read forward) is optimal. No random access, no file locking.
- **Crash resilience**: Files pre-allocated at creation. `append()` flushes per snapshot. If job is killed (HPC wall-time), all flushed snapshots are intact.
- **Julia interop**: Julia's `NPZ.jl` or `Mmap.mmap()` reads `.npy` files natively. A Julia solver can write the same directory layout.
- **Large files**: 64-bit systems handle files >2GB. A 3+1D 64^3 x 500 snapshot field = 1 GB per `.npy` file, well within limits.
- **Virtual machines**: Standard filesystem operations. No special kernel modules or shared memory needed.

## Future Extensions

1. **Zarr backend**: If compression or cloud storage is needed, a `ZarrWriter` can implement the same `append()`/`close()` API. Measurement code is unchanged because `SimulationData` provides the abstraction layer.
2. **Constraint field exclusion**: `SnapshotWriter` can skip constraint fields (`time_derivative_order == 0`) that are re-derivable from dynamical fields, saving ~20% storage for systems like coupled_proca.
3. **Parallel writing**: For domain-decomposed simulations, each rank writes its own subdirectory; a merge step concatenates along spatial axes.

## References

- [PDEBench data format (HDF5)](https://github.com/pdebench/PDEBench) — similar PDE simulation storage problem, chose HDF5
- [mmap vs Zarr/HDF5 performance](https://pythonspeed.com/articles/mmap-vs-zarr-hdf5/) — mmap is optimal when access aligns with storage layout
- [numpy.lib.format.open_memmap](https://numpy.org/doc/stable/reference/generated/numpy.lib.format.open_memmap.html) — creates proper `.npy` with header
- [py-pde CallbackTracker](https://py-pde.readthedocs.io/) — `CallbackTracker(func, interrupts=interval)` for streaming snapshots
