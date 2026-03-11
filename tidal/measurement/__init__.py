"""Measurement and analysis tools for coupled-field simulations.

Provides energy computation, conversion probability, mixing length extraction,
spectral decomposition, and energy conservation diagnostics for simulation
output from the TIDAL Lagrangian-to-PDE pipeline.

Typical usage::

    from tidal.measurement import SimulationData, compute_conversion_probability
    from tidal.measurement import compute_mixing_length, compute_mixing_spectrum

    data = SimulationData.from_storage(storage, spec, grid, params)
    # Or from a snapshot directory (memory-mapped, O(1) RAM):
    # data = SimulationData.load("output_dir", spec)
    result = compute_conversion_probability(data, "phi_0", "chi_0")
    mixing = compute_mixing_length(result)
    print(f"L_mix = {mixing.mixing_length:.4f} +/- {mixing.mixing_length_uncertainty:.4f}")
"""

from __future__ import annotations

from tidal.measurement._asymptotic import (
    AsymptoticConversionResult as AsymptoticConversionResult,
)
from tidal.measurement._asymptotic import (
    compute_asymptotic_conversion as compute_asymptotic_conversion,
)
from tidal.measurement._conversion import (
    ConversionResult as ConversionResult,
)
from tidal.measurement._conversion import (
    compute_conversion_probability as compute_conversion_probability,
)
from tidal.measurement._conversion import (
    compute_group_conversion as compute_group_conversion,
)
from tidal.measurement._diagnostics import (
    EnergyDiagnostics as EnergyDiagnostics,
)
from tidal.measurement._diagnostics import (
    check_energy_conservation as check_energy_conservation,
)
from tidal.measurement._diagnostics import (
    summarize as summarize,
)
from tidal.measurement._dispersion import (
    DispersionResult as DispersionResult,
)
from tidal.measurement._dispersion import (
    compute_dispersion as compute_dispersion,
)
from tidal.measurement._effective_mass import (
    EffectiveMassResult as EffectiveMassResult,
)
from tidal.measurement._effective_mass import (
    compute_effective_mass as compute_effective_mass,
)
from tidal.measurement._energy import (
    FieldEnergy as FieldEnergy,
)
from tidal.measurement._energy import (
    SystemEnergy as SystemEnergy,
)
from tidal.measurement._energy import (
    compute_energy_timeseries as compute_energy_timeseries,
)
from tidal.measurement._energy import (
    compute_system_energy as compute_system_energy,
)
from tidal.measurement._io import SimulationData as SimulationData
from tidal.measurement._mixing import (
    MixingResult as MixingResult,
)
from tidal.measurement._mixing import (
    MixingSpectrum as MixingSpectrum,
)
from tidal.measurement._mixing import (
    SpectralPeak as SpectralPeak,
)
from tidal.measurement._mixing import (
    compute_mixing_length as compute_mixing_length,
)
from tidal.measurement._mixing import (
    compute_mixing_spectrum as compute_mixing_spectrum,
)
from tidal.measurement._resonance import (
    ResonanceResult as ResonanceResult,
)
from tidal.measurement._resonance import (
    compute_resonance_analysis as compute_resonance_analysis,
)
from tidal.measurement._spectral import (
    SpectralSnapshot as SpectralSnapshot,
)
from tidal.measurement._spectral import (
    compute_mode_amplitudes as compute_mode_amplitudes,
)
from tidal.measurement._spectral import (
    compute_spectral_energy as compute_spectral_energy,
)
from tidal.measurement._spectral import (
    compute_spectrum as compute_spectrum,
)
from tidal.measurement._spectral_conversion import (
    SpectralConversion as SpectralConversion,
)
from tidal.measurement._spectral_conversion import (
    compute_group_spectral_conversion as compute_group_spectral_conversion,
)
from tidal.measurement._spectral_conversion import (
    compute_spectral_conversion as compute_spectral_conversion,
)
from tidal.measurement._velocity import (
    VelocityMismatchResult as VelocityMismatchResult,
)
from tidal.measurement._velocity import (
    VelocityResult as VelocityResult,
)
from tidal.measurement._velocity import (
    compute_velocities as compute_velocities,
)
from tidal.measurement._velocity import (
    compute_velocity_mismatch as compute_velocity_mismatch,
)
from tidal.measurement._writer import (
    SnapshotWriter as SnapshotWriter,
)
from tidal.measurement._writer import (
    compute_snapshot_count as compute_snapshot_count,
)
from tidal.measurement._writer import (
    create_snapshot_callback as create_snapshot_callback,
)

__all__ = [
    "AsymptoticConversionResult",
    "ConversionResult",
    "DispersionResult",
    "EffectiveMassResult",
    "EnergyDiagnostics",
    "FieldEnergy",
    "MixingResult",
    "MixingSpectrum",
    "ResonanceResult",
    "SimulationData",
    "SnapshotWriter",
    "SpectralConversion",
    "SpectralPeak",
    "SpectralSnapshot",
    "SystemEnergy",
    "VelocityMismatchResult",
    "VelocityResult",
    "check_energy_conservation",
    "compute_asymptotic_conversion",
    "compute_conversion_probability",
    "compute_dispersion",
    "compute_effective_mass",
    "compute_energy_timeseries",
    "compute_group_conversion",
    "compute_group_spectral_conversion",
    "compute_mixing_length",
    "compute_mixing_spectrum",
    "compute_mode_amplitudes",
    "compute_resonance_analysis",
    "compute_snapshot_count",
    "compute_spectral_conversion",
    "compute_spectral_energy",
    "compute_spectrum",
    "compute_system_energy",
    "compute_velocities",
    "compute_velocity_mismatch",
    "create_snapshot_callback",
    "summarize",
]
