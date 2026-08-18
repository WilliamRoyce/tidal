---
name: spec
description: Read and check equation specification JSONs — coefficients, sign consistency, families, and diffs between derivations. Use when inspecting a spec, comparing a re-derivation, or checking whether a component's signs are consistent, instead of reading the JSON directly.
---

# Reading equation specifications

A spec JSON is ~96,000 tokens; these commands answer a question in a few hundred, and
the answer comes from a vetted accessor rather than inference by eye. Six confidently
wrong diagnoses came from ad-hoc `json.load` scans (GH #401) — that is what this
replaces.

## Start here for corpus-level questions

`tests/data/spec_semantics.txt` already records families, index structure and proven
sign conflicts for every spec in `examples/data/`. **Read it rather than scanning.**
Regenerate with `python -m scripts.spec_semantics_report`; a test keeps it current.

## Recipes

**"Is this component's sign wrong?"** — do not compare raw coefficients across
components; that is misreading #1. Use the check that normalizes by the kinetic
coefficient:

```bash
tidal validate SPEC          # reports proven sibling sign conflicts; exit 1 if any
```

**"What is this coefficient?"**

```bash
tidal inspect SPEC --coefficient 'h_5:identity(h_5)'
```

Reports the summed coefficient, the kinetic divisor, the proven sign with its deciding
tactic, and every other place the quantity is recorded — separating parts of it from
redundant re-encodings from the *related but distinct* Hamiltonian term.

**"Did re-derivation change the physics?"**

```bash
tidal inspect OLD.json --diff NEW.json      # exit 1 = real change, 0 = none
```

Distinguishes a real change from an equation merely rescaled on both sides. A naive diff
reports the latter as three separate "fixes".

**"What does this equation say?"**

```bash
tidal inspect SPEC --equation h_5           # also accepts 'a_0,a_1' or 'all'
tidal inspect SPEC --detail summary         # proven signs only, whole spec, ~4k tokens
```

**"Which components belong together?"**

```bash
tidal inspect SPEC --families
```

Grouped by `tensor_head` and classified by `tensor_indices` — never by the numeric
suffix. Index 0 is the temporal component of a rank-1 field, but rank-3 torsion has
components like `t_13 = [2, 0, 2]` where that reading is meaningless.

## What to trust

- A sign verdict says **`unknown` rather than guessing**. Only ~8% of coefficients have
  a provable sign; the rest are free sweep parameters whose sign genuinely is not
  knowable. Supply `--assume-positive` / `--assume-nonzero` when you have physical
  grounds — `kappa` cannot vanish, but `xi` and `b5` genuinely reach zero.
- Every verdict names the **tactic** that decided it, so it can be audited without
  re-deriving.
- Prefer the **text output over `--json`**: measured, `--json` costs ~3x more and buys
  nothing unless a script is parsing it.

## Reading the equation listing

```text
[-kappa^(-2)] d2_t(h_5) =
    + [(B0) + (-2*B0^3*rho)] gradient_x(a_1)
    + [-kappa^(-2)] laplacian_x(h_5)   eff>0
```

The kinetic coefficient sits on the LHS, so the line *is* the equation. Each
`(operator, field)` key appears once with its terms summed. `eps:0,1` marks a key
spanning perturbative orders.

`eff>0` is the proven sign of the coefficient **after** dividing by the LHS kinetic
coefficient — **not** the sign of the bracket beside it. Above, `laplacian_x(h_5)`
prints a negative `[-kappa^(-2)]` yet is `eff>0`, because the LHS carries the same
factor and the two cancel: a standard wave operator. Reading the bracket's sign as the
physical sign is the mistake; that is why the marker is prefixed.

## Python API

`tidal.symbolic.spec_query` — `effective_coefficient`, `field_families`,
`coefficient_provenance`, `diff_systems`, `sibling_sign_conflicts`.
`tidal.symbolic.sign_algebra` — the sound sign/ratio decisions underneath.
