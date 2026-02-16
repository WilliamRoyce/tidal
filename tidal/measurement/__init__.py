"""Measurement and analysis tools for coupled-field simulations.

Provides energy computation, conversion probability, mixing length extraction,
spectral decomposition, and energy conservation diagnostics for simulation
output from the TIDAL Lagrangian-to-PDE pipeline.

Typical usage::

    from tidal.measurement import SimulationData, compute_conversion_probability
    from tidal.measurement import compute_mixing_length, compute_mixing_spectrum

    data = SimulationData.from_storage(storage, spec, grid, params)
    result = compute_conversion_probability(data, "phi_0", "chi_0")
    mixing = compute_mixing_length(result)
    print(f"L_mix = {mixing.mixing_length:.4f} +/- {mixing.mixing_length_uncertainty:.4f}")
"""

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
    compute_field_energy as compute_field_energy,
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

__all__ = [
    "ConversionResult",
    "DispersionResult",
    "EnergyDiagnostics",
    "FieldEnergy",
    "MixingResult",
    "MixingSpectrum",
    "SimulationData",
    "SpectralConversion",
    "SpectralPeak",
    "SpectralSnapshot",
    "SystemEnergy",
    "check_energy_conservation",
    "compute_conversion_probability",
    "compute_dispersion",
    "compute_energy_timeseries",
    "compute_field_energy",
    "compute_group_conversion",
    "compute_group_spectral_conversion",
    "compute_mixing_length",
    "compute_mixing_spectrum",
    "compute_mode_amplitudes",
    "compute_spectral_conversion",
    "compute_spectral_energy",
    "compute_spectrum",
    "compute_system_energy",
    "summarize",
]
