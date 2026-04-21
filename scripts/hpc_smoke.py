"""Minimal smoke test — verify the tidal Python stack imports on a CSD3 compute node."""

import matplotlib as mpl
import numpy as np
import scipy
import sksundae

import tidal

print(f"numpy={np.__version__}")
print(f"scipy={scipy.__version__}")
print(f"matplotlib={mpl.__version__}")
print(f"sksundae={sksundae.__version__}")
print(f"tidal imported from {tidal.__file__}")
print("SMOKE OK")
