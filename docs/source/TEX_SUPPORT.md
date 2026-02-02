# TeX Support Configuration

By default, the development container does not install TeX to keep rebuild times fast. However, TeX support can be enabled for high-quality LaTeX-rendered plots.

## Enabling TeX Support

To enable TeX installation in the development container:

1. **Temporary**: Set the environment variable before rebuilding:
   ```bash
   export INSTALL_TEX=true
   ```
   Then rebuild the container.

2. **Permanent**: Edit `.devcontainer/devcontainer.json` and change:
   ```json
   "containerEnv": {
     "INSTALL_TEX": "true"
   }
   ```

## Using Plots Without TeX

The plotting functionality automatically falls back to standard matplotlib when TeX is not available:

```python
from torsion_gertsenshtein import has_tex_support
from torsion_gertsenshtein.plot_pgf import enable_pgf

# This will work with or without TeX
enable_pgf("xelatex")  # Automatically falls back if TeX unavailable

# Check if TeX is available
if has_tex_support():
    print("High-quality LaTeX rendering available")
else:
    print("Using standard matplotlib fonts")
```

## What TeX Provides

- **With TeX**: Vector-based LaTeX-rendered text, mathematical symbols, and fonts
- **Without TeX**: Standard matplotlib fonts with similar styling but bitmap text rendering

## Package Size Impact

Installing TeX adds approximately 2-3 GB to the container and increases rebuild time by several minutes.