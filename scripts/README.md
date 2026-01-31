# Scripts

Utility scripts for project setup and maintenance.

## Wolfram Engine Setup

The project uses Wolfram Engine for symbolic computation (xCoba tensor algebra).

### Prerequisites

1. **Download the Wolfram Engine installer** from [wolfram.com/engine](https://www.wolfram.com/engine/)
2. Place the downloaded `.sh` file in the `third_party/` directory
3. The installer should be named like: `WolframEngine_14.3.0_LIN.sh`

### Installation

#### Automatic (Dev Container)

If you place the installer in `third_party/` before building the dev container, it will be installed automatically during container creation.

#### Manual Installation

```bash
# Run the installation script (requires sudo)
sudo ./scripts/install-wolfram-engine.sh

# Or if installer is already downloaded:
sudo ./scripts/install-wolfram-engine.sh --skip-download
```

### Activation

After installation, you must activate Wolfram Engine with a free Wolfram ID:

```bash
# Interactive activation (recommended)
./scripts/activate-wolfram.sh

# Or directly:
wolframscript -activate
```

You will be prompted for your Wolfram ID credentials. If you don't have one, create a free account at [account.wolfram.com](https://account.wolfram.com/).

### Verification

Check that Wolfram Engine is working:

```bash
# Check activation status
./scripts/activate-wolfram.sh --check

# Or test directly:
wolframscript -code "1+1"
# Should output: 2

# Test symbolic computation:
wolframscript -code "Integrate[x^2, {x, 0, 1}]"
# Should output: 1/3
```

### Troubleshooting

**"Wolfram Engine is not activated"**

- Run `wolframscript -activate` and enter your Wolfram ID credentials

**"wolframscript not found"**

- Ensure installation completed successfully
- Check that `/usr/local/bin` is in your PATH

**License limitations**

- Free Wolfram Engine license has 2GB memory limit
- For larger computations, consider a commercial license

### Files

- `install-wolfram-engine.sh` - Downloads and installs Wolfram Engine
- `activate-wolfram.sh` - Helps with license activation

### Environment Variables

| Variable              | Default                            | Description                       |
| --------------------- | ---------------------------------- | --------------------------------- |
| `WOLFRAM_VERSION`     | `14.3.0`                           | Wolfram Engine version to install |
| `WOLFRAM_INSTALL_DIR` | `/usr/local/Wolfram/WolframEngine` | Installation directory            |
