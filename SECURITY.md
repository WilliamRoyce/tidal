# Security Policy

## Supported Versions

The following versions of this project are currently supported with security updates:

| Version | Supported          |
| ------- | ------------------ |
| 0.45.x   | :white_check_mark: |
| < 0.45   | :x:                |

## Scope

This security policy covers:

- **Python Package** (`tidal/`)
  - PDE builders and simulation code
  - JSON parsing and validation
  - Initial condition generators

- **Wolfram Scripts** (`tidal/wolfram/`, `*.wls`)
  - xAct/Mathematica symbolic computation pipeline
  - JSON export routines

- **Development Container** (`.devcontainer/`)
  - Setup scripts and configuration
  - Wolfram Engine installation

## Reporting a Vulnerability

If you discover a security vulnerability, please report it via **GitHub Security Advisories**:

1. Go to the [Security tab](https://github.com/WilliamRoyce/torsion-gertsenshtein/security)
2. Click "Report a vulnerability"
3. Fill out the form with details about the vulnerability

### What to Include

Please include the following information in your report:

- Description of the vulnerability
- Steps to reproduce the issue
- Affected components (Python, Wolfram, dev container)
- Potential impact
- Any suggested fixes (if available)

### Response Timeline

- **Initial Response**: Within 48 hours of report submission
- **Status Update**: Within 7 days with assessment and planned fix timeline
- **Resolution**: Depends on severity and complexity

## Security Best Practices

When contributing to this project:

- Do not commit sensitive data (API keys, credentials, private data)
- Review `eval()` usage in `pde_builder.py` - it's sandboxed with `__builtins__={}` for safety
- Validate all user inputs, especially JSON schemas
- Use the provided dev container for reproducible, isolated environments
- Keep dependencies up to date via `uv` lock file

## Out of Scope

The following are explicitly out of scope for security reports:

- Theoretical physics errors (use regular issue tracker)
- Numerical precision/stability issues (use regular issue tracker)
- Performance optimizations (use regular issue tracker)
- Documentation typos or clarity (use regular issue tracker)

Thank you for helping keep this project secure!
