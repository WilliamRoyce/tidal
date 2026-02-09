"""Allow ``python -m tidal.cli`` invocation."""

from __future__ import annotations

import sys

from tidal.cli import main

if __name__ == "__main__":
    sys.exit(main())
