"""Allow ``python -m torsion_gertsenshtein.cli`` invocation."""

from __future__ import annotations

import sys

from torsion_gertsenshtein.cli import main

if __name__ == "__main__":
    sys.exit(main())
