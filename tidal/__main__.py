"""Allow ``python -m tidal`` invocation."""

from __future__ import annotations

import sys

from tidal.cli import main

sys.exit(main())
