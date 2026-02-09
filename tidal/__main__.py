"""Allow ``python -m tidal`` invocation."""

import sys

from tidal.cli import main

sys.exit(main())
