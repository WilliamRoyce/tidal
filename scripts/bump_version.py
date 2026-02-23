#!/usr/bin/env python3
"""Version management script for TIDAL.

Updates version numbers across:
- pyproject.toml
- tidal/__init__.py
- CITATION.cff (version + date)
- uv.lock (via uv lock regeneration)

Note: docs/source/conf.py reads the version dynamically via
importlib.metadata.version("tidal"), so it does not need updating.

Usage:
    python scripts/bump_version.py              # Interactive mode
    python scripts/bump_version.py 0.3.0        # Direct mode
    python scripts/bump_version.py 0.3.0 --commit  # With git commit
    python scripts/bump_version.py 0.3.0 --dry-run # Preview changes
"""

import argparse
import re
import shutil
import subprocess  # noqa: S404  # Needed for git/uv commands
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

# Optional dependencies for safer file handling
try:
    import tomli  # type: ignore[import-not-found]
    import tomli_w  # type: ignore[import-not-found]

    has_tomli_w = True
except ImportError:
    has_tomli_w = False

try:
    from ruamel.yaml import YAML  # type: ignore[import-not-found]

    has_ruamel_yaml = True
except ImportError:
    has_ruamel_yaml = False


class VersionBumpError(Exception):
    """Raised when version bump fails."""


class VersionBumper:
    """Handles atomic version updates across project files."""

    def __init__(
        self,
        new_version: str,
        *,
        dry_run: bool = False,
        allow_dirty: bool = False,
        commit: bool = False,
    ) -> None:
        """Initialize the version bumper.

        Parameters
        ----------
        new_version : str
            The new version to set (semantic versioning format)
        dry_run : bool
            If True, show what would change without modifying files
        allow_dirty : bool
            If True, allow uncommitted changes in git
        commit : bool
            If True, create a git commit after updating

        """
        self.new_version = new_version
        self.dry_run = dry_run
        self.allow_dirty = allow_dirty
        self.commit = commit
        self.root = Path(__file__).parent.parent
        self.backups: list[Path] = []

        # Files to update
        # Note: docs/source/conf.py reads version dynamically via
        # importlib.metadata, so it does not need explicit updating.
        self.pyproject_toml = self.root / "pyproject.toml"
        self.init_py = self.root / "tidal" / "__init__.py"
        self.citation_cff = self.root / "CITATION.cff"

    def run(self) -> None:
        """Execute the three-phase update process."""
        try:
            print(f"Version Bump: {self.get_current_version()} → {self.new_version}")
            print("=" * 60)
            print()

            self.validate()
            self.update_files()
            self.verify_and_commit()

            print()
            print("✓ Version bump complete!")
            if not self.dry_run and not self.commit:
                print()
                print("Next steps:")
                print("  1. Review changes: git diff")
                print(
                    f'  2. Commit changes: git add -A && git commit -m "chore: bump version to {self.new_version}"'
                )
                print(f"  3. Create tag: git tag v{self.new_version}")
                print("  4. Push: git push && git push --tags")

        except Exception as e:  # noqa: BLE001  # Intentional catch-all for rollback
            print()
            print(f"✗ Error: {e}", file=sys.stderr)
            if not self.dry_run:
                self.rollback()
            sys.exit(1)

    def validate(self) -> None:
        """Phase 1: Validate inputs and preconditions.

        Raises
        ------
        VersionBumpError
            If version format is invalid, required files are missing,
            git working directory has uncommitted changes, or other
            preconditions are not met.

        """
        print("Phase 1: Validation")

        # Validate version format
        if not self._is_valid_version(self.new_version):
            msg = (
                f"Invalid version format: {self.new_version}\n"
                "Expected semantic versioning (e.g., 0.3.0, 1.0.0-alpha.1)"
            )
            raise VersionBumpError(msg)
        print("  ✓ Version format valid (semantic versioning)")

        # Check all files exist
        for file_path in [self.pyproject_toml, self.init_py, self.citation_cff]:
            if not file_path.exists():
                msg = f"Required file not found: {file_path}\nEnsure you are running this script from the project root and that all files are present."
                raise VersionBumpError(msg)
        print("  ✓ All target files exist")

        # Check git status
        if (
            not self.allow_dirty
            and not self.dry_run
            and self._has_uncommitted_changes()
        ):
            msg = (
                "Git working directory has uncommitted changes.\n"
                "Commit or stash changes first, or use --allow-dirty flag."
            )
            raise VersionBumpError(msg)
        print(
            "  ✓ Git working directory clean"
            if not self.allow_dirty
            else "  ⚠ Skipping git check (--allow-dirty)"
        )

        # Report current version inconsistencies
        current_versions = self._get_all_current_versions()
        unique_versions = set(current_versions.values())
        if len(unique_versions) > 1:
            print("  ⚠ Version mismatch detected:")
            for file_name, version in current_versions.items():
                print(f"    - {file_name}: {version}")
        else:
            print(f"  ✓ Current version consistent: {next(iter(unique_versions))}")

        print()

    def update_files(self) -> None:
        """Phase 2: Update all files atomically."""
        print("Phase 2: Update Files")

        if self.dry_run:
            print("  [DRY RUN] Showing what would be updated:")
        else:
            print("  ✓ Creating backups")
            self._create_backups()

        # Update pyproject.toml
        if self.dry_run:
            print("  [DRY RUN] Would update pyproject.toml")
        else:
            self._update_pyproject_toml()
            print("  ✓ Updated pyproject.toml")

        # Update __init__.py
        if self.dry_run:
            print("  [DRY RUN] Would update __init__.py")
        else:
            self._update_init_py()
            print("  ✓ Updated __init__.py")

        # Update CITATION.cff
        today = datetime.now().strftime("%Y-%m-%d")  # noqa: DTZ005  # Local date for citation file
        if self.dry_run:
            print(f"  [DRY RUN] Would update CITATION.cff (version + date: {today})")
        else:
            self._update_citation_cff(today)
            print(f"  ✓ Updated CITATION.cff (version + date: {today})")

        # Regenerate uv.lock
        if self.dry_run:
            print("  [DRY RUN] Would regenerate uv.lock (via uv lock)")
        else:
            self._regenerate_uv_lock()
            print("  ✓ Regenerated uv.lock (via uv lock)")

        print()

    def verify_and_commit(self) -> None:
        """Phase 3: Verify updates and optionally commit.

        Raises
        ------
        VersionBumpError
            If post-update verification fails or git commit creation fails.

        """
        print("Phase 3: Verification")

        if self.dry_run:
            print("  [DRY RUN] Skipping verification")
            print()
            return

        # Verify all versions match
        current_versions = self._get_all_current_versions()
        for file_name, version in current_versions.items():
            if version != self.new_version:
                msg = f"Post-update verification failed: {file_name} has {version}, expected {self.new_version}"
                raise VersionBumpError(msg)
        print(f"  ✓ All files updated to {self.new_version}")

        # Show changes summary
        self._show_changes_summary(current_versions)

        # Clean up backups
        self._clean_backups()
        print("  ✓ Cleaned up backups")

        # Create git commit if requested
        if self.commit:
            self._create_git_commit()
            print("  ✓ Created git commit")

    def rollback(self) -> None:
        """Restore backups on failure."""
        if not self.backups:
            return

        print()
        print("Rolling back changes...")
        for backup_path in self.backups:
            original_path = backup_path.with_suffix("")
            if backup_path.exists():
                shutil.copy2(backup_path, original_path)
                backup_path.unlink()
                print(f"  ✓ Restored {original_path.name}")

    @staticmethod
    def _is_valid_version(version: str) -> bool:
        """Check if version follows semantic versioning."""
        # Semantic versioning pattern: X.Y.Z or X.Y.Z-suffix
        pattern = r"^\d+\.\d+\.\d+(?:-[a-zA-Z0-9.]+)?$"
        return bool(re.match(pattern, version))

    def _has_uncommitted_changes(self) -> bool:
        """Check if git has uncommitted changes."""
        try:
            result = subprocess.run(
                ["git", "status", "--porcelain"],
                cwd=self.root,
                capture_output=True,
                text=True,
                check=True,
            )
            return bool(result.stdout.strip())
        except (subprocess.CalledProcessError, FileNotFoundError):
            # If git is not available, assume clean
            return False

    def get_current_version(self) -> str:
        """Get the current version from pyproject.toml."""
        versions = self._get_all_current_versions()
        return versions.get("pyproject.toml", "unknown")

    def _get_all_current_versions(self) -> dict[str, str]:
        """Read current versions from all files."""
        versions: dict[str, str] = {}

        # pyproject.toml
        content = self.pyproject_toml.read_text()
        match = re.search(r'^version\s*=\s*"([^"]+)"', content, re.MULTILINE)
        if match:
            versions["pyproject.toml"] = match.group(1)

        # __init__.py
        content = self.init_py.read_text()
        match = re.search(r'^__version__\s*=\s*"([^"]+)"', content, re.MULTILINE)
        if match:
            versions["__init__.py"] = match.group(1)

        # CITATION.cff
        content = self.citation_cff.read_text()
        match = re.search(r"^version:\s*(.+)$", content, re.MULTILINE)
        if match:
            versions["CITATION.cff"] = match.group(1).strip()

        return versions

    def _create_backups(self) -> None:
        """Create .bak copies of all files before modification."""
        for file_path in [self.pyproject_toml, self.init_py, self.citation_cff]:
            backup_path = file_path.with_suffix(file_path.suffix + ".bak")
            shutil.copy2(file_path, backup_path)
            self.backups.append(backup_path)

    def _clean_backups(self) -> None:
        """Remove backup files."""
        for backup_path in self.backups:
            if backup_path.exists():
                backup_path.unlink()
        self.backups.clear()

    def _update_pyproject_toml(self) -> None:
        """Update version in pyproject.toml."""
        if has_tomli_w:
            content = self.pyproject_toml.read_text()
            data: dict[str, Any] = tomli.loads(content)  # type: ignore[possibly-undefined]
            data["project"]["version"] = self.new_version
            self.pyproject_toml.write_text(tomli_w.dumps(data))  # type: ignore[possibly-undefined]
        else:
            # Fallback to regex
            content = self.pyproject_toml.read_text()
            updated = re.sub(
                r'^(version\s*=\s*)"[^"]+"',
                rf'\1"{self.new_version}"',
                content,
                flags=re.MULTILINE,
            )
            self.pyproject_toml.write_text(updated)

    def _update_init_py(self) -> None:
        """Update __version__ in __init__.py."""
        content = self.init_py.read_text()
        updated = re.sub(
            r'^(__version__\s*=\s*)"[^"]+"',
            rf'\1"{self.new_version}"',
            content,
            flags=re.MULTILINE,
        )
        self.init_py.write_text(updated)

    def _update_citation_cff(self, date: str) -> None:
        """Update version and date in CITATION.cff."""
        if has_ruamel_yaml:
            # Use YAML parser (safer, preserves formatting)
            yaml = YAML()  # type: ignore[possibly-undefined]
            yaml.preserve_quotes = True
            yaml.default_flow_style = False
            data: dict[str, Any] = yaml.load(self.citation_cff)  # type: ignore[assignment]
            data["version"] = self.new_version
            data["date-released"] = date
            yaml.dump(data, self.citation_cff)
        else:
            # Fallback to regex
            content = self.citation_cff.read_text()
            updated = re.sub(
                r"^version:\s*.+$",
                f"version: {self.new_version}",
                content,
                flags=re.MULTILINE,
            )
            updated = re.sub(
                r"^date-released:\s*.+$",
                f"date-released: {date}",
                updated,
                flags=re.MULTILINE,
            )
            self.citation_cff.write_text(updated)

    def _regenerate_uv_lock(self) -> None:
        """Regenerate uv.lock file via uv lock command.

        Raises
        ------
        VersionBumpError
            If uv lock command fails or uv is not installed.

        """
        try:
            subprocess.run(
                ["uv", "lock"],
                cwd=self.root,
                check=True,
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as e:
            msg = f"Failed to regenerate uv.lock: {e.stderr}"
            raise VersionBumpError(msg) from e
        except FileNotFoundError as e:
            msg = "uv command not found. Install uv: https://github.com/astral-sh/uv"
            raise VersionBumpError(msg) from e

    @staticmethod
    def _show_changes_summary(versions: dict[str, str]) -> None:
        """Show summary of what changed."""
        print()
        print("  Changes:")
        for file_name in versions:
            print(f"    - {file_name}")

    def _create_git_commit(self) -> None:
        """Create a git commit with the version update.

        Raises
        ------
        VersionBumpError
            If git commit creation fails.

        """
        try:
            # Stage all changes
            subprocess.run(
                ["git", "add", "-A"],
                cwd=self.root,
                check=True,
                capture_output=True,
            )

            # Create commit
            commit_msg = (
                f"chore: bump version to {self.new_version}\n\n"
                f"Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>"
            )
            subprocess.run(
                ["git", "commit", "-m", commit_msg],
                cwd=self.root,
                check=True,
                capture_output=True,
            )

        except subprocess.CalledProcessError as e:
            msg = f"Failed to create git commit: {e}"
            raise VersionBumpError(msg) from e


def main() -> None:
    """Run the version bump script."""
    parser = argparse.ArgumentParser(
        description="Update version numbers across project files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/bump_version.py                # Interactive mode
  python scripts/bump_version.py 0.3.0          # Direct mode
  python scripts/bump_version.py 0.3.0 --commit # With git commit
  python scripts/bump_version.py 0.3.0 --dry-run # Preview changes
        """,
    )

    parser.add_argument(
        "version",
        nargs="?",
        help="New version number (semantic versioning format, e.g., 0.3.0)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would change without modifying files",
    )
    parser.add_argument(
        "--allow-dirty",
        action="store_true",
        help="Allow uncommitted changes in git",
    )
    parser.add_argument(
        "--commit",
        action="store_true",
        help="Create a git commit after updating",
    )

    args = parser.parse_args()

    # Interactive mode if no version provided
    if args.version is None:
        print("Interactive Version Bump")
        print("=" * 60)
        print()
        current = VersionBumper("0.0.0", dry_run=True).get_current_version()
        print(f"Current version: {current}")
        print()
        version = input("Enter new version (e.g., 0.3.0): ").strip()
        if not version:
            print("Error: Version required")
            sys.exit(1)
        args.version = version

    # Run version bump
    bumper = VersionBumper(
        args.version,
        dry_run=args.dry_run,
        allow_dirty=args.allow_dirty,
        commit=args.commit,
    )
    bumper.run()


if __name__ == "__main__":
    main()
