"""Tests for version bumping script."""

# Import the module - adjust path if needed
import sys
import tempfile
from collections.abc import Generator
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

import scripts.bump_version  # pyright: ignore[reportMissingTypeStubs]
from scripts.bump_version import (  # pyright: ignore[reportMissingTypeStubs]
    VersionBumper,
    VersionBumpError,
)


class TestVersionValidation:
    """Test version format validation."""

    def test_valid_versions(self) -> None:
        """Test that valid semantic versions are accepted."""
        valid_versions = [
            "0.0.0",
            "0.1.0",
            "1.0.0",
            "1.2.3",
            "10.20.30",
            "1.0.0-alpha",
            "1.0.0-alpha.1",
            "1.0.0-beta.2",
            "1.0.0-rc.1",
        ]
        bumper = VersionBumper("0.0.0", dry_run=True)
        for version in valid_versions:
            assert bumper._is_valid_version(version), f"Should accept {version}"

    def test_invalid_versions(self) -> None:
        """Test that invalid versions are rejected."""
        invalid_versions = [
            "1",
            "1.0",
            "1.0.0.0",
            "v1.0.0",
            "1.0.0 alpha",
            "abc",
            "1.0.x",
            "",
        ]
        bumper = VersionBumper("0.0.0", dry_run=True)
        for version in invalid_versions:
            assert not bumper._is_valid_version(version), f"Should reject {version}"


class TestFileUpdates:
    """Test file update functions with temporary files."""

    @pytest.fixture
    def temp_dir(self) -> Generator[Path, None, None]:
        """Create a temporary directory for testing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def mock_project(self, temp_dir: Path) -> Path:
        """Create a mock project structure."""
        # Create directories
        (temp_dir / "tidal").mkdir()

        # Create pyproject.toml
        pyproject_content = """[project]
name = "tidal"
version = "0.1.0"
description = "Test project"
"""
        (temp_dir / "pyproject.toml").write_text(pyproject_content)

        # Create __init__.py
        init_content = '''"""Test package."""

__version__ = "0.1.0"
'''
        (temp_dir / "tidal" / "__init__.py").write_text(init_content)

        # Create CITATION.cff
        citation_content = """cff-version: 1.2.0
title: "Test Project"
version: 0.1.0
date-released: 2026-01-01
"""
        (temp_dir / "CITATION.cff").write_text(citation_content)

        return temp_dir

    def test_get_current_versions(
        self, mock_project: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Test reading current versions from files."""
        # Monkeypatch the root directory
        bumper = VersionBumper("0.2.0", dry_run=True)
        monkeypatch.setattr(bumper, "root", mock_project)
        monkeypatch.setattr(bumper, "pyproject_toml", mock_project / "pyproject.toml")
        monkeypatch.setattr(
            bumper, "init_py", mock_project / "tidal" / "__init__.py"
        )
        monkeypatch.setattr(bumper, "citation_cff", mock_project / "CITATION.cff")

        versions = bumper._get_all_current_versions()

        assert versions["pyproject.toml"] == "0.1.0"
        assert versions["__init__.py"] == "0.1.0"
        assert versions["CITATION.cff"] == "0.1.0"

    def test_update_init_py(
        self, mock_project: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Test updating __init__.py version."""
        bumper = VersionBumper("0.2.0", dry_run=False)
        init_path = mock_project / "tidal" / "__init__.py"
        monkeypatch.setattr(bumper, "init_py", init_path)

        bumper._update_init_py()

        content = init_path.read_text()
        assert '__version__ = "0.2.0"' in content

    def test_update_pyproject_toml_regex(
        self, mock_project: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Test updating pyproject.toml with regex (no tomli_w)."""
        bumper = VersionBumper("0.2.0", dry_run=False)
        pyproject_path = mock_project / "pyproject.toml"
        monkeypatch.setattr(bumper, "pyproject_toml", pyproject_path)

        # Force regex path by disabling tomli_w

        original_has_tomli: bool = scripts.bump_version.has_tomli_w
        scripts.bump_version.has_tomli_w = False

        try:
            bumper._update_pyproject_toml()
            content = pyproject_path.read_text()
            assert 'version = "0.2.0"' in content
        finally:
            scripts.bump_version.has_tomli_w = original_has_tomli

    def test_update_citation_cff_regex(
        self, mock_project: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Test updating CITATION.cff with regex (no ruamel.yaml)."""
        bumper = VersionBumper("0.2.0", dry_run=False)
        citation_path = mock_project / "CITATION.cff"
        monkeypatch.setattr(bumper, "citation_cff", citation_path)

        # Force regex path by disabling ruamel.yaml

        original_has_ruamel: bool = scripts.bump_version.has_ruamel_yaml
        scripts.bump_version.has_ruamel_yaml = False

        try:
            bumper._update_citation_cff("2026-02-06")
            content = citation_path.read_text()
            assert "version: 0.2.0" in content
            # Date might have quotes depending on parser
            assert (
                "date-released: 2026-02-06" in content
                or "date-released: '2026-02-06'" in content
            )
        finally:
            scripts.bump_version.has_ruamel_yaml = original_has_ruamel


class TestBackupRestore:
    """Test backup and restore mechanism."""

    @pytest.fixture
    def temp_dir(self) -> Generator[Path, None, None]:
        """Create a temporary directory for testing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_backup_creation(
        self, temp_dir: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Test that backups are created correctly."""
        # Create test files
        test_file = temp_dir / "test.txt"
        test_file.write_text("original content")

        bumper = VersionBumper("0.2.0", dry_run=False)
        monkeypatch.setattr(bumper, "pyproject_toml", test_file)
        monkeypatch.setattr(bumper, "init_py", test_file)
        monkeypatch.setattr(bumper, "citation_cff", test_file)

        bumper._create_backups()

        # Check backup exists
        backup_file = test_file.with_suffix(".txt.bak")
        assert backup_file.exists()
        assert backup_file.read_text() == "original content"

    def test_backup_restore(
        self, temp_dir: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Test that backups are restored on rollback."""
        # Create test files
        test_file = temp_dir / "test.txt"
        test_file.write_text("original content")

        bumper = VersionBumper("0.2.0", dry_run=False)
        monkeypatch.setattr(bumper, "pyproject_toml", test_file)
        monkeypatch.setattr(bumper, "init_py", test_file)
        monkeypatch.setattr(bumper, "citation_cff", test_file)

        bumper._create_backups()

        # Modify the original file
        test_file.write_text("modified content")

        # Rollback
        bumper.rollback()

        # Check original content restored
        assert test_file.read_text() == "original content"

        # Check backup removed
        backup_file = test_file.with_suffix(".txt.bak")
        assert not backup_file.exists()


class TestErrorConditions:
    """Test error handling."""

    def test_invalid_version_raises_error(self) -> None:
        """Test that invalid version raises error during validation."""
        bumper = VersionBumper("invalid-version", dry_run=False, allow_dirty=True)

        with pytest.raises(VersionBumpError, match="Invalid version format"):
            bumper.validate()

    def test_missing_file_raises_error(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Test that missing file raises error during validation."""
        bumper = VersionBumper("0.2.0", dry_run=False, allow_dirty=True)
        monkeypatch.setattr(bumper, "root", tmp_path)
        monkeypatch.setattr(bumper, "pyproject_toml", tmp_path / "nonexistent.toml")
        monkeypatch.setattr(bumper, "init_py", tmp_path / "nonexistent.py")
        monkeypatch.setattr(bumper, "citation_cff", tmp_path / "nonexistent.cff")

        with pytest.raises(VersionBumpError, match="Required file not found"):
            bumper.validate()
