"""Tests for version bumping script."""

from __future__ import annotations

# Import the module - adjust path if needed
import sys
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from collections.abc import Generator

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


class TestComputeNextVersion:
    """Test automatic version computation."""

    def test_patch_bump(self) -> None:
        """Test patch version increment."""
        assert VersionBumper.compute_next_version("0.5.0", "patch") == "0.5.1"
        assert VersionBumper.compute_next_version("0.5.3", "patch") == "0.5.4"
        assert VersionBumper.compute_next_version("1.2.9", "patch") == "1.2.10"

    def test_minor_bump(self) -> None:
        """Test minor version increment resets patch to 0."""
        assert VersionBumper.compute_next_version("0.5.0", "minor") == "0.6.0"
        assert VersionBumper.compute_next_version("0.5.3", "minor") == "0.6.0"
        assert VersionBumper.compute_next_version("0.9.0", "minor") == "0.10.0"
        assert VersionBumper.compute_next_version("0.99.5", "minor") == "0.100.0"

    def test_strips_prerelease(self) -> None:
        """Test that pre-release suffixes are stripped before bumping."""
        assert VersionBumper.compute_next_version("1.0.0-alpha.1", "patch") == "1.0.1"
        assert VersionBumper.compute_next_version("0.5.0-beta.2", "minor") == "0.6.0"

    def test_invalid_level_raises(self) -> None:
        """Test that unsupported bump level raises ValueError."""
        with pytest.raises(ValueError, match="Unsupported bump level"):
            VersionBumper.compute_next_version("0.5.0", "major")

    def test_high_minor_stays_same_major(self) -> None:
        """Test that minor bumps never change the major version."""
        # 0.9.0 → 0.10.0, not 1.0.0
        result = VersionBumper.compute_next_version("0.9.0", "minor")
        assert result == "0.10.0"
        assert result.split(".")[0] == "0"


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
        (temp_dir / "docs").mkdir()

        # Create pyproject.toml
        pyproject_content = """[project]
name = "tidal"
version = "0.1.0"
description = "Test project"
"""
        (temp_dir / "pyproject.toml").write_text(pyproject_content)

        # Create CITATION.cff
        citation_content = """cff-version: 1.2.0
title: "Test Project"
version: 0.1.0
date-released: 2026-01-01
"""
        (temp_dir / "CITATION.cff").write_text(citation_content)

        # Create docs/NEXT_PHASES.md
        next_phases_content = """# Next Phases
**Version:** 0.1.0 | **Tests:** 100 collected
"""
        (temp_dir / "docs" / "NEXT_PHASES.md").write_text(next_phases_content)

        # Create docs/ROADMAP.md
        roadmap_content = """# Roadmap
**Current Version:** 0.1.0
**Previous Milestones:** 0.0.1 delivered initial
"""
        (temp_dir / "docs" / "ROADMAP.md").write_text(roadmap_content)

        # Create SECURITY.md
        security_content = """# Security
| Version | Supported          |
| ------- | ------------------ |
| 0.1.x   | :white_check_mark: |
| < 0.1   | :x:                |
"""
        (temp_dir / "SECURITY.md").write_text(security_content)

        return temp_dir

    def test_get_current_versions(
        self,
        mock_project: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Test reading current versions from files."""
        bumper = VersionBumper("0.2.0", dry_run=True)
        monkeypatch.setattr(bumper, "root", mock_project)
        monkeypatch.setattr(bumper, "pyproject_toml", mock_project / "pyproject.toml")
        monkeypatch.setattr(bumper, "citation_cff", mock_project / "CITATION.cff")

        versions = bumper._get_all_current_versions()

        assert versions["pyproject.toml"] == "0.1.0"
        assert versions["CITATION.cff"] == "0.1.0"
        # __init__.py is no longer tracked (reads version dynamically)
        assert "__init__.py" not in versions

    def test_update_pyproject_toml_regex(
        self,
        mock_project: Path,
        monkeypatch: pytest.MonkeyPatch,
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
        self,
        mock_project: Path,
        monkeypatch: pytest.MonkeyPatch,
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

    def test_update_next_phases_md(
        self,
        mock_project: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Test updating version in docs/NEXT_PHASES.md."""
        bumper = VersionBumper("0.2.0", dry_run=False)
        next_phases_path = mock_project / "docs" / "NEXT_PHASES.md"
        monkeypatch.setattr(bumper, "next_phases_md", next_phases_path)

        bumper._update_next_phases_md()
        content = next_phases_path.read_text()
        assert "**Version:** 0.2.0" in content
        # Should NOT touch the Tests count
        assert "**Tests:** 100 collected" in content

    def test_update_roadmap_md(
        self,
        mock_project: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Test updating current version in docs/ROADMAP.md."""
        bumper = VersionBumper("0.2.0", dry_run=False)
        roadmap_path = mock_project / "docs" / "ROADMAP.md"
        monkeypatch.setattr(bumper, "roadmap_md", roadmap_path)

        bumper._update_roadmap_md()
        content = roadmap_path.read_text()
        assert "**Current Version:** 0.2.0" in content
        # Should NOT touch previous milestones history
        assert "0.0.1 delivered initial" in content

    def test_update_security_md(
        self,
        mock_project: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Test updating supported version table in SECURITY.md."""
        bumper = VersionBumper("0.6.0", dry_run=False)
        security_path = mock_project / "SECURITY.md"
        monkeypatch.setattr(bumper, "security_md", security_path)

        bumper._update_security_md()
        content = security_path.read_text()
        assert "0.6.x" in content
        assert "< 0.6" in content
        assert ":white_check_mark:" in content
        assert ":x:" in content


class TestBackupRestore:
    """Test backup and restore mechanism."""

    @pytest.fixture
    def temp_dir(self) -> Generator[Path, None, None]:
        """Create a temporary directory for testing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_backup_creation(
        self,
        temp_dir: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Test that backups are created correctly."""
        # Create test files
        test_file = temp_dir / "test.txt"
        test_file.write_text("original content")

        bumper = VersionBumper("0.2.0", dry_run=False)
        monkeypatch.setattr(bumper, "pyproject_toml", test_file)
        monkeypatch.setattr(bumper, "citation_cff", test_file)
        # Set optional files to non-existent paths (should be skipped)
        monkeypatch.setattr(bumper, "next_phases_md", temp_dir / "nonexistent1.md")
        monkeypatch.setattr(bumper, "roadmap_md", temp_dir / "nonexistent2.md")
        monkeypatch.setattr(bumper, "security_md", temp_dir / "nonexistent3.md")

        bumper._create_backups()

        # Check backup exists
        backup_file = test_file.with_suffix(".txt.bak")
        assert backup_file.exists()
        assert backup_file.read_text() == "original content"

    def test_backup_restore(
        self,
        temp_dir: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Test that backups are restored on rollback."""
        # Create test files
        test_file = temp_dir / "test.txt"
        test_file.write_text("original content")

        bumper = VersionBumper("0.2.0", dry_run=False)
        monkeypatch.setattr(bumper, "pyproject_toml", test_file)
        monkeypatch.setattr(bumper, "citation_cff", test_file)
        monkeypatch.setattr(bumper, "next_phases_md", temp_dir / "nonexistent1.md")
        monkeypatch.setattr(bumper, "roadmap_md", temp_dir / "nonexistent2.md")
        monkeypatch.setattr(bumper, "security_md", temp_dir / "nonexistent3.md")

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
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Test that missing file raises error during validation."""
        bumper = VersionBumper("0.2.0", dry_run=False, allow_dirty=True)
        monkeypatch.setattr(bumper, "root", tmp_path)
        monkeypatch.setattr(bumper, "pyproject_toml", tmp_path / "nonexistent.toml")
        monkeypatch.setattr(bumper, "citation_cff", tmp_path / "nonexistent.cff")

        with pytest.raises(VersionBumpError, match="Required file not found"):
            bumper.validate()
