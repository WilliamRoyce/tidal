# GitHub Issue Tracking Setup

This document explains the issue tracking system created for the TIDAL project based on comprehensive codebase analysis.

## Overview

> **Update (v0.3.1):** Since this document was created, 17 of 25 identified issues have been resolved. Phases 1, 3, and most of Phase 4 are complete. The project has also been renamed to TIDAL, gained a CLI tool, professional branding, and 743 Python tests. See [ROADMAP.md](ROADMAP.md) for current status of all issues.

**Date Created:** February 2026
**Total Issues Identified:** 25 (17 resolved, 8 remaining)
**Source:** Comprehensive exploration of Wolfram pipeline, Python simulation code, tests, examples, documentation, and CI/CD infrastructure

## What Was Created

### 1. Issue Templates (`.github/ISSUE_TEMPLATE/`)

Standardized templates for contributors:

- **[bug_report.md](.github/ISSUE_TEMPLATE/bug_report.md)** - Template for bug reports
- **[feature_request.md](.github/ISSUE_TEMPLATE/feature_request.md)** - Template for feature requests
- **[config.yml](.github/ISSUE_TEMPLATE/config.yml)** - Issue template configuration

### 2. Issue Creation Script (`scripts/create_github_issues.py`)

Python script containing all 25 identified issues, ready to be created on GitHub.

**Features:**
- Complete issue definitions with descriptions, labels, and priorities
- Dry-run mode to preview issues before creation
- Automatic label assignment
- Detailed problem descriptions with code examples

### 3. Project Documentation

- **[CONTRIBUTING.md](../CONTRIBUTING.md)** - Updated contributor guidelines (already existed)
- **[docs/ROADMAP.md](ROADMAP.md)** - Public roadmap with prioritized phases
- **[This file]** - Setup and usage instructions

## Issue Summary by Priority

### 🔴 Critical (6 issues) — 4 resolved
**Security, correctness, or major user impact**

1. ✅ Replace Assertions with Explicit Error Handling
2. ✅ Strengthen Mathematica Expression Evaluation Security
3. 🔄 Add Wolfram Tests to GitHub Actions CI
4. ✅ Support Rank-3+ Tensor Decomposition
5. ✅ Add 3+1D Spacetime Examples
6. 🔄 Implement Automatic Gauge Fixing

### 🟠 High Priority (5 issues) — 3 resolved
**Important features or significant gaps**

7. 🔄 Add Animation Module Test Coverage
8. ⚡ Add Code Coverage Reporting to CI (partial: runs in CI, codecov.io pending)
9. ✅ Validate Grid Dimensions During PDE Construction
10. ✅ Improve Coefficient Resolution Performance
11. ✅ Document JSON Schema with Detailed Guide

### 🟡 Medium Priority (7 issues) — 4 resolved
**Improvements that enhance robustness**

12. 🔄 Add Tests for Observers, Profiling, Runners Modules
13. ✅ Handle Mixed Time-Space Cross-Derivatives Properly
14. ✅ Expand _mathematica_to_python Function Set
15. ✅ Add Non-Cartesian Coordinate System Examples
16. ✅ Add Convergence and Stability Stress Tests
17. 🔄 Create Architecture Diagrams
18. 🔄 Add Full Pipeline Validation to CI

### 🟢 Low Priority (7 issues) — 2 resolved
**Nice-to-have enhancements**

19. ✅ Refactor ContainsTimeDerivative and IsMixedTimeSpaceDerivative
20. 🔄 Add Parameter Sweep Examples
21. 🔄 Add Python 3.12+ Testing to CI Matrix
22. 🔄 Add PGF Export Module Tests
23. 🔄 Cache Coordinate Symbols and Chart Dimension
24. 🔄 Add Debugging and Performance Tuning Guides
25. ✅ Support Elliptic PDE Solving (Constraint Equations)

## How to Create Issues on GitHub

### Option 1: Using the Script (Recommended)

**Prerequisites:**
- GitHub CLI (`gh`) installed and authenticated
- Repository access

**Step 1: Install GitHub CLI**
```bash
# On Ubuntu/Debian
sudo apt install gh

# On macOS
brew install gh

# Authenticate
gh auth login
```

**Step 2: Preview Issues (Dry Run)**
```bash
python scripts/create_github_issues.py --dry-run
```

This shows what will be created without actually creating issues.

**Step 3: Create All Issues**
```bash
python scripts/create_github_issues.py
```

This will create all 25 issues with proper labels and descriptions.

### Option 2: Manual Creation

If you prefer to create issues manually or selectively:

1. Navigate to: https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/new
2. Open `scripts/create_github_issues.py` in an editor
3. Find the issue you want to create in the `ISSUES` list
4. Copy the title, body, and labels
5. Paste into the GitHub issue form
6. Submit

### Option 3: Selective Creation

To create only specific issues, edit `scripts/create_github_issues.py`:

```python
# At the bottom of ISSUES list, keep only the ones you want
ISSUES = [
    ISSUES[0],  # Issue #1
    ISSUES[6],  # Issue #7
    # etc.
]
```

Then run the script normally.

## Labels to Create First

Before running the script, ensure these labels exist in the repository:

### Priority Labels
- `priority: critical` (color: #d73a4a)
- `priority: high` (color: #e99695)
- `priority: medium` (color: #fbca04)
- `priority: low` (color: #0e8a16)

### Type Labels
- `bug` (color: #d73a4a)
- `enhancement` (color: #a2eeef)
- `documentation` (color: #0075ca)
- `testing` (color: #0e8a16)
- `ci-cd` (color: #bfdadc)
- `performance` (color: #e99695)
- `security` (color: #ee0701)
- `refactoring` (color: #fef2c0)
- `validation` (color: #1d76db)

### Component Labels
- `wolfram` (color: #fbca04)
- `python` (color: #f9d0c4)
- `examples` (color: #c5def5)
- `animation` (color: #bfdadc)

**To create labels via CLI:**
```bash
gh label create "priority: critical" --color d73a4a
gh label create "priority: high" --color e99695
gh label create "priority: medium" --color fbca04
gh label create "priority: low" --color 0e8a16
# ... etc for all labels
```

## Implementation Phases

See [ROADMAP.md](ROADMAP.md) for detailed phase breakdown:

1. **Phase 1:** ✅ COMPLETE — assertions, security, validation
2. **Phase 2:** 🔄 PRIMARY REMAINING FOCUS — coverage, animation tests, Wolfram CI
3. **Phase 3:** ✅ COMPLETE — 3+1D examples, JSON schema guide, coefficient resolution
4. **Phase 4:** ✅ MOSTLY COMPLETE — rank-3 tensors, mixed derivatives, curvilinear. Remaining: gauge fixing
5. **Phase 5:** Partially complete — derivative refactor done, elliptic solving done

## Managing Issues

### Project Board (Recommended)

Create a GitHub Project board to track progress:

1. Go to: https://github.com/WilliamRoyce/torsion-gertsenshtein/projects
2. Create new project: "Pipeline Improvements"
3. Add columns: `Backlog`, `Phase 1`, `Phase 2`, `Phase 3`, `Phase 4`, `Phase 5`, `In Progress`, `Done`
4. Add issues to appropriate phase columns

### Milestones

Create milestones for each phase:

```bash
gh milestone create "Phase 1: Critical Fixes" --due-date 2026-03-01
gh milestone create "Phase 2: Testing & CI" --due-date 2026-04-01
gh milestone create "Phase 3: Features & Docs" --due-date 2026-05-01
# etc.
```

Then assign issues to milestones:

```bash
gh issue edit <issue-number> --milestone "Phase 1: Critical Fixes"
```

## Integration with Development

### When Working on Issues

1. **Claim the issue:** Comment "I'm working on this"
2. **Create a branch:** `git checkout -b fix/issue-<number>-description`
3. **Follow CONTRIBUTING.md:** Ensure tests, docs, and code quality
4. **Reference in commits:** Use `Fixes #<number>` or `Closes #<number>` in PR description
5. **Link PR to issue:** GitHub will auto-link when you use keywords

### PR Review Checklist

Before marking issue as complete, ensure:
- [ ] All acceptance criteria met
- [ ] Tests added/updated
- [ ] Documentation updated
- [ ] CHANGELOG.md updated
- [ ] CI passes
- [ ] Code reviewed and approved

## Tracking Progress

### View by Priority
```bash
gh issue list --label "priority: critical"
gh issue list --label "priority: high"
```

### View by Phase
```bash
gh issue list --milestone "Phase 1: Critical Fixes"
```

### View by Type
```bash
gh issue list --label "bug"
gh issue list --label "enhancement"
```

## Updating the Roadmap

As issues are completed or priorities change:

1. Update [ROADMAP.md](ROADMAP.md) with current status
2. Move issues between phases in project board
3. Adjust milestone due dates if needed
4. Communicate changes in project discussions

## Questions or Issues?

- **Script problems:** Open issue with `ci-cd` label
- **Priority disputes:** Discuss in issue comments or project discussions
- **New issues:** Follow templates in `.github/ISSUE_TEMPLATE/`

---

## Appendix: Complete Issue List

See `scripts/create_github_issues.py` for full issue definitions with:
- Detailed problem descriptions
- Code examples
- Implementation notes
- Impact assessments
- File locations

**To view all issues without creating:**
```bash
python scripts/create_github_issues.py --dry-run | less
```

---

*This tracking system is based on comprehensive codebase exploration identifying 25 improvement areas. The project is in excellent shape (v0.3.1, 743 Python tests + ~108 Wolfram tests, 17/25 issues resolved). Remaining work is primarily testing/CI infrastructure and gauge automation.*
