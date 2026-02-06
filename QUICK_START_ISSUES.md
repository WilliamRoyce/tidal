# Quick Start: Creating GitHub Issues

## Prerequisites ✅
- [x] GitHub CLI installed (gh 2.86.0)
- [ ] Authenticated with GitHub
- [ ] Labels created
- [ ] Issues created

## Commands (Copy-Paste Ready)

### 1. Authenticate
```bash
gh auth login
```
Follow prompts → Choose GitHub.com → HTTPS → Browser authentication

### 2. Create Labels
```bash
./scripts/create_github_labels.sh
```

### 3. Preview Issues (Optional)
```bash
python scripts/create_github_issues.py --dry-run | less
```

### 4. Create All Issues
```bash
python scripts/create_github_issues.py
```

---

## What Gets Created

**25 Issues** across 5 phases:

### 🔴 Critical (6)
1. Replace assertions with explicit errors
2. Strengthen eval() validation (NaN/Inf)
3. Add Wolfram tests to CI
4. Support rank-3+ tensors
5. Add 3+1D examples
6. Implement automatic gauge fixing

### 🟠 High (5)
7. Animation module tests (0 coverage)
8. Code coverage reporting
9. Grid dimension validation at construction
10. Coefficient resolution performance
11. JSON schema documentation

### 🟡 Medium (7)
12-18. Testing, features, validation, CI improvements

### 🟢 Low (7)
19-25. Code quality, optimization, future-proofing

---

## After Issues Are Created

1. **Set up project board** (optional):
   - Go to: https://github.com/WilliamRoyce/torsion-gertsenshtein/projects
   - Create: "Pipeline Improvements"
   - Add columns: Phase 1, Phase 2, Phase 3, Phase 4, Phase 5, In Progress, Done

2. **Create milestones** (optional):
   ```bash
   gh milestone create "Phase 1: Critical Fixes" --due-date 2026-03-01
   gh milestone create "Phase 2: Testing & CI" --due-date 2026-04-01
   ```

3. **Start with Phase 1** - Critical fixes (Issues #1, #2, #9)

---

## Documentation

| File | Purpose |
|------|---------|
| [ISSUE_TRACKING_SUMMARY.md](ISSUE_TRACKING_SUMMARY.md) | Quick reference, all 25 issues |
| [docs/ISSUE_TRACKING_SETUP.md](docs/ISSUE_TRACKING_SETUP.md) | Complete setup guide |
| [docs/ROADMAP.md](docs/ROADMAP.md) | 5-phase plan with timeline |
| [scripts/create_github_issues.py](scripts/create_github_issues.py) | Full issue definitions |

---

## Troubleshooting

**"gh: command not found"** → Restart terminal or run `source ~/.bashrc`

**"gh auth: not authenticated"** → Run `gh auth login` first

**"label already exists"** → Safe to ignore (script uses `--force`)

**Preview issues taking too long** → Use `head -200` instead of `less`:
```bash
python scripts/create_github_issues.py --dry-run | head -200
```

---

**Status:** GitHub CLI installed ✅ | Ready to authenticate and create issues!
