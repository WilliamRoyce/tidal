#!/bin/bash
# Script to create GitHub labels for issue tracking
# Requires: gh CLI installed and authenticated

set -e

echo "Creating GitHub labels for torsion-gertsenshtein..."
echo ""

# Priority labels
echo "📌 Creating priority labels..."
gh label create "priority: critical" --color "d73a4a" --description "🔴 Security, correctness, or major user impact" --force || true
gh label create "priority: high" --color "e99695" --description "🟠 Important features or significant gaps" --force || true
gh label create "priority: medium" --color "fbca04" --description "🟡 Improvements that enhance robustness" --force || true
gh label create "priority: low" --color "0e8a16" --description "🟢 Nice-to-have enhancements" --force || true

# Type labels
echo "📝 Creating type labels..."
gh label create "bug" --color "d73a4a" --description "Something isn't working" --force || true
gh label create "enhancement" --color "a2eeef" --description "New feature or request" --force || true
gh label create "documentation" --color "0075ca" --description "Improvements or additions to documentation" --force || true
gh label create "testing" --color "0e8a16" --description "Test additions or improvements" --force || true
gh label create "ci-cd" --color "bfdadc" --description "Continuous integration/deployment" --force || true
gh label create "performance" --color "e99695" --description "Performance optimization" --force || true
gh label create "security" --color "ee0701" --description "Security vulnerability" --force || true
gh label create "refactoring" --color "fef2c0" --description "Code quality improvement" --force || true
gh label create "validation" --color "1d76db" --description "Input/output validation" --force || true
gh label create "type-safety" --color "1d76db" --description "Type checking and safety" --force || true

# Component labels
echo "🔧 Creating component labels..."
gh label create "wolfram" --color "fbca04" --description "Mathematica/xAct pipeline" --force || true
gh label create "python" --color "f9d0c4" --description "Python simulation code" --force || true
gh label create "examples" --color "c5def5" --description "Example additions/improvements" --force || true
gh label create "animation" --color "bfdadc" --description "Visualization features" --force || true
gh label create "tensors" --color "fbca04" --description "Tensor field operations" --force || true
gh label create "gauge-theory" --color "fbca04" --description "Gauge theory features" --force || true
gh label create "3d" --color "c5def5" --description "3D spacetime features" --force || true
gh label create "json-schema" --color "0075ca" --description "JSON schema and validation" --force || true
gh label create "derivatives" --color "fbca04" --description "Derivative operators" --force || true
gh label create "mathematica" --color "fbca04" --description "Mathematica expression handling" --force || true
gh label create "coordinates" --color "c5def5" --description "Coordinate systems" --force || true
gh label create "convergence" --color "0e8a16" --description "Convergence and stability" --force || true
gh label create "architecture" --color "0075ca" --description "System architecture" --force || true
gh label create "code-quality" --color "fef2c0" --description "Code quality improvements" --force || true
gh label create "future-proofing" --color "bfdadc" --description "Future compatibility" --force || true
gh label create "visualization" --color "bfdadc" --description "Visualization and plotting" --force || true
gh label create "elliptic-pde" --color "a2eeef" --description "Elliptic PDE support" --force || true
gh label create "optimization" --color "e99695" --description "Performance optimization" --force || true

# Special labels
echo "⭐ Creating special labels..."
gh label create "roadmap" --color "0052cc" --description "Part of project roadmap" --force || true
gh label create "good first issue" --color "7057ff" --description "Good for newcomers" --force || true
gh label create "help wanted" --color "008672" --description "Extra attention is needed" --force || true

echo ""
echo "✅ Done! Created/updated all labels."
echo ""
echo "Next steps:"
echo "1. Run: python scripts/create_github_issues.py --dry-run (preview)"
echo "2. Run: python scripts/create_github_issues.py (create issues)"
