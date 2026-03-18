## Description

<!-- Provide a brief summary of the changes -->

## Related Issue

<!-- Link to the issue this PR addresses, if applicable -->
Fixes #(issue number)

## Type of Change

<!-- Mark the relevant option with an "x" -->

- [ ] Bug fix (non-breaking change that fixes an issue)
- [ ] New feature (non-breaking change that adds functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to not work as expected)
- [ ] Documentation update
- [ ] Physics validation / equation correction
- [ ] Refactoring (no functional changes)

## Changes Made

<!-- Describe the changes in detail -->

-
-
-

## Physics Context

<!-- If relevant, describe the physics/equations involved -->

-

## Testing

<!-- Describe the tests you ran and their results -->

- [ ] All existing tests pass (`uv run pytest tests/`)
- [ ] New tests added for new functionality
- [ ] Wolfram tests pass (if applicable)
- [ ] Manual testing performed

## Checklist

<!-- Ensure all items are checked before submitting -->

- [ ] Code follows project conventions (see CONTRIBUTING.md)
- [ ] **No hardcoded physics equations** - all PDEs derive from Lagrangian
- [ ] Type hints added to all functions (`-> None` for tests)
- [ ] Docstrings added/updated (NumPy style)
- [ ] Documentation updated (docstrings, relevant `docs/*.md` if needed)
- [ ] Ruff linting passes (`uv run ruff check`)
- [ ] Pyright type checking passes (`uv run pyright`)
- [ ] Commit messages are descriptive
- [ ] I have read and follow the [Code of Conduct](../CODE_OF_CONDUCT.md)

## Additional Notes

<!-- Any additional information reviewers should know -->
