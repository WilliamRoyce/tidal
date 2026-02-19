# Documentation Directory

This directory contains living documentation that captures the project's evolution, design decisions, and hard-won lessons.

## 📚 Documentation Files

### [MEMORY.md](MEMORY.md)

**The main reference guide** - Start here for architecture, patterns, and critical decisions.

**Contents:**

- Critical architecture decisions (multi-field coupling, dimension handling)
- xAct/Wolfram patterns and best practices
- Pipeline module integration patterns
- Known issues and solutions
- Example implementations and their key features
- Quick reference (file locations, key functions)

**Use when:** Planning new features, debugging architectural issues, onboarding to the codebase

---

### [troubleshooting.md](troubleshooting.md)

**Error encyclopedia** - Look here when something breaks.

**Contents:**

- Common Wolfram/xAct errors and fixes
- Python/py-pde issues (operators, grids, state)
- Debugging techniques for both sides
- Verification checklists after changes

**Use when:** Debugging errors, after hitting a new issue (add it!), before making risky changes

---

### [chern-simons-notes.md](chern-simons-notes.md)

**Example-specific deep dive** - Reference for complex implementation patterns.

**Contents:**

- Physics background for Chern-Simons theory
- Implementation status (symbolic vs manual approaches)
- Wolfram hybrid approach details
- JSON structure for topological terms
- Python simulation specifics
- Future automation roadmap

**Use as template for:** New complex examples, topological theories, gauge theories, any case requiring hybrid symbolic+manual approaches

---

## 🔄 Maintenance Philosophy

These documents are **living references** that must evolve with the codebase:

1. **Update immediately** when solving non-trivial bugs
2. **Add patterns** after implementing new features
3. **Refine sections** when better approaches are discovered
4. **Create new example notes** for complex implementations
5. **Prune obsolete info** when code changes eliminate old issues

## 🎯 Document Relationships

```
MEMORY.md (Architecture & Patterns)
    ↓
    ├─→ troubleshooting.md (Error Solutions)
    └─→ chern-simons-notes.md (Example Details)
         └─→ [future: yang-mills-notes.md, etc.]
```

- **MEMORY.md** provides the "what and why" (design decisions, patterns)
- **troubleshooting.md** provides the "when it breaks" (error resolution)
- **Example notes** provide the "how for this case" (implementation specifics)

## 📝 Creating New Documentation

### When to create a new example-notes file:

Create `docs/{example-name}-notes.md` when:

- The example requires hybrid symbolic/manual approaches
- Special tensor structures need careful handling (epsilon, field strength, etc.)
- The implementation has interesting physics or mathematical subtleties
- You want to document a roadmap for future automation

Use `chern-simons-notes.md` as a template.

### What belongs where:

| Content Type          | Destination                                 |
| --------------------- | ------------------------------------------- |
| Architecture decision | MEMORY.md → Critical Architecture Decisions |
| General xAct pattern  | MEMORY.md → xAct/Wolfram Patterns           |
| Error you solved      | troubleshooting.md → Appropriate section    |
| Example physics       | {example}-notes.md → Physics Background     |
| Example-specific code | {example}-notes.md → Implementation Pattern |
| New file location     | MEMORY.md → Quick Reference                 |
| New function purpose  | MEMORY.md → Quick Reference                 |

---

## 🚀 Quick Start

**New to the project?** Read in this order:

1. Project README.md (root)
2. MEMORY.md (architecture overview)
3. Run an example: `examples/scalar_field/` or `examples/coupled_scalars/`
4. Skim troubleshooting.md to know what to watch for

**Implementing something new?**

1. Check MEMORY.md for existing patterns
2. Look at similar examples (scalar → coupled scalars → electromagnetic → chern-simons)
3. Add your findings back to these docs when done

**Hit an error?**

1. Search troubleshooting.md for symptoms
2. If not found, debug and **add your solution** to troubleshooting.md
3. If it reveals an architectural insight, update MEMORY.md too

---

Last updated: 2026-02-09 (Documentation comprehensive review and update)
