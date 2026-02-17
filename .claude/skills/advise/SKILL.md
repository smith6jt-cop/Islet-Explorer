---
name: advise
description: Search the Skills Registry for relevant experiments before starting new work
---

# /advise — Skills Registry Lookup

Search the Skills Registry for relevant experiments before starting new work.

## Steps

1. Read the user's goal from the conversation context
2. Search `Skills_Registry/plugins/` for related skills by:
   - Scanning `description` fields in each `plugin.json` (look for matching trigger conditions)
   - Reading `SKILL.md` files in matching plugins for detailed findings
   - Check all categories: `general/`, `scientific/`, `trading/`, `kintsugi/`, `maxfuse/`, `training/`
3. For each relevant skill found, summarize:
   - What was tried and what worked (from "Verified Workflow" section)
   - What failed and why (from "Failed Attempts" table — the most valuable section)
   - Recommended parameters (from "Final Parameters" section)
   - Environment/version notes
4. If no relevant skills found:
   - Inform the user that no prior experiments match their goal
   - Suggest running `/retrospective` after completing their task to capture learnings

## Where to Search

The Skills Registry is a git submodule at `Skills_Registry/`:

```
Skills_Registry/plugins/
  ├── general/              # Cross-project Python/dev skills
  ├── scientific/           # Scientific computing patterns
  ├── trading/              # Alpaca Trading system skills
  ├── kintsugi/             # KINTSUGI-specific skills
  ├── maxfuse/              # MaxFuse imaging integration
  ├── training/             # Training experiment documentation
  └── templates/            # Skill templates
```

Each plugin has:
- `.claude-plugin/plugin.json` — metadata with trigger conditions in `description`
- `skills/<name>/SKILL.md` — detailed documentation with verified workflows and failed attempts

You can also consult `Skills_Registry/marketplace.json` for a quick index of all plugins.
