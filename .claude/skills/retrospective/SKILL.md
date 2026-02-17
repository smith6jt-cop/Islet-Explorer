---
name: retrospective
description: Save learnings from the current session as a new skill in the Skills Registry
---

# /retrospective — Save Session Learnings

Capture learnings from the current session as a new skill in the Skills Registry.

## Steps

1. **Summarize key findings** from the conversation:
   - What was the goal?
   - What approaches worked?
   - What failed and why? (most valuable)
   - What were the final parameters/configuration?

2. **Create a new plugin** using the template:
   - Copy `Skills_Registry/templates/experiment-skill-template/` to the appropriate category:
     - `plugins/general/` — cross-project Python/dev skills
     - `plugins/scientific/` — scientific computing patterns
     - `plugins/trading/` — Alpaca Trading system skills
     - `plugins/kintsugi/` — KINTSUGI-specific skills
   - Rename `TEMPLATE_NAME` to a descriptive kebab-case name

3. **Fill in plugin.json** (`.claude-plugin/plugin.json`):
   ```json
   {
     "name": "your-skill-name",
     "version": "1.0.0",
     "description": "Specific trigger conditions: (1) scenario one, (2) scenario two. Verified on [environment].",
     "author": { "name": "smith6jt" },
     "skills": "./skills",
     "repository": "https://github.com/smith-cop/Skills_Registry"
   }
   ```

4. **Fill in SKILL.md** (`skills/<name>/SKILL.md`) with YAML frontmatter and sections:
   - **Experiment Overview** — date, goal, environment, status
   - **Context** — problem description
   - **Verified Workflow** — step-by-step with exact commands/code
   - **Failed Attempts (Critical)** — table: approach, result, lesson learned
   - **Final Parameters** — copy-pasteable configurations
   - **Key Insights** — bullet points of important learnings

5. **Create a branch and open a PR** to the Skills Registry:
   ```bash
   cd Skills_Registry
   git checkout -b skill/your-skill-name
   git add plugins/category/your-skill-name/
   git commit -m "Add your-skill-name skill: brief description"
   git push -u origin skill/your-skill-name
   gh pr create --title "Add skill: your-skill-name" --body "Description of what was learned"
   ```

## Rules

- Every skill needs specific trigger conditions in the description (not vague advice)
- The "Failed Attempts" table is the most valuable section — always include it
- Include exact hyperparameters and configurations, not general guidance
- Document the environment (software versions, hardware) where verified
- Skills should be specific enough to be actionable but general enough to be reusable
