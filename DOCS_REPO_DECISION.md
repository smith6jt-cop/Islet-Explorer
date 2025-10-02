# Documentation Repository Strategy

## Decision: Keep Documentation in Main Repository ✅

### Summary

**Keep documentation with the main codebase.** The source files are small (308 KB), and this is the industry-standard approach for most projects.

## Size Analysis

```
docs/source/      308 KB   ✅ Include in git (Markdown/RST source files)
docs/build/        15 MB   ❌ Exclude from git (generated HTML)
docs_env/         110 MB   ❌ Exclude from git (Python virtual environment)
```

**What gets committed:** Only 69 lightweight text files (~308 KB total)

## Updated `.gitignore`

Changed from:
```
docs/              # ❌ Excluded everything
```

To:
```
docs/build/        # ✅ Only exclude generated files
docs/doctrees/     # ✅ Only exclude Sphinx cache
docs_env/          # ✅ Only exclude virtual environment
```

## Why This Approach is Best

### ✅ Advantages of Single Repo

1. **Tiny footprint** - 308 KB of Markdown is negligible
2. **Industry standard** - 95% of projects keep docs with code
3. **Atomic commits** - Code and docs update together
4. **Version sync** - Docs always match the code version
5. **Easier maintenance** - One repo to manage, not two
6. **Simpler CI/CD** - Read the Docs pulls from one place
7. **Better discoverability** - Users find docs immediately

### ❌ Why NOT Separate Repo

- **Extra complexity** - Need to sync versions between repos
- **Version drift risk** - Docs can get out of sync with code
- **Contributor friction** - PRs need updates in two places
- **CI/CD overhead** - Need to coordinate builds across repos
- **No size benefit** - Source files are already tiny

## Examples from Popular Projects

**Single Repo (Most Common):**
- Django (Python) - docs in `/docs`
- React (JavaScript) - docs in `/docs`
- Kubernetes (Go) - docs in `/docs`
- TensorFlow (Python) - docs in `/docs`
- Vue.js (JavaScript) - docs in `/docs`

**Separate Repo (Rare, for specific reasons):**
- Linux kernel - Massive project, many subsystems
- LLVM - Multiple interrelated projects
- (Usually only for multi-repo projects)

## What Gets Committed

**Included (69 files, 308 KB):**
```
docs/
├── Makefile                    # Build commands
├── README.md                   # How to build docs
├── build_docs.sh               # Automated build script
├── requirements.txt            # Python dependencies for docs
└── source/
    ├── conf.py                 # Sphinx configuration
    ├── index.rst               # Main landing page
    └── [67 Markdown files]     # All documentation content
```

**Excluded (automatically by .gitignore):**
```
docs/build/                     # Generated HTML (15 MB)
docs/doctrees/                  # Sphinx cache
docs_env/                       # Python virtual environment (110 MB)
```

## Read the Docs Integration

When you push to GitHub and connect Read the Docs:

1. **RTD clones your repo** - Gets `docs/source/` files only
2. **RTD creates its own venv** - Installs from `requirements.txt`
3. **RTD runs Sphinx** - Generates fresh HTML on their servers
4. **RTD hosts the result** - At `yourproject.readthedocs.io`

**Your local `docs/build/` is never used by Read the Docs!**

## Developer Workflow

```bash
# Clone repo - gets code + docs source
git clone https://github.com/smith6jt-cop/Islet-Explorer.git
cd Islet-Explorer

# Build docs locally (optional)
cd docs
./build_docs.sh

# View locally (optional)
cd build/html
python3 -m http.server 8000

# Commit changes to code + docs together
git add app/ docs/source/
git commit -m "Add new feature + documentation"
git push

# Read the Docs auto-rebuilds (no manual steps!)
```

## Comparison: Single vs Separate Repos

| Aspect | Single Repo (✅ Recommended) | Separate Repo |
|--------|------------------------------|---------------|
| **Complexity** | Simple | Complex |
| **Version sync** | Automatic | Manual |
| **Commits** | Atomic (code+docs) | Split across repos |
| **PR workflow** | Single PR | Two PRs |
| **CI/CD** | One build | Two builds |
| **Discoverability** | Immediate | Need to find other repo |
| **Size impact** | 308 KB | None (but unnecessary) |
| **Maintenance** | Low | Higher |
| **Industry practice** | ✅ Standard | Rare |

## When to Use Separate Repo

Only consider separate docs repo if:
- ❌ Docs are >100 MB of binary assets (not your case - 308 KB text)
- ❌ You have a multi-repo project (e.g., microservices)
- ❌ Docs are generated from multiple source repos
- ❌ You have dedicated docs team with different release cycle

**None of these apply to Islet-Explorer!**

## Decision Summary

✅ **Keep documentation in main repo**
- Adds only 308 KB of text files
- Standard industry practice
- Simpler maintenance
- Better version control
- Easier for contributors

❌ **Don't use separate repo**
- No benefit (docs are tiny)
- Adds unnecessary complexity
- Risk of version drift

## Action Items

- [x] Update `.gitignore` to include `docs/source/` but exclude `docs/build/`
- [x] Stage documentation files for commit
- [ ] Commit documentation to main branch
- [ ] Push to GitHub
- [ ] Connect to Read the Docs
- [ ] Verify RTD builds successfully

## Commands to Complete Setup

```bash
# Check what will be committed (should be ~69 files)
git status docs/ --short

# Commit documentation
git commit -m "Add Sphinx documentation with Read the Docs theme"

# Push to GitHub
git push origin master

# Then: Go to readthedocs.org and import repository
```

---

**Recommendation: Proceed with single-repo approach.** ✅
