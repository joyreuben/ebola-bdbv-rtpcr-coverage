# GitHub workflow

How this project got onto GitHub, and how to keep tracking changes from here.

## One-time setup (done 2026-08-25)

1. Installed GitHub CLI into the conda base environment (not tied to any one
   project env, since it's a general tool):
   ```
   conda install -n base -c conda-forge gh -y
   ```
2. Authenticated via device-flow login (opens a one-time code + URL, approved
   in the browser — no password or SSH key needed):
   ```
   gh auth login --hostname github.com --git-protocol https --web
   ```
   Logged in as `joyreuben`. Check status any time with `gh auth status`.
3. Added `.gitignore` (excludes `__pycache__/`, editor/OS cruft, and the local
   `.claude/` session-state folder — none of that belongs in version control).
4. Initialized the repo and made the first commit:
   ```
   git init
   git branch -m main
   git add -A
   git commit -m "..."
   ```
5. Created the GitHub repo and pushed in one step:
   ```
   gh repo create ebola-bdbv-rtpcr-coverage --public --source=. --remote=origin --push
   ```
   Repo: https://github.com/joyreuben/ebola-bdbv-rtpcr-coverage

Note: `gh` was installed via conda, so it's on PATH inside a shell that has
run `conda init` / activated an environment. If a fresh terminal reports
`gh: command not found`, either open a shell where conda's `base` env is on
PATH, or run it via the full path:
`/home/joy/miniconda3/bin/gh`.

## Day-to-day: committing your progress

Run from `~/ebola project/`:

```bash
git status                 # see what changed
git diff                   # review the actual changes
git add <files>             # stage specific files (avoid `git add -A` if unsure —
                             # check status first so you don't commit something unintended)
git commit -m "short description of what and why"
git push
```

Suggested cadence: commit at the end of each step (e.g. after fixing Step 5,
after running Step 6 thermodynamics, after generating each figure) rather
than one giant commit at the end. Small commits make it easy to see, in
`git log`, exactly when a result changed and why — useful when writing the
methods section and when responding to reviewer questions later.

## Suggested commit messages for the known next steps

Done: Step 5 phantom-10 fix, Step 6 thermodynamics (see `git log` for the
actual messages used). Remaining:

- `Add coverage heatmap and binding-site schematic figures`
- `Add phylogeny annotated with mismatch lineages`
- `Draft results write-up`

## Checking history

```bash
git log --oneline          # compact history
git log -p -- step5_combine_assays.py   # full history of one file's changes
git show <commit>           # see exactly what a given commit changed
```

## Before every commit: check for secrets

This repo is **public**. Before staging anything new, glance at `git status`
and, if a file looks unfamiliar, open it — don't commit API keys, tokens, or
anything not meant to be public. The source genome data itself is fine (it's
OPEN-use-terms Pathoplexus data, already public).

## Resolved documentation gap (was: known gap)

`CLAUDE.md` used to reference a `diagnose_phantom10.py` script that was
never present in the working directory. The Step 5 skip-fraction fix was
applied directly instead, verified working (Marburg control now 0/77
LIKELY DETECTED), and `CLAUDE.md` was updated to drop the reference. See
`git log -- CLAUDE.md` for the update commit.
