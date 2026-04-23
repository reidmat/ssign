# `scripts/` — one-off analysis scripts

Scripts here perform single-purpose analyses that produced results or figures
in the ssign paper. They are **not** part of the ssign pipeline itself — those
live in `src/ssign_app/` and are importable from the installed package.

## When to put a script here vs in `src/ssign_app/`

- **Here (`scripts/`):** one-off paper-specific analysis, figure-regeneration
  wrappers, dataset-curation helpers. Not imported by other code. Runnable
  standalone.
- **In `src/ssign_app/`:** reusable modules, CLI entry points, code that any
  ssign user would invoke as part of the pipeline.

If a script becomes generally useful to multiple analyses, promote it into
`src/ssign_app/` as a module.

## Expected scripts (to be populated in Phase 8)

- `fetch_weights.sh` — download model weights from Zenodo.
- `fetch_databases.sh` — download reference databases from Zenodo / RDS.
- Per-figure regeneration scripts referenced from `figures/`.
- Validation-set curation scripts.
- Benchmark driver scripts.

## Running

Every script should be runnable from the repo root:

```bash
cd /path/to/ssign
python scripts/<script>.py [args]
# or
bash scripts/<script>.sh
```

and should document its inputs, outputs, and any required environment at the
top of the file.

## Style

- Python scripts: use the same deps as the main package; no separate venv.
- Shell scripts: POSIX-compatible, `set -euo pipefail` at the top.
- Document non-obvious assumptions in comments.
- Record random seeds for any stochastic analysis.

## Status

Directory scaffold only; content populated as paper analyses are finalised.
