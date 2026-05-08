#!/usr/bin/env bash
# tests/run_tutorial.sh — verify the E. coli K-12 tutorial in docs/tutorials/first_run.md
#
# Runs the exact commands from the tutorial against a real K-12 genome,
# checks expected outputs exist + are non-empty, and spot-checks that the
# T2SS gsp operon shows up in the secretion-systems section. Use this as
# a smoke test before tagging a release, or after touching anything in
# the run-pipeline path.
#
# Usage:
#     bash tests/run_tutorial.sh            # ~30 min, internet required
#     bash tests/run_tutorial.sh --force    # re-run even if output exists
#
# Requires: ssign on PATH (pip install ssign), wget, gunzip, python3.
# Internet required for: K-12 GenBank download + DTU webserver calls.

set -euo pipefail

WORK_DIR="${HOME}/ssign-tutorial-test"
OUTDIR="${WORK_DIR}/ecoli_results"
GENBANK="${WORK_DIR}/ecoli_k12.gbff"
GENBANK_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gbff.gz"

FORCE=0
[[ "${1:-}" == "--force" ]] && FORCE=1

# ── Sanity checks ────────────────────────────────────────────────────────
command -v ssign >/dev/null || { echo "FAIL: ssign not on PATH. Install with 'pip install ssign'." >&2; exit 2; }
command -v wget  >/dev/null || { echo "FAIL: wget not installed." >&2; exit 2; }
command -v python3 >/dev/null || { echo "FAIL: python3 not installed." >&2; exit 2; }

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# ── Download K-12 (~13 MB) ──────────────────────────────────────────────
if [[ ! -f "$GENBANK" ]]; then
    echo ">> Downloading E. coli K-12 GenBank from NCBI..."
    wget -q "$GENBANK_URL" -O ecoli_k12.gbff.gz
    gunzip -f ecoli_k12.gbff.gz
    [[ -f "$GENBANK" ]] || { echo "FAIL: download did not produce $GENBANK" >&2; exit 1; }
fi
echo ">> K-12 GenBank: $(du -h "$GENBANK" | cut -f1)"

# ── Clean output ─────────────────────────────────────────────────────────
if [[ -d "$OUTDIR" ]]; then
    if (( FORCE )); then
        echo ">> --force: clearing $OUTDIR"
        rm -rf "$OUTDIR"
    else
        echo "FAIL: $OUTDIR exists. Pass --force to re-run, or delete manually." >&2
        exit 1
    fi
fi

# ── Run ssign (the exact tutorial command) ───────────────────────────────
echo ">> Running ssign — DLP + SignalP via DTU webserver, ~30 min total..."
START=$(date +%s)
ssign run "$GENBANK" --outdir "$OUTDIR" \
    --use-input-annotations \
    --signalp-mode remote \
    --deeplocpro-mode remote \
    --skip-blastp \
    --skip-eggnog \
    --skip-hhsuite \
    --skip-interproscan \
    --skip-plmblast
ELAPSED=$(( $(date +%s) - START ))
echo ">> ssign finished in ${ELAPSED}s"

# ── Verify outputs ───────────────────────────────────────────────────────
echo ">> Checking output files..."
EXPECTED=(
    "$OUTDIR/ecoli_k12_results.csv"
    "$OUTDIR/ecoli_k12_results_raw.csv"
    "$OUTDIR/ecoli_k12_summary.txt"
    "$OUTDIR/figures/ecoli_k12"
)
for f in "${EXPECTED[@]}"; do
    if [[ ! -e "$f" ]]; then
        echo "FAIL: missing $f" >&2
        exit 1
    fi
    if [[ -f "$f" && ! -s "$f" ]]; then
        echo "FAIL: $f is empty" >&2
        exit 1
    fi
done
echo "   all expected outputs present + non-empty"

# ── Spot-check: T2SS should appear in section 2 of the chunked CSV ──────
echo ">> Spot-checking K-12 biology..."
python3 - "$OUTDIR/ecoli_k12_results.csv" <<'PY'
import sys
import pandas as pd

path = sys.argv[1]
with open(path) as fh:
    raw = fh.read()

# Section 2 is between "# Section 2:" and "# Section 3:"; chunked CSV
# uses these markers as separators (output_files.md § Layout).
if "# Section 2" not in raw:
    sys.exit("FAIL: no Section 2 header in chunked CSV")

# Read the whole CSV and find rows whose ss_type column contains T2SS.
df = pd.read_csv(path, comment="#", skip_blank_lines=True)
if "ss_type" not in df.columns:
    sys.exit(f"FAIL: ss_type column missing. Columns: {list(df.columns)}")
t2ss_rows = df[df["ss_type"].astype(str).str.contains("T2SS", na=False)]
if t2ss_rows.empty:
    sys.exit("FAIL: no T2SS rows in K-12 results — gsp operon not detected?")

print(f"   T2SS rows found: {len(t2ss_rows)}")
PY

echo
echo "=========================================="
echo "TUTORIAL VERIFICATION PASSED in ${ELAPSED}s"
echo "=========================================="
echo "Outputs: $OUTDIR"
