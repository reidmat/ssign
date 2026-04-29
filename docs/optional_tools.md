# Optional Tool Installation Guide

ssign works out of the box with `pip install ssign`. The tools below extend its capabilities and can be installed at any time.

---

## DeepSecE — Additional Secreted Protein Prediction

DeepSecE predicts which secretion system type a protein is secreted by, providing cross-validation against DeepLocPro. Download is approximately **7.3 GB** (PyTorch + ESM protein language model).

**Important:** Install inside your ssign virtual environment.

### All platforms

```bash
# Activate your ssign environment first
source ~/.ssign-env/bin/activate    # Linux/macOS
# or on Windows WSL: source ~/.ssign-env/bin/activate

pip install ssign[deepsece]
```

ssign will automatically detect DeepSecE on next launch.

---

## Bakta — Higher-Quality Gene Prediction

If your input is raw FASTA contigs (not GenBank), ssign uses Pyrodigal by default. Bakta provides richer functional annotation. Requires a **~2 GB** database download.

### All platforms

```bash
source ~/.ssign-env/bin/activate

pip install ssign[bakta]
bakta_db download --output ~/bakta_db --type light
```

Then enable Bakta in the ssign GUI and enter the database path (`~/bakta_db`).

---

## BLAST+ — Ortholog Grouping Across Genomes

When analysing multiple genomes, BLAST+ enables ortholog grouping to identify shared secreted proteins. BLAST+ is **not pip-installable** and must be installed through your system package manager. Download is approximately **200 MB**.

### Linux (Ubuntu / Debian)

```bash
sudo apt install ncbi-blast+
```

### macOS

```bash
brew install blast
```

### Conda (cross-platform)

```bash
conda install -c bioconda blast
```

ssign will automatically detect BLAST+ on next launch.

---

## DeepLocPro (Local Mode) — Faster Localization Prediction

By default, ssign runs DeepLocPro via the free DTU cloud API. Local mode is faster but requires a DTU academic license and approximately **5 GB** of model files. GPU recommended.

### All platforms

1. Register for a free academic license at [DTU DeepLocPro](https://services.healthtech.dtu.dk/services/DeepLocPro-1.0/)
2. Download and install following DTU's instructions
3. In ssign, select "Local install" for DeepLocPro and enter the install path

---

## SignalP 6.0 (Local Mode) — Faster Signal Peptide Prediction

By default, ssign runs SignalP via the free DTU cloud API. Local mode is faster but requires a DTU academic license and approximately **1.5 GB** download.

### Why a separate environment

SignalP 6.0 pins **Python ≤ 3.10** and **PyTorch < 2.0**, while ssign itself runs on Python 3.11+ with PyTorch 2.x. Installing SignalP into your ssign env will downgrade PyTorch and break DeepSecE / DeepLocPro / PLM-Effector / pLM-BLAST. Install SignalP into its own env and point ssign at the binary via `--signalp-path`.

### Install (Linux/macOS)

```bash
# 1. Register and download from DTU (requires academic email):
#    https://services.healthtech.dtu.dk/services/SignalP-6.0/
#    Pick the "fast" model variant (signalp-6.0i.fast.tar.gz, ~1.5 GB).

# 2. Create a dedicated Python 3.10 env (using micromamba, conda, or venv).
#    Example with micromamba:
micromamba create -n signalp6 -c conda-forge python=3.10 pip "numpy<2" -y
micromamba activate signalp6

# 3. Pre-install CPU-only torch<2 BEFORE the package, so torch 2.x isn't
#    pulled in transitively (also keeps the wheel small at ~150 MB):
pip install "torch<2.0" --index-url https://download.pytorch.org/whl/cpu

# 4. Extract and install:
tar xzf ~/Downloads/signalp-6.0i.fast.tar.gz -C ~/build
cd ~/build/signalp6_fast
pip install ./signalp-6-package/

# 5. Copy the model weights into the installed package directory:
SIGNALP_DIR=$(python -c "import signalp, os; print(os.path.dirname(signalp.__file__))")
cp -r signalp-6-package/models/* "$SIGNALP_DIR/model_weights/"

# 6. Verify:
signalp6 --version
```

### Wire ssign to the local install

Pass the directory containing the `signalp6` console script (typically `<env>/bin`) via `--signalp-path` on the CLI, or set it in the GUI's local-install field.

```bash
# Example for the integration test:
SSIGN_SIGNALP_PATH=~/micromamba/envs/signalp6/bin pytest -m integration \
    tests/integration/test_run_signalp_integration.py::TestLocal
```

ssign invokes SignalP with `--organism other --mode fast --format txt` (gram-negative bacteria use the `other` group in v6 — v5's `gram-` was removed).

---

## pLM-BLAST — Remote-Homology Search via ProtT5 Embeddings

Not on PyPI; clone from GitHub and point ssign at the script.

```bash
git clone https://github.com/labstructbioinf/pLM-BLAST.git ~/pLM-BLAST
export SSIGN_PLMBLAST_SCRIPT=~/pLM-BLAST/scripts/plmblast.py
```

Pre-built ECOD70 database (~21 GB compressed, ~24 GB extracted):

```bash
mkdir -p ~/pLM-BLAST/db && cd ~/pLM-BLAST/db
wget http://ftp.tuebingen.mpg.de/ebio/protevo/toolkit/databases/plmblast_dbs/ecod70db_20240417.tar.gz
tar -xzf ecod70db_20240417.tar.gz && rm ecod70db_20240417.tar.gz
export SSIGN_ECOD70_DB=~/pLM-BLAST/db/ECOD70
```

GPU strongly recommended: ProtT5 embedding takes ~5-10 sec per
500-aa protein on CPU vs ~0.1 sec on a modern GPU. A whole-genome
run (~5,000 proteins) on CPU is ~10+ hours just for embedding.

---

## PLM-Effector — Five-Type Secretion Effector Prediction

Vendored under CC-BY 3.0 in `src/ssign_app/scripts/plm_effector/`. Code
ships with ssign; only the weights need to be downloaded separately.

```bash
# Download trained models (~1.7 GB)
wget https://www.mgc.ac.cn/PLM-Effector/download/sourcecode.zip
unzip sourcecode.zip
mv sourcecode/trained_models /path/to/plm_effector_weights/

# Download pretrained PLMs from HuggingFace (~17 GB total)
hf download Rostlab/prot_t5_xl_uniref50 \
    --include "config.json" "spiece.model" "tokenizer_config.json" \
    "special_tokens_map.json" "pytorch_model.bin" \
    --local-dir /path/to/plm_effector_weights/transformers_pretrained/prot_t5_xl_uniref50
hf download facebook/esm1b_t33_650M_UR50S \
    --local-dir /path/to/plm_effector_weights/transformers_pretrained/esm1b_t33_650M_UR50S
hf download facebook/esm2_t33_650M_UR50D \
    --local-dir /path/to/plm_effector_weights/transformers_pretrained/esm2_t33_650M_UR50D
hf download Rostlab/prot_bert \
    --local-dir /path/to/plm_effector_weights/transformers_pretrained/prot_bert

export SSIGN_PLM_EFFECTOR_WEIGHTS=/path/to/plm_effector_weights
```

CUDA GPU required for practical runtime: ~100x slower on CPU. ssign's
PLM-Effector test hard-skips on no-GPU systems.

---

## Environment constraints (extended tier)

A few Python deps are version-pinned because of upstream API changes
that break the bioinformatics tools:

| Package | Pin | Reason |
|---|---|---|
| `transformers` | `>=4.38,<5.0` | 5.0 removed `batch_encode_plus` (used by pLM-BLAST + PLM-Effector tokenizers) |
| `numpy` | `>=1.26,<2.0` | 2.0 removed `np.issubsctype` (used by pLM-BLAST's alignment code) |
| `protobuf` | any | Required by ProtT5's SentencePiece tokenizer at load time |
| `mkl`, `mkl-service` | any | pLM-BLAST's `plmblast.py` imports these directly |

These are captured in the `extended` and `full` extras in
`pyproject.toml`, so `pip install ssign[extended]` resolves to the
correct versions automatically. Listed here as a reference if you're
debugging an install or assembling a Conda env from scratch.

---

## HH-suite — Profile-vs-Profile Remote Homology Search

`hhsearch` and `hhblits` (Steinegger / Söding labs) detect remote
evolutionary relationships by comparing HMM-vs-HMM profiles. ssign
uses HH-suite to annotate substrate proteins with Pfam domain
families and PDB structural homologs.

### Install the binaries

```bash
# Bioconda (recommended)
conda install -c bioconda hhsuite

# Or build from source: https://github.com/soedinglab/hh-suite
```

### Download the databases

ssign needs three precomputed databases. Two canonical mirrors host
them; **prefer Tübingen** (fresher, recommended by Söding lab issue
#382) and fall back to GWDG only if Tübingen is unreachable.

```bash
# Pfam — domain families. Tübingen has v38 (2024+), GWDG has v35 (2021).
mkdir -p $HHSUITE_DBS && cd $HHSUITE_DBS
wget http://ftp.tuebingen.mpg.de/pub/ebio/protevo/toolkit/databases/hhsuite_dbs/PfamA_v38_2.tar.gz
tar -xzf PfamA_v38_2.tar.gz && rm PfamA_v38_2.tar.gz

# PDB70 — structural homology. Tübingen has 2026-02 build, GWDG has 2022-03.
wget http://ftp.tuebingen.mpg.de/pub/ebio/protevo/toolkit/databases/hhsuite_dbs/pdb70_from_mmcif_2026-02-20.tar.gz
tar -xzf pdb70_from_mmcif_2026-02-20.tar.gz && rm pdb70_from_mmcif_2026-02-20.tar.gz

# UniRef30 — clustered UniProt for hhblits MSA generation. Only at GWDG.
wget https://wwwuser.gwdg.de/~compbiol/uniclust/2023_02/UniRef30_2023_02_hhsuite.tar.gz
tar -xzf UniRef30_2023_02_hhsuite.tar.gz && rm UniRef30_2023_02_hhsuite.tar.gz
```

Total disk after extraction: ~55 GB. Tell ssign where to find them:

```bash
export SSIGN_HHSUITE_PFAM=$HHSUITE_DBS/pfam
export SSIGN_HHSUITE_PDB70=$HHSUITE_DBS/pdb70_from_mmcif_2026-02-20
export SSIGN_HHSUITE_UNICLUST=$HHSUITE_DBS/UniRef30_2023_02
```

### Mirror caveats

The Söding lab is in maintenance mode: per [hh-suite issue
#382](https://github.com/soedinglab/hh-suite/issues/382), no new
funding for support. Files at GWDG are stable but stale. Tübingen
mirror has fresher builds. If both go down for an extended period,
ssign's Phase 7a HPC test (where DBs are pre-staged on CX3) is the
fallback path.

`curl -I` (without `-L`) reports these mirrors as 404 because the
GWDG host issues a 302 redirect to a sibling domain. **Always probe
with `curl -IL`** (capital L = follow redirects) to verify.
