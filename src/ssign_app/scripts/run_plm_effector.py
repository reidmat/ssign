#!/usr/bin/env python3
"""Run PLM-Effector secreted-protein prediction and convert output to ssign format.

PLM-Effector (Zheng et al. 2026, Briefings in Bioinformatics 27(2) bbag143)
predicts bacterial secretion effectors using a two-layer ensemble over
pretrained protein language models (ESM-1b, ESM-2, ProtT5) plus per-type
trained neural networks and an XGBoost stacking meta-model. Vendored from
https://github.com/zhengdd0422/PLM-Effector under CC-BY 3.0 into
`src/ssign_app/scripts/plm_effector/` (see that package's `LICENSE` file
and `__init__.py` for attribution details).

Usage:
    run_plm_effector.py --input proteins.faa \\
        --weights-dir /path/to/plm_effector_weights \\
        --effector-type T1SE \\
        --out predictions.tsv \\
        [--device cuda|cpu] [--batch-size 5]

Weights directory layout expected:
    {weights-dir}/transformers_pretrained/
        esm1b_t33_650M_UR50S/
        esm2_t33_650M_UR50D/
        prot_t5_xl_uniref50/
        prot_bert/                 (only needed for T4SE)
    {weights-dir}/trained_models/
        {TYPE}_model{i}_fold{f}.pth
        {TYPE}_XGB_stackingmeta_model.json

Output format (TSV):
    seq_id  model1  model2  ...  stacking  passes_threshold  effector_type
"""

from __future__ import annotations

import argparse
import logging
import os
import sys

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


_VALID_EFFECTOR_TYPES = ("T1SE", "T2SE", "T3SE", "T4SE", "T6SE")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run PLM-Effector and convert output to ssign format",
    )
    parser.add_argument("--input", required=True, help="Input protein FASTA")
    parser.add_argument(
        "--weights-dir",
        required=True,
        help="Directory containing PLM-Effector weights (see docstring for layout)",
    )
    parser.add_argument(
        "--effector-type",
        required=True,
        choices=_VALID_EFFECTOR_TYPES,
        help="Secretion system effector type to predict",
    )
    parser.add_argument(
        "--out", required=True, help="Output predictions TSV (ssign format)"
    )
    parser.add_argument(
        "--device",
        default="cuda",
        choices=["cuda", "cpu"],
        help="Inference device (default: cuda; cpu works but is ~100x slower)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=5,
        help="PLM forward-pass batch size (drop to 1-2 if VRAM-constrained)",
    )
    args = parser.parse_args()

    # FRAGILE: PLM-Effector requires torch + transformers + xgboost +
    # NumPy/BioPython, plus the pretrained PLM weights (~15 GB on disk).
    # If this breaks:
    #   - pip install ssign[extended]    # installs torch, transformers, xgboost
    #   - scripts/fetch_weights.sh       # fetches pretrained PLM + trained models
    try:
        from ssign_app.scripts.plm_effector import predict
    except ImportError as e:
        print(
            f"ERROR: PLM-Effector dependencies not available: {e}\n"
            f"  Common causes:\n"
            f"    - transformers / xgboost / torch are not installed\n"
            f"  How to fix:\n"
            f"    - pip install ssign[extended]\n"
            f"    - Or install manually: pip install transformers xgboost torch",
            file=sys.stderr,
        )
        return 2

    if not os.path.exists(args.input):
        print(f"ERROR: Input FASTA not found: {args.input}", file=sys.stderr)
        return 2
    if not os.path.isdir(args.weights_dir):
        print(
            f"ERROR: Weights directory not found: {args.weights_dir}\n"
            f"  Fetch with scripts/fetch_weights.sh or mirror from "
            f"http://www.mgc.ac.cn/PLM-Effector/downloads.html",
            file=sys.stderr,
        )
        return 2

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)

    try:
        n_positive = predict(
            proteins_fasta=args.input,
            weights_dir=args.weights_dir,
            effector_type=args.effector_type,
            out_path=args.out,
            device=args.device,
            batch_size=args.batch_size,
        )
    except RuntimeError as e:
        print(f"ERROR: PLM-Effector prediction failed: {e}", file=sys.stderr)
        return 1
    except FileNotFoundError as e:
        print(
            f"ERROR: required PLM-Effector weight file missing:\n  {e}\n"
            f"  See --weights-dir layout in the run_plm_effector.py docstring.",
            file=sys.stderr,
        )
        return 2

    logger.info(
        "Done: %d proteins flagged as %s substrates (full table at %s)",
        n_positive,
        args.effector_type,
        args.out,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
