"""Top-level PLM-Effector prediction API.

Wraps `feature_extraction.extract_all_features` and `ensemble.run_ensemble`
into a single `predict()` call. Writes results as a TSV matching the
upstream output format plus an explicit `passes_threshold` column so
downstream ssign code does not need to know each effector type's
threshold.
"""

from __future__ import annotations

import logging
import os

import numpy as np

# torch / ensemble / feature_extraction are imported lazily inside the
# heavy functions so `write_predictions_tsv` (pure-numpy, unit-testable)
# can be imported without the ML deps installed.


logger = logging.getLogger(__name__)

_VALID_EFFECTOR_TYPES = ("T1SE", "T2SE", "T3SE", "T4SE", "T6SE")


def _resolve_device(device):
    """Accept None / "cuda" / "cpu" / torch.device and return a torch.device."""
    import torch

    if isinstance(device, torch.device):
        return device
    if device is None or device == "cuda":
        if not torch.cuda.is_available():
            raise RuntimeError(
                "PLM-Effector requires a CUDA GPU (default) but torch reports "
                "no CUDA device available. Set device='cpu' to force CPU "
                "inference (very slow, ~1-2 min per protein)."
            )
        return torch.device("cuda")
    if device == "cpu":
        return torch.device("cpu")
    raise ValueError(f"Unrecognised device: {device!r}")


def predict(
    proteins_fasta: str,
    weights_dir: str,
    effector_type: str,
    out_path: str,
    device=None,
    batch_size: int = 5,
) -> int:
    """Run PLM-Effector on a protein FASTA and write predictions to TSV.

    Args:
        proteins_fasta: input protein FASTA.
        weights_dir: directory containing
            `transformers_pretrained/<model>/` subdirs (HuggingFace PLMs)
            and `trained_models/<effector_type>_*.pth` + the
            `<effector_type>_XGB_stackingmeta_model.json` stacking model.
        effector_type: one of `T1SE`, `T2SE`, `T3SE`, `T4SE`, `T6SE`.
        out_path: output TSV path.
        device: None (auto-pick cuda), "cuda", "cpu", or a torch.device.
            CPU works but is very slow.
        batch_size: PLM forward-pass batch size. Drop to 1-2 if VRAM-constrained.

    Returns:
        The number of proteins flagged as positive (passing threshold).
    """
    if effector_type not in _VALID_EFFECTOR_TYPES:
        raise ValueError(
            f"effector_type must be one of {_VALID_EFFECTOR_TYPES}; got {effector_type!r}"
        )
    if not os.path.exists(proteins_fasta):
        raise FileNotFoundError(f"Input FASTA not found: {proteins_fasta}")
    if not os.path.isdir(weights_dir):
        raise FileNotFoundError(
            f"PLM-Effector weights directory not found: {weights_dir}\n"
            f"  Fetch with scripts/fetch_weights.sh or mirror manually from "
            f"http://www.mgc.ac.cn/PLM-Effector/downloads.html"
        )

    from .ensemble import run_ensemble
    from .feature_extraction import extract_all_features

    torch_device = _resolve_device(device)
    logger.info(
        "PLM-Effector: extracting features for %s on %s (batch_size=%d)",
        effector_type,
        torch_device,
        batch_size,
    )
    features = extract_all_features(
        proteins_fasta=proteins_fasta,
        effector_type=effector_type,
        weights_dir=weights_dir,
        device=torch_device,
        batch_size=batch_size,
    )

    logger.info("PLM-Effector: running ensemble for %s", effector_type)
    seq_ids, stacked, final_probs, passes = run_ensemble(
        features=features,
        weights_dir=weights_dir,
        effector_type=effector_type,
        device=torch_device,
    )

    n_base = stacked.shape[1]
    write_predictions_tsv(
        out_path, seq_ids, stacked, final_probs, passes, effector_type, n_base
    )

    n_positive = int(passes.sum())
    logger.info(
        "PLM-Effector: wrote %d predictions (%d passing threshold) to %s",
        len(seq_ids),
        n_positive,
        out_path,
    )
    return n_positive


def write_predictions_tsv(
    out_path: str,
    seq_ids: np.ndarray,
    stacked: np.ndarray,
    final_probs: np.ndarray,
    passes: np.ndarray,
    effector_type: str,
    n_base_models: int,
) -> None:
    """Write predictions as TSV: `seq_id\\tmodel1..N\\tstacking\\tpasses_threshold\\teffector_type`.

    One row per input protein. Upstream's format only emitted rows above
    threshold; we emit all proteins and let downstream code filter so the
    output shape is predictable regardless of positive-hit count.
    """
    headers = (
        ["seq_id"]
        + [f"model{i + 1}" for i in range(n_base_models)]
        + ["stacking", "passes_threshold", "effector_type"]
    )
    with open(out_path, "w") as f:
        f.write("\t".join(headers) + "\n")
        for i, sid in enumerate(seq_ids):
            sid_str = sid.decode() if isinstance(sid, bytes) else str(sid)
            # Strip leading '>' if present (upstream parser kept it)
            if sid_str.startswith(">"):
                sid_str = sid_str[1:].split()[0]
            row = (
                [sid_str]
                + [f"{stacked[i, j]:.4f}" for j in range(n_base_models)]
                + [
                    f"{final_probs[i]:.4f}",
                    "1" if passes[i] else "0",
                    effector_type,
                ]
            )
            f.write("\t".join(row) + "\n")
