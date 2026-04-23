# `models/` — trained model weights used by ssign

This directory is a **placeholder**. Model weights are too large to version in
Git, so they are hosted externally and fetched at install time.

## What goes here at runtime

| Model                | Size     | Used by                                   |
| -------------------- | -------- | ----------------------------------------- |
| DeepSecE checkpoint  | ~500 MB  | `run_deepsece.py`                         |
| ESM-1b               | ~7 GB    | DeepSecE feature extraction               |
| ESM-2                | ~3 GB    | PLM-Effector feature extraction           |
| ProtT5               | ~2.5 GB  | PLM-Effector feature extraction           |
| PLM-Effector weights | ~200 MB  | `run_plm_effector.py`                     |
| DeepLocPro weights   | ~2 GB    | `run_deeplocpro.py`                       |
| SignalP 6.0 weights  | variable | `run_signalp.py` (subject to DTU license) |

## Getting the weights

At v1.0.0 release, a fetch script will pull everything from our Zenodo
deposit:

```bash
bash scripts/fetch_weights.sh
```

The **DeepSecE checkpoint** is currently hosted on an unreliable SJTU server;
mirroring to Zenodo before release is part of the longevity mitigation stack
(see project plan).

## DTU academic-licensed models

`SignalP 6.0` and `DeepLocPro` are under DTU academic licenses that may
restrict redistribution. If the Docker bundle cannot ship them directly,
install docs will cover DTU registration + separate install.

## Integrity checking

```bash
cd models && sha256sum -c checksums.sha256
```

## What must **not** go here

- User-trained or user-fine-tuned models (those are out-of-scope for ssign).
- Raw training data (that's in `data/` if anywhere).
- Model outputs (those are pipeline `results/`).
