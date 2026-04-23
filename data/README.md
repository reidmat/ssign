# `data/` — reference databases and raw inputs

This directory is a **placeholder**. The reference databases used by ssign are
too large (~1.5 TB combined) to distribute through Git, so they are hosted
externally and fetched at install time.

## What goes here at runtime

After a full install, this directory (or a path you pass to ssign via config)
will contain:

| Database                                 | Tier            | Size    | Purpose                                               |
| ---------------------------------------- | --------------- | ------- | ----------------------------------------------------- |
| NCBI taxdump (`nodes.dmp` + `names.dmp`) | base            | ~1.5 GB | Local taxonomy resolution (replaces NCBI E-utilities) |
| Bakta (light)                            | base            | ~2 GB   | Gene annotation (light variant)                       |
| Bakta (full)                             | full            | ~84 GB  | Whole-genome gene annotation                          |
| EggNOG-mapper                            | extended / full | ~50 GB  | Orthologous groups, COGs, KEGG                        |
| HH-suite (Pfam + PDB70)                  | extended        | ~20 GB  | Pfam + PDB70 profile HMM databases                    |
| HH-suite (+ UniRef30)                    | full            | +~25 GB | Adds UniRef30 for deep remote-homology MSA generation |
| InterProScan                             | extended / full | ~24 GB  | Domain / motif annotation                             |
| pLM-BLAST (ECOD70)                       | extended / full | ~20 GB  | Embedding-based remote homology                       |
| BLAST NR                                 | full            | ~390 GB | Primary homology search                               |

## Getting the databases

At v1.0.0 release, ssign will ship a fetch script:

```bash
bash scripts/fetch_databases.sh             # full ~1.5 TB bundle
bash scripts/fetch_databases.sh --minimal   # ~20 GB subset for testing
```

The databases are (will be) mirrored to **Zenodo** (DOI TBD at release). If a
deposit exceeds Zenodo's per-record size limit, it will be split across
multiple linked records, or hosted on Imperial's Research Data Store with a
persistent URL and checksum manifest.

## Provenance — original data used in the paper

Analysis inputs used to generate figures and validation sets for the ssign
paper (Reid et al., _in preparation_) are listed here. Raw data is not
redistributed via this repo; each dataset has a persistent source.

- _List to be populated in Phase 3–8 as analysis is finalised._

## Integrity checking

Every file delivered via `fetch_databases.sh` has an SHA-256 checksum recorded
alongside it. To verify after download:

```bash
cd data && sha256sum -c checksums.sha256
```

## What must **not** go here

- User-uploaded genomes or proteins for a specific run (those live in the
  user's chosen output directory).
- Any cleaned, preprocessed, or derived data (that belongs in `results/`).
- Model weights (those live in `models/`).
