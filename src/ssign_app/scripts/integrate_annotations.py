#!/usr/bin/env python3
"""Integrate all annotation sources into master CSV.

Left-join merges on locus_tag, adds GBFF annotation from gene_info,
protein sequences, and computes annotation tool hit counts.

Annotation tools counted: BLASTp, HHpred (Pfam), HHpred (PDB),
InterProScan, ProtParam, GBFF (original genome annotation).

SignalP is a secretion prediction tool (not an annotation tool) and
is already included in the substrates table via cross_validate.
"""

import argparse
import logging
import os

import pandas as pd
from Bio import SeqIO

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# Annotation tool columns — used for n_tools_with_hits counting
# SignalP is NOT here (it's a prediction tool, not annotation)
TOOL_HIT_COLUMNS = {
    'BLASTp': 'blastp_hit_description',
    'HHpred_Pfam': 'pfam_top1_description',
    'HHpred_PDB': 'pdb_top1_description',
    'InterProScan': 'interpro_domains',
    'Foldseek': 'foldseek_hit_description',
    'pLM-BLAST': 'ecod70_top1_description',
    'GBFF': 'gbff_annotation',
}


def _compute_tool_counts(df):
    """Add n_tools_with_hits and annotation_tools columns."""
    def _count_hits(row):
        hits = []
        for tool, col in TOOL_HIT_COLUMNS.items():
            if col not in df.columns:
                continue
            val = row.get(col, '')
            if pd.isna(val) or not str(val).strip():
                continue
            # For GBFF, skip generic annotations
            if tool == 'GBFF' and str(val).strip().lower() in (
                'hypothetical protein', 'uncharacterized protein', ''):
                continue
            hits.append(tool)
        return hits

    hit_lists = df.apply(_count_hits, axis=1)
    df['n_tools_with_hits'] = hit_lists.apply(len)
    df['annotation_tools'] = hit_lists.apply(lambda x: ','.join(x) if x else '')
    return df


def main():
    parser = argparse.ArgumentParser(description="Integrate annotations")
    parser.add_argument("--substrates-filtered", required=True)
    parser.add_argument("--substrates-all", required=True)
    parser.add_argument("--annotations", nargs='*', default=[])
    parser.add_argument("--gene-info", default="",
                        help="Gene info TSV from extract_proteins (for GBFF annotation)")
    parser.add_argument("--proteins", default="",
                        help="Proteins FASTA (to include sequences in output)")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # Load base substrate table
    df = pd.read_csv(args.substrates_filtered, sep='\t')
    logger.info(f"Base: {len(df)} filtered substrates for {args.sample}")

    # Add GBFF annotation from gene_info (original genome annotations)
    if args.gene_info and os.path.exists(args.gene_info):
        try:
            gi = pd.read_csv(args.gene_info, sep='\t')
            if 'locus_tag' in gi.columns and 'product' in gi.columns:
                gi_ann = gi[['locus_tag', 'product']].copy()
                gi_ann = gi_ann.rename(columns={'product': 'gbff_annotation'})
                gi_ann = gi_ann.drop_duplicates(subset='locus_tag')
                before = len(df)
                df = df.merge(gi_ann, on='locus_tag', how='left')
                assert len(df) == before
                logger.info(f"Added GBFF annotations from gene_info")
        except Exception as e:
            logger.warning(f"Failed to add GBFF annotations: {e}")

    # Add protein sequences
    if args.proteins and os.path.exists(args.proteins):
        try:
            seqs = {}
            for rec in SeqIO.parse(args.proteins, 'fasta'):
                seqs[rec.id] = str(rec.seq)
            if seqs and 'locus_tag' in df.columns:
                df['sequence'] = df['locus_tag'].map(seqs)
                df['aa_length'] = df['sequence'].apply(
                    lambda s: len(s) if pd.notna(s) else 0)
                logger.info(f"Added sequences for {df['sequence'].notna().sum()} proteins")
        except Exception as e:
            logger.warning(f"Failed to add sequences: {e}")

    # Left-join each annotation file
    for ann_file in args.annotations:
        if not os.path.exists(ann_file):
            continue

        try:
            ann_df = pd.read_csv(ann_file, sep=None, engine='python')

            # Determine join column
            join_col = 'locus_tag' if 'locus_tag' in ann_df.columns else None
            if not join_col:
                for col in ann_df.columns:
                    if 'protein' in col.lower() or 'id' == col.lower():
                        join_col = col
                        break

            if join_col and join_col in df.columns:
                # Avoid duplicate columns
                overlap = set(ann_df.columns) & set(df.columns) - {join_col}
                if overlap:
                    ann_df = ann_df.drop(columns=list(overlap))

                before = len(df)
                df = df.merge(ann_df, on=join_col, how='left')
                assert len(df) == before, (
                    f"Row count changed after merging {ann_file}: "
                    f"{before} -> {len(df)}"
                )
                logger.info(f"Merged {ann_file}: +{len(ann_df.columns)-1} columns")
            else:
                logger.warning(f"No join column found in {ann_file}, skipping")

        except Exception as e:
            logger.warning(f"Failed to merge {ann_file}: {e}")

    # Compute tool hit counts
    df = _compute_tool_counts(df)

    # Write output
    df.to_csv(args.output, index=False)
    logger.info(f"Wrote {len(df)} rows x {len(df.columns)} columns to {args.output}")


if __name__ == '__main__':
    main()
