"""Regression + parser tests for run_blastp.py.

Pins the NCBI " >" stitle-split behaviour: stitle "hemolysin >ID hypothetical
protein" must reduce to "hemolysin" before the EXCLUDE_TERMS check, otherwise
real hits get wrongly dropped.
"""

import os
import sys

import pytest

SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src", "ssign_app", "scripts"))
sys.path.insert(0, SCRIPTS_DIR)

from _helpers import make_blast_outfmt_row  # noqa: E402
from run_blastp import EXCLUDE_TERMS, filter_hits, parse_blast_tabular  # noqa: E402


@pytest.mark.parametrize(
    "stitle, expected_desc",
    [
        ("hemolysin BL21-like protein >WP_X.1 hypothetical protein", "hemolysin BL21-like protein"),
        ("hemolysin >ID_2 something >ID_3 hypothetical protein", "hemolysin"),
        ("hemolysin", "hemolysin"),
    ],
)
def test_ncbi_title_split_takes_first_segment(stitle, expected_desc):
    hits = parse_blast_tabular(make_blast_outfmt_row(stitle=stitle))
    assert hits["GENE_0001"]["blastp_hit_description"] == expected_desc


def test_description_truncated_to_200_chars():
    long_desc = "x" * 500
    hits = parse_blast_tabular(
        make_blast_outfmt_row(stitle=f"{long_desc} >ID hypothetical protein"),
    )
    assert len(hits["GENE_0001"]["blastp_hit_description"]) == 200


class TestFilterRespectsSplitDescription:
    """Without the split, a real hit whose stitle ends in
    "...>X hypothetical protein" would be wrongly excluded — the unsplit
    string contains "hypothetical protein", an EXCLUDE_TERM."""

    def test_real_hit_kept_when_redundant_tail_says_hypothetical(self):
        hits = parse_blast_tabular(
            make_blast_outfmt_row(stitle="hemolysin >WP_X.1 hypothetical protein"),
        )
        assert "GENE_0001" in filter_hits(hits, min_pident=80, min_qcov=80)

    def test_truly_hypothetical_hit_dropped(self):
        hits = parse_blast_tabular(
            make_blast_outfmt_row(stitle="hypothetical protein >WP_X.1 hemolysin"),
        )
        assert filter_hits(hits, min_pident=80, min_qcov=80) == {}


class TestBitscoreSort:
    """parse_blast_tabular re-sorts hits per query by bitscore — BLAST's row
    order is DB-traversal order, not bitscore order, so without the sort the
    wrong best hit gets reported."""

    def test_best_bitscore_wins(self):
        rows = "\n".join(
            [
                make_blast_outfmt_row(qseqid="Q", sseqid="LOWER", bitscore=100.0),
                make_blast_outfmt_row(qseqid="Q", sseqid="HIGHER", bitscore=900.0),
                make_blast_outfmt_row(qseqid="Q", sseqid="MID", bitscore=500.0),
            ]
        )
        hits = parse_blast_tabular(rows)
        assert hits["Q"]["blastp_hit_accession"] == "HIGHER"


class TestParserResilience:
    def test_empty_output_returns_empty_dict(self):
        assert parse_blast_tabular("") == {}

    def test_blank_lines_ignored(self):
        rows = "\n".join(
            [
                "",
                make_blast_outfmt_row(qseqid="Q1"),
                "",
                make_blast_outfmt_row(qseqid="Q2"),
                "",
            ]
        )
        assert set(parse_blast_tabular(rows).keys()) == {"Q1", "Q2"}

    def test_short_rows_skipped(self):
        # Fewer than the 15 required fields → silently dropped (corrupt output)
        assert parse_blast_tabular("Q1\tWP_X.1\t95.0\t100") == {}


@pytest.mark.parametrize("term", EXCLUDE_TERMS)
def test_every_exclude_term_filters_out(term):
    hits = parse_blast_tabular(make_blast_outfmt_row(stitle=term))
    assert filter_hits(hits, min_pident=0, min_qcov=0) == {}
