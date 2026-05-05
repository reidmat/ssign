"""Regression + parser tests for run_deepsece.py.

Pins T3SS preservation through parse_deepsece_output — downstream code
(cross_validate_predictions, proximity_analysis) drops T3SS calls
intentionally, but the wrapper itself must not silently filter them, or
the audit trail vanishes. Downstream T3SS guards are covered by
test_proximity_analysis.py and test_cross_validate_predictions.py.
"""

import os
import sys

# Production module
SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src", "ssign_app", "scripts"))
sys.path.insert(0, SCRIPTS_DIR)

from _helpers import DSE_RAW_FIELDS, make_dse_raw_row, write_tsv  # noqa: E402
from run_deepsece import PREDICTED_LABELS, SS_MAP, parse_deepsece_output  # noqa: E402


def _write_dse_csv(tmp_dir, rows):
    return write_tsv(os.path.join(tmp_dir, "dse.csv"), DSE_RAW_FIELDS, rows, delimiter=",")


class TestT3ssRowsPreserved:
    def test_t3ss_row_present_after_parse(self, tmp_dir):
        path = _write_dse_csv(
            tmp_dir,
            [
                make_dse_raw_row("GENE_0001", ss_type="T1SS", max_prob=0.95, T1_prob=0.95),
                make_dse_raw_row("GENE_0002", ss_type="T3SS", max_prob=0.91, T3_prob=0.91),
                make_dse_raw_row("GENE_0003", ss_type="Non-secreted", nonsec_prob=0.99),
            ],
        )
        by_locus = {e["locus_tag"]: e for e in parse_deepsece_output(path)}
        assert by_locus["GENE_0002"]["dse_ss_type"] == "T3SS"
        assert by_locus["GENE_0002"]["dse_max_prob"] == "0.91"

    def test_all_ss_type_calls_pass_through(self, tmp_dir):
        rows = [
            make_dse_raw_row(f"GENE_{i:04d}", ss_type=ss_type, max_prob=0.9)
            for i, ss_type in enumerate(["Non-secreted", "T1SS", "T2SS", "T3SS", "T4SS", "T6SS"])
        ]
        entries = parse_deepsece_output(_write_dse_csv(tmp_dir, rows))
        assert len(entries) == 6
        assert {e["dse_ss_type"] for e in entries} == {
            "Non-secreted",
            "T1SS",
            "T2SS",
            "T3SS",
            "T4SS",
            "T6SS",
        }


class TestColumnMapping:
    def test_protein_id_renamed_to_locus_tag(self, tmp_dir):
        path = _write_dse_csv(tmp_dir, [make_dse_raw_row("BIMENO_04457", ss_type="T1SS")])
        entry = parse_deepsece_output(path)[0]
        assert entry["locus_tag"] == "BIMENO_04457"
        assert "protein_id" not in entry

    def test_t_prob_columns_renamed_with_dse_prefix(self, tmp_dir):
        path = _write_dse_csv(tmp_dir, [make_dse_raw_row("GENE_0001", T3_prob=0.7, T6_prob=0.3)])
        entry = parse_deepsece_output(path)[0]
        assert entry["dse_T3_prob"] == "0.7"
        assert entry["dse_T6_prob"] == "0.3"


class TestModelLabelContract:
    """If the upstream model retrains and label order shifts, every downstream
    T3SS guard breaks silently. Pin the contract."""

    def test_predicted_labels_position_3_is_t3ss(self):
        assert PREDICTED_LABELS[3] == "III"
        assert SS_MAP["III"] == "T3SS"

    def test_ss_map_covers_every_predicted_label(self):
        for label in PREDICTED_LABELS:
            assert label in SS_MAP, f"PREDICTED_LABEL {label!r} missing from SS_MAP"

    def test_ss_map_t_indices_match_t_prob_columns(self):
        # dse_T3_prob means "T3SS probability" downstream — pin the alignment.
        assert SS_MAP["I"] == "T1SS"
        assert SS_MAP["II"] == "T2SS"
        assert SS_MAP["III"] == "T3SS"
        assert SS_MAP["IV"] == "T4SS"
        assert SS_MAP["VI"] == "T6SS"
