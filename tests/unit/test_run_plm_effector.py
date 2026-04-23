"""Unit tests for the PLM-Effector vendored package and its ssign wrapper.

These tests exercise only the pure-Python preprocessing + output-writing
paths. The actual prediction path requires a GPU, ~15 GB of pretrained PLM
weights, plus `transformers` and `xgboost` — that's covered in
tests/integration/test_run_plm_effector_integration.py.
"""

import os
import sys

import numpy as np

SCRIPTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "src", "ssign_app", "scripts")
)
sys.path.insert(0, SCRIPTS_DIR)

from plm_effector.utils import (  # noqa: E402
    _normalise_sequence,
    read_fasta_for_prediction,
    read_fasta_for_prediction_terminal,
)


class TestNormaliseSequence:
    def test_bert_like_inserts_spaces(self):
        assert _normalise_sequence("MKTLLL", "ProtBert") == "M K T L L L"
        assert _normalise_sequence("MKTLLL", "ProtT5") == "M K T L L L"

    def test_esm_keeps_raw(self):
        assert _normalise_sequence("MKTLLL", "esm1") == "MKTLLL"
        assert _normalise_sequence("MKTLLL", "esm2_t33") == "MKTLLL"

    def test_rare_residues_replaced_with_X(self):
        # U, Z, O, B are rare; tokenizers don't know them reliably
        assert _normalise_sequence("MKUZOB", "esm1") == "MKXXXX"
        assert "U" not in _normalise_sequence("MKU", "esm1")


class TestReadFastaForPrediction:
    def _write_fasta(self, tmp_dir, content):
        path = os.path.join(tmp_dir, "sample.fasta")
        with open(path, "w") as f:
            f.write(content)
        return path

    def test_extracts_ids_and_sequences(self, tmp_dir):
        path = self._write_fasta(
            tmp_dir,
            ">GENE_A\nMKTLLLTLLCAFSVAQA\n>GENE_B\nMFVFLVLLPLVSSQ\n",
        )
        ids, seqs = read_fasta_for_prediction(path, model_type="esm1")
        assert ids == [">GENE_A", ">GENE_B"]
        assert seqs == ["MKTLLLTLLCAFSVAQA", "MFVFLVLLPLVSSQ"]

    def test_bert_inserts_spaces_between_residues(self, tmp_dir):
        path = self._write_fasta(tmp_dir, ">G\nMKT\n")
        _, seqs = read_fasta_for_prediction(path, model_type="ProtBert")
        assert seqs == ["M K T"]

    def test_rare_residues_x_substituted(self, tmp_dir):
        path = self._write_fasta(tmp_dir, ">G\nMKUZ\n")
        _, seqs = read_fasta_for_prediction(path, model_type="esm1")
        assert seqs == ["MKXX"]


class TestReadFastaTerminal:
    def test_nterminal_truncates_from_start(self, tmp_dir):
        seq = "M" + "A" * 3000  # 3001 residues
        path = os.path.join(tmp_dir, "long.fasta")
        with open(path, "w") as f:
            f.write(f">LONG\n{seq}\n")

        _, seqs = read_fasta_for_prediction_terminal(
            path, model_type="esm1", terminal="Nterminal", maxlen=1022
        )
        assert len(seqs) == 1
        assert len(seqs[0]) == 1022
        assert seqs[0][0] == "M"  # First residue preserved

    def test_cterminal_truncates_from_end(self, tmp_dir):
        seq = "A" * 3000 + "E"  # Last residue is E
        path = os.path.join(tmp_dir, "long.fasta")
        with open(path, "w") as f:
            f.write(f">LONG\n{seq}\n")

        _, seqs = read_fasta_for_prediction_terminal(
            path, model_type="esm1", terminal="Cterminal", maxlen=1022
        )
        assert len(seqs) == 1
        assert len(seqs[0]) == 1022
        assert seqs[0][-1] == "E"  # Last residue preserved

    def test_short_sequence_passes_through_unchanged(self, tmp_dir):
        path = os.path.join(tmp_dir, "short.fasta")
        with open(path, "w") as f:
            f.write(">SHORT\nMKTLLL\n")

        _, seqs = read_fasta_for_prediction_terminal(
            path, model_type="esm1", terminal="Nterminal", maxlen=1022
        )
        assert seqs == ["MKTLLL"]


class TestWritePredictionsTsv:
    def test_writes_all_proteins_with_threshold_flag(self, tmp_dir):
        from plm_effector.predict_api import write_predictions_tsv

        seq_ids = np.array([">GENE_A", ">GENE_B desc", "GENE_C"])
        stacked = np.array(
            [
                [0.9, 0.8, 0.85, 0.7, 0.6, 0.65],
                [0.1, 0.2, 0.15, 0.3, 0.4, 0.35],
                [0.5, 0.55, 0.6, 0.65, 0.7, 0.75],
            ]
        )
        final_probs = np.array([0.95, 0.15, 0.55])
        passes = np.array([True, False, True])

        out = os.path.join(tmp_dir, "preds.tsv")
        write_predictions_tsv(out, seq_ids, stacked, final_probs, passes, "T1SE", 6)

        with open(out) as f:
            lines = f.readlines()

        # Header + 3 data rows
        assert len(lines) == 4

        header = lines[0].strip().split("\t")
        assert header[0] == "seq_id"
        assert header[1:7] == [f"model{i}" for i in range(1, 7)]
        assert header[7] == "stacking"
        assert header[8] == "passes_threshold"
        assert header[9] == "effector_type"

        row_a = lines[1].strip().split("\t")
        assert row_a[0] == "GENE_A"  # leading ">" stripped
        assert row_a[8] == "1"  # passes
        assert row_a[9] == "T1SE"

        row_b = lines[2].strip().split("\t")
        assert row_b[0] == "GENE_B"  # description after first token dropped
        assert row_b[8] == "0"

        row_c = lines[3].strip().split("\t")
        assert row_c[0] == "GENE_C"  # no leading ">" to strip
        assert row_c[8] == "1"

    def test_t4se_writes_eight_base_models(self, tmp_dir):
        from plm_effector.predict_api import write_predictions_tsv

        stacked = np.array([[0.5] * 8])
        seq_ids = np.array([">X"])
        write_predictions_tsv(
            os.path.join(tmp_dir, "out.tsv"),
            seq_ids,
            stacked,
            np.array([0.7]),
            np.array([True]),
            "T4SE",
            8,
        )
        with open(os.path.join(tmp_dir, "out.tsv")) as f:
            header = f.readline().strip().split("\t")
        assert header[1:9] == [f"model{i}" for i in range(1, 9)]
