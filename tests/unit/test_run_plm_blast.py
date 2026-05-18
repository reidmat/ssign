"""Unit tests for run_plm_blast.py — CSV parser.

Does not exercise the `plmblast.py` subprocess itself (that requires the
pLM-BLAST install plus a ~20 GB ECOD70 database). Tests target the
pure-Python helper that parses pLM-BLAST's CSV output.
"""

import os

from run_plm_blast import (
    _use_cuda_for_embedding,
    load_substrate_ids,
    parse_plmblast_csv,
    write_substrates_only_fasta,
)

# Representative pLM-BLAST CSV fixture. Real output columns verified on
# first integration run — swap in the actual captured output as a fixture
# once we run pLM-BLAST against the T1SS fixture on CX3.
_PLMBLAST_FIXTURE = (
    "qid,sid,score,qstart,qend,tstart,tend\n"
    "GENE_00001,ecod_1a2bA1,0.912,10,145,5,140\n"
    "GENE_00001,ecod_3c4dB2,0.871,20,150,10,140\n"
    "GENE_00003,ecod_5e6fC1,0.789,1,200,1,199\n"
)


class TestParsePlmblastCsv:
    def test_extracts_all_hits(self, tmp_dir):
        path = os.path.join(tmp_dir, "plm_blast.csv")
        with open(path, "w") as f:
            f.write(_PLMBLAST_FIXTURE)

        entries = parse_plmblast_csv(path)
        assert len(entries) == 3

    def test_multiple_hits_per_protein_preserved(self, tmp_dir):
        path = os.path.join(tmp_dir, "plm_blast.csv")
        with open(path, "w") as f:
            f.write(_PLMBLAST_FIXTURE)

        entries = parse_plmblast_csv(path)
        hits_for_gene1 = [e for e in entries if e["protein_id"] == "GENE_00001"]
        assert len(hits_for_gene1) == 2

    def test_captures_all_expected_fields(self, tmp_dir):
        path = os.path.join(tmp_dir, "plm_blast.csv")
        with open(path, "w") as f:
            f.write(_PLMBLAST_FIXTURE)

        entries = parse_plmblast_csv(path)
        first = next(e for e in entries if e["protein_id"] == "GENE_00001")
        assert first["target_id"] == "ecod_1a2bA1"
        assert first["score"] == "0.912"
        assert first["qstart"] == "10"
        assert first["qend"] == "145"
        assert first["tstart"] == "5"
        assert first["tend"] == "140"

    def test_empty_qid_row_skipped(self, tmp_dir):
        path = os.path.join(tmp_dir, "plm_blast.csv")
        with open(path, "w") as f:
            f.write(
                "qid,sid,score,qstart,qend,tstart,tend\n,ecod_empty,0.5,1,10,1,10\nGENE_X,ecod_real,0.9,1,100,1,100\n"
            )

        entries = parse_plmblast_csv(path)
        assert len(entries) == 1
        assert entries[0]["protein_id"] == "GENE_X"

    def test_empty_csv_returns_empty_list(self, tmp_dir):
        path = os.path.join(tmp_dir, "plm_blast.csv")
        with open(path, "w") as f:
            f.write("qid,sid,score,qstart,qend,tstart,tend\n")

        assert parse_plmblast_csv(path) == []

    def test_output_shape_is_stable(self, tmp_dir):
        """Every entry has the six expected fields for downstream consumers."""
        path = os.path.join(tmp_dir, "plm_blast.csv")
        with open(path, "w") as f:
            f.write(_PLMBLAST_FIXTURE)

        entries = parse_plmblast_csv(path)
        required = {
            "protein_id",
            "target_id",
            "score",
            "qstart",
            "qend",
            "tstart",
            "tend",
        }
        for e in entries:
            assert required <= set(e.keys())


class TestSubstrateFiltering:
    """Runs before pLM-BLAST itself; verifies we only search substrates, not whole genome."""

    def test_load_substrate_ids_basic(self, tmp_dir):
        path = os.path.join(tmp_dir, "substrates.tsv")
        with open(path, "w") as f:
            f.write("locus_tag\tsample\n")
            f.write("G1\tsample1\n")
            f.write("G2\tsample1\n")
        assert load_substrate_ids(path) == {"G1", "G2"}

    def test_load_substrate_ids_skips_empty_rows(self, tmp_dir):
        path = os.path.join(tmp_dir, "substrates.tsv")
        with open(path, "w") as f:
            f.write("locus_tag\n")
            f.write("G1\n")
            f.write("\n")
            f.write("G2\n")
        assert load_substrate_ids(path) == {"G1", "G2"}

    def test_write_substrates_only_fasta_filters(self, tmp_dir):
        src = os.path.join(tmp_dir, "proteins.faa")
        with open(src, "w") as f:
            f.write(">G1\nMKT\n>G2\nMFV\n>G3\nMQK\n")
        out = os.path.join(tmp_dir, "substrates.faa")
        n = write_substrates_only_fasta(src, {"G1", "G3"}, out)
        assert n == 2
        with open(out) as f:
            body = f.read()
        assert ">G1" in body and ">G3" in body
        assert ">G2" not in body

    def test_write_substrates_only_fasta_missing_ids_skipped(self, tmp_dir):
        src = os.path.join(tmp_dir, "proteins.faa")
        with open(src, "w") as f:
            f.write(">G1\nMKT\n")
        out = os.path.join(tmp_dir, "substrates.faa")
        n = write_substrates_only_fasta(src, {"G1", "NOT_IN_FASTA"}, out)
        assert n == 1


class TestUseCudaForEmbedding:
    """Regression coverage for the GPU auto-detect that controls --cuda on embeddings.py.

    Pre-fix, the wrapper never passed --cuda regardless of GPU availability,
    causing ProtT5 embedding to run on CPU on full-tier GPU nodes (~100x
    slower). These tests pin the detection contract so a future refactor
    doesn't silently regress it.
    """

    def test_force_cpu_env_var_overrides_detection(self, monkeypatch):
        monkeypatch.setenv("SSIGN_PLMBLAST_FORCE_CPU", "1")
        assert _use_cuda_for_embedding() is False

    def test_force_cpu_accepts_alternate_truthy_values(self, monkeypatch):
        for v in ("true", "TRUE", "yes", "Yes"):
            monkeypatch.setenv("SSIGN_PLMBLAST_FORCE_CPU", v)
            assert _use_cuda_for_embedding() is False, f"failed for value {v!r}"

    def test_force_cpu_falsy_value_does_not_override(self, monkeypatch):
        # An empty value (or "0", "false") should NOT force CPU — auto-detect wins.
        monkeypatch.setenv("SSIGN_PLMBLAST_FORCE_CPU", "0")
        # Result depends on whether torch + CUDA is available; we just verify
        # the env-var didn't short-circuit to False the way "1" does.
        try:
            import torch

            expected = bool(torch.cuda.is_available())
        except ImportError:
            expected = False
        assert _use_cuda_for_embedding() is expected

    def test_missing_torch_falls_back_to_cpu(self, monkeypatch):
        import builtins

        real_import = builtins.__import__

        def fake_import(name, *args, **kwargs):
            if name == "torch":
                raise ImportError("torch not installed (test stub)")
            return real_import(name, *args, **kwargs)

        monkeypatch.delenv("SSIGN_PLMBLAST_FORCE_CPU", raising=False)
        monkeypatch.setattr(builtins, "__import__", fake_import)
        assert _use_cuda_for_embedding() is False
