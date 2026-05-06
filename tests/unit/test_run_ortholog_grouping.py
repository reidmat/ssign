"""Tests for run_ortholog_grouping.py.

Three pure-Python surfaces here:

1. `read_fasta_simple` — header-only-first-token FASTA reader (no
   BioPython dep). Multi-line sequences must concatenate.
2. `cluster_union_find` — single-linkage clustering via path-compressed
   union-find. Hits between known IDs unify; transitivity must hold;
   IDs not mentioned in any hit must end up as their own singleton.
3. `compute_group_stats` — per-group `OG_NNN` IDs by descending size,
   mean within-group %identity (bidirectional hit lookup), singleton
   default identity = 100.

The `run_local_blast` subprocess path requires NCBI BLAST+ on PATH and
is exercised by tests/integration/.
"""

import os
import sys

SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src", "ssign_app", "scripts"))
sys.path.insert(0, SCRIPTS_DIR)

from run_ortholog_grouping import (  # noqa: E402
    cluster_union_find,
    compute_group_stats,
    read_fasta_simple,
)

from ssign_app.scripts.ssign_lib.fasta_io import write_fasta  # noqa: E402

# ---------------------------------------------------------------------------
# read_fasta_simple
# ---------------------------------------------------------------------------


class TestReadFastaSimple:
    def test_single_sequence(self, tmp_dir):
        path = os.path.join(tmp_dir, "x.fasta")
        write_fasta({"GENE_001 some annotation": "MKTLLLTLLCAFSV"}, path)
        assert read_fasta_simple(path) == {"GENE_001": "MKTLLLTLLCAFSV"}

    def test_multi_line_sequence_concatenated(self, tmp_dir):
        # write_fasta wraps at line_width=80; supply a 140-char seq so the
        # output spans two lines and the parser's concatenation logic gets
        # exercised.
        seq = "MKTLLLTLLCAFSV" * 10
        path = os.path.join(tmp_dir, "x.fasta")
        write_fasta({"GENE_001": seq}, path)
        assert read_fasta_simple(path) == {"GENE_001": seq}

    def test_id_is_first_whitespace_token(self, tmp_dir):
        path = os.path.join(tmp_dir, "x.fasta")
        write_fasta(
            {"GENE_001 [Escherichia coli] hypothetical protein": "MKT"},
            path,
        )
        assert "GENE_001" in read_fasta_simple(path)

    def test_multiple_sequences(self, tmp_dir):
        path = os.path.join(tmp_dir, "x.fasta")
        write_fasta({"A": "AAAA", "B": "BBBB", "C": "CCCC"}, path)
        assert read_fasta_simple(path) == {"A": "AAAA", "B": "BBBB", "C": "CCCC"}

    def test_blank_lines_ignored(self, tmp_dir):
        # write_fasta won't emit blank lines — write raw to test the parser's
        # blank-line handling.
        path = os.path.join(tmp_dir, "x.fasta")
        with open(path, "w") as f:
            f.write(">A\n\nMKT\n\n>B\nGGG\n")
        assert read_fasta_simple(path) == {"A": "MKT", "B": "GGG"}

    def test_empty_file_returns_empty(self, tmp_dir):
        path = os.path.join(tmp_dir, "empty.fasta")
        open(path, "w").close()
        assert read_fasta_simple(path) == {}

    def test_orphan_lines_before_first_header_ignored(self, tmp_dir):
        # Garbage prefix before any `>` — must be silently dropped
        path = os.path.join(tmp_dir, "x.fasta")
        with open(path, "w") as f:
            f.write("garbage line\n>A\nMKT\n")
        assert read_fasta_simple(path) == {"A": "MKT"}


# ---------------------------------------------------------------------------
# cluster_union_find
# ---------------------------------------------------------------------------


class TestClusterUnionFind:
    def test_no_hits_yields_all_singletons(self):
        groups = cluster_union_find(hits=[], all_protein_ids={"A", "B", "C"})
        # Every protein its own group
        sizes = sorted(len(s) for s in groups.values())
        assert sizes == [1, 1, 1]

    def test_pairwise_hit_unifies_two(self):
        groups = cluster_union_find(
            hits=[("A", "B", 80.0, 80.0)],
            all_protein_ids={"A", "B", "C"},
        )
        # 2 groups: {A,B} and {C}
        sizes = sorted(len(s) for s in groups.values())
        assert sizes == [1, 2]

    def test_transitive_closure(self):
        # A-B and B-C → all three in one group
        groups = cluster_union_find(
            hits=[("A", "B", 80.0, 80.0), ("B", "C", 80.0, 80.0)],
            all_protein_ids={"A", "B", "C"},
        )
        assert len(groups) == 1
        assert next(iter(groups.values())) == {"A", "B", "C"}

    def test_two_disconnected_components(self):
        groups = cluster_union_find(
            hits=[
                ("A", "B", 80.0, 80.0),
                ("C", "D", 80.0, 80.0),
            ],
            all_protein_ids={"A", "B", "C", "D"},
        )
        member_sets = sorted([frozenset(s) for s in groups.values()], key=lambda s: sorted(s)[0])
        assert member_sets == [frozenset({"A", "B"}), frozenset({"C", "D"})]

    def test_hits_referencing_unknown_ids_skipped(self):
        # "X" is not in all_protein_ids → the hit should silently no-op
        groups = cluster_union_find(
            hits=[("A", "X", 80.0, 80.0)],
            all_protein_ids={"A", "B"},
        )
        # Both A and B remain as singletons
        sizes = sorted(len(s) for s in groups.values())
        assert sizes == [1, 1]

    def test_self_hit_idempotent(self):
        # A union(A, A) shouldn't change cluster structure
        groups = cluster_union_find(
            hits=[("A", "A", 100.0, 100.0)],
            all_protein_ids={"A", "B"},
        )
        sizes = sorted(len(s) for s in groups.values())
        assert sizes == [1, 1]

    def test_chain_of_five(self):
        # A-B-C-D-E chain → one group of 5
        hits = [
            ("A", "B", 80.0, 80.0),
            ("B", "C", 80.0, 80.0),
            ("C", "D", 80.0, 80.0),
            ("D", "E", 80.0, 80.0),
        ]
        groups = cluster_union_find(hits=hits, all_protein_ids={"A", "B", "C", "D", "E"})
        assert len(groups) == 1
        assert next(iter(groups.values())) == {"A", "B", "C", "D", "E"}


# ---------------------------------------------------------------------------
# compute_group_stats
# ---------------------------------------------------------------------------


class TestComputeGroupStats:
    def test_groups_ordered_by_size_desc(self):
        # 3-member group first (OG_001), then two singletons
        groups = {
            "rep1": {"A", "B", "C"},
            "rep2": {"D"},
            "rep3": {"E"},
        }
        stats = compute_group_stats(
            groups,
            hits=[("A", "B", 90.0, 90.0), ("B", "C", 95.0, 95.0)],
            all_protein_ids={"A", "B", "C", "D", "E"},
        )
        # First entry is the 3-member group
        assert stats[0]["ortholog_group"] == "OG_001"
        assert stats[0]["n_members"] == 3
        # Subsequent entries are singletons (in some order, both n_members=1)
        assert all(s["n_members"] == 1 for s in stats[1:])

    def test_singleton_default_identity_100(self):
        groups = {"rep1": {"A"}}
        stats = compute_group_stats(
            groups,
            hits=[],
            all_protein_ids={"A"},
        )
        assert stats[0]["mean_pident"] == 100.0
        assert stats[0]["n_members"] == 1

    def test_mean_identity_within_group(self):
        # Group {A, B, C} with within-group hits at 80 and 90% → mean 85.0
        groups = {"rep1": {"A", "B", "C"}}
        stats = compute_group_stats(
            groups,
            hits=[
                ("A", "B", 80.0, 80.0),
                ("B", "C", 90.0, 90.0),
                # Cross-group hit (irrelevant — A and B both in same group)
            ],
            all_protein_ids={"A", "B", "C"},
        )
        assert stats[0]["mean_pident"] == 85.0

    def test_bidirectional_hit_lookup(self):
        # Hit recorded only as (A, B); compute_group_stats must also find (B, A)
        groups = {"rep1": {"A", "B"}}
        stats = compute_group_stats(
            groups,
            hits=[("A", "B", 75.0, 80.0)],
            all_protein_ids={"A", "B"},
        )
        assert stats[0]["mean_pident"] == 75.0

    def test_members_field_sorted_and_semicolon_joined(self):
        groups = {"rep1": {"C", "A", "B"}}
        stats = compute_group_stats(
            groups,
            hits=[],
            all_protein_ids={"A", "B", "C"},
        )
        assert stats[0]["members"] == "A;B;C"

    def test_og_id_zero_padded(self):
        # 3-digit zero-padding contract pinned: OG_001, OG_002, ...
        groups = {f"rep{i}": {f"P{i}"} for i in range(1, 6)}
        stats = compute_group_stats(
            groups,
            hits=[],
            all_protein_ids={f"P{i}" for i in range(1, 6)},
        )
        og_ids = [s["ortholog_group"] for s in stats]
        assert og_ids == ["OG_001", "OG_002", "OG_003", "OG_004", "OG_005"]

    def test_no_hits_in_multi_member_group_falls_back_to_100(self):
        # group {A, B} but no hits between them — defensively defaults to 100
        # rather than crashing on division-by-zero
        groups = {"rep1": {"A", "B"}}
        stats = compute_group_stats(
            groups,
            hits=[],
            all_protein_ids={"A", "B"},
        )
        assert stats[0]["mean_pident"] == 100.0
