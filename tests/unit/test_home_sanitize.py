"""Unit tests for Home._sanitize_sample_id.

The sample_id derived from a user-uploaded filename ends up joined into
output paths via `os.path.join(outdir, f"{sid}_results.csv")`. A naive
"strip suffix" approach lets a maliciously-crafted upload write outside
outdir. These tests pin the sanitisation rules so a future refactor
can't silently weaken them.

Security goal: for any user input `s`, `os.path.join(outdir, sid + "_x")`
must resolve to a path strictly inside outdir.
"""

import os

import pytest

from ssign_app.Home import _sanitize_sample_id


class TestSanitizeSampleId:
    def test_normal_filename_passes_through(self):
        assert _sanitize_sample_id("xanthomonas_campestris") == "xanthomonas_campestris"

    def test_dots_dashes_underscores_kept(self):
        assert _sanitize_sample_id("strain.A-1_genomic") == "strain.A-1_genomic"

    @pytest.mark.parametrize(
        "evil",
        [
            "../../../etc/passwd",
            "..\\..\\windows\\system32",
            "/absolute/path/to/file",
            "C:\\Users\\victim\\file",
            "innocent/../../etc/passwd",
            "..",
            "../",
            "../..",
        ],
    )
    def test_path_traversal_neutralised(self, evil):
        sid = _sanitize_sample_id(evil)
        # No directory separators of any flavour
        assert "/" not in sid
        assert "\\" not in sid
        # Joining must stay inside outdir — i.e. abspath(join(outdir, sid))
        # is a child of abspath(outdir).
        outdir = "/tmp/ssign_test_outdir"
        joined = os.path.abspath(os.path.join(outdir, f"{sid}_results.csv"))
        assert joined.startswith(os.path.abspath(outdir) + os.sep)

    def test_leading_dots_stripped(self):
        # Prevents both `..` (parent dir) and `.hidden` (hidden file).
        assert _sanitize_sample_id("..secret") == "secret"
        assert _sanitize_sample_id(".bashrc") == "bashrc"

    def test_runs_of_unsafe_chars_collapse_to_single_underscore(self):
        assert _sanitize_sample_id("a$$$b") == "a_b"
        assert _sanitize_sample_id("a   b") == "a_b"

    def test_empty_or_pure_unsafe_falls_back_to_sample(self):
        assert _sanitize_sample_id("") == "sample"
        assert _sanitize_sample_id("///") == "sample"
        assert _sanitize_sample_id("...") == "sample"

    def test_unicode_chars_replaced(self):
        # Non-ASCII becomes _; existing underscores in the input are safe-
        # set chars and pass through, so adjacent unicode + underscore can
        # produce a double underscore. Fine functionally — the goal is
        # ASCII-safe, not minimum-underscore.
        assert _sanitize_sample_id("escherichia_coli_strain_α") == "escherichia_coli_strain__"

    def test_basename_strips_directories(self):
        # Even with no special chars, a basename with directory components
        # gets stripped of the directory part.
        result = _sanitize_sample_id("subdir/genome.gbff")
        assert result.startswith("genome")
        assert "/" not in result
