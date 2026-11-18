import pytest

from equivalence import assert_identical, cached_reference, run_predict


@pytest.mark.parametrize("nproc", [1, 2, 8])
def test_fixture_matches_upstream(nproc, fork_bin, upstream_bin, fixture_dir, ref_root, tmp_path):
    ref = cached_reference(upstream_bin, fixture_dir, ref_root)
    out = run_predict(fork_bin, fixture_dir, tmp_path / f"p{nproc}", nproc)
    assert_identical(ref, out)


@pytest.mark.slow
def test_chr22_matches_upstream(fork_bin, upstream_bin, chr22_dir, ref_root, tmp_path):
    ref = cached_reference(upstream_bin, chr22_dir, ref_root)
    out = run_predict(fork_bin, chr22_dir, tmp_path / "p8", 8)
    assert_identical(ref, out)
