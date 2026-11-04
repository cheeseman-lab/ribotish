"""Edge cases of `predict -p N`: tiny inputs, early exit, cached TIS background."""

import shutil
import subprocess
from pathlib import Path

import pytest

from equivalence import cached_reference

TIMEOUT = 900  # a hang (lost pool result, zombie worker) fails instead of blocking the suite


def _bams(fixture_dir):
    return {p.name[:3]: p for p in sorted(Path(fixture_dir).glob("72?_*.bam"))}


def _predict(fork_bin, out_dir, tis, chx, gtf, fa, nproc, extra=()):
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(fork_bin),
        "predict",
        "-t",
        ",".join(map(str, tis)),
        "-b",
        ",".join(map(str, chx)),
        "-g",
        str(gtf),
        "-f",
        str(fa),
        "-o",
        str(out_dir / "out.txt"),
        "--minaalen",
        "3",
        "--alt",
        "--seq",
        "--aaseq",
        "-p",
        str(nproc),
        *extra,
    ]
    r = subprocess.run(cmd, cwd=out_dir, capture_output=True, text=True, timeout=TIMEOUT)
    return r.returncode, r.stdout + r.stderr


def _small_gtf(mini_dir, dest, ngenes, genes=None):
    """First `ngenes` genes of the mini GTF, or of `genes` if given (fewer genes than workers)."""
    keep, out = set(), []
    for line in (Path(mini_dir) / "chr22.gtf").read_text().splitlines(keepends=True):
        if line.startswith("#") or 'gene_id "' not in line:
            continue
        gid = line.split('gene_id "', 1)[1].split('"', 1)[0]
        if genes is not None and gid not in genes:
            continue
        if gid not in keep and len(keep) == ngenes:
            continue
        keep.add(gid)
        out.append(line)
    dest.write_text("".join(out))
    return dest


def _called_genes(out_txt):
    """Gene ids with calls in a predict out.txt, in file order."""
    lines = Path(out_txt).read_text().splitlines()[1:]
    return list(dict.fromkeys(line.split("\t", 1)[0] for line in lines))


def test_more_workers_than_genes(fork_bin, upstream_bin, mini_dir, ref_root, tmp_path):
    bg = cached_reference(upstream_bin, mini_dir, ref_root) / "bg.txt"
    called = _called_genes(bg.parent / "out.txt")[:3]
    gtf = _small_gtf(mini_dir, tmp_path / "small.gtf", 3, set(called))
    b = _bams(mini_dir)
    tis, chx = [b["727"], b["729"]], [b["726"], b["728"]]
    fa, extra = mini_dir / "chr22.fa", ["-s", str(bg)]
    ref = tmp_path / "upstream_p1"
    rc, log = _predict(upstream_bin, ref, tis, chx, gtf, fa, 1, extra)
    assert rc == 0, log[-2000:]
    for name in ("out.txt", "out_all.txt"):
        assert len((ref / name).read_text().splitlines()) > 1, f"{name}: no data rows"
    for n in (1, 2, 16):
        out = tmp_path / f"p{n}"
        rc, log = _predict(fork_bin, out, tis, chx, gtf, fa, n, extra)
        assert rc == 0, log[-2000:]
        for name in ("out.txt", "out_all.txt"):
            assert (ref / name).read_bytes() == (out / name).read_bytes(), (n, name)


@pytest.mark.parametrize("nproc", [1, 2])
def test_too_few_background_samples_exits(fork_bin, mini_dir, tmp_path, nproc):
    # Background estimation calls exit(1) inside the child process when -p > 1; must not hang.
    gtf = _small_gtf(mini_dir, tmp_path / "small.gtf", 3)
    b = _bams(mini_dir)
    rc, log = _predict(
        fork_bin,
        tmp_path / "out",
        [b["727"], b["729"]],
        [b["726"], b["728"]],
        gtf,
        mini_dir / "chr22.fa",
        nproc,
        ["-e", str(tmp_path / "out" / "bg.txt")],
    )
    assert rc == 1, log[-2000:]
    assert "Too few TIS background samples" in log


def test_cached_background_warning(fork_bin, upstream_bin, mini_dir, ref_root, tmp_path):
    bg = cached_reference(upstream_bin, mini_dir, ref_root) / "bg.txt"
    src = _bams(mini_dir)["727"]
    data = tmp_path / "data"
    data.mkdir()
    for suffix in ("", ".bai", ".para.py"):
        (data / (src.name + suffix)).symlink_to(Path(str(src) + suffix))
    tis = data / src.name
    cached = Path(str(tis) + ".bgest.txt")
    shutil.copyfile(bg, cached)
    chx = [_bams(mini_dir)["726"]]
    gtf, fa = mini_dir / "chr22.gtf", mini_dir / "chr22.fa"

    rc, log = _predict(fork_bin, tmp_path / "implicit", [tis], chx, gtf, fa, 2)
    assert rc == 0, log[-2000:]
    assert f"WARNING: loading cached TIS background {cached}" in log

    # -e does not stop the cache from loading (upstream behaviour); the warning must still print.
    rc, log = _predict(
        fork_bin, tmp_path / "with_e", [tis], chx, gtf, fa, 2, ["-e", str(tmp_path / "bg_e.txt")]
    )
    assert rc == 0, log[-2000:]
    assert f"WARNING: loading cached TIS background {cached}" in log

    rc, log = _predict(fork_bin, tmp_path / "explicit", [tis], chx, gtf, fa, 2, ["-s", str(cached)])
    assert rc == 0, log[-2000:]
    assert "WARNING" not in log
    for name in ("out.txt", "out_all.txt"):
        assert (tmp_path / "implicit" / name).read_bytes() == (
            tmp_path / "explicit" / name
        ).read_bytes()
