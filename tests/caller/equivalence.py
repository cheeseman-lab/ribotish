"""Helpers to run ribotish predict and compare outputs byte-for-byte."""

import hashlib
import shutil
import subprocess
from pathlib import Path

OUTPUTS = ("out.txt", "out_all.txt", "bg.txt")


def load_fixture(fixture_dir):
    """TIS BAMs, CHX BAMs, GTF and FASTA of a fixture directory.

    A fixture lists its BAMs in fixture.tsv (role, bam) and uses ../region.gtf and ../region.fa
    (scripts/make_fixture.py). A directory without fixture.tsv has the local chr22 layout:
    72?_*.bam, chr22.gtf and chr22.fa.
    """
    d = Path(fixture_dir).resolve()
    if (d / "fixture.tsv").exists():
        rows = [ln.split("\t") for ln in (d / "fixture.tsv").read_text().splitlines()[1:] if ln]
        tis = [d / b for role, b in rows if role == "tis"]
        chx = [d / b for role, b in rows if role == "chx"]
        return tis, chx, d.parent / "region.gtf", d.parent / "region.fa"
    bams = {p.name[:3]: p for p in sorted(d.glob("72?_*.bam"))}
    return [bams["727"], bams["729"]], [bams["726"], bams["728"]], d / "chr22.gtf", d / "chr22.fa"


def run_predict(ribotish_bin, fixture_dir, out_dir, nproc):
    """Run `ribotish predict` on a fixture dir; return out_dir."""
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    tis_bams, chx_bams, gtf, fa = load_fixture(fixture_dir)
    tis = ",".join(map(str, tis_bams))
    chx = ",".join(map(str, chx_bams))
    cmd = [
        str(ribotish_bin),
        "predict",
        "-t",
        tis,
        "-b",
        chx,
        "-g",
        str(gtf),
        "-f",
        str(fa),
        "-o",
        str(out_dir / "out.txt"),
        "-e",
        str(out_dir / "bg.txt"),
        "--minaalen",
        "3",
        "--alt",
        "--seq",
        "--aaseq",
        "-p",
        str(nproc),
    ]
    log = out_dir / "predict.log"
    with open(log, "w") as fh:
        r = subprocess.run(cmd, cwd=out_dir, stdout=fh, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        raise RuntimeError(f"predict failed (rc={r.returncode}); tail:\n" + log.read_text()[-2000:])
    return out_dir


def _stamp(upstream_bin, fixture_dir):
    h = hashlib.sha256()
    _, _, gtf, fa = load_fixture(fixture_dir)
    for p in [*sorted(Path(fixture_dir).iterdir()), gtf, fa]:
        st = p.stat()  # follows symlinks (chr22.fa)
        h.update(f"{p.name}:{st.st_size}:{int(st.st_mtime)}\n".encode())
    ver = subprocess.run([str(upstream_bin), "--version"], capture_output=True, text=True)
    h.update(f"{Path(upstream_bin).resolve()}\n{ver.stdout}{ver.stderr}".encode())
    return h.hexdigest()


def cached_reference(upstream_bin, fixture_dir, ref_root):
    """Upstream -p 1 outputs for a fixture, cached under ref_root/<fixture name>/.

    The cache is rebuilt when the fixture files or the upstream binary change (.stamp).
    """
    ref = Path(ref_root) / Path(fixture_dir).name
    stamp = _stamp(upstream_bin, fixture_dir)
    sf = ref / ".stamp"
    ok = sf.exists() and sf.read_text() == stamp and all((ref / f).exists() for f in OUTPUTS)
    if not ok:
        if ref.exists():
            shutil.rmtree(ref)
        run_predict(upstream_bin, fixture_dir, ref, 1)
        sf.write_text(stamp)
    return ref


def assert_identical(a_dir, b_dir):
    """Assert out.txt, out_all.txt and bg.txt are byte-identical in both dirs."""
    bad = []
    for name in OUTPUTS:
        a, b = Path(a_dir) / name, Path(b_dir) / name
        if not a.exists() or not b.exists():
            bad.append(f"{name}: missing ({a.exists()=}, {b.exists()=})")
        elif a.stat().st_size == 0 or b.stat().st_size == 0:
            bad.append(f"{name}: empty output ({a.stat().st_size=}, {b.stat().st_size=})")
        elif a.read_bytes() != b.read_bytes():
            la, lb = a.read_text().splitlines(), b.read_text().splitlines()
            for i in range(max(len(la), len(lb))):
                x = la[i] if i < len(la) else "<EOF>"
                y = lb[i] if i < len(lb) else "<EOF>"
                if x != y:
                    bad.append(f"{name}: line {i + 1} differs\n  a: {x}\n  b: {y}")
                    break
            else:
                bad.append(f"{name}: differs (line endings/trailing bytes)")
    assert not bad, "outputs differ:\n" + "\n".join(bad)
