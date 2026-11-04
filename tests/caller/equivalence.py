"""Helpers to run ribotish predict and compare outputs byte-for-byte."""

import hashlib
import shutil
import subprocess
from pathlib import Path

OUTPUTS = ("out.txt", "out_all.txt", "bg.txt")


def run_predict(ribotish_bin, fixture_dir, out_dir, nproc):
    """Run `ribotish predict` on a fixture dir; return out_dir."""
    fixture_dir = Path(fixture_dir).resolve()
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    bams = {p.name[:3]: p for p in sorted(fixture_dir.glob("72?_*.bam"))}
    tis = ",".join(str(bams[k]) for k in ("727", "729"))
    chx = ",".join(str(bams[k]) for k in ("726", "728"))
    cmd = [
        str(ribotish_bin),
        "predict",
        "-t",
        tis,
        "-b",
        chx,
        "-g",
        str(fixture_dir / "chr22.gtf"),
        "-f",
        str(fixture_dir / "chr22.fa"),
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
    for p in sorted(Path(fixture_dir).iterdir()):
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
