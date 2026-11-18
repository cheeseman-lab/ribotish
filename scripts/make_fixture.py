"""Build a small, self-contained Ribo-TISH fixture from one genomic region.

The region is cut out of the genome and becomes its own contig, so the fixture does not need the
full chromosome. Coordinates in the BAMs and the GTF are shifted by the region start. A read is kept
only when it lies fully inside the region. A gene is kept only when all its rows lie inside it. BAM
headers keep only @HD and @SQ, and read names are replaced by numbers. The equivalence test runs the
fork and pinned upstream on the same fixture, so the cut does not affect what it checks.

Usage:
    python scripts/make_fixture.py OUT_DIR --fasta GENOME.fa --gtf ANNOTATION.gtf
        --region chr22:19000000-21000000 --tis A.bam[,B.bam] --chx C.bam[,D.bam]
        [--para-dir DIR] [--name LABEL ...]

Each BAM needs its ``ribotish quality`` offsets: ``<bam>.para.py`` next to it, or
``<para-dir>/<library>.para.py`` where the library is the BAM name before ``_Aligned``. OUT_DIR gets
``region.fa``, ``region.gtf``, one BAM, index and ``.para.py`` per library, and ``fixture.tsv``
(role, file). Needs pysam and samtools on PATH.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

import pysam


def parse_region(text):
    chrom, span = text.split(":")
    start, end = (int(x.replace(",", "")) for x in span.split("-"))
    return chrom, start, end


def write_fasta(fasta, chrom, start, end, out):
    seq = pysam.FastaFile(fasta).fetch(chrom, start, end)
    with open(out, "w") as fh:
        fh.write(f">{chrom}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i : i + 60] + "\n")
    pysam.faidx(str(out))
    return len(seq)


def write_gtf(gtf, chrom, start, end, out):
    rows, span = [], {}
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] != chrom:
                continue
            m = re.search(r'gene_id "([^"]+)"', f[8])
            if not m:
                continue
            g, a, b = m.group(1), int(f[3]), int(f[4])
            lo, hi = span.get(g, (a, b))
            span[g] = (min(lo, a), max(hi, b))
            rows.append((g, f))
    # GTF is 1-based closed; the region is 0-based half-open
    keep = {g for g, (lo, hi) in span.items() if lo > start and hi <= end}
    with open(out, "w") as fh:
        for g, f in rows:
            if g in keep:
                f[3], f[4] = str(int(f[3]) - start), str(int(f[4]) - start)
                fh.write("\t".join(f) + "\n")
    return len(keep)


def write_bam(bam, chrom, start, end, length, out):
    header = {"HD": {"VN": "1.4", "SO": "coordinate"}, "SQ": [{"SN": chrom, "LN": length}]}
    n = 0
    with pysam.AlignmentFile(bam) as src, pysam.AlignmentFile(out, "wb", header=header) as dst:
        for r in src.fetch(chrom, start, end):
            if r.reference_start < start or r.reference_end > end:
                continue
            n += 1
            d = r.to_dict()
            d["name"], d["ref_name"] = str(n), chrom
            d["ref_pos"] = str(r.reference_start - start + 1)
            a = pysam.AlignedSegment.from_dict(d, dst.header)
            dst.write(a)
    pysam.index(str(out))
    return n


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("out_dir", type=Path)
    p.add_argument("--fasta", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--region", required=True)
    p.add_argument("--tis", required=True, help="comma-separated initiation BAMs")
    p.add_argument("--chx", required=True, help="comma-separated elongation BAMs")
    p.add_argument("--para-dir", type=Path, help="offset files, when not next to the BAMs")
    p.add_argument("--name", nargs="*", default=None, help="output names, TIS then CHX order")
    a = p.parse_args(argv)

    chrom, start, end = parse_region(a.region)
    out = a.out_dir
    out.mkdir(parents=True, exist_ok=True)
    length = write_fasta(a.fasta, chrom, start, end, out / "region.fa")
    genes = write_gtf(a.gtf, chrom, start, end, out / "region.gtf")
    print(f"{a.region}: {length} nt, {genes} genes", file=sys.stderr)

    libs = [("tis", b) for b in a.tis.split(",")] + [("chx", b) for b in a.chx.split(",")]
    names = a.name or [Path(b).name.split("_")[0] for _, b in libs]
    if len(names) != len(libs):
        p.error("--name needs one name per BAM")
    rows = []
    for (role, bam), name in zip(libs, names):
        dst = out / f"{name}.bam"
        n = write_bam(bam, chrom, start, end, length, dst)
        para = Path(bam + ".para.py")
        if not para.exists() and a.para_dir:
            para = a.para_dir / (Path(bam).name.split("_Aligned")[0] + ".para.py")
        offdict = [ln for ln in para.read_text().splitlines() if ln.startswith("offdict")]
        (out / f"{name}.bam.para.py").write_text(offdict[-1] + "\n")
        rows.append(f"{role}\t{dst.name}")
        print(f"{name} ({role}): {n} reads", file=sys.stderr)
    (out / "fixture.tsv").write_text("role\tbam\n" + "\n".join(rows) + "\n")
    subprocess.run(["du", "-sh", str(out)], check=True)


if __name__ == "__main__":
    main()
