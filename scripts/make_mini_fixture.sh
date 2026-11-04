#!/bin/bash
# Build work/mini/ from work/chr22: BAMs restricted to a 2 Mb region, genes fully inside.
set -euo pipefail
cd "$(dirname "$0")/.."
SRC=work/chr22
OUT=work/mini
REG=chr22:19000000-21000000
S=19000000
E=21000000
rm -rf "$OUT"; mkdir -p "$OUT"
for bam in "$SRC"/72[6-9]_*.bam; do
  b=$(basename "$bam")
  samtools view -b "$bam" "$REG" -o "$OUT/$b"
  samtools index "$OUT/$b"
  cp "$bam.para.py" "$OUT/$b.para.py"
done
# genes fully inside the region: keep all rows whose gene_id belongs to such a gene
python - "$SRC/chr22.gtf" "$OUT/chr22.gtf" $S $E <<'PY'
import re, sys
src, dst, s, e = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
rows, span = [], {}
for line in open(src):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    m = re.search(r'gene_id "([^"]+)"', f[8])
    g = m.group(1) if m else None
    rows.append((g, line))
    if g:
        a, b = int(f[3]), int(f[4])
        lo, hi = span.get(g, (a, b))
        span[g] = (min(lo, a), max(hi, b))
keep = {g for g, (lo, hi) in span.items() if lo >= s and hi <= e}
with open(dst, "w") as out:
    for g, line in rows:
        if g in keep:
            out.write(line)
print(len(keep), "genes kept")
PY
ln -s "$(readlink -f $SRC/chr22.fa)" "$OUT/chr22.fa"
ln -s "$(readlink -f $SRC/chr22.fa.fai)" "$OUT/chr22.fa.fai"
rm -rf work/ref/mini
du -sh "$OUT"
