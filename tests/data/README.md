# Test fixtures

`chr22_19M/` contains the region chr22:19,000,000-21,000,000 of GRCh38, cut out as its own
2,000,000 nt contig named `chr22`. Position 1 of the contig is chr22:19,000,001.

| File | Contents |
|---|---|
| `region.fa` | GRCh38 primary assembly, the region only |
| `region.gtf` | GENCODE v50 primary assembly annotation. 130 genes, each with all of its rows inside the region; coordinates shifted by -19,000,000 |
| `<fixture>/fixture.tsv` | role (`tis` or `chx`) and BAM of each library |
| `<fixture>/<library>.bam`, `.bai` | reads that lie fully inside the region, shifted by -19,000,000 |
| `<fixture>/<library>.bam.para.py` | the library's P-site offsets (`ribotish quality`, curated) |

Built by `scripts/make_fixture.py`:
- The BAM headers keep only `@HD` and `@SQ`.
- Read names are replaced by numbers.
- Alignments, flags and tags are unchanged.

The reads were aligned with STAR 2.7.1a to GRCh38 + GENCODE v50, with
`--outFilterMultimapNmax 5 --outReadsUnmapped Fastx --outFilterType BySJout --outSAMattributes All
--outSAMtype BAM SortedByCoordinate`, after per-source adapter trimming.

## Fixtures

| Fixture | Library | Role | Source | Reads |
|---|---|---|---|---|
| `u2os` | `727_S24_L002` | tis (harringtonine), rep 1 | U2OS, Cheeseman lab (Whitehead Institute), unpublished | 21,982 |
| `u2os` | `729_S26_L002` | tis (harringtonine), rep 2 | same | 17,317 |
| `u2os` | `726_S23_L002` | chx, rep 1 | same | 71,033 |
| `u2os` | `728_S25_L002` | chx, rep 2 | same | 62,005 |
| `chen2020_ipsc` | `SRR9113062` | tis (harringtonine) | human iPSC, Chen et al. 2020, GEO GSE131650 (GSM3791723) | 93,180 |
| `chen2020_ipsc` | `SRR9113064` | chx, rep 1 | same study (GSM3791725) | 46,610 |
| `chen2020_ipsc` | `SRR9113065` | chx, rep 2 | same study (GSM3791726) | 43,433 |

Chen J, Brunner AD, Cogan JZ, et al. Pervasive functional translation of noncanonical human open
reading frames. *Science* 367, 1140-1146 (2020).

The fixtures are for testing this software only. The equivalence test runs the fork and pinned
upstream Ribo-TISH 0.2.8 on each fixture and requires byte-identical output, so the cut changes
nothing that the test checks.
