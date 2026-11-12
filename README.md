# No ALT, no patch — adaptive-sampling capture-reference emulation

Code and data for the in-silico part of **"No alt, no patch: ALT contigs in the capture reference
deplete MHC coverage under nanopore adaptive sampling."**

Oxford Nanopore adaptive sampling (AS) decides in real time whether to keep or eject each nascent
read by aligning its first few hundred bases against a **capture reference**. A read ejected at that
point is physically lost. A read that is kept is sequenced to full length and is then re-aligned
downstream on the analysis reference, whatever the capture reference said about it.

Two questions must therefore be kept apart:

| | What decides it | Reversible? | Metric here |
|---|---|---|---|
| **Capture** | capture reference + target BED + decision rule | no — the molecule is gone | `keep_rate` |
| **Placement** | alignment reference | yes — re-align | `placement_ok` |

## What this repo contains

```
code/as_emulate.py           emulator: several capture references x decision rules x best_n x prefix length
code/patch_divergence.py     per-contig alignment to the primary region, 10-kb divergence windows
code/run_p14.sh, run_p14_patches.sh   hg38.p14 drivers (gene level; per-patch, Fig. 5)
code/fig5_patches.py         manuscript Fig. 5
code/make_alt_gene_bed.py    genes of the primary assembly that have an ALT / HLA / fix-patch counterpart
code/as_altmapping_demo.py   legacy (V13 Figs. 2 and 5): primary-hit rule only
code/build_t2t_patch_and_run.sh  legacy T2T + patch driver (see caveats below)
data/CKD_SU_V6_no_alt.bed    study capture panel, primary-only
data/CKD_SU_V6_plus_alt.bed  same panel + the ALT intervals of an HLA-aware panel
data/loci_mhc.bed            MHC loci + a non-MHC chr6 control used for the emulation
results/emu_mhc_*            emulation output (summary TSV, figure)
run.sh                       one-command reproduction from hs38DH.fa
```

No human sequencing data are distributed. Reads are simulated from the reference (ground truth by
construction), or taken from a public whole-genome ONT BAM with `--reads-bam`.

## Decision rules (`code/as_emulate.py`)

MinKNOW's built-in AS algorithm, used for the clinical runs, is not public. readfish (Payne et al.
2021) is. We therefore never call one rule "the" AS rule:

| rule | keep the read when … | corresponds to |
|---|---|---|
| `primary` | the **primary** alignment overlaps a target interval | hypothesis for MinKNOW |
| `readfish` | `single_on` or `multi_on`: the end coordinate (`r_en` on +, `r_st` on −) of **any** alignment returned by mappy is in a target, as in readfish `make_decision()`; no MAPQ filter | readfish, documented TOML (`multi_on = "stop_receiving"`) |
| `mapqN` | primary in target **and** MAPQ ≥ N | stricter hypothesis |

`readfish` depends on how many secondary alignments mappy returns, so `--best-n` takes a list.

The predicted relative depth is `k + (1 − k) · L_ej / L_full`. Here k is the keep-rate, L_ej the length
of an ejected stub and L_full the full read length. Unmapped prefixes count as ejected.

## Key results (simulated, 1,000 molecules per locus, 400-b prefixes, 5 % error)

Keep-rate by capture reference and rule (`results/emu_mhc_summary.tsv`):

| locus | noalt | hs38DH primary | hs38DH readfish (best_n 5 / 20) | hs38DH mapq1 | ALT only, primary | HLA only, primary |
|---|---|---|---|---|---|---|
| CFB | 1.00 | **0.23** | 0.94 / 1.00 | 0.06 | 0.23 | 1.00 |
| C2–CFB block | 1.00 | 0.25 | 0.94 / 1.00 | 0.03 | 0.25 | 1.00 |
| C4A–C4B | 1.00 (mapq1: 0.24) | 0.34 | 0.96 / 1.00 | 0.06 | 0.34 | 1.00 |
| HLA-A | 1.00 | 0.38 | 0.76 / 0.94 | 0.28 | 0.91 | 0.39 |
| HLA-DRB1 | 1.00 | 0.12 | 0.75 / 1.00 | 0.01 | 1.00 | 0.13 |
| chr6 control | 1.00 | 1.00 | 0.99 | 1.00 | 1.00 | 1.00 |

What this shows:

1. **Class III (CFB, C2, C4) is depleted by the 7 MHC ALT haplotypes only.** The 525 HLA contigs
   have no effect there. HLA-DRB1 is depleted by the HLA contigs only, and HLA-A mostly by them.
2. **readfish's documented rule would not reproduce the depletion.** The read keeps a secondary
   alignment on primary chr6, which gives `multi_on`, so it is kept. With best_n = 20 the read is
   rescued completely.
3. **Only the `primary` rule fits both observations in the cohort:**
   - Under hs38DH, CFB reads are ejected: in the cohort 75 % of CFB reads are shorter than 1 kb,
     against 77 % of prefixes ejected in the emulation. This rules out readfish-like `multi_on`.
   - Under the alt-free reference, C4A–C4B reads are *not* ejected: median 6.1 kb, 4 % below 1 kb
     (Additional file 5), although 76 % of C4 prefixes have MAPQ 0 on the alt-free reference
     because C4A and C4B are paralogous. This rules out a MAPQ gate.
4. Predicted depth under `primary` is 0.31, against 0.22–0.23 observed. It is the same order of
   magnitude, not an exact match.

## Genome-wide extension

```bash
# loci = protein-coding genes overlapping an ALT/HLA/fix counterpart present in the capture reference
python3 code/make_alt_gene_bed.py --gtf gencode.v44.annotation.gtf.gz \
    --alt-file hs38DH.fa.alt [--locations hg38.p14.fixLocations.txt.gz] --fai hs38DH.fa.fai \
    [--genes clinical_genes.txt] --flank 2000 --out loci_alt_genes.bed
# capture references in full; target BED = the same genes
python3 code/as_emulate.py --ref noalt=GRCh38_no_alt.fa --ref hs38DH=hs38DH.fa [--ref hg38p14=hg38.p14.fa] \
    --target-bed loci_alt_genes.bed --loci-bed loci_alt_genes.bed --simulate 200 \
    --prefix-len 400 --best-n 5,20 --threads 16 --index-dir idx --out results/emu_genome
```

Whole-genome indexes need about 8–10 GB of RAM each. `--reads-bam` with a public ONT WGS
(e.g. HG002) replaces simulated reads with real haplotypes.

## "No patch": UCSC hg38.p14 (`code/run_p14.sh`)

```bash
bash code/run_p14.sh hg38.p14.fa.gz fixLocations.txt.gz altLocations.txt.gz gencode.v38.annotation.gtf.gz hs38DH.fa.fai
```

The run splits hg38.p14 into the primary assembly, the 261 GRCh38 ALT contigs, 166 GRC fix patches and
90 novel patches. A patch is any `*_alt` contig absent from hs38DH. It then emulates capture against
the primary assembly alone, the primary plus each class, and the full p14. The loci are the 1,673
protein-coding genes (GENCODE v38, ± 2 kb) that overlap a counterpart contig, each with 100
molecules. The target BED is the same genes in primary coordinates. Results:
`results/emu_p14_summary.tsv` and `results/p14_gene_keep.tsv` (one row per gene; reads simulated
from the primary sequence).

Number of genes with keep-rate < 0.5 under the primary-hit rule:

| capture reference | genes < 0.5 | genes < 0.8 |
|---|---|---|
| primary only | 1 | 19 |
| + 261 ALT | 341 | 721 |
| + 166 fix patches | 242 | 708 |
| + 90 novel patches | 63 | 234 |
| hg38.p14 (all) | 680 | 1,504 |

- **Fix patches deplete as much as ALT contigs do.** A fix patch is a near-identical copy of a
  primary region, so the primary hit falls on it about half the time: the median keep-rate is 0.53
  for the 666 genes under a fix patch only.
- **Genes on the study panel (CKD_SU_V6) that would be hit under hg38.p14**, with keep-rate under the
  patch class alone:
  - *ADAMTS13* (fix, 0.44);
  - *CFH* 0.64, *CFHR2* 0.66, *CFHR4* 0.54, *CFHR5* 0.60 and *F13B* 0.47 (novel patch over the RCA
    cluster);
  - *ACTN4* (0.50), *WDR73*, *PTEN*, *SDHB*, *NLRP3*;
  - in addition to the MHC.
- **The readfish rule (best_n 5) keeps almost every gene**: median keep-rate 1.00 under hg38.p14, and
  every panel gene listed above at 0.99-1.00. The 36 genes still below 0.8 lie in multi-copy regions
  (the chr19 LRC/KIR, the MHC, the chrY TSPY array), plausibly because more near-equal copies compete
  than the five alignments mappy reports (best_n 20 not run on p14).
- Reads simulated from the patch sequence (`@patch` loci) map patch coordinates by offset from the
  primary region. This is approximate (patch orientation and indels are ignored), so use those
  values for direction only.

### General or patch-specific? (`code/run_p14_patches.sh`, manuscript Fig. 5)

```bash
bash code/run_p14_patches.sh      # after run_p14.sh; about 25 min on 8 cores
```

Every ALT contig, fix patch and novel patch is aligned to its primary region (`code/patch_divergence.py`:
minimap2 asm20 through mappy, in 50-kb chunks) and cut into 10-kb windows of local divergence
(mismatches + indel events per kb). Each window is sampled with 20 molecules from the **primary
haplotype** and 20 from the **contig haplotype**. Capture is emulated against the primary assembly plus
the contig's own class; the target BED is the primary regions of all contigs (17,255 windows, 507
contigs). `code/fig5_patches.py` draws the figure and writes `results/fig5_patches_windows.tsv`,
`results/fig5_patches_contigs.tsv` and `results/fig5_patches_numbers.txt`.

The answer is **both**:

| differences per kb | windows | primary haplotype | contig haplotype | readfish, contig haplotype |
|---|---|---|---|---|
| 0 | 8,204 | 0.46 | 0.46 | 0.97 |
| 1–2 | 2,112 | 0.53 | 0.23 | 0.90 |
| 2–5 | 2,314 | 0.59 | 0.16 | 0.85 |
| > 20 | 254 | 0.59 | 0.03 | 0.23 |

- **General part.** Where the contig copies the primary sequence (48 % of windows, 63 % of fix-patch
  windows), the primary hit falls on either copy, and about half of the reads are kept whichever
  haplotype they come from.
- **Sequence-dependent part.** Where the sequences differ, reads from the contig haplotype go to the
  contig and are ejected (Spearman ρ = −0.69 against divergence). For a fix patch, the contig haplotype
  is the corrected sequence. Reads from the primary haplotype are kept more often (0.59 at 2–5 /kb), but
  not always.
- Per-contig medians (primary / contig haplotype): ALT 0.52 / 0.38, fix 0.51 / 0.48, novel 0.53 / 0.43.
- **Panel genes, contig haplotype:** *CFH* 0.33, *CFHR2* 0.40, *CFHR4* 0.34, *CFHR5* 0.44, *F13B* 0.44,
  *ADAMTS13* 0.51.
- **readfish does not rescue divergent contig haplotypes** (0.23 above 20 /kb): they have no on-target
  secondary alignment to fall back on.
- Against the primary assembly alone, contig-haplotype reads are kept at 0.95–0.99, or 0.83 for the most
  divergent fix windows. This is the control.
- Gene values in Fig. 5c average 10-kb windows and differ slightly from the gene ± 2 kb values of
  `p14_gene_keep.tsv`.
- One patch (chr22_ML143378v1_fix) crashed mappy with 50-kb chunks and was re-aligned with 10-kb chunks.
  Five contigs yielded no window of at least 2 kb.

## Caveats on the legacy T2T / CFHR3-CFHR1 demonstration (V13 Fig. 5)

The legacy run reports 0 % of CFHR3-CFHR1 reads at their "true locus" on a CFHR3-CFHR1-masked
reference. This is a **placement** metric, and it is 0 by construction because the locus does not
exist in that reference. The capture keep-rate is 0.91 against 0.99. The reads land on CFH, CFHR2 and
CFHR4, which are in the BED, so they are kept, sequenced in full and placed correctly once
re-aligned on GRCh38. A missing segment therefore costs little at capture when its paralogues are
targeted. Its cost is at alignment (Hamza et al. 2023), and that cost is reversible.

The T2T + patch driver targets every contig of the capture reference in full, so its keep-rate is
not informative. Its "~14 % ejected at MAPQ 0" depends on a MAPQ gate that readfish does not apply,
and that the C4 data argue MinKNOW does not apply either. The reference used is masked GRCh38, not
CHM13. To test a real T2T capture reference, use `as_emulate.py --ref chm13=… --ref-bed chm13=<BED
in CHM13 coordinates>` and read `keep_rate`, not `placement_ok`.

## Citation

- Hamza A, et al. Eur J Hum Genet. 2023;31:730-732. doi:10.1038/s41431-023-01350-8.
- Payne A, et al. Readfish. Nat Biotechnol. 2021;39:442-450. Li H. Minimap2. Bioinformatics 2018.

Functional Nephrogenomics, Sorbonne Université / AP-HP.
