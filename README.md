# No ALT, no patch — and no missing target

Code and data to reproduce the in-silico demonstration behind the adaptive-sampling (AS) capture
reference study **"No alt, no patch: ALT contigs in the capture reference deplete MHC coverage under
nanopore adaptive sampling."**

Oxford Nanopore adaptive sampling decides, in real time, whether to keep or eject each nascent read
by aligning its first bases against a **capture reference** and keeping the molecule only when the
best alignment falls in a target-BED region. That decision is made **before** the read is fully
sequenced, so any read ejected here is physically lost and cannot be recovered by any downstream
re-alignment. The choice of capture reference therefore matters, in two symmetric ways:

| Failure mode | Locus | What happens to a target read | Signature |
|---|---|---|---|
| **Surplus copies** (ALT/HLA contigs present) | MHC class III (CFB, C2) | best hit lands on an off-target ALT/HLA contig -> **ejected** | coverage hole |
| **Missing segment** (target absent from the reference) | RCA / Factor-H cluster (CFHR3-CFHR1) | no true target -> **mis-maps onto paralogues** | locus destroyed, keep-rate stays high |

**The principle is symmetric: the AS capture reference must be a primary, alt-free assembly that
contains every target locus and its common structural alleles - neither surplus copies nor missing
segments. No ALT, no patch, and no missing target.**

## What this repo contains

```
code/as_altmapping_demo.py   emulates the readfish keep/eject decision, read by read
data/CKD_SU_V6_no_alt.bed    the study capture panel - primary-only (23 contigs, no ALT/HLA)
data/CKD_SU_V6_plus_alt.bed  the same panel + the ALT-contig intervals of an HLA-aware panel
results/fig2_mhc_alt.*       Figure 2 - MHC/CFB, ALT-excess ejection
results/fig3_rca_cfhr_t2t.*  Figure 3 - CFHR3-CFHR1, segment-absence (T2T-like) mis-mapping
run.sh                       one-command reproduction from a capture reference FASTA
```

No human sequencing data are distributed here: the demonstration uses **simulated** reads from the
target sequence (ground truth by construction) and public reference genomes. Real-read confirmation
for the CFHR3-CFHR1 case is possible from SRA **PRJNA916871** (Hamza et al. 2023).

## How the emulation works (`code/as_altmapping_demo.py`)

For each read prefix (default 400 bp, the order MinKNOW sees first) it runs `minimap2` (`map-ont`,
via `mappy`) against two references and applies the readfish keep-rule (best alignment overlaps a
target-BED interval):

- `--noalt-ref` : the alt-free primary assembly (the recommended reference)
- `--alt-ref`   : the ALT-bearing (hs38DH) or segment-deleted (T2T-like) reference under test
- `--mapq-min`  : 0 = vanilla readfish; >0 emulates a decision that rejects ambiguous (MAPQ-0) multi-mappers
- `--simulate N` / `--reads-bam` / `--reads-fastq` : read source (simulate draws from `--noalt-ref`)

Outputs: per-read TSV, a summary (keep-rate, best-hit destination, MAPQ), and a **locus-integrity**
metric for simulated reads (`--misplace-kb`) - the fraction of reads that reach their true locus,
which is the metric that matters at segmental-duplication loci where a broad BED masks mis-mapping
as apparent "keep".

## Reproduce

Dependencies: Python 3 with `mappy`, `pysam`, `matplotlib`, plus a capture reference FASTA
(`hs38DH.fa`, the GRCh38 full analysis set with decoy + HLA; e.g. from `bwa.kit`).

```bash
bash run.sh /path/to/hs38DH.fa
```

`run.sh` builds the reduced references (chr6 +/- ALT/HLA for the MHC case; chr1 intact vs
CFHR3-CFHR1-masked for the RCA case), runs every condition, and writes figures/TSVs into `results/`.

### Key results (simulated n = 3000, 400-bp prefixes)

**MHC / CFB (Figure 2)** - capture BED is primary-only (`CKD_SU_V6_no_alt.bed`):

| capture reference | keep-rate |
|---|---|
| alt-free (chr6) | 1.00 |
| ALT-bearing (chr6 + ALT/HLA) | **0.23** - 77% of prefixes' best hit lands on an off-target MHC ALT haplotype (chr6_GL000250-256v2_alt) |

Adding those ALT intervals to the BED does **not** rescue the reads (they scatter across all seven
MHC alt haplotypes); targeting every ALT contig in full restores 1.00 under a permissive rule but
collapses again (0.05) once a confident-MAPQ rule is applied - i.e. patching is fragile.

**CFHR3-CFHR1 / T2T (Figure 3)** - capture BED is the RCA/Factor-H cluster (chr1:196,598,161-197,040,250):

| capture reference | keep-rate | reads placed at the true CFHR3-CFHR1 locus |
|---|---|---|
| intact GRCh38 primary | 0.99 | **96%** |
| CFHR3-CFHR1-deleted (T2T-like) | 0.92 | **0%** |

On a T2T-like reference the keep-rate looks benign (0.92) but **0%** of the segment's reads reach
their true locus - they mis-map onto the CFHR paralogues, reproducing the uninterpretable "aberrant
patterns" of Hamza et al. 2023 at the capture stage.

## Citation

- Companion methods paper: *No alt, no patch: ALT contigs in the capture reference deplete MHC
  coverage under nanopore adaptive sampling* (Georget, Yousfi, Bensouna, Mesnard).
- Hamza A, El-Sissy C, Yousfi N, et al. *The absence of CFHR3 and CFHR1 genes from the T2T-CHM13
  assembly can limit the molecular diagnosis of complement-related diseases.* Eur J Hum Genet.
  2023;31(7):730-732. doi:10.1038/s41431-023-01350-8.
- Adaptive sampling: Payne A, et al. *Readfish enables targeted nanopore sequencing of gigabase-sized
  genomes.* Nat Biotechnol. 2021;39:442-450. Aligner: Li H. *Minimap2.* Bioinformatics 2018.

Functional Nephrogenomics, Sorbonne Universite / AP-HP.
