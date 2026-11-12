#!/bin/bash
# Per-patch "no patch" emulation (manuscript Fig. 5). Run after code/run_p14.sh (needs work/p14).
#
# 1. every ALT / fix / novel contig of hg38.p14 is aligned to its primary region and cut into 10-kb windows
#    of local divergence (code/patch_divergence.py);
# 2. each window is sampled twice, from the primary haplotype and from the contig haplotype, and the
#    capture decision is emulated against the primary assembly alone, the primary assembly plus each
#    contig class, and the full hg38.p14 (target BED = the primary regions of all contigs);
# 3. code/fig5_patches.py draws Fig. 5 and writes the per-window and per-contig tables.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
P14="${P14:-$HERE/work/p14}"; WORK="${WORK:-$HERE/work/p14div}"; mkdir -p "$WORK" "$HERE/results"
THREADS="${THREADS:-8}"; N="${N:-20}"

python3 "$HERE/code/patch_divergence.py" --fa "$P14/p14_full.fa" \
    --locations "$P14/fixLocations.txt.gz" --locations "$P14/altLocations.txt.gz" \
    --orig-alt "$P14/orig_alt.txt" --threads "$THREADS" --out "$WORK/div"

python3 "$HERE/code/as_emulate.py" \
    --ref primary="$P14/p14_primary.fa" --ref plus_ALT="$P14/p14_primary_plus_alt.fa" \
    --ref plus_fix="$P14/p14_primary_plus_fix.fa" --ref plus_novel="$P14/p14_primary_plus_novel.fa" \
    --ref hg38.p14="$P14/p14_full.fa" \
    --source-ref "$P14/p14_full.fa" --target-bed "$WORK/div_target.bed" --loci-bed "$WORK/div_loci.bed" \
    --simulate "$N" --prefix-len 400 --best-n 5 --threads "$THREADS" --no-figure --out "$WORK/emu"

python3 "$HERE/code/fig5_patches.py" --windows "$WORK/div_windows.tsv" --segments "$WORK/div_segments.tsv" \
    --emu "$WORK/emu_summary.tsv" --genes "$HERE/results/p14_gene_keep.tsv" --out "$HERE/results/fig5_patches"
