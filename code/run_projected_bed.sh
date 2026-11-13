#!/bin/bash
# "Projected BED": would adding the contig coordinates of each target to the BED remove the loss?
# Run after code/run_p14.sh and code/run_p14_patches.sh (needs work/p14 and work/p14div).
#
# 1. code/project_bed.py lifts the primary BED onto every ALT / fix / novel contig through the 10-kb alignment
#    windows of patch_divergence.py (panel BED: all classes, and ALT only; Fig. 5 target BED: all classes);
# 2. panel run: the CKD panel genes that have a counterpart contig, sampled from the primary and the contig
#    haplotype, decided on the same alignments against the standard, clinical (+ALT HLA), projected-ALT and
#    projected BEDs (as_emulate.py --extra-bed);
# 3. windows run: the Fig. 5 windows, standard versus projected target BED;
# 4. code/summarize_projbed.py writes results/projbed_*.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
P="${P14:-$HERE/work/p14}"; D="${DIV:-$HERE/work/p14div}"; W="${WORK:-$HERE/work/projbed}"
mkdir -p "$W" "$HERE/results"; THREADS="${THREADS:-8}"
PB="$HERE/code/project_bed.py"; PANEL="$HERE/data/CKD_SU_V6_no_alt.bed"

# panel genes (GENCODE, +-2 kb) that overlap the panel and have a counterpart contig
python3 - "$PANEL" "$HERE/results/p14_genes.bed" "$W/genes_primary.bed" <<'EOF'
import collections, sys
panel = collections.defaultdict(list)
for l in open(sys.argv[1]):
    f = l.split("\t"); panel[f[0]].append((int(f[1]), int(f[2])))
with open(sys.argv[3], "w") as o:
    for l in open(sys.argv[2]):
        if l.startswith("#"): continue
        f = l.rstrip("\n").split("\t"); c, s, e = f[0], int(f[1]), int(f[2])
        if any(a < e and b > s for a, b in panel.get(c, [])):
            o.write(f"{c}\t{s}\t{e}\t{f[3]}|src=primary\n")
EOF
sed 's/|src=primary//' "$W/genes_primary.bed" > "$W/genes_plain.bed"
python3 "$PB" --bed "$W/genes_plain.bed" --windows "$D/div_windows.tsv" --loci --out "$W/genes_contig.bed"
cat "$W/genes_primary.bed" "$W/genes_contig.bed" > "$W/loci.bed"
python3 "$PB" --bed "$PANEL" --windows "$D/div_windows.tsv" --with-primary --out "$W/panel_projected.bed"
python3 "$PB" --bed "$PANEL" --windows "$D/div_windows.tsv" --classes alt --with-primary \
    --out "$W/panel_projected_altonly.bed"
python3 "$PB" --bed "$D/div_target.bed" --windows "$D/div_windows.tsv" --with-primary \
    --out "$W/div_target_projected.bed"
cp "$W/panel_projected.bed" "$HERE/results/CKD_SU_V6_projected.bed"

python3 "$HERE/code/as_emulate.py" \
    --ref primary="$P/p14_primary.fa" --ref plus_ALT="$P/p14_primary_plus_alt.fa" --ref hg38.p14="$P/p14_full.fa" \
    --source-ref "$P/p14_full.fa" --target-bed "$PANEL" \
    --extra-bed clinical="$HERE/data/CKD_SU_V6_plus_alt.bed" \
    --extra-bed projALT="$W/panel_projected_altonly.bed" --extra-bed proj="$W/panel_projected.bed" \
    --loci-bed "$W/loci.bed" --simulate 50 --prefix-len 400 --best-n 5 --threads "$THREADS" --no-figure \
    --out "$W/panel" > "$W/panel.log" 2>&1

python3 "$HERE/code/as_emulate.py" \
    --ref plus_ALT="$P/p14_primary_plus_alt.fa" --ref plus_fix="$P/p14_primary_plus_fix.fa" \
    --ref plus_novel="$P/p14_primary_plus_novel.fa" --ref hg38.p14="$P/p14_full.fa" \
    --source-ref "$P/p14_full.fa" --target-bed "$D/div_target.bed" --extra-bed proj="$W/div_target_projected.bed" \
    --loci-bed "$D/div_loci.bed" --simulate 20 --prefix-len 400 --best-n 5 --threads "$THREADS" --no-figure \
    --out "$W/win" > "$W/win.log" 2>&1

python3 "$HERE/code/summarize_projbed.py" --panel "$W/panel_summary.tsv" --windows "$W/win_summary.tsv" \
    --div "$D/div_windows.tsv" --out "$HERE/results/projbed"
