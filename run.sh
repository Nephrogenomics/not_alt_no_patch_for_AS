#!/bin/bash
# Reproduce the two adaptive-sampling capture-reference demonstrations from a capture reference FASTA.
# Usage: bash run.sh /path/to/hs38DH.fa
# Needs: python3 with mappy, pysam, matplotlib. Writes figures/TSVs into results/.
set -euo pipefail
REF="${1:?usage: bash run.sh /path/to/hs38DH.fa}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="${WORK:-$HERE/work}"; mkdir -p "$WORK" "$HERE/results"
DEMO="$HERE/code/as_altmapping_demo.py"
[ -f "$REF.fai" ] || python3 -c "import pysam,sys; pysam.faidx(sys.argv[1])" "$REF"

echo "[refs] building reduced references (chr6 +/- ALT for MHC; chr1 intact vs CFHR3-CFHR1-masked for RCA)..."
python3 - "$REF" "$WORK" <<'PY'
import pysam, sys, re
ref, work = sys.argv[1], sys.argv[2]
fa = pysam.FastaFile(ref)
def write(path, pairs):
    with open(path, "w") as o:
        for name, seq in pairs:
            o.write(">"+name+"\n")
            for i in range(0, len(seq), 60): o.write(seq[i:i+60]+"\n")
    pysam.faidx(path)
    print("  wrote", path)
# MHC: chr6 + ALT/HLA/decoy  vs  chr6 only
alt = [c for c in fa.references if c=="chr6" or "_alt" in c.lower() or c.lower().startswith("hla-") or "decoy" in c.lower()]
write(f"{work}/mhc_ALT.fa",   [(c, fa.fetch(c)) for c in alt])
write(f"{work}/mhc_noalt.fa", [("chr6", fa.fetch("chr6"))])
# RCA: chr1 intact  vs  chr1 with CFHR3-CFHR1 masked to N (T2T-like deletion)
chr1 = list(fa.fetch("chr1")); a, b = 196745000, 196805000
for i in range(a, b): chr1[i] = "N"
write(f"{work}/chr1_RCApresent.fa", [("chr1", fa.fetch("chr1"))])
write(f"{work}/chr1_RCAdel.fa",     [("chr1", "".join(chr1))])
PY

printf 'chr1\t196598161\t197040250\n' > "$WORK/RCA_cluster.bed"
V6="$HERE/data/CKD_SU_V6_no_alt.bed"; V6ALT="$HERE/data/CKD_SU_V6_plus_alt.bed"
MHCREG=chr6:31887785-31962084; RCAREG=chr1:196745000-196805000

echo "[fig2] MHC / CFB - ALT-excess ejection"
python3 "$DEMO" --alt-ref "$WORK/mhc_ALT.fa" --noalt-ref "$WORK/mhc_noalt.fa" --target-bed "$V6" \
    --simulate 3000 --region $MHCREG --out "$HERE/results/fig2_mhc_alt"
echo "[fig2b] partial ALT-in-BED (fragile patch)"
python3 "$DEMO" --alt-ref "$WORK/mhc_ALT.fa" --noalt-ref "$WORK/mhc_noalt.fa" --target-bed "$V6ALT" \
    --simulate 3000 --region $MHCREG --mapq-min 1 --out "$HERE/results/fig2_mhc_altbed_mq1"

echo "[fig3] CFHR3-CFHR1 / T2T - segment-absence mis-mapping"
python3 "$DEMO" --alt-ref "$WORK/chr1_RCAdel.fa" --noalt-ref "$WORK/chr1_RCApresent.fa" --target-bed "$WORK/RCA_cluster.bed" \
    --simulate 3000 --region $RCAREG --out "$HERE/results/fig3_rca_cfhr_t2t"

echo "[done] figures + TSVs in $HERE/results/"
