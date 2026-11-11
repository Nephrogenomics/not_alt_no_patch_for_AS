#!/bin/bash
# Reproduce the adaptive-sampling capture-reference emulations from a capture reference FASTA.
# Usage: bash run.sh /path/to/hs38DH.fa [--legacy]
# Needs: python3 with mappy, pysam, matplotlib, numpy. Writes figures/TSVs into results/.
#   default : code/as_emulate.py — 4 capture references x 3 decision rules x best_n x prefix length
#   --legacy: also re-run the V13 figures (code/as_altmapping_demo.py; primary-hit rule only)
set -euo pipefail
REF="${1:?usage: bash run.sh /path/to/hs38DH.fa [--legacy]}"
LEGACY="${2:-}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="${WORK:-$HERE/work}"; mkdir -p "$WORK" "$HERE/results"
THREADS="${THREADS:-8}"
[ -f "$REF.fai" ] || python3 -c "import pysam,sys; pysam.faidx(sys.argv[1])" "$REF"

echo "[refs] building reduced references (chr6 alone, + ALT, + HLA, + both; decoys kept with ALT/HLA)..."
python3 - "$REF" "$WORK" <<'PY'
import pysam, sys, os
ref, work = sys.argv[1], sys.argv[2]
fa = pysam.FastaFile(ref)
def write(path, names):
    if os.path.exists(path + ".fai"): print("  exists", path); return
    with open(path, "w") as o:
        for c in names:
            s = fa.fetch(c); o.write(">" + c + "\n")
            for i in range(0, len(s), 60): o.write(s[i:i + 60] + "\n")
    pysam.faidx(path); print("  wrote", path, len(names), "contigs")
alt = [c for c in fa.references if c.endswith("_alt")]
hla = [c for c in fa.references if c.startswith("HLA-")]
dec = [c for c in fa.references if "decoy" in c.lower()]
write(f"{work}/mhc_noalt.fa",    ["chr6"])
write(f"{work}/mhc_ALT_HLA.fa",  ["chr6"] + alt + hla + dec)
write(f"{work}/mhc_ALT_only.fa", ["chr6"] + alt + dec)
write(f"{work}/mhc_HLA_only.fa", ["chr6"] + hla + dec)
# legacy RCA references (V13 Fig. 5)
chr1 = list(fa.fetch("chr1")); a, b = 196745000, 196805000
for i in range(a, b): chr1[i] = "N"
for name, seq in (("chr1_RCApresent", fa.fetch("chr1")), ("chr1_RCAdel", "".join(chr1))):
    p = f"{work}/{name}.fa"
    if not os.path.exists(p + ".fai"):
        with open(p, "w") as o:
            o.write(">chr1\n")
            for i in range(0, len(seq), 60): o.write(seq[i:i + 60] + "\n")
        pysam.faidx(p)
PY

V6="$HERE/data/CKD_SU_V6_no_alt.bed"
echo "[emu] MHC loci: 4 capture references x rules {primary, readfish, mapq1} x best_n {5,20} x prefix {200,400,800}"
python3 "$HERE/code/as_emulate.py" \
    --ref noalt="$WORK/mhc_noalt.fa" --ref hs38DH="$WORK/mhc_ALT_HLA.fa" \
    --ref ALT_only="$WORK/mhc_ALT_only.fa" --ref HLA_only="$WORK/mhc_HLA_only.fa" \
    --target-bed "$V6" --loci-bed "$HERE/data/loci_mhc.bed" --simulate 1000 \
    --prefix-len 200,400,800 --best-n 5,20 --threads "$THREADS" --index-dir "$WORK/idx" \
    --out "$HERE/results/emu_mhc"

if [ "$LEGACY" = "--legacy" ]; then
  DEMO="$HERE/code/as_altmapping_demo.py"
  echo "[legacy] V13 Fig. 2 (primary-hit rule only)"
  python3 "$DEMO" --alt-ref "$WORK/mhc_ALT_HLA.fa" --noalt-ref "$WORK/mhc_noalt.fa" --target-bed "$V6" \
      --simulate 3000 --region chr6:31887785-31962084 --out "$HERE/results/fig2_mhc_alt"
  printf 'chr1\t196598161\t197040250\n' > "$WORK/RCA_cluster.bed"
  echo "[legacy] V13 Fig. 5 (CFHR3-CFHR1 masked) — see README: placement metric, not a capture loss"
  python3 "$DEMO" --alt-ref "$WORK/chr1_RCAdel.fa" --noalt-ref "$WORK/chr1_RCApresent.fa" \
      --target-bed "$WORK/RCA_cluster.bed" --simulate 3000 --region chr1:196745000-196805000 \
      --out "$HERE/results/fig3_rca_cfhr_t2t"
fi
echo "[done] results in $HERE/results/"
