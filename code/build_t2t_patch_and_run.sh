#!/bin/bash
# T2T + patch scenario: can patching CFHR3-CFHR1 back into a T2T capture reference (with a BED
# targeting it) rescue capture of a CARRIER's CFHR3-CFHR1 reads? T2T-CHM13 lacks these genes
# (Hamza et al. 2023), so this tests the pragmatic mitigation of adding the segment as a patch contig.
#
# Usage: bash build_t2t_patch_and_run.sh /path/hs38DH.fa /path/chm13_rca.fa
#   arg1 = GRCh38 primary FASTA (source of the CFHR3-CFHR1 patch and of the simulated carrier reads)
#   arg2 = T2T-CHM13 RCA/Factor-H region FASTA (lacks CFHR3-CFHR1), or T2T chr1
# Needs python3 + mappy + pysam. Writes results/t2tpatch_mq{0,1}_perread.tsv and a breakdown.
set -euo pipefail
HG38="${1:?usage: build_t2t_patch_and_run.sh hs38DH.fa chm13_rca.fa}"
T2T="${2:?usage: build_t2t_patch_and_run.sh hs38DH.fa chm13_rca.fa}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="${WORK:-$HERE/../work_t2tpatch}"; mkdir -p "$WORK" "$HERE/../results"
DEMO="$HERE/as_altmapping_demo.py"
[ -f "$HG38.fai" ] || python3 -c "import pysam,sys; pysam.faidx(sys.argv[1])" "$HG38"

python3 - "$HG38" "$T2T" "$WORK" <<'PY'
import pysam, sys
hg38, t2t, work = sys.argv[1], sys.argv[2], sys.argv[3]
patch = pysam.FastaFile(hg38).fetch("chr1", 196745000, 196805000)   # hg38 CFHR3-CFHR1 = the patch
tr = pysam.FastaFile(t2t)
cap = f"{work}/t2t_plus_patch.fa"
with open(cap, "w") as o:
    for c in tr.references:
        s = tr.fetch(c); o.write(">T2T_"+c+"\n")
        for i in range(0, len(s), 60): o.write(s[i:i+60]+"\n")
    o.write(">CFHR3_CFHR1_patch\n")
    for i in range(0, len(patch), 60): o.write(patch[i:i+60]+"\n")
pysam.faidx(cap)
ff = pysam.FastaFile(cap)
with open(f"{work}/t2t_plus_patch.bed", "w") as b:
    for c, l in zip(ff.references, ff.lengths): b.write(f"{c}\t0\t{l}\n")   # target the whole capture ref
c1 = f"{work}/hg38_chr1.fa"                                                 # read source (carrier)
with open(c1, "w") as o:
    s = pysam.FastaFile(hg38).fetch("chr1"); o.write(">chr1\n")
    for i in range(0, len(s), 60): o.write(s[i:i+60]+"\n")
pysam.faidx(c1)
print("built capture reference (T2T RCA + CFHR3-CFHR1 patch):", ff.references)
PY

CAP="$WORK/t2t_plus_patch.fa"; BED="$WORK/t2t_plus_patch.bed"; SRC="$WORK/hg38_chr1.fa"
for MQ in 0 1; do
  python3 "$DEMO" --alt-ref "$CAP" --noalt-ref "$SRC" --target-bed "$BED" \
    --simulate 3000 --region chr1:196745000-196805000 --mapq-min "$MQ" --out "$HERE/../results/t2tpatch_mq$MQ"
done
echo ""
echo "=== Where do CFHR3-CFHR1 reads land on the T2T+patch capture reference? ==="
for f in t2tpatch_mq0 t2tpatch_mq1; do
  echo "--- $f (readfish keep-rule; mq0 = permissive, mq1 = confident-mapping) ---"
  awk -F'\t' 'NR>1{tot++; c=$7; m=$9;
    if(c=="CFHR3_CFHR1_patch"){patch++; if(m==0)p0++}
    else if($6=="unmapped"){un++} else {par++; if(m==0)pa0++}}
    END{printf "  patch (correct): %d/%d (%.0f%%) MAPQ0=%d | T2T paralogue: %d (%.0f%%) MAPQ0=%d | unmapped=%d\n",
        patch,tot,100*patch/tot,p0,par,100*par/tot,pa0,un}' "$HERE/../results/${f}_perread.tsv"
done
echo ""
echo "NOTE: ignore the tool's [locus-integrity] line for this run (cross-contig coordinate artifact);"
echo "      the patch-vs-paralogue breakdown above is the correct metric."
