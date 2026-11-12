#!/bin/bash
# "No patch": emulate adaptive sampling against UCSC hg38.p14 (primary + 261 GRCh38 ALT + GRC fix and
# novel patches) and against the same assembly with each contig class added alone.
#
# Usage: bash code/run_p14.sh <hg38.p14.fa.gz> <fixLocations.txt.gz> <altLocations.txt.gz> <gencode.gtf.gz> <hs38DH.fa.fai>
#   the hs38DH .fai lists the 261 original GRCh38 ALT contigs; any other *_alt contig of p14 is a novel patch
# Outputs: results/p14_genes.bed, results/emu_p14_summary.tsv (+ work/p14/*.fa reduced references)
#
# Loci = protein-coding genes (+/- 2 kb) of the primary assembly overlapped by an ALT, fix or novel patch
# contig, plus non-patched panel controls. Each gene is sampled from the primary sequence; genes under a
# fix or novel patch are also sampled from the patch sequence ("@patch"), i.e. from an individual who
# carries the patched (for fix patches: corrected) sequence. Target BED = the gene loci on the primary.
set -euo pipefail
FA="${1:?hg38.p14.fa.gz}"; FIX="${2:?fixLocations}"; ALT="${3:?altLocations}"; GTF="${4:?gtf}"; HSFAI="${5:?hs38DH.fa.fai}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="${WORK:-$HERE/work/p14}"; mkdir -p "$WORK" "$HERE/results"
THREADS="${THREADS:-8}"; N="${N:-100}"

echo "[refs] splitting hg38.p14 by contig class"
cut -f1 "$HSFAI" | grep '_alt$' > "$WORK/orig_alt.txt"
python3 - "$FA" "$WORK" <<'PY'
import gzip, sys, os, re
fa, work = sys.argv[1], sys.argv[2]
ORIG = set(open(f"{work}/orig_alt.txt").read().split())
def cls(n):
    if n.endswith("_fix"): return "fix"
    if n.endswith("_alt"):
        return "alt" if n in ORIG else "novel"
    return "primary"
outs = {k: open(f"{work}/p14_{k}.fa", "w") for k in ("primary", "alt", "fix", "novel")}
full = open(f"{work}/p14_full.fa", "w")
counts = {}
w = None
with gzip.open(fa, "rt") as fh:
    for ln in fh:
        if ln.startswith(">"):
            name = ln[1:].split()[0]; c = cls(name); counts[c] = counts.get(c, 0) + 1
            w = [full, outs[c]]
        for o in w: o.write(ln)
for o in list(outs.values()) + [full]: o.close()
print("  contigs:", counts)
PY
# reference = primary + one class (the class file holds only its own contigs)
for c in alt fix novel; do cat "$WORK/p14_primary.fa" "$WORK/p14_$c.fa" > "$WORK/p14_primary_plus_$c.fa"; rm "$WORK/p14_$c.fa"; done
python3 -c "import pysam,sys; [pysam.faidx(f) for f in sys.argv[1:]]" "$WORK/p14_full.fa" "$WORK/p14_primary.fa"

echo "[loci] genes overlapped by ALT / fix / novel contigs"
python3 "$HERE/code/make_alt_gene_bed.py" --gtf "$GTF" --locations "$FIX" --locations "$ALT" \
    --fai "$WORK/p14_full.fa.fai" --flank 2000 --out "$HERE/results/p14_genes.bed"
python3 - "$HERE/results/p14_genes.bed" "$WORK/p14_full.fa.fai" "$FIX" "$ALT" "$WORK" <<'PY'
import sys, gzip
genes, fai, fixf, altf, work = sys.argv[1:]
L = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai)}
loc = {}
for f in (fixf, altf):
    for ln in gzip.open(f, "rt"):
        x = ln.rstrip("\n").split("\t")
        if x[0].isdigit(): x = x[1:]
        if x[3] in L: loc.setdefault(x[3], (x[0], int(x[1]), int(x[2])))
ORIG = set(open(f"{work}/orig_alt.txt").read().split())
def cls(n):
    if n.endswith("_fix"): return "fix"
    return "alt" if n in ORIG else "novel"
tgt = open(f"{work}/target.bed", "w"); loci = open(f"{work}/loci.bed", "w")
seen = set()
for ln in open(genes):
    if ln.startswith("#"): continue
    c, a, b, g, n, classes, cps, _ = ln.rstrip("\n").split("\t")
    a, b = int(a), int(b)
    if (c, g) in seen: continue
    seen.add((c, g))
    cl = sorted({cls(x) for x in cps.split(",") if x and x != "..."})
    tag = "+".join(cl)
    tgt.write(f"{c}\t{a}\t{b}\t{g}\n")
    loci.write(f"{c}\t{a}\t{b}\t{g}|{tag}|primary\n")
    for x in cps.split(","):                      # first fix, else first novel patch: sample from it
        if x in loc and cls(x) in ("fix", "novel"):
            pc, ps, pe = loc[x]
            s, e = max(0, a - ps), min(L[x], b - ps)
            if e - s > 1000:
                loci.write(f"{x}\t{s}\t{e}\t{g}|{tag}|patch\n")
            break
for c, a, b, g in (("chr1", 207752038, 207795516, "CD46"), ("chr19", 6677704, 6730562, "C3"),
                   ("chr6", 31946247, 31952084, "CFB")):
    tgt.write(f"{c}\t{a}\t{b}\t{g}\n"); loci.write(f"{c}\t{a}\t{b}\t{g}|control|primary\n")
tgt.close(); loci.close()
print("  loci:", sum(1 for _ in open(f"{work}/loci.bed")))
PY

echo "[emu] 5 capture references, rules primary / readfish / mapq1, prefix 400"
python3 "$HERE/code/as_emulate.py" \
    --ref primary="$WORK/p14_primary.fa" --ref plus_ALT="$WORK/p14_primary_plus_alt.fa" \
    --ref plus_fix="$WORK/p14_primary_plus_fix.fa" --ref plus_novel="$WORK/p14_primary_plus_novel.fa" \
    --ref hg38.p14="$WORK/p14_full.fa" \
    --source-ref "$WORK/p14_full.fa" --target-bed "$WORK/target.bed" --loci-bed "$WORK/loci.bed" \
    --simulate "$N" --prefix-len 400 --best-n 5 --threads "$THREADS" --no-figure \
    --out "$HERE/results/emu_p14"
echo "[done] results/emu_p14_summary.tsv"
