#!/usr/bin/env python3
"""Build the loci BED for the genome-wide emulation: genes of the primary assembly that have an
ALT, HLA or patch counterpart in a capture reference.

Inputs
  --gtf          GENCODE / Ensembl GTF (.gz accepted), primary-assembly coordinates, 'chr' names
  --alt-file     bwakit hs38DH.fa.alt (SAM records of each ALT/HLA contig aligned to the primary)
  --locations    UCSC altLocations.txt(.gz) / fixLocations.txt(.gz), or any BED-like
                 'chrom start end name' table of primary regions with a counterpart contig; repeatable
  --fai          .fai of the capture reference: keep only counterpart contigs it actually contains
  --genes        optional list of clinically relevant gene symbols (one per line, first column)
                 -> column 'in_list'; with --only-list, output is restricted to them
Output (stdout or --out): chrom start end gene n_counterparts counterpart_classes counterparts in_list
Use it as --loci-bed (and, with --flank, as --target-bed) of as_emulate.py.
"""
import argparse, gzip, re, sys
from collections import defaultdict

CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def contig_class(c):
    c = c.lower()
    return ("hla" if c.startswith("hla-") else "alt" if c.endswith("_alt") else
            "fix" if c.endswith("_fix") else "decoy" if "_decoy" in c else "other")


def regions_from_alt(path):
    out = []
    with opener(path) as fh:
        for ln in fh:
            if ln.startswith("@"):
                continue
            f = ln.split("\t")
            if len(f) < 6 or f[2] == "*" or f[5] == "*":
                continue
            span = sum(int(n) for n, op in CIG.findall(f[5]) if op in "MDN=X")
            st = int(f[3]) - 1
            out.append((f[2], st, st + span, f[0]))
    return out


def regions_from_locations(path):
    out = []
    with opener(path) as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) >= 5 and f[0].isdigit():      # UCSC table with leading 'bin' column
                f = f[1:]
            out.append((f[0], int(f[1]), int(f[2]), (f[3] if len(f) > 3 else "?").split(":")[0]))
    return out


def genes_from_gtf(path, biotypes):
    out = []
    with opener(path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            a = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
            bt = a.get("gene_type", a.get("gene_biotype", ""))
            if biotypes and bt not in biotypes:
                continue
            chrom = f[0] if f[0].startswith("chr") else "chr" + f[0]
            out.append((chrom, int(f[3]) - 1, int(f[4]), a.get("gene_name", a.get("gene_id", "?"))))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--alt-file", action="append", default=[])
    ap.add_argument("--locations", action="append", default=[])
    ap.add_argument("--fai")
    ap.add_argument("--genes")
    ap.add_argument("--only-list", action="store_true")
    ap.add_argument("--biotype", default="protein_coding", help="comma list; '' = all")
    ap.add_argument("--flank", type=int, default=0)
    ap.add_argument("--out")
    args = ap.parse_args()

    regs = []
    for p in args.alt_file:
        regs += regions_from_alt(p)
    for p in args.locations:
        regs += regions_from_locations(p)
    if args.fai:
        present = {ln.split("\t")[0] for ln in open(args.fai)}
        regs = [r for r in regs if r[3] in present]
    if not regs:
        sys.exit("no counterpart regions (check --alt-file/--locations/--fai)")
    by = defaultdict(list)
    for c, a, b, n in regs:
        by[c].append((a, b, n))
    wanted = None
    if args.genes:
        wanted = {ln.split()[0].strip() for ln in open(args.genes) if ln.strip() and not ln.startswith("#")}
    bts = {x for x in args.biotype.split(",") if x}
    o = open(args.out, "w") if args.out else sys.stdout
    o.write("#chrom\tstart\tend\tgene\tn_counterparts\tcounterpart_classes\tcounterparts\tin_list\n")
    n = 0
    for c, a, b, g in genes_from_gtf(args.gtf, bts):
        hits = sorted({nm for (x, y, nm) in by.get(c, ()) if x < b and y > a})
        inl = "NA" if wanted is None else ("yes" if g in wanted else "no")
        if not hits or (args.only_list and inl != "yes"):
            continue
        cls = ",".join(sorted({contig_class(h) for h in hits}))
        o.write(f"{c}\t{max(0, a - args.flank)}\t{b + args.flank}\t{g}\t{len(hits)}\t{cls}\t"
                f"{','.join(hits) if len(hits) <= 20 else ','.join(hits[:20]) + ',...'}\t{inl}\n")
        n += 1
    print(f"{n} genes with a counterpart contig ({len(regs)} counterpart regions)", file=sys.stderr)


if __name__ == "__main__":
    main()
