#!/usr/bin/env python3
"""Convert an NCBI GFF3 for GMI1000 into the 'cleaned' layout used for nf-core/rnaseq.

Layout (reverse-engineered from GMI1000_cleaned.gff):
  protein-coding: gene -> mRNA (transcript-X) -> exon (exon-transcript-X) + CDS (CDS-<locus>)
  tRNA / rRNA:    gene -> tRNA|rRNA (<type>-X) -> exon (exon-<type>-X)
  everything else (region, repeats, mobile elements, parentless misc_RNA, ...) is dropped.
Gene IDs whose name occurs more than once get a __<seqid>_<start>_<end>_<plus|minus> suffix.

usage: clean_ncbi_gff.py raw.gff3 [--existing-ids cleaned.gff] > out.gff
"""
import re
import sys
from collections import Counter, defaultdict
from urllib.parse import quote, unquote


def parse_attrs(s):
    out = {}
    for kv in s.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = unquote(v)
    return out


def enc(v):
    return quote(v.replace(";", ","), safe="',:+-._~*")


def fmt(attrs):
    return ";".join(f"{k}={enc(v)}" for k, v in attrs if v is not None)


def main():
    args = sys.argv[1:]
    existing = set()
    if "--existing-ids" in args:
        i = args.index("--existing-ids")
        for line in open(args[i + 1]):
            m = re.search(r"(?:^|\t|;)ID=([^;\n]+)", line)
            if m:
                existing.add(m.group(1))
        del args[i : i + 2]
    raw = args[0]

    genes, children = [], defaultdict(list)
    for line in open(raw):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        a = parse_attrs(f[8])
        if f[2] == "gene":
            genes.append((f, a))
        elif "Parent" in a:
            children[a["Parent"]].append((f, a))

    def base(gid):
        return re.sub(r"-\d+$", "", gid)

    counts = Counter(base(a["ID"]) for _, a in genes)
    out = []
    for f, a in genes:
        seqid, start, end, strand = f[0], f[3], f[4], f[6]
        gid = a["ID"]
        if counts[base(gid)] > 1 or gid in existing:
            gid = f"{base(gid)}__{seqid}_{start}_{end}_{'plus' if strand == '+' else 'minus'}"
        name = gid[len("gene-"):]
        col = lambda t, s, e: [seqid, "GenBank", t, s, e, ".", strand, "."]
        common = [("gene", a.get("gene")), ("locus_tag", a.get("locus_tag")), ("old_locus_tag", a.get("old_locus_tag"))]
        kids = children.get(a["ID"], [])
        cds = [k for k in kids if k[0][2] == "CDS"]
        rna = [k for k in kids if k[0][2] in ("tRNA", "rRNA")]
        if cds:
            tid = f"transcript-{name}"
            out.append(col("gene", start, end) + [fmt([("ID", gid)] + common)])
            out.append(col("mRNA", start, end) + [fmt([("ID", tid), ("Parent", gid)] + common)])
            out.append(col("exon", start, end) + [fmt([("ID", f"exon-{tid}"), ("Parent", tid)])])
            for cf, ca in cds:
                cid = f"CDS-{a.get('locus_tag') or name}"
                out.append(col("CDS", cf[3], cf[4]) + [fmt([("ID", cid)] + common + [("product", ca.get("product")), ("protein_id", ca.get("protein_id")), ("Parent", tid)])])
        elif rna:
            rf, ra = rna[0]
            rid = f"{rf[2]}-{name}"
            out.append(col("gene", start, end) + [fmt([("ID", gid)] + common)])
            out.append(col(rf[2], rf[3], rf[4]) + [fmt([("ID", rid)] + common + [("product", ra.get("product")), ("Parent", gid)])])
            out.append(col("exon", rf[3], rf[4]) + [fmt([("ID", f"exon-{rid}"), ("Parent", rid)])])
        else:
            print(f"WARN: dropped gene without CDS/RNA child: {a['ID']}", file=sys.stderr)

    print("##gff-version 3")
    for r in out:
        print("\t".join(r))


if __name__ == "__main__":
    main()
