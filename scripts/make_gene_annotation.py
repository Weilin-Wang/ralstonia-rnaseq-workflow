#!/usr/bin/env python3
"""Gene annotation table (locus tag, gene name, product, replicon, coordinates) from the cleaned GFF.

usage: make_gene_annotation.py cleaned.gff out.tsv
"""
import sys
from urllib.parse import unquote

REPLICONS = {"AL646052.1": "chromosome", "AL646053.1": "megaplasmid"}


def records(path):
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        attrs = {k: unquote(v) for k, v in (kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)}
        yield f, attrs


def main():
    gff, out = sys.argv[1:3]
    genes, tx2gene, products = {}, {}, {}
    for f, a in records(gff):
        if f[2] == "gene":
            genes[a["ID"]] = dict(gene_id=a["ID"], locus_tag=a.get("locus_tag", ""), gene_name=a.get("gene", ""),
                                  old_locus_tag=a.get("old_locus_tag", ""), replicon=REPLICONS.get(f[0], f[0]),
                                  start=f[3], end=f[4], strand=f[6], product="")
        elif f[2] in ("mRNA", "tRNA", "rRNA"):
            tx2gene[a["ID"]] = a["Parent"]
        if f[2] in ("CDS", "tRNA", "rRNA"):
            products[a["Parent"]] = a.get("product", "")
    for parent, prod in products.items():
        gid = tx2gene.get(parent, parent)
        if gid in genes:
            genes[gid]["product"] = prod

    cols = ["gene_id", "locus_tag", "gene_name", "old_locus_tag", "replicon", "start", "end", "strand", "product"]
    with open(out, "w") as o:
        o.write("\t".join(cols) + "\n")
        for g in genes.values():
            o.write("\t".join(g[c].replace("\t", " ") for c in cols) + "\n")
    print(f"{len(genes)} genes written to {out}")


if __name__ == "__main__":
    main()
