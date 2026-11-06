#!/usr/bin/env bash
# Build the GMI1000 reference (chromosome AL646052.1 + megaplasmid AL646053.1) from NCBI
# and convert the annotation into the layout nf-core/rnaseq + RSEM/kallisto accept
# (every gene -> transcript -> exon; see scripts/clean_ncbi_gff.py).
#
# Output: inputs/GMI1000_chr_megaplasmid.fasta, inputs/GMI1000_chr_megaplasmid_cleaned.gff
set -euo pipefail
cd "$(dirname "$0")/.."

DL=ref_download
mkdir -p "$DL" inputs

for acc in AL646052.1 AL646053.1; do
    curl -sf "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$acc&rettype=fasta&retmode=text" -o "$DL/$acc.fasta"
    sleep 1
    curl -sf "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=$acc" -o "$DL/$acc.gff3"
    sleep 1
done

python3 scripts/clean_ncbi_gff.py "$DL/AL646052.1.gff3" > "$DL/chromosome_cleaned.gff"
# megaplasmid IDs that clash with chromosome IDs get a __<seqid>_<start>_<end>_<strand> suffix
python3 scripts/clean_ncbi_gff.py "$DL/AL646053.1.gff3" --existing-ids "$DL/chromosome_cleaned.gff" > "$DL/megaplasmid_cleaned.gff"

cat "$DL/AL646052.1.fasta" "$DL/AL646053.1.fasta" > inputs/GMI1000_chr_megaplasmid.fasta
{ cat "$DL/chromosome_cleaned.gff"; grep -v '^#' "$DL/megaplasmid_cleaned.gff"; } > inputs/GMI1000_chr_megaplasmid_cleaned.gff

echo "Sequences: $(grep -c '>' inputs/GMI1000_chr_megaplasmid.fasta)"
echo "Genes:     $(awk -F'\t' '$3=="gene"' inputs/GMI1000_chr_megaplasmid_cleaned.gff | wc -l)"
