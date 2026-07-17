#!/usr/bin/env bash

# Script to download FASTQ files from ENA
# Integrated longitudinal analysis of SARS-CoV-2 within-host evolution during prolonged infection
# Alen Suljič (alen.suljic@mf.uni-lj.si)

set -euo pipefail

CSV="${1:-../data/samples.csv}"
OUTDIR="${2:-fastq_data}"

mkdir -p "$OUTDIR"

tail -n +2 "$CSV" |
while IFS=',' read -r sample accession; do

    sample="${sample//$'\r'/}"
    accession="${accession//$'\r'/}"

    echo "Downloading ${sample} (${accession})"

    urls=$(
    curl -fsSL \
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${accession}&result=read_run&fields=fastq_ftp&format=tsv" |
    tail -n +2 |
    cut -f2 |
    tr -d '\r'
    )
    
    IFS=';' read -ra fastq_urls <<< "$urls"

    for url in "${fastq_urls[@]}"; do
        filename="${url##*/}"

        case "$filename" in
            *_1.fastq.gz)
                output="${OUTDIR}/${sample}_R1.fastq.gz"
                ;;
            *_2.fastq.gz)
                output="${OUTDIR}/${sample}_R2.fastq.gz"
                ;;
            *)
                output="${OUTDIR}/${sample}.fastq.gz"
                ;;
        esac

        wget -c "https://${url}" -O "$output"
    done

done