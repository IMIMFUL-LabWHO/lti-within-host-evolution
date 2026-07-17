#!/bin/bash
set -euo pipefail

# Load configuration
source config.sh

# Create necessary directories
mkdir -p trimmed qc mappings stats consensus variants results logs

# Function to log messages with timestamps
log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

# Function to generate sample list
generate_sample_list() {
    log "Generating sample list from input directory: ${INPUT_DIR}"

    # Check if input directory exists
    if [[ ! -d "${INPUT_DIR}" ]]; then
        log "Error: Input directory ${INPUT_DIR} not found."
        exit 1
    fi

    # Generate sample list
    if ls "${INPUT_DIR}"/*.fastq.gz &> /dev/null; then
        ls "${INPUT_DIR}" | cut -d "_" -f 1 | sort -V | uniq > samples
        log "Sample list generated and saved to 'samples' file."
    else
        log "Error: No FASTQ files found in ${INPUT_DIR}."
        exit 1
    fi

    # Check if samples file is empty
    if [[ ! -s samples ]]; then
        log "Error: No samples detected. Check input files and naming convention."
        exit 1
    fi
}

# Function to trim reads
trim_reads() {
    log "Trimming reads"
    for i in $(cat samples); do
        sample="${i}"
        log "Processing sample: ${sample}"

        fastp \
            --in1 "${INPUT_DIR}/${sample}_R1.fastq.gz" \
            --in2 "${INPUT_DIR}/${sample}_R2.fastq.gz" \
            --out1 trimmed/${sample}_trim_R1.fastq.gz \
            --out2 trimmed/${sample}_trim_R2.fastq.gz \
            --cut_front \
            --cut_tail \
            --cut_window_size 4 \
            --cut_mean_quality $qqp \
            --qualified_quality_phred $qqp \
            --length_required $lr \
            --correction \
            --detect_adapter_for_pe \
            --trim_poly_x \
            --trim_poly_g \
            --html qc/${sample}.fastp.html \
            --json qc/${sample}.fastp.json \
            --thread $thr
    done
    log "Read trimming completed."
}

# Function to align reads
align_reads() {
    log "Indexing reference sequence"
    bwa index reference/"${REFERENCE}".fasta
    samtools faidx reference/"${REFERENCE}".fasta
    samtools dict reference/"${REFERENCE}".fasta > reference/"${REFERENCE}".dict

    log "Mapping reads to reference"
    for i in $(cat samples); do
        sample="${i}"
        log "Aligning sample: ${sample}"

        bwa mem -t $thr -k 17 -R "@RG\tID:1\tSM:${sample}\tPL:ILLUMINA" reference/"${REFERENCE}".fasta trimmed/${sample}_trim_R1.fastq.gz trimmed/${sample}_trim_R2.fastq.gz | \
            samtools view -uhS -F4 -@$thr - | \
            samtools sort -@$thr -n - | \
            ivar trim -b reference/"${REFERENCE}".primer.bed -m 0 -q 0 -e | \
            samtools sort -@$thr -n - | \
            samtools fixmate -@$thr -m - mappings/${sample}_PE_fixmate.bam
        samtools sort -@$thr mappings/${sample}_PE_fixmate.bam | \
            samtools markdup -@$thr -s - mappings/${sample}_final.bam 
        samtools index mappings/${sample}_final.bam 
    done
    log "Alignment completed."
}

# Function to calculate mapping statistics
calculate_mapping_stats() {
    log "Calculating mapping statistics"

    for i in $(cat samples); do
        sample="${i}"
        log "Calculating mapping statistics for: ${sample}"

        echo -e "${sample}" > stats/${sample}_mapstats_name.log
        samtools flagstat mappings/${sample}_final.bam > stats/${sample}_allstats.log
        samtools coverage -q $qqp -Q $qqp mappings/${sample}_final.bam > stats/${sample}_coverage.log
        samtools depth -aa -H -q $qqp -Q $qqp mappings/${sample}_final.bam -o stats/${sample}.covdepth
        cat stats/${sample}_mapstats_name.log stats/${sample}_allstats.log stats/${sample}_coverage.log > stats/${sample}_stats.log
        rm  stats/${sample}_allstats.log stats/${sample}_mapstats_name.log stats/${sample}_coverage.log
    done

    log "Mapping statistics calculated."

}

# Function to generate mapping statistics report
generate_mapping_stats_report() {
    log "Generating mapping statistics report"

    echo -e "sample\trname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\tprimary_mapped\tr1_nreads\tr2_nreads\ttotal_reads" > results/mapstats.tsv

    for i in $(cat samples); do
        sample="${i}"
        log "Generating mapping statistics report for sample: ${sample}"

        first_line=$(head -n 1 stats/${sample}_stats.log)
        nineteenth_line=$(sed -n '19p' stats/${sample}_stats.log)
        primary_mapped=$(grep 'with itself and mate mapped' stats/${sample}_stats.log | awk '{print $1}')
        r1_nreads=$(gunzip -c trimmed/${sample}_trim_R1.fastq.gz | awk '{s++}END{print s/4}')
        r2_nreads=$(gunzip -c trimmed/${sample}_trim_R2.fastq.gz | awk '{s++}END{print s/4}')
        total_reads=$(($r1_nreads + $r2_nreads))
        echo -e "$first_line\t$nineteenth_line\t$primary_mapped\t$r1_nreads\t$r2_nreads\t$total_reads" >> results/mapstats.tsv
    done
    log "Mapping statistics generated."
    }

# Function to generate consensus sequences
generate_consensus() {
    log "Generating consensus sequence"

    for i in $(cat samples); do
        sample="${i}"
        log "Generating consensus for sample: ${sample}"

        samtools mpileup -aa -A -d 0 -E -f reference/"${REFERENCE}".fasta -q $qqp -Q $qqp mappings/${sample}_final.bam | \
            ivar consensus -p consensus/${sample} -t 0.5 -q $qqp -m $dc

        seqtk seq consensus/${sample}.fa | \
            tr "?RYSWKMBDHVN.ryswkmbdhvn" "N" | \
            tr "-" "N" | \
            tr [:lower:] [:upper:] | grep -iv ">" - > consensus/tmp.seq

        echo ">"$i > consensus/${sample}.fa
        cat consensus/${sample}.fa consensus/tmp.seq  > consensus/${sample}.fasta
        find ./consensus -type f ! -name '*.fasta' -delete
    done

    cat consensus/*.fasta > results/consensus_sequences.fasta
    log "Consensus sequences generated."
}

# Function to generate variant calls
call_variants() {
    log "Calling variants"

    for i in $(cat samples); do
        sample="${i}"
        log "Calling variants for sample: ${sample}"

        samtools mpileup -aa -A -d 0 -E -f reference/"${REFERENCE}".fasta -q $qqp -Q $qqp mappings/${sample}_final.bam | \
            ivar variants -p variants/${sample} -t 0.01 -q $qqp -m 10 -r reference/"${REFERENCE}".fasta -g reference/"${REFERENCE}".gff3
    done
    log "Variant calls generated."
}

# Function for data transformation
transform_data() {
    log "Transforming and consolidating results"

    # Add sample names to variant files
    for file in variants/*.tsv; do
        filename_noext=$(basename -- "$file")
        filename_noext="${filename_noext%.*}"
        awk -v filename="$filename_noext" 'BEGIN{OFS="\t"} {$1=filename; print}' "$file" > "${file}.tmp"
        mv "${file}.tmp" "$file"
    done

    # Add sample names to coverage files
    for file in stats/*.covdepth; do
        filename_noext=$(basename -- "$file")
        filename_noext="${filename_noext%.*}"
        awk -v filename="$filename_noext" 'BEGIN{OFS="\t"} {$1=filename; print}' "$file" > "${file}.tmp"
        mv "${file}.tmp" "$file"
    done

    # Merge variant annotation files
    awk 'FNR==1 && NR!=1{next;}{print}' variants/*.tsv > results/sleek_variants.tsv

    # Merge coverage files
    awk 'FNR==1 && NR!=1{next;}{print}' stats/*.covdepth > results/coverage.tsv

    # Convert files to CSV format
    tr '\t' ',' < results/coverage.tsv > results/coverage.csv
    tr '\t' ',' < results/sleek_variants.tsv > results/sleek_variants.csv
    tr '\t' ',' < results/mapstats.tsv > results/mapstats.csv

    # Cleanup intermediate files
    rm results/coverage.tsv
    rm results/sleek_variants.tsv
    rm results/mapstats.tsv
}

# Main pipeline
log "Starting pipeline"
generate_sample_list
trim_reads
align_reads
calculate_mapping_stats
generate_mapping_stats_report
generate_consensus
call_variants
transform_data
log "Pipeline completed successfully"