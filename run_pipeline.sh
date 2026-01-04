#!/usr/bin/env bash
set -euo pipefail

# Parameters
THREADS=24
SAMPLES=("A" "B") 

# 0. QC
mkdir -p qc
for S in "${SAMPLES[@]}"; do
  fastp \
    -i reads/${S}.fastq.gz \
    -o qc/${S}.fq.gz \
    -w ${THREADS}
done

# 1. Co-assembly
megahit \
  -r qc/A.fq.gz,qc/B.fq.gz \
  -o assembly \
  -t ${THREADS} \
  --min-contig-len 1000

CONTIGS=assembly/final.contigs.fa

# 2. ORF prediction
prodigal \
  -i ${CONTIGS} \
  -a proteins.faa \
  -p meta

# 3. HMM search (STRICT)
mkdir -p hmm_out
for G in hzsA hzsB hzsC hdh; do
  hmmsearch \
    --cpu ${THREADS} \
    -E 1e-20 --domE 1e-20 \
    hmm/${G}.hmm \
    proteins.faa \
    > hmm_out/${G}.out
done

# 4. filter results
mkdir -p filtered
for G in hzsA hzsB hzsC hdh; do
  grep -v "^#" hmm_out/${G}.out \
  | awk '$6 >= 300 {print $1}' \
  | sort -u \
  > filtered/${G}_proteins.txt
done

# 5. contig assembly
for G in hzsA hzsB hzsC hdh; do
  sed 's/_.*//' filtered/${G}_proteins.txt \
  | sort -u > filtered/${G}_contigs.txt
done

# 6. anammox contig filter
comm -12 filtered/hzsA_contigs.txt filtered/hzsB_contigs.txt \
| comm -12 - filtered/hzsC_contigs.txt \
| comm -12 - filtered/hdh_contigs.txt \
> results/anammox_contigs.txt

# 7. mapping
bowtie2-build ${CONTIGS} contigs

for S in "${SAMPLES[@]}"; do
  bowtie2 -x contigs \
    -U qc/${S}.fq.gz \
    --threads ${THREADS} \
  | samtools sort -@ 8 -o ${S}.bam
  samtools index ${S}.bam
done

# 8. quantification (TPM / RPKM)
coverm contig \
  --bam-files A.bam B.bam \
  --methods mean covered_fraction \
  --output-file results/coverm_all.tsv

awk 'NR==FNR {a[$1]; next} NR==1 || ($1 in a)' \
  results/anammox_contigs.txt \
  results/coverm_all.tsv \
  > results/anammox_coverm.tsv

# 9. plotting
Rscript plot_anammox_abundance.R

echo "Pipeline finished successfully."
