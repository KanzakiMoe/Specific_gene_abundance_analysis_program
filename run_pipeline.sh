#!/usr/bin/env bash
set -euo pipefail
set -x

# Parameters
THREADS=24
SAMPLES=("A" "B" "C")

# Cache root
CACHE=cache

QC=${CACHE}/qc
ASSEMBLY=${CACHE}/assembly
HMM=${CACHE}/hmm
FILTERED=${CACHE}/filtered
MAPPING=${CACHE}/mapping
LOGS=${CACHE}/logs

mkdir -p ${QC} ${ASSEMBLY} ${HMM} ${FILTERED} ${MAPPING} ${LOGS} results

# 0. QC
for S in "${SAMPLES[@]}"; do
  fastp \
    -i reads/${S}.fastq.gz \
    -o ${QC}/${S}.fq.gz \
    -w ${THREADS} \
    > ${LOGS}/fastp_${S}.log 2>&1
done

# 1. Co-assembly
rm -rf ${ASSEMBLY}

megahit \
  -r ${QC}/A.fq.gz,${QC}/B.fq.gz \
  -o ${ASSEMBLY} \
  -t ${THREADS} \
  --min-contig-len 1000 \
  > ${LOGS}/megahit.log 2>&1

CONTIGS=${ASSEMBLY}/final.contigs.fa

# 2. ORF prediction
prodigal \
  -i ${CONTIGS} \
  -a proteins.faa \
  -p meta

# 3. HMM search (STRICT)
for G in hzsA hzsB hzsC hdh; do
  hmmsearch \
    --cpu ${THREADS} \
    -E 1e-1 --domE 1e-1 \
    --tblout ${HMM}/${G}.tbl \
    hmm/${G}.hmm \
    proteins.faa \
    > ${HMM}/${G}.log 2>&1
done

# 4. filter HMM hits by score
for G in hzsA hzsB hzsC hdh; do
  awk '$1 !~ /^#/ && $6 >= 50 {print $1}' \
    ${HMM}/${G}.tbl \
    | sort -u > ${FILTERED}/${G}_proteins.txt
done

# 5. contig assembly
for G in hzsA hzsB hzsC hdh; do
  sed 's/_.*//' ${FILTERED}/${G}_proteins.txt \
    | sort -u > ${FILTERED}/${G}_contigs.txt
done

# 6. anammox contig filter
cat \
  ${FILTERED}/hzsA_contigs.txt \
  ${FILTERED}/hzsB_contigs.txt \
  ${FILTERED}/hzsC_contigs.txt \
  ${FILTERED}/hdh_contigs.txt \
| sort | uniq -c | awk '$1 >= 1 {print $2}' \
> results/anammox_contigs.txt

# 7. mapping
bowtie2-build ${CONTIGS} ${MAPPING}/contigs

for S in "${SAMPLES[@]}"; do
  bowtie2 -x ${MAPPING}/contigs \
    -U ${QC}/${S}.fq.gz \
    --threads ${THREADS} \
  | samtools sort -@ 8 -o ${MAPPING}/${S}.bam

  samtools index ${MAPPING}/${S}.bam
done > ${LOGS}/mapping.log 2>&1

# 8. quantification
# temp fix bam naming issue
coverm contig \
  --bam-files ${MAPPING}/A.bam ${MAPPING}/B.bam ${MAPPING}/C.bam \
  --methods mean covered_fraction \
  --output-file results/coverm_all.tsv \
  --threads ${THREADS}

# Generate mapping summary
echo -e "sample\tmapped\ttotal" > results/bam_mapped_summary.tsv
for sample in "${SAMPLES[@]}"; do
    total=$(samtools view -c "cache/mapping/${sample}.bam")
    mapped=$(samtools view -c -F 4 "cache/mapping/${sample}.bam")
    echo -e "${sample}\t${mapped}\t${total}" >> results/bam_mapped_summary.tsv
done


awk 'NR==FNR {a[$1]; next} NR==1 || ($1 in a)' \
  results/anammox_contigs.txt \
  results/coverm_all.tsv \
  > results/anammox_coverm.tsv

# 9. plotting
Rscript plot_anammox_abundance.R

echo "Pipeline finished successfully."
