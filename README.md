# Specific_gene_abundance_analysis_program
Specific ANAMMOX gene HzsA/B/C and Hdh abundance analysis program storage for paper
# How to use
Requirements: Windows 11 25H2 (other systems didn't test), Docker desktop, R, WSL
1. Prepare the reference genome using MAFFT & HMMER 
2. Put .hmm on ./hmm folder
3. Put .fastq.gz sequences on ./reads folder
4. Set CPU threads and sample IDs in run_pipeline.sh
5. Run run.ps1
