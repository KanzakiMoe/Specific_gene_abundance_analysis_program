# Specific_gene_abundance_analysis_program
Specific ANAMMOX gene HzsA/B/C and Hdh abundance analysis program storage for paper

# How to use
Requirements: Windows 11 25H2 (other systems didn't test), Docker desktop, R, WSL
1. Prepare the reference genome using MAFFT & HMMER 
2. Put .hmm on ./hmm folder
3. Put .fastq.gz sequences on ./reads folder
4. Set CPU threads and sample IDs in run_pipeline.sh
5. Run run.ps1

# References
1. Aroney, S. T., Newell, R. J., Nissen, J. N., Camargo, A. P., Tyson, G. W., & Woodcroft, B. J. (2025). CoverM: read alignment statistics for metagenomics. Bioinformatics, 41(4), btaf147.
2. Chen, S. (2025). Fastp 1.0: An ultra‐fast all‐round tool for FASTQ data quality control and preprocessing. Imeta, 4(5), e70078.
3. Hmmer.org. (2023). HMMER: biosequence analysis using profile hidden Markov models. HMMER 3.4. Web. http://hmmer.org/.
4. Hyatt, D., Chen, G. L., LoCascio, P. F., Land, M. L., Larimer, F. W., & Hauser, L. J. (2010). Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC bioinformatics, 11(1), 119.
5. Katoh, K., Rozewicki, J., & Yamada, K. D. (2019). MAFFT online service: multiple sequence alignment, interactive sequence choice and visualization. Briefings in bioinformatics, 20(4), 1160-1166.
6. Kuraku, S., Zmasek, C. M., Nishimura, O., & Katoh, K. (2013). aLeaves facilitates on-demand exploration of metazoan gene family trees on MAFFT sequence alignment server with enhanced interactivity. Nucleic acids research, 41(W1), W22-W28.
7. Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), 357-359.
8. Li, D., Liu, C-M., Luo, R., Sadakane, K., and Lam, T-W., (2015) MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics, doi: 10.1093/bioinformatics/btv033 [PMID: 25609793].
9. Li, D., Luo, R., Liu, C.M., Leung, C.M., Ting, H.F., Sadakane, K., Yamashita, H. and Lam, T.W., 2016. MEGAHIT v1.0: A Fast and Scalable Metagenome Assembler driven by Advanced Methodologies and Community Practices. Methods.
10. Merkel, D. (2014). Docker: lightweight linux containers for consistent development and deployment. Linux j, 239(2), 2.
11. Wickham, H. (2016). Data analysis. In ggplot2: elegant graphics for data analysis (pp. 189-201). Cham: Springer international publishing.
