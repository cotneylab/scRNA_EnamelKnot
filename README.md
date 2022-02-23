# Dental Disease Gene Prioritization
Welcome to our github respository of analyses on the developing mouse tooth!

This respository contains the scripts used on ChIPseq data of E13.5 incisors, bulk RNAseq of E12.5-E17.5 molars and E12.5 incisors, scRNA-seq of E14 molars, and compiled human GWAS data for dental and craniofacial phenotypes and diseases.

For exploration of the data itself, please visit our laboratory website (https://cotney.research.uchc.edu/scrna-mouse-molar/).

To understand our processes, follow this general outline:
1. GREGOR_Caries.sh details our GWAS enrichment analyses using GREGOR. 
2. TSE_generation_motifs_ucscfiles.sh contains our generation of tooth-specific enhancers, as well as our criteria for identifying the transcription factor binding motifs enriched in these regions (using HOMER). It also generates the bed files used to format VISTA enhancers prior to the analysis in #2. 
3. TSE_VISTA_enrichment.Rscript includes the analysis used to identify tooth enhancers (and tooth-specific enhancers) which were previously valided by VISTA. 
4. TSE_motif_enrichment_figures.Rscript generates enrichment values for the motifs identifed in #1 using HOMER. 
5. TSE_GREAT_GO_DO_Enrichment.Rscript outlines how we leveraged GREAT and clusterProfiler to predict TSE target genes, and determine the enriched gene ontologies of these genes. note that this script filters TSE target genes to include only genes which we included in our annotation for subsequent bulk and single cell RNA-seq analyses. 
6. kallisto_kb_alignments.sh details our kallisto and kallisto-bustools alignments of publicly available bulk and single cell transcriptomic data obtained from GEO. 
7. WGCNA.Rscript contains detailed WGCNA scripts including all settings used, as well as our formatting of the eigengene trajectory over PC1 for all samples. 
8. kallisto_Seurat_singlecell.Rscript imports kallisto-bustools matrices into R and generates Seurat objects for each sample. 
9. Seurat_tooth_scRNA.Rscript details our single cell analysis of E14 molar samples published by Hallikas et al 2021, including identification of the putative enamel knot and its transcriptomic signature. 
10. GeneListPrioritization.Rscript demonstrates how we integrated all these analyses to prioritize genes based on their likelihood to contribute to dental disease. 
