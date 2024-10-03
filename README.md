Comparative Analysis of High Versus No Mutations in the TP53 Gene in Prostate Cancer (PRAD) Using the TCGA Pipeline

Raj Habib

Thesis Statement: The primary objective of this project is to compare the high versus no mutations in the TP53 gene in prostate cancer (PRAD) using the TCGA pipeline for data isolation and filtering.
Materials and Methods: The project utilized the Genomic Data Commons (GDC) portal to isolate the TCGA, PRAD, TP53, high mutations, and no mutations datasets. The key steps involved in the pipeline were:
1.	Selection of the gene of interest.
2.	Application of filters on the GDC to isolate the relevant datasets.
3.	Following the TCGA pipeline to combine and filter the datasets to create the samDF2 dataframe.
4.	Grouping the data based on high versus no mutations.
Filters we used
●	TCGA
●	TCGA-PRAD
●	TP53 (Gene)
○	High Mutation
○	No Mutation
The measurements used to filter loci are detailed in the pipeline. The software used for this project was R Studio, with the following packages: ggplot2, ggrepel, tidyselect, TCGAbiolinks, DESeq2, and readr.
Supplementary Files: The supplementary files attached to this report include:
1.	highmutationsgdc.tsv: A TSV file containing all the high mutations in TP53 in PRAD.
2.	nomutationsgdc.tsv: A TSV file containing no mutations in TP53 in PRAD.
3.	samDF2: A combined and filtered dataset of the highmutationsgdc.tsv and nomutationsgdc.tsv files.
Group 2 originally had too many samples for simple analysis, so we used 50 random samples to use in downstream analysis.
Results:
●	Samples per grouping: 70, 156
●	Number of genes in final analysis: 29036
●	Number of genes detected to be significantly differentially expressed: 10,818
●	Most significantly differentially expressed gene: ENSG00000123407.4
Analysis
The primary objective of this project was to compare the high versus no mutations in the TP53 gene in prostate cancer (PRAD). The most significantly differentially expressed gene was ENSG00000123407.4. Due to high level of filtering, our PCA plot shows that all of our data is correlated and fits the criteria to be used for analysis
The volcano plot provided represents a comprehensive overview of gene expression data, encompassing a total of 29,306 gene features. The x-axis, delineating the log2(Fold Change), reveals the extent of upregulation or downregulation of genes, with a range extending from -2.5 to 2.5. Genes situated on the right exhibit upregulation, while those on the left are downregulated. The y-axis, indicating the -log10(adjusted p-value), scales from 0 to 6, with higher values signifying greater statistical significance. The color coding of the data points, corresponds to the log10(Mean) expression levels, providing a visual gradient of gene expression intensity. The most critical insights are drawn from the outliers that punctuate the plot’s upper extremes, flagging genes with significant fold changes and low p-values, which are prime candidates for further biological investigation. This plot serves as a pivotal tool for identifying key genes that may play crucial roles in the biological processes under study. 
The supplementary files attached to this report include:
●	highmutationsgdc.tsv: A TSV file containing all the high mutations in TP53 in PRAD.
●	nomutationsgdc.tsv: A TSV file containing no mutations in TP53 in PRAD.
●	Project1.csv (SamDf2): A combined and filtered dataset of the highmutationsgdc.tsv and nomutationsgdc.tsv files.

●	Table_GeneData.csv: All the genes we used and extra information
