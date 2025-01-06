# CLL-1 as a cellular immunotherapy target in juvenile myelomonocytic leukemia

# Description
This repository contains codes specifically designed for the data analysis performed in Werner et al, 2025. The repository provides scripts for analyzing and recreating figures pertaining to the publication. Note this version cited in the paper will not have the updated PMID but later versions will be updated. 

* PMID: not available
* PMCID: not available
* DOI: not available
  
# Features
## Bulk RNAseq and scRNAseq: 
  + There are two R markdowns to generate all figures pertaining to the RNAseq and scRNAseq for the manuscript.
  + There is a single auxiliary script that will provide all the custom functions to run the above markdowns. 
  + Moreover, there is a single config file that provides the paths to all the necessary input files. 
    + Please note that due to data sensitivity and licensing restrictions, we cannot provide all the files mentioned in this manuscript. Users must manually download these files using the source references provided here or in the manuscript. For example, all raw data must be downloaded as cited in the original manuscript. Preprocessed data can be requested from the Stieglitz Lab.
# Data
The dataset used for the Werner et al study can be accessed from dbGaP: phs002504.v3.p3.

# Minimum Requirements
  * Download and install R from here. R (4.1.2)
  * Download the dataset dbGaP: phs002504.v3.p3 and all associated inputs. 
  * Main R packages
 
    + edgeR (3.28.1)
    + limma package(3.42.2)
    + Seurat (4.3.0.1)
    + Harmony (0.10)

  <details>
  <summary>Click to expand Bioinformatics and Genomics packages!</summary>
  
  ### Genomic Data Annotation:
  1. `library(biomaRt)`: Tools for BioMart databases (like Ensembl).
  2. `library(BSgenome)`: Infrastructure for Bioconductor packages using large-scale genomic or other data.
  3. `library(org.Hs.eg.db)`: Mapping information for human genes.
  4. `library(GenomicFeatures)`: Tools for making and manipulating transcript centric annotations.

  ### Genomic Data Analysis (Omics):
  1. `library(DESeq2)`: Differential gene expression analysis based on the negative binomial distribution.
  2. `library(edgeR)`: Empirical analysis of digital gene expression data in R.
  3. `library(GenomicRanges)`: Representations and manipulations of genomic intervals and variables defined along a genome.
  4. `library(GSVA)`: Gene set variation analysis for microarray and RNA-seq data.
  5. `library(Gviz)`: Plotting data and annotation information along genomic coordinates.
  6. `library(pathview)`: Plots pathway maps and overlays experimental data.
  7. `library(ggbio)`: Visualization tools for genomic data.
 
  
  ### Heatmaps and Clustering:
  1. `library(clusterProfiler)`: Statistical analysis and visualization of functional profiles for genes and gene clusters.
  2. `library(ComplexHeatmap)`: Making complex heatmaps.
  3. `library(d3heatmap)`: Interactive heatmaps.
  4. `library(dendextend)`: Extending R's dendrogram functionality.
  5. `library(dendroextras)`: Extra functions to cut, label and colour dendrogram clusters.
  6. `library(parallelDist)`: Parallel distance matrix computation.
  
  ### Visualization:
  1. `library(corrplot)`: Visualization of a correlation matrix.
  2. `library(factoextra)`: Extract and visualize the results of multivariate data analyses.
  3. `library(ggdendro)`: Create dendrograms using ggplot.
  4. `library(ggplot2)`: An implementation of the Grammar of Graphics.
  5. `library(ggplotify)`: Convert plot function call to 'ggplot' objects.
  6. `library(ggpubr)`: 'ggplot2' based publication ready plots.
  7. `library(ggpval)`: Annotate statistical significance onto 'ggplot' objects.
  8. `library(ggrepel)`: Automatically position non-overlapping text labels with 'ggplot2'.
  9. `library(gplots)`: Various R programming tools for plotting data.
  10. `library(gridExtra)`: Miscellaneous functions for "grid" graphics.
  11. `library(kableExtra)`: Build complex HTML or 'LaTeX' tables using 'kable()' and pipe syntax.
  12. `library(patchwork)`: The composer of ggplots.
  13. `library(RColorBrewer)`: ColorBrewer palettes.
  14. `library(VennDiagram)`: Generate high-resolution Venn and Euler plots.
  15. `library(Vennerable)`: Venn and Euler area-proportional diagrams.
  16. `library(wesanderson)`: Wes Anderson color palettes.
  17. `library(igraph)`: Network analysis and visualization.
  18. `library (ggbeeswarm)` # Beeswarm plots helper
  19. `library(forestplot)` # forest plot helper, mostly use in meta-analysis
  20. `library (ggridges)` # Ridgeline plots 
  21. `library(cowplot)` # functions to align plots and arrange them into complex compound figures
  
  ### Statistical Analysis:
  1. `library(FactoMineR)`: An R package for multivariate analysis.
  2. `library(fgsea)`: Fast gene set enrichment analysis.
  3. `library(MASS)`: Functions and datasets to support Venables and Ripley's MASS.
  4. `library(matrixStats)`: Functions that apply to rows and columns of matrices (and to vectors).
  5. `library(PerformanceAnalytics)`: Econometric tools for performance and risk analysis.
  6. `library(psych)`: Procedures for psychological, psychometric, and personality research.
  7. `library(survival)`: Survival analysis.
  8. `library(survminer)`: Drawing survival curves using 'ggplot2'.
  9. `library(vegan)`: Community Ecology Package.
  10. `library(scales)`: Scale functions for visualization.
  11. `library(Rtsne)`: T-distributed stochastic neighbor embedding using a Barnes-Hut implementation.
  12. `library(umap)`: Uniform Manifold Approximation and Projection.
  
  ### Data Manipulation:
  1. `library(data.table)`: Extension of `data.frame`.
  2. `library(dplyr)`: A grammar of data manipulation.
  3. `library(DT)`: A wrapper of the JavaScript library 'DataTables'.
  4. `library(forcats)`: Tools for working with categorical variables (factors).
  5. `library(plyr)`: Tools for splitting, applying and combining data.
  6. `library(reshape)`: Flexibly reshape data.
  7. `library(stringr)`: Simple, consistent wrappers for common string operations.
  8. `library(tidyr)`: Easily tidy data with 'spread()' and 'gather()' functions.

  ### Document Generation and Reporting:
  1. `library(knitr)`: A general-purpose tool for dynamic report generation in R.
  2. `library(pander)`: An R Pandoc writer.
  3. `library(stargazer)` # LATEX, HTML and ASCII tables from R statistical output
  
  ### File I/O:
  1. `library(openxlsx)`: Read, write and edit XLSX files.
  

  
</details>

## Additional public datasets

| Description | Type     | Source | Link |
| ---------------------------- | ------------------ | ----------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Pathway                      | database           | MSigDB                  | [GSEA](https://www.gsea-msigdb.org/gsea/msigdb)                                                                                                                                                                                                                                                                                                                                                     |
| Normal RNAseq             | Normal Counts  | Genotype-Tissue Expression (GTEx)              | [Portal](https://gtexportal.org/home/)             |
| scRNAseq Normal                   | Chen et al, 2022    | PMID: 34864916                 | [PUBMED](https://pubmed.ncbi.nlm.nih.gov/34864916/)                                                                                                                                                                                     |
| RNAseq index                 | Genomic index      | GENCODE                 | [INDEX](ftps://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz)                                                                                                                                                                                                                     |
| RNAseq GTF                   | Gene Annotation    | GENCODE                 | [GTF](ftps://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz)|
| Annotation of cell types     |  Ianevski et al, 2022 | ScType | [Github](https://github.com/IanevskiAleksandr/sc-type) |
| Cell Surface Markers                | The Human Protein Atlas, Bausch-Fluck et al, 2015 and Cancer Surfaceome Atlas | Multiple |  [The Human Protein Atlas](https://www.proteinatlas.org/),<br> [PMID: 25894527](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0121314 ), <br> [Cancer Surfaceome Atlas, Hu et al, 2021](http://fcgportal.org/TCSA)
