
# RNA-seq Enrichment Analysis

This R package provides a pipeline for analyzing RNA-seq count data, including filtering, differential expression analysis, and enrichment analysis. 

### collaborators: 
- Sandra Desatová
- Zozan Ismail

  
## Features

- **Filtering**: Removes lowly expressed genes using `edgeR::filterByExpr`.
- **Differential Expression Analysis**: Identifies differentially expressed genes (DEGs) using `DESeq2`.
- **Enrichment Analysis**:
  - Performs Gene Ontology (GO) and KEGG pathway enrichment using `clusterProfiler`.
  - Converts gene SYMBOLs to ENTREZ IDs using `org.Hs.eg.db`.
- **Visualization**:
  - Generates dot plots for enriched pathways.
  - Creates GO term similarity tree plots.
- **Export**:
  - Saves results in user-friendly formats (Excel files and PNG plots).

## Installation

To use this package, you need to clone the repository and install the dependencies.

```r
# Clone the repository
git clone https://github.com/your-repo-name.git

# Install required R packages
install.packages(c("edgeR", "DESeq2", "clusterProfiler", "org.Hs.eg.db", "openxlsx", "enrichplot", "ggplot2"))
