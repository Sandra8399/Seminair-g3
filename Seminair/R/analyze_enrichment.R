#' Perform Enrichment Analysis and Visualization
#'
#' This function performs Gene Ontology (GO) and KEGG pathway enrichment analysis
#' for significant differentially expressed genes (DEGs) and generates results
#' as Excel files and plots, including a tree plot for GO term similarity.
#'
#' @param significant_degs A data frame of significant DEGs (row names as SYMBOL IDs).
#' @param output_dir Directory where enrichment results and plots will be saved. Defaults to "output".
#' @return A list containing:
#'   \itemize{
#'     \item `go_results`: Results of GO enrichment analysis.
#'     \item `kegg_results`: Results of KEGG pathway enrichment analysis.
#'   }
#' @examples
#' \dontrun{
#' enrichment_results <- analyze_enrichment(results$significant_degs, "output")
#' }
#' @export
analyze_enrichment <- function(significant_degs, output_dir = "output") {
  if (!dir.exists(output_dir)) dir.create(output_dir)

  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(openxlsx)
  library(enrichplot)
  library(ggplot2)  # Load ggplot2 for ggsave

  # Convert gene SYMBOLs to ENTREZ IDs
  gene_conversion <- bitr(rownames(significant_degs),
                          fromType = "SYMBOL",
                          toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db)
  entrez_ids <- gene_conversion$ENTREZID

  # GO enrichment analysis
  go_results <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                         ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)

  # KEGG enrichment analysis
  kegg_results <- enrichKEGG(gene = entrez_ids, organism = "hsa", pvalueCutoff = 0.05)

  # Save results to Excel
  wb <- createWorkbook()
  addWorksheet(wb, "GO_Enrichment")
  writeData(wb, "GO_Enrichment", go_results@result)
  addWorksheet(wb, "KEGG_Enrichment")
  writeData(wb, "KEGG_Enrichment", kegg_results@result)
  saveWorkbook(wb, file.path(output_dir, "Enrichment_Results.xlsx"), overwrite = TRUE)

  # Generate and save dot plots
  dotplot(go_results, showCategory = 10, title = "Top 10 GO Enrichment Terms")
  ggsave(file.path(output_dir, "Top10_GO_Enrichment.png"))

  dotplot(kegg_results, showCategory = 10, title = "Top 10 KEGG Pathways")
  ggsave(file.path(output_dir, "Top10_KEGG_Pathways.png"))

  # Generate and save tree plot
  if (!is.null(go_results) && nrow(go_results@result) > 1) {
    go_termsim <- pairwise_termsim(go_results)  # Calculate term similarity
    treeplot(go_termsim, showCategory = 10, title = "GO Term Similarity Tree")
    ggsave(file.path(output_dir, "GO_Term_Similarity_Tree.png"))
  }

  # Return results
  list(go_results = go_results, kegg_results = kegg_results)
}

# Usage
enrichment_output_dir <- "enrichment_output"
enrichment_results <- analyze_enrichment(results$significant_degs, enrichment_output_dir)
