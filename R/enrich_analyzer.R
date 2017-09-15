#' All kinds of enrichment Analysis
#'
#' @param res result object from DESeq2 results
#' @param pThreshold p value ajust
#' @param goOnt one of "BP", "MF", "CC" or "GO"
#' @param goOrg GO orgnism database
#' @param goKeytype GO keytype
#' @param keggOrg  KEGG orgnism
#' @param keggKeytype KEGG keytype
#' @param pAdjustMethod p value ajust method
#'
#' @importFrom clusterProfiler enrichGO
#' @importFrom clusterProfiler gseGO
#' @importFrom clusterProfiler enrichKEGG
#' @importFrom clusterProfiler gseKEGG
#' @return
#' @export enrich_analyzer
#'
#' @examples
enrich_analyzer <- function(res,
                            pThreshold,
                            goOrg,
                            goKeytype = "TAIR",
                            goOnt = "MF",
                            keggOrg = "ath",
                            keggKeytype = "kegg",
                            pAdjustMethod = "BH",
                            ...) {
  # Enrichment Analysis
  res_sub <- subset(res, padj < pThreshold)
  genes <- rownames(res_sub)

  ## GO Enrichment Analysis
  ego <- enrichGO(
    gene = genes,
    keyType = goKeytype,
    OrgDb = goOrg,
    ont = goOnt,
    pAdjustMethod = pAdjustMethod
  )

  # KEGG Enrichment Analysis
  ekegg <- enrichKEGG(
    gene = genes,
    organism = keggOrg,
    keyType = keggKeytype,
    pAdjustMethod  = pAdjustMethod
  )

  # GSEA
  ## prepare the gene list
  genelist <- res$log2FoldChange
  names(genelist) <- rownames(res)
  genelist <- sort(genelist, decreasing = TRUE)

  gsego <- gseGO(
    geneList = genelist,
    ont = goOnt,
    OrgDb = goOrg,
    keyType = goKeytype,
    pAdjustMethod = pAdjustMethod
  )

  gsekegg <- gseKEGG(
    geneList = genelist,
    organism = keggOrg ,
    keyType = keggKeytype,
    pAdjustMethod = pAdjustMethod,
    verbose = TRUE
  )
  # combine all enrichment anaylsis results
  enrich_list <- list(
    go_enrich = ego,
    go_gsea = gsego,
    kegg_enrich = ekegg,
    kegg_gsea = gsekegg
  )
  return(enrich_list)
}
