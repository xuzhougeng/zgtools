#' Auto Gene Expression Analysis
#'
#'
#' @inheritParams file2dds
#' @param orgDb orgDb, will be used for GO enrichment analysis
#' @param lfcThreashold a non-negative value which specifies a log2 fold change
#'   threshold.
#' @param pAdjustMethod the method to use for adjusting p-values, see ?p.adjust
#' @param pThreshold a numeric to subset results of DESeq by pajust, default is
#'   0.1
#' @param goOnt one of "MF", "BP", and "CC" subontologies. BP: biological process
#'   CC: cellular component MF: molecular function
#' @param goKeytype GO keytype
#' @param keggOrg  KEGG orgnism
#' @param keggKeytype KEGG keytype
#' @param ... additional parameter
#'
#' @return Gene Expression Report
#' @importFrom DESeq2 DESeqDataSetFromTximport
#' @importFrom DESeq2 DESeq
#' @importFrom DESeq2 results
#' @importFrom Matrix rowSums
#' @importFrom magrittr %<>%
#' @importFrom stringr str_split
#' @export auto_ge_analyzer
#'
#' @author Xu Zhougeng
#' @examples
auto_ge_analyzer <- function(filepath,
                             type,
                             colData,
                             design,
                             txDb,
                             orgDb,
                             keyType = "GENEID",
                             txName = "TXNAME",
                             lfcThreashold = 0,
                             pAdjustMethod  = "BH",
                             pThreshold = 0.05,
                             goOnt = "MF",
                             goKeytype = "TAIR",
                             keggOrg = "ath",
                             keggKeytype = "kegg",
                             ...) {

  raw_dds <- file2dds(files, type, col_data, design = ~ treat_time, txDb = txdb)
  dds <- raw_dds[rowSums(counts(raw_dds)) > 1,]

  # differential expresssion analysis
  dds <- DESeq(dds)

  pAdjustMethod %<>% match.arg(
    c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
  ))

  res <- results(dds,
                 lfcThreashold = lfcThreashold,
                 pAdjustMethod = pAdjustMethod)
  # Enrichment Analysis of GO and KEGG
  enrich_results <- enrich_analyzer(
    res = res,
    pThreshold = pThreshold,
    goOrg = orgDb,
    goKeytype = goKeytype,
    goOnt = goOnt,
    keggOrg = keggOrg,
    keggKeytype = keggKeytype,
    pAdjustMethod = pAdjustMethod
  )

  output = list(result = res, enrichment = enrich_results)
  return(output)

}




