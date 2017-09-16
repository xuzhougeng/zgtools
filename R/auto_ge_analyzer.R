#' Auto Gene Expression Analysis
#'
#'
#' @param files file path of count matrix
#' @param sampleInfo sample infomation
#' @param design experiment design, i.e. design = ~ condition
#' @param txDb txDb, will be to import count matrix
#' @param orgDb orgDb, will be used for GO enrichment analysis
#' @param lfcThreashold a non-negative value which specifies a log2 fold change
#'   threshold.
#' @param pAdjustMethod the method to use for adjusting p-values, see ?p.adjust
#' @param pThreshold a numeric to subset results of DESeq by pajust, default is
#'   0.1
#' @param ont one of "MF", "BP", and "CC" subontologies. BP: biological process
#'   CC: cellular component MF: molecular function
#' @param res result object from DESeq2 results
#' @param pThreshold p value ajust
#' @param goOnt one of "BP", "MF", "CC" or "GO"
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
auto_ge_analyzer <- function(files,
                             sampleInfo,
                             design,
                             txDb,
                             orgDb,
                             lfcThreashold = 0,
                             pAdjustMethod  = "BH",
                             pThreshold = 0.05,
                             goOnt = "MF",
                             goKeytype = "TAIR",
                             keggOrg = "ath",
                             keggKeytype = "kegg",
                             ...) {
  txi <- files2txi(files = files,
                   sampleInfo = sampleInfo,
                   txDb = txDb)
  dds <- DESeqDataSetFromTximport(txi, sample_info, design)
  dds <- dds[rowSums(counts(dds)) > 1,]
  # Explorary data analysis plot
  intgroup = unlist(str_split(as.character(design)[2], ' '))
  eda_plot(dds, intgroup = intgroup[!intgroup == '+'])

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



#' @importFrom readr read_tsv
#' @importFrom tximport tximport
#' @importFrom AnnotationDbi keys
#' @importFrom AnnotationDbi select
files2txi <- function(files, sampleInfo, txDb) {
  # test if files are exists
  if (!all(file.exists(files))) {
    stop("some files are not exitst! ")
  }
  # import data, build a
  k <- keys(txDb, keytype = "GENEID")
  df <- select(txdb,
               keys = k,
               keytype = 'GENEID',
               columns = 'TXNAME')
  tx2gene <- df[, 2:1]
  # suppress the message from tximport during the load files
  suppressMessages(txi <-
                     tximport(files, type = "salmon", tx2gene = tx2gene))
  return(txi)
}
