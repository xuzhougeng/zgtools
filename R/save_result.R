
#' Annotate genes
#' 
#'
#' @param result result from DESeq2 results
#' @param path which path to store data 
#' @param org organism database, such as org.Hs.eg.db
#'
#' @return
#' @export
#'
#' @examples
save_result <- function(result, path, org){
  ID <- rownames(result)
  SYMBOL <- AnnotationDbi::mapIds(org, ID, column = "SYMBOL", keytype ="TAIR", multiVals = "first" )
  res <- as.data.frame(result) %>% 
    dplyr::mutate(ID = rownames(result), SYMBOL = SYMBOL) %>% 
    dplyr::select(ID, SYMBOL, baseMean:padj) %>% 
    dplyr::arrange(padj)
  write.csv(res, path)
}