#' @title Building DESeqDataSet
#'
#' @description  \code{file2dds}  build \code{DESeqDataSet} base on  a vector
#'   input of sample file path, its counter and experiment design for use of
#'   \code{\link[DESeq]{DESeq2}}
#'
#' @param filepath a \code{vector} of file path of gene count
#' @param type method to count gene, one of HTSeq, Salmon, Sailfish,
#'   kallisto
#' @param colData a \code{data.frame} containning sample information with
#'   columns of sample name, treatment
#' @param design a \code{formula} which expresses how the counts for each gene
#'   depend on the variables in \code{colData}. Many R \code{formula} are valid,
#'   e.g, ~ group, ~ group + condition, ~genotype + treatment + genotype:
#'   treatment.
#' @param txDb a \code{\link[TxDb]{GenomicFeatures}} object, which could be
#'   downloaded from \code{\link[AnnotationHub-objects]{AnnotationHub}}.
#' @param keyType
#' @param txName
#' @param ... other paramter from \code{\link[tximport]{tximport}}
#'
#' @return DESeqDataSet object
#' @importFrom magrittr %<>%
#' @importFrom readr read_tsv
#'
#' @examples
#' @export
file2dds <- function(filepath, type, colData, design, txDb=NULL,
                     keyType = "GENEID", txName = 'TXNAME', ...) {
  # test if files are exists
  if (!all(file.exists(filepath))) {
    stop("some files are not exitst! ")
  }

  # file type decide
  type <- tolower(type)
  type %<>% match.arg(choices = c("htseq", "salmon", "sailfish","kallisto"))

  # build DESeqDataSet
  if (type == "htseq"){
    print("unfinished")
    return(NULL)
  } else{
    # tx2gene
    if ( is.null(txDb) ) stop("txdb not exists")
    k <- AnnotationDbi::keys(txDb, keytype = keyType)
    df <- AnnotationDbi::select(txDb,
                 keys = k,
                 keytype = keyType,
                 columns = txName)
    tx2gene <- df[, 2:1]
    suppressMessages(txi <- tximport::tximport(files = filepath,
                              type = type,
                              tx2gene = tx2gene,
                              importer = function(x) read_tsv(x)))
    dds <- DESeq2::DESeqDataSetFromTximport(txi, colData, design )
    return(dds)

  }


}


