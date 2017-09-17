#' All kinds of enrichment Analysis
#'
#' @param res result object from DESeq2 results
#' @param pThreshold p value ajust
#' @param goOnt one of "BP", "MF", "CC" or "GO"
#' @param goOrg GO orgnism database
#' @param goKeytype GO keytype
#' @param keggOrg  KEGG orgnism
#' @param keggKeytype KEGG keytype
#' @param pAdjustMethod p value adjust method
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
  message("running enrichment analysis")
  res_sub <- AnnotationDbi::subset(res, padj <pThreshold)
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


#' Enrichment Analysis Results Output
#'
#' Write the enrichment analysis results into csv files
#'
#' @param enrichList a list return from enrich_analyzer
#' @param dirname directory name for storage of results
#' @param orgDb orgniams database for gnen annotation
#' @param column the columns or kinds of things that can be retrieved from the
#'   database. As with keys, all possible columns are returned by using the
#'   columns method.
#' @param keyType the keytype that matches the keys used. For the select
#'   methods, this is used to indicate the kind of ID being used with the keys
#'   argument. For the keys method this is used to indicate which kind of keys
#'   are desired from keys
#'
#' @return
#' @export
#'
#' @examples
write_enrich <- function(enrichList, dirname, orgDb=NULL, column =NULL, keyType = NULL){
  # create directory
  if (!dir.exists(dirname)) dir.create(dirname)

  ## write DESeq2  raw and sigificant resutls
  filepath_x = file.path(dirname, "DESeqResults.csv")
  x <- file(filepath_x, open = "at" )
  filepath_y = file.path(dirname, "DESeqResultsSignificant.csv")
  y <- file(filepath_y, open = "at" )
  ## meta data
  metaData <- S4Vectors::elementMetadata(enrichList$result)
  ## raw result data
  rawData <- as.data.frame(enrichList$result)
  rawData <- rawData[order(rawData$padj), ]
  ## sigificant resutls
  signData <- as.data.frame(
    AnnotationDbi::subset(enrichList$result,
                          padj < 0.05 & abs(log2FoldChange) > 2))
  signData <- signData[order(signData$padj), ]

  if ( ! (! isS4(orgDb) & is.na(column) &  is.na(keyType)) ){
    rawData$symbol <- AnnotationDbi::mapIds(org, keys = rownames(rawData),
                                            column = column, keytype = keyType)
    signData$symbol <- AnnotationDbi::mapIds(org, keys = rownames(signData),
                                            column = column, keytype = keyType)
  }
  write.table(metaData, file = x, sep= ",", fileEncoding = "UTF-8")
  write.table(rawData, file =x, append = TRUE, sep =",", fileEncoding = "UTF-8")
  write.table(metaData, file = y, sep= ",", fileEncoding = "UTF-8")
  write.table(signData, file =y, append = TRUE, sep =",", fileEncoding = "UTF-8")
  close(x)
  close(y)

  # write the enrichmnent results
  for (elist in enrichList$enrichments){
    if (class(elist)[1] == 'enrichResult'){
      # over-repressent analysis(ORA)
      organism <- paste(strsplit(elist@organism," ")[[1]], collapse = "_")
      ontology <- elist@ontology

      filepath  <- file.path(dirname, paste0(organism,"_",ontology, "_enrichResult", ".csv"))

      message(filepath)

      write.table(as.data.frame(elist), file = filepath, sep = ",")
    } else{
      # gene set enrichment analysis(GESA)
      organism <- paste(strsplit(elist@organism," ")[[1]], collapse = "_")
      settype <- elist@setType

      filepath <- file.path(dirname , paste0(organism, "_", settype,"_gseaResult", ".csv"))

      message(filepath)

      write.table(as.data.frame(elist), file = filepath, sep = ",")
    }
  }
}




