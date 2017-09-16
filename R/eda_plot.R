#' Exploratory Data Analysis Plot
#'
#'
#' @param dds DESeqDataSet object
#' @param blind logical, whether to blind the transformation to the experimental design.
#' @param intgroup interesting groups: a character vector of names to use for grouping
#' @return graphics
#'
#' @importFrom DESeq2 rlog
#' @importFrom DESeq2 vst
#' @importFrom DESeq2 estimateSizeFactors
#' @importFrom DESeq2 plotPCA
#' @importFrom SummarizedExperiment assay
#' @importFrom magrittr %>%
#' @importFrom dplyr as_data_frame
#' @importFrom dplyr bind_rows
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_hex
#' @importFrom ggplot2 coord_fixed
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggtitle
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#'
#' @export eda_plot
#'
#' @examples
eda_plot <- function(dds, intgroup, blind = FALSE) {
  # begin drawing
  # The rlog and variance stabilizing transformations
  rld <- rlog(dds, blind = blind)
  vsd <- vst(dds, blind = blind)
  dds <- estimateSizeFactors(dds)
  df <- bind_rows(
    as_data_frame(log2(counts(dds, normalized = TRUE)[, 1:2] + 1)) %>%
      mutate(transformation = "log2(x+1)"),
    as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
    as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst")
  )
  colnames(df)[1:2] <- c("x", "y")
  p1 <-
    ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) + coord_fixed() +
    facet_grid(. ~ transformation)

  # distance heatmap
  samplesDists <- dist(t(assay(rld)))
  sampleDistMatrix <- as.matrix(samplesDists)
  rownames(sampleDistMatrix) <- colnames(dds)
  colnames(sampleDistMatrix) <- colnames(dds)
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  pheatmap(
    sampleDistMatrix,
    clustering_distance_rows = samplesDists,
    clustering_distance_cols = samplesDists,
    col = colors
  )
  # PCA
  p2 <-
    plotPCA(rld, intgroup = intgroup) + ggtitle("PCA plot using rlog-transformed values")
  print(p2)

  # # MDS
  # mds <- as.data.frame(colData(rld))  %>%
  #     cbind(cmdscale(sampleDistMatrix))
  # p3 <- ggplot(mds, aes(x = `1`, y = `2`, color = dex, shape = cell)) +
  #     geom_point(size = 3) + coord_fixed()
  # print(p3)

  dev.off()
}
