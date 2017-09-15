#' Title
#'
#' @param seqA
#' @param seqB
#' @param template
#' @param winSize
#' @param slide
#' @param ...
#'
#' @return
#' @export windowCompare
#' @importFrom tidyr gather
#'
windowCompare <- function(seqA, seqB, template = NULL, winSize = 1000, slide = 500, ...) {
    contigs <- template[, 1]
    win_list <- list()
    for (contig in contigs) {
        contigs_len <- template[template[, 1] == contig, 2]
        win_list[[contig]] <- make_wins(contigs_len, winSize, slide)
    }
    win_list

}




make_wins <- function(seq_len, winSize = 1000, slide = 500) {
    if (seq_len <= winSize + slide) {
        win.start <- 1
    } else {
        win.start <- seq(1, seq_len, slide)
    }
    win.end <- win.start + winSize
    win.mid <- win.start + slide/2
    data_frame(starts = win.start, ends = win.end, means = win.mid)
}




windows_cal_ <- function(tbl_wds) {
    win_start <- tbl_wds[1]
    win_end <- tbl_wds[2]
    win_mid <- win_start + win_end/2
    win_df <- filter(tbls, DP > DP_low & DP < DP_up & quals > qual) %>% filter(pos >= win_start & pos <
        win_end)
    win_maf_mean <- mean(win_df$MAF)
    loci_num <- nrow(win_df)
    win_sd <- sd(win_df$MAF)
    return(tibble(win_mid, win_maf_mean, win_sd, loci_num))
}


