#!/usr/bin/env Rscript

usage <- function() {
  exe <- commandArgs(trailingOnly = FALSE)[1]
  cat("Usage:", exe, "<output-file> [<grups-rs-results-dir> ...]", "\n")
}


#' Get the average SVM probability of every Most_Likely_rel
#' @param dir a GRUPS-rs output results directory (must contain a .result and
#'        a .probs file)
#' @return a merged dataframe.
merge_results <- function(dir) {
  results <- list.files(dir, pattern = ".result", full.names = TRUE)[1]
  probs   <- list.files(dir, pattern = ".probs",  full.names = TRUE)[1]
  merge(
    read.table(results, header = TRUE, sep = "\t"),
    read.table(probs,   header = TRUE, sep = "\t")
  )
}
#' Get the average SVM probability of every Most_Likely_rel
#' @param merged a merged dataframe containing the .results and .probs of a
#'        GRUPS-rs run.
#' @param exclude a vector of relationships to exclude from the computation
#'        (e.g.: c("Self", "Unrelated"))
#' @return an average SVM probability
get_average_max_prob <- function(merged, exclude = c()) {
  merged <- merged[!(merged$Most_Likely_rel %in% exclude), ]

  max_svm_probs <- apply(merged, 1, function(row) {
    mlr <- row[["Most_Likely_rel"]]
    mlr <- sub("^([0-9])", "X\\1", mlr, perl = TRUE)
    mlr <- gsub("[-+]", ".", mlr, perl = TRUE)

    as.numeric(row[[mlr]])
  })

  mean(max_svm_probs)
}

#' Create a list of all possible combinations from the elements of a vector
#' @param x a vector
#' @return a list of vectors, each containing a possible combination of the
#'         elements of x
combinations <- function(x) {
  do.call("c", lapply(seq_along(x), function(i) combn(x, i, list)))
}

#' Get the common investigated relationship between a set of merged GRUPS-rs
#' results.
#' @param dfs a list of merged GRUPS-rs output results dataframes.
#' @return a vector of common relatedness levels.
get_common_rels <- function(dfs) {
  Reduce(intersect, lapply(dfs, function(x) colnames(x)[10:ncol(x)]))

}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 2 || !(all(dir.exists(args[2:length(args)])))) {
    usage()
    stop()
  }

  output_file <- args[1]
  input_dirs <- args[2:length(args)]

  merged_dfs <- lapply(input_dirs, FUN = merge_results)


  output         <- data.frame(dir = input_dirs)
  output$max.all <- sapply(merged_dfs, FUN = get_average_max_prob)


  common_rels <- get_common_rels(merged_dfs)

  for (exclude_combination in combinations(common_rels)) {
    label <- paste0(
      "max.exclude(",
      paste(exclude_combination, collapse = "-"),
      ")"
    )
    output[[label]] <- sapply(
      X   = merged_dfs,
      FUN = get_average_max_prob,
      exclude = exclude_combination
    )
  }

  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  write.table(
    x         = output,
    file      = output_file,
    quote     = FALSE,
    sep       = "\t",
    row.names = FALSE)

}