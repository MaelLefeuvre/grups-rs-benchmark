#!/usr/bin/env Rscript
require(ggplot2)
require(caret)
require(viridis)
require(plyr)
require(svglite)

.KINSHIP_LEVELS <- c("Self", "First", "Second", "Third", "Unrelated")

.msg <- function(..., end = "\n", shift = 0){
  if (!is.character(end)) stop("Invalid end character:", end)
  cat(rep(" ", times = shift), ..., end)
}

.load_data_frame <- function(path=NULL, header = TRUE, sep = "\t", ...) {
  if (is.null(path)) stop("function .load_data_frame() must be given a path")
  .msg("Loading file:", path)
  read.table(path, header = header, sep = sep, ...)
}

#' Convert a vector of character into a predefined set of factorized levels
.sanitize_kinship_levels <- function(x) {
  patterns <- c(
    "Twins_or_self", 
    "Siblings|Parent_child|First_degree",
    "Half-siblings|Avuncular|GPGC|Second_degree",
    "Cousins|Third_degree"
  )

  for (i in seq_along(patterns)) {
    x <- sub(patterns[i], .KINSHIP_LEVELS[i], x, perl = TRUE)
  }
  return(x)
}

#' Load a prediction file with column <Ind1-id> <Ind2>
load_expected_results_file <- function(path = NULL){
  .load_data_frame(path)
}

load_grups_results_file <- function(path = NULL){
  df                 <- .load_data_frame(path)
  df$Most_Likely_rel <- .sanitize_kinship_levels(df$Most_Likely_rel)
  df
}

.str_extract <- function(pattern, x, global = FALSE, perl = TRUE, ...) {
  regex_func <- ifelse(global, gregexpr, regexpr)
  regmatches(x, regex_func(pattern, x, perl = perl, ...))
}

#' parse the pedigree id from a column of strings
#' @param x a vector of character
parse_pedigree_id <- function(x, pattern = "ped[0-9]+_", perl = TRUE) {
  ids <- .str_extract(pattern, x)
  id  <- unique(ids)
  if (length(id) > 1)
    stop("Invalid file: column contains multiple pedigree id.")
  return(id)
}

check_results <- function(expect_df, results_df) {
  ped_id <- parse_pedigree_id(results_df$Pair_name)

  expect_df$Pair_name <- with(expect_df,
    paste0(ped_id, pair_1, "-", ped_id, pair_2)
  )

  merged_df <- merge(
    x    = expect_df,
    y    = results_df,
    by.x = "Pair_name",
    by.y = "Pair_name"
  )

  output <- data.frame(
    pair = merged_df$Pair_name,
    true = factor(merged_df$relatedness,     levels = .KINSHIP_LEVELS),
    pred = factor(merged_df$Most_Likely_rel, levels = .KINSHIP_LEVELS),
    overlap = results_df$Corr.Overlap

  )
  output$TP   <- output$true == output$pred
  output$diff <- as.numeric(output$true) - as.numeric(output$pred)
  return(output)
}

bind_predictions <- function(expect_df, prediction_paths) {
  predictions <- lapply(prediction_paths, FUN=function(path){
    .msg("Checking results in", path)
    id           <- .str_extract("(?<=[.])[0-9]+(?=[.]result$)", path)
    coverage     <- .str_extract(paste0("(?<=[-])[0-9]+(?=[.]", id, "[.]result$)"), path)
    out          <- check_results(expected_df, load_grups_results_file(path))
    out$coverage <- coverage
    out$id       <- id 
    out
  })

  .msg("Binding predictions")
  do.call("rbind", predictions)
}

get_cm <- function(pred, true) {
  caret::confusionMatrix(data = pred, reference = true)
}

cms_to_barplot <- function(cms, class, text.size = 11) {
  df <- plyr::ldply(
    lapply(cms, FUN=function(cm) {
      df <-as.data.frame(cm$table);
      df[df$Reference == class, ]
    })
  )
  df$diff <- factor(as.numeric(df$Prediction) - as.numeric(df$Reference))
  df$.id  <- factor(df$.id, levels = as.numeric(unique(df$.id)))
  df$text.col <- ifelse(df$Prediction %in% c("Unrelated", "Third"), "white", "black")

  accuracies <- plyr::ldply(lapply(
    cms, FUN = function(cm) {
      round(cm$byClass[paste("Class:", class), "Balanced Accuracy"], 2)
    }
  ))

  geom_text_size <- (text.size - 2) / (14/5)

  ggplot2::ggplot(data = df) +
  ggplot2::geom_bar(
    aes(x = .id, y = Freq, fill = Prediction),
    position = position_fill(reverse = TRUE),
    stat     = "identity"
  ) +
  ggplot2::geom_text(
    aes(
      x      = .id,
      y      = Freq,
      label  = ifelse(Freq == 0, "", Freq),
      color  = text.col,
      family = "sans"
    ),
    position = position_fill(vjust = 0.5, reverse = TRUE),
    size     = geom_text_size,
  ) +
  scale_color_manual(values = c("black", "white")) + 
  viridis::scale_fill_viridis(discrete = TRUE, option = "D", direction = -1) +
  guides(color = "none") + 
  ggplot2::geom_text(
    data = accuracies,
    aes(x = .id, y = 1.1, label = V1),
    hjust  = 1,
    size   = geom_text_size,
    family = "sans"
  ) +
  ggplot2::coord_flip(clip = "off", ylim = c(0, 1)) +
  ggplot2::scale_y_reverse(expand = c(0, 0), breaks = seq(0, 1, 0.1)) +
  ggplot2::theme_classic() +
  ggplot2::xlab("Overlap") +
  ggplot2::ylab("Frequency") +
  ggplot2::ggtitle(class) +
  ggplot2::theme(
    legend.position = "bottom",
    plot.margin = unit(c(1, 4, 1, 1),  "lines"),
    axis.text.x = element_text(size = text.size, family = "sans"),
    axis.text.y = element_text(size = text.size, family = "sans"),
    plot.title  = element_text(size = text.size + 2, family = "sans", face = "bold")
  )

}

# --------------------------------------------------------------------------- #
# ---- Main

usage <- function(){
  exe <- commandArgs(trailingOnly = FALSE)[0]
  cat(exe, ": <output-dir> <expected-results.tsv> [grups-rs-output.result ...] \n")
}

if (sys.nframe() == 0) {
  args             <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 2 || any(c("--help", "-h") %in% args)) {
    usage()
  }

  output_dir       <- args[1]
  expected_path    <- args[2]
  prediction_paths <- args[3:length(args)]

  expected_df <- load_expected_results_file(expected_path)
  predictions <- bind_predictions(expected_df, prediction_paths)

  # compute the overall overlap mean
  average_overlap <- aggregate(
    overlap~coverage,
    data = predictions,
    FUN  = function(x) round(as.numeric(mean(x)))
  )
  apply(average_overlap, 1 , function(x) {
    predictions$overlap[
      predictions$coverage == x[["coverage"]]
    ] <<- as.numeric(x[["overlap"]])
  })

  avg_accuracy <- aggregate(TP~coverage + true, predictions, FUN = mean)


  coverages <- sort(as.numeric(unique(predictions$coverage)))
  overlaps  <- sort(as.numeric(unique(predictions$overlap)))
  cms       <- lapply(setNames(overlaps, overlaps), FUN = function(cov){
    with(predictions,
      get_cm(pred = pred[overlap == cov], true = true[overlap == cov])
    )
  })

  print(cms)

  for (rel in unique(predictions$true)) {
    ggsave(
      filename = paste0(output_dir, "/accuracy-barplot-", rel, ".png"),
      plot     = cms_to_barplot(cms, rel),
      width    = 17,
      height   = 7,
      units    = c("cm"),
      dpi      = 1200
    )
  }
}
