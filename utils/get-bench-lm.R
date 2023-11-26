#!/usr/bin/env Rscript

if (!require(plyr)) {
  install.packages("plyr")
  library(plyr)
}

if (!require(dplyr)) {
    install.packages("dplyr")
    library(plyr)
}

usage <- function(){
  exe <- sub("--file=", "", commandArgs()[4])
  cat("Usage:", exe, "<benchmarks-dir> [output-file]", '\n')
}

.str_extract <- function(pattern, x, global = FALSE, perl = TRUE, ...) {
  regex_func <- ifelse(global, gregexpr, regexpr)
  regmatches(x, regex_func(pattern, x, perl = perl, ...))
}

sci <- function(x){
  formatC(x, format = "e", digits = 2)
}

if (sys.nframe() == 0) {

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) stop(usage())

  bench_files <- list.files(args[1], pattern="*.0.log", recursive = TRUE, full.names = TRUE)

  bench_data <- do.call("rbind", lapply(bench_files, function(x) {
    df      <- read.table(x, header = TRUE)
    cov     <- .str_extract("/[0-9]+/", x) %>% gsub("/", "", .)
    samples <- .str_extract("(?<=[-])[0-9]+(?=.0.log)", x)
    df$cov <- as.numeric(cov)
    df$samples <- as.numeric(samples)
    df 
  }))

  head(bench_data)

  hms_lm <- lm(s ~ 0 + cov * samples, bench_data)
  rss_lm <- lm(max_rss ~ cov * samples, bench_data)

  print(summary(hms_lm))
  print(summary(rss_lm))

  print(paste0(
    "runtime(seconds) = ",
    sci(hms_lm$coefficients["cov"]), "g", " + ",
    sci(hms_lm$coefficients["samples"]), "π", " + ",
    sci(hms_lm$coefficients['cov:samples']), "πg"
  ))

  print(paste0(
    "memory(KiB) = ",
    sci(rss_lm$coefficients['(Intercept)']), " + ",
    sci(rss_lm$coefficients["cov"]), "g", " + ",
    sci(rss_lm$coefficients["samples"]), "π", " + ",
    sci(rss_lm$coefficients['cov:samples']), "πg"
  ))

  if (length(args>1)) {
    output_file <- args[2]
    write.table(
      x         = bench_data,
      file      = output_file,
      quote     = FALSE,
      sep       = "\t",
      row.names = FALSE
    )
  }
}