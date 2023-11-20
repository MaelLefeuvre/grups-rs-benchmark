#!/usr/bin/env Rscript
# --------------------------------------------------------------------------- #
# ---- Dependencies
library(lubridate)
library(plotly)
library(reshape2)
library(stringr)
library(htmlwidgets)
# --------------------------------------------------------------------------- #
#       Functions

seconds_to_datetime <- function(s) {
  MINUTE   <- 60
  HOUR     <- 60 * 60

  hours    <- floor(s / HOUR)
  minutes <-  floor((s - (hours * HOUR)) / MINUTE)
  seconds <-  floor(s - (hours * HOUR) - (minutes * MINUTE))

  pad <- function(x) ifelse(x < 10, paste0("0", x), as.character(x))

  paste0("2017-01-01 ", pad(hours), ":", pad(minutes), ":", pad(seconds))
  #return `2017-01-01 ${h}:${pad(m)}:${pad(s)}`
}

# Compute an average from a dataframe column
compute_metric_avg <- function(df, metric, ci = 0.95) {
    if (metric == "h.m.s") {
        values <- lubridate::period_to_seconds(lubridate::hms(df[[metric]]))
    } else {
        values <- df[[metric]]
    }
    avg <- mean(values)

    eps <- qnorm((1 - ci) / 2)

    ci <- eps * (sd(values) / sqrt(length(values)))

    list(avg = avg, ci = ci)
}

# Apply 'compute_metric_avg' to a list of dataframes
get_bench_average <- function(bench, metric, ci = 0.95) {
    values <- lapply(bench, FUN = compute_metric_avg, metric = metric)
    values
}

# Return a number of comparisons, given a number of individuals.
pair_to_comparisons <- function(n) {
    (n * (n - 1)) / 2
}

view_plotly <- function(widget) {

  # Generate random file name
  temp <- paste("plots/plotly", "html", sep = ".")

  # Save. Note, leaving selfcontained=TRUE created files that froze my browser
  htmlwidgets::saveWidget(widget, temp, selfcontained = TRUE)

  # Launch with desired application
  system(sprintf("firefox file://%s", paste0(getwd(), "/", temp)))

  # Return file name if it's needed for any other purpose
  temp
}

plot_2D_metric <- function(
  df,
  metric_label = "Average Runtime",
  log_scale    = TRUE,
  supp_trace   = NULL
) {

  if (log_scale) {
      log_args <- list(type = "log", dtick = 0.30102999566)
  } else {
      log_args <- list()
  }

  fig <- plot_ly(
      type    = "scatter",
      mode    = "marker+lines",
      y       = df$value[which(df$L3 == "avg")],
      x       = df$L1[which(df$L3 == "avg")],
      error_y = list(array = df$value[which(df$L3 == "ci")])
  )

  if (length(supp_trace) != 0) {
    fig <- fig %>%  add_trace(y = supp_trace, color = "red")
  }

  fig <- fig %>% layout(
      title = paste(metric_label, "as a function of depth"),
      xaxis = c(
        list(title = list(text = "SNPs")),
        log_args,
        rangemode = "tozero"
      ),
      yaxis = c(
        list(title = list(text = metric_label)),
        log_args,
        rangemode = "tozero"
      )
  )
}

parse_2D_plot_bench_df <- function(
  benches,
  by = "coverage",
  set_other,
  metric = "h.m.s"
) {

  if (!(by %in% c("coverage", "samples"))) {
    return - 1
  }
  metric_values <- lapply(benches, FUN = get_bench_average, metric = metric)
  plot_data     <- reshape2::melt(metric_values)

  other_col <- ifelse(by == "coverage", "L2", "L1")
  plot_data[which(plot_data[[other_col]] == set_other), ]
}

plot_3D_metric <- function(
  benches,
  metric       = "h.m.s",
  metric_label = "Average Runtime",
  log_scale    = TRUE
) {
  metric_values <- lapply(benches, FUN = get_bench_average, metric = metric)
  plot_data <- reshape2::melt(metric_values)

  # Exclude confidence intervals. This is dirty...
  plot_data <- plot_data[which(plot_data$L3 == "avg"), ]
  plot_x <- unique(as.numeric(plot_data$L1))                           # SNPs
  plot_y <- pair_to_comparisons(unique(as.numeric(plot_data$L2) + 1))  # Comparisons # nolint
  plot_z <- matrix(plot_data$value, c(length(plot_y), length(plot_x)))  # Metric

  # ---- Rescale to minutes if doing average runtime.
  if (metric == "h.m.s") plot_z <- plot_z / 60

  # ---- Rescale axis to log scale if requested.
  if (log_scale) {
      log_args <- list(type = "log", dtick = 0.30102999566)
  } else {
      log_args <- list()
  }

  axisfont <- list(size = 16)
  tickfont <- list(size = 14)

  fig <- plotly::plot_ly(
    type     = "surface",
    x        = plot_x,
    y        = plot_y,
    z        = plot_z,
    contours = list(
      x = list(
        show  = TRUE,
        start = min(plot_x),
        end   = max(plot_x),
        color = "white"
      ),
      y = list(
        show  = TRUE,
        start = min(plot_y),
        end   = max(plot_y),
        color = "white"
      ),
      z = list(
        show  = TRUE,
        start = floor(min(plot_z) - 3600) * 1000,
        end   = floor(min(plot_z) - 3600) * 1000,
        color = "white"
      )
    )
  ) %>%
  plotly::layout(
    title = paste(
      metric_label,
      "as a function of depth and number of comparisons."
    ),
    scene = list(
      xaxis = c(log_args,
        list(
          title     = list(text = "SNPs", font = axisfont),
          range     = c(0, 1200000),
          tickfont  = tickfont,
          rangemode = "tozero"
        )
      ),
      yaxis = c(log_args,
        list(
          title     = list(text = "Comparisons", font = axisfont),
          range     = c(0, 55),
          tickfont  = tickfont,
          rangemode = "tozero"
        )
      ),
      zaxis = c(log_args,
        list(
          title      = list(text = metric_label, font = axisfont),
          tickfont   = tickfont,
          colorbar   = list(title = "Average Runtime (minutes)")
        )
      )
    )
  ) %>%
  plotly::config(
      editable             = TRUE,
      displaylogo          = FALSE,
      scrollZoom           = TRUE,
      toImageButtonOptions = list(
        format   = "svg",
        filename = "benchmark"
      )
    )
  fig
}

# ---------------------------------------------------------------------------- #
# ---- Main

usage <- function() {
  script_name <- commandArgs()[0]
  cat("Usage:", script_name, "<benchmarks dir> <output dir> <pattern>\n")
  cat("Example:", script_name, "benchmarks/grups-rs-bench/1X/",
    "./plots/", "'.*[.]0[.]log$'", '\n'
  )
}

main <- function() {
  # ---- Parse arguments
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 3) {
    usage()
    stop("Invalid number of arguments")
  }
  dir  <- args[1]
  out  <- ifelse(grepl("/$", args[2]), args[2], paste0(args[2], "/"))
  snps <- as.numeric(unique(stringr::str_extract(list.files(dir), "[0-9]{2,}")))
  snps <- sort(snps)

  # ---- Load benchmark files
  benches <- list.files(
    path       = dir,
    pattern    = args[3],
    full.names = TRUE,
    recursive  = TRUE
  )

  # ---- Separate bench files across number of samples
  benches <- lapply(snps, FUN = function(x) benches[grep(x, benches)])

  # ---- Subdivide bench files across number of samples
  benches        <- lapply(benches,
    FUN = function(cov) {
      dfs <- lapply(cov, FUN = function(file) {
        read.table(file, header = TRUE)
      })

      names(dfs) <- sapply(cov, FUN = function(x) {
        gsub(".log", "", tail(strsplit(split = "-", x)[[1]], n = 1))
      })
    dfs[order(as.numeric(names(dfs)))]  # Return with numeric ordering.
  })

  # benches is now a 2D list with ['coverage']['n-samples']
  # Ex: benches["4494"]["4"]  == benchmark for coverage 4494 SNPs & 4 samples
  names(benches) <- snps

  # Plot Runtime h.m.s and RSS 3D scatterplot
  runtime_hms <- plot_3D_metric(
    benches      = benches,
    metric       = "h.m.s",
    metric_label = "Average Runtime (minutes)",
    log_scale    = FALSE
  )
  runtime_rss <- plot_3D_metric(
    benches      = benches,
    metric       = "max_rss",
    metric_label = "Max. Resident Set Size (KiB)",
    log_scale    = FALSE
  )

  # ---- Export as html
  if (!dir.exists(out)) dir.create(out)
  htmlwidgets::saveWidget(
    widget        = runtime_hms,
    file          = paste0(out, "runtime-hms-3d-scatterplot.html"),
    selfcontained = TRUE
  )
  htmlwidgets::saveWidget(
    widget        = runtime_rss,
    file          = paste0(out, "runtime-rss-3d-scatterplot.html"),
    selfcontained = TRUE
  )
  ?htmlwidgets::saveWidget
}

if (!interactive()) {
  main()
}
