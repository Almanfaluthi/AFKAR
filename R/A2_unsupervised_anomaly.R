#' @title Clinical Anomaly and Outlier Detection (Unsupervised)
#'
#' @description
#' Identifies multivariate anomalies in clinical or administrative datasets
#' (e.g., medication errors, abnormal lab results, or billing fraud)
#' using Mahalanobis distance. Includes dual visualizations to demonstrate
#' the mathematical cutoff process and the scatter distribution.
#'
#' @param data A data.frame containing the clinical or hospital management data.
#' @param features Character vector specifying the numeric columns to use for anomaly detection.
#' @param p_threshold Numeric. The Chi-Square p-value threshold to flag an anomaly (default: 0.001).
#' @param plot_result Logical. If TRUE, generates a 1x2 panel plot highlighting anomalies.
#' @param colors Character vector of length 2 for plot colors c("Normal", "Anomaly").
#' @param reverse_palette Logical. If TRUE, reverses the order of the color palette.
#' @param save_plot Logical. If TRUE, saves the generated plot as a high-resolution 4K PNG.
#' @param save_prefix Character. Optional prefix string for the saved plot filename.
#' @param footer Logical. If TRUE, adds a reproducible generation footer to the plot.
#' @param cex Numeric. Text and symbol scaling factor.
#' @param font Numeric. Font style (1 = plain, 2 = bold, 3 = italic, 4 = bold italic).
#' @param family Character. Font family to use for plotting (e.g., "sans", "serif").
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{data_scored}: The original data frame appended with distance, p-value, and Anomaly_Flag.
#'   \item \code{anomaly_indices}: Row indices of detected anomalies.
#'   \item \code{summary}: A numeric vector summarizing normal vs anomalous counts.
#' }
#'
#' @importFrom stats cov mahalanobis pchisq qchisq complete.cases
#' @importFrom graphics plot points legend mtext par abline
#' @importFrom grDevices png dev.copy dev.off
#'
#' @export
#'
#' @examples
#' # Run anomaly detection with process visualization
#' result <- A2_unsupervised_anomaly(
#'   data = df_A2,
#'   features = c("Dose_mg", "Frequency_per_day"),
#'   p_threshold = 0.001,
#'   plot_result = TRUE,
#'   save_plot = FALSE
#' )

A2_unsupervised_anomaly <- function(data,
                                    features,
                                    p_threshold = 0.001,
                                    plot_result = TRUE,
                                    colors = c("#1f77b4", "#d62728"),
                                    reverse_palette = FALSE,
                                    save_plot = FALSE,
                                    save_prefix = "",
                                    footer = TRUE,
                                    cex = 1,
                                    font = 1,
                                    family = "sans") {

  # 1. Robust Error Handling & Input Validation
  stopifnot("Error: 'data' must be a data.frame." = is.data.frame(data))
  stopifnot("Error: 'features' must be a character vector." = is.character(features))

  missing_cols <- setdiff(features, names(data))
  if (length(missing_cols) > 0) {
    stop(paste("Error: The following features are not found in data:", paste(missing_cols, collapse = ", ")))
  }

  for (f in features) {
    if (!is.numeric(data[[f]])) {
      stop(paste("Error: Feature '", f, "' must be numeric.", sep = ""))
    }
  }

  # 2. Data Preparation
  df_subset <- data[, features, drop = FALSE]
  comp_cases <- stats::complete.cases(df_subset)

  if (sum(comp_cases) < length(features) + 1) {
    stop("Error: Not enough complete observations to calculate covariance matrix.")
  }
  df_clean <- df_subset[comp_cases, , drop = FALSE]

  # 3. Core Algorithm: Mahalanobis Distance
  center <- colMeans(df_clean)
  cov_mat <- stats::cov(df_clean)

  if (det(cov_mat) == 0 || is.na(det(cov_mat))) {
    stop("Error: Covariance matrix is singular. Features might be collinear.")
  }

  m_dist <- stats::mahalanobis(x = df_clean, center = center, cov = cov_mat)
  p_values <- stats::pchisq(m_dist, df = length(features), lower.tail = FALSE)

  # Hitung batas matematis Chi-Square untuk garis bantu di plot
  chi_cutoff <- stats::qchisq(p_threshold, df = length(features), lower.tail = FALSE)

  flags <- ifelse(p_values < p_threshold, "Anomaly", "Normal")

  # Merge back to original data structure
  data_scored <- data
  data_scored[["Mahalanobis_Dist"]] <- NA
  data_scored[["P_Value"]] <- NA
  data_scored[["Anomaly_Flag"]] <- NA

  data_scored[comp_cases, "Mahalanobis_Dist"] <- m_dist
  data_scored[comp_cases, "P_Value"] <- p_values
  data_scored[comp_cases, "Anomaly_Flag"] <- flags

  anomaly_idx <- which(data_scored[["Anomaly_Flag"]] == "Anomaly")

  # 4. Strict Visualization (Plotting & Export)
  if (plot_result) {
    if (reverse_palette) {
      colors <- rev(colors)
    }

    # Simpan par lama untuk dikembalikan nanti
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))

    # KUNCI PERBAIKAN: mfrow (1x2 panel) dan oma (Outer Margin Area)
    # mar = margin dalam (bawah, kiri, atas, kanan)
    # oma = margin luar global (memberi ruang 3 baris di bawah untuk footer)
    graphics::par(mfrow = c(1, 2),
                  mar = c(5, 4, 4, 2) + 0.1,
                  oma = c(3, 0, 0, 0),
                  family = family, font = font, cex = cex)

    # ================= PLOT 1: PROSES DISTRIBUSI =================
    # Mengurutkan titik dari jarak terkecil hingga terbesar
    sorted_idx <- order(m_dist)
    sorted_dist <- m_dist[sorted_idx]
    sorted_flags <- flags[sorted_idx]

    plot_colors_1 <- ifelse(sorted_flags == "Anomaly", colors[2], colors[1])
    plot_pch_1 <- ifelse(sorted_flags == "Anomaly", 19, 1)

    graphics::plot(sorted_dist,
                   col = plot_colors_1,
                   pch = plot_pch_1,
                   xlab = "Ranked Observations (Normal to Extreme)",
                   ylab = "Mahalanobis Distance",
                   main = "Anomaly Detection Process")

    # Tambahkan patahan/batas (Cut-off)
    graphics::abline(h = chi_cutoff, col = colors[2], lty = 2, lwd = 2)
    graphics::legend("topleft",
                     legend = paste("Chi-Square Cutoff (p =", p_threshold, ")"),
                     col = colors[2], lty = 2, lwd = 2, bty = "n", cex = 0.8)

    # ================= PLOT 2: SCATTER PLOT =================
    if (length(features) >= 2) {
      x_var <- df_clean[[features[1]]]
      y_var <- df_clean[[features[2]]]
      xlab_text <- features[1]
      ylab_text <- features[2]
    } else {
      x_var <- seq_along(df_clean[[features[1]]])
      y_var <- df_clean[[features[1]]]
      xlab_text <- "Index"
      ylab_text <- features[1]
    }

    plot_colors_2 <- ifelse(flags == "Anomaly", colors[2], colors[1])
    plot_pch_2 <- ifelse(flags == "Anomaly", 19, 1)

    graphics::plot(x_var, y_var,
                   col = plot_colors_2,
                   pch = plot_pch_2,
                   xlab = xlab_text,
                   ylab = ylab_text,
                   main = "Clinical Feature Scatter Space",
                   sub = paste("Anomalies Identified: ", length(anomaly_idx)))

    graphics::legend("topright",
                     legend = c("Normal", "Anomaly"),
                     col = c(colors[1], colors[2]),
                     pch = c(1, 19),
                     bty = "n")

    # ================= GLOBAL FOOTER (TIDAK AKAN TABRAKAN) =================
    if (footer) {
      footer_text <- paste0("Generated by R-Studio (", R.version.string, ") on ", format(Sys.time(), "%B %d, %Y at %H:%M:%S"))
      # outer = TRUE melempar teks ini ke kanvas terluar (oma) yang sudah kita siapkan
      graphics::mtext(footer_text, side = 1, line = 1, outer = TRUE, adj = 0.5, cex = 0.8)
    }

    # ================= HIGH-RES 4K EXPORT =================
    if (save_plot) {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      clean_prefix <- ifelse(nchar(save_prefix) > 0, paste0("_", save_prefix), "")
      filename <- paste0("A2_unsupervised_anomaly_", timestamp, clean_prefix, ".png")

      # Karena kita menggunakan panel 1x2, rasio 16:9 (3840x2160) ini akan sangat cantik!
      grDevices::dev.copy(grDevices::png, filename = filename, width = 3840, height = 2160, res = 300)
      grDevices::dev.off()

      message(paste("High-resolution 4K plot saved as:", filename))
    }
  }

  # 5. Return List
  summary_stats <- c(
    Total_Records = nrow(data),
    Evaluated = sum(comp_cases),
    Normal = sum(flags == "Normal", na.rm = TRUE),
    Anomalies = length(anomaly_idx)
  )

  result_list <- list(
    data_scored = data_scored,
    anomaly_indices = anomaly_idx,
    summary = summary_stats
  )

  return(result_list)
}
