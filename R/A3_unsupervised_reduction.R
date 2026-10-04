#' @title Dimensionality Reduction (Unsupervised: PCA)
#'
#' @description
#' Performs Principal Component Analysis (PCA) to reduce high-dimensional
#' clinical data, surveys (e.g., PSQI, MBI), or biomarker panels into
#' a lower-dimensional space. Includes dual visualizations: a Scree Plot
#' showing variance explained, and a 2D Scatter Plot of PC1 vs PC2.
#'
#' @param data A data.frame containing the clinical dataset.
#' @param features Character vector specifying the numeric columns to reduce.
#' @param scale Logical. If TRUE (default), scales variables to have unit variance before PCA.
#' @param plot_result Logical. If TRUE, generates a 1x2 panel plot.
#' @param colors Character vector of length 2 for plot colors c("Line/Point", "Highlight").
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
#'   \item \code{pca_model}: The complete prcomp object.
#'   \item \code{data_reduced}: Original data appended with PC1 and PC2 scores.
#'   \item \code{variance_summary}: A matrix showing variance explained by each principal component.
#' }
#'
#' @importFrom stats prcomp complete.cases
#' @importFrom graphics par plot lines points mtext abline
#' @importFrom grDevices png dev.copy dev.off
#'
#' @export
#'
#' @examples
#' # Run Dimensionality Reduction
#' res_pca <- A3_unsupervised_reduction(
#'   data = df_A3,
#'   features = paste0("Q", 1:10),
#'   scale = TRUE,
#'   plot_result = TRUE,
#'   save_plot = FALSE
#' )
#'
#' # View the variance summary
#' print(res_pca$variance_summary[, 1:3]) # First 3 PCs

A3_unsupervised_reduction <- function(data,
                                      features,
                                      scale = TRUE,
                                      plot_result = TRUE,
                                      colors = c("#1f77b4", "#ff7f0e"),
                                      reverse_palette = FALSE,
                                      save_plot = FALSE,
                                      save_prefix = "",
                                      footer = TRUE,
                                      cex = 1,
                                      font = 1,
                                      family = "sans") {

  # 1. Input Validation
  stopifnot("Error: 'data' must be a data.frame." = is.data.frame(data))
  stopifnot("Error: 'features' must be a character vector." = is.character(features))

  missing_cols <- setdiff(features, names(data))
  if (length(missing_cols) > 0) {
    stop(paste("Error: Features not found in data:", paste(missing_cols, collapse = ", ")))
  }

  for (f in features) {
    if (!is.numeric(data[[f]])) {
      stop(paste("Error: Feature '", f, "' must be numeric for PCA.", sep = ""))
    }
  }

  if (length(features) < 3) {
    warning("Medical Insight: PCA is most useful for 3 or more variables. For 2 variables, a simple scatter plot usually suffices.")
  }

  # 2. Data Preparation
  df_subset <- data[, features, drop = FALSE]
  comp_cases <- stats::complete.cases(df_subset)

  if (sum(comp_cases) < 3) {
    stop("Error: Insufficient complete observations. Clean your dataset first.")
  }
  df_clean <- df_subset[comp_cases, , drop = FALSE]

  # 3. Core Algorithm: Principal Component Analysis (PCA)
  pca_res <- stats::prcomp(df_clean, center = TRUE, scale. = scale)

  # Extract Variance Information
  pca_summary <- summary(pca_res)
  var_explained <- pca_summary$importance[2, ] * 100 # Proportion of Variance
  cum_var_explained <- pca_summary$importance[3, ] * 100 # Cumulative Proportion

  # Merge PC1 and PC2 back into the original dataset
  data_reduced <- data
  data_reduced[["PC1"]] <- NA
  data_reduced[["PC2"]] <- NA

  if (ncol(pca_res$x) >= 2) {
    data_reduced[comp_cases, "PC1"] <- pca_res$x[, 1]
    data_reduced[comp_cases, "PC2"] <- pca_res$x[, 2]
  } else {
    data_reduced[comp_cases, "PC1"] <- pca_res$x[, 1]
    warning("Only 1 Principal Component was generated.")
  }

  # 4. Strict Visualization (1x2 Panel Plot)
  if (plot_result && ncol(pca_res$x) >= 2) {

    if (reverse_palette) {
      colors <- rev(colors)
    }

    # Save old par
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))

    # KUNCI: mfrow 1x2, mar untuk dalam, oma untuk footer luar (mencegah overlap)
    graphics::par(mfrow = c(1, 2),
                  mar = c(5, 5, 4, 2) + 0.1,
                  oma = c(4, 0, 0, 0),
                  family = family, font = font, cex = cex)

    # ================= PLOT 1: SCREE PLOT (Variance Explained) =================
    num_pcs_to_plot <- min(10, length(var_explained)) # Plot max top 10 PCs
    x_coords <- 1:num_pcs_to_plot

    graphics::plot(x_coords, var_explained[1:num_pcs_to_plot],
                   type = "b",
                   pch = 19,
                   col = colors[1],
                   lwd = 2,
                   xaxt = "n",
                   ylim = c(0, max(var_explained) + 10),
                   xlab = "Principal Component (PC)",
                   ylab = "Variance Explained (%)",
                   main = "Scree Plot: Dimensionality Process")

    graphics::axis(1, at = x_coords, labels = paste0("PC", x_coords))

    # Add cumulative variance text above the points
    for (i in x_coords) {
      graphics::text(i, var_explained[i] + 4,
                     labels = paste0(round(cum_var_explained[i], 1), "%"),
                     col = "darkgray", cex = 0.8)
    }

    # ================= PLOT 2: PCA SCATTER SPACE =================
    pc1_var <- round(var_explained[1], 1)
    pc2_var <- round(var_explained[2], 1)

    graphics::plot(pca_res$x[, 1], pca_res$x[, 2],
                   pch = 21,
                   bg = colors[2], # Fill color for points
                   col = "white",  # Border color for points
                   cex = 1.5,
                   xlab = paste0("PC1 (", pc1_var, "% Variance)"),
                   ylab = paste0("PC2 (", pc2_var, "% Variance)"),
                   main = "2D Patient/Record Mapping")

    # Add subtle axis lines through the origin (0,0)
    graphics::abline(h = 0, v = 0, col = "gray70", lty = 3, lwd = 1.5)

    # ================= GLOBAL FOOTER =================
    if (footer) {
      footer_text <- paste0("Generated by R-Studio (", R.version.string, ") on ", format(Sys.time(), "%B %d, %Y at %H:%M:%S"))
      graphics::mtext(footer_text, side = 1, line = 2, outer = TRUE, adj = 0.5, cex = 0.8, col = "dimgray")
    }

    # ================= HIGH-RES 4K EXPORT =================
    if (save_plot) {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      clean_prefix <- ifelse(nchar(save_prefix) > 0, paste0("_", save_prefix), "")
      filename <- paste0("A3_unsupervised_reduction_", timestamp, clean_prefix, ".png")

      grDevices::dev.copy(grDevices::png, filename = filename, width = 3840, height = 2160, res = 300)
      grDevices::dev.off()

      message(paste("High-resolution 4K plot saved as:", filename))
    }
  }

  # 5. Return List
  result_list <- list(
    pca_model = pca_res,
    data_reduced = data_reduced,
    variance_summary = pca_summary$importance
  )

  return(result_list)
}
