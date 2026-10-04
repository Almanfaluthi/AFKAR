#' @title Time-Series Trajectory Clustering (Unsupervised DTW)
#'
#' @description
#' Groups patients based on the shape of their clinical time-series curves
#' (e.g., ICU vital signs or lab values over 7 days) using a native Base R
#' implementation of Dynamic Time Warping (DTW). It calculates the DTW distance
#' matrix and applies hierarchical clustering to find hidden trajectory patterns.
#'
#' @param data A data.frame containing the clinical dataset.
#' @param time_features Character vector of column names representing the sequential time points (e.g., Day_1 to Day_7).
#' @param num_clusters Integer. The number of trajectory clusters to identify (k value).
#' @param plot_result Logical. If TRUE, generates a 1x2 panel plot (Dendrogram & Trajectory Curves).
#' @param colors Character vector specifying the colors for different clusters.
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
#'   \item \code{dtw_matrix}: The symmetric DTW distance matrix.
#'   \item \code{cluster_assignments}: A named vector of cluster assignments for each patient.
#'   \item \code{hclust_obj}: The hierarchical clustering object.
#'   \item \code{data_scored}: Original data appended with the Trajectory_Cluster column.
#' }
#'
#' @importFrom stats complete.cases as.dist hclust cutree rect.hclust
#' @importFrom graphics par plot lines legend mtext
#' @importFrom grDevices colorRampPalette png dev.copy dev.off
#'
#' @export
#'
#' @examples
#' # Run Unsupervised Time-Series Clustering
#' res_ts <- A7_unsupervised_timeseries(
#'   data = df_A7,
#'   time_features = paste0("Day_", 1:5),
#'   num_clusters = 3,
#'   plot_result = TRUE,
#'   save_plot = FALSE
#' )
#'
#' # View the cluster assignment for the first 5 patients
#' print(head(res_ts$data_scored[, c("Patient_ID", "Trajectory_Cluster")]))

A7_unsupervised_timeseries <- function(data,
                                       time_features,
                                       num_clusters = 3,
                                       plot_result = TRUE,
                                       colors = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"),
                                       reverse_palette = FALSE,
                                       save_plot = FALSE,
                                       save_prefix = "",
                                       footer = TRUE,
                                       cex = 1,
                                       font = 1,
                                       family = "sans") {

  # 1. Input Validation
  stopifnot("Error: 'data' must be a data.frame." = is.data.frame(data))
  stopifnot("Error: 'time_features' must be a character vector." = is.character(time_features))

  missing_cols <- setdiff(time_features, names(data))
  if (length(missing_cols) > 0) {
    stop(paste("Error: Time features not found in data:", paste(missing_cols, collapse = ", ")))
  }

  for (f in time_features) {
    if (!is.numeric(data[[f]])) {
      stop(paste("Error: Time feature '", f, "' must be numeric.", sep = ""))
    }
  }

  n_timepoints <- length(time_features)
  if (n_timepoints < 3) {
    stop("Error: Time-series clustering requires at least 3 sequential time points.")
  }

  # 2. Data Preparation
  df_subset <- data[, time_features, drop = FALSE]
  comp_cases <- stats::complete.cases(df_subset)
  df_clean <- df_subset[comp_cases, , drop = FALSE]
  n_cases <- nrow(df_clean)

  if (n_cases < num_clusters) {
    stop("Error: Number of complete cases is less than the requested number of clusters.")
  }

  if (n_cases > 200) {
    message("Medical Insight: N > 200 detected. Native Base R DTW calculation might take a moment. Please wait as it maps the trajectory matrices...")
  }

  # 3. Core Algorithm: Native Dynamic Time Warping (DTW) Helper
  dtw_dist_native <- function(x, y) {
    nx <- length(x)
    ny <- length(y)
    cost_mat <- matrix(Inf, nx + 1, ny + 1)
    cost_mat[1, 1] <- 0

    for (i in seq_len(nx)) {
      for (j in seq_len(ny)) {
        cost <- abs(x[i] - y[j])
        cost_mat[i + 1, j + 1] <- cost + min(cost_mat[i, j + 1],     # Insertion
                                             cost_mat[i + 1, j],     # Deletion
                                             cost_mat[i, j])         # Match
      }
    }
    return(cost_mat[nx + 1, ny + 1])
  }

  # Construct Distance Matrix O(N^2)
  dist_mat <- matrix(0, nrow = n_cases, ncol = n_cases)
  rownames(dist_mat) <- colnames(dist_mat) <- 1:n_cases

  for (i in 1:(n_cases - 1)) {
    for (j in (i + 1):n_cases) {
      d <- dtw_dist_native(as.numeric(df_clean[i, ]), as.numeric(df_clean[j, ]))
      dist_mat[i, j] <- d
      dist_mat[j, i] <- d # Make it symmetric
    }
  }

  # 4. Hierarchical Clustering
  dtw_dist_obj <- stats::as.dist(dist_mat)
  hc_obj <- stats::hclust(dtw_dist_obj, method = "ward.D2")

  clusters <- stats::cutree(hc_obj, k = num_clusters)

  # Merge back to original data
  data_scored <- data
  data_scored[["Trajectory_Cluster"]] <- NA
  data_scored[comp_cases, "Trajectory_Cluster"] <- paste("Cluster", clusters)

  # 5. Strict Visualization (1x2 Panel Plot)
  if (plot_result) {

    if (reverse_palette) {
      colors <- rev(colors)
    }

    color_pal <- grDevices::colorRampPalette(colors)(num_clusters)

    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))

    # Setup global layout with oma for header/footer
    graphics::par(mfrow = c(1, 2),
                  oma = c(4, 0, 3, 0),
                  family = family, font = font, cex = cex)

    # ================= PLOT 1: TRAJECTORY DENDROGRAM =================
    graphics::par(mar = c(3, 4, 2, 2) + 0.1)

    graphics::plot(hc_obj,
                   labels = FALSE, # Hide labels to prevent clutter on large clinical data
                   hang = -1,
                   main = paste("DTW Trajectory Dendrogram (k =", num_clusters, ")"),
                   xlab = "Patients (Clustered by Curve Shape)",
                   sub = "",
                   ylab = "DTW Distance",
                   cex = 0.9)

    stats::rect.hclust(hc_obj, k = num_clusters, border = color_pal)

    # ================= PLOT 2: SPAGHETTI CURVE PLOT =================
    graphics::par(mar = c(4, 4, 2, 2) + 0.1)

    # Setup empty plot
    y_min <- min(df_clean, na.rm = TRUE)
    y_max <- max(df_clean, na.rm = TRUE)
    x_range <- 1:n_timepoints

    graphics::plot(x_range, type = "n",
                   ylim = c(y_min, y_max + (y_max - y_min) * 0.1), # Add headroom for legend
                   xaxt = "n",
                   xlab = "Sequential Time Points",
                   ylab = "Clinical Value",
                   main = "Clinical Trajectory Space")

    graphics::axis(1, at = x_range, labels = time_features, las = 2, cex.axis = 0.8)

    # Draw individual lines (semi-transparent)
    for (i in 1:n_cases) {
      c_idx <- clusters[i]
      # Create hex color with ~40% opacity (alpha = "66")
      trans_col <- paste0(color_pal[c_idx], "66")
      graphics::lines(x_range, as.numeric(df_clean[i, ]), col = trans_col, lwd = 1)
    }

    # Calculate and draw Centroid (Mean) lines for each cluster
    for (k in 1:num_clusters) {
      cluster_data <- df_clean[clusters == k, , drop = FALSE]
      centroid_line <- colMeans(cluster_data, na.rm = TRUE)

      # Thick opaque line for the average trajectory
      graphics::lines(x_range, centroid_line, col = color_pal[k], lwd = 4, lty = 1)

      # Add subtle points to the centroid
      graphics::points(x_range, centroid_line, col = color_pal[k], pch = 21, bg = "white", cex = 1.2, lwd = 2)
    }

    graphics::legend("topleft",
                     legend = paste("Cluster", 1:num_clusters),
                     col = color_pal,
                     lwd = 4,
                     bty = "n",
                     cex = 0.8)

    # ================= GLOBAL TITLES & FOOTER =================
    graphics::mtext("Unsupervised Time-Series Trajectory Clustering", outer = TRUE, side = 3, line = 1, cex = 1.3, font = 2)

    if (footer) {
      footer_text <- paste0("Generated by R-Studio (", R.version.string, ") on ", format(Sys.time(), "%B %d, %Y at %H:%M:%S"))
      graphics::mtext(footer_text, side = 1, line = 2, outer = TRUE, adj = 0.5, cex = 0.8, col = "dimgray")
    }

    # ================= HIGH-RES 4K EXPORT =================
    if (save_plot) {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      clean_prefix <- ifelse(nchar(save_prefix) > 0, paste0("_", save_prefix), "")
      filename <- paste0("A7_unsupervised_timeseries_", timestamp, clean_prefix, ".png")

      grDevices::dev.copy(grDevices::png, filename = filename, width = 3840, height = 2160, res = 300)
      grDevices::dev.off()

      message(paste("High-resolution 4K plot saved as:", filename))
    }
  }

  # 6. Return Data
  result_list <- list(
    dtw_matrix = dist_mat,
    cluster_assignments = clusters,
    hclust_obj = hc_obj,
    data_scored = data_scored
  )

  return(result_list)
}
