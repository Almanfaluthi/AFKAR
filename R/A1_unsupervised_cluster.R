#' @title AFKAR Unsupervised Machine Learning for Clinical Clustering
#'
#' @description
#' Performs unsupervised machine learning (K-Means and Hierarchical Clustering)
#' strictly using Base R. Evaluates cluster tendency, determines optimal k via
#' the Elbow method, and generates two sequential 4K-ready 1x2 dashboards:
#' 1) Cluster Tendency Profile (Elbow Curve & Dendrogram).
#' 2) Clinical Phenotype Profile (PCA Scatter Space & Feature Z-Score Barplot).
#'
#' @param data A data.frame containing the dataset.
#' @param cols A character vector of column names to be used for clustering (must be numeric).
#' @param k_max Integer. Maximum number of clusters to evaluate (default: 10).
#' @param force_k Integer or NULL. Overrides automatic k detection to force a specific number of clusters.
#' @param colors Character vector of colors for plotting clusters.
#' @param reverse_palette Logical. If TRUE, reverses the order of the color palette.
#' @param plot_result Logical. If TRUE, generates the dual 1x2 panel dashboards.
#' @param save_plot Logical. If TRUE, saves the plots as two high-resolution 4K PNGs.
#' @param save_prefix Character. A prefix to add to the saved plot filenames.
#' @param footer Logical. If TRUE, adds a reproducible generation footer to the plots.
#' @param cex Numeric. Text and symbol scaling factor.
#' @param font Numeric. Font style (1 = plain, 2 = bold, 3 = italic, 4 = bold italic).
#' @param family Character. Font family to use for plotting (e.g., "sans", "serif").
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{data_clustered}: The original dataset with a new 'Cluster' column.
#'   \item \code{cluster_summary}: A data.frame showing the mean of each variable per cluster.
#'   \item \code{optimal_k}: The calculated or forced optimal number of clusters.
#'   \item \code{hopkins_stat}: The estimated Hopkins statistic for cluster tendency.
#' }
#'
#' @importFrom stats dist runif kmeans hclust cutree prcomp aggregate complete.cases rect.hclust
#' @importFrom graphics par plot lines points mtext polygon legend barplot axis
#' @importFrom grDevices png dev.copy dev.off adjustcolor chull
#'
#' @export
#'
#' @examples
#' # Run Unsupervised Clustering
#' res_cluster <- A1_unsupervised_cluster(
#'   data = df_A1,
#'   cols = c("Age", "BMI", "Fasting_Glucose", "Total_Cholesterol"),
#'   k_max = 8,
#'   plot_result = TRUE,
#'   save_plot = FALSE
#' )
#'
#' # View clinical cluster summary
#' print(res_cluster$cluster_summary)

A1_unsupervised_cluster <- function(data,
                                    cols,
                                    k_max = 10,
                                    force_k = NULL,
                                    colors = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33"),
                                    reverse_palette = FALSE,
                                    plot_result = TRUE,
                                    save_plot = FALSE,
                                    save_prefix = "digimed",
                                    footer = TRUE,
                                    cex = 1,
                                    font = 1,
                                    family = "sans") {

  # 1. ROBUST ERROR HANDLING & VALIDATION
  stopifnot("Error: 'data' must be a data.frame." = is.data.frame(data))
  stopifnot("Error: 'cols' must be a character vector." = is.character(cols))

  missing_cols <- setdiff(cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste("Error: Missing columns in data:", paste(missing_cols, collapse = ", ")))
  }

  df_ml <- data[, cols, drop = FALSE]

  for (col in cols) {
    if (!is.numeric(df_ml[[col]])) {
      stop(paste("Error: Column", col, "is not numeric. Clustering requires continuous features."))
    }
  }

  comp_cases <- stats::complete.cases(df_ml)
  if (sum(!comp_cases) > 0) {
    warning(paste("Medical Insight:", sum(!comp_cases), "rows with missing values were omitted."))
  }

  df_ml <- df_ml[comp_cases, , drop = FALSE]
  original_data <- data[comp_cases, , drop = FALSE]

  # 2. DATA SCALING
  df_scaled <- scale(df_ml)
  n_rows <- nrow(df_scaled)
  n_cols <- ncol(df_scaled)

  if (n_rows < 5) stop("Error: Dataset is too small for meaningful clinical clustering.")
  actual_k_max <- min(k_max, n_rows - 1)

  # 3. CLUSTER TENDENCY (Hopkins Statistic Estimation)
  m <- min(50, floor(n_rows * 0.1))
  if (m < 2) m <- 2
  sample_idx <- sample(1:n_rows, m)

  rand_data <- matrix(
    stats::runif(m * n_cols,
                 min = apply(df_scaled, 2, min),
                 max = apply(df_scaled, 2, max)),
    nrow = m
  )

  dist_matrix <- as.matrix(stats::dist(rbind(df_scaled, rand_data)))

  W <- numeric(m)
  U <- numeric(m)
  for (i in 1:m) {
    real_idx <- sample_idx[i]
    dists_real <- dist_matrix[real_idx, 1:n_rows][-real_idx]
    W[i] <- min(dists_real)

    rand_idx <- n_rows + i
    dists_rand <- dist_matrix[rand_idx, 1:n_rows]
    U[i] <- min(dists_rand)
  }
  hopkins_stat <- sum(U) / (sum(U) + sum(W))

  # 4. OPTIMAL K DETERMINATION (Elbow Method via WSS)
  wss <- numeric(actual_k_max)
  for (i in 1:actual_k_max) {
    set.seed(123)
    km_temp <- stats::kmeans(df_scaled, centers = i, nstart = 25)
    wss[i] <- km_temp$tot.withinss
  }

  if (is.null(force_k)) {
    x1 <- 1; y1 <- wss[1]
    x2 <- actual_k_max; y2 <- wss[actual_k_max]

    distances <- numeric(actual_k_max)
    for (i in 1:actual_k_max) {
      x0 <- i; y0 <- wss[i]
      numerator <- abs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1)
      denominator <- sqrt((y2 - y1)^2 + (x2 - x1)^2)
      distances[i] <- numerator / denominator
    }
    optimal_k <- which.max(distances)
  } else {
    optimal_k <- force_k
  }

  # 5. FINAL CLUSTERING EXECUTION
  set.seed(123)
  final_km <- stats::kmeans(df_scaled, centers = optimal_k, nstart = 25)
  clusters <- final_km$cluster

  cluster_summary <- stats::aggregate(df_ml, by = list(Cluster = clusters), FUN = mean)
  original_data[["Cluster"]] <- as.factor(clusters)

  # 6. STRICT 4K VISUALIZATION (DUAL DASHBOARDS)
  if (plot_result) {
    if (reverse_palette) colors <- rev(colors)

    # Ensure enough colors
    pal_len <- length(colors)
    color_pal <- if(optimal_k > pal_len) grDevices::colorRampPalette(colors)(optimal_k) else colors[1:optimal_k]

    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))

    footer_text <- paste0("Generated by R-Studio (", R.version.string, ") on ", format(Sys.time(), "%B %d, %Y at %H:%M:%S"))
    algo_meta <- paste0("Base R K-Means & H-Clust | Total N: ", n_rows, " | Hopkins Stat: ", round(hopkins_stat, 3))

    # =========================================================================
    # GAMBAR 1: CLUSTER TENDENCY PROFILE (ELBOW & DENDROGRAM)
    # =========================================================================
    graphics::par(mfrow = c(1, 2), oma = c(4, 0, 4, 0), family = family, font = font, cex = cex)

    # Panel 1.1: WSS Elbow
    graphics::par(mar = c(5, 5, 2, 2) + 0.1)
    graphics::plot(1:actual_k_max, wss, type = "b", pch = 19, col = "#2C3E50", lwd = 2,
                   xlab = "Number of Clusters (k)", ylab = "Total Within Sum of Squares",
                   main = "Mathematical Optimal K (Elbow Method)", las = 1)
    graphics::points(optimal_k, wss[optimal_k], col = "#E74C3C", pch = 19, cex = 2)
    graphics::lines(c(optimal_k, optimal_k), c(0, wss[optimal_k]), col = "#E74C3C", lty = 2, lwd = 2)
    graphics::legend("topright", legend = paste("Optimal k =", optimal_k), col = "#E74C3C", pch = 19, bty = "n")

    # Panel 1.2: Dendrogram
    dist_hc <- stats::dist(df_scaled, method = "euclidean")
    hc <- stats::hclust(dist_hc, method = "ward.D2")

    graphics::par(mar = c(5, 4, 2, 2) + 0.1)
    graphics::plot(hc, hang = -1, labels = FALSE, main = "Hierarchical Cluster Architecture",
                   xlab = "Patients / Observations", ylab = "Euclidean Distance", sub = "")
    stats::rect.hclust(hc, k = optimal_k, border = color_pal)

    header_1 <- "Unsupervised Cluster Tendency Profile"
    graphics::mtext(header_1, outer = TRUE, side = 3, line = 2, cex = 1.3, font = 2)
    graphics::mtext(algo_meta, outer = TRUE, side = 3, line = 0.5, cex = 0.9, font = 1, col = "darkblue")
    if (footer) graphics::mtext(footer_text, side = 1, line = 2, outer = TRUE, adj = 0.5, cex = 0.8, col = "dimgray")

    if (save_plot) {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      clean_prefix <- ifelse(nchar(save_prefix) > 0, paste0("_", save_prefix), "")
      filename1 <- paste0("A1_cluster_", timestamp, clean_prefix, "_Part1_Tendency.png")
      grDevices::dev.copy(grDevices::png, filename = filename1, width = 3840, height = 2160, res = 300)
      grDevices::dev.off()
      message(paste("High-resolution 4K plot (Part 1) saved as:", filename1))
    }

    # =========================================================================
    # GAMBAR 2: CLINICAL PHENOTYPE PROFILE (PCA HULLS & FEATURE Z-SCORES)
    # =========================================================================
    graphics::par(mfrow = c(1, 2), oma = c(4, 0, 4, 0), family = family, font = font, cex = cex)

    # Panel 2.1: PCA Scatter Space with Convex Hulls
    pca_res <- stats::prcomp(df_scaled, scale. = FALSE)
    pca_data <- pca_res$x[, 1:2]
    var_exp <- round((pca_res$sdev^2 / sum(pca_res$sdev^2))[1:2] * 100, 1)

    graphics::par(mar = c(5, 5, 2, 2) + 0.1)
    graphics::plot(pca_data, col = color_pal[clusters], pch = 21, bg = color_pal[clusters],
                   xlab = paste0("Principal Component 1 (", var_exp[1], "%)"),
                   ylab = paste0("Principal Component 2 (", var_exp[2], "%)"),
                   main = "Clinical Phenotype Space (PCA)", las = 1)

    for (i in 1:optimal_k) {
      pts <- pca_data[clusters == i, , drop = FALSE]
      if (nrow(pts) >= 3) {
        hull_idx <- grDevices::chull(pts)
        hull_idx <- c(hull_idx, hull_idx[1])
        graphics::polygon(pts[hull_idx, ], border = color_pal[i], lwd = 2,
                          col = grDevices::adjustcolor(color_pal[i], alpha.f = 0.2))
      }
    }

    graphics::legend("topright", legend = paste("Cluster", 1:optimal_k),
                     fill = grDevices::adjustcolor(color_pal, alpha.f = 0.5),
                     border = color_pal, bty = "n", cex = 0.9)

    # Panel 2.2: Clinical Feature Profile (Z-Score Barplot)
    # Plotting the Z-scores (scaled centers) to see what characterizes each cluster
    graphics::par(mar = c(5, 7, 2, 2) + 0.1) # Extra left margin for variable names

    z_centers <- final_km$centers
    rownames(z_centers) <- paste("Cluster", 1:optimal_k)

    # Transpose for grouped barplot: Features on X/Y axis, grouped by clusters
    bp <- graphics::barplot(z_centers, beside = TRUE, horiz = TRUE, las = 1,
                            col = color_pal, border = "white",
                            xlab = "Standardized Mean (Z-Score)",
                            main = "Phenotype Characteristics per Cluster")

    graphics::abline(v = 0, lty = 2, col = "gray40", lwd = 2)
    graphics::legend("bottomright", legend = paste("Cluster", 1:optimal_k),
                     fill = color_pal, border = "white", bty = "n", cex = 0.8)

    header_3 <- "Clinical Phenotype Discovery Profile"
    graphics::mtext(header_3, outer = TRUE, side = 3, line = 2, cex = 1.3, font = 2)
    graphics::mtext(algo_meta, outer = TRUE, side = 3, line = 0.5, cex = 0.9, font = 1, col = "darkblue")
    if (footer) graphics::mtext(footer_text, side = 1, line = 2, outer = TRUE, adj = 0.5, cex = 0.8, col = "dimgray")

    if (save_plot) {
      filename2 <- paste0("A1_cluster_", timestamp, clean_prefix, "_Part2_Phenotypes.png")
      grDevices::dev.copy(grDevices::png, filename = filename2, width = 3840, height = 2160, res = 300)
      grDevices::dev.off()
      message(paste("High-resolution 4K plot (Part 2) saved as:", filename2))
    }
  }

  return(list(
    data_clustered = original_data,
    cluster_summary = cluster_summary,
    optimal_k = optimal_k,
    hopkins_stat = hopkins_stat
  ))
}
