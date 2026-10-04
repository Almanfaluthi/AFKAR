#' @title AFKAR Unsupervised Machine Learning for Clinical Clustering
#' @description
#' Performs unsupervised machine learning (K-Means and Hierarchical Clustering)
#' strictly using Base R. It automatically evaluates cluster tendency,
#' determines the optimal number of clusters using the mathematical Elbow method,
#' and generates beautiful, CRAN-ready visualizations (WSS Curve, Dendrogram,
#' and PCA-based Cluster Plot with Convex Hulls).
#'
#' @param data A data.frame containing the dataset.
#' @param cols A character vector of column names to be used for clustering.
#' Must be numeric columns.
#' @param k_max Integer. The maximum number of clusters to evaluate for the
#' Elbow method (default: 10).
#' @param force_k Integer or NULL. If provided, overrides the automatic optimal
#' k detection and forces the algorithm to use this number of clusters.
#' @param palette Character vector of colors for plotting.
#' @param reverse_palette Logical. If TRUE, reverses the order of the color palette.
#' @param footer Logical. If TRUE, adds a reproducibility footer to the plots.
#' @param save_plot Logical. If TRUE, saves the plots as 4K resolution PNGs.
#' @param save_prefix Character. A prefix to add to the saved plot filenames.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{data_clustered}: The original dataset with a new 'cluster' column.
#'   \item \code{cluster_summary}: A data.frame showing the mean of each variable per cluster.
#'   \item \code{optimal_k}: The calculated or forced optimal number of clusters.
#'   \item \code{hopkins_stat}: The estimated Hopkins statistic for cluster tendency.
#' }
#'
#' @importFrom stats dist runif kmeans hclust cutree prcomp aggregate complete.cases rect.hclust
#' @importFrom graphics par plot lines points mtext polygon legend
#' @importFrom grDevices dev.copy dev.off png adjustcolor chull
#'
#' @export
#'
#' @examples
#' # Use built-in iris dataset for clinical-like simulation
#' df <- iris
#' # Run unsupervised clustering on numeric columns
#' res <- A1_unsupervised_cluster(
#'   data = df_A1,
#'   cols = c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"),
#'   save_plot = TRUE
#' )
#' head(res$data_clustered)
#' print(res$cluster_summary)
A1_unsupervised_cluster <- function(data,
                                    cols,
                                    k_max = 10,
                                    force_k = NULL,
                                    palette = c("#E41A1C", "#377EB8", "#4DAF4A",
                                                "#984EA3", "#FF7F00", "#FFFF33"),
                                    reverse_palette = FALSE,
                                    footer = TRUE,
                                    save_plot = TRUE,
                                    save_prefix = "digimed") {

  # 1. ROBUST ERROR HANDLING & VALIDATION
  stopifnot("Input 'data' must be a data.frame." = is.data.frame(data))
  stopifnot("Input 'cols' must be a character vector of column names." = is.character(cols))

  missing_cols <- setdiff(cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste("The following columns are not found in the data:", paste(missing_cols, collapse = ", ")))
  }

  # Extract and clean data
  df_ml <- data[, cols, drop = FALSE]

  # Ensure all selected columns are numeric
  for (col in cols) {
    if (!is.numeric(df_ml[[col]])) {
      stop(paste("Column", col, "is not numeric. Clustering requires numeric features."))
    }
  }

  # Remove NAs strictly for the clustering subset
  complete_cases_idx <- stats::complete.cases(df_ml)
  if (sum(!complete_cases_idx) > 0) {
    warning(paste(sum(!complete_cases_idx), "rows with missing values were removed prior to modeling."))
  }

  df_ml <- df_ml[complete_cases_idx, , drop = FALSE]
  original_data <- data[complete_cases_idx, , drop = FALSE]

  # 2. DATA SCALING
  df_scaled <- scale(df_ml)
  n_rows <- nrow(df_scaled)
  n_cols <- ncol(df_scaled)

  # Prevent errors on tiny datasets
  if (n_rows < 3) stop("Dataset is too small for clustering.")
  actual_k_max <- min(k_max, n_rows - 1)

  # 3. CLUSTER TENDENCY (Custom Hopkins Statistic Estimation)
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
    # Nearest real point to sampled real point
    dists_real <- dist_matrix[real_idx, 1:n_rows][-real_idx]
    W[i] <- min(dists_real)

    # Nearest real point to random point
    rand_idx <- n_rows + i
    dists_rand <- dist_matrix[rand_idx, 1:n_rows]
    U[i] <- min(dists_rand)
  }
  hopkins_stat <- sum(U) / (sum(U) + sum(W))

  # 4. OPTIMAL K DETERMINATION (Mathematical Elbow Method via WSS)
  wss <- numeric(actual_k_max)
  for (i in 1:actual_k_max) {
    set.seed(123) # For reproducibility
    km_temp <- stats::kmeans(df_scaled, centers = i, nstart = 25)
    wss[i] <- km_temp$tot.withinss
  }

  # Calculate perpendicular distance to find the elbow if force_k is NULL
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

  # Aggregation (Summary)
  cluster_summary <- stats::aggregate(df_ml, by = list(Cluster = clusters), FUN = mean)

  # Bind cluster to original dataset
  original_data[["cluster"]] <- as.factor(clusters)

  # 6. STRICT VISUALIZATION RULES (AFKAR Guidelines)
  if (reverse_palette) palette <- rev(palette)

  # Save current par parameters to restore later
  old_par <- graphics::par(no.readonly = TRUE)
  timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")

  footer_text <- paste0("Generated by R-Studio (", R.version.string,
                        ") on ", format(Sys.time(), "%B %d, %Y at %H:%M:%S"))

  # Helper to plot footer
  add_footer <- function() {
    if (footer) {
      graphics::mtext(footer_text, side = 1, line = 4, adj = 0.5, cex = 0.75, col = "darkgray")
    }
  }

  # --- PLOT 1: WSS / ELBOW ---
  graphics::par(mar = c(6, 4, 4, 2) + 0.1)
  graphics::plot(1:actual_k_max, wss, type = "b", pch = 19, col = "#2C3E50",
                 xlab = "Number of Clusters (k)", ylab = "Total Within Sum of Squares",
                 main = "Optimal K Determination (Elbow Method)",
                 las = 1, cex.main = 1.2, cex.lab = 1.1)
  graphics::points(optimal_k, wss[optimal_k], col = "#E74C3C", pch = 19, cex = 2)
  graphics::lines(c(optimal_k, optimal_k), c(0, wss[optimal_k]), col = "#E74C3C", lty = 2)
  add_footer()

  if (save_plot) {
    filename_1 <- paste0(save_prefix, "_", timestamp_str, "_1_elbow.png")
    grDevices::dev.copy(grDevices::png, filename = filename_1, width = 3840, height = 2160, res = 300)
    grDevices::dev.off()
  }

  # --- PLOT 2: DENDROGRAM ---
  dist_hc <- stats::dist(df_scaled, method = "euclidean")
  hc <- stats::hclust(dist_hc, method = "ward.D2")

  graphics::par(mar = c(6, 4, 4, 2) + 0.1)
  graphics::plot(hc, hang = -1, cex = 0.6, main = "Hierarchical Clustering Dendrogram",
                 xlab = "Observations", ylab = "Height", sub = "")
  rect_colors <- rep(palette, length.out = optimal_k)
  stats::rect.hclust(hc, k = optimal_k, border = rect_colors)
  add_footer()

  if (save_plot) {
    filename_2 <- paste0(save_prefix, "_", timestamp_str, "_2_dendro.png")
    grDevices::dev.copy(grDevices::png, filename = filename_2, width = 3840, height = 2160, res = 300)
    grDevices::dev.off()
  }

  # --- PLOT 3: PCA CLUSTER PLOT ---
  pca_res <- stats::prcomp(df_scaled, scale. = FALSE)
  pca_data <- pca_res$x[, 1:2]

  var_explained <- round((pca_res$sdev^2 / sum(pca_res$sdev^2))[1:2] * 100, 1)

  graphics::par(mar = c(6, 4, 4, 2) + 0.1)
  graphics::plot(pca_data, col = palette[clusters], pch = 19,
                 xlab = paste0("Principal Component 1 (", var_explained[1], "%)"),
                 ylab = paste0("Principal Component 2 (", var_explained[2], "%)"),
                 main = "Cluster Plot (PCA Dimensionality Reduction)",
                 las = 1)

  # Draw Elegant Convex Hulls
  for (i in 1:optimal_k) {
    pts <- pca_data[clusters == i, , drop = FALSE]
    if (nrow(pts) >= 3) {
      hull_idx <- grDevices::chull(pts)
      hull_idx <- c(hull_idx, hull_idx[1]) # close polygon
      graphics::polygon(pts[hull_idx, ], border = palette[i],
                        col = grDevices::adjustcolor(palette[i], alpha.f = 0.2))
    }
  }

  graphics::legend("topright", legend = paste("Cluster", 1:optimal_k),
                   fill = grDevices::adjustcolor(palette[1:optimal_k], alpha.f = 0.5),
                   border = palette[1:optimal_k], bty = "n", cex = 0.9)
  add_footer()

  if (save_plot) {
    filename_3 <- paste0(save_prefix, "_", timestamp_str, "_3_pca.png")
    grDevices::dev.copy(grDevices::png, filename = filename_3, width = 3840, height = 2160, res = 300)
    grDevices::dev.off()
  }

  # Restore graphics parameters
  graphics::par(old_par)

  # 7. RETURN OBJECTS
  message(sprintf("[AFKAR ML] Unsupervised clustering completed successfully."))
  message(sprintf("Estimated Hopkins Statistic: %.3f", hopkins_stat))
  message(sprintf("Optimal Clusters (k): %d", optimal_k))

  return(list(
    data_clustered = original_data,
    cluster_summary = cluster_summary,
    optimal_k = optimal_k,
    hopkins_stat = hopkins_stat
  ))
}
