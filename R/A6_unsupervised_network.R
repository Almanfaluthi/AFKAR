#' @title Network Community Detection for Syndemic Architecture (Unsupervised)
#'
#' @description
#' Identifies and visualizes disease communities/syndemics using native Base R.
#' It calculates the Jaccard similarity index for co-occurrences of clinical
#' conditions (e.g., ICD-10 codes) and applies hierarchical clustering to
#' detect hidden communities. Features a dual-panel visualization: a Community
#' Dendrogram and a Circular Syndemic Web Graph.
#'
#' @param data A data.frame of binary numeric columns (1 = Present, 0 = Absent) representing patients and diseases/conditions.
#' @param features Character vector specifying the column names (diseases) to include in the network.
#' @param min_weight Numeric. Minimum Jaccard similarity (0 to 1) required to draw a connecting line (edge) between two diseases.
#' @param num_communities Integer. The number of disease communities to detect (k value for clustering).
#' @param plot_result Logical. If TRUE, generates a 1x2 panel plot.
#' @param colors Character vector specifying the colors for different communities. Must be at least length 2.
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
#'   \item \code{jaccard_matrix}: The symmetric Jaccard similarity matrix of diseases.
#'   \item \code{communities}: A named integer vector showing the assigned community for each disease.
#'   \item \code{hclust_obj}: The hierarchical clustering object.
#' }
#'
#' @importFrom stats complete.cases as.dist hclust cutree rect.hclust
#' @importFrom graphics par plot segments points text mtext
#' @importFrom grDevices colorRampPalette png dev.copy dev.off
#'
#' @export
#'
#' @examples
#' # Run Unsupervised Network Detection
#' dis_network <- A6_unsupervised_network(
#'   data = df_A6,
#'   features = names(df_A6),
#'   min_weight = 0.15,
#'   num_communities = 2,
#'   plot_result = TRUE,
#'   save_plot = FALSE
#' )
#'
#' # Print discovered communities
#' print(dis_network$communities)

A6_unsupervised_network <- function(data,
                                    features,
                                    min_weight = 0.1,
                                    num_communities = 3,
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
  stopifnot("Error: 'features' must be a character vector." = is.character(features))

  missing_cols <- setdiff(features, names(data))
  if (length(missing_cols) > 0) {
    stop(paste("Error: Features not found in data:", paste(missing_cols, collapse = ", ")))
  }

  for (f in features) {
    if (!is.numeric(data[[f]]) || !all(data[[f]] %in% c(0, 1, NA))) {
      stop(paste("Error: Feature '", f, "' must be binary numeric (0 and 1) for Jaccard indexing.", sep = ""))
    }
  }

  n_feat <- length(features)
  if (n_feat < 3) {
    stop("Error: Network architecture requires at least 3 clinical features.")
  }

  if (num_communities > n_feat - 1) {
    stop("Error: 'num_communities' must be less than the total number of features.")
  }

  # 2. Data Preparation
  df_subset <- data[, features, drop = FALSE]
  comp_cases <- stats::complete.cases(df_subset)
  df_clean <- df_subset[comp_cases, , drop = FALSE]

  if (nrow(df_clean) < 10) {
    warning("Medical Insight: Dataset is extremely small. Network linkages might be statistically unstable.")
  }

  # 3. Core Algorithm: Jaccard Similarity Matrix
  # Custom native implementation to ensure zero-dependency
  jaccard_mat <- matrix(0, nrow = n_feat, ncol = n_feat)
  rownames(jaccard_mat) <- colnames(jaccard_mat) <- features

  for (i in 1:n_feat) {
    for (j in 1:n_feat) {
      if (i == j) {
        jaccard_mat[i, j] <- 1
      } else {
        vec_i <- df_clean[[features[i]]]
        vec_j <- df_clean[[features[j]]]

        intersection_count <- sum(vec_i == 1 & vec_j == 1)
        union_count <- sum(vec_i == 1 | vec_j == 1)

        if (union_count == 0) {
          jaccard_mat[i, j] <- 0
        } else {
          jaccard_mat[i, j] <- intersection_count / union_count
        }
      }
    }
  }

  # 4. Community Detection (Hierarchical Clustering)
  # Convert similarity to distance (1 - similarity)
  dist_mat <- stats::as.dist(1 - jaccard_mat)
  hc_obj <- stats::hclust(dist_mat, method = "ward.D2")

  # Extract communities
  community_vector <- stats::cutree(hc_obj, k = num_communities)

  # Prepare Colors
  if (reverse_palette) {
    colors <- rev(colors)
  }

  color_pal <- grDevices::colorRampPalette(colors)(num_communities)
  node_colors <- color_pal[community_vector]

  # 5. Strict Visualization (1x2 Panel Plot)
  if (plot_result) {

    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))

    # Global layout: oma (Outer Margin Area) for global footer
    graphics::par(mfrow = c(1, 2),
                  oma = c(4, 0, 3, 0),
                  family = family, font = font, cex = cex)

    # ================= PLOT 1: SYNDEMIC DENDROGRAM =================
    graphics::par(mar = c(3, 4, 2, 2) + 0.1)

    graphics::plot(hc_obj,
                   hang = -1,
                   main = paste("Disease Communities (k =", num_communities, ")"),
                   xlab = "",
                   sub = "",
                   ylab = "Jaccard Distance",
                   cex = 0.9)

    # Strict namespace for rect.hclust
    stats::rect.hclust(hc_obj, k = num_communities, border = color_pal)

    # ================= PLOT 2: CIRCULAR NETWORK WEB =================
    graphics::par(mar = c(1, 1, 2, 1) + 0.1)

    # Native Circle Layout Mathematics (Trigonometry)
    theta <- seq(0, 2 * pi, length.out = n_feat + 1)[-(n_feat + 1)]
    # Align nodes by reordering theta based on clustering order to prevent edge crossing clutter
    ord <- hc_obj$order
    x_coords <- numeric(n_feat)
    y_coords <- numeric(n_feat)

    for (i in 1:n_feat) {
      idx <- which(ord == i)
      x_coords[i] <- cos(theta[idx])
      y_coords[i] <- sin(theta[idx])
    }

    # Empty Plot Canvas
    graphics::plot(1, 1, type = "n",
                   xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5),
                   axes = FALSE, xlab = "", ylab = "",
                   main = "Syndemic Web Architecture")

    # Draw Edges (Lines)
    for (i in 1:(n_feat - 1)) {
      for (j in (i + 1):n_feat) {
        weight <- jaccard_mat[i, j]
        if (weight >= min_weight) {
          # Dynamic line thickness based on Jaccard strength
          lwd_scaled <- weight * 10
          graphics::segments(x0 = x_coords[i], y0 = y_coords[i],
                             x1 = x_coords[j], y1 = y_coords[j],
                             col = "#d3d3d380", # Semi-transparent light gray
                             lwd = lwd_scaled)
        }
      }
    }

    # Draw Nodes (Points)
    graphics::points(x_coords, y_coords,
                     pch = 21,
                     bg = node_colors,
                     col = "white", # Clean border
                     cex = 4)

    # Draw Text Labels radially expanded
    label_multiplier <- 1.25
    graphics::text(x = x_coords * label_multiplier,
                   y = y_coords * label_multiplier,
                   labels = features,
                   cex = 0.8,
                   font = 2,
                   col = "gray10")

    # ================= GLOBAL TITLES & FOOTER =================
    graphics::mtext("Unsupervised Network Community Detection", outer = TRUE, side = 3, line = 1, cex = 1.3, font = 2)

    if (footer) {
      footer_text <- paste0("Generated by R-Studio (", R.version.string, ") on ", format(Sys.time(), "%B %d, %Y at %H:%M:%S"))
      graphics::mtext(footer_text, side = 1, line = 2, outer = TRUE, adj = 0.5, cex = 0.8, col = "dimgray")
    }

    # ================= HIGH-RES 4K EXPORT =================
    if (save_plot) {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      clean_prefix <- ifelse(nchar(save_prefix) > 0, paste0("_", save_prefix), "")
      filename <- paste0("A6_unsupervised_network_", timestamp, clean_prefix, ".png")

      grDevices::dev.copy(grDevices::png, filename = filename, width = 3840, height = 2160, res = 300)
      grDevices::dev.off()

      message(paste("High-resolution 4K plot saved as:", filename))
    }
  }

  # 6. Return Data
  result_list <- list(
    jaccard_matrix = jaccard_mat,
    communities = community_vector,
    hclust_obj = hc_obj
  )

  return(result_list)
}
