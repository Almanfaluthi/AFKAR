#' @title Association Rule Mining for Clinical Pathways (Unsupervised)
#'
#' @description
#' Discovers pairwise clinical association rules (e.g., Comorbidities -> Treatments)
#' using a native Base R implementation. Calculates Support, Confidence, and Lift
#' for binary clinical datasets (Electronic Medical Records/EMR).
#' Includes dual visualizations: a Support-Confidence Scatter Space and
#' a Top Rules Barplot.
#'
#' @param data A data.frame of binary numeric columns (1 = Present, 0 = Absent) representing patients (rows) and clinical items (columns).
#' @param min_support Numeric. Minimum support threshold (proportion of total cases, 0 to 1).
#' @param min_confidence Numeric. Minimum confidence threshold (0 to 1).
#' @param top_n Integer. Number of top rules (sorted by Lift) to display in the barplot.
#' @param plot_result Logical. If TRUE, generates a 1x2 panel plot.
#' @param colors Character vector of length 2 for the low and high color gradient (default: c("#ADD8E6", "#00008B")).
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
#'   \item \code{rules}: A data.frame of discovered rules sorted by Lift.
#'   \item \code{summary}: A numeric vector summarizing the mining process.
#' }
#'
#' @importFrom graphics par plot points mtext barplot text
#' @importFrom grDevices colorRampPalette png dev.copy dev.off
#' @importFrom stats complete.cases
#'
#' @export
#'
#' @examples
#' # Run Association Mining
#' res_assoc <- A4_unsupervised_association(
#'   data = df_A4,
#'   min_support = 0.05,
#'   min_confidence = 0.50,
#'   top_n = 5,
#'   plot_result = TRUE,
#'   save_plot = FALSE
#' )
#'
#' # View the strongest rules
#' print(head(res_assoc$rules))

A4_unsupervised_association <- function(data,
                                        min_support = 0.1,
                                        min_confidence = 0.5,
                                        top_n = 10,
                                        plot_result = TRUE,
                                        colors = c("#ADD8E6", "#00008B"),
                                        reverse_palette = FALSE,
                                        save_plot = FALSE,
                                        save_prefix = "",
                                        footer = TRUE,
                                        cex = 1,
                                        font = 1,
                                        family = "sans") {

  # 1. Input Validation
  stopifnot("Error: 'data' must be a data.frame." = is.data.frame(data))

  # Check if all data is binary numeric (0 or 1)
  for (col in names(data)) {
    if (!is.numeric(data[[col]]) || !all(data[[col]] %in% c(0, 1, NA))) {
      stop(paste("Error: Column '", col, "' must be binary numeric (0 and 1 only) for association mining.", sep = ""))
    }
  }

  # 2. Data Preparation
  comp_cases <- stats::complete.cases(data)
  df_clean <- data[comp_cases, , drop = FALSE]
  n_total <- nrow(df_clean)

  if (n_total == 0) {
    stop("Error: No complete cases found. Please handle missing values.")
  }

  # 3. Core Algorithm: Native Pairwise Rule Mining (A => B)
  items <- names(df_clean)
  n_items <- length(items)

  # Pre-calculate item supports to optimize loop
  item_supports <- colSums(df_clean) / n_total

  rule_list <- list()
  counter <- 1

  for (i in seq_len(n_items)) {
    for (j in seq_len(n_items)) {
      if (i != j) {
        supp_A <- item_supports[i]
        supp_B <- item_supports[j]

        # Calculate joint probability P(A & B)
        joint_cases <- sum(df_clean[[items[i]]] == 1 & df_clean[[items[j]]] == 1)
        supp_AB <- joint_cases / n_total

        # Apply min_support threshold
        if (supp_AB >= min_support) {
          conf_AB <- supp_AB / supp_A

          # Apply min_confidence threshold
          if (conf_AB >= min_confidence) {
            lift_AB <- conf_AB / supp_B

            rule_list[[counter]] <- data.frame(
              Rule = paste(items[i], "=>", items[j]),
              Antecedent = items[i],
              Consequent = items[j],
              Support = supp_AB,
              Confidence = conf_AB,
              Lift = lift_AB,
              stringsAsFactors = FALSE
            )
            counter <- counter + 1
          }
        }
      }
    }
  }

  # Combine rules and sort by Lift (descending)
  if (length(rule_list) == 0) {
    warning("Medical Insight: No rules found meeting the minimum support and confidence thresholds. Try lowering the parameters.")
    return(NULL)
  }

  rules_df <- do.call(rbind, rule_list)
  rules_df <- rules_df[order(-rules_df$Lift), ]
  rownames(rules_df) <- NULL

  # 4. Strict Visualization (1x2 Panel Plot)
  if (plot_result) {

    if (reverse_palette) {
      colors <- rev(colors)
    }

    # Save old par
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))

    # Setup global layout with oma for footer
    graphics::par(mfrow = c(1, 2),
                  oma = c(4, 0, 0, 0),
                  family = family, font = font, cex = cex)

    # ================= PLOT 1: SCATTER SPACE =================
    # Adjust margin for Plot 1
    graphics::par(mar = c(5, 5, 4, 2) + 0.1)

    # Create color gradient based on Lift
    color_pal <- grDevices::colorRampPalette(c(colors[1], colors[2]))(100)

    # Normalize Lift to 1-100 index for color mapping
    min_lift <- min(rules_df$Lift)
    max_lift <- max(rules_df$Lift)

    if (max_lift == min_lift) {
      lift_indices <- rep(50, nrow(rules_df))
    } else {
      lift_indices <- as.integer(((rules_df$Lift - min_lift) / (max_lift - min_lift)) * 99) + 1
    }

    point_colors <- color_pal[lift_indices]

    graphics::plot(rules_df$Support, rules_df$Confidence,
                   pch = 21,
                   bg = point_colors,
                   col = "darkgray",
                   cex = 1.5 + (rules_df$Lift / max(rules_df$Lift)), # Size scaled by lift
                   xlab = "Support (Case Frequency)",
                   ylab = "Confidence (Probability A => B)",
                   main = "Rule Scatter Space",
                   sub = paste("Total Rules Found:", nrow(rules_df)))

    # ================= PLOT 2: TOP RULES BARPLOT =================
    # Adjust margin for Plot 2 (Larger left margin for long rule names)
    graphics::par(mar = c(5, 9, 4, 2) + 0.1)

    n_plot <- min(top_n, nrow(rules_df))
    top_rules <- rules_df[1:n_plot, ]

    # Reverse order so the highest is at the top of the barplot
    top_rules <- top_rules[order(top_rules$Lift, decreasing = FALSE), ]

    bp <- graphics::barplot(top_rules$Confidence,
                            names.arg = top_rules$Rule,
                            horiz = TRUE,
                            las = 1,          # Horizontal labels
                            col = colors[2],
                            border = "white",
                            xlim = c(0, max(top_rules$Confidence) * 1.2),
                            xlab = "Confidence",
                            main = paste("Top", n_plot, "Clinical Pathways\n(Ranked by Lift)"))

    # Add text labels inside/next to bars showing Confidence %
    graphics::text(x = top_rules$Confidence,
                   y = bp,
                   labels = paste0(round(top_rules$Confidence * 100, 1), "%"),
                   pos = 4,
                   cex = 0.8,
                   col = "black")

    # ================= GLOBAL FOOTER =================
    if (footer) {
      footer_text <- paste0("Generated by R-Studio (", R.version.string, ") on ", format(Sys.time(), "%B %d, %Y at %H:%M:%S"))
      graphics::mtext(footer_text, side = 1, line = 2, outer = TRUE, adj = 0.5, cex = 0.8, col = "dimgray")
    }

    # ================= HIGH-RES 4K EXPORT =================
    if (save_plot) {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      clean_prefix <- ifelse(nchar(save_prefix) > 0, paste0("_", save_prefix), "")
      filename <- paste0("A4_unsupervised_association_", timestamp, clean_prefix, ".png")

      grDevices::dev.copy(grDevices::png, filename = filename, width = 3840, height = 2160, res = 300)
      grDevices::dev.off()

      message(paste("High-resolution 4K plot saved as:", filename))
    }
  }

  # 5. Return List
  summary_stats <- c(
    Total_Patients = n_total,
    Items_Evaluated = n_items,
    Rules_Found = nrow(rules_df)
  )

  result_list <- list(
    rules = rules_df,
    summary = summary_stats
  )

  return(result_list)
}
