#' @title Automated Feature Selection (Native BIC-Driven RFE)
#' @description Zero-dependency Recursive Feature Elimination using Bayesian Information Criterion (BIC).
#' Strips away clinical noise and expensive variables to isolate the 'Elite Biomarkers'
#' (top N parameters) that maximize predictive accuracy while minimizing diagnostic cost.
#' @param data A data.frame containing all candidate clinical predictors.
#' @param target_var Character. The name of the binary outcome variable (0/1).
#' @param max_features Numeric. The maximum number of elite biomarkers to retain (default 5).
#' @param split_ratio Numeric. Proportion of data for training (default 0.7).
#' @param reverse_palette Logical. If TRUE, reverses the plot color palette.
#' @param save_plot Logical. If TRUE, exports 4K high-resolution plots.
#' @param save_prefix Character. Custom prefix for exported plot filenames.
#' @param footer Logical. If TRUE, adds reproducibility footer to plots.
#' @param seed Numeric. Seed for reproducibility.
#' @return A list containing Elite Biomarker names, models (Full vs Elite), and metrics.
#' @export
#' @examples
#' # Run Automated Feature Selection (Keep Top 3)
#' result <- B8_supervised_biomarker(
#'   data = df_B8,
#'   target_var = "disease_status",
#'   max_features = 3,
#'   save_plot = FALSE
#' )
B8_supervised_biomarker <- function(data,
                                    target_var,
                                    max_features = 5,
                                    split_ratio = 0.7,
                                    reverse_palette = FALSE,
                                    save_plot = FALSE,
                                    save_prefix = "",
                                    footer = TRUE,
                                    seed = 42) {

  # 1. Input Validation
  stopifnot("Input 'data' must be a data.frame." = is.data.frame(data))
  stopifnot("Target variable not found in data." = target_var %in% names(data))

  set.seed(seed)

  y <- data[[target_var]]
  if (!is.numeric(y) || length(unique(stats::na.omit(y))) != 2) {
    stop("Target variable must be binary (numeric 0 and 1).")
  }
  if (!all(unique(stats::na.omit(y)) %in% c(0, 1))) {
    stop("Target variable values must be exactly 0 and 1.")
  }

  data <- data[stats::complete.cases(data), ]
  n_total <- nrow(data)
  n_total_vars <- ncol(data) - 1

  # 2. Stratified Train/Test Split
  idx_1 <- which(data[[target_var]] == 1)
  idx_0 <- which(data[[target_var]] == 0)

  if(length(idx_1) < 2 || length(idx_0) < 2) {
    stop("Insufficient class distribution for stratified split.")
  }

  train_1_idx <- sample(idx_1, size = floor(length(idx_1) * split_ratio))
  train_0_idx <- sample(idx_0, size = floor(length(idx_0) * split_ratio))

  train_data <- data[c(train_0_idx, train_1_idx), ]
  test_data <- data[-c(train_0_idx, train_1_idx), ]
  actual <- test_data[[target_var]]

  # 3. Phase I: Train The "Expensive" Full Model
  formula_full <- stats::as.formula(paste(target_var, "~ ."))
  mod_full <- suppressWarnings(stats::glm(formula_full, data = train_data, family = stats::binomial(link = "logit")))

  # 4. Phase II: Native BIC-Driven Recursive Feature Elimination (RFE)
  # penalty k = log(n) enforces strict BIC instead of loose AIC
  penalty_bic <- log(nrow(train_data))
  mod_bic <- suppressWarnings(stats::step(mod_full, direction = "both", k = penalty_bic, trace = 0))

  # Extract selected variables from BIC survival
  surviving_vars <- attr(stats::terms(mod_bic), "term.labels")

  # Failsafe: If BIC strips everything, take the strongest single variable from Full Model
  if (length(surviving_vars) == 0) {
    z_full <- abs(summary(mod_full)$coefficients[-1, "z value", drop = FALSE])
    z_full[is.na(z_full)] <- 0
    surviving_vars <- rownames(z_full)[which.max(z_full)]
  }

  # 5. Phase III: Elite Biomarker Selection & Ranking (Top N)
  # Rank the surviving variables by absolute Wald Z-Statistic (Standardized Importance)
  z_scores <- abs(summary(mod_bic)$coefficients[-1, "z value"])
  z_scores[is.na(z_scores)] <- 0
  names(z_scores) <- surviving_vars

  z_scores_sorted <- sort(z_scores, decreasing = TRUE)
  elite_vars <- names(z_scores_sorted)[1:min(max_features, length(z_scores_sorted))]

  # 6. Phase IV: Train The "Cost-Effective" Elite Model
  formula_elite <- stats::as.formula(paste(target_var, "~", paste(elite_vars, collapse = " + ")))
  mod_elite <- suppressWarnings(stats::glm(formula_elite, data = train_data, family = stats::binomial(link = "logit")))

  # 7. Predictions & AUC Calculation Helper
  prob_full <- stats::predict(mod_full, test_data, type = "response")
  prob_elite <- stats::predict(mod_elite, test_data, type = "response")

  calc_auc <- function(prob, act) {
    ord <- order(prob, decreasing = TRUE)
    act_ord <- act[ord]
    tpr <- cumsum(act_ord == 1) / sum(act_ord == 1)
    fpr <- cumsum(act_ord == 0) / sum(act_ord == 0)
    tpr <- c(0, tpr)
    fpr <- c(0, fpr)
    auc <- sum(diff(fpr) * (tpr[-1] + tpr[-length(tpr)]) / 2)
    return(list(AUC = auc, tpr = tpr, fpr = fpr))
  }

  res_full <- calc_auc(prob_full, actual)
  res_elite <- calc_auc(prob_elite, actual)

  # Elite Odds Ratios
  coef_vals <- stats::coef(mod_elite)[elite_vars]
  ci_vals <- suppressMessages(stats::confint.default(mod_elite))[elite_vars, , drop = FALSE]

  or_data <- data.frame(
    Predictor = elite_vars,
    OR = round(exp(coef_vals), 2),
    CI_Lower = round(exp(ci_vals[, 1]), 2),
    CI_Upper = round(exp(ci_vals[, 2]), 2),
    Importance_Z = z_scores_sorted[elite_vars],
    stringsAsFactors = FALSE
  )
  or_data$Significance <- ifelse(or_data$CI_Lower <= 1.00 & or_data$CI_Upper >= 1.00, "ns", "sig")

  # 8. Visualization Setup
  base_colors <- c("#95A5A6", "#E74C3C", "#18BC9C", "#2C3E50", "#F39C12")
  if (reverse_palette) base_colors <- rev(base_colors)
  footer_text <- paste0("Generated by R-Studio (", R.version.string, ") on ", format(Sys.time(), "%B %d, %Y at %H:%M:%S"))
  timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")

  # ==========================================
  # PLOT 1: Feature Elimination & Cost-Efficiency
  # ==========================================
  # Initialize OMA for centered footer
  graphics::par(mfrow = c(1, 2), mar = c(7, 5, 4, 2) + 0.1, oma = c(3, 0, 0, 0))

  # Subplot 1A: Cost-Efficiency ROC (Full vs Elite)
  graphics::plot(res_full$fpr, res_full$tpr, type = "l", col = grDevices::adjustcolor(base_colors[1], alpha.f = 0.6), lwd = 3,
                 main = "Cost-Efficiency Diagnostic ROC",
                 xlab = "False Positive Rate (1 - Specificity)",
                 ylab = "True Positive Rate (Sensitivity)",
                 xlim = c(0, 1), ylim = c(0, 1), bty = "n")

  graphics::lines(res_elite$fpr, res_elite$tpr, col = base_colors[2], lwd = 3.5)
  graphics::abline(a = 0, b = 1, lty = 2, col = "grey60")

  legend_full <- sprintf("Full Model (%d Vars) - AUC: %.3f", n_total_vars, res_full$AUC)
  legend_elite <- sprintf("Elite Model (%d Vars) - AUC: %.3f", length(elite_vars), res_elite$AUC)

  graphics::legend("bottomright", legend = c(legend_full, legend_elite),
                   col = c(grDevices::adjustcolor(base_colors[1], alpha.f = 0.6), base_colors[2]),
                   lwd = c(3, 3.5), bty = "n", cex = 0.85)

  # Subplot 1B: Biomarker Importance Ranking (Z-Score)
  # Show all surviving variables from BIC (capped at top 10 for visualization)
  plot_z <- sort(z_scores, decreasing = FALSE)
  if(length(plot_z) > 10) plot_z <- plot_z[(length(plot_z)-9):length(plot_z)]

  bar_cols <- ifelse(names(plot_z) %in% elite_vars, base_colors[2], base_colors[1])

  bp <- graphics::barplot(plot_z, horiz = TRUE, col = bar_cols, border = NA,
                          main = "Biomarker Importance Ranking (Wald-Z)",
                          xlab = "Standardized Importance Score",
                          las = 1, cex.names = 0.7, xlim = c(0, max(plot_z) * 1.2))

  graphics::text(x = plot_z, y = bp, labels = round(plot_z, 1), pos = 4, cex = 0.75, col = "black")

  graphics::legend("bottom", legend = c("Elite Biomarker Retained", "Discarded / Lower Rank"),
                   fill = c(base_colors[2], base_colors[1]), bty = "n",
                   horiz = TRUE, inset = c(0, -0.35), xpd = TRUE, cex = 0.9)

  # ABSOLUTE CENTERED FOOTER
  if (footer) graphics::mtext(footer_text, side = 1, line = 1, outer = TRUE, adj = 0.5, cex = 0.7, col = "grey40")

  if (save_plot) {
    file_name1 <- paste0("B8_Topology_", timestamp_str, "_", save_prefix, ".png")
    grDevices::dev.copy(grDevices::png, filename = file_name1, width = 3840, height = 2160, res = 300)
    grDevices::dev.off()
    message("Saved: ", file_name1)
  }

  # ==========================================
  # PLOT 2: Elite Biomarkers Clinical Profile
  # ==========================================
  graphics::par(mfrow = c(1, 2), mar = c(7, 5, 4, 2) + 0.1, oma = c(3, 0, 0, 0))

  # Subplot 2A: Elite Model Probability Density
  dens_0 <- stats::density(prob_elite[actual == 0])
  dens_1 <- stats::density(prob_elite[actual == 1])

  graphics::plot(dens_0, main = "Elite Model Prediction Density",
                 xlab = "Predicted Diagnostic Probability", ylab = "Density",
                 col = base_colors[4], lwd = 2,
                 xlim = c(0, 1), ylim = c(0, max(c(dens_0$y, dens_1$y)) * 1.2), bty = "n")
  graphics::polygon(dens_0, col = grDevices::adjustcolor(base_colors[4], alpha.f = 0.3), border = NA)
  graphics::lines(dens_1, col = base_colors[2], lwd = 2)
  graphics::polygon(dens_1, col = grDevices::adjustcolor(base_colors[2], alpha.f = 0.3), border = NA)
  graphics::abline(v = 0.5, lty = 2, col = "grey50")

  graphics::legend("bottom", legend = c("Control (Actual)", "Disease (Actual)"),
                   fill = c(grDevices::adjustcolor(base_colors[4], alpha.f = 0.3),
                            grDevices::adjustcolor(base_colors[2], alpha.f = 0.3)),
                   bty = "n", horiz = TRUE, inset = c(0, -0.35), xpd = TRUE, cex = 0.9)

  # Subplot 2B: Pillar-Aligned Forest Plot (Odds Ratios)
  n_vars <- nrow(or_data)
  y_pos <- n_vars:1

  point_colors <- ifelse(or_data$Significance == "ns", "#B0B0B0",
                         ifelse(or_data$OR > 1.0, base_colors[2], base_colors[3]))
  font_weights <- ifelse(or_data$Significance == "ns", 1, 2)

  # PILLAR ALIGNMENT ALGORITHM
  x_range <- max(or_data$CI_Upper) - min(or_data$CI_Lower)
  if (x_range == 0) x_range <- 1

  # 60% Left Padding, 50% Right Padding
  x_min_plot <- min(min(or_data$CI_Lower), 0.5) - (x_range * 0.6)
  x_max_plot <- max(max(or_data$CI_Upper), 1.5) + (x_range * 0.5)

  graphics::plot(or_data$OR, y_pos, pch = 16, cex = 1.5, col = point_colors,
                 xlim = c(x_min_plot, x_max_plot), ylim = c(0.5, n_vars + 1),
                 yaxt = "n", xlab = "Odds Ratio (95% CI)", ylab = "",
                 main = "Elite Biomarker Impacts\n(Pillar-Aligned Clinical Boundary)", bty = "n")
  graphics::abline(v = 1.0, lty = 2, col = "black", lwd = 1.5)

  # Render Pillar-Aligned Elements
  for (i in 1:n_vars) {
    graphics::segments(x0 = or_data$CI_Lower[i], y0 = y_pos[i],
                       x1 = or_data$CI_Upper[i], y1 = y_pos[i],
                       col = point_colors[i], lwd = 2)

    # Pillar 1 (Left Aligned Names)
    graphics::text(x = x_min_plot, y = y_pos[i],
                   labels = or_data$Predictor[i], pos = 4, xpd = TRUE,
                   col = point_colors[i], font = font_weights[i], cex = 0.9)

    # Pillar 2 (Right Aligned Values)
    graphics::text(x = x_max_plot, y = y_pos[i],
                   labels = sprintf("%.2f [%.2f, %.2f]", or_data$OR[i], or_data$CI_Lower[i], or_data$CI_Upper[i]),
                   pos = 2, cex = 0.8, col = point_colors[i], xpd = TRUE)
  }

  # ABSOLUTE CENTERED FOOTER (Plot 2)
  if (footer) graphics::mtext(footer_text, side = 1, line = 1, outer = TRUE, adj = 0.5, cex = 0.7, col = "grey40")

  if (save_plot) {
    file_name2 <- paste0("B8_Diagnostics_", timestamp_str, "_", save_prefix, ".png")
    grDevices::dev.copy(grDevices::png, filename = file_name2, width = 3840, height = 2160, res = 300)
    grDevices::dev.off()
    message("Saved: ", file_name2)
  }

  # Reset par & oma
  graphics::par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0))

  # 9. Output Consolidation
  results <- list(
    Initial_Variables = n_total_vars,
    Elite_Biomarkers = elite_vars,
    AUC_Full_Model = res_full$AUC,
    AUC_Elite_Model = res_elite$AUC,
    Elite_OR_Table = or_data,
    Models = list(
      Full = mod_full,
      Elite = mod_elite
    )
  )

  return(invisible(results))
}
