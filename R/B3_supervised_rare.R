#' @title Imbalanced Event Handler (Rare Disease Prediction)
#' @description Handles highly imbalanced medical datasets (e.g., rare diseases, fatal events)
#' using Native R SMOTE-Lite, Random Oversampling (ROS), or Undersampling (RUS).
#' Integrates strict clinical boundary checks for predictor significance and directional color coding.
#' @param data A data.frame containing clinical/epidemiological data.
#' @param target_var Character. The name of the binary outcome variable (0 = Control, 1 = Rare Event).
#' @param method Character. Balancing method: "SMOTE" (default), "ROS", or "RUS".
#' @param split_ratio Numeric. Proportion of data for training (default 0.7).
#' @param reverse_palette Logical. If TRUE, reverses the plot color palette.
#' @param save_plot Logical. If TRUE, exports 4K high-resolution plots.
#' @param save_prefix Character. Custom prefix for exported plot filenames.
#' @param footer Logical. If TRUE, adds reproducibility footer to plots.
#' @param seed Numeric. Seed for reproducibility.
#' @return A list containing model summaries, performance metrics, and balanced training data.
#' @export
#' @examples
#' # Run Model
#' result <- B3_supervised_rare(
#'   data = df_B3,
#'   target_var = "anaphylaxis",
#'   method = "SMOTE",
#'   save_plot = FALSE
#' )
B3_supervised_rare <- function(data,
                               target_var,
                               method = c("SMOTE", "ROS", "RUS"),
                               split_ratio = 0.7,
                               reverse_palette = FALSE,
                               save_plot = FALSE,
                               save_prefix = "",
                               footer = TRUE,
                               seed = 42) {

  # 1. Input Validation
  stopifnot("Input 'data' must be a data.frame." = is.data.frame(data))
  stopifnot("Target variable not found in data." = target_var %in% names(data))

  method <- match.arg(method)
  set.seed(seed)

  # Ensure target is binary 0/1 numeric
  y <- data[[target_var]]
  if (!is.numeric(y) || length(unique(stats::na.omit(y))) != 2) {
    stop("Target variable must be binary (numeric 0 and 1) representing Control (0) and Rare Event (1).")
  }
  if (!all(unique(stats::na.omit(y)) %in% c(0, 1))) {
    stop("Target variable values must be exactly 0 and 1.")
  }

  # Clean NA from data
  data <- data[stats::complete.cases(data), ]
  n_total <- nrow(data)

  # 2. Train/Test Split (Stratified-like by random sampling per class)
  idx_1 <- which(data[[target_var]] == 1)
  idx_0 <- which(data[[target_var]] == 0)

  train_1_idx <- sample(idx_1, size = floor(length(idx_1) * split_ratio))
  train_0_idx <- sample(idx_0, size = floor(length(idx_0) * split_ratio))

  test_1_idx <- setdiff(idx_1, train_1_idx)
  test_0_idx <- setdiff(idx_0, train_0_idx)

  train_data <- data[c(train_0_idx, train_1_idx), ]
  test_data <- data[c(test_0_idx, test_1_idx), ]

  # 3. Balancing Engine (Native R) - Applied ONLY to Train Data
  t_idx_1 <- which(train_data[[target_var]] == 1)
  t_idx_0 <- which(train_data[[target_var]] == 0)
  n_maj <- length(t_idx_0)
  n_min <- length(t_idx_1)

  if (n_min == 0) stop("Training data has no rare events. Increase data size.")
  if (n_min > n_maj) warning("Rare event is actually the majority in this dataset.")

  balanced_train <- train_data

  if (method == "ROS") {
    syn_idx <- sample(t_idx_1, size = (n_maj - n_min), replace = TRUE)
    balanced_train <- rbind(train_data, train_data[syn_idx, ])
  } else if (method == "RUS") {
    keep_idx <- sample(t_idx_0, size = n_min, replace = FALSE)
    balanced_train <- rbind(train_data[t_idx_1, ], train_data[keep_idx, ])
  } else if (method == "SMOTE") {
    n_synthetic <- n_maj - n_min
    if (n_synthetic > 0) {
      syn_data <- data.frame(matrix(ncol = ncol(train_data), nrow = n_synthetic))
      colnames(syn_data) <- colnames(train_data)

      pair1_idx <- sample(t_idx_1, size = n_synthetic, replace = TRUE)
      pair2_idx <- sample(t_idx_1, size = n_synthetic, replace = TRUE)

      for (col in names(train_data)) {
        if (is.numeric(train_data[[col]])) {
          gap <- stats::runif(n_synthetic, 0, 1)
          syn_data[[col]] <- train_data[[col]][pair1_idx] + gap * (train_data[[col]][pair2_idx] - train_data[[col]][pair1_idx])
        } else {
          pick <- ifelse(stats::runif(n_synthetic) > 0.5, pair1_idx, pair2_idx)
          syn_data[[col]] <- train_data[[col]][pick]
        }
      }
      balanced_train <- rbind(train_data, syn_data)
    }
  }

  # 4. Modeling (Logistic Regression Base R)
  formula_obj <- stats::as.formula(paste(target_var, "~ ."))
  model <- suppressWarnings(stats::glm(formula_obj, data = balanced_train, family = stats::binomial(link = "logit")))

  # 5. Prediction and Metrics on Unseen Test Data
  test_probs <- stats::predict(model, newdata = test_data, type = "response")
  test_preds <- ifelse(test_probs > 0.5, 1, 0)
  actual <- test_data[[target_var]]

  tp <- sum(test_preds == 1 & actual == 1)
  tn <- sum(test_preds == 0 & actual == 0)
  fp <- sum(test_preds == 1 & actual == 0)
  fn <- sum(test_preds == 0 & actual == 1)

  sensitivity <- ifelse((tp + fn) == 0, 0, tp / (tp + fn))
  specificity <- ifelse((tn + fp) == 0, 0, tn / (tn + fp))
  accuracy <- (tp + tn) / length(actual)

  # Fast Base R AUC Calculation
  pred_ord <- order(test_probs, decreasing = TRUE)
  sorted_act <- actual[pred_ord]
  tpr_vec <- cumsum(sorted_act == 1) / sum(sorted_act == 1)
  fpr_vec <- cumsum(sorted_act == 0) / sum(sorted_act == 0)
  tpr_vec <- c(0, tpr_vec)
  fpr_vec <- c(0, fpr_vec)
  auc_val <- sum(diff(fpr_vec) * (tpr_vec[-1] + tpr_vec[-length(tpr_vec)]) / 2)

  # 6. Strict Clinical Boundary Check for Odds Ratios
  coef_vals <- stats::coef(model)[-1] # Remove intercept
  ci_vals <- suppressMessages(stats::confint.default(model))[-1, , drop = FALSE]

  or_data <- data.frame(
    Predictor = names(coef_vals),
    OR = round(exp(coef_vals), 2),
    CI_Lower = round(exp(ci_vals[, 1]), 2),
    CI_Upper = round(exp(ci_vals[, 2]), 2),
    stringsAsFactors = FALSE
  )

  # Core Logic: If CI crosses 1.00, force label "ns"
  or_data$Significance <- ifelse(or_data$CI_Lower <= 1.00 & or_data$CI_Upper >= 1.00, "ns", "sig")

  # 7. Visualization Setup
  base_colors <- c("#2C3E50", "#E74C3C", "#18BC9C", "#F39C12") # 2 = Red (Risk), 3 = Green (Protective)
  if (reverse_palette) base_colors <- rev(base_colors)
  footer_text <- paste0("Generated by R-Studio (", R.version.string, ") on ", format(Sys.time(), "%B %d, %Y at %H:%M:%S"))

  timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")

  # ==========================================
  # PLOT 1: Screening & Data Topology Profile
  # ==========================================
  graphics::par(mar = c(8, 5, 4, 2) + 0.1, mfrow = c(1, 2))

  # Subplot 1A: Imbalance Profile (Barplot)
  orig_counts <- table(train_data[[target_var]])
  bal_counts <- table(balanced_train[[target_var]])
  plot_mat <- rbind(orig_counts, bal_counts)
  rownames(plot_mat) <- c("Original", paste("Balanced", method))

  bp <- graphics::barplot(plot_mat, beside = TRUE, col = base_colors[1:2],
                          main = "Class Distribution Profile",
                          ylab = "Number of Cases",
                          xlab = "",
                          ylim = c(0, max(plot_mat) * 1.2),
                          border = NA)
  graphics::text(x = bp, y = plot_mat + (max(plot_mat) * 0.05), labels = plot_mat, cex = 0.9)

  graphics::legend("bottom", legend = rownames(plot_mat), fill = base_colors[1:2],
                   bty = "n", horiz = TRUE, inset = c(0, -0.3), xpd = TRUE, cex = 1)

  # Subplot 1B: Probability Density of Predictions
  dens_0 <- stats::density(test_probs[actual == 0])
  dens_1 <- stats::density(test_probs[actual == 1])

  graphics::plot(dens_0, main = "Probability Density (Test Set)",
                 xlab = "",
                 ylab = "Density", col = base_colors[1], lwd = 2,
                 xlim = c(0, 1), ylim = c(0, max(c(dens_0$y, dens_1$y))),
                 bty = "n")
  graphics::polygon(dens_0, col = grDevices::adjustcolor(base_colors[1], alpha.f = 0.3), border = NA)
  graphics::lines(dens_1, col = base_colors[2], lwd = 2)
  graphics::polygon(dens_1, col = grDevices::adjustcolor(base_colors[2], alpha.f = 0.3), border = NA)
  graphics::abline(v = 0.5, lty = 2, col = "grey50")

  graphics::legend("bottom", legend = c("Control", "Rare Event"),
                   fill = grDevices::adjustcolor(base_colors[1:2], alpha.f = 0.3),
                   bty = "n", horiz = TRUE, inset = c(0, -0.3), xpd = TRUE, cex = 1)

  if (footer) graphics::mtext(footer_text, side = 1, line = 6.5, adj = 0, cex = 0.7, col = "grey40", outer = FALSE)

  # Export Plot 1
  if (save_plot) {
    file_name1 <- paste0("B3_Topology_", timestamp_str, "_", save_prefix, ".png")
    grDevices::dev.copy(grDevices::png, filename = file_name1, width = 3840, height = 2160, res = 300)
    grDevices::dev.off()
    message("Saved: ", file_name1)
  }

  # ==========================================
  # PLOT 2: Diagnostic Performance (ROC & Forest)
  # ==========================================
  graphics::par(mfrow = c(1, 2), mar = c(7, 4, 4, 2) + 0.1)

  # Subplot 2A: ROC Curve
  graphics::plot(fpr_vec, tpr_vec, type = "l", col = base_colors[3], lwd = 3,
                 main = sprintf("ROC Curve (AUC = %.3f)", auc_val),
                 xlab = "False Positive Rate (1 - Specificity)",
                 ylab = "True Positive Rate (Sensitivity)",
                 xlim = c(0, 1), ylim = c(0, 1), bty = "n")
  graphics::polygon(c(fpr_vec, 1, 0), c(tpr_vec, 0, 0),
                    col = grDevices::adjustcolor(base_colors[3], alpha.f = 0.2), border = NA)
  graphics::abline(a = 0, b = 1, lty = 2, col = "grey60")

  metrics_txt <- sprintf("Sensitivity: %.1f%%\nSpecificity: %.1f%%\nAccuracy: %.1f%%",
                         sensitivity * 100, specificity * 100, accuracy * 100)
  graphics::legend("bottomright", legend = metrics_txt, bty = "n", text.col = "black", cex = 0.9)

  # Subplot 2B: Directional Clinical Boundary Forest Plot
  n_vars <- nrow(or_data)
  y_pos <- n_vars:1

  # Directional Color Coding Logic
  point_colors <- ifelse(or_data$Significance == "ns", "#B0B0B0",
                         ifelse(or_data$OR > 1.0, base_colors[2], base_colors[3]))
  font_weights <- ifelse(or_data$Significance == "ns", 1, 2)

  graphics::plot(or_data$OR, y_pos, pch = 16, cex = 1.5, col = point_colors,
                 xlim = c(min(or_data$CI_Lower, 0.5), max(or_data$CI_Upper) + 1),
                 ylim = c(0.5, n_vars + 1),
                 yaxt = "n", xlab = "Odds Ratio (95% CI)", ylab = "",
                 main = "Predictor Odds Ratios\n(Directional Boundary)",
                 bty = "n")
  graphics::abline(v = 1.0, lty = 2, col = "black", lwd = 1.5)

  # Draw directional error bars and labels
  for (i in 1:n_vars) {
    graphics::segments(x0 = or_data$CI_Lower[i], y0 = y_pos[i],
                       x1 = or_data$CI_Upper[i], y1 = y_pos[i],
                       col = point_colors[i], lwd = 2)
    graphics::text(x = min(or_data$CI_Lower, 0.5), y = y_pos[i],
                   labels = or_data$Predictor[i], pos = 2, xpd = TRUE,
                   col = point_colors[i], font = font_weights[i])
    graphics::text(x = max(or_data$CI_Upper) + 0.1, y = y_pos[i],
                   labels = sprintf("%.2f [%.2f-%.2f]", or_data$OR[i], or_data$CI_Lower[i], or_data$CI_Upper[i]),
                   pos = 4, cex = 0.8, col = point_colors[i])
  }

  if (footer) graphics::mtext(footer_text, side = 1, line = 5.5, adj = 0.5, cex = 0.7, col = "grey40")

  # Export Plot 2
  if (save_plot) {
    file_name2 <- paste0("B3_Diagnostics_", timestamp_str, "_", save_prefix, ".png")
    grDevices::dev.copy(grDevices::png, filename = file_name2, width = 3840, height = 2160, res = 300)
    grDevices::dev.off()
    message("Saved: ", file_name2)
  }

  graphics::par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

  # 8. Return Objects
  results <- list(
    Balancing_Method = method,
    Original_Dist = orig_counts,
    Balanced_Dist = bal_counts,
    Metrics = data.frame(
      Accuracy = accuracy,
      Sensitivity = sensitivity,
      Specificity = specificity,
      AUC = auc_val
    ),
    Strict_OR_Table = or_data,
    Model = model
  )

  return(invisible(results))
}
