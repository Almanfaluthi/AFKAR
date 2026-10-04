#' @title Diagnostic Screening Classifier Platform (Supervised)
#'
#' @description
#' A robust supervised learning platform using Logistic Regression to screen
#' for a specific disease. Generates two sequential, high-resolution 1x2 dashboards:
#' 1) Clinical Screening Profile (Probability Density & Confusion Metrics).
#' 2) Diagnostic Performance Profile (Test-Set ROC & Predictor Forest Plot).
#' Features strict clinical CI validation to prevent visual rounding errors.
#'
#' @param data A data.frame containing the clinical dataset.
#' @param target Character. The name of the binary target variable (0 = Negative, 1 = Positive).
#' @param features Character vector specifying the predictor columns.
#' @param split_ratio Numeric. Proportion of data used for training (default: 0.7 for 70% Train / 30% Test).
#' @param prob_threshold Numeric. Probability threshold to classify as positive (default: 0.5).
#' @param plot_result Logical. If TRUE, generates two sequential 1x2 panel dashboards.
#' @param colors Character vector of length 5 specifying colors for plots (Negative, Positive, Accuracy, Sens, Spec).
#' @param reverse_palette Logical. If TRUE, reverses the order of the color palette.
#' @param save_plot Logical. If TRUE, saves the generated plots as two separate high-resolution 4K PNGs.
#' @param save_prefix Character. Optional prefix string for the saved plot filenames.
#' @param footer Logical. If TRUE, adds a reproducible generation footer to the plots.
#' @param cex Numeric. Text and symbol scaling factor.
#' @param font Numeric. Font style (1 = plain, 2 = bold, 3 = italic, 4 = bold italic).
#' @param family Character. Font family to use for plotting (e.g., "sans", "serif").
#'
#' @return A list containing the model, odds ratios, performance metrics, and split datasets.
#'
#' @importFrom stats glm as.formula predict coef confint.default complete.cases density
#' @importFrom graphics par plot lines abline segments points text mtext axis legend polygon barplot
#' @importFrom grDevices png dev.copy dev.off
#'
#' @export
#'
#' @examples
#' res_screen <- B1_supervised_screening(
#'   data = df_B1,
#'   target = "DHF_Diagnosis",
#'   features = c("Age", "Temperature", "Platelets_10k", "Hematocrit"),
#'   split_ratio = 0.7,
#'   prob_threshold = 0.4,
#'   plot_result = TRUE,
#'   save_plot = FALSE
#' )

B1_supervised_screening <- function(data,
                                    target,
                                    features,
                                    split_ratio = 0.7,
                                    prob_threshold = 0.5,
                                    plot_result = TRUE,
                                    colors = c("#377eb8", "#e41a1c", "#4daf4a", "#ff7f00", "#984ea3"),
                                    reverse_palette = FALSE,
                                    save_plot = FALSE,
                                    save_prefix = "",
                                    footer = TRUE,
                                    cex = 1,
                                    font = 1,
                                    family = "sans") {

  # 1. Strict Input Validation
  stopifnot("Error: 'data' must be a data.frame." = is.data.frame(data))
  stopifnot("Error: 'target' must be a single character string." = length(target) == 1 && is.character(target))
  stopifnot("Error: 'split_ratio' must be between 0.1 and 0.9." = split_ratio > 0 && split_ratio < 1)

  if (!(target %in% names(data))) stop(paste("Error: Target '", target, "' not found.", sep = ""))
  if (!all(data[[target]] %in% c(0, 1, NA))) stop("Error: Target must be strictly binary (0 and 1).")

  missing_cols <- setdiff(features, names(data))
  if (length(missing_cols) > 0) stop(paste("Error: Features missing:", paste(missing_cols, collapse = ", ")))

  # 2. Data Preparation & Train/Test Split
  df_subset <- data[, c(target, features), drop = FALSE]
  comp_cases <- stats::complete.cases(df_subset)
  df_clean <- df_subset[comp_cases, , drop = FALSE]

  n_total <- nrow(df_clean)
  n_train <- floor(split_ratio * n_total)
  n_test <- n_total - n_train

  if (n_train < length(features) * 10) warning("Medical Insight: Small training sample size risk.")

  train_idx <- sample(seq_len(n_total), size = n_train)
  train_data <- df_clean[train_idx, ]
  test_data <- df_clean[-train_idx, ]

  # 3. Algorithm: Logistic Regression (Training)
  algo_name <- "Logistic Regression"
  formula_str <- paste(target, "~", paste(features, collapse = " + "))
  model <- stats::glm(stats::as.formula(formula_str), data = train_data, family = "binomial")

  # 4. Predictions on TEST SET (Validation)
  pred_probs <- stats::predict(model, newdata = test_data, type = "response")
  pred_class <- ifelse(pred_probs >= prob_threshold, 1, 0)

  # 5. Extract Odds Ratios & Generate Asterisks (WITH STRICT CI VISUAL CHECK)
  coef_est <- stats::coef(model)
  ci_est <- stats::confint.default(model)
  p_vals <- summary(model)$coefficients[, 4]

  get_stars <- function(p, ci_low, ci_high) {
    sapply(seq_along(p), function(i) {
      if (is.na(p[i])) return("ns")

      # CLINICAL FIX: If the plotted (rounded) CI touches or crosses 1.00, force 'ns'
      if (round(ci_low[i], 2) <= 1.00 && round(ci_high[i], 2) >= 1.00) {
        return("ns")
      }

      if (p[i] < 0.001) return("***")
      if (p[i] < 0.01) return("**")
      if (p[i] < 0.05) return("*")
      return("ns")
    })
  }

  or_df <- data.frame(
    Predictor = names(coef_est),
    Odds_Ratio = exp(coef_est),
    CI_Lower = exp(ci_est[, 1]),
    CI_Upper = exp(ci_est[, 2]),
    P_Value = p_vals,
    stringsAsFactors = FALSE
  )

  or_df$Significance <- get_stars(or_df$P_Value, or_df$CI_Lower, or_df$CI_Upper)
  or_df$Label <- sprintf("%.2f [%.2f-%.2f] %s", or_df$Odds_Ratio, or_df$CI_Lower, or_df$CI_Upper, or_df$Significance)

  or_plot_df <- or_df[-1, , drop = FALSE] # Remove Intercept

  # 6. Evaluate Performance on Test Set
  calc_roc <- function(probs, labels) {
    ord <- order(probs, decreasing = TRUE)
    p <- probs[ord]
    l <- labels[ord]
    tpr <- cumsum(l) / sum(l)
    fpr <- cumsum(1 - l) / sum(1 - l)
    auc <- sum(diff(fpr) * (tpr[-1] + tpr[-length(tpr)]) / 2)
    return(list(fpr = c(0, fpr), tpr = c(0, tpr), auc = auc))
  }

  roc_data <- calc_roc(pred_probs, test_data[[target]])

  tp <- sum(pred_class == 1 & test_data[[target]] == 1)
  tn <- sum(pred_class == 0 & test_data[[target]] == 0)
  fp <- sum(pred_class == 1 & test_data[[target]] == 0)
  fn <- sum(pred_class == 0 & test_data[[target]] == 1)

  acc <- (tp + tn) / n_test
  sens <- ifelse((tp + fn) == 0, 0, tp / (tp + fn))
  spec <- ifelse((tn + fp) == 0, 0, tn / (tn + fp))
  ppv <- ifelse((tp + fp) == 0, 0, tp / (tp + fp))
  npv <- ifelse((tn + fn) == 0, 0, tn / (tn + fn))

  perf_metrics <- c(Total_N = n_total, Train_N = n_train, Test_N = n_test,
                    Accuracy = acc, Sensitivity = sens, Specificity = spec,
                    PPV = ppv, NPV = npv, AUC = roc_data$auc)

  # 7. Strict Sequential Visualization (Two 1x2 Panels)
  if (plot_result) {
    if (reverse_palette) colors <- rev(colors)

    col_neg <- colors[1] # Blue (Negative)
    col_pos <- colors[2] # Red (Positive)
    col_nsig <- "#999999" # Gray for non-significant OR
    metric_cols <- colors[1:5]

    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))

    footer_text <- paste0("Generated by R-Studio (", R.version.string, ") on ", format(Sys.time(), "%B %d, %Y at %H:%M:%S"))

    # =========================================================================
    # GAMBAR 1: CLINICAL SCREENING PROFILE (DENSITY & METRICS)
    # =========================================================================
    graphics::par(mfrow = c(1, 2), oma = c(4, 0, 4, 0), family = family, font = font, cex = cex)

    # Panel 1.1: Probability Density Distribution
    graphics::par(mar = c(5, 5, 2, 2) + 0.1)

    dens_0 <- stats::density(pred_probs[test_data[[target]] == 0], from = 0, to = 1)
    dens_1 <- stats::density(pred_probs[test_data[[target]] == 1], from = 0, to = 1)
    y_max <- max(c(dens_0$y, dens_1$y))

    graphics::plot(dens_0, col = col_neg, lwd = 3, main = "Probability Separation (Test Set)",
                   xlab = "Predicted Probability of Disease", ylab = "Patient Density",
                   xlim = c(0, 1), ylim = c(0, y_max * 1.2), las = 1)

    graphics::polygon(dens_0, col = paste0(col_neg, "40"), border = NA)
    graphics::lines(dens_1, col = col_pos, lwd = 3)
    graphics::polygon(dens_1, col = paste0(col_pos, "40"), border = NA)

    graphics::abline(v = prob_threshold, col = "black", lty = 2, lwd = 2)
    graphics::legend("topright", legend = c("Negative Class", "Positive Class", paste("Cut-off:", prob_threshold)),
                     col = c(col_neg, col_pos, "black"), lwd = c(3, 3, 2), lty = c(1, 1, 2), bty = "n", cex = 0.9)

    # Panel 1.2: Clinical Metrics Barplot
    graphics::par(mar = c(5, 5, 2, 2) + 0.1)
    metrics_to_plot <- c(acc, sens, spec, ppv, npv) * 100
    metric_names <- c("Accuracy", "Sensitivity", "Specificity", "PPV", "NPV")

    bp <- graphics::barplot(metrics_to_plot, names.arg = metric_names, col = metric_cols,
                            border = "white", ylim = c(0, 115), las = 1,
                            xlab = "Screening Metrics", ylab = "Percentage (%)",
                            main = paste("Model Evaluation at Threshold", prob_threshold))

    for (i in seq_along(metrics_to_plot)) {
      val <- round(metrics_to_plot[i], 1)
      graphics::text(x = bp[i], y = ifelse(val < 15, val + 8, val - 8),
                     labels = paste0(val, "%"), col = ifelse(val < 15, "black", "white"), font = 2, cex = 0.9)
    }

    header_1 <- "Clinical Screening Profile Dashboard"
    header_2 <- paste0("Algorithm: ", algo_name, " | Total N: ", n_total, " (Train: ", n_train, ", Validation: ", n_test, ")")

    graphics::mtext(header_1, outer = TRUE, side = 3, line = 2, cex = 1.3, font = 2)
    graphics::mtext(header_2, outer = TRUE, side = 3, line = 0.5, cex = 0.9, font = 1, col = "darkblue")
    if (footer) graphics::mtext(footer_text, side = 1, line = 2, outer = TRUE, adj = 0.5, cex = 0.8, col = "dimgray")

    if (save_plot) {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      clean_prefix <- ifelse(nchar(save_prefix) > 0, paste0("_", save_prefix), "")
      filename1 <- paste0("B1_supervised_screening_", timestamp, clean_prefix, "_Part1_Clinical_Profile.png")
      grDevices::dev.copy(grDevices::png, filename = filename1, width = 3840, height = 2160, res = 300)
      grDevices::dev.off()
      message(paste("High-resolution 4K plot (Part 1 - Clinical Profile) saved as:", filename1))
    }

    # =========================================================================
    # GAMBAR 2: DIAGNOSTIC PERFORMANCE PROFILE (ROC & FOREST PLOT)
    # =========================================================================
    graphics::par(mfrow = c(1, 2), oma = c(4, 0, 4, 0), family = family, font = font, cex = cex)

    # Panel 2.1: ROC Curve (Test Set)
    graphics::par(mar = c(5, 5, 2, 2) + 0.1)
    graphics::plot(roc_data$fpr, roc_data$tpr, type = "l", lwd = 3, col = col_pos,
                   xlab = "False Positive Rate (1 - Specificity)",
                   ylab = "True Positive Rate (Sensitivity)",
                   main = "ROC Curve (Validation / Test Set)", las = 1)
    graphics::abline(a = 0, b = 1, lty = 2, col = "gray50")
    legend_text <- c(paste("Test AUC:", round(perf_metrics["AUC"], 3)))
    graphics::legend("bottomright", legend = legend_text, bty = "n", cex = 1.2, text.col = "black")

    # Panel 2.2: Forest Plot w/ CI Text (Train Set)
    graphics::par(mar = c(5, 8, 2, 2) + 0.1)
    or_plot_df <- or_plot_df[order(nrow(or_plot_df):1), ]
    y_coords <- 1:nrow(or_plot_df)

    x_min <- min(c(0.1, min(or_plot_df$CI_Lower)))
    x_max <- max(c(2, max(or_plot_df$CI_Upper))) * 5 # Ruang ekstra

    # Tautkan warna murni dari status Signifikansi yang sudah diverifikasi batas CI-nya
    is_sig <- or_plot_df$Significance != "ns"
    sig_col <- ifelse(is_sig, col_pos, col_nsig)

    graphics::plot(or_plot_df$Odds_Ratio, y_coords, type = "n",
                   xlim = c(x_min, x_max), ylim = c(0.5, nrow(or_plot_df) + 0.5),
                   yaxt = "n", log = "x", xlab = "Odds Ratio (Log Scale)", ylab = "",
                   main = "Predictor Clinical Impact (Train Set)")

    graphics::abline(v = 1, lty = 2, col = "gray50", lwd = 2)
    graphics::axis(2, at = y_coords, labels = or_plot_df$Predictor, las = 1, cex.axis = 0.9)

    for (i in 1:nrow(or_plot_df)) {
      graphics::segments(x0 = or_plot_df$CI_Lower[i], y0 = y_coords[i],
                         x1 = or_plot_df$CI_Upper[i], y1 = y_coords[i], col = sig_col[i], lwd = 2)
      graphics::points(or_plot_df$Odds_Ratio[i], y_coords[i], pch = 15, col = sig_col[i], cex = 1.5)
      graphics::text(x = or_plot_df$CI_Upper[i] * 1.2, y = y_coords[i],
                     labels = or_plot_df$Label[i], pos = 4, cex = 0.9, col = sig_col[i])
    }

    header_3 <- "Diagnostic Performance Profile Dashboard"
    graphics::mtext(header_3, outer = TRUE, side = 3, line = 2, cex = 1.3, font = 2)
    graphics::mtext(header_2, outer = TRUE, side = 3, line = 0.5, cex = 0.9, font = 1, col = "darkblue")
    if (footer) graphics::mtext(footer_text, side = 1, line = 2, outer = TRUE, adj = 0.5, cex = 0.8, col = "dimgray")

    if (save_plot) {
      filename2 <- paste0("B1_supervised_screening_", timestamp, clean_prefix, "_Part2_Diagnostic_Performance.png")
      grDevices::dev.copy(grDevices::png, filename = filename2, width = 3840, height = 2160, res = 300)
      grDevices::dev.off()
      message(paste("High-resolution 4K plot (Part 2 - Diagnostic Performance) saved as:", filename2))
    }
  }

  result_list <- list(
    model = model,
    odds_ratios = or_df,
    performance = perf_metrics,
    data_train = train_data,
    data_test = test_data
  )

  return(result_list)
}
