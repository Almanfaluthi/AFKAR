#' @title Emergency Triage & Deterioration Predictor (Supervised)
#'
#' @description
#' Trains a clinical early warning model (Logistic Regression) to predict
#' rapid deterioration (e.g., ICU admission). Features automatic Train/Test
#' splitting and generates two separate dual-panel visualizations:
#' 1) Diagnostic Performance (ROC & Predictor Impact)
#' 2) Triage Calibration (Risk Separation & Stratification).
#'
#' @param data A data.frame containing the clinical triage dataset.
#' @param target Character. The binary target variable (0 = Stable, 1 = Deterioration).
#' @param features Character vector specifying the predictor columns.
#' @param split_ratio Numeric. Proportion of data used for training (default: 0.7).
#' @param risk_thresholds Numeric vector of length 2 defining the probability cutoffs for Moderate and High Risk (default: c(0.2, 0.5)).
#' @param plot_result Logical. If TRUE, generates two sequential 1x2 panel dashboards.
#' @param colors Character vector of length 5 for plotting elements (ROC, Signif_OR/High_Risk, Low_Risk, Mod_Risk, Non_Signif_OR).
#' @param reverse_palette Logical. If TRUE, reverses the order of the color palette.
#' @param save_plot Logical. If TRUE, saves the generated plots as two separate high-resolution 4K PNGs.
#' @param save_prefix Character. Optional prefix string for the saved plot filenames.
#' @param footer Logical. If TRUE, adds a reproducible generation footer to the plots.
#' @param cex Numeric. Text and symbol scaling factor.
#' @param font Numeric. Font style (1 = plain, 2 = bold, 3 = italic, 4 = bold italic).
#' @param family Character. Font family to use for plotting (e.g., "sans", "serif").
#'
#' @return A list containing the model, odds ratios, triage calibration, performance metrics, and split datasets.
#'
#' @importFrom stats glm as.formula predict coef confint.default complete.cases
#' @importFrom graphics par boxplot barplot text mtext abline legend axis segments points
#' @importFrom grDevices png dev.copy dev.off
#'
#' @export
#'
#' @examples
#' res_triage <- B2_supervised_triage(
#'   data = df_B2,
#'   target = "ICU_Admission",
#'   features = c("Heart_Rate", "Resp_Rate", "Systolic_BP", "Temperature"),
#'   split_ratio = 0.7,
#'   plot_result = TRUE,
#'   save_plot = FALSE
#' )

B2_supervised_triage <- function(data,
                                 target,
                                 features,
                                 split_ratio = 0.7,
                                 risk_thresholds = c(0.2, 0.5),
                                 plot_result = TRUE,
                                 colors = c("#377eb8", "#e41a1c", "#4daf4a", "#ff7f00", "#999999"),
                                 reverse_palette = FALSE,
                                 save_plot = FALSE,
                                 save_prefix = "",
                                 footer = TRUE,
                                 cex = 1,
                                 font = 1,
                                 family = "sans") {

  # 1. Strict Input Validation
  stopifnot("Error: 'data' must be a data.frame." = is.data.frame(data))
  stopifnot("Error: 'target' must be a single string." = length(target) == 1 && is.character(target))
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

  # 3. Core Algorithm: Logistic Regression (Training)
  algo_name <- "EWS Probabilistic Regression"
  formula_str <- paste(target, "~", paste(features, collapse = " + "))
  model <- stats::glm(stats::as.formula(formula_str), data = train_data, family = "binomial")

  # 4. Extract Odds Ratios (Train Set)
  coef_est <- stats::coef(model)
  ci_est <- stats::confint.default(model)
  p_vals <- summary(model)$coefficients[, 4]

  get_stars <- function(p) {
    sapply(p, function(x) {
      if (is.na(x)) return("")
      if (x < 0.001) return("***")
      if (x < 0.01) return("**")
      if (x < 0.05) return("*")
      return("ns")
    })
  }

  or_df <- data.frame(
    Predictor = names(coef_est), Odds_Ratio = exp(coef_est),
    CI_Lower = exp(ci_est[, 1]), CI_Upper = exp(ci_est[, 2]),
    P_Value = p_vals, Significance = get_stars(p_vals), stringsAsFactors = FALSE
  )

  or_df$Label <- sprintf("%.2f [%.2f-%.2f] %s", or_df$Odds_Ratio, or_df$CI_Lower, or_df$CI_Upper, or_df$Significance)
  or_plot_df <- or_df[-1, , drop = FALSE] # Remove Intercept

  # 5. Triage Predictions on TEST SET (Validation)
  pred_probs <- stats::predict(model, newdata = test_data, type = "response")

  risk_tiers <- rep("1_Low_Risk", n_test)
  risk_tiers[pred_probs >= risk_thresholds[1]] <- "2_Moderate_Risk"
  risk_tiers[pred_probs >= risk_thresholds[2]] <- "3_High_Risk"

  test_data$Predicted_Probability <- pred_probs
  test_data$Risk_Tier <- risk_tiers

  # 6. Performance & Triage Calibration (Test Set)
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
  auc_val <- roc_data$auc

  tier_summary <- data.frame(
    Tier = c("1_Low_Risk", "2_Moderate_Risk", "3_High_Risk"),
    Total_Patients = 0, Actual_Deterioration = 0, Deterioration_Rate = 0, stringsAsFactors = FALSE
  )

  for (i in 1:3) {
    tier_name <- tier_summary$Tier[i]
    tier_pts <- test_data[test_data$Risk_Tier == tier_name, ]
    tier_summary$Total_Patients[i] <- nrow(tier_pts)
    if (nrow(tier_pts) > 0) {
      tier_summary$Actual_Deterioration[i] <- sum(tier_pts[[target]] == 1)
      tier_summary$Deterioration_Rate[i] <- tier_summary$Actual_Deterioration[i] / tier_summary$Total_Patients[i]
    }
  }

  perf_metrics <- c(Total_N = n_total, Train_N = n_train, Test_N = n_test, Test_Set_AUC = auc_val)

  # 7. Strict Sequential Visualization (Two 1x2 Panels)
  if (plot_result) {
    if (reverse_palette) colors <- rev(colors)

    col_roc <- colors[1]
    col_sig <- colors[2]
    col_nsig <- colors[5]
    col_triage <- c(colors[3], colors[4], colors[2]) # Low, Mod, High

    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))

    # =========================================================================
    # GAMBAR 1: DIAGNOSTIC PERFORMANCE PROFILE (ROC & FOREST PLOT)
    # =========================================================================
    graphics::par(mfrow = c(1, 2), oma = c(4, 0, 4, 0), family = family, font = font, cex = cex)

    # Panel 1.1: ROC Curve
    graphics::par(mar = c(5, 5, 2, 2) + 0.1)
    graphics::plot(roc_data$fpr, roc_data$tpr, type = "l", lwd = 3, col = col_roc,
                   xlab = "False Positive Rate (1 - Specificity)", ylab = "True Positive Rate (Sensitivity)",
                   main = "ROC Curve (Validation / Test Set)")
    graphics::abline(a = 0, b = 1, lty = 2, col = "gray50")
    graphics::legend("bottomright", legend = paste("Test AUC:", round(auc_val, 3)), bty = "n", cex = 1.1)

    # Panel 1.2: Forest Plot
    graphics::par(mar = c(5, 8, 2, 2) + 0.1)
    or_plot_df <- or_plot_df[order(nrow(or_plot_df):1), ]
    y_coords <- 1:nrow(or_plot_df)

    x_min <- min(c(0.1, min(or_plot_df$CI_Lower)))
    x_max <- max(c(2, max(or_plot_df$CI_Upper))) * 5 # Ruang ekstra untuk teks

    point_cols <- ifelse(or_plot_df$P_Value < 0.05, col_sig, col_nsig)

    graphics::plot(or_plot_df$Odds_Ratio, y_coords, type = "n",
                   xlim = c(x_min, x_max), ylim = c(0.5, nrow(or_plot_df) + 0.5),
                   yaxt = "n", log = "x", xlab = "Odds Ratio (Log Scale)", ylab = "",
                   main = "Predictor Impact (Train Set)")
    graphics::abline(v = 1, lty = 2, col = "gray50", lwd = 2)
    graphics::axis(2, at = y_coords, labels = or_plot_df$Predictor, las = 1, cex.axis = 0.9)

    for (i in 1:nrow(or_plot_df)) {
      graphics::segments(x0 = or_plot_df$CI_Lower[i], y0 = y_coords[i],
                         x1 = or_plot_df$CI_Upper[i], y1 = y_coords[i], col = point_cols[i], lwd = 2)
      graphics::points(or_plot_df$Odds_Ratio[i], y_coords[i], pch = 15, col = point_cols[i], cex = 1.5)
      graphics::text(x = or_plot_df$CI_Upper[i] * 1.2, y = y_coords[i],
                     labels = or_plot_df$Label[i], pos = 4, cex = 0.85, col = point_cols[i])
    }

    header_1 <- "Diagnostic Performance Profile"
    header_2 <- paste0("Algorithm: ", algo_name, " | Total N: ", n_total, " (Train: ", n_train, ", Valid: ", n_test, ")")
    graphics::mtext(header_1, outer = TRUE, side = 3, line = 2, cex = 1.3, font = 2)
    graphics::mtext(header_2, outer = TRUE, side = 3, line = 0.5, cex = 0.9, font = 1, col = "darkblue")

    if (footer) {
      footer_text <- paste0("Generated by R-Studio (", R.version.string, ") on ", format(Sys.time(), "%B %d, %Y at %H:%M:%S"))
      graphics::mtext(footer_text, side = 1, line = 2, outer = TRUE, adj = 0.5, cex = 0.8, col = "dimgray")
    }

    if (save_plot) {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      clean_prefix <- ifelse(nchar(save_prefix) > 0, paste0("_", save_prefix), "")
      filename1 <- paste0("B2_supervised_triage_", timestamp, clean_prefix, "_Part1_Diagnostic.png")
      grDevices::dev.copy(grDevices::png, filename = filename1, width = 3840, height = 2160, res = 300)
      grDevices::dev.off()
      message(paste("High-resolution 4K plot (Part 1 - Diagnostics) saved as:", filename1))
    }

    # =========================================================================
    # GAMBAR 2: EMERGENCY TRIAGE CALIBRATION (BOXPLOT & BARPLOT)
    # =========================================================================
    # Reset par untuk gambar kedua yang terpisah
    graphics::par(mfrow = c(1, 2), oma = c(4, 0, 4, 0), family = family, font = font, cex = cex)

    # Panel 2.1: Risk Separation Boxplot
    graphics::par(mar = c(5, 5, 2, 2) + 0.1)
    graphics::boxplot(test_data$Predicted_Probability ~ test_data[[target]],
                      col = c(col_roc, col_sig),
                      xlab = "Actual Outcome (0=Stable, 1=ICU)", ylab = "Probability of Deterioration",
                      main = "Risk Separation (Test Set)", las = 1, ylim = c(0, 1))

    graphics::abline(h = risk_thresholds[1], col = col_triage[2], lty = 2, lwd = 2)
    graphics::abline(h = risk_thresholds[2], col = col_triage[3], lty = 2, lwd = 2)
    graphics::legend("topleft", legend = c("Mod Risk Threshold", "High Risk Threshold"),
                     col = c(col_triage[2], col_triage[3]), lty = 2, lwd = 2, bty = "n", cex = 0.8)

    # Panel 2.2: Triage Calibration Barplot
    graphics::par(mar = c(5, 5, 2, 2) + 0.1)
    bp <- graphics::barplot(tier_summary$Deterioration_Rate * 100,
                            names.arg = c("Low Risk", "Moderate Risk", "High Risk"),
                            col = col_triage, border = "white", ylim = c(0, 110),
                            xlab = "Triage Stratification", ylab = "Actual Deterioration Rate (%)",
                            main = "EWS Calibration Validation (Test Set)")

    for (i in 1:3) {
      rate_pct <- round(tier_summary$Deterioration_Rate[i] * 100, 1)
      n_pts_tier <- tier_summary$Total_Patients[i]
      label_txt <- paste0(rate_pct, "%\n(n=", n_pts_tier, ")")

      graphics::text(x = bp[i], y = ifelse(rate_pct < 15, rate_pct + 10, rate_pct - 10),
                     labels = label_txt, col = ifelse(rate_pct < 15, "black", "white"), font = 2, cex = 1)
    }

    header_3 <- "Emergency Triage Calibration Dashboard"
    header_4 <- paste0("Algorithm: ", algo_name, " | Validation Test Set (N = ", n_test, ")")
    graphics::mtext(header_3, outer = TRUE, side = 3, line = 2, cex = 1.3, font = 2)
    graphics::mtext(header_4, outer = TRUE, side = 3, line = 0.5, cex = 0.9, font = 1, col = "darkblue")

    if (footer) {
      graphics::mtext(footer_text, side = 1, line = 2, outer = TRUE, adj = 0.5, cex = 0.8, col = "dimgray")
    }

    if (save_plot) {
      filename2 <- paste0("B2_supervised_triage_", timestamp, clean_prefix, "_Part2_Calibration.png")
      grDevices::dev.copy(grDevices::png, filename = filename2, width = 3840, height = 2160, res = 300)
      grDevices::dev.off()
      message(paste("High-resolution 4K plot (Part 2 - Calibration) saved as:", filename2))
    }
  }

  result_list <- list(
    model = model,
    odds_ratios = or_df,
    triage_calibration = tier_summary,
    performance = perf_metrics,
    data_train = train_data,
    data_test = test_data
  )

  return(result_list)
}
