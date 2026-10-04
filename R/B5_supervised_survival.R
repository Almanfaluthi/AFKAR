#' @title Machine Learning Survival Analysis (Prognostic Forecaster)
#' @description Zero-dependency discrete-time survival engine to predict exact time-to-event
#' (e.g., cancer relapse) using a Poisson-Hazard exponential model. Calculates native Kaplan-Meier.
#' Integrates strict clinical boundary checks (CI crossing 1.0) for Hazard Ratio significance.
#' @param data A data.frame containing clinical follow-up data.
#' @param time_var Character. The name of the time-to-event numeric variable (e.g., months).
#' @param status_var Character. The name of the event status variable (1 = Event/Relapse, 0 = Censored).
#' @param split_ratio Numeric. Proportion of data for training (default 0.7).
#' @param reverse_palette Logical. If TRUE, reverses the plot color palette.
#' @param save_plot Logical. If TRUE, exports 4K high-resolution plots.
#' @param save_prefix Character. Custom prefix for exported plot filenames.
#' @param footer Logical. If TRUE, adds reproducibility footer to plots.
#' @param seed Numeric. Seed for reproducibility.
#' @return A list containing Hazard Ratios, C-Index metrics, and predicted median survival times.
#' @export
#' @examples
#' # Run Native Survival Engine
#' result <- B5_supervised_survival(
#'   data = df_B5,
#'   time_var = "time_to_relapse",
#'   status_var = "relapse_status",
#'   save_plot = FALSE
#' )
B5_supervised_survival <- function(data,
                                   time_var,
                                   status_var,
                                   split_ratio = 0.7,
                                   reverse_palette = FALSE,
                                   save_plot = FALSE,
                                   save_prefix = "",
                                   footer = TRUE,
                                   seed = 42) {

  # 1. Input Validation
  stopifnot("Input 'data' must be a data.frame." = is.data.frame(data))
  stopifnot("Time variable not found." = time_var %in% names(data))
  stopifnot("Status variable not found." = status_var %in% names(data))

  set.seed(seed)

  # Check Status Binary & Time > 0
  if (!all(unique(stats::na.omit(data[[status_var]])) %in% c(0, 1))) {
    stop("Status variable must be binary (1 = Event, 0 = Censored).")
  }
  if (any(stats::na.omit(data[[time_var]]) <= 0)) {
    data <- data[data[[time_var]] > 0, ]
    warning("Observations with time <= 0 removed.")
  }

  data <- data[stats::complete.cases(data), ]
  n_total <- nrow(data)

  # Inject log-time uniformly prior to splitting to ensure availability
  data$log_t <- log(data[[time_var]])

  # 2. Train/Test Split
  train_idx <- sample(seq_len(n_total), size = floor(n_total * split_ratio))
  train_data <- data[train_idx, ]
  test_data <- data[-train_idx, ]

  # 3. Base R Poisson-Hazard Modeling (Cox Approximation)
  predictors <- setdiff(names(train_data), c(time_var, status_var, "log_t"))

  form_str <- paste(status_var, "~", paste(predictors, collapse = " + "), "+ offset(log_t)")
  formula_obj <- stats::as.formula(form_str)

  model <- suppressWarnings(stats::glm(formula_obj, data = train_data, family = stats::poisson(link = "log")))

  # 4. Prediction Engine (Exact Time-to-Event)
  # CLINICAL MATH FIX: We must neutralize the time offset (log_t = 0) to extract pure log-hazard rate (X * Beta)
  test_data_pred <- test_data
  test_data_pred$log_t <- 0

  test_log_lambda <- stats::predict(model, newdata = test_data_pred)
  test_lambda <- exp(test_log_lambda)

  # Median Survival Time = -ln(0.5) / Lambda (for exponential distribution framework)
  pred_median_time <- -log(0.5) / test_lambda

  actual_time <- test_data[[time_var]]
  actual_status <- test_data[[status_var]]

  # 5. Base R Native C-Index (Concordance) Calculation
  n_test <- length(actual_time)
  concordant <- 0
  discordant <- 0

  for (i in seq_len(n_test - 1)) {
    for (j in (i + 1):n_test) {
      if (actual_status[i] == 1 && actual_time[i] < actual_time[j]) {
        if (pred_median_time[i] < pred_median_time[j]) concordant <- concordant + 1 else discordant <- discordant + 1
      } else if (actual_status[j] == 1 && actual_time[j] < actual_time[i]) {
        if (pred_median_time[j] < pred_median_time[i]) concordant <- concordant + 1 else discordant <- discordant + 1
      }
    }
  }
  c_index <- ifelse((concordant + discordant) == 0, 0.5, concordant / (concordant + discordant))

  # 6. Strict Clinical Boundary Check for Hazard Ratios
  coef_vals <- stats::coef(model)[predictors] # Exclude intercept
  ci_vals <- suppressMessages(stats::confint.default(model))[predictors, , drop = FALSE]

  hr_data <- data.frame(
    Predictor = predictors,
    HR = round(exp(coef_vals), 2),
    CI_Lower = round(exp(ci_vals[, 1]), 2),
    CI_Upper = round(exp(ci_vals[, 2]), 2),
    stringsAsFactors = FALSE
  )

  # Core Logic: If CI crosses 1.00, force label "ns"
  hr_data$Significance <- ifelse(hr_data$CI_Lower <= 1.00 & hr_data$CI_Upper >= 1.00, "ns", "sig")

  # 7. Native Kaplan-Meier Calculation
  ord <- order(train_data[[time_var]])
  t_sorted <- train_data[[time_var]][ord]
  s_sorted <- train_data[[status_var]][ord]

  n_at_risk <- length(t_sorted):1
  haz_step <- s_sorted / n_at_risk
  surv_prob <- cumprod(1 - haz_step)
  km_df <- data.frame(time = c(0, t_sorted), surv = c(1, surv_prob))

  # 8. Visualization Setup
  base_colors <- c("#2C3E50", "#E74C3C", "#18BC9C", "#F39C12") # 2 = Red (Hazardous), 3 = Green (Protective)
  if (reverse_palette) base_colors <- rev(base_colors)
  footer_text <- paste0("Generated by R-Studio (", R.version.string, ") on ", format(Sys.time(), "%B %d, %Y at %H:%M:%S"))
  timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")

  # ==========================================
  # PLOT 1: Survival Topology & Time Map
  # ==========================================
  graphics::par(mar = c(8, 5, 4, 2) + 0.1, mfrow = c(1, 2))

  # Subplot 1A: Native Kaplan-Meier Curve
  graphics::plot(km_df$time, km_df$surv, type = "s", col = base_colors[1], lwd = 2,
                 main = "Baseline (Kaplan-Meier)",
                 xlab = "", ylab = "Probabilitas Survival",
                 ylim = c(0, 1), bty = "n")

  censor_idx <- which(s_sorted == 0)
  if (length(censor_idx) > 0) {
    graphics::points(t_sorted[censor_idx], surv_prob[censor_idx], pch = 3, col = base_colors[4], cex = 0.6)
  }

  graphics::legend("bottom", legend = c("Survival Curve", "Censored Event"),
                   col = c(base_colors[1], base_colors[4]), lty = c(1, NA), pch = c(NA, 3),
                   bty = "n", horiz = TRUE, inset = c(0, -0.3), xpd = TRUE, cex = 1)

  # Subplot 1B: Predicted vs Actual Event Time (Scatter for Event = 1 only)
  event_idx <- which(actual_status == 1)

  if (length(event_idx) > 0) {
    graphics::plot(actual_time[event_idx], pred_median_time[event_idx], pch = 16,
                   col = grDevices::adjustcolor(base_colors[2], alpha.f = 0.6),
                   main = "Relaps Precision (Test Set)",
                   xlab = "", ylab = "Prediksi Waktu (Bulan/Hari)",
                   xlim = c(0, max(actual_time)), ylim = c(0, max(pred_median_time)), bty = "n")
    graphics::abline(a = 0, b = 1, col = base_colors[1], lwd = 2, lty = 2)

    graphics::legend("bottom", legend = c("Relaps Cases", "Ideal Line(1:1)"),
                     col = c(grDevices::adjustcolor(base_colors[2], alpha.f = 0.6), base_colors[1]),
                     pch = c(16, NA), lty = c(NA, 2), lwd = 2,
                     bty = "n", horiz = TRUE, inset = c(0, -0.3), xpd = TRUE, cex = 1)
  } else {
    graphics::plot(1, type = "n", main = "No Events in Test Set", xlab = "", ylab = "", bty = "n")
  }

  if (footer) graphics::mtext(footer_text, side = 1, line = 6.5, adj = 0, cex = 0.7, col = "grey40", outer = FALSE)

  if (save_plot) {
    file_name1 <- paste0("B5_Topology_", timestamp_str, "_", save_prefix, ".png")
    grDevices::dev.copy(grDevices::png, filename = file_name1, width = 3840, height = 2160, res = 300)
    grDevices::dev.off()
    message("Saved: ", file_name1)
  }

  # ==========================================
  # PLOT 2: Diagnostic Performance & Prognostic Forest
  # ==========================================
  graphics::par(mfrow = c(1, 2), mar = c(7, 4, 4, 2) + 0.1)

  # Subplot 2A: Hazard Density (High Risk vs Low Risk)
  median_hr <- stats::median(test_lambda)
  risk_group <- ifelse(test_lambda > median_hr, "High Risk", "Low Risk")
  dens_high <- stats::density(pred_median_time[risk_group == "High Risk"])
  dens_low <- stats::density(pred_median_time[risk_group == "Low Risk"])

  graphics::plot(dens_high, main = paste0("Model C-Index: ", round(c_index, 3)),
                 xlab = "Relaps Prediction Distribution", ylab = "Densitas",
                 col = base_colors[2], lwd = 2, bty = "n",
                 xlim = c(0, max(c(dens_high$x, dens_low$x))),
                 ylim = c(0, max(c(dens_high$y, dens_low$y)) * 1.1))
  graphics::polygon(dens_high, col = grDevices::adjustcolor(base_colors[2], alpha.f = 0.3), border = NA)
  graphics::lines(dens_low, col = base_colors[3], lwd = 2)
  graphics::polygon(dens_low, col = grDevices::adjustcolor(base_colors[3], alpha.f = 0.3), border = NA)

  graphics::legend("topright", legend = c("High Risk", "Low Risk"),
                   fill = c(grDevices::adjustcolor(base_colors[2], alpha.f = 0.3),
                            grDevices::adjustcolor(base_colors[3], alpha.f = 0.3)),
                   bty = "n", cex = 0.9)

  # Subplot 2B: Pillar-Aligned Forest Plot (Directional Hazard Ratio)
  n_vars <- nrow(hr_data)
  y_pos <- n_vars:1

  point_colors <- ifelse(hr_data$Significance == "ns", "#B0B0B0",
                         ifelse(hr_data$HR > 1.0, base_colors[2], base_colors[3]))
  font_weights <- ifelse(hr_data$Significance == "ns", 1, 2)

  # PILLAR-ALIGNMENT FIX: Mencegah teks bertabrakan dengan CI ekstrim
  x_range <- max(hr_data$CI_Upper) - min(hr_data$CI_Lower)
  if (x_range == 0) x_range <- 1

  x_min_plot <- min(min(hr_data$CI_Lower), 0.5) - (x_range * 0.6)
  x_max_plot <- max(max(hr_data$CI_Upper), 1.5) + (x_range * 0.5)

  graphics::plot(hr_data$HR, y_pos, pch = 16, cex = 1.5, col = point_colors,
                 xlim = c(x_min_plot, x_max_plot), ylim = c(0.5, n_vars + 1),
                 yaxt = "n", xlab = "Hazard Ratio (95% CI)", ylab = "",
                 main = "Prognostic Impacts\n(Pillar-Aligned HR Boundary)", bty = "n")
  graphics::abline(v = 1.0, lty = 2, col = "black", lwd = 1.5)

  for (i in 1:n_vars) {
    graphics::segments(x0 = hr_data$CI_Lower[i], y0 = y_pos[i],
                       x1 = hr_data$CI_Upper[i], y1 = y_pos[i],
                       col = point_colors[i], lwd = 2)

    # Nama Variabel (Rata Kiri Ekstrim di x_min_plot)
    graphics::text(x = x_min_plot, y = y_pos[i],
                   labels = hr_data$Predictor[i], pos = 4, xpd = TRUE,
                   col = point_colors[i], font = font_weights[i], cex = 0.9)

    # Teks Nilai HR & CI (Rata Kanan Ekstrim di x_max_plot)
    graphics::text(x = x_max_plot, y = y_pos[i],
                   labels = sprintf("%.2f [%.2f, %.2f]", hr_data$HR[i], hr_data$CI_Lower[i], hr_data$CI_Upper[i]),
                   pos = 2, cex = 0.8, col = point_colors[i], xpd = TRUE)
  }

  if (footer) graphics::mtext(footer_text, side = 1, line = 5.5, adj = 0.5, cex = 0.7, col = "grey40")

  if (save_plot) {
    file_name2 <- paste0("B5_Diagnostics_", timestamp_str, "_", save_prefix, ".png")
    grDevices::dev.copy(grDevices::png, filename = file_name2, width = 3840, height = 2160, res = 300)
    grDevices::dev.off()
    message("Saved: ", file_name2)
  }

  graphics::par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

  # 9. Return Objects
  results <- list(
    Metrics = data.frame(
      C_Index = c_index
    ),
    Strict_HR_Table = hr_data,
    Model = model
  )

  return(invisible(results))
}
