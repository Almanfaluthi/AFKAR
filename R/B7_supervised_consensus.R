#' @title Automated Ensemble Benchmark (Multi-Algorithm Consensus)
#' @description Zero-dependency ensemble learning engine mimicking a medical
#' multidisciplinary team. It deploys 5 distinct statistical models (Logit,
#' Probit, C-LogLog, Cauchit, Stepwise AIC) and aggregates their predictions
#' via Majority Voting to create a highly robust 'Super Model'.
#' @param data A data.frame containing clinical predictors.
#' @param target_var Character. The name of the binary outcome variable (0/1).
#' @param split_ratio Numeric. Proportion of data for training (default 0.7).
#' @param reverse_palette Logical. If TRUE, reverses the plot color palette.
#' @param save_plot Logical. If TRUE, exports 4K high-resolution plots.
#' @param save_prefix Character. Custom prefix for exported plot filenames.
#' @param footer Logical. If TRUE, adds reproducibility footer to plots.
#' @param seed Numeric. Seed for reproducibility.
#' @return A list containing individual model performances, Super Model metrics,
#' and comparison tables.
#' @export
#' @examples
#' # Run Multi-Algorithm Consensus
#' result <- B7_supervised_consensus(
#'   data = df_B7,
#'   target_var = "sepsis_status",
#'   save_plot = FALSE
#' )
B7_supervised_consensus <- function(data,
                                    target_var,
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

  # 2. Stratified Train/Test Split
  idx_1 <- which(data[[target_var]] == 1)
  idx_0 <- which(data[[target_var]] == 0)

  if(length(idx_1) < 2 || length(idx_0) < 2) {
    stop("Insufficient occurrences of either Class 0 or Class 1 to perform stratified split.")
  }

  train_1_idx <- sample(idx_1, size = floor(length(idx_1) * split_ratio))
  train_0_idx <- sample(idx_0, size = floor(length(idx_0) * split_ratio))

  train_data <- data[c(train_0_idx, train_1_idx), ]
  test_data <- data[-c(train_0_idx, train_1_idx), ]

  actual <- test_data[[target_var]]

  # 3. Deploy The 5 "Specialists" (Native Base R Algorithms)
  formula_obj <- stats::as.formula(paste(target_var, "~ ."))

  mod_logit   <- suppressWarnings(stats::glm(formula_obj, data = train_data, family = stats::binomial(link = "logit")))
  mod_probit  <- suppressWarnings(stats::glm(formula_obj, data = train_data, family = stats::binomial(link = "probit")))
  mod_cloglog <- suppressWarnings(stats::glm(formula_obj, data = train_data, family = stats::binomial(link = "cloglog")))
  mod_cauchit <- suppressWarnings(stats::glm(formula_obj, data = train_data, family = stats::binomial(link = "cauchit")))
  mod_step    <- suppressWarnings(stats::step(mod_logit, trace = 0, direction = "both"))

  # 4. Independent Predictions on Test Set
  prob_logit   <- stats::predict(mod_logit, test_data, type = "response")
  prob_probit  <- stats::predict(mod_probit, test_data, type = "response")
  prob_cloglog <- stats::predict(mod_cloglog, test_data, type = "response")
  prob_cauchit <- stats::predict(mod_cauchit, test_data, type = "response")
  prob_step    <- stats::predict(mod_step, test_data, type = "response")

  # 5. Super Model Construction via Majority Voting
  votes <- (prob_logit > 0.5) + (prob_probit > 0.5) + (prob_cloglog > 0.5) +
    (prob_cauchit > 0.5) + (prob_step > 0.5)

  super_pred <- ifelse(votes >= 3, 1, 0)
  super_prob <- (prob_logit + prob_probit + prob_cloglog + prob_cauchit + prob_step) / 5

  # 6. Performance Evaluation Helper
  eval_metrics <- function(prob, pred_class, act) {
    tp <- sum(pred_class == 1 & act == 1)
    tn <- sum(pred_class == 0 & act == 0)
    fp <- sum(pred_class == 1 & act == 0)
    fn <- sum(pred_class == 0 & act == 1)

    sens <- ifelse((tp + fn) == 0, 0, tp / (tp + fn))
    spec <- ifelse((tn + fp) == 0, 0, tn / (tn + fp))

    ord <- order(prob, decreasing = TRUE)
    act_ord <- act[ord]
    tpr <- cumsum(act_ord == 1) / sum(act_ord == 1)
    fpr <- cumsum(act_ord == 0) / sum(act_ord == 0)
    tpr <- c(0, tpr)
    fpr <- c(0, fpr)
    auc <- sum(diff(fpr) * (tpr[-1] + tpr[-length(tpr)]) / 2)

    return(list(Sens = sens, Spec = spec, AUC = auc, tpr = tpr, fpr = fpr))
  }

  # 7. Compute Metrics for All Models
  res_logit   <- eval_metrics(prob_logit, ifelse(prob_logit > 0.5, 1, 0), actual)
  res_probit  <- eval_metrics(prob_probit, ifelse(prob_probit > 0.5, 1, 0), actual)
  res_cloglog <- eval_metrics(prob_cloglog, ifelse(prob_cloglog > 0.5, 1, 0), actual)
  res_cauchit <- eval_metrics(prob_cauchit, ifelse(prob_cauchit > 0.5, 1, 0), actual)
  res_step    <- eval_metrics(prob_step, ifelse(prob_step > 0.5, 1, 0), actual)
  res_super   <- eval_metrics(super_prob, super_pred, actual)

  bench_table <- data.frame(
    Algorithm = c("Logit", "Probit", "C-LogLog", "Cauchit", "StepAIC", "SUPER MODEL"),
    Sensitivity = c(res_logit$Sens, res_probit$Sens, res_cloglog$Sens, res_cauchit$Sens, res_step$Sens, res_super$Sens),
    Specificity = c(res_logit$Spec, res_probit$Spec, res_cloglog$Spec, res_cauchit$Spec, res_step$Spec, res_super$Spec),
    AUC = c(res_logit$AUC, res_probit$AUC, res_cloglog$AUC, res_cauchit$AUC, res_step$AUC, res_super$AUC)
  )

  # 8. Visualization Setup
  base_colors <- c("#95A5A6", "#7F8C8D", "#BDC3C7", "#34495E", "#2C3E50")
  super_color <- "#E74C3C"
  if (reverse_palette) {
    base_colors <- rev(base_colors)
    super_color <- "#18BC9C"
  }

  footer_text <- paste0("Generated by R-Studio (", R.version.string, ") on ", format(Sys.time(), "%B %d, %Y at %H:%M:%S"))
  timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")

  # ==========================================
  # PLOT 1: Consensus Topology (ROC & Density)
  # ==========================================
  # OMA (Outer Margin Area) activated: c(bottom=3, left=0, top=0, right=0)
  graphics::par(mfrow = c(1, 2), mar = c(7, 5, 4, 2) + 0.1, oma = c(3, 0, 0, 0))

  # Subplot 1A: ROC Overlay
  graphics::plot(res_logit$fpr, res_logit$tpr, type = "l", col = grDevices::adjustcolor(base_colors[1], alpha.f = 0.5), lwd = 1.5,
                 main = "Multi-Algorithm Consensus (ROC Overlay)",
                 xlab = "False Positive Rate (1 - Specificity)",
                 ylab = "True Positive Rate (Sensitivity)",
                 xlim = c(0, 1), ylim = c(0, 1), bty = "n")

  graphics::lines(res_probit$fpr, res_probit$tpr, col = grDevices::adjustcolor(base_colors[2], alpha.f = 0.5), lwd = 1.5)
  graphics::lines(res_cloglog$fpr, res_cloglog$tpr, col = grDevices::adjustcolor(base_colors[3], alpha.f = 0.5), lwd = 1.5)
  graphics::lines(res_cauchit$fpr, res_cauchit$tpr, col = grDevices::adjustcolor(base_colors[4], alpha.f = 0.5), lwd = 1.5)
  graphics::lines(res_step$fpr, res_step$tpr, col = grDevices::adjustcolor(base_colors[5], alpha.f = 0.5), lwd = 1.5)

  graphics::lines(res_super$fpr, res_super$tpr, col = super_color, lwd = 3.5)
  graphics::abline(a = 0, b = 1, lty = 2, col = "grey60")

  legend_names <- c(paste0("Logit (AUC: ", round(res_logit$AUC, 3), ")"),
                    paste0("Probit (AUC: ", round(res_probit$AUC, 3), ")"),
                    paste0("Super Model (AUC: ", round(res_super$AUC, 3), ")"))

  graphics::legend("bottom", legend = legend_names,
                   col = c(base_colors[1], base_colors[2], super_color),
                   lwd = c(1.5, 1.5, 3.5), bty = "n",
                   horiz = TRUE, inset = c(0, -0.35), xpd = TRUE, cex = 0.85)

  # Subplot 1B: Super Model Probability Density
  dens_0 <- stats::density(super_prob[actual == 0])
  dens_1 <- stats::density(super_prob[actual == 1])

  graphics::plot(dens_0, main = "Super Model Probability Density",
                 xlab = "",
                 ylab = "Density", col = base_colors[4], lwd = 2,
                 xlim = c(0, 1), ylim = c(0, max(c(dens_0$y, dens_1$y)) * 1.2),
                 bty = "n")
  graphics::polygon(dens_0, col = grDevices::adjustcolor(base_colors[4], alpha.f = 0.3), border = NA)
  graphics::lines(dens_1, col = super_color, lwd = 2)
  graphics::polygon(dens_1, col = grDevices::adjustcolor(super_color, alpha.f = 0.3), border = NA)
  graphics::abline(v = 0.5, lty = 2, col = "grey50")

  graphics::legend("bottom", legend = c("Control (Actual)", "Event (Actual)"),
                   fill = c(grDevices::adjustcolor(base_colors[4], alpha.f = 0.3),
                            grDevices::adjustcolor(super_color, alpha.f = 0.3)),
                   bty = "n", horiz = TRUE, inset = c(0, -0.35), xpd = TRUE, cex = 1)

  # ABSOLUTE CENTERED FOOTER: outer = TRUE forces it into OMA zone, adj = 0.5 centers it to canvas
  if (footer) graphics::mtext(footer_text, side = 1, line = 1, outer = TRUE, adj = 0.5, cex = 0.7, col = "grey40")

  # Export Visual 1
  if (save_plot) {
    file_name1 <- paste0("B7_Topology_", timestamp_str, "_", save_prefix, ".png")
    grDevices::dev.copy(grDevices::png, filename = file_name1, width = 3840, height = 2160, res = 300)
    grDevices::dev.off()
    message("Saved: ", file_name1)
  }

  # ==========================================
  # PLOT 2: Comparative Performance Benchmark
  # ==========================================
  graphics::par(mfrow = c(1, 2), mar = c(7, 5, 4, 2) + 0.1, oma = c(3, 0, 0, 0))

  # Subplot 2A: Barplot Benchmark
  mat_plot <- t(as.matrix(bench_table[, c("AUC", "Sensitivity")]))
  colnames(mat_plot) <- bench_table$Algorithm

  bp <- graphics::barplot(mat_plot, beside = TRUE,
                          col = c(base_colors[4], super_color),
                          main = "Automated Benchmark Metrics (Test Set)",
                          ylab = "Performance Score (0.0 - 1.0)",
                          xlab = "",
                          ylim = c(0, 1.25), border = NA,
                          las = 2, cex.names = 0.8)

  graphics::text(x = bp, y = mat_plot + 0.05,
                 labels = sprintf("%.2f", mat_plot),
                 cex = 0.75, col = "black")

  graphics::rect(xleft = bp[1, 6] - 1, ybottom = 0, xright = bp[2, 6] + 1, ytop = 1.25,
                 col = grDevices::adjustcolor(super_color, alpha.f = 0.1), border = NA)

  graphics::legend("bottom", legend = c("AUC Score", "Sensitivity"),
                   fill = c(base_colors[4], super_color), bty = "n",
                   horiz = TRUE, inset = c(0, -0.35), xpd = TRUE, cex = 1)

  # Subplot 2B: Voting Agreement Boxplot
  graphics::boxplot(votes ~ actual, col = c(grDevices::adjustcolor(base_colors[4], 0.5), grDevices::adjustcolor(super_color, 0.5)),
                    main = "Voting Agreement by Clinical Status",
                    ylab = "Number of 'Event' Votes (0 to 5)",
                    xlab = "",
                    names = c("Control (Actual)", "Event (Actual)"),
                    border = c(base_colors[4], super_color),
                    outline = FALSE, bty = "n")

  graphics::points(jitter(actual + 1, factor = 1.2), votes, pch = 16,
                   col = grDevices::adjustcolor("black", alpha.f = 0.15), cex = 0.6)

  graphics::abline(h = 2.5, lty = 2, col = "red", lwd = 1.5)

  graphics::legend("bottom", legend = "Threshold: >= 3 Votes",
                   col = "red", lty = 2, lwd = 1.5, bty = "n",
                   horiz = TRUE, inset = c(0, -0.35), xpd = TRUE, cex = 1)

  # ABSOLUTE CENTERED FOOTER (Plot 2)
  if (footer) graphics::mtext(footer_text, side = 1, line = 1, outer = TRUE, adj = 0.5, cex = 0.7, col = "grey40")

  # Export Visual 2
  if (save_plot) {
    file_name2 <- paste0("B7_Diagnostics_", timestamp_str, "_", save_prefix, ".png")
    grDevices::dev.copy(grDevices::png, filename = file_name2, width = 3840, height = 2160, res = 300)
    grDevices::dev.off()
    message("Saved: ", file_name2)
  }

  # Reset graphic parameters & turn off OMA
  graphics::par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0))

  # 9. Return Unified Object
  results <- list(
    Benchmark_Table = bench_table,
    Super_Model_AUC = res_super$AUC,
    Models = list(
      Logit = mod_logit,
      Probit = mod_probit,
      CLogLog = mod_cloglog,
      Cauchit = mod_cauchit,
      StepAIC = mod_step
    )
  )

  return(invisible(results))
}
