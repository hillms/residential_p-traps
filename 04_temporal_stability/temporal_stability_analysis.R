################################################################################
# Temporal Stability and Predictive Modeling Analysis
# Description: Analyzes community stability over time using multiple metrics
#              including CV, MSSD, autocorrelation, GAMs, and ARIMA models
################################################################################

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(mgcv)      # For GAMs (v1.9)
library(forecast)  # For ARIMA/autocorrelation (v8.21)

################################################################################
# LOAD DATA
################################################################################

# Load beta diversity from baseline (created in beta_diversity.R)
load("beta_div_graph_dfs.RData")

# The beta_div_dfs object should contain:
# - sink_id: Sink A or B
# - time_point: Week number
# - distances: Bray-Curtis dissimilarity from baseline (week 1)

cat("========================================\n")
cat("Temporal Stability Analysis\n")
cat("Sink Microbiome Communities\n")
cat("========================================\n\n")

# Clean up sink_id if needed
beta_div_dfs$sink_id <- gsub("Sink ", "", beta_div_dfs$sink_id)

# Sort by sink and time
df <- beta_div_dfs %>% arrange(sink_id, time_point)

cat("Data summary:\n")
cat("Total samples:", nrow(df), "\n")
print(table(df$sink_id))

################################################################################
# 1. STABILITY METRICS
################################################################################

cat("\n========================================\n")
cat("1. STABILITY METRICS\n")
cat("========================================\n\n")

# Calculate multiple stability metrics
# Exclude timepoint 1 (distance = 0 by definition)
stability_metrics <- df %>%
  filter(time_point > 1) %>%
  group_by(sink_id) %>%
  summarise(
    n_timepoints = n(),
    mean_beta = mean(distances),
    sd_beta = sd(distances),
    cv_beta = sd(distances) / mean(distances) * 100,  # Coefficient of variation
    min_beta = min(distances),
    max_beta = max(distances),
    range_beta = max(distances) - min(distances),
    .groups = "drop"
  )

cat("Stability Metrics (excluding t=1):\n\n")
print(as.data.frame(stability_metrics))

cat("\nInterpretation:\n")
cat("- Lower CV (coefficient of variation) indicates more stable community\n")
cat("- Lower range indicates less fluctuation over time\n\n")

# Mean Squared Successive Difference (MSSD)
# Measures volatility/smoothness of trajectory
mssd_calc <- df %>%
  filter(time_point > 1) %>%
  arrange(sink_id, time_point) %>%
  group_by(sink_id) %>%
  mutate(diff_sq = (distances - lag(distances))^2) %>%
  summarise(
    mssd = mean(diff_sq, na.rm = TRUE),
    .groups = "drop"
  )

cat("Mean Squared Successive Difference (MSSD):\n")
cat("(Lower values = more stable/predictable trajectory)\n\n")
print(as.data.frame(mssd_calc))

################################################################################
# 2. AUTOCORRELATION ANALYSIS
################################################################################

cat("\n========================================\n")
cat("2. AUTOCORRELATION ANALYSIS\n")
cat("========================================\n\n")

# Calculate autocorrelation for each sink
acf_results <- list()

for (sink in c("A", "B")) {
  sink_data <- df %>%
    filter(sink_id == sink, time_point > 1) %>%
    arrange(time_point)
  
  cat("Sink", sink, ":\n")
  
  # Create time series object
  ts_data <- ts(sink_data$distances)
  
  # Calculate ACF
  acf_result <- acf(ts_data, lag.max = 10, plot = FALSE)
  acf_results[[sink]] <- acf_result
  
  cat("  Autocorrelation coefficients (lag 1-5):\n")
  for (i in 2:min(6, length(acf_result$acf))) {
    cat(sprintf("    Lag %d: %.3f\n", i-1, acf_result$acf[i]))
  }
  
  # Ljung-Box test for autocorrelation
  lb_test <- Box.test(ts_data, lag = 5, type = "Ljung-Box")
  cat(sprintf("  Ljung-Box test (lag 5): χ² = %.2f, p = %.4f\n",
              lb_test$statistic, lb_test$p.value))
  cat("  Interpretation:", ifelse(lb_test$p.value < 0.05,
                                   "Significant autocorrelation (predictable)",
                                   "No significant autocorrelation"), "\n\n")
}

# Save ACF plots
pdf("autocorrelation_plots.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))

for (sink in c("A", "B")) {
  sink_data <- df %>%
    filter(sink_id == sink, time_point > 1) %>%
    arrange(time_point)
  
  ts_data <- ts(sink_data$distances)
  acf(ts_data, main = paste("Sink", sink, "- Autocorrelation"), lag.max = 15)
}

dev.off()

################################################################################
# 3. GENERALIZED ADDITIVE MODELS (GAMs)
################################################################################

cat("\n========================================\n")
cat("3. GENERALIZED ADDITIVE MODELS (GAMs)\n")
cat("========================================\n\n")

# Fit GAM for each sink
gam_results <- list()

for (sink in c("A", "B")) {
  sink_data <- df %>%
    filter(sink_id == sink) %>%
    arrange(time_point)
  
  # Fit GAM with smooth term for time
  # REML: Restricted Maximum Likelihood estimation
  # k = 10: Maximum basis dimension for smooth
  gam_model <- gam(distances ~ s(time_point, k = 10), 
                   data = sink_data, 
                   method = "REML")
  
  gam_results[[sink]] <- gam_model
  
  cat("Sink", sink, "GAM Summary:\n")
  cat("  Deviance explained:", round(summary(gam_model)$dev.expl * 100, 1), "%\n")
  cat("  Effective df for smooth:", round(summary(gam_model)$s.table[, "edf"], 2), "\n")
  cat("  Smooth term p-value:", format(summary(gam_model)$s.table[, "p-value"], scientific = TRUE), "\n\n")
}

# Compare models
cat("Model Comparison:\n")
cat("  Sink A - AIC:", round(AIC(gam_results$A), 2), "\n")
cat("  Sink B - AIC:", round(AIC(gam_results$B), 2), "\n")
cat("  (Lower AIC = better fit; not directly comparable between sinks)\n\n")

# Create prediction data for plotting
pred_data <- expand.grid(
  time_point = seq(1, max(df$time_point), length.out = 100)
)

# Generate predictions with confidence intervals
gam_predictions <- data.frame()

for (sink in c("A", "B")) {
  pred <- predict(gam_results[[sink]], newdata = pred_data, se.fit = TRUE)
  
  temp_df <- data.frame(
    time_point = pred_data$time_point,
    fit = pred$fit,
    se = pred$se.fit,
    lower = pred$fit - 1.96 * pred$se.fit,
    upper = pred$fit + 1.96 * pred$se.fit,
    sink_id = sink
  )
  
  gam_predictions <- rbind(gam_predictions, temp_df)
}

# Plot GAM results
p_gam <- ggplot() +
  geom_ribbon(data = gam_predictions, 
              aes(x = time_point, ymin = lower, ymax = upper, fill = sink_id), 
              alpha = 0.3) +
  geom_line(data = gam_predictions, 
            aes(x = time_point, y = fit, color = sink_id), 
            linewidth = 1) +
  geom_point(data = df, 
             aes(x = time_point, y = distances, color = sink_id), 
             alpha = 0.5) +
  scale_color_manual(values = c("A" = "#4575b4", "B" = "#2ca25f"), name = "Sink") +
  scale_fill_manual(values = c("A" = "#4575b4", "B" = "#2ca25f"), name = "Sink") +
  labs(
    title = "GAM Smoothed Trends in Beta Diversity Over Time",
    subtitle = "Shaded regions show 95% confidence intervals",
    x = "Time (weeks)",
    y = "Beta Diversity from Baseline"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave("gam_temporal_trends.pdf", p_gam, width = 10, height = 6)

################################################################################
# 4. AUTOREGRESSIVE PREDICTIVE MODELING (ARIMA)
################################################################################

cat("\n========================================\n")
cat("4. AUTOREGRESSIVE PREDICTIVE MODELING\n")
cat("========================================\n\n")

ar_results <- list()

for (sink in c("A", "B")) {
  sink_data <- df %>%
    filter(sink_id == sink, time_point > 1) %>%
    arrange(time_point)
  
  ts_data <- ts(sink_data$distances)
  
  # Fit auto ARIMA
  # Automatic model selection based on AIC
  ar_model <- auto.arima(ts_data, 
                         seasonal = FALSE, 
                         stepwise = FALSE, 
                         approximation = FALSE)
  
  ar_results[[sink]] <- ar_model
  
  cat("Sink", sink, "- Best ARIMA model:", 
      paste0("ARIMA(", paste(arimaorder(ar_model), collapse = ","), ")"), "\n")
  cat("  AIC:", round(AIC(ar_model), 2), "\n")
  cat("  BIC:", round(BIC(ar_model), 2), "\n")
  
  # Model accuracy (in-sample)
  acc <- accuracy(ar_model)
  cat("  RMSE:", round(acc[1, "RMSE"], 4), "\n")
  cat("  MAE:", round(acc[1, "MAE"], 4), "\n\n")
}

################################################################################
# 5. LEAVE-FUTURE-OUT CROSS-VALIDATION
################################################################################

cat("\n========================================\n")
cat("5. LEAVE-FUTURE-OUT CROSS-VALIDATION\n")
cat("========================================\n\n")

cv_results <- data.frame()

for (sink in c("A", "B")) {
  sink_data <- df %>%
    filter(sink_id == sink, time_point > 1) %>%
    arrange(time_point)
  
  n <- nrow(sink_data)
  train_size <- floor(n * 0.7)  # Use 70% for training
  
  errors <- c()
  
  # Rolling window prediction
  for (i in (train_size + 1):n) {
    train_data <- ts(sink_data$distances[1:(i-1)])
    test_value <- sink_data$distances[i]
    
    # Fit AR model on training data
    tryCatch({
      model <- auto.arima(train_data, seasonal = FALSE, max.p = 3, max.q = 3)
      pred <- forecast(model, h = 1)$mean[1]
      errors <- c(errors, (test_value - pred)^2)
    }, error = function(e) {
      # Skip if model fitting fails
    })
  }
  
  rmse <- sqrt(mean(errors, na.rm = TRUE))
  
  cv_results <- rbind(cv_results, data.frame(
    sink_id = sink,
    cv_rmse = rmse,
    n_predictions = length(errors)
  ))
  
  cat("Sink", sink, ":\n")
  cat("  Cross-validation RMSE:", round(rmse, 4), "\n")
  cat("  Number of predictions:", length(errors), "\n\n")
}

cat("Interpretation:\n")
cat("Lower CV-RMSE indicates more predictable community dynamics.\n")

################################################################################
# 6. SUMMARY COMPARISON
################################################################################

cat("\n========================================\n")
cat("6. SUMMARY COMPARISON\n")
cat("========================================\n\n")

summary_table <- stability_metrics %>%
  left_join(mssd_calc, by = "sink_id") %>%
  left_join(cv_results, by = "sink_id") %>%
  select(sink_id, mean_beta, cv_beta, range_beta, mssd, cv_rmse)

colnames(summary_table) <- c("Sink", "Mean Beta Div", "CV (%)", "Range", "MSSD", "CV-RMSE")

cat("Summary Table:\n\n")
print(as.data.frame(summary_table), row.names = FALSE)

cat("\n\nKey Findings:\n")
cat("-------------\n")

# Compare metrics
if (summary_table$`CV (%)`[summary_table$Sink == "A"] < 
    summary_table$`CV (%)`[summary_table$Sink == "B"]) {
  cat("- Sink A has LOWER coefficient of variation (more stable)\n")
} else {
  cat("- Sink B has LOWER coefficient of variation (more stable)\n")
}

if (summary_table$MSSD[summary_table$Sink == "A"] < 
    summary_table$MSSD[summary_table$Sink == "B"]) {
  cat("- Sink A has LOWER MSSD (smoother trajectory)\n")
} else {
  cat("- Sink B has LOWER MSSD (smoother trajectory)\n")
}

if (summary_table$`CV-RMSE`[summary_table$Sink == "A"] < 
    summary_table$`CV-RMSE`[summary_table$Sink == "B"]) {
  cat("- Sink A has LOWER prediction error (more predictable)\n")
} else {
  cat("- Sink B has LOWER prediction error (more predictable)\n")
}

################################################################################
# SAVE RESULTS
################################################################################

save(stability_metrics, mssd_calc, cv_results, gam_results, ar_results, summary_table,
     file = "temporal_stability_results.RData")

write.csv(summary_table, "stability_summary.csv", row.names = FALSE)
write.csv(stability_metrics, "stability_metrics_detailed.csv", row.names = FALSE)

cat("\n========================================\n")
cat("Output files saved:\n")
cat("  - autocorrelation_plots.pdf\n")
cat("  - gam_temporal_trends.pdf\n")
cat("  - stability_summary.csv\n")
cat("  - stability_metrics_detailed.csv\n")
cat("  - temporal_stability_results.RData\n")
cat("========================================\n")

################################################################################
# NOTES:
# - CV (Coefficient of Variation): Relative stability metric
# - MSSD (Mean Squared Successive Difference): Trajectory smoothness
# - ACF (Autocorrelation Function): Temporal dependence structure
# - Ljung-Box test: Formal test for autocorrelation
# - GAMs (Generalized Additive Models): Non-linear temporal trends (mgcv v1.9)
# - ARIMA: Autoregressive integrated moving average models (forecast v8.21)
# - Cross-validation: 70% training, rolling window for out-of-sample prediction
# - All analyses use Bray-Curtis dissimilarity from baseline (week 1)
################################################################################
