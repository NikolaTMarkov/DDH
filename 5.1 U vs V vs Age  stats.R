#───────────────────────────────────────────────────────────────
# Statistical Test: Regress U and V against Age with FDR Correction
#───────────────────────────────────────────────────────────────

# Create the data frame with canonical variates and age
cv_age_df <- data.frame(
  Age = NoddiCog_imputed$Age,
  U = U,  # Imaging canonical variates
  V = V   # Cognitive canonical variates
)

# Rename columns to U1-U7 and V1-V7
colnames(cv_age_df) <- c("Age", paste0("U", 1:ncol(U)), paste0("V", 1:ncol(V)))

# Simple function that works with explicit namespace
run_age_regressions_simple <- function(data) {
  results <- list()
  variate_names <- c(paste0("U", 1:7), paste0("V", 1:7))
  
  for(cv in variate_names) {
    # Use stats::lm to be explicit
    model <- stats::lm(as.formula(paste(cv, "~ Age")), data = data)
    model_summary <- summary(model)
    
    # Extract coefficients safely
    coefs <- model_summary$coefficients
    age_row <- which(rownames(coefs) == "Age")
    
    if(length(age_row) > 0) {
      results[[cv]] <- data.frame(
        Variable = cv,
        Beta = stats::coef(model)[["Age"]],
        SE = coefs[age_row, "Std. Error"],
        t_value = coefs[age_row, "t value"],
        p_value = coefs[age_row, "Pr(>|t|)"],
        R_squared = model_summary$r.squared
      )
    }
  }
  
  return(do.call(rbind, results))
}

# Run the analysis
all_results <- run_age_regressions_simple(cv_age_df)

# Apply FDR correction
all_results$p_adj_fdr <- p.adjust(all_results$p_value, method = "fdr")

# Calculate confidence intervals
all_results$CI_lower <- all_results$Beta - 1.96 * all_results$SE
all_results$CI_upper <- all_results$Beta + 1.96 * all_results$SE

# Add significance stars
all_results$significance <- cut(all_results$p_adj_fdr,
                               breaks = c(0, 0.001, 0.01, 0.05, 1),
                               labels = c("***", "**", "*", ""),
                               include.lowest = TRUE)

# Order by significance
all_results <- all_results[order(all_results$p_adj_fdr), ]

# Format results
formatted_results <- all_results %>%
  mutate(
    `Beta [95% CI]` = paste0(round(Beta, 3), " [", round(CI_lower, 3), ", ", round(CI_upper, 3), "]"),
    `p (FDR)` = ifelse(p_adj_fdr < 0.001, "<0.001", as.character(round(p_adj_fdr, 4)))
  ) %>%
  select(Variable, `Beta [95% CI]`, p_value, `p (FDR)`, R_squared, significance)

# Print results
cat("Age Regression Results for Canonical Variates (FDR-corrected)\n")
cat("───────────────────────────────────────────────────────────────\n")
print(formatted_results, row.names = FALSE)

# Create plot
age_effect_plot <- ggplot(all_results, aes(x = reorder(Variable, -p_adj_fdr), y = Beta, 
                                          color = p_adj_fdr < 0.05)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Age Effects on Canonical Variates",
       x = "Canonical Variate", y = "Beta Coefficient") +
  theme(legend.position = "none")

print(age_effect_plot)

# Summary
cat("\nSignificant variates (FDR < 0.05):\n")
significant <- all_results[all_results$p_adj_fdr < 0.05, ]
if(nrow(significant) > 0) {
  print(significant[, c("Variable", "Beta", "p_adj_fdr")], row.names = FALSE)
} else {
  cat("No significant age effects after FDR correction.\n")
}

cat("\nNumber of significant variates:", sum(all_results$p_adj_fdr < 0.05), "out of", nrow(all_results), "\n")