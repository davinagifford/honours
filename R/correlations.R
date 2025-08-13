


# use the lm_fit_north and lm_fit_south from the calc_cci_south.r and calc_cci_north.r files

# Extract residuals
residuals_north <- resid(lm_fit_north)

residuals_south <- resid(lm_fit_south)

# Correlate the residuals
residuals_correlation <- cor(residuals_north, residuals_south, method = "pearson")
residuals_cor_test <- cor.test(residuals_north, residuals_south, method = "pearson")

# Output results
print(paste("Pearson correlation of residuals:", round(residuals_correlation, 3)))
print(residuals_cor_test)

# Get fitted values from both models
fitted_north <- fitted(lm_fit_north)   # Predicted North from South
fitted_south <- fitted(lm_fit_south) # Predicted South from North


# Run Pearson correlation
fitted_correlation <- cor(fitted_south, fitted_north, method = "pearson")
fitted_cor_test <- cor.test(fitted_south, fitted_north, method = "pearson")

# Output results
print(paste("Pearson correlation of fitted values:", round(fitted_correlation, 3)))
print(fitted_cor_test)
