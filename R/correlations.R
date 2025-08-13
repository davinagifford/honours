


# use the lm_fit_north and lm_fit_south from the calc_cci_south.r and calc_cci_north.r files

# Extract residuals
residuals_north <- resid(lm_fit_north)

residuals_south <- resid(lm_fit_south)

# Compare means and variances
mean(residuals_north); mean(residuals_south)
var(residuals_north); var(residuals_south)

# Plot distributions
hist(residuals_north, col="blue", main="Residual North", xlab="Residual Value")
hist(residuals_south, col="red", main="Residual South", xlab="Residual Value")


# Get fitted values from both models
fitted_north <- fitted(lm_fit_north)   # Predicted North from South
fitted_south <- fitted(lm_fit_south) # Predicted South from North


# Compare means and variances
mean(fitted_north); mean(fitted_south)
var(fitted_north); var(fitted_south)

# Plot distributions
hist(fitted_north, col="blue", main="Fitted North", xlab="Fitted Value")
hist(fitted_south, col="red", main="Fitted South", xlab="Fitted Value")


boxplot(fitted_north, fitted_south, names = c("North", "South"), main = "Fitted Value Comparison")
