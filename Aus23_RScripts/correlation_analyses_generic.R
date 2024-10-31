# Load required libraries
library(dplyr)
library(stringr)
library(tidyr)
library(readr)
library(corrplot)
library(ggplot2)

library(writexl)

# Load and prepare datasets
combined_df <- your_dataset #pooled
# Ensure all factor variables are converted to character
combined_df <- combined_df %>% mutate(across(where(is.factor), as.character))


# Filter out rows with missing numeric values
combined_df_cleaned <- combined_df %>%
  filter(across(where(is.numeric), ~ !is.na(.)))

# Define a function to perform correlation analysis
correlation_analysis <- function(data, response_var, predictor_var, filter_conditions = NULL) {
  
  # Apply filter conditions if specified
  if (!is.null(filter_conditions)) {
    data <- data %>% filter(!!rlang::parse_expr(filter_conditions))
  }
  
  # Convert variables to numeric
  data[[response_var]] <- as.numeric(data[[response_var]])
  data[[predictor_var]] <- as.numeric(data[[predictor_var]])
  
  # Check for NA values post-conversion
  if (any(is.na(data[[response_var]])) || any(is.na(data[[predictor_var]]))) {
    warning("Conversion to numeric resulted in NAs, please check your data.")
  }
  
  # Perform correlation test
  cor_test <- cor.test(data[[response_var]], data[[predictor_var]])
  
  return(cor_test)
}

# Define response and predictor variables for the correlation analysis
# Each response var and predictor var will be paired e.g. cor.test(response_var, predictor_var, data = your_dataset) will be run for each possible combination of response and predictor variables
hypotheses <- expand_grid(
  response_var = c("species_richness", "species_richness_native", "species_richness_exotic", 
                   "global_functional_diversity", "mean_pairwise_dist", "PSV", 
                   "shannon_diversity_abundance"),
  predictor_var = c("canopy_cover", "canopy_relief_ratio", "foliar_height_diversity", 
                    "mean_canopy_height", "median_canopy_height", "mean_outer_canopy_height", 
                    "min_canopy_height", "max_canopy_height", "std_outer_canopy_height", 
                    "rh98", "mean_nn_dist", "mean_fhd_grid"),
  filter_conditions = c(NA)
)

# Store correlation results in a list
results <- list()

# Loop through each hypothesis and perform the correlation analysis
for (i in 1:nrow(hypotheses)) {
  result <- correlation_analysis(
    data = combined_df_cleaned,
    response_var = hypotheses$response_var[i],
    predictor_var = hypotheses$predictor_var[i],
    filter_conditions = hypotheses$filter_conditions[i]
  )
  results[[i]] <- result
}

# Create a table of correlation results
result_table <- data.frame(
  ResponseVariable = hypotheses$response_var,
  PredictorVariable = hypotheses$predictor_var,
  CorrelationCoefficient = sapply(results, function(x) x$estimate),
  PValue = sapply(results, function(x) x$p.value)
)

# Display results as a table
kable(result_table, caption = "Correlation Analysis Results")


