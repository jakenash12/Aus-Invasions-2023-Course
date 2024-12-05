#BB
#10/23/24
#Looking at how soil stuff changes with the biolog data
#Mostly just biolog stuff in here.

library(ggplot2)
library(tidyr)

theme_set(theme_bw())


ds <- read.csv("Aus23_allData_19Nov24.csv")

ds_long <- ds %>%
  pivot_longer(cols = colnames(ds)[5:36], names_to = "Substrate", values_to = "Biolog")  
#Change the above col numbers to reflect whatever we end up using
#AMF

library(RColorBrewer)

plot_filtered_slopes <- function(data, x_var, y_var, group_var, slope_threshold = 0.2) {
  # Convert input variable names to symbols for dplyr
  x_var <- rlang::sym(x_var)
  y_var <- rlang::sym(y_var)
  group_var <- rlang::sym(group_var)
  
  # Fit linear models for each group and extract the slope
  slope_data <- data %>%
    group_by(!!group_var) %>%
    do(model = lm(!!y_var ~ !!x_var, data = .)) %>%
    mutate(slope = coef(model)[2]) %>%
    # Filter slopes based on the threshold
    filter(abs(slope) > slope_threshold) %>%
    mutate(slope_category = ifelse(slope > 0, "Positive", "Negative"))  # Classify slopes
  
  # Merge the slope data back into the original data
  data_with_slopes <- data %>%
    left_join(slope_data %>% select(!!group_var, slope), by = as.character(rlang::as_name(group_var)))
  
  # Filter the data to include only those with valid slopes
  data_filtered <- data_with_slopes %>%
    filter(!is.na(slope))
  
  # Create the plot
  plot <- ggplot(data_filtered, aes(x = !!x_var, y = !!y_var, color = !!group_var)) +
    geom_point(show.legend = TRUE) +  # Points colored by Substrate
    geom_smooth(
      method = "lm",
      se = FALSE,
      aes(group = !!group_var, fill = slope),  # Fill by slope for the smooth line
      linewidth = 0.5
    ) +
    scale_color_viridis_d(name = "Substrate") +  # Use viridis for substrates
    scale_fill_gradient2(
      low = "red",     # Color for negative slopes
      mid = "grey",    # Color for slopes around zero
      high = "blue",   # Color for positive slopes
      midpoint = 0,    # Midpoint for the gradient
      name = NULL      # Remove legend title for slope
    ) +
    labs(x = as.character(x_var), y = as.character(y_var)) +
    guides(fill = "none")  # Remove the legend for slope
  
  # Create a dynamic filename based on the x variable name
  filename <- paste0("Plots/", as.character(x_var), "_biolog.png")
  
  # Save the plot
  ggsave(filename, plot, width = 10, height = 6)  # Adjust width and height as needed
  
  return(plot)  # Optionally return the plot
}

plot_filtered_slopes(ds_long, "perc_P", "Biolog", "Substrate", slope_threshold = 32)
#Need to togglethe slope_threshold depending on the values of the plot

colnames(ds_long)


plot_filtered_slopes <- function(data, x_vars, y_var, group_var) {
  # Convert input variable names to symbols for dplyr
  y_var <- rlang::sym(y_var)
  group_var <- rlang::sym(group_var)
  
  # Iterate over each x variable
  plots <- lapply(x_vars, function(x_var) {
    x_var <- rlang::sym(x_var)  # Convert the x variable name to a symbol
    
    # Fit linear models for each group and extract the slope
    slope_data <- data %>%
      group_by(!!group_var) %>%
      do(model = lm(!!y_var ~ !!x_var, data = .)) %>%
      mutate(slope = coef(model)[2]) %>%
      ungroup() %>%
      # Get the 2 largest and 2 smallest slopes
      filter(slope %in% c(sort(unique(slope), decreasing = TRUE)[1:2], sort(unique(slope))[1:2])) %>%
      mutate(slope_category = ifelse(slope > 0, "Positive", "Negative"))  # Classify slopes
    
    # Merge the slope data back into the original data
    data_with_slopes <- data %>%
      left_join(slope_data %>% select(!!group_var, slope), by = as.character(rlang::as_name(group_var)))
    
    # Filter the data to include only those with valid slopes
    data_filtered <- data_with_slopes %>%
      filter(!is.na(slope))
    
    # Create the plot
    plot <- ggplot(data_filtered, aes(x = !!x_var, y = !!y_var, color = !!group_var)) +
      geom_point(show.legend = TRUE) +  # Points colored by Substrate
      geom_smooth(
        method = "lm",
        se = FALSE,
        aes(group = !!group_var, fill = slope),  # Fill by slope for the smooth line
        linewidth = 0.5
      ) +
      scale_color_viridis_d(name = "Substrate") +  # Use viridis for substrates
      scale_fill_gradient2(
        low = "red",     # Color for negative slopes
        mid = "grey",    # Color for slopes around zero
        high = "blue",   # Color for positive slopes
        midpoint = 0,    # Midpoint for the gradient
        name = NULL      # Remove legend title for slope
      ) +
      labs(x = as.character(x_var), y = as.character(y_var)) +
      guides(fill = "none") +  # Remove the legend for slope
      facet_wrap(~ tree_spp) +  # Facet by tree_spp
      theme_bw()  # Optional theme for better appearance
    
    # Create a dynamic filename based on the x variable name
    filename <- paste0("Plots/", as.character(x_var), "_biolog.png")
    
    # Save the plot
    ggsave(filename, plot, width = 10, height = 6)  # Adjust width and height as needed
    
    return(plot)  # Return the plot
  })
  
  return(plots)  # Return a list of plots
}

column_names <- colnames(ds_long)[7:73]  # Get column names from 7 to 73

plots <- plot_filtered_slopes(ds_long, column_names, "Biolog", "Substrate")
