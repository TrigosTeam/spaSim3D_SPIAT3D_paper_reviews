# Code to plot SD1 analysis, showing only cell pair A/B (ignoring cell pairs A/A, B/A, B/B).
# Set seed ----
set.seed(678999821)
# Libraries ----
library(cowplot)
library(ggplot2)
library(S4Vectors)
library(stringr)
library(dplyr)

# Read data and set up ----
setwd("~/R/SD1_data")

metric_df_list <- readRDS("SD1_metric_df_list.RDS")

parameters_df <- readRDS("SD1_parameters_df.RDS")

# Functions ----

# Get plot_df by sampling a random slice
get_plot_df_for_random_slice <- function(metric_df_list, 
                                         metrics,
                                         parameters_df) {
  
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change slice column to numeric for safety
    metric_df$slice <- as.numeric(metric_df$slice)
    
    if (metric %in% c("EBSAC", "EBP_AUC")) {
      # For EBSAC and EBP_AUC, assume pair is the same as cell_types for consistency
      metric_df$pair <- gsub(',', '/', metric_df$cell_types)
    }
    else if (metric %in% c("ANE_AUC")) {
      # For ANE_AUC, assume pair is the same as target_cell_type for consistency (as target is of form A,B already)
      metric_df$pair <- gsub(',', '/', metric_df$target)
    }
    else {
      # Add reference-target column
      metric_df$pair <- paste(metric_df$reference, metric_df$target, sep = "/")
    }
    
    metric_df$metric <- metric  
    
    # Select 3D values, when slice == 0, then remove from data frame
    metric_values3D <- metric_df[[metric]][metric_df[["slice"]] == 0]
    metric_df <- metric_df[metric_df[["slice"]] != 0, ]
    
    # Select a random 2D slice
    metric_df <- metric_df %>%
      group_by(simulation, pair) %>%
      slice_sample(n = 1) %>%   # pick one random slice per simulation
      ungroup() %>%
      select(simulation, pair, metric, value2D = .data[[metric]])
    
    # Combine 3D and 2D values
    metric_df[["value3D"]] <- metric_values3D
    
    # Merge metric_df and parameters_df
    parameters_df$simulation <- seq(nrow(parameters_df))
    metric_df <- metric_df %>% left_join(parameters_df[, c("simulation", "arrangement", "shape")], by = "simulation")
    
    # Get structure column
    metric_df$structure <- paste(metric_df$arrangement, metric_df$shape, sep = "_")
    
    metric_df$metric <- metric  
    
    # Subset for A/B pair
    metric_df <- metric_df[metric_df$pair == "A/B", ]
    
    
    # Add metric_df to plot_df
    plot_df <- rbind(plot_df, metric_df[ , c("value3D", "value2D", "pair", "metric", "structure")])
    
  }
  return(plot_df)
}

# Plot 2D vs 3D for each metric and pair for a random slice
plot_2D_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot <- function(plot_df,
                                                                               metrics) {
  
  # Update metrics (remove _AUC from specific metrics)
  metrics <- sub("_AUC$", "", metrics)
  plot_df$metric <- sub("_AUC$", "", plot_df$metric)
  
  # For nicer tick labels
  sci_clean_threshold <- function(x) {
    
    sapply(x, function(v) {
      if (is.na(v)) return("")
      
      if (abs(v) < 1000) {
        return(as.character(v))
      }
      
      # scientific notation with 1 decimal place
      s <- formatC(v, format = "e", digits = 1)   # e.g. "1.5e+03"
      
      # remove "+" in exponent
      s <- gsub("e\\+", "e", s)
      
      # split mantissa and exponent
      parts <- strsplit(s, "e")[[1]]
      mant <- parts[1]
      exp  <- parts[2]
      
      # remove trailing .0 (so 1.0e3 → 1e3)
      mant <- sub("\\.0$", "", mant)
      
      # remove leading zeros in exponent
      exp <- sub("^0+", "", exp)
      
      paste0(mant, "e", exp)
    })
  }
  
  # Factor metrics
  plot_df$metric <- factor(plot_df$metric, metrics)
  
  # Get correlation
  corr_df <- plot_df %>% 
    group_by(metric) %>% 
    summarise(
      corr = cor(value3D, value2D, method = "spearman", use = "complete.obs"), 
      .groups = "drop")
  
  write.csv(corr_df, "~/R/values_from_figures/random_slice_corr_df.csv")
  
  fig <- ggplot(plot_df, aes(x = value3D, y = value2D)) +
    geom_point(alpha = 0.25, color = "#0062c5", size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                color = "#bb0036", linewidth = 1) +
    # labs(
    #   title = "Scatterplots showing 2D vs 3D, for each metric, a random slice and cell pair A/B",
    #   x = "3D value",
    #   y = "2D value"
    # ) +
    labs(
      x = "",
      y = ""
    ) +
    
    facet_wrap(~ interaction(metric), scales = "free", ncol = length(metrics)) +
    
    scale_x_continuous(n.breaks = 3, labels = sci_clean_threshold) + 
    scale_y_continuous(n.breaks = 3, labels = sci_clean_threshold) +
    
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      
      # Axis tick labels
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      
      # Axis titles
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      
      # Plot title
      plot.title = element_text(size = 18),
      
      # Facet strip titles
      strip.text = element_text(size = 15)
    ) +
    geom_text(
      data = corr_df,
      aes(x = -Inf, y = Inf, label = paste0("r: ", round(corr, 3))),
      inherit.aes = FALSE,
      hjust = -0.1, 
      vjust = 1.4,
      size = 5,   # r-value font size
      color = "black"
    )
  
  
  return(fig)
}

# Plot percentage difference vs 3D for each metric and pair for a random slice
plot_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot <- function(plot_df,
                                                                                                  metrics) {
  
  # Update metrics (remove _AUC from specific metrics)
  metrics <- sub("_AUC$", "", metrics)
  plot_df$metric <- sub("_AUC$", "", plot_df$metric)
  
  # For nicer tick labels
  sci_clean_threshold <- function(x) {
    
    x[!(x %in% range(x, na.rm = T))] <- NA
    
    sapply(x, function(v) {
      if (is.na(v)) return("")
      
      if (abs(v) < 1000) {
        return(as.character(v))
      }
      
      # scientific notation with 1 decimal place
      s <- formatC(v, format = "e", digits = 1)   # e.g. "1.5e+03"
      
      # remove "+" in exponent
      s <- gsub("e\\+", "e", s)
      
      # split mantissa and exponent
      parts <- strsplit(s, "e")[[1]]
      mant <- parts[1]
      exp  <- parts[2]
      
      # remove trailing .0 (so 1.0e3 → 1e3)
      mant <- sub("\\.0$", "", mant)
      
      # remove leading zeros in exponent
      exp <- sub("^0+", "", exp)
      
      paste0(mant, "e", exp)
    })
  }
  
  # Factor for metric
  plot_df$metric <- factor(plot_df$metric, metrics)
  
  # Get error column
  plot_df$error <- (plot_df$value2D - plot_df$value3D) / plot_df$value3D * 100
  
  fig <- ggplot(plot_df, aes(x = value3D, y = error)) +
    geom_point(alpha = 0.25, color = "#0062c5", size = 1) +
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) +
    labs(
      x = "3D value",
      y = "Percentage difference (%)"
    ) +
    facet_wrap(~ interaction(metric), scales = "free", ncol = length(metrics)) +
    scale_x_continuous(n.breaks = 5, labels = sci_clean_threshold) +
    scale_y_continuous(limits = c(-1000, 1000), n.breaks = 5, labels = sci_clean_threshold) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x  = element_text(size = 11),
      axis.text.y  = element_text(size = 11),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      strip.text   = element_text(size = 11)
    )
  
  
  return(fig)
}

# Plot normalised percentage difference vs 3D for each metric and pair for a random slice
plot_normalised_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot <- function(plot_df,
                                                                                                             metrics) {
  
  # Update metrics (remove _AUC from specific metrics)
  metrics <- sub("_AUC$", "", metrics)
  plot_df$metric <- sub("_AUC$", "", plot_df$metric)
  
  # For nicer tick labels
  sci_clean_threshold <- function(x) {
    
    x[!(x %in% range(x, na.rm = T))] <- NA
    
    sapply(x, function(v) {
      if (is.na(v)) return("")
      
      if (abs(v) < 1000) {
        return(as.character(v))
      }
      
      # scientific notation with 1 decimal place
      s <- formatC(v, format = "e", digits = 1)   # e.g. "1.5e+03"
      
      # remove "+" in exponent
      s <- gsub("e\\+", "e", s)
      
      # split mantissa and exponent
      parts <- strsplit(s, "e")[[1]]
      mant <- parts[1]
      exp  <- parts[2]
      
      # remove trailing .0 (so 1.0e3 → 1e3)
      mant <- sub("\\.0$", "", mant)
      
      # remove leading zeros in exponent
      exp <- sub("^0+", "", exp)
      
      paste0(mant, "e", exp)
    })
  }
  
  # Factor for metric
  plot_df$metric <- factor(plot_df$metric, metrics)
  
  # Get error column
  plot_df <- plot_df %>%
    group_by(metric) %>%
    mutate(
      error = (value2D - value3D) / mean(value3D, na.rm = TRUE) * 100
    ) %>%
    ungroup()
  
  fig <- ggplot(plot_df, aes(x = value3D, y = error)) +
    geom_point(alpha = 0.25, color = "#0062c5", size = 1) +
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) +
    labs(
      x = "3D value",
      y = "Normalised percentage difference (%)"
    ) +
    facet_wrap(~ interaction(metric), scales = "free", ncol = length(metrics)) +
    scale_x_continuous(n.breaks = 5, labels = sci_clean_threshold) +
    scale_y_continuous(n.breaks = 5, labels = sci_clean_threshold) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x  = element_text(size = 11),
      axis.text.y  = element_text(size = 11),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      strip.text   = element_text(size = 11)
    )
  
  
  return(fig)
}

# Plot SMAPE vs 3D for each metric and pair for a random slice
plot_SMAPE_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot <- function(plot_df,
                                                                                  metrics) {
  
  # Update metrics (remove _AUC from specific metrics)
  metrics <- sub("_AUC$", "", metrics)
  plot_df$metric <- sub("_AUC$", "", plot_df$metric)
  
  # For nicer tick labels
  sci_clean_threshold <- function(x) {
    
    x[!(x %in% range(x, na.rm = T))] <- NA
    
    sapply(x, function(v) {
      if (is.na(v)) return("")
      
      if (abs(v) < 1000) {
        return(as.character(v))
      }
      
      # scientific notation with 1 decimal place
      s <- formatC(v, format = "e", digits = 1)   # e.g. "1.5e+03"
      
      # remove "+" in exponent
      s <- gsub("e\\+", "e", s)
      
      # split mantissa and exponent
      parts <- strsplit(s, "e")[[1]]
      mant <- parts[1]
      exp  <- parts[2]
      
      # remove trailing .0 (so 1.0e3 → 1e3)
      mant <- sub("\\.0$", "", mant)
      
      # remove leading zeros in exponent
      exp <- sub("^0+", "", exp)
      
      paste0(mant, "e", exp)
    })
  }
  
  # Factor for metric
  plot_df$metric <- factor(plot_df$metric, metrics)
  
  # Get error column
  plot_df$error <- abs(plot_df$value2D - plot_df$value3D) / ((abs(plot_df$value3D) + abs(plot_df$value2D)) / 2) * 100
  
  fig <- ggplot(plot_df, aes(x = value3D, y = error)) +
    geom_point(alpha = 0.25, color = "#0062c5", size = 1) +
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) +
    labs(
      x = "3D value",
      y = "SMAPE (%)"
    ) +
    facet_wrap(~ interaction(metric), scales = "free", ncol = length(metrics)) +
    scale_x_continuous(n.breaks = 5, labels = sci_clean_threshold) +
    scale_y_continuous(limits = c(0, 200), n.breaks = 5, labels = sci_clean_threshold) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x  = element_text(size = 11),
      axis.text.y  = element_text(size = 11),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      strip.text   = element_text(size = 11)
    )
  
  
  return(fig)
}


# Running the functions ----

metrics <- c("AMD",
             "ANC_AUC", "ACIN_AUC", "ANE_AUC",
             "MS_AUC", "NMS_AUC",
             "CK_AUC", "CL_AUC", "CG_AUC",
             "COO_AUC",
             "PBP_AUC", "EBP_AUC", "PBSAC", "EBSAC")



plot_df <- get_plot_df_for_random_slice(metric_df_list,
                                        metrics,
                                        parameters_df)



# 2D vs 3D
setwd("~/R/plots/SD1")
fig_2D_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot <- plot_2D_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot(plot_df,
                                                                                                                                        metrics)
pdf("fig_2D_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot.pdf", width = 30, height = 4)
print(fig_2D_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot)
dev.off()


# Percentage difference vs 3D
setwd("~/R/plots/SD1")
fig_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot <- plot_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot(plot_df,
                                                                                                                                                                              metrics)
pdf("fig_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot.pdf", width = 17.5, height = 2.5)
print(fig_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot)
dev.off()



# Normalised percentage difference vs 3D
setwd("~/R/plots/SD1")
fig_normalised_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot <- plot_normalised_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot(plot_df,
                                                                                                                                                                                                    metrics)
pdf("fig_normalised_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot.pdf", width = 17.5, height = 2.5)
print(fig_normalised_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot)
dev.off()



# SMAPE vs 3D
setwd("~/R/plots/SD1")
fig_SMAPE_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot <- plot_SMAPE_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot(plot_df,
                                                                                                                                              metrics)
pdf("fig_SMAPE_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot.pdf", width = 17.5, height = 2.5)
print(fig_SMAPE_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot)
dev.off()