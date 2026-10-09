# Code to generate SD1 and to perform SPIAT-3D analyis on SD1.
# Set seed for consistency
set.seed(678999821)

library(dplyr)
library(ggplot2)

# Functions -----
# Generate data frame containing values of important parameters
generate_random_parameters <- function(
    n_simulations,
    bg_prop_A_range,
    bg_prop_B_range,
    E_radius_x_range,
    E_radius_y_range,
    E_radius_z_range,
    N_width_range,
    cluster_prop_A_range,
    ring_width_factor_range,
    cluster1_x_coord_range
) {
  
  # For each parameter, get 'n_simulations' values between specified min and max of each parameter  
  parameters_df <- data.frame(
    bg_prop_A = runif(n_simulations, bg_prop_A_range["min"], bg_prop_A_range["max"]),
    bg_prop_B = runif(n_simulations, bg_prop_B_range["min"], bg_prop_B_range["max"]),
    E_radius_x = runif(n_simulations, E_radius_x_range["min"], E_radius_x_range["max"]),
    E_radius_y = runif(n_simulations, E_radius_y_range["min"], E_radius_y_range["max"]),
    E_radius_z = runif(n_simulations, E_radius_z_range["min"], E_radius_z_range["max"]),
    N_width = runif(n_simulations, N_width_range["min"], N_width_range["max"]),
    cluster_prop_A = runif(n_simulations, cluster_prop_A_range["min"], cluster_prop_A_range["max"]),
    ring_width_factor = runif(n_simulations, ring_width_factor_range["min"], ring_width_factor_range["max"]),
    cluster1_x_coord = runif(n_simulations, cluster1_x_coord_range["min"], cluster1_x_coord_range["max"])
  )
  
  # Randomly choose which arrangement and shape to use
  arrangements <- c("mixed", "ringed", "separated")
  shapes <- c("ellipsoid", "network")
  
  parameters_df$arrangement <- sample(arrangements, size = nrow(parameters_df), replace = TRUE)
  parameters_df$shape <- sample(shapes, size = nrow(parameters_df), replace = TRUE)
  
  return(parameters_df)
}

# Analyse SD1 simulations
analyse_SD1_simulations <- function(parameters_df) {
  
  # Set defined parameters/values
  cell_types <- c('A', 'B')
  n_simulations <- nrow(parameters_df)
  
  radii <- seq(20, 400, 5)
  radii_colnames <- paste("r", radii, sep = "")
  
  n_splits <- 10
  thresholds <- seq(0.01, 1, 0.01)
  thresholds_colnames <- paste("t", thresholds, sep = "")
  
  # Generate simulation metadata using parameters data frame from 'generate_random_parameters' function
  generate_simulation_metadata <- function(parameters_df, simulation_index) {
    
    # Constant metadata for simulations
    simulation_metadata <- spe_metadata_background_template("random")
    simulation_metadata$background$n_cells <- 30000
    simulation_metadata$background$length <- 600
    simulation_metadata$background$width <- 600
    simulation_metadata$background$height <- 300
    simulation_metadata$background$minimum_distance_between_cells <- 10
    simulation_metadata$background$cell_types <- c("A", "B", "O") # Cell proportions will change later
    
    N_radius <- 125
    N_n_edges <- 20
    
    mixed_cluster_cell_types <- c("A", "B")
    mixed_cluster_centre_loc <- c(300, 300, 150)
    
    ringed_cluster_cell_type <- "A"
    ringed_cluster_cell_prop <- 1
    ringed_ring_cell_type <- "B"
    ringed_ring_cell_prop <- 1
    ringed_cluster_centre_loc <- c(300, 300, 150)
    
    separated_cluster1_cell_type <- "A"
    separated_cluster1_cell_prop <- 1
    separated_cluster1_y_z_centre_loc <- c(300, 150) # Append x-coord to it later
    separated_cluster2_shape <- "sphere"
    separated_cluster2_cell_type <- "B"
    separated_cluster2_cell_prop <- 1
    separated_cluster2_centre_loc <- c(450, 300, 150)
    separated_cluster2_radius <- 100
    
    # Alter background cell proportions
    simulation_metadata$background$cell_proportions <- c(parameters_df$bg_prop_A[simulation_index], 
                                                         parameters_df$bg_prop_B[simulation_index],
                                                         1 - parameters_df$bg_prop_A[simulation_index] - parameters_df$bg_prop_B[simulation_index]) # prop(O) = 1 - prop(A) - prop(B)
    
    shape <- parameters_df[simulation_index, "shape"]
    arrangement <- parameters_df[simulation_index, "arrangement"]
    
    # Determine what type of cluster to add to current metadata
    if (arrangement %in% c("mixed", "separated")) {
      simulation_metadata <- spe_metadata_cluster_template("regular", shape, simulation_metadata)
    }
    else if (arrangement == "ringed") {
      simulation_metadata <- spe_metadata_cluster_template("ring", shape, simulation_metadata)
    }
    
    # Modify shape parameters
    if (shape == "ellipsoid") {
      simulation_metadata$cluster_1$radii <- c(parameters_df$E_radius_x[simulation_index], parameters_df$E_radius_y[simulation_index], parameters_df$E_radius_z[simulation_index])
      simulation_metadata$cluster_1$axes_rotation <- c(runif(1, 0, 180), runif(1, 0, 180), runif(1, 0, 180)) # Choose random rotation angles
    }
    
    else if (shape == "network") {
      simulation_metadata$cluster_1$n_edges <- N_n_edges
      simulation_metadata$cluster_1$width <- parameters_df$N_width[simulation_index]
      simulation_metadata$cluster_1$radius <- N_radius
    }
    
    # Modify arrangemetn parameters
    if (arrangement == "mixed") {
      simulation_metadata$cluster_1$cluster_cell_types <- mixed_cluster_cell_types
      simulation_metadata$cluster_1$cluster_cell_proportions <- c(parameters_df$cluster_prop_A[simulation_index], 1 - parameters_df$cluster_prop_A[simulation_index])
      simulation_metadata$cluster_1$centre_loc <- mixed_cluster_centre_loc
    }
    else if (arrangement == "ringed") {
      simulation_metadata$cluster_1$cluster_cell_types <- ringed_cluster_cell_type
      simulation_metadata$cluster_1$cluster_cell_proportions <- ringed_cluster_cell_prop
      simulation_metadata$cluster_1$centre_loc <- ringed_cluster_centre_loc
      simulation_metadata$cluster_1$ring_cell_types <- ringed_ring_cell_type
      simulation_metadata$cluster_1$ring_cell_proportions <- ringed_ring_cell_prop
      simulation_metadata$cluster_1$ring_width <- parameters_df$ring_width_factor[simulation_index] * ifelse(shape == "network", 
                                                                                                             parameters_df$N_width[simulation_index],
                                                                                                             apply(parameters_df[, c("E_radius_x", "E_radius_y", "E_radius_z")], 1, mean))
    }
    else if (arrangement == "separated") {
      simulation_metadata$cluster_1$cluster_cell_types <- separated_cluster1_cell_type
      simulation_metadata$cluster_1$cluster_cell_proportions <- separated_cluster1_cell_prop
      simulation_metadata$cluster_1$centre_loc <- c(parameters_df$cluster1_x_coord[simulation_index], separated_cluster1_y_z_centre_loc)
      
      simulation_metadata <- spe_metadata_cluster_template("regular", "sphere", simulation_metadata)
      simulation_metadata$cluster_2$cluster_cell_types <- separated_cluster2_cell_type
      simulation_metadata$cluster_2$cluster_cell_proportions <- separated_cluster2_cell_prop
      simulation_metadata$cluster_2$centre_loc <- separated_cluster2_centre_loc
      simulation_metadata$cluster_2$radius <- separated_cluster2_radius
    }
    
    return(simulation_metadata)
  }
  
  # Function to create empty metric df list
  create_empty_metric_df_list <- function(
    cell_types,
    n_simulations,
    radii_colnames,
    thresholds_colnames
  ) {
    n_cell_type_combinations <- length(cell_types)^2
    
    # Define AMD data frames as well as constants
    # n_slices + 1, to include the 3D value
    AMD_df_colnames <- c("simulation", "reference", "target", "AMD")
    AMD_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(AMD_df_colnames)))
    colnames(AMD_df) <- AMD_df_colnames
    
    # Define MS, NMS, ANC, ACIN, FCO, ANE data frames as well as constants
    radii_colnames <- paste("r", radii, sep = "")
    
    MS_df_colnames <- c("simulation", "reference", "target", radii_colnames)
    MS_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(MS_df_colnames)))
    colnames(MS_df) <- MS_df_colnames
    
    NMS_df <- ANC_df <- ANE_df <- ACIN_df <- FCO_df <- CK_df <- CL_df <- CG_df <- MS_df
    
    # Define SAC and prevalence data frames as well as constants
    thresholds_colnames <- paste("t", thresholds, sep = "")
    
    PBSAC_df_colnames <- c("simulation", "reference", "target", "PBSAC")
    PBSAC_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(PBSAC_df_colnames)))
    colnames(PBSAC_df) <- PBSAC_df_colnames
    
    PBP_df_colnames <- c("simulation", "reference", "target", thresholds_colnames)
    PBP_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(PBP_df_colnames)))
    colnames(PBP_df) <- PBP_df_colnames
    
    EBSAC_df_colnames <- c("simulation", "cell_types", "EBSAC")
    EBSAC_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(EBSAC_df_colnames)))
    colnames(EBSAC_df) <- EBSAC_df_colnames
    
    EBP_df_colnames <- c("simulation", "cell_types", thresholds_colnames)
    EBP_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(EBP_df_colnames)))
    colnames(EBP_df) <- EBP_df_colnames
    
    
    # Add all to list:
    metric_df_list <- list(AMD = AMD_df,
                           MS = MS_df,
                           NMS = NMS_df,
                           ACIN = ACIN_df,
                           ANE = ANE_df,
                           ANC = ANC_df,
                           FCO = FCO_df,
                           CK = CK_df,
                           CL = CL_df,
                           CG = CG_df,
                           PBSAC = PBSAC_df,
                           PBP = PBP_df,
                           EBSAC = EBSAC_df,
                           EBP = EBP_df)
    
    return(metric_df_list)
  }
  
  # Function analyse spe in 2D
  analyse_simulation2D <- function(spe, 
                                   cell_types, 
                                   radii, 
                                   thresholds, 
                                   n_splits, 
                                   simulation_index, 
                                   metric_df_list) {
    n_cell_type_combinations <- length(cell_types)^2
    radii_colnames <- paste("r", radii, sep = "")
    thresholds_colnames <- paste("t", thresholds, sep = "")
    
    index <- n_cell_type_combinations * (simulation_index - 1) + 1 
    
    minimum_distance_data <- calculate_minimum_distances_between_cell_types2D(spe,
                                                                              cell_types,
                                                                              feature_colname = "Cell.Type",
                                                                              show_summary = F,
                                                                              plot_image = F)
    
    minimum_distance_data_summary <- summarise_distances_between_cell_types2D(minimum_distance_data)
    
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "simulation"] <- simulation_index
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "reference"] <- minimum_distance_data_summary$reference
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "target"] <- minimum_distance_data_summary$target
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "AMD"] <- minimum_distance_data_summary$mean
    
    # Need a new index for 'pair' data which increments after each reference-target pair
    pair_index <- n_cell_type_combinations * (simulation_index - 1) + 1 
    
    for (reference_cell_type in cell_types) {
      gradient_data <- calculate_all_gradient_cc_metrics2D(spe,
                                                           reference_cell_type,
                                                           cell_types,
                                                           radii,
                                                           feature_colname = "Cell.Type")
      for (target_cell_type in cell_types) {
        print(paste(reference_cell_type, target_cell_type, sep = "/"))
        metric_df_list[["ANC"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["ACIN"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["FCO"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["CK"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["CL"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["CG"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["MS"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["NMS"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["ANE"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, 
                                                                                         paste(reference_cell_type, target_cell_type, sep = ","))
        
        
        if (is.null(gradient_data)) {
          metric_df_list[["ANC"]][pair_index, radii_colnames] <- NA
          metric_df_list[["FCO"]][pair_index, radii_colnames] <- NA
          metric_df_list[["CK"]][pair_index, radii_colnames] <- NA
          metric_df_list[["CL"]][pair_index, radii_colnames] <- NA
          metric_df_list[["CG"]][pair_index, radii_colnames] <- NA
          metric_df_list[["ACIN"]][pair_index, radii_colnames] <- NA
          metric_df_list[["MS"]][pair_index, radii_colnames] <- NA
          metric_df_list[["NMS"]][pair_index, radii_colnames] <- NA
          metric_df_list[["ANE"]][pair_index, radii_colnames] <- NA
        }
        else {
          metric_df_list[["ANC"]][pair_index, radii_colnames] <- gradient_data[["neighbourhood_counts"]][[target_cell_type]]
          metric_df_list[["FCO"]][pair_index, radii_colnames] <- gradient_data[["fast_co_occurrence"]][[target_cell_type]]
          metric_df_list[["CK"]][pair_index, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]]
          metric_df_list[["CL"]][pair_index, radii_colnames] <- gradient_data[["cross_L"]][[target_cell_type]]
          metric_df_list[["CG"]][pair_index, radii_colnames] <- gradient_data[["cross_G"]][[target_cell_type]][["observed_cross_G"]]
          
          if (reference_cell_type != target_cell_type) {
            metric_df_list[["ACIN"]][pair_index, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
            metric_df_list[["MS"]][pair_index, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$mixing_score
            metric_df_list[["NMS"]][pair_index, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$normalised_mixing_score
            metric_df_list[["ANE"]][pair_index, radii_colnames] <- gradient_data[["neighbourhood_entropy"]][[target_cell_type]]
          }
        }
        if (reference_cell_type == target_cell_type) {
          metric_df_list[["ACIN"]][pair_index, radii_colnames] <- Inf
          metric_df_list[["MS"]][pair_index, radii_colnames] <- Inf
          metric_df_list[["NMS"]][pair_index, radii_colnames] <- Inf
          metric_df_list[["ANE"]][pair_index, radii_colnames] <- Inf
        }
        
        # Spatial heterogeneity metrics
        metric_df_list[["PBSAC"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["PBP"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        
        metric_df_list[["EBSAC"]][pair_index, c("simulation", "cell_types")] <- c(simulation_index, paste(reference_cell_type, target_cell_type, sep = ","))
        metric_df_list[["EBP"]][pair_index, c("simulation", "cell_types")] <- c(simulation_index, paste(reference_cell_type, target_cell_type, sep = ","))
        
        if (reference_cell_type != target_cell_type) {
          proportion_grid_metrics <- calculate_cell_proportion_grid_metrics2D(spe, 
                                                                              n_splits,
                                                                              reference_cell_type, 
                                                                              target_cell_type,
                                                                              feature_colname = "Cell.Type",
                                                                              plot_image = F)
          
          if (is.null(proportion_grid_metrics)) {
            metric_df_list[["PBSAC"]][pair_index, "PBSAC"] <- NA
            metric_df_list[["PBP"]][pair_index, thresholds_colnames] <- NA
          }
          else {
            PBSAC <- calculate_spatial_autocorrelation2D(proportion_grid_metrics, 
                                                         "proportion",
                                                         weight_method = "queen")
            
            PBP_df <- calculate_prevalence_gradient2D(proportion_grid_metrics,
                                                      "proportion",
                                                      show_AUC = F)
            
            
            metric_df_list[["PBSAC"]][pair_index, "PBSAC"] <- PBSAC
            metric_df_list[["PBP"]][pair_index, thresholds_colnames] <- PBP_df$prevalence
          } 
          
          entropy_grid_metrics <- calculate_entropy_grid_metrics2D(spe, 
                                                                   n_splits,
                                                                   c(reference_cell_type, target_cell_type), 
                                                                   feature_colname = "Cell.Type",
                                                                   plot_image = F)
          
          if (is.null(entropy_grid_metrics)) {
            metric_df_list[["EBSAC"]][pair_index, "EBSAC"] <- NA
            metric_df_list[["EBP"]][pair_index, thresholds_colnames] <- NA
          }
          else {
            EBSAC <- calculate_spatial_autocorrelation2D(entropy_grid_metrics, 
                                                         "entropy",
                                                         weight_method = "queen")
            
            EBP_df <- calculate_prevalence_gradient2D(entropy_grid_metrics,
                                                      "entropy",
                                                      show_AUC = F)
            
            metric_df_list[["EBSAC"]][pair_index, "EBSAC"] <- EBSAC
            metric_df_list[["EBP"]][pair_index, thresholds_colnames] <- EBP_df$prevalence
          }    
        }
        else {
          metric_df_list[["PBSAC"]][pair_index, "PBSAC"] <- Inf
          metric_df_list[["PBP"]][pair_index, thresholds_colnames] <- Inf
          metric_df_list[["EBSAC"]][pair_index, "EBSAC"] <- Inf
          metric_df_list[["EBP"]][pair_index, thresholds_colnames] <- Inf
        }
        
        pair_index <- pair_index + 1
      }
    }
    return(metric_df_list)
  }
  
  # Define metric df lists
  metric_df_list <- create_empty_metric_df_list(cell_types, n_simulations, radii_colnames, thresholds_colnames)
  
  # Define list containing x,y,z,cell type for each spe
  spe_list <- list()
  
  for (simulation_index in seq_len(n_simulations)) {
    print(simulation_index)
    # Simulate spes in 3D, and get 2D slices
    simulation_metadata <- generate_simulation_metadata(parameters_df, simulation_index)
    
    spe3D <- simulate_spe_metadata3D(simulation_metadata, plot_image = F)
    
    spe2D <- spe3D[ , 145 < SpatialExperiment::spatialCoords(spe3D)[ , "Cell.Z.Position"] & SpatialExperiment::spatialCoords(spe3D)[ , "Cell.Z.Position"] < 155]
    
    spe_list[[paste("sim", simulation_index, sep = "_")]] <- data.frame(SpatialExperiment::spatialCoords(spe2D), "Cell.Type" = spe2D[["Cell.Type"]])
    
    # Analyse spes
    metric_df_list <- analyse_simulation2D(spe2D,
                                           cell_types,
                                           radii,
                                           thresholds,
                                           n_splits,
                                           simulation_index,
                                           metric_df_list)
  }
  
  ## Turn gradient metrics into AUC and add to metric_df list
  get_AUC_for_radii_gradient_metrics <- function(y) {
    x <- radii
    h <- diff(x)[1]
    n <- length(x)
    
    AUC <- (h / 2) * (y[1] + 2 * sum(y[2:(n - 1)]) + y[n])
    
    return(AUC)
  }
  
  add_AUC_for_radii_gradient_metrics_to_metric_df_list <- function(metric_df_list) {
    gradient_radii_metrics <- c("MS", "NMS", "ACIN", "ANE", "ANC", "FCO", "CK", "CL", "CG")
    
    for (metric in gradient_radii_metrics) {
      metric_AUC_name <- paste(metric, "AUC", sep = "_")
      
      if (metric %in% c("MS", "NMS", "ANC", "ACIN", "ANE", "FCO", "CK", "CL", "CG")) {
        subset_colnames <- c("simulation", "reference", "target", metric_AUC_name)
      }
      else {
        stop("Unknown metric?")
      }
      
      df <- metric_df_list[[metric]]
      df[[metric_AUC_name]] <- apply(df[ , radii_colnames], 1, get_AUC_for_radii_gradient_metrics)
      df <- df[ , subset_colnames]
      metric_df_list[[metric_AUC_name]] <- df
    }
    
    # PBP_AUC 3D
    PBP_df <- metric_df_list[["PBP"]]
    PBP_df$PBP_AUC <- apply(PBP_df[ , thresholds_colnames], 1, function(y) {
      sum(diff(thresholds) * (head(y, -1) + tail(y, -1)) / 2)
    })
    PBP_AUC_df <- PBP_df[ , c("simulation", "reference", "target", "PBP_AUC")]
    metric_df_list[["PBP_AUC"]] <- PBP_AUC_df
    
    
    # EBP_AUC 3D
    EBP_df <- metric_df_list[["EBP"]]
    EBP_df$EBP_AUC <- apply(EBP_df[ , thresholds_colnames], 1, function(y) {
      sum(diff(thresholds) * (head(y, -1) + tail(y, -1)) / 2)
    })
    EBP_AUC_df <- EBP_df[ , c("simulation", "cell_types", "EBP_AUC")]
    metric_df_list[["EBP_AUC"]] <- EBP_AUC_df
    
    return(metric_df_list)
  }
  
  metric_df_list <- add_AUC_for_radii_gradient_metrics_to_metric_df_list(metric_df_list)
  
  return(list("metric_df_list2D" = metric_df_list, "spe_list" = spe_list))

}

# Running the functions ----
parameters_df <- generate_random_parameters(
  n_simulations = 1,
  bg_prop_A_range = c("min" = 0, "max" = 0.10),
  bg_prop_B_range = c("min" = 0, "max" = 0.10),
  E_radius_x_range = c("min" = 75, "max" = 125),
  E_radius_y_range = c("min" = 75, "max" = 125),
  E_radius_z_range = c("min" = 75, "max" = 125),
  N_width_range = c("min" = 25, "max" = 35),
  cluster_prop_A_range = c("min" = 0.3, "max" = 0.7),
  ring_width_factor_range = c("min" = 0.1, "max" = 0.2) ,
  cluster1_x_coord_range = c("min" = 125, "max" = 175)
)

SD1_metric_df_list <- analyse_SD1_simulations(parameters_df)
saveRDS(SD1_metric_df_list, "~/R/non_interactive_slurm/SD1_metric_df_list.RDS")
write.csv(SD1_metric_df_list$spe_list$sim_1, "~/R/non_interactive_slurm/SD1_dfs/spe1.csv")


SD1_metric_df_list <- readRDS("~/R/non_interactive_slurm/SD1_metric_df_list.RDS")



# Updated SPIAT function to match squidpy ----
calculate_all_single_radius_cc_metrics2D <- function(spe,
                                                     reference_cell_type,
                                                     target_cell_types,
                                                     radius,
                                                     feature_colname) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spe[[feature_colname]][trimws(spe[[feature_colname]]) == ""]) > 0) {
    stop("spe cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (!(is.character(reference_cell_type) && length(reference_cell_type) == 1)) {
    stop("`reference_cell_type` is not a character.")
  }
  if (!is.character(target_cell_types)) {
    stop("`target_cell_types` is not a character vector.")
  }
  if (!(is.numeric(radius) && length(radius) == 1 && radius > 0)) {
    stop(paste(radius, " is not a positive numeric."))
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  
  ## For reference_cell_type, check it is found in the spe object
  if (!(reference_cell_type %in% spe[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spe object"))
    return(NULL)
  }
  ## For target_cell_types, check they are found in the spe object
  unknown_cell_types <- setdiff(target_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
  }
  
  # Define result
  result <- list("neighbourhood_counts" = list(),
                 "cells_in_neighbourhood" = list(),
                 "neighbourhood_entropy" = list(),
                 "mixing_score" = list(),
                 "cross_K" = list(),
                 "cross_L" = list(),
                 "cross_G" = list(),
                 "fast_co_occurrence" = list())
  
  # Define other constants
  mixing_score_df_colnames <- c("ref_cell_type",
                                "tar_cell_type",
                                "n_ref_cells",
                                "n_tar_cells",
                                "n_ref_tar_interactions",
                                "n_ref_ref_interactions",
                                "mixing_score",
                                "normalised_mixing_score")
  
  cross_K_df_colnames <- c("reference",
                           "expected",
                           target_cell_types)
  cross_G_df_colnames <- c("observed_cross_G",
                           "expected_cross_G")
  fast_co_occurrence_df_colnames <- c("reference",
                                      target_cell_types)
  
  # Get rough dimensions of window for cross_K
  spe_coords <- data.frame(SpatialExperiment::spatialCoords(spe))
  length <- (max(spe_coords$Cell.X.Position) - min(spe_coords$Cell.X.Position))
  width  <- (max(spe_coords$Cell.Y.Position) - min(spe_coords$Cell.Y.Position))
  
  # Get area of the window the cells are in
  area <- length * width
  
  # All single radius cc metrics stem from calculate_neighbourhood_entropy2D function
  neighbourhood_entropy_df <- calculate_neighbourhood_entropy2D(spe,
                                                                reference_cell_type,
                                                                target_cell_types,
                                                                radius,
                                                                feature_colname)
  
  ## Neighbourhood counts ----------
  result[["neighbourhood_counts"]] <- neighbourhood_entropy_df[ , c("ref_cell_id", target_cell_types)]
  
  ## Cells in neighbourhood ----------
  result[["cells_in_neighbourhood"]] <- neighbourhood_entropy_df[ , c("ref_cell_id", paste(target_cell_types, "_prop", sep = ""))]
  
  ## Neighbourhood_entropy --------------
  result[["neighbourhood_entropy"]] <- neighbourhood_entropy_df[ , c("ref_cell_id", paste(target_cell_types, "_entropy", sep = ""))]
  
  ## Mixing score -----------------
  for (target_cell_type in target_cell_types) {
    mixing_score_df <- data.frame(matrix(nrow = 1, ncol = length(mixing_score_df_colnames)))
    colnames(mixing_score_df) <- mixing_score_df_colnames
    mixing_score_df$ref_cell_type <- reference_cell_type
    
    # No need to fill in mixing_score_df if the reference and target cell is the same
    if (reference_cell_type != target_cell_type) {
      mixing_score_df$tar_cell_type <- target_cell_type
      mixing_score_df$n_ref_cells <- sum(spe[[feature_colname]] == reference_cell_type)
      mixing_score_df$n_tar_cells <- sum(spe[[feature_colname]] == target_cell_type)
      mixing_score_df$n_ref_tar_interactions <- sum(neighbourhood_entropy_df[[target_cell_type]])
      mixing_score_df$n_ref_ref_interactions <- sum(neighbourhood_entropy_df[[reference_cell_type]])
      mixing_score_df$mixing_score <- mixing_score_df$n_ref_tar_interactions / (0.5 * mixing_score_df$n_ref_ref_interactions)
      mixing_score_df$normalised_mixing_score <- 0.5 * mixing_score_df$mixing_score * mixing_score_df$n_ref_cells / mixing_score_df$n_tar_cell
      if (is.infinite(mixing_score_df$mixing_score)) mixing_score_df$mixing_score <- NA
      if (is.infinite(mixing_score_df$normalised_mixing_score)) mixing_score_df$normalised_mixing_score <- NA
      result[["mixing_score"]][[target_cell_type]] <- mixing_score_df
    }
  }
  
  ## Cross_K ---------------------
  cross_K_df <- data.frame(matrix(nrow = 1, ncol = length(cross_K_df_colnames)))
  colnames(cross_K_df) <- cross_K_df_colnames
  cross_K_df$reference <- reference_cell_type
  cross_K_df$expected <- pi * radius^2
  
  for (target_cell_type in target_cell_types) {
    cross_K_df[[target_cell_type]] <- ((area * sum(neighbourhood_entropy_df[[target_cell_type]])) / length(spe[[feature_colname]]) ^ 2)
  }
  result[["cross_K"]] <- cross_K_df
  
  ## Cross_L ---------------------
  cross_L_df <- cross_K_df
  cross_L_df[ , c("expected", target_cell_types)] <- (cross_L_df[ , c("expected", target_cell_types)] / (pi)) ^ (1/2)
  result[["cross_L"]] <- cross_L_df
  
  ## Cross_G ---------------------
  for (target_cell_type in target_cell_types) {
    cross_G_df <- data.frame(matrix(nrow = 1, ncol = length(cross_G_df_colnames)))
    colnames(cross_G_df) <- cross_G_df_colnames
    
    reference_target_interactions <- neighbourhood_entropy_df[[target_cell_type]]
    n_target_cells <- sum(spe[[feature_colname]] == target_cell_type)
    target_cell_type_intensity <- n_target_cells / area
    observed_cross_G <- sum(reference_target_interactions != 0) / length(reference_target_interactions)
    expected_cross_G <- 1 - exp(-1 * target_cell_type_intensity * pi * radius^2)
    
    cross_G_df$observed_cross_G <- observed_cross_G
    cross_G_df$expected_cross_G <- expected_cross_G
    result[["cross_G"]][[target_cell_type]] <- cross_G_df
  }
  
  
  ## Fast_co_occurrence ---------------
  all_cell_types <- unique(spe[[feature_colname]])
  
  neighbourhood_counts_df <- calculate_neighbourhood_counts2D(spe,
                                                              reference_cell_type,
                                                              all_cell_types,
                                                              radius,
                                                              feature_colname)
  
  neighbourhood_counts_df$total <- rowSums(neighbourhood_counts_df[, -1], na.rm = TRUE)
  
  fast_co_occurrence_df <- data.frame(matrix(nrow = 1, ncol = length(fast_co_occurrence_df_colnames)))
  colnames(fast_co_occurrence_df) <- fast_co_occurrence_df_colnames
  fast_co_occurrence_df$reference <- reference_cell_type
  
  # Get total number of cells in radius around reference cell type
  n_cells_in_reference_cell_type_radius <- sum(neighbourhood_counts_df$total)
  
  # Get total number of cells in radius around all cell types
  n_cells_in_all_cell_type_radius <- dbscan::frNN(SpatialExperiment::spatialCoords(spe),
                                                  eps = radius,
                                                  query = SpatialExperiment::spatialCoords(spe),
                                                  sort = FALSE)
  n_cells_in_all_cell_type_radius <- sum(rapply(n_cells_in_all_cell_type_radius$id, base::length))
  
  for (target_cell_type in target_cell_types) {
    
    # Account for case where not enough target cells in the data
    if (sum(spe[[feature_colname]] == target_cell_type) <= 1) {
      fast_co_occurrence_df[[target_cell_type]] <- NA
      next
    }
    
    # Get total number of target cells in radius around reference cell type
    n_target_cells_in_reference_cell_type_radius <- sum(neighbourhood_counts_df[[target_cell_type]])
    
    # Get proportion of target cells in radius around reference cell type
    target_cell_type_proportion_in_reference_cell_type_radius <- n_target_cells_in_reference_cell_type_radius / n_cells_in_reference_cell_type_radius
    
    # Get total number of target cells in radius around all cell types
    n_target_cells_in_all_cell_type_radius <- dbscan::frNN(SpatialExperiment::spatialCoords(spe)[spe[[feature_colname]] == target_cell_type, ],
                                                           eps = radius,
                                                           query = SpatialExperiment::spatialCoords(spe),
                                                           sort = FALSE)
    n_target_cells_in_all_cell_type_radius <- sum(rapply(n_target_cells_in_all_cell_type_radius$id, base::length))
    
    # Get proportion of target cells in radius around all cell types
    target_cell_type_proportion_in_all_cell_type_radius <- n_target_cells_in_all_cell_type_radius / n_cells_in_all_cell_type_radius
    
    # Get co-occurence value for target cell type
    target_cell_type_co_occurrence <- target_cell_type_proportion_in_reference_cell_type_radius / target_cell_type_proportion_in_all_cell_type_radius
    
    # Add to result data frame
    fast_co_occurrence_df[[target_cell_type]] <- target_cell_type_co_occurrence
  }
  result[["fast_co_occurrence"]] <- fast_co_occurrence_df
  
  return(result)
}
# Plot spatial data----
spe2D <- read.csv("~/R/non_interactive_slurm/SD1_dfs/spe1.csv")

ggplot(spe2D,
       aes(Cell.X.Position, Cell.Y.Position, color = Cell.Type)) +
  geom_point(size = 4) +
  theme_bw() +
  labs(x = "x", y = "y", color = "Cell Type") +
  theme(
    # X and Y axis titles
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    
    # X and Y axis tick labels
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    
    # Legend title
    legend.title = element_text(size = 14),
    
    # Legend labels
    legend.text = element_text(size = 12)
  )

# Comparing SPIAT and squidpy -----
SD1_metric_df_list <- readRDS("~/R/non_interactive_slurm/SD1_metric_df_list.RDS")
spiat_co_df <- SD1_metric_df_list$metric_df_list2D$FCO
spiat_L_df <- SD1_metric_df_list$metric_df_list2D$CL
spiat_G_df <- SD1_metric_df_list$metric_df_list2D$CG


squidpy_L_stat <- read.csv("~/R/non_interactive_slurm/L_stat.csv")
colnames(squidpy_L_stat) <- c("radius", "reference", "value")
squidpy_L_stat$target <- squidpy_L_stat$reference

squidpy_G_stat <- read.csv("~/R/non_interactive_slurm/G_stat.csv")
colnames(squidpy_G_stat) <- c("radius", "reference", "value")
squidpy_G_stat$target <- squidpy_G_stat$reference

squidpy_co_stat <- read.csv("~/R/non_interactive_slurm/co_stat.csv")



compare_SPIAT_squidpy_co_occurrence <- function(spiat_df, squidpy_df, reference) {
  squidpy_df <- squidpy_df[squidpy_df$reference == reference, ]
  squidpy_df <- squidpy_df[squidpy_df$target != "O", ]
  colnames(squidpy_df)[colnames(squidpy_df) == "distance"] <- "radius"
  squidpy_df$tool <- "squidpy"
  
  spiat_df <- spiat_df[spiat_df$reference == reference, ]
  spiat_df$simulation <- NULL
  
  spiat_df <- reshape2::melt(spiat_df)
  spiat_df$variable <- as.numeric(gsub("r", "", spiat_df$variable))
  colnames(spiat_df)[colnames(spiat_df) == "variable"] <- "radius"
  spiat_df$tool <- "SPIAT"
  
  spiat_squidpy_co_df <- rbind(squidpy_df, spiat_df)
  spiat_squidpy_co_df$target_tool <- paste(spiat_squidpy_co_df$target, spiat_squidpy_co_df$tool, sep = "-")
  
  fig <- ggplot(spiat_squidpy_co_df, aes(radius, value, color = target_tool)) +
    geom_line() +
    theme_bw() +
    labs(y = "Co-occurrence", color = "") +
    scale_color_manual(
      labels = c(
        "A-SPIAT" = paste(reference, "/A using SPIAT", sep = ""),
        "A-squidpy" = paste(reference, "/A using squidpy", sep = ""),
        "B-SPIAT" = paste(reference, "/B using SPIAT", sep = ""),
        "B-squidpy" = paste(reference, "/B using squidpy", sep = "")
      ),
      values = c(
        "A-SPIAT" = "#48bbff",
        "A-squidpy" = "#bb0036",
        "B-SPIAT" = "#0062c5",
        "B-squidpy" = "#770026"
      )
    )
  
  methods::show(fig)
}

compare_SPIAT_squidpy_co_occurrence(spiat_co_df, squidpy_co_stat, "A")
compare_SPIAT_squidpy_co_occurrence(spiat_co_df, squidpy_co_stat, "B")




compare_SPIAT_squidpy_ripley_function <- function(spiat_df, squidpy_df) {
  squidpy_df <- squidpy_df[squidpy_df$target != "O", ]
  colnames(squidpy_df)[colnames(squidpy_df) == "distance"] <- "radius"
  squidpy_df$tool <- "squidpy"
  
  spiat_df$simulation <- NULL
  
  spiat_df <- spiat_df[spiat_df$reference == spiat_df$target, ]
  spiat_df <- reshape2::melt(spiat_df)
  spiat_df$variable <- as.numeric(gsub("r", "", spiat_df$variable))
  colnames(spiat_df)[colnames(spiat_df) == "variable"] <- "radius"
  spiat_df$tool <- "SPIAT"
  
  spiat_squidpy_combined_df <- rbind(squidpy_df, spiat_df)
  spiat_squidpy_combined_df$target_tool <- paste(spiat_squidpy_combined_df$target, spiat_squidpy_combined_df$tool, sep = "-")
  
  fig <- ggplot(spiat_squidpy_combined_df, aes(radius, value, color = target_tool)) +
    geom_line() +
    theme_bw() +
    labs(y = "Ripley's L function value", color = "") +
    scale_color_manual(
      labels = c(
        "A-SPIAT" = "A/A using SPIAT",
        "A-squidpy" = "A/A using squidpy",
        "B-SPIAT" = "B/B using SPIAT",
        "B-squidpy" = "B/B using squidpy"
      ),
      values = c(
        "A-SPIAT" = "#48bbff",
        "A-squidpy" = "#bb0036",
        "B-SPIAT" = "#0062c5",
        "B-squidpy" = "#770026"
      )
    )
  
  methods::show(fig)
}

compare_SPIAT_squidpy_ripley_function(spiat_L_df, squidpy_L_stat)
compare_SPIAT_squidpy_ripley_function(spiat_G_df, squidpy_G_stat)
