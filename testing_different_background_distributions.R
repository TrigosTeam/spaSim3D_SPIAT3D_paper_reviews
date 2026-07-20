# Set seed for consistency
set.seed(678999821)
library(SpatialExperiment)
library(dplyr)

# Functions -----
check_input_parameters <- function(input_parameters) {
  
  input_parameter_names <- names(input_parameters)
  
  check_value <- 0
  
  is_equal_with_tolerance <- function(x, y, tolerance = 1e-6) {
    abs(x - y) <= tolerance
  }
  
  for (input_parameter_name in input_parameter_names) {
    input_parameter <- input_parameters[[input_parameter_name]]
    
    # spe
    if (input_parameter_name == "spe" && class(input_parameter) != "SpatialExperiment") {
      check_value <- 1
      break
    }
    # Positive integer
    if (input_parameter_name %in% c("n_cells", "n_edges") && !(is.integer(input_parameter) && length(input_parameter) == 1 || (is.numeric(input_parameter) && length(input_parameter) == 1 && input_parameter > 0 && input_parameter%%1 == 0))) {
      check_value <- 2
      break
    }  
    # Positive numeric
    if (input_parameter_name %in% c("length", "width", "height", "radius", "x_radius", "y_radius", "z_radius", "ring_width", "inner_ring_width", "outer_ring_width") && !(is.numeric(input_parameter) && length(input_parameter) == 1 && input_parameter > 0)) {
      check_value <- 3
      break
    }
    # Non-negative numeric
    if (input_parameter_name %in% c("minimum_distance_between_cells") && !(is.numeric(input_parameter) && length(input_parameter) == 1 && input_parameter >= 0)) {
      check_value <- 4
      break
    }
    # Numeric between 0 and 1
    if (input_parameter_name %in% c("jitter_proportion") && !(is.numeric(input_parameter) && length(input_parameter) == 1 && input_parameter >= 0 && input_parameter <= 1)) {
      check_value <- 5
      break
    }
    # Character
    if (input_parameter_name %in% c("background_cell_type") && !(is.character(input_parameter)) && length(input_parameter) == 1) {
      check_value <- 6
      break
    }
    # Logical
    if (input_parameter_name %in% c("plot_image") && !(is.logical(input_parameter)) && length(input_parameter) == 1) {
      check_value <- 7
      break
    }
    # Character vector
    if (input_parameter_name %in% c("cell_types", "cluster_cell_types", "ring_cell_types", "inner_ring_cell_types", "outer_ring_cell_types") && 
        !(is.character(input_parameter))) {
      check_value <- 8
      break
    }
    # Numeric vector
    if (input_parameter_name %in% c("cell_proportions", "cluster_cell_proportions", "ring_cell_proportions", "inner_ring_cell_proportions", "outer_ring_cell_proportions") && 
        !(is.numeric(input_parameter))) {
      check_value <- 9
      break
    }
    # Numeric vector contains values between 0 and 1
    if (input_parameter_name %in% c("cell_proportions", "cluster_cell_proportions", "ring_cell_proportions", "inner_ring_cell_proportions", "outer_ring_cell_proportions") && 
        sum(input_parameter < 0 | input_parameter > 1) != 0) {
      check_value <- 10
      break
    }
    # Numeric vector contains values that sum to 1
    if (input_parameter_name %in% c("cell_proportions", "cluster_cell_proportions", "ring_cell_proportions", "inner_ring_cell_proportions", "outer_ring_cell_proportions") && 
        !is_equal_with_tolerance(sum(input_parameter), 1)) {
      check_value <- 11
      break
    }
    # Numeric vector of length 3
    if (input_parameter_name %in% c("centre_loc", "start_loc", "end_loc") && !(is.numeric(input_parameter) && length(input_parameter) == 3)) {
      check_value <- 12
      break
    }
    # Numeric
    if (input_parameter_name %in% c("y_z_rotation", "x_z_rotation", "x_y_rotation") && !(is.numeric(input_parameter) && length(input_parameter) == 1)) {
      check_value <- 13
      break
    }
  }
  
  # Two vectors match in length
  if (check_value == 0) {
    pairs <- data.frame(name1 = c("cell_types", "cluster_cell_types", "ring_cell_types", "inner_ring_cell_types", "outer_ring_cell_types"),
                        name2 = c("cell_proportions", "cluster_cell_proportions", "ring_cell_proportions", "inner_ring_cell_proportions", "outer_ring_cell_proportions"))  
    
    for (i in seq(nrow(pairs))) {
      name1 <- pairs[["name1"]][i]
      name2 <- pairs[["name2"]][i]
      if (name1 %in% input_parameter_names && name2 %in% input_parameter_names && length(input_parameters[[name1]]) != length(input_parameters[[name2]])) {
        check_value <- 14
        input_parameter_name <- c(name1, name2)
        break
      }
    }
  }
  
  # If check_value equals 0, all inputs are valid.
  if (check_value == 0) {
    return(TRUE)
  }
  # At least one input is not valid, return the first invalid input. 
  else {
    return(list(input_parameter_name = input_parameter_name, check_value = check_value)) 
  }
}

input_parameter_error_message <- function(input_parameter_check_value) {
  
  input_parameter_name <- input_parameter_check_value[[1]]
  check_value <- input_parameter_check_value[[2]]
  
  error_message <- switch(check_value,
                          "1" = paste(input_parameter_name, "is not a SpatialExperiment object."),
                          "2" = paste(input_parameter_name, "is not a positive integer."),
                          "3" = paste(input_parameter_name, "is not a positive numeric."),
                          "4" = paste(input_parameter_name, "is not a non-negative numeric."),
                          "5" = paste(input_parameter_name, "is not a numeric between 0 and 1."),
                          "6" = paste(input_parameter_name, "is not a character."),
                          "7" = paste(input_parameter_name, "is not a logical (TRUE or FALSE)."),
                          "8" = paste(input_parameter_name, "is not a character vector."),
                          "9" = paste(input_parameter_name, "is not a numeric vector."),
                          "10" = paste(input_parameter_name, "cannot be negative or greater than 1."),
                          "11" = paste(input_parameter_name, "does not sum to 1."),
                          "12" = paste(input_parameter_name, "is not a numeric vector of length 3."),
                          "13" = paste(input_parameter_name, "is not a numeric."),
                          "14" = paste(input_parameter_name[1], "and", input_parameter_name[2], "do not match in length."))
  
  return(error_message)
}

# Generate data frame containing values of important parameters
generate_random_parameters <- function(
    n_simulations_per_parameter,
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
  
  parameters <- c("bg_prop_A", 
                  "bg_prop_B", 
                  "E_radius_x", # E_radius_x, E_radius_y, E_radius_z count as 1
                  "N_width", 
                  "cluster_prop_A", 
                  "ring_width_factor", 
                  "cluster1_x_coord")
  n_parameters <- length(parameters) 
  n_simulations <- n_simulations_per_parameter * n_parameters
  
  # Set parameter values to be the average of the min and max set in the parameter range, bg_props are an exception and will be equal to 0
  parameters_df <- data.frame(
    bg_prop_A = 0,
    bg_prop_B = 0,
    E_radius_x = mean(c(E_radius_x_range["min"], E_radius_x_range["max"])),
    E_radius_y = mean(c(E_radius_y_range["min"], E_radius_y_range["max"])),
    E_radius_z = mean(c(E_radius_z_range["min"], E_radius_z_range["max"])),
    N_width = mean(c(N_width_range["min"], N_width_range["max"])),
    cluster_prop_A = mean(c(cluster_prop_A_range["min"], cluster_prop_A_range["max"])),
    ring_width_factor = mean(c(ring_width_factor_range["min"], ring_width_factor_range["max"])),
    cluster1_x_coord = mean(c(cluster1_x_coord_range["min"], cluster1_x_coord_range["max"])),
    variable_parameter = rep(parameters, each = n_simulations_per_parameter)
  )
  
  # Put all parameter ranges into a list
  parameter_ranges <- list(
    bg_prop_A = bg_prop_A_range,
    bg_prop_B = bg_prop_B_range,
    E_radius_x = E_radius_x_range,
    E_radius_y = E_radius_y_range,
    E_radius_z = E_radius_z_range,
    N_width = N_width_range,
    cluster_prop_A = cluster_prop_A_range,
    ring_width_factor = ring_width_factor_range,
    cluster1_x_coord = cluster1_x_coord_range
  )
  
  # For the current variable parameter, vary the parameter...
  for (i in seq_len(nrow(parameters_df))) {
    
    curr_variable_parameter <- parameters_df[["variable_parameter"]][i]
    
    if (curr_variable_parameter == "E_radius_x") {
      parameters_df[i, "E_radius_x"] <- runif(1, E_radius_x_range[["min"]], E_radius_x_range[["max"]])
      parameters_df[i, "E_radius_y"] <- runif(1, E_radius_y_range[["min"]], E_radius_y_range[["max"]])
      parameters_df[i, "E_radius_z"] <- runif(1, E_radius_z_range[["min"]], E_radius_z_range[["max"]])
    }
    else {
      parameter_min <- parameter_ranges[[curr_variable_parameter]][["min"]]
      parameter_max  <- parameter_ranges[[curr_variable_parameter]][["max"]]
      parameters_df[i, curr_variable_parameter] <- runif(1, parameter_min, parameter_max)
    }
  }
  
  return(parameters_df)
}

# Generate updated parameters data frame to be specific for arrangement and shape
generate_updated_parameters_for_arrangement_and_shape <- function(
    n_simulations_per_parameter,
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
  
  arrangements <- c("mixed", "ringed", "separated")
  shapes <- c("ellipsoid", "network")
  
  updated_parameter_df <- data.frame()
  
  for (arrangement in arrangements) {
    for (shape in shapes) {
      parameters_df <- generate_random_parameters(
        n_simulations_per_parameter = n_simulations_per_parameter,
        bg_prop_A_range = bg_prop_A_range,
        bg_prop_B_range = bg_prop_B_range,
        E_radius_x_range = E_radius_x_range,
        E_radius_y_range = E_radius_y_range,
        E_radius_z_range = E_radius_z_range,
        N_width_range = N_width_range,
        cluster_prop_A_range = cluster_prop_A_range,
        ring_width_factor_range = ring_width_factor_range,
        cluster1_x_coord_range = cluster1_x_coord_range
      )   
      
      parameters_df$arrangement <- arrangement
      parameters_df$shape <- shape
      
      if (arrangement == "mixed") {
        parameters_df <- parameters_df[!(parameters_df[["variable_parameter"]] %in% c("ring_width_factor", "cluster1_x_coord")), ] 
      }
      else if (arrangement == "ringed") {
        parameters_df <- parameters_df[!(parameters_df[["variable_parameter"]] %in% c("cluster_prop_A", "cluster1_x_coord")), ] 
      }
      else if (arrangement == "separated") {
        parameters_df <- parameters_df[!(parameters_df[["variable_parameter"]] %in% c("cluster_prop_A", "ring_width_factor")), ] 
      }
      
      if (shape == "ellipsoid") {
        parameters_df <- parameters_df[!(parameters_df[["variable_parameter"]] %in% c("N_width")), ] 
      }
      else if (shape == "network") {
        parameters_df <- parameters_df[!(parameters_df[["variable_parameter"]] %in% c("E_radius_x")), ] 
      }
      updated_parameter_df <- rbind(updated_parameter_df, parameters_df)
    }
  }
  return(updated_parameter_df)
}


# Analyse SD2 simulations
analyse_SD2_simulations <- function(parameters_df, background_type) {
  
  # Set defined parameters/values
  cell_types <- c('A', 'B')
  n_simulations <- nrow(parameters_df)
  
  radii <- seq(20, 100, 10)
  radii_colnames <- paste("r", radii, sep = "")
  
  n_splits <- 10
  thresholds <- seq(0.01, 1, 0.01)
  thresholds_colnames <- paste("t", thresholds, sep = "")
  
  # z-coords for slices
  bottom_z_coord_of_slices <- c(145, 175, 205)
  top_z_coord_of_slices <- bottom_z_coord_of_slices + 10
  n_slices <- length(bottom_z_coord_of_slices)
  
  # Generate simulation metadata using parameters data frame from 'generate_random_parameters' function
  generate_simulation_metadata <- function(parameters_df, simulation_index) {
    
    # Constant metadata for simulations
    simulation_metadata <- spe_metadata_background_template(background_type)
    simulation_metadata$background$n_cells <- 30000
    simulation_metadata$background$length <- 600
    simulation_metadata$background$width <- 600
    simulation_metadata$background$height <- 300
    simulation_metadata$background$cell_types <- c("A", "B", "O") # Cell proportions will change later
    
    if (background_type == "hardcore") simulation_metadata$background$minimum_distance_between_cells <- 10
    if (background_type == "ordered") simulation_metadata$background$jitter_proportion <- 0
    
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
  
  # Function to get all slices from spe
  get_all_spe_slices <- function(spe, 
                                 bottom_z_coord_of_slices, 
                                 top_z_coord_of_slices) {
    
    spe_slices <- list()
    
    if (length(bottom_z_coord_of_slices) != length(top_z_coord_of_slices)) stop("Lengths of bottom_z_coords_of_slices and top_z_coords_of_slices should be equal.")
    number_of_slices <- length(bottom_z_coord_of_slices)
    
    for (i in seq(number_of_slices)) {
      bottom_z_coord_of_slice <- bottom_z_coord_of_slices[i]
      top_z_coord_of_slice <- top_z_coord_of_slices[i]
      z_coords_of_cells_in_spe <- SpatialExperiment::spatialCoords(spe)[ , "Cell.Z.Position"]
      spe_for_slice <- spe[, bottom_z_coord_of_slice < z_coords_of_cells_in_spe & z_coords_of_cells_in_spe < top_z_coord_of_slice]
      SpatialExperiment::spatialCoords(spe_for_slice) <- SpatialExperiment::spatialCoords(spe_for_slice)[ , c("Cell.X.Position", "Cell.Y.Position")]
      
      spe_slices[[i]] <- spe_for_slice
    }
    
    return(spe_slices)
  }
  
  # Function to create empty metric df list
  create_empty_metric_df_list3D <- function(
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
    
    # Define MS, NMS, ANC, ACIN, COO, ANE data frames as well as constants
    radii_colnames <- paste("r", radii, sep = "")
    
    MS_df_colnames <- c("simulation", "reference", "target", radii_colnames)
    MS_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(MS_df_colnames)))
    colnames(MS_df) <- MS_df_colnames
    
    NMS_df <- ANC_df <- ANE_df <- ACIN_df <- COO_df <- CK_df <- CL_df <- CG_df <- MS_df
    
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
                           COO = COO_df,
                           CK = CK_df,
                           CL = CL_df,
                           CG = CG_df,
                           PBSAC = PBSAC_df,
                           PBP = PBP_df,
                           EBSAC = EBSAC_df,
                           EBP = EBP_df)
    
    return(metric_df_list)
  }
  
  create_empty_metric_df_list2D <- function(
    cell_types,
    n_simulations,
    n_slices,
    radii_colnames,
    thresholds_colnames
  ) {
    n_cell_type_combinations <- length(cell_types)^2
    
    # Define AMD data frames as well as constants
    AMD_df_colnames <- c("simulation", "slice", "reference", "target", "AMD")
    AMD_df <- data.frame(matrix(nrow = n_simulations * n_slices * n_cell_type_combinations, ncol = length(AMD_df_colnames)))
    colnames(AMD_df) <- AMD_df_colnames
    
    # Define MS, NMS, ANC, ACIN, COO, ANE data frames as well as constants
    radii_colnames <- paste("r", radii, sep = "")
    
    MS_df_colnames <- c("simulation", "slice", "reference", "target", radii_colnames)
    MS_df <- data.frame(matrix(nrow = n_simulations * n_slices * n_cell_type_combinations, ncol = length(MS_df_colnames)))
    colnames(MS_df) <- MS_df_colnames
    
    NMS_df <- ANC_df <- ANE_df <- ACIN_df <- COO_df <- CK_df <- CL_df <- CG_df <- MS_df
    
    # Define SAC and prevalence data frames as well as constants
    thresholds_colnames <- paste("t", thresholds, sep = "")
    
    PBSAC_df_colnames <- c("simulation", "slice", "reference", "target", "PBSAC")
    PBSAC_df <- data.frame(matrix(nrow = n_simulations * n_slices * n_cell_type_combinations, ncol = length(PBSAC_df_colnames)))
    colnames(PBSAC_df) <- PBSAC_df_colnames
    
    PBP_df_colnames <- c("simulation", "slice", "reference", "target", thresholds_colnames)
    PBP_df <- data.frame(matrix(nrow = n_simulations * n_slices * n_cell_type_combinations, ncol = length(PBP_df_colnames)))
    colnames(PBP_df) <- PBP_df_colnames
    
    EBSAC_df_colnames <- c("simulation", "slice", "cell_types", "EBSAC")
    EBSAC_df <- data.frame(matrix(nrow = n_simulations * n_slices * n_cell_type_combinations, ncol = length(EBSAC_df_colnames)))
    colnames(EBSAC_df) <- EBSAC_df_colnames
    
    EBP_df_colnames <- c("simulation", "slice", "cell_types", thresholds_colnames)
    EBP_df <- data.frame(matrix(nrow = n_simulations * n_slices * n_cell_type_combinations, ncol = length(EBP_df_colnames)))
    colnames(EBP_df) <- EBP_df_colnames
    
    
    # Add all to list:
    metric_df_list <- list(AMD = AMD_df,
                           MS = MS_df,
                           NMS = NMS_df,
                           ACIN = ACIN_df,
                           ANE = ANE_df,
                           ANC = ANC_df,
                           CK = CK_df,
                           CL = CL_df,
                           CG = CG_df,
                           COO = COO_df,
                           PBSAC = PBSAC_df,
                           PBP = PBP_df,
                           EBSAC = EBSAC_df,
                           EBP = EBP_df)
    
    return(metric_df_list)
  }
  
  
  # Function analyse spe in 3D
  analyse_simulation3D <- function(spe, 
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
    
    minimum_distance_data <- calculate_minimum_distances_between_cell_types3D(spe,
                                                                              cell_types,
                                                                              feature_colname = "Cell.Type",
                                                                              show_summary = F,
                                                                              plot_image = F)
    
    minimum_distance_data_summary <- summarise_distances_between_cell_types3D(minimum_distance_data)
    
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "simulation"] <- simulation_index
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "reference"] <- minimum_distance_data_summary$reference
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "target"] <- minimum_distance_data_summary$target
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "AMD"] <- minimum_distance_data_summary$mean
    
    # Need a new index for 'pair' data which increments after each reference-target pair
    pair_index <- n_cell_type_combinations * (simulation_index - 1) + 1 
    
    for (reference_cell_type in cell_types) {
      gradient_data <- calculate_all_gradient_cc_metrics3D(spe,
                                                           reference_cell_type,
                                                           cell_types,
                                                           radii,
                                                           feature_colname = "Cell.Type",
                                                           plot_image = F)
      
      for (target_cell_type in cell_types) {
        print(paste(reference_cell_type, target_cell_type, sep = "/"))
        metric_df_list[["ANC"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["ACIN"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["COO"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["CK"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["CL"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["CG"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["MS"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["NMS"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["ANE"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, 
                                                                                         paste(reference_cell_type, target_cell_type, sep = ","))
        
        
        if (is.null(gradient_data)) {
          metric_df_list[["ANC"]][pair_index, radii_colnames] <- NA
          metric_df_list[["COO"]][pair_index, radii_colnames] <- NA
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
          metric_df_list[["COO"]][pair_index, radii_colnames] <- gradient_data[["co_occurrence"]][[target_cell_type]]
          metric_df_list[["CK"]][pair_index, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]] - gradient_data[["cross_K"]][["expected"]]
          metric_df_list[["CL"]][pair_index, radii_colnames] <- gradient_data[["cross_L"]][[target_cell_type]] - gradient_data[["cross_L"]][["expected"]]
          metric_df_list[["CG"]][pair_index, radii_colnames] <- gradient_data[["cross_G"]][[target_cell_type]][["observed_cross_G"]] - gradient_data[["cross_G"]][[target_cell_type]][["expected_cross_G"]]
          
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
          proportion_grid_metrics <- calculate_cell_proportion_grid_metrics3D(spe, 
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
            PBSAC <- calculate_spatial_autocorrelation3D(proportion_grid_metrics, 
                                                         "proportion",
                                                         weight_method = "queen")
            
            PBP_df <- calculate_prevalence_gradient3D(proportion_grid_metrics,
                                                      "proportion",
                                                      show_AUC = F,
                                                      plot_image = F)
            
            
            metric_df_list[["PBSAC"]][pair_index, "PBSAC"] <- PBSAC
            metric_df_list[["PBP"]][pair_index, thresholds_colnames] <- PBP_df$prevalence
          } 
          
          entropy_grid_metrics <- calculate_entropy_grid_metrics3D(spe, 
                                                                   n_splits,
                                                                   c(reference_cell_type, target_cell_type), 
                                                                   feature_colname = "Cell.Type",
                                                                   plot_image = F)
          
          if (is.null(entropy_grid_metrics)) {
            metric_df_list[["EBSAC"]][pair_index, "EBSAC"] <- NA
            metric_df_list[["EBP"]][pair_index, thresholds_colnames] <- NA
          }
          else {
            EBSAC <- calculate_spatial_autocorrelation3D(entropy_grid_metrics, 
                                                         "entropy",
                                                         weight_method = "queen")
            
            EBP_df <- calculate_prevalence_gradient3D(entropy_grid_metrics,
                                                      "entropy",
                                                      show_AUC = F,
                                                      plot_image = F)
            
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
  
  # Function analyse spe in 2D
  analyse_simulation2D <- function(spe_slices, 
                                   cell_types, 
                                   radii, 
                                   thresholds, 
                                   n_splits, 
                                   simulation_index, 
                                   metric_df_list) {
    n_cell_type_combinations <- length(cell_types)^2
    n_slices <- length(spe_slices)
    radii_colnames <- paste("r", radii, sep = "")
    thresholds_colnames <- paste("t", thresholds, sep = "")
    
    for (slice_index in seq_len(n_slices)) {
      spe <- spe_slices[[slice_index]]
      
      index <- n_cell_type_combinations * (simulation_index - 1) * n_slices +  n_cell_type_combinations * (slice_index - 1) + 1 
      
      minimum_distance_data <- calculate_minimum_distances_between_cell_types2D(spe,
                                                                                cell_types,
                                                                                feature_colname = "Cell.Type",
                                                                                show_summary = F,
                                                                                plot_image = F)
      
      minimum_distance_data_summary <- summarise_distances_between_cell_types2D(minimum_distance_data)
      
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "simulation"] <- simulation_index
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "slice"] <- slice_index
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "reference"] <- minimum_distance_data_summary$reference
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "target"] <- minimum_distance_data_summary$target
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "AMD"] <- minimum_distance_data_summary$mean
      
      # Need a new index for 'pair' data which increments after each reference-target pair
      pair_index <- n_cell_type_combinations * (simulation_index - 1) * n_slices +  n_cell_type_combinations * (slice_index - 1) + 1 
      
      for (reference_cell_type in cell_types) {
        gradient_data <- calculate_all_gradient_cc_metrics2D(spe,
                                                             reference_cell_type,
                                                             cell_types,
                                                             radii,
                                                             feature_colname = "Cell.Type",
                                                             plot_image = F)
        
        for (target_cell_type in cell_types) {
          print(paste(reference_cell_type, target_cell_type, sep = "/"))
          metric_df_list[["ANC"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["ACIN"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["COO"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["CK"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["CL"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["CG"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["MS"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["NMS"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["ANE"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, 
                                                                                                    paste(reference_cell_type, target_cell_type, sep = ","))
          
          
          if (is.null(gradient_data)) {
            metric_df_list[["ANC"]][pair_index, radii_colnames] <- NA
            metric_df_list[["COO"]][pair_index, radii_colnames] <- NA
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
            metric_df_list[["COO"]][pair_index, radii_colnames] <- gradient_data[["co_occurrence"]][[target_cell_type]]
            metric_df_list[["CK"]][pair_index, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]] - gradient_data[["cross_K"]][["expected"]]
            metric_df_list[["CL"]][pair_index, radii_colnames] <- gradient_data[["cross_L"]][[target_cell_type]] - gradient_data[["cross_L"]][["expected"]]
            metric_df_list[["CG"]][pair_index, radii_colnames] <- gradient_data[["cross_G"]][[target_cell_type]][["observed_cross_G"]] - gradient_data[["cross_G"]][[target_cell_type]][["expected_cross_G"]]
            
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
          metric_df_list[["PBSAC"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["PBP"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          
          metric_df_list[["EBSAC"]][pair_index, c("simulation", "slice", "cell_types")] <- c(simulation_index, slice_index, paste(reference_cell_type, target_cell_type, sep = ","))
          metric_df_list[["EBP"]][pair_index, c("simulation", "slice", "cell_types")] <- c(simulation_index, slice_index, paste(reference_cell_type, target_cell_type, sep = ","))
          
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
                                                        show_AUC = F,
                                                        plot_image = F)
              
              
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
                                                        show_AUC = F,
                                                        plot_image = F)
              
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
    }
    return(metric_df_list)
  }
  
  # Define metric df lists
  metric_df_list3D <- create_empty_metric_df_list3D(cell_types, n_simulations, radii_colnames, thresholds_colnames)
  metric_df_list2D <- create_empty_metric_df_list2D(cell_types, n_simulations, n_slices, radii_colnames, thresholds_colnames)
  
  for (simulation_index in seq_len(n_simulations)) {
    print(simulation_index)
    # Simulate spes in 3D, and get 2D slices
    simulation_metadata <- generate_simulation_metadata(parameters_df, simulation_index)
    
    spe3D <- simulate_spe_metadata3D(simulation_metadata, plot_image = F)
    
    spe_slices2D <- get_all_spe_slices(spe3D, bottom_z_coord_of_slices, top_z_coord_of_slices)
    
    # Analyse spes
    metric_df_list3D <- analyse_simulation3D(spe3D,
                                             cell_types,
                                             radii,
                                             thresholds,
                                             n_splits,
                                             simulation_index,
                                             metric_df_list3D)
    
    metric_df_list2D <- analyse_simulation2D(spe_slices2D,
                                             cell_types,
                                             radii,
                                             thresholds,
                                             n_splits,
                                             simulation_index,
                                             metric_df_list2D)
  }
  
  ## Turn gradient metrics into AUC and add to metric_df list
  get_AUC_for_radii_gradient_metrics <- function(y) {
    x <- radii
    h <- diff(x)[1]
    n <- length(x)
    
    AUC <- (h / 2) * (y[1] + 2 * sum(y[2:(n - 1)]) + y[n])
    
    return(AUC)
  }
  
  add_AUC_for_radii_gradient_metrics_to_metric_df_list3D <- function(metric_df_list) {
    gradient_radii_metrics <- c("MS", "NMS", "ACIN", "ANE", "ANC", "COO", "CK", "CL", "CG")
    
    for (metric in gradient_radii_metrics) {
      metric_AUC_name <- paste(metric, "AUC", sep = "_")
      
      if (metric %in% c("MS", "NMS", "ANC", "ACIN", "ANE", "COO", "CK", "CL", "CG")) {
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
  add_AUC_for_radii_gradient_metrics_to_metric_df_list2D <- function(metric_df_list) {
    gradient_radii_metrics <- c("MS", "NMS", "ACIN", "ANE", "ANC", "COO", "CK", "CL", "CG")
    
    for (metric in gradient_radii_metrics) {
      metric_AUC_name <- paste(metric, "AUC", sep = "_")
      
      if (metric %in% c("MS", "NMS", "ANC", "ACIN", "ANE", "COO", "CK", "CL", "CG")) {
        subset_colnames <- c("simulation", "slice", "reference", "target", metric_AUC_name)
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
    PBP_AUC_df <- PBP_df[ , c("simulation", "slice", "reference", "target", "PBP_AUC")]
    metric_df_list[["PBP_AUC"]] <- PBP_AUC_df
    
    # EBP_AUC 3D
    EBP_df <- metric_df_list[["EBP"]]
    EBP_df$EBP_AUC <- apply(EBP_df[ , thresholds_colnames], 1, function(y) {
      sum(diff(thresholds) * (head(y, -1) + tail(y, -1)) / 2)
    })
    EBP_AUC_df <- EBP_df[ , c("simulation", "slice", "cell_types", "EBP_AUC")]
    metric_df_list[["EBP_AUC"]] <- EBP_AUC_df
    
    return(metric_df_list)
  }
  
  metric_df_list3D <- add_AUC_for_radii_gradient_metrics_to_metric_df_list3D(metric_df_list3D)
  metric_df_list2D <- add_AUC_for_radii_gradient_metrics_to_metric_df_list2D(metric_df_list2D)
  
  # Combine 3D and 2D metric df lists
  metric_df_list_combined <- list()
  
  metrics <- names(metric_df_list3D)
  for (metric in metrics) {
    metric_df_list3D[[metric]][["slice"]] <- 0 # Set slice index for 3D to be 0
    
    metric_df_combined <- rbind(metric_df_list3D[[metric]], metric_df_list2D[[metric]])
    
    metric_df_list_combined[[metric]] <- metric_df_combined
  }
  
  return(metric_df_list_combined)
}


# Updated functions to generate different background distributions -----
simulate_spe_metadata3D <- function(spe_metadata,
                                    plot_image = TRUE,
                                    plot_cell_types = NULL,
                                    plot_colours = NULL) {
  
  # First element should contain background metadata
  bg_metadata <- spe_metadata[[1]]
  if (!(is.character(bg_metadata$background_type) && length(bg_metadata$background_type) == 1)) {
    stop("background_type parameter found in the metadata background list is not a character.")
  }
  
  if (bg_metadata$background_type == "uniform") {
    spe <- simulate_uniform_background_cells3D(bg_metadata$n_cells,
                                               bg_metadata$length,
                                               bg_metadata$width,
                                               bg_metadata$height,
                                               plot_image = F)
  }
  else if (bg_metadata$background_type == "poisson") {
    spe <- simulate_poisson_background_cells3D(bg_metadata$n_cells,
                                               bg_metadata$length,
                                               bg_metadata$width,
                                               bg_metadata$height,
                                               plot_image = F)
  }
  else if (bg_metadata$background_type == "hardcore") {
    spe <- simulate_hardcore_background_cells3D(bg_metadata$n_cells,
                                                bg_metadata$length,
                                                bg_metadata$width,
                                                bg_metadata$height,
                                                bg_metadata$minimum_distance_between_cells,
                                                plot_image = F)
  }
  else if (bg_metadata$background_type == "ordered") {
    spe <- simulate_ordered_background_cells3D(bg_metadata$n_cells,
                                               bg_metadata$length,
                                               bg_metadata$width,
                                               bg_metadata$height,
                                               bg_metadata$jitter_proportion,
                                               plot_image = F)
  }
  else {
    stop("background_type parameter is invalid.")
  }
  # Apply background mixing
  spe <- simulate_mixing3D(spe,
                           bg_metadata$cell_types,
                           bg_metadata$cell_proportions,
                           plot_image = F)
  
  ### If there is only background metadata, we are done
  if (length(spe_metadata) == 1) {
    
    # Plot
    if (plot_image) {
      fig <- plot_cells3D(spe,
                          plot_cell_types,
                          plot_colours)
      methods::show(fig)
    }
    
    return(spe)
  }
  
  ### All other elements should help to simulate clusters
  for (i in 2:length(spe_metadata)) {
    cluster_metadata <- spe_metadata[[i]]
    
    if (!(is.character(cluster_metadata$cluster_type) && length(cluster_metadata$cluster_type) == 1)) {
      stop(paste("cluster_type parameter found in the metadata cluster list", i,"is not a character."))
    }
    
    if (cluster_metadata$cluster_type == "regular") {
      spe <- simulate_clusters3D(spe, list(cluster_metadata), plot_image = F)
    }
    else if (cluster_metadata$cluster_type == "ring") {
      spe <- simulate_rings3D(spe, list(cluster_metadata), plot_image = F)
    }
    else if (cluster_metadata$cluster_type == "double ring") {
      spe <- simulate_double_rings3D(spe, list(cluster_metadata), plot_image = F)
    }
    else {
      stop("cluster_type parameter must be either 'regular', 'ring' or 'double ring'.")
    }
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe,
                        plot_cell_types,
                        plot_colours)
    methods::show(fig)
  }
  
  return(spe)
}


simulate_uniform_background_cells3D <- function(n_cells,
                                                length,
                                                width,
                                                height,
                                                background_cell_type = "Others",
                                                plot_image = TRUE) {
  
  # Check input parameters
  input_parameters <- list("n_cells" = n_cells,
                           "length" = length,
                           "width" = width,
                           "height" = height,
                           "background_cell_type" = background_cell_type,
                           "plot_image" = plot_image)
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  spatial_df <- data.frame("Cell.X.Position" = runif(n_cells, 0, length),
                           "Cell.Y.Position" = runif(n_cells, 0, width),
                           "Cell.Z.Position" = runif(n_cells, 0, height))
  
  # Add Cell.Type and Cell.ID
  spatial_df$Cell.Type <- background_cell_type
  spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")
  
  # Get metadata
  background_metadata <- list("background_type" = "uniform",
                              "n_cells" = n_cells,
                              "length" = length,
                              "width" = width,
                              "height" = height,
                              "cell_types" = background_cell_type,
                              "cell_proportions" = 1)
  simulation_metadata <- list(background = background_metadata)
  
  ## Convert data frame to spe object
  spe <- SpatialExperiment::SpatialExperiment(
    assay = matrix(data = NA, nrow = 0, ncol = nrow(spatial_df)),
    colData = spatial_df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = list(simulation = simulation_metadata))
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe,
                        background_cell_type,
                        "lightgray")
    methods::show(fig)
  }
  
  return(spe)
}

simulate_poisson_background_cells3D <- function(n_cells,
                                                length,
                                                width,
                                                height,
                                                background_cell_type = "Others",
                                                plot_image = TRUE) {
  
  # Update n_cells using poisson distribution
  n_cells <- rpois(1, n_cells)
  
  # Check input parameters
  input_parameters <- list("n_cells" = n_cells,
                           "length" = length,
                           "width" = width,
                           "height" = height,
                           "background_cell_type" = background_cell_type,
                           "plot_image" = plot_image)
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  spatial_df <- data.frame("Cell.X.Position" = runif(n_cells, 0, length),
                           "Cell.Y.Position" = runif(n_cells, 0, width),
                           "Cell.Z.Position" = runif(n_cells, 0, height))
  
  # Add Cell.Type and Cell.ID
  spatial_df$Cell.Type <- background_cell_type
  spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")
  
  # Get metadata
  background_metadata <- list("background_type" = "poisson",
                              "n_cells" = n_cells,
                              "length" = length,
                              "width" = width,
                              "height" = height,
                              "cell_types" = background_cell_type,
                              "cell_proportions" = 1)
  simulation_metadata <- list(background = background_metadata)
  
  ## Convert data frame to spe object
  spe <- SpatialExperiment::SpatialExperiment(
    assay = matrix(data = NA, nrow = 0, ncol = nrow(spatial_df)),
    colData = spatial_df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = list(simulation = simulation_metadata))
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe,
                        background_cell_type,
                        "lightgray")
    methods::show(fig)
  }
  
  return(spe)
}

simulate_hardcore_background_cells3D <- function(n_cells,
                                                 length,
                                                 width,
                                                 height,
                                                 minimum_distance_between_cells,
                                                 background_cell_type = "Others",
                                                 plot_image = TRUE) {
  
  # Check input parameters
  input_parameters <- list("n_cells" = n_cells,
                           "length" = length,
                           "width" = width,
                           "height" = height,
                           "minimum_distance_between_cells" = minimum_distance_between_cells,
                           "background_cell_type" = background_cell_type,
                           "plot_image" = plot_image)
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Need to over-sample as cells which are too close will be removed later
  n_cells_inflated <- n_cells * 2
  
  spatial_df <- data.frame("Cell.X.Position" = runif(n_cells_inflated, 0, length),
                           "Cell.Y.Position" = runif(n_cells_inflated, 0, width),
                           "Cell.Z.Position" = runif(n_cells_inflated, 0, height))
  
  # Add integer rownames to data frame - each cell is labelled by an integer
  rownames(spatial_df) <- seq(nrow(spatial_df))
  
  ### Check if all other cells are to close to the current cell
  # Use frNN function: for each point, get all points within min_d of it
  spatial_df_distances <- dbscan::frNN(spatial_df,
                                       eps = minimum_distance_between_cells,
                                       query = NULL,
                                       sort = FALSE)
  
  # For each cell, get all other cells which were within 'minimum_distance_between_cells'
  spatial_df_distances_ids <- spatial_df_distances$id
  
  # Filter out zero length cells
  spatial_df_distances_ids <- Filter(function(x) length(x) != 0, spatial_df_distances_ids)
  
  # Get integer labels for the remaining cells
  spatial_df_distances_ids_names <- as.integer(names(spatial_df_distances_ids))
  
  # Determine which cells should be chosen from spatial_df
  cells_chosen <- rep(T, nrow(spatial_df))
  for (i in seq_len(length(spatial_df_distances_ids))) {
    cells_too_close <- spatial_df_distances_ids[[i]]
    
    if (cells_chosen[spatial_df_distances_ids_names[i]]) cells_chosen[cells_too_close] <- F
  }
  
  spatial_df <- spatial_df[cells_chosen, ]
  
  # If number of cells remaining is still higher than n_cells, randomly subset n_cells
  if (nrow(spatial_df) > n_cells) {
    spatial_df <- dplyr::sample_n(spatial_df, n_cells)
  }
  else {
    warning("The number of cells in your object could not reach the desired number of cells, consider lowering the number of cells, the minimum distance between cells, or increasing length, width or height.")
  }
  
  # Add Cell.Type and Cell.ID
  spatial_df$Cell.Type <- background_cell_type
  spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")
  
  # Get metadata
  background_metadata <- list("background_type" = "hardcore",
                              "n_cells" = n_cells,
                              "length" = length,
                              "width" = width,
                              "height" = height,
                              "minimum_distance_between_cells" = minimum_distance_between_cells,
                              "cell_types" = background_cell_type,
                              "cell_proportions" = 1)
  simulation_metadata <- list(background = background_metadata)
  
  ## Convert data frame to spe object
  spe <- SpatialExperiment::SpatialExperiment(
    assay = matrix(data = NA, nrow = 0, ncol = nrow(spatial_df)),
    colData = spatial_df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = list(simulation = simulation_metadata))
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe,
                        background_cell_type,
                        "lightgray")
    methods::show(fig)
  }
  
  return(spe)
}

simulate_ordered_background_cells3D <- function(n_cells,
                                                length,
                                                width,
                                                height,
                                                jitter_proportion,
                                                background_cell_type = "Others",
                                                plot_image = TRUE) {
  
  # Check input parameters
  input_parameters <- list("n_cells" = n_cells,
                           "length" = length,
                           "width" = width,
                           "height" = height,
                           "jitter_proportion" = jitter_proportion,
                           "background_cell_type" = background_cell_type,
                           "plot_image" = plot_image)
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Obtain distance between each point using MAGIC formula
  d_cells <- ((sqrt(2) * length * width * height)/n_cells)^(1/3)
  
  # Get distance between rows, columns and layers using 'd_cells'
  d_rows <- d_cells
  d_cols <- (sqrt(3) / 2) * d_cells
  d_lays <- (sqrt(6) / 3) * d_cells
  
  # Get number of rows, columns and layers
  n_rows <- round(length / d_rows)
  n_cols <- round(width / d_cols)
  n_lays <- round(height / d_lays)
  
  # Step 0. Assume points are on a 3D rectangular grid
  rows <- rep(seq(n_rows), n_cols * n_lays) * d_rows
  cols <- rep(rep(seq(n_cols), each = n_rows), n_lays) * d_cols
  lays <- rep(seq(n_lays), each = n_rows * n_cols) * d_lays
  
  # Step 1. For every odd sheet, every even row shifts by d_cells/2 right
  if (n_cols %% 2 == 0) {
    shift <- rep(c(rep(0, n_rows), rep(d_cells/2, n_rows)), n_cols/2)
  }
  else {
    shift <- c(rep(c(rep(0, n_rows), rep(d_cells/2, n_rows)), n_cols/2), rep(0, n_rows))
  }
  rows <- rows + c(shift, rep(0, n_rows * n_cols)) # Shift each even row by d_cells/2 right
  
  # Step 2. For every even sheet, odd rows shift d_cells/2 right, all rows shift d_cells/(2*sqrt(3)) up
  if (n_cols %% 2 == 0) {
    shift <- rep(c(rep(d_cells/2, n_rows), rep(0, n_rows)), n_cols/2)
  }
  else {
    shift <- c(rep(c(rep(d_cells/2, n_rows), rep(0, n_rows)), n_cols/2), rep(d_cells/2, n_rows))
  }
  rows <- rows + c(rep(0, n_rows * n_cols), shift) # Shift each odd row by d_cells/2 right
  cols <- cols + rep(c(0, d_cells/(2 * sqrt(3))), each = n_rows * n_cols) # Shift all rows by d_cells/(2*sqrt(3)) up
  
  # Get total number of cells (should be roughly equal to n_cells)
  n_total <- n_rows * n_cols * n_lays
  
  # Add randomness to the location of the cells
  jitter <- jitter_proportion * d_cells # Jitter is proportional to distance between points in hexagonal grid
  jitter_row <- runif(n_total, -jitter, jitter)
  jitter_col <- runif(n_total, -jitter, jitter)
  jitter_lay <- runif(n_total, -jitter, jitter)
  
  rows <- rows + jitter_row
  cols <- cols + jitter_col
  lays <- lays + jitter_lay
  
  # Put data into data frame
  df <- data.frame("Cell.X.Position" = rows,
                   "Cell.Y.Position" = cols,
                   "Cell.Z.Position" = lays,
                   "Cell.Type" = background_cell_type)
  df$Cell.ID <- paste("Cell", seq(nrow(df)), sep = "_")
  
  # Get metadata
  background_metadata <- list("background_type" = "ordered",
                              "n_cells" = n_cells,
                              "length" = length,
                              "width" = width,
                              "height" = height,
                              "jitter_proportion" = jitter_proportion,
                              "cell_types" = background_cell_type,
                              "cell_proportions" = 1)
  simulation_metadata <- list(background = background_metadata)
  
  ## Convert data frame to spe object
  spe <- SpatialExperiment::SpatialExperiment(
    assay = matrix(data = NA, nrow = 0, ncol = nrow(df)),
    colData = df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = list(simulation = simulation_metadata))
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe,
                        background_cell_type,
                        "lightgray")
    methods::show(fig)
  }
  
  return(spe)
}

spe_metadata_background_template <- function(background_type,
                                             original_spe_metadata = NULL) {
  
  if (background_type == "hardcore") {
    background_metadata <- list(background = list(background_type = "hardcore",
                                                  n_cells = 20000,
                                                  length = 600,
                                                  width = 600,
                                                  height = 300,
                                                  minimum_distance_between_cells = 10,
                                                  cell_types = c("Tumour", "Others"),
                                                  cell_proportions = c(0.05, 0.95)))
  }
  else if (background_type == "ordered") {
    background_metadata <- list(background = list(background_type = "ordered",
                                                  n_cells = 20000,
                                                  length = 600,
                                                  width = 300,
                                                  height = 300,
                                                  jitter_proportion = 0.25,
                                                  cell_types = c("Immune", "Others"),
                                                  cell_proportions = c(0.05, 0.95)))
  }
  else if (background_type == "uniform") {
    background_metadata <- list(background = list(background_type = "uniform",
                                                  n_cells = 20000,
                                                  length = 600,
                                                  width = 600,
                                                  height = 300,
                                                  cell_types = c("Tumour", "Others"),
                                                  cell_proportions = c(0.05, 0.95)))
  }
  else if (background_type == "poisson") {
    background_metadata <- list(background = list(background_type = "poisson",
                                                  n_cells = 20000,
                                                  length = 600,
                                                  width = 300,
                                                  height = 300,
                                                  cell_types = c("Immune", "Others"),
                                                  cell_proportions = c(0.05, 0.95)))
  }
  else {
    stop("background_type parameter is invalid.")
  }
  
  
  # If original_spe_metadata input is not null, replace its background metadata with new background metadata
  if (!is.null(original_spe_metadata) && !is.null(original_spe_metadata[["background"]])) {
    original_spe_metadata[["background"]] <- background_metadata
    return(original_spe_metadata)
  }
  else if (!is.null(original_spe_metadata) && is.null(original_spe_metadata[["background"]])) {
    original_spe_metadata <- c(background_metadata, original_spe_metadata)
    return(original_spe_metadata)
  }
  
  # Else, just return the background_metadata
  return(background_metadata)
}

# Running the functions ------------
updated_parameters_df <- generate_updated_parameters_for_arrangement_and_shape(
  n_simulations_per_parameter = 10,
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


poisson_metric_df_list <- analyse_SD2_simulations(updated_parameters_df, "poisson")
ordered_metric_df_list <- analyse_SD2_simulations(updated_parameters_df, "ordered")
uniform_metric_df_list <- analyse_SD2_simulations(updated_parameters_df, "uniform")
hardcore_metric_df_list <- analyse_SD2_simulations(updated_parameters_df, "hardcore")

