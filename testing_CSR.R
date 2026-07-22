# Updated functions to generate different background distributions -----
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

# Generate background
poisson_bg <- spe_metadata_background_template("poisson")
poisson_spe <- simulate_spe_metadata3D(poisson_bg, plot_image = F)

hardcore_bg <- spe_metadata_background_template("hardcore")
hardcore_spe <- simulate_spe_metadata3D(hardcore_bg, plot_image = F)

ordered_bg <- spe_metadata_background_template("ordered")
ordered_spe <- simulate_spe_metadata3D(ordered_bg, plot_image = F)

uniform_bg <- spe_metadata_background_template("uniform")
uniform_spe <- simulate_spe_metadata3D(uniform_bg, plot_image = F)

# Check if it follows CSR
library(spaSim3D)
library(SPIAT3D)
library(SpatialExperiment)
library(FNN)
library(ggplot2)
 
check_CSR_with_nnd3D <- function(spe) {
  # Using nearest neighbour distance
  # If clustered: nn distance is smaller than expected
  # If CSR: nn distance matches expected
  # If regular: nn distance is larger than expected
  
  coords <- spatialCoords(spe)
  colnames(coords) <- c("x", "y", "z")
  
  # Perform nearest neighbour test on spe
  nn <- get.knn(coords, k = 1)
  nn.dist <- nn$nn.dist
  mean.obs <- mean(nn.dist)
  
  
  # Perform nearest neighbour test on 999 simulated CSR distributions
  n <- nrow(coords)
  
  xrange <- range(coords[ , "x"])
  yrange <- range(coords[ , "y"])
  zrange <- range(coords[ , "z"])
  
  B <- 999
  
  mean.sim <- replicate(B, {
    sim <- cbind(
      runif(n, xrange[1], xrange[2]),
      runif(n, yrange[1], yrange[2]),
      runif(n, zrange[1], zrange[2])
    )
    mean(get.knn(sim, k = 1)$nn.dist)
  })
  
  # Rank of observed statistic among simulations
  all.stats <- c(mean.sim, mean.obs)
  
  r <- rank(all.stats, ties.method = "random")[length(all.stats)]
  
  p.low  <- r / (B + 1)
  p.high <- (B + 2 - r) / (B + 1)
  
  p.value <- min(1, 2 * min(p.low, p.high))
  
  return(p.value)
}

check_CSR_with_ripleyK3D <- function(spe) {

  # Make background cells only contain one cell type
  spe$Cell.Type <- "Others"
  
  # Perform ripleyK on spe
  K.obs <- SPIAT3D::calculate_cross_K_gradient3D(spe = spe, 
                                                 reference_cell_type = "Others", 
                                                 target_cell_type = "Others", 
                                                 radii = seq(1, 150, length.out = 11), 
                                                 feature_colname = "Cell.Type", 
                                                 plot_image = F)
  
  
  # Perform ripleyK test on 999 simulated CSR distributions
  B <- 999
  
  K.sim <- vector("list", B)
  
  for (i in seq_len(B)) {
    
    sim <- data.frame(
      x = runif(n, xrange[1], xrange[2]),
      y = runif(n, yrange[1], yrange[2]),
      z = runif(n, zrange[1], zrange[2]),
      Cell.Type = "Others"
    )
    
    sim_spe <- SpatialExperiment::SpatialExperiment(
      assay = matrix(data = NA, nrow = 0, ncol = nrow(sim)),
      colData = sim,
      spatialCoordsNames = c("x", "y", "z"))
    
    K.sim[[i]] <-   SPIAT3D::calculate_cross_K_gradient3D(spe = sim_spe, 
                                                          reference_cell_type = "Others", 
                                                          target_cell_type = "Others", 
                                                          radii = seq(0, 150, length.out = 51), 
                                                          feature_colname = "Cell.Type", 
                                                          plot_image = F)
  }
  
  K.mat <- sapply(K.sim, function(x) x$Others)
  
  lower <- apply(K.mat, 1, quantile, probs = 0.025)
  upper <- apply(K.mat, 1, quantile, probs = 0.975)
  
  envelope <- data.frame(
    radius = K.obs$radius,
    observed = K.obs$Others,
    lower = lower,
    upper = upper
  )
  
  fig <- ggplot(envelope, aes(radius)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                fill = "grey80") +
    geom_line(aes(y = observed),
              colour = "red",
              linewidth = 1)
  
  methods::show(fig)
}

poisson_p_value <- check_CSR_with_nnd3D(poisson_spe)
hardcore_p_value <- check_CSR_with_nnd3D(hardcore_spe)
ordered_p_value <- check_CSR_with_nnd3D(ordered_spe)
uniform_p_value <- check_CSR_with_nnd3D(uniform_spe)
