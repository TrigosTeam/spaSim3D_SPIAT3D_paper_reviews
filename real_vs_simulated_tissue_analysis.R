library(spaSim3D)
library(SPIAT3D)
# Functions ----
convert_df_to_spe <- function(df) {
  spe <- SpatialExperiment::SpatialExperiment(
    assay = matrix(data = NA, nrow = 0, ncol = nrow(df)),
    colData = df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"))
  
  return(spe)
}

analyse_3D_data_real_and_sim <- function(
    data3D_real_spe,
    data3D_sim_spe,
    cell_types,
    radii = seq(20, 100, 10),
    n_splits = 10,
    thresholds = seq(0.01, 1, 0.01)
) {
  
  n_cell_type_combinations <- length(cell_types)^2
  radii_colnames <- paste("r", radii, sep = "")
  thresholds_colnames <- paste("t", thresholds, sep = "")
  
  create_empty_metric_df_list <- function(
    cell_types,
    radii_colnames,
    thresholds_colnames
  ) {
    n_cell_type_combinations <- length(cell_types)^2
    
    # Define AMD data frames as well as constants
    AMD_df_colnames <- c("reference", "target", "AMD")
    AMD_df <- data.frame(matrix(nrow = 2 * n_cell_type_combinations, ncol = length(AMD_df_colnames)))
    colnames(AMD_df) <- AMD_df_colnames
    
    # Define MS, NMS, ANC, ACIN, CKR, CLR, CGR, COO, ANE data frames as well as constants
    radii_colnames <- paste("r", radii, sep = "")
    
    MS_df_colnames <- c("reference", "target", radii_colnames)
    MS_df <- data.frame(matrix(nrow = 2 * n_cell_type_combinations, ncol = length(MS_df_colnames)))
    colnames(MS_df) <- MS_df_colnames
    
    NMS_df <- ANC_df <- ANE_df <- ACIN_df <- CKR_df <- CLR_df <- COO_df <- CGR_df <- CK_df <- CL_df <- CG_df <- MS_df
    
    # Define SAC and prevalence data frames as well as constants
    thresholds_colnames <- paste("t", thresholds, sep = "")
    
    PBSAC_df_colnames <- c("reference", "target", "PBSAC")
    PBSAC_df <- data.frame(matrix(nrow = 2 * n_cell_type_combinations, ncol = length(PBSAC_df_colnames)))
    colnames(PBSAC_df) <- PBSAC_df_colnames
    
    PBP_df_colnames <- c("reference", "target", thresholds_colnames)
    PBP_df <- data.frame(matrix(nrow = 2 * n_cell_type_combinations, ncol = length(PBP_df_colnames)))
    colnames(PBP_df) <- PBP_df_colnames
    
    EBSAC_df_colnames <- c("cell_types", "EBSAC")
    EBSAC_df <- data.frame(matrix(nrow = 2 * n_cell_type_combinations, ncol = length(EBSAC_df_colnames)))
    colnames(EBSAC_df) <- EBSAC_df_colnames
    
    EBP_df_colnames <- c("cell_types", thresholds_colnames)
    EBP_df <- data.frame(matrix(nrow = 2 * n_cell_type_combinations, ncol = length(EBP_df_colnames)))
    colnames(EBP_df) <- EBP_df_colnames
    
    
    # Add all to list:
    metric_df_list <- list(AMD = AMD_df,
                           MS = MS_df,
                           NMS = NMS_df,
                           ACIN = ACIN_df,
                           ANE = ANE_df,
                           ANC = ANC_df,
                           CKR = CKR_df,
                           CLR = CLR_df,
                           CGR = CGR_df,
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
  
  metric_df_list <- create_empty_metric_df_list(
    cell_types,
    radii_colnames,
    thresholds_colnames
  )
  
  
  for (i in c(1, 2)) {
    print(i)
    
    if (i == 1) {
      spe <- data3D_real_spe
    }
    else if (i == 2) {
      spe <- data3D_sim_spe
    }
    
    ### 3D analysis -----------------------------
    index <- n_cell_type_combinations * (i - 1) + 1 
    
    minimum_distance_data <- calculate_minimum_distances_between_cell_types3D(spe,
                                                                              cell_types,
                                                                              feature_colname = "Cell.Type",
                                                                              show_summary = F,
                                                                              plot_image = F)
    
    minimum_distance_data_summary <- summarise_distances_between_cell_types3D(minimum_distance_data)
    
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "reference"] <- minimum_distance_data_summary$reference
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "target"] <- minimum_distance_data_summary$target
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "AMD"] <- minimum_distance_data_summary$mean
    
    # Need a new index for gradient-based data which increments after each target cell type
    pair_index <- n_cell_type_combinations * (i - 1) + 1 
    
    for (reference_cell_type in cell_types) {
      gradient_data <- calculate_all_gradient_cc_metrics3D(spe,
                                                           reference_cell_type,
                                                           cell_types,
                                                           radii,
                                                           feature_colname = "Cell.Type",
                                                           plot_image = F)
      
      for (target_cell_type in cell_types) {
        print(paste(reference_cell_type, target_cell_type, sep = "/"))
        metric_df_list[["ANC"]][pair_index, c("reference", "target")] <- c(reference_cell_type, target_cell_type)
        metric_df_list[["ACIN"]][pair_index, c("reference", "target")] <- c(reference_cell_type, target_cell_type)
        metric_df_list[["CKR"]][pair_index, c("reference", "target")] <- c(reference_cell_type, target_cell_type)
        metric_df_list[["CLR"]][pair_index, c("reference", "target")] <- c(reference_cell_type, target_cell_type)
        metric_df_list[["COO"]][pair_index, c("reference", "target")] <- c(reference_cell_type, target_cell_type)
        metric_df_list[["CGR"]][pair_index, c("reference", "target")] <- c(reference_cell_type, target_cell_type)
        metric_df_list[["CK"]][pair_index, c("reference", "target")] <- c(reference_cell_type, target_cell_type)
        metric_df_list[["CL"]][pair_index, c("reference", "target")] <- c(reference_cell_type, target_cell_type)
        metric_df_list[["CG"]][pair_index, c("reference", "target")] <- c(reference_cell_type, target_cell_type)
        metric_df_list[["MS"]][pair_index, c("reference", "target")] <- c(reference_cell_type, target_cell_type)
        metric_df_list[["NMS"]][pair_index, c("reference", "target")] <- c(reference_cell_type, target_cell_type)
        metric_df_list[["ANE"]][pair_index, c("reference", "target")] <- c(reference_cell_type, 
                                                                           paste(reference_cell_type, target_cell_type, sep = ","))
        
        
        if (is.null(gradient_data)) {
          metric_df_list[["ANC"]][pair_index, radii_colnames] <- NA
          metric_df_list[["CKR"]][pair_index, radii_colnames] <- NA
          metric_df_list[["CLR"]][pair_index, radii_colnames] <- NA
          metric_df_list[["COO"]][pair_index, radii_colnames] <- NA
          metric_df_list[["CGR"]][pair_index, radii_colnames] <- NA
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
          metric_df_list[["CKR"]][pair_index, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]] / gradient_data[["cross_K"]][["expected"]]
          metric_df_list[["CLR"]][pair_index, radii_colnames] <- gradient_data[["cross_L"]][[target_cell_type]] / gradient_data[["cross_L"]][["expected"]]
          metric_df_list[["COO"]][pair_index, radii_colnames] <- gradient_data[["co_occurrence"]][[target_cell_type]]
          metric_df_list[["CGR"]][pair_index, radii_colnames] <- gradient_data[["cross_G"]][[target_cell_type]][["observed_cross_G"]] / gradient_data[["cross_G"]][[target_cell_type]][["expected_cross_G"]]
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
        metric_df_list[["PBSAC"]][pair_index, c("reference", "target")] <- c(reference_cell_type, target_cell_type)
        metric_df_list[["PBP"]][pair_index, c("reference", "target")] <- c(reference_cell_type, target_cell_type)
        
        metric_df_list[["EBSAC"]][pair_index, c("cell_types")] <- c(paste(reference_cell_type, target_cell_type, sep = ","))
        metric_df_list[["EBP"]][pair_index, c("cell_types")] <- c(paste(reference_cell_type, target_cell_type, sep = ","))
        
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
  }
  return(metric_df_list)
}

### 1. Download and simulate mouse hypothalamus ----

# Original
setwd("/home/dle/R/data3D/merfish_squidpy")
merfish_squidpy_df <- read.csv("merfish_squidpy_df.csv")

merfish_squidpy_df$Cell.X.Position <- merfish_squidpy_df$Cell.X.Position * 1.8
merfish_squidpy_df$Cell.Y.Position <- merfish_squidpy_df$Cell.Y.Position * 1.8
merfish_squidpy_df$Cell.Z.Position <- ((as.integer(merfish_squidpy_df$Cell.Z.Position) + 29) / 100)

mouse_hypothalamus_real_spe <- convert_df_to_spe(merfish_squidpy_df)

SPIAT3D::plot_cells3D(mouse_hypothalamus_real_spe, feature_colname = "Cell.Type")


# using spaSim-3D
mouse_hypothalamus_md <- spe_metadata_background_template("random")
mouse_hypothalamus_md <- spe_metadata_cluster_template("regular", "cylinder", mouse_hypothalamus_md)
mouse_hypothalamus_md <- spe_metadata_cluster_template("regular", "cylinder", mouse_hypothalamus_md)
mouse_hypothalamus_md <- spe_metadata_cluster_template("regular", "cylinder", mouse_hypothalamus_md)
mouse_hypothalamus_md <- spe_metadata_cluster_template("regular", "cylinder", mouse_hypothalamus_md)
mouse_hypothalamus_md <- spe_metadata_cluster_template("regular", "cylinder", mouse_hypothalamus_md)
mouse_hypothalamus_md <- spe_metadata_cluster_template("regular", "cylinder", mouse_hypothalamus_md)
mouse_hypothalamus_md <- spe_metadata_cluster_template("regular", "network", mouse_hypothalamus_md)
mouse_hypothalamus_md <- spe_metadata_cluster_template("regular", "ellipsoid", mouse_hypothalamus_md)


mouse_hypothalamus_md$background$n_cells <- 70000
mouse_hypothalamus_md$background$length <- 100
mouse_hypothalamus_md$background$width <- 100
mouse_hypothalamus_md$background$height <- 50
mouse_hypothalamus_md$background$minimum_distance_between_cells <- 0
mouse_hypothalamus_md$background$cell_types <- 
  c("Ambiguous", "Astrocyte", "Endothelial", 
    "Ependymal", "Excitatory", "Inhibitory")
mouse_hypothalamus_md$background$cell_proportions <-
  c(0.125, 0.125, 0,
    0, 0.25, 0.5)

mouse_hypothalamus_md$cluster_1$cluster_cell_types <- "Ependymal"
mouse_hypothalamus_md$cluster_1$cluster_cell_proportions <- 1
mouse_hypothalamus_md$cluster_1$radius <- 3
mouse_hypothalamus_md$cluster_1$start_loc <- c(10, 50, 47)
mouse_hypothalamus_md$cluster_1$end_loc <- c(100, 50, 47)

mouse_hypothalamus_md$cluster_2$cluster_cell_types <- "Ependymal"
mouse_hypothalamus_md$cluster_2$cluster_cell_proportions <- 1
mouse_hypothalamus_md$cluster_2$radius <- 3
mouse_hypothalamus_md$cluster_2$start_loc <- c(10, 50, 44)
mouse_hypothalamus_md$cluster_2$end_loc <- c(90, 50, 44)

mouse_hypothalamus_md$cluster_3$cluster_cell_types <- "Ependymal"
mouse_hypothalamus_md$cluster_3$cluster_cell_proportions <- 1
mouse_hypothalamus_md$cluster_3$radius <- 3
mouse_hypothalamus_md$cluster_3$start_loc <- c(10, 50, 41)
mouse_hypothalamus_md$cluster_3$end_loc <- c(90, 50, 41)

mouse_hypothalamus_md$cluster_4$cluster_cell_types <- "Ependymal"
mouse_hypothalamus_md$cluster_4$cluster_cell_proportions <- 1
mouse_hypothalamus_md$cluster_4$radius <- 3
mouse_hypothalamus_md$cluster_4$start_loc <- c(10, 50, 38)
mouse_hypothalamus_md$cluster_4$end_loc <- c(100, 50, 38)

mouse_hypothalamus_md$cluster_5$cluster_cell_types <- "Ependymal"
mouse_hypothalamus_md$cluster_5$cluster_cell_proportions <- 1
mouse_hypothalamus_md$cluster_5$radius <- 3
mouse_hypothalamus_md$cluster_5$start_loc <- c(10, 50, 35)
mouse_hypothalamus_md$cluster_5$end_loc <- c(90, 50, 35)

mouse_hypothalamus_md$cluster_6$cluster_cell_types <- "Ependymal"
mouse_hypothalamus_md$cluster_6$cluster_cell_proportions <- 1
mouse_hypothalamus_md$cluster_6$radius <- 3
mouse_hypothalamus_md$cluster_6$start_loc <- c(10, 50, 32)
mouse_hypothalamus_md$cluster_6$end_loc <- c(90, 50, 32)


mouse_hypothalamus_md$cluster_7$cluster_cell_types <- "Endothelial"
mouse_hypothalamus_md$cluster_7$cluster_cell_proportions <- 1
mouse_hypothalamus_md$cluster_7$n_edges <- 30
mouse_hypothalamus_md$cluster_7$width <- 6
mouse_hypothalamus_md$cluster_7$centre_loc <- c(80, 50, 0)
mouse_hypothalamus_md$cluster_7$radius <- 40

mouse_hypothalamus_md$cluster_8$cluster_cell_types <- 
  c("Ambiguous", "Astrocyte", "Endothelial", 
    "Ependymal", "Excitatory", "Inhibitory")
mouse_hypothalamus_md$cluster_8$cluster_cell_proportions <-
  c(0.125, 0, 0,
    0, 0.75, 0.125)

mouse_hypothalamus_md$cluster_8$radii <- c(20, 30, 20)
mouse_hypothalamus_md$cluster_8$centre_loc <- c(0, 50, 40)
mouse_hypothalamus_md$cluster_8$axes_rotation <- c(0, 0, 0)

mouse_hypothalamus_sim_spe <- simulate_spe_metadata3D(mouse_hypothalamus_md, plot_image = F)

# Update coords to match real data
SpatialExperiment::spatialCoords(mouse_hypothalamus_sim_spe)[ , "Cell.X.Position"] <- SpatialExperiment::spatialCoords(mouse_hypothalamus_sim_spe)[ , "Cell.X.Position"] * 1.8 / 100
SpatialExperiment::spatialCoords(mouse_hypothalamus_sim_spe)[ , "Cell.Y.Position"] <- SpatialExperiment::spatialCoords(mouse_hypothalamus_sim_spe)[ , "Cell.Y.Position"] * 1.8 / 100
SpatialExperiment::spatialCoords(mouse_hypothalamus_sim_spe)[ , "Cell.Z.Position"] <- SpatialExperiment::spatialCoords(mouse_hypothalamus_sim_spe)[ , "Cell.Z.Position"] / 100

SPIAT3D::plot_cells3D(mouse_hypothalamus_sim_spe, feature_colname = "Cell.Type")


### 2. Download and simulate Drosphilia Egg (E16_18h) ----
# Original
setwd("/home/dle/R/data3D/stomics_fly")
E16_18h_df <- read.csv("E16-18h_df.csv")
E16_18h_df$Cell.X.Position <- E16_18h_df$Cell.X.Position * 10 + 265
E16_18h_df$Cell.Y.Position <- E16_18h_df$Cell.Y.Position * 10 + 160
E16_18h_df$Cell.Z.Position <- E16_18h_df$Cell.Z.Position * 10

drosphilia_egg_real_spe <- convert_df_to_spe(E16_18h_df)

SPIAT3D::plot_cells3D(drosphilia_egg_real_spe, 
                      feature_colname = "Cell.Type",
                      plot_cell_types = 
                        c("salivary gland", "epidermis", "CNS", "carcass", "fat body",
                          "muscle", "trachea", "midgut", "hemolymph", "foregut"),
                      plot_colours = c(
                        "#33135b",
                        "#73ec81",
                        "#9437a8",
                        "#b8db50",
                        "#0062c5",
                        "#f77e3b",
                        "#48bbff",
                        "#bb0036",
                        "#4deeac",
                        "#770026"
                      ))


# using spasim_3D
egg_md <- spe_metadata_background_template("ordered")
egg_md <- spe_metadata_cluster_template("double ring", "ellipsoid", egg_md)
egg_md <- spe_metadata_cluster_template("regular", "sphere", egg_md)
egg_md <- spe_metadata_cluster_template("regular", "sphere", egg_md)
egg_md <- spe_metadata_cluster_template("regular", "sphere", egg_md)
egg_md <- spe_metadata_cluster_template("regular", "ellipsoid", egg_md)
egg_md <- spe_metadata_cluster_template("regular", "ellipsoid", egg_md)
egg_md <- spe_metadata_cluster_template("regular", "ellipsoid", egg_md)
egg_md <- spe_metadata_cluster_template("regular", "sphere", egg_md)

egg_md$background$n_cells <- 15000
egg_md$background$length <- 100
egg_md$background$width <- 70
egg_md$background$height <- 40
egg_md$background$jitter_proportion <- 0.1
egg_md$background$cell_types <- 
  c("salivary gland", "epidermis", "CNS", "carcass", "fat body",
    "muscle", "trachea", "midgut", "hemolymph", "foregut",
    "other")
egg_md$background$cell_proportions <-
  c(0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    1)

egg_md$cluster_1$cluster_cell_types <- "midgut"
egg_md$cluster_1$cluster_cell_proportions <- 1
egg_md$cluster_1$centre_loc <- c(50, 50, 25)
egg_md$cluster_1$radii <- c(34, 14, 8)
egg_md$cluster_1$axes_rotation <- c(0, 0, 0)
egg_md$cluster_1$inner_ring_cell_types <- "fat body"
egg_md$cluster_1$inner_ring_cell_proportions <- 1
egg_md$cluster_1$inner_ring_width <- 4
egg_md$cluster_1$outer_ring_cell_types <- c("epidermis", "carcass", "muscle")
egg_md$cluster_1$outer_ring_cell_proportions <- c(0.8, 0.1, 0.1)
egg_md$cluster_1$outer_ring_width <- 4

egg_md$cluster_2$cluster_cell_types <- "muscle"
egg_md$cluster_2$cluster_cell_proportions <- 1
egg_md$cluster_2$radius <- 14
egg_md$cluster_2$centre_loc <- c(30, 50, 25)

egg_md$cluster_3$cluster_cell_types <- "salivary gland"
egg_md$cluster_3$cluster_cell_proportions <- 1
egg_md$cluster_3$radius <- 10
egg_md$cluster_3$centre_loc <- c(30, 60, 25)

egg_md$cluster_4$cluster_cell_types <- "midgut"
egg_md$cluster_4$cluster_cell_proportions <- 1
egg_md$cluster_4$radius <- 8
egg_md$cluster_4$centre_loc <- c(70, 50, 40)

egg_md$cluster_5$cluster_cell_types <- "carcass"
egg_md$cluster_5$cluster_cell_proportions <- 1
egg_md$cluster_5$centre_loc <- c(50, 65, 25)
egg_md$cluster_5$radii <- c(32, 5, 8)
egg_md$cluster_5$axes_rotation <- c(0, 0, 0)

egg_md$cluster_6$cluster_cell_types <- "carcass"
egg_md$cluster_6$cluster_cell_proportions <- 1
egg_md$cluster_6$centre_loc <- c(50, 35, 25)
egg_md$cluster_6$radii <- c(32, 5, 8)
egg_md$cluster_6$axes_rotation <- c(0, 0, 0)

egg_md$cluster_7$cluster_cell_types <- c("trachea")
egg_md$cluster_7$cluster_cell_proportions <- c(1)
egg_md$cluster_7$centre_loc <- c(75, 60, 32)
egg_md$cluster_7$radii <- c(5, 12, 6)
egg_md$cluster_7$axes_rotation <- c(45, 0, 0)

egg_md$cluster_8$cluster_cell_types <- "foregut"
egg_md$cluster_8$cluster_cell_proportions <- 1
egg_md$cluster_8$radius <- 6
egg_md$cluster_8$centre_loc <- c(80, 50, 32)

drosphilia_egg_sim_spe <- simulate_spe_metadata3D(egg_md, plot_image = F)

# Update coords to match real data
SpatialExperiment::spatialCoords(drosphilia_egg_sim_spe)[ , "Cell.X.Position"] <- SpatialExperiment::spatialCoords(drosphilia_egg_sim_spe)[ , "Cell.X.Position"] * 5
SpatialExperiment::spatialCoords(drosphilia_egg_sim_spe)[ , "Cell.Y.Position"] <- SpatialExperiment::spatialCoords(drosphilia_egg_sim_spe)[ , "Cell.Y.Position"] * 4.5
SpatialExperiment::spatialCoords(drosphilia_egg_sim_spe)[ , "Cell.Z.Position"] <- SpatialExperiment::spatialCoords(drosphilia_egg_sim_spe)[ , "Cell.Z.Position"] * 2

plot_cells3D(drosphilia_egg_sim_spe, 
             feature_colname = "Cell.Type",
             plot_cell_types = 
               c("salivary gland", "epidermis", "CNS", "carcass", "fat body",
                 "muscle", "trachea", "midgut", "hemolymph", "foregut"),
             plot_colours = c(
               "#33135b",
               "#73ec81",
               "#9437a8",
               "#b8db50",
               "#0062c5",
               "#f77e3b",
               "#48bbff",
               "#bb0036",
               "#4deeac",
               "#770026"
             ))



### 3. Download and simulate Mouse embryo (E11.5h) ----
# Original
setwd("/home/dle/R/data3D/spateo")
mouse_E11.5_embryo <- read.csv("mouse_E11.5_embryo.csv")
mouse_embryo_cells <- names(sort(table(mouse_E11.5_embryo$Cell.Type), decreasing = TRUE))[1:12]

mouse_E11.5_embryo$Cell.X.Position <- mouse_E11.5_embryo$Cell.X.Position - 580
mouse_E11.5_embryo$Cell.Y.Position <- mouse_E11.5_embryo$Cell.Y.Position - 38
mouse_E11.5_embryo$Cell.Z.Position <- mouse_E11.5_embryo$Cell.Z.Position - 305

mouse_embryo_real_spe <- convert_df_to_spe(mouse_E11.5_embryo)

SPIAT3D::plot_cells3D(mouse_embryo_real_spe, 
                      feature_colname = "Cell.Type",
                      plot_cell_types = mouse_embryo_cells,
                      plot_colours = c(
                        "#33135b",
                        "#73ec81",
                        "#9437a8",
                        "#b8db50",
                        "#0062c5",
                        "#f77e3b",
                        "#48bbff",
                        "#bb0036",
                        "#4deeac",
                        "#770026",
                        "#007128",
                        "#d99dff"
                      ))


# using spaSim-3D
mouse_md <- spe_metadata_background_template("random")
mouse_md <- spe_metadata_cluster_template("regular", "sphere", mouse_md)
mouse_md <- spe_metadata_cluster_template("regular", "ellipsoid", mouse_md)
mouse_md <- spe_metadata_cluster_template("regular", "ellipsoid", mouse_md)
mouse_md <- spe_metadata_cluster_template("regular", "sphere", mouse_md)
mouse_md <- spe_metadata_cluster_template("regular", "sphere", mouse_md)
mouse_md <- spe_metadata_cluster_template("regular", "ellipsoid", mouse_md)
mouse_md <- spe_metadata_cluster_template("regular", "ellipsoid", mouse_md)
mouse_md <- spe_metadata_cluster_template("regular", "cylinder", mouse_md)
mouse_md <- spe_metadata_cluster_template("regular", "cylinder", mouse_md)
mouse_md <- spe_metadata_cluster_template("regular", "cylinder", mouse_md)


mouse_md$background$n_cells <- 7000000
mouse_md$background$length <- 100
mouse_md$background$width <- 100
mouse_md$background$height <- 50
mouse_md$background$minimum_distance_between_cells <- 0
mouse_md$background$cell_types <- "other"
mouse_md$background$cell_proportions <- 1

mouse_md$cluster_1$cluster_cell_types <- c("Glutamatergic neurons", "Cranial mesenchyme", "Cranial mesoderm")
mouse_md$cluster_1$cluster_cell_proportions <- c(0.3, 0.4, 0.3)
mouse_md$cluster_1$radius <- 25
mouse_md$cluster_1$centre_loc <- c(25, 25, 25)

mouse_md$cluster_2$cluster_cell_types <- c("Sclerotome", "Fibroblasts", "Glutamatergic neurons", "Lateral plate and intermediate mesoderm")
mouse_md$cluster_2$cluster_cell_proportions <- c(0.4, 0.3, 0.2, 0.1)
mouse_md$cluster_2$centre_loc <- c(60, 45, 25)
mouse_md$cluster_2$radii <- c(45, 20, 20)
mouse_md$cluster_2$axes_rotation <- c(-10, 0, 0)

mouse_md$cluster_3$cluster_cell_types <- c("Sclerotome", "Fibroblasts", "Glutamatergic neurons", "Lateral plate and intermediate mesoderm")
mouse_md$cluster_3$cluster_cell_proportions <- c(0.6, 0.20, 0.15, 0.05)
mouse_md$cluster_3$centre_loc <- c(80, 20, 25)
mouse_md$cluster_3$radii <- c(25, 15, 20)
mouse_md$cluster_3$axes_rotation <- c(38, 0, 0)

mouse_md$cluster_4$cluster_cell_types <- c("Limb progenitors and lateral plate mesenchyme")
mouse_md$cluster_4$cluster_cell_proportions <- c(1)
mouse_md$cluster_4$radius <- 9
mouse_md$cluster_4$centre_loc <- c(60, 45, 40)

mouse_md$cluster_5$cluster_cell_types <- c("Limb progenitors and lateral plate mesenchyme")
mouse_md$cluster_5$cluster_cell_proportions <- c(1)
mouse_md$cluster_5$radius <- 10
mouse_md$cluster_5$centre_loc <- c(75, 15, 40)

mouse_md$cluster_6$cluster_cell_types <- c("Facial mesenchyme")
mouse_md$cluster_6$cluster_cell_proportions <- c(1)
mouse_md$cluster_6$centre_loc <- c(35, 25, 30)
mouse_md$cluster_6$radii <- c(10, 18, 15)
mouse_md$cluster_6$axes_rotation <- c(0, 0, 0)

mouse_md$cluster_7$cluster_cell_types <- c("Telencephalon neuroectoderm")
mouse_md$cluster_7$cluster_cell_proportions <- c(1)
mouse_md$cluster_7$centre_loc <- c(25, 10, 25)
mouse_md$cluster_7$radii <- c(18, 10, 15)
mouse_md$cluster_7$axes_rotation <- c(-60, 0, 0)

mouse_md$cluster_8$cluster_cell_types <- c("Spinal cord neuroectoderm")
mouse_md$cluster_8$cluster_cell_proportions <- c(1)
mouse_md$cluster_8$radius <- 5
mouse_md$cluster_8$start_loc <- c(20, 60, 25)
mouse_md$cluster_8$end_loc <- c(65, 64, 25)


mouse_md$cluster_9$cluster_cell_types <- c("Spinal cord neuroectoderm")
mouse_md$cluster_9$cluster_cell_proportions <- c(1)
mouse_md$cluster_9$radius <- 5
mouse_md$cluster_9$start_loc <- c(65, 64, 25)
mouse_md$cluster_9$end_loc <- c(100, 42, 25)

mouse_md$cluster_10$cluster_cell_types <- c("Spinal cord neuroectoderm")
mouse_md$cluster_10$cluster_cell_proportions <- c(1)
mouse_md$cluster_10$radius <- 5
mouse_md$cluster_10$start_loc <- c(100, 40, 25)
mouse_md$cluster_10$end_loc <- c(100, 25, 25)

mouse_embryo_sim_spe <- simulate_spe_metadata3D(mouse_md, plot_image = F)
SpatialExperiment::spatialCoords(mouse_embryo_sim_spe)[ , "Cell.X.Position"] <- SpatialExperiment::spatialCoords(mouse_embryo_sim_spe)[ , "Cell.X.Position"] * 13
SpatialExperiment::spatialCoords(mouse_embryo_sim_spe)[ , "Cell.Y.Position"] <- SpatialExperiment::spatialCoords(mouse_embryo_sim_spe)[ , "Cell.Y.Position"] * 16
SpatialExperiment::spatialCoords(mouse_embryo_sim_spe)[ , "Cell.Z.Position"] <- SpatialExperiment::spatialCoords(mouse_embryo_sim_spe)[ , "Cell.Z.Position"] * 10

plot_cells3D(mouse_embryo_sim_spe,
             feature_colname = "Cell.Type",
             plot_cell_types = mouse_embryo_cells,
             plot_colours = c(
               "#33135b",
               "#73ec81",
               "#9437a8",
               "#b8db50",
               "#0062c5",
               "#f77e3b",
               "#48bbff",
               "#bb0036",
               "#4deeac",
               "#770026",
               "#007128",
               "#d99dff"
             ))



### 4. Analysing and comparing datasets ------

# Mouse hypothalamus
mouse_hypothalamus_real_spe$Cell.Type[mouse_hypothalamus_real_spe$Cell.Type %in% c("Endothelial 2")] <- "Endothelial"
mouse_hypothalamus_cell_types <- intersect(mouse_hypothalamus_real_spe$Cell.Type, mouse_hypothalamus_sim_spe$Cell.Type)

mouse_hypothalamus_metric_df_list <- 
  analyse_3D_data_real_and_sim(
    mouse_hypothalamus_real_spe,
    mouse_hypothalamus_sim_spe,
    mouse_hypothalamus_cell_types
  )

# Drosphilia Egg

# Mouse embryo
