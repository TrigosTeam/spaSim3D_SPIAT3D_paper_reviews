### spaSim3D testing
## parameters:
# - number of cells in simulation
# - size of clusters
# - shape of clusters
# - number of clusters
# - complexity of networks

plot_results_n_cells <- function(results, n_cells) {
  
  
  # Combine results
  df <- map2_dfr(results, n_cells, ~{
    tibble(
      size = .y,
      median_time = as.numeric(.x$median),      # seconds
      mem_alloc = as.numeric(.x$mem_alloc)      # bytes
    )
  })
  
  # Runtime plot
  p_runtime <- ggplot(df, aes(size, median_time)) +
    geom_line() +
    geom_point(size = 3) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      x = "Number of cells",
      y = "Median runtime (s)"
    ) +
    theme_minimal()
  
  # Memory plot
  p_memory <- ggplot(df, aes(size, mem_alloc / 1024^2)) +
    geom_line() +
    geom_point(size = 3) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      x = "Number of cells",
      y = "Memory allocated (MB)"
    ) +
    theme_minimal()
  
  # Side-by-side
  fig <- plot_grid(
    p_runtime,
    p_memory,
    ncol = 2,
    align = "hv"
  )
  
  methods::show(fig)
  
  return(df)
}

# Libraries ----
library(bench)
library(dplyr)
library(purrr)
library(ggplot2)
library(cowplot)
library(spaSim3D)
library(SPIAT3D)

# Testing simulate_random_background_cells3D function with increasing number of cells ----
n_cells <- c(1e3, 1e4, 1e5, 1e6, 1e7)

results <- lapply(n_cells, function(n) {
  
  mark(
    spaSim3D::simulate_random_background_cells3D(n_cells = n,
                                                 length = 100,
                                                 width = 100,
                                                 height = 100,
                                                 minimum_distance_between_cells = (5 / (2.5) ^ (log(n, 10) - 2)),
                                                 plot_image = FALSE),
    iterations = 10,
    check = F
  )
  
})

plot_results_n_cells(results, n_cells)

# Testing simulate_ordered_background_cells3D function with increasing number of cells ----
n_cells <- c(1e3, 1e4, 1e5, 1e6, 1e7)

results <- lapply(n_cells, function(n) {
  
  mark(
    spaSim3D::simulate_ordered_background_cells3D(n_cells = n,
                                                  length = 100,
                                                  width = 100,
                                                  height = 100,
                                                  jitter_proportion = 0,
                                                  plot_image = FALSE),
    iterations = 10,
    check = F
  )
  
})

plot_results_n_cells(results, n_cells)

# Testing simulating ellipsoid clusters with increasing number of cells in the background ----
n_cells <- c(1e3, 1e4, 1e5, 1e6, 1e7)

background_spes <- list()
i <- 1
for (n in n_cells) {
  background_spes[[i]] <- spaSim3D::simulate_random_background_cells3D(n_cells = n,
                                                                       length = 100,
                                                                       width = 100,
                                                                       height = 100,
                                                                       minimum_distance_between_cells = 0,
                                                                       plot_image = FALSE)
  i <- i + 1
}


indices <- seq_len(length(n_cells))

results <- lapply(indices, function(index) {
  
  mark(
    spaSim3D::simulate_clusters3D(background_spes[[index]],
                                  cluster_properties_list = list(
                                    C1 = list(
                                      shape = "ellipsoid",
                                      cluster_cell_types = c("A", "B", "O"),
                                      cluster_cell_proportions = c(0.65, 0.3, 0.05),
                                      radii = c(15, 20, 25),
                                      centre_loc = c(70, 70, 70),
                                      axes_rotation = c(0, 0, 0)
                                    )
                                  ),
                                  plot_image = FALSE),
    iterations = 10,
    check = F
  )
  
})

plot_results_n_cells(results, n_cells)



# Testing simulating network clusters with increasing number of cells in the background ----
n_cells <- c(1e3, 1e4, 1e5, 1e6, 1e7)

background_spes <- list()
i <- 1
for (n in n_cells) {
  background_spes[[i]] <- spaSim3D::simulate_random_background_cells3D(n_cells = n,
                                                                       length = 100,
                                                                       width = 100,
                                                                       height = 100,
                                                                       minimum_distance_between_cells = 0,
                                                                       plot_image = FALSE)
  i <- i + 1
}

indices <- seq_len(length(n_cells))

results <- lapply(indices, function(index) {
  
  mark(
    spaSim3D::simulate_clusters3D(background_spes[[index]],
                                  cluster_properties_list = list(
                                    C1 = list(
                                      shape = "network",
                                      cluster_cell_types = c("Immune"),
                                      cluster_cell_proportions = c(1),
                                      n_edges = 15,
                                      width = 8,
                                      centre_loc = c(50, 50, 50), # Rough centre of network cluster
                                      radius = 40 # Rough radius spanned by the network cluster
                                    )
                                  ),
                                  plot_image = FALSE),
    iterations = 10,
    check = F
  )
  
})

plot_results_n_cells(results, n_cells)

# Testing simulating ellipsoid clusters with increasing size of clusters ----
background_spe <- spaSim3D::simulate_random_background_cells3D(n_cells = 1e6,
                                                               length = 100,
                                                               width = 100,
                                                               height = 100,
                                                               minimum_distance_between_cells = 0,
                                                               plot_image = FALSE)



ellipsoid_sizes <- list(
  c(10, 15, 20),
  c(15, 20, 25),
  c(20, 25, 30),
  c(25, 30, 35),
  c(35, 40, 45)
)


results <- lapply(ellipsoid_sizes, function(size) {
  
  mark(
    spaSim3D::simulate_clusters3D(background_spe,
                                  cluster_properties_list = list(
                                    C1 = list(
                                      shape = "ellipsoid",
                                      cluster_cell_types = c("A", "B", "O"),
                                      cluster_cell_proportions = c(0.65, 0.3, 0.05),
                                      radii = size,
                                      centre_loc = c(50, 50, 50),
                                      axes_rotation = c(0, 0, 0)
                                    )
                                  ),
                                  plot_image = FALSE),
    iterations = 10,
    check = F
  )
  
})

# Combine results
df <- map2_dfr(results, seq_len(length(ellipsoid_sizes)), ~{
  tibble(
    size = .y,
    median_time = as.numeric(.x$median),      # seconds
    mem_alloc = as.numeric(.x$mem_alloc)      # bytes
  )
})

# Runtime plot
p_runtime <- ggplot(df, aes(size, median_time)) +
  geom_line() +
  geom_point(size = 3) +
  labs(
    x = "Ellipsoid size",
    y = "Median runtime (s)"
  ) +
  theme_minimal()

# Memory plot
p_memory <- ggplot(df, aes(size, mem_alloc / 1024^2)) +
  geom_line() +
  geom_point(size = 3) +
  labs(
    x = "Ellipsoid size",
    y = "Memory allocated (MB)"
  ) +
  theme_minimal()

# Side-by-side
fig <- plot_grid(
  p_runtime,
  p_memory,
  ncol = 2,
  align = "hv"
)

methods::show(fig)

# Testing simulating networks with increasing complexity of networks ----
background_spe <- spaSim3D::simulate_random_background_cells3D(n_cells = 1e6,
                                                               length = 100,
                                                               width = 100,
                                                               height = 100,
                                                               minimum_distance_between_cells = 0,
                                                               plot_image = FALSE)



network_n_edges <- c(10, 15, 20, 25, 30)

results <- lapply(network_n_edges, function(n) {
  
  mark(
    spaSim3D::simulate_clusters3D(background_spe,
                                  cluster_properties_list = list(
                                    C1 = list(
                                      shape = "network",
                                      cluster_cell_types = c("Immune"),
                                      cluster_cell_proportions = c(1),
                                      n_edges = n,
                                      width = 8,
                                      centre_loc = c(50, 50, 50), # Rough centre of network cluster
                                      radius = 40 # Rough radius spanned by the network cluster
                                    )
                                  ),
                                  plot_image = FALSE),
    iterations = 10,
    check = F
  )
  
})

# Combine results
df <- map2_dfr(results, network_n_edges, ~{
  tibble(
    network_n_edges = .y,
    median_time = as.numeric(.x$median),      # seconds
    mem_alloc = as.numeric(.x$mem_alloc)      # bytes
  )
})

# Runtime plot
p_runtime <- ggplot(df, aes(network_n_edges, median_time)) +
  geom_line() +
  geom_point(size = 3) +
  labs(
    x = "Network edges",
    y = "Median runtime (s)"
  ) +
  theme_minimal()

# Memory plot
p_memory <- ggplot(df, aes(network_n_edges, mem_alloc / 1024^2)) +
  geom_line() +
  geom_point(size = 3) +
  labs(
    x = "Network edges",
    y = "Memory allocated (MB)"
  ) +
  theme_minimal()

# Side-by-side
fig <- plot_grid(
  p_runtime,
  p_memory,
  ncol = 2,
  align = "hv"
)

methods::show(fig)
# Testing simulating ellipsoid clusters with increasing number of clusters ----
background_metadata <- spaSim3D::spe_metadata_background_template("random")
background_metadata$background$n_cells <- 1e6
background_metadata$background$length <- 100
background_metadata$background$width <- 100
background_metadata$background$height <- 100
background_metadata$background$minimum_distance_between_cells <- 0
background_metadata$background$cell_types <- "O"
background_metadata$background$cell_proportions <- 1

n_clusters <- seq(5, 25, 5)

results <- lapply(n_clusters, function(n) {
  
  for (i in seq_len(n)) {
    if (i == 1) {
      cluster_metadata <- spaSim3D::spe_metadata_cluster_template("regular",
                                                                  "ellipsoid", 
                                                                  background_metadata
      )
    }
    else {
      cluster_metadata <- spaSim3D::spe_metadata_cluster_template("regular",
                                                                  "ellipsoid", 
                                                                  cluster_metadata
      )
      cluster_metadata[[i + 1]]$radii <- runif(3, 15, 30)
      cluster_metadata[[i + 1]]$centre_loc <- runif(3, 25, 75)
      cluster_metadata[[i + 1]]$axes_rotation <- runif(3, 0, 180)
    }
  }
  
  mark(
    spaSim3D::simulate_spe_metadata3D(cluster_metadata,
                                      plot_image = FALSE),
    iterations = 10,
    check = F
  )
  
})

# Combine results
df <- map2_dfr(results, n_clusters, ~{
  tibble(
    n_clusters = .y,
    median_time = as.numeric(.x$median),      # seconds
    mem_alloc = as.numeric(.x$mem_alloc)      # bytes
  )
})

# Runtime plot
p_runtime <- ggplot(df, aes(n_clusters, median_time)) +
  geom_line() +
  geom_point(size = 3) +
  labs(
    x = "Number of clusters",
    y = "Median runtime (s)"
  ) +
  theme_minimal()

# Memory plot
p_memory <- ggplot(df, aes(n_clusters, mem_alloc / 1024^2)) +
  geom_line() +
  geom_point(size = 3) +
  labs(
    x = "Number of clusters",
    y = "Memory allocated (MB)"
  ) +
  theme_minimal()

# Side-by-side
fig <- plot_grid(
  p_runtime,
  p_memory,
  ncol = 2,
  align = "hv"
)

methods::show(fig)

# Testing simulating empty background vs mixed background vs ME vs RE vs SE vs MN vs RN vs SN ----
background_metadata <- spaSim3D::spe_metadata_background_template("random")
background_metadata$background$n_cells <- 1e6
background_metadata$background$length <- 100
background_metadata$background$width <- 100
background_metadata$background$height <- 100
background_metadata$background$minimum_distance_between_cells <- 0
background_metadata$background$cell_types <- "O"
background_metadata$background$cell_proportions <- 1

tissue_structures <- 
  c(
    "empty_background",
    "mixed_background",
    "mixed_ellipsoid",
    "ringed_ellipsoid",
    "separated_ellipsoid",
    "mixed_network",
    "ringed_network",
    "separated_network"
  )


results <- lapply(tissue_structures, function(tissue_structure) {
  
  if (tissue_structure == "empty_background") {
    cluster_metadata <- background_metadata
  }
  else if (tissue_structure == "mixed_background") {
    cluster_metadata <- background_metadata
    
    cluster_metadata$background$cell_types <- c("A", "B", "O")
    cluster_metadata$background$cell_proportions <- c(0.025, 0.025, 0.95)
    
  }
  else if (tissue_structure == "mixed_ellipsoid") {
    cluster_metadata <- spaSim3D::spe_metadata_cluster_template("regular",
                                                                "ellipsoid", 
                                                                background_metadata)

    cluster_metadata$cluster_1$cluster_cell_types <- c("A", "B")
    cluster_metadata$cluster_1$cluster_cell_proportions <- c(0.5, 0.5)  
    cluster_metadata$cluster_1$radii <- c(20, 30, 40)
    cluster_metadata$cluster_1$centre_loc <- c(50, 50, 50)
    cluster_metadata$cluster_1$axes_rotation <- runif(3, 0, 180)
    
  }
  else if (tissue_structure == "ringed_ellipsoid") {
    cluster_metadata <- spaSim3D::spe_metadata_cluster_template("ring",
                                                                "ellipsoid", 
                                                                background_metadata)
    
    cluster_metadata$cluster_1$cluster_cell_types <- "A"
    cluster_metadata$cluster_1$cluster_cell_proportions <- 1
    cluster_metadata$cluster_1$radii <- c(25, 50, 75)
    cluster_metadata$cluster_1$centre_loc <- c(50, 50, 50)
    cluster_metadata$cluster_1$axes_rotation <- runif(3, 0, 180)
    cluster_metadata$cluster_1$ring_cell_types <- "B"
    cluster_metadata$cluster_1$ring_cell_proportions <- 1
    cluster_metadata$cluster_1$ring_width <- 10
     
  }
  else if (tissue_structure == "separated_ellipsoid") {
    cluster_metadata <- spaSim3D::spe_metadata_cluster_template("regular",
                                                                "ellipsoid", 
                                                                background_metadata)
    
    cluster_metadata <- spaSim3D::spe_metadata_cluster_template("regular",
                                                                "ellipsoid", 
                                                                cluster_metadata)
    
    
  }
  else if (tissue_structure == "mixed_network") {
    cluster_metadata <- spaSim3D::spe_metadata_cluster_template("regular",
                                                                "network", 
                                                                background_metadata)
    
    cluster_metadata$cluster_1$cluster_cell_types <- c("A", "B")
    cluster_metadata$cluster_1$cluster_cell_proportions <- c(0.5, 0.5)  
    cluster_metadata$cluster_1$radius <- 50
    cluster_metadata$cluster_1$centre_loc <- c(50, 50, 50)
    cluster_metadata$cluster_1$width <- 5
  }
  else if (tissue_structure == "ringed_network") {
    cluster_metadata <- spaSim3D::spe_metadata_cluster_template("ring",
                                                                "network", 
                                                                background_metadata)
    
    cluster_metadata$cluster_1$cluster_cell_types <- "A"
    cluster_metadata$cluster_1$cluster_cell_proportions <- 1
    cluster_metadata$cluster_1$radius <- 50
    cluster_metadata$cluster_1$centre_loc <- c(50, 50, 50)
    cluster_metadata$cluster_1$width <- 5
    cluster_metadata$cluster_1$ring_cell_types <- "B"
    cluster_metadata$cluster_1$ring_cell_proportions <- 1
    cluster_metadata$cluster_1$ring_width <- 1
    
  }
  else if (tissue_structure == "separated_network") {
    cluster_metadata <- spaSim3D::spe_metadata_cluster_template("regular",
                                                                "network", 
                                                                background_metadata)
    
    cluster_metadata <- spaSim3D::spe_metadata_cluster_template("regular",
                                                                "network", 
                                                                cluster_metadata)
  }
  else {
    stop("tissue structure not found")
  }
  
  mark(
    spaSim3D::simulate_spe_metadata3D(cluster_metadata,
                                      plot_image = FALSE),
    iterations = 10,
    check = F
  )
  
})

# Combine results
df <- map2_dfr(results, tissue_structures, ~{
  tibble(
    tissue_structures = .y,
    median_time = as.numeric(.x$median),      # seconds
    mem_alloc = as.numeric(.x$mem_alloc)      # bytes
  )
})

# Runtime plot
p_runtime <- ggplot(df, aes(tissue_structures, median_time)) +
  geom_point(size = 3) +
  labs(
    x = "Tissue structure",
    y = "Median runtime (s)"
  ) +
  theme_minimal()

# Memory plot
p_memory <- ggplot(df, aes(tissue_structures, mem_alloc / 1024^2)) +
  geom_line() +
  geom_point(size = 3) +
  labs(
    x = "Tissue structure",
    y = "Memory allocated (MB)"
  ) +
  theme_minimal()

# Side-by-side
fig <- plot_grid(
  p_runtime,
  p_memory,
  ncol = 2,
  align = "hv"
)

methods::show(fig)


