### spaSim3D testing
## parameters:
# - number of cells in simulation
# - size of clusters
# - shape of clusters
# - number of clusters
# - complexity of networks

### SPIAT3D testing
## parameters:
# - number of cells in 3D data
# - number of radii for gradient-based metrics
# - number of cell comparisons

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

# Testing simulate_clusters3D function (ellipsoids) with increasing number of cells in the background ----
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



# Testing simulate_clusters3D function (network) with increasing number of cells in the background ----
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

# Testing simulate_clusters3D function (ellipsoids) with increasing size of clusters ----
background_spe <- spaSim3D::simulate_random_background_cells3D(n_cells = 1e5,
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

# Testing simulate_clusters3D function (networks) with increasing complexity of networks ----
background_spe <- spaSim3D::simulate_random_background_cells3D(n_cells = 1e5,
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
# Testing simulate_clusters3D function with increasing number of clusters ----