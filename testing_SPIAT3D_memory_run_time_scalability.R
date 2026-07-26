### SPIAT3D testing
## parameters:
# - number of cells in 3D data
# - number of radii for gradient-based metrics
# - number of splits for grid-based metrics
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
library(spaSim3D)
library(SPIAT3D)


# Testing calculate_minimum_distances_between_cell_types3D with increasing number of cells ----
n_cells <- c(1e3, 1e4, 1e5, 1e6, 1e7)

background_metadata <- spaSim3D::spe_metadata_background_template("random")
background_metadata$background$length <- 100
background_metadata$background$width <- 100
background_metadata$background$height <- 100
background_metadata$background$minimum_distance_between_cells <- 0
background_metadata$background$cell_types <- c("A", "B", "O")
background_metadata$background$cell_proportions <- c(0.4, 0.4, 0.2)

results <- lapply(n_cells, function(n) {
  background_metadata$background$n_cells <- n
  
  background_spe <- spaSim3D::simulate_spe_metadata3D(background_metadata,
                                                      plot_image = FALSE)
  
  mark(
    SPIAT3D::calculate_minimum_distances_between_cell_types3D(background_spe,
                                                              cell_types_of_interest = NULL,
                                                              feature_colname = "Cell.Type",
                                                              show_summary = FALSE,
                                                              plot_image = FALSE),
    iterations = 10,
    check = F
  )
  
})

plot_results_n_cells(results, n_cells)

# Testing calculate_neighbourhood_counts_gradient3D with increasing number of cells ----
n_cells <- c(1e3, 1e4, 1e5, 1e6, 1e7)

background_metadata <- spaSim3D::spe_metadata_background_template("random")
background_metadata$background$length <- 100
background_metadata$background$width <- 100
background_metadata$background$height <- 100
background_metadata$background$minimum_distance_between_cells <- 0
background_metadata$background$cell_types <- c("A", "B", "O")
background_metadata$background$cell_proportions <- c(0.4, 0.4, 0.2)

results <- lapply(n_cells, function(n) {
  background_metadata$background$n_cells <- n
  
  background_spe <- spaSim3D::simulate_spe_metadata3D(background_metadata,
                                                      plot_image = FALSE)
  
  mark(
    SPIAT3D::calculate_neighbourhood_counts_gradient3D(background_spe,
                                                       reference_cell_type = "A",
                                                       target_cell_types = c("A", "B"),
                                                       radii = seq(10, 100, 10),
                                                       feature_colname = "Cell.Type",
                                                       plot_image = FALSE),
    iterations = 10,
    check = F
  )
  
})

plot_results_n_cells(results, n_cells)

# Testing calculate_co_occurrence_gradient3D with increasing number of cells ----
n_cells <- c(1e3, 1e4, 1e5, 1e6, 1e7)

background_metadata <- spaSim3D::spe_metadata_background_template("random")
background_metadata$background$length <- 100
background_metadata$background$width <- 100
background_metadata$background$height <- 100
background_metadata$background$minimum_distance_between_cells <- 0
background_metadata$background$cell_types <- c("A", "B", "O")
background_metadata$background$cell_proportions <- c(0.4, 0.4, 0.2)

results <- lapply(n_cells, function(n) {
  background_metadata$background$n_cells <- n
  
  background_spe <- spaSim3D::simulate_spe_metadata3D(background_metadata,
                                                      plot_image = FALSE)
  
  mark(
    SPIAT3D::calculate_co_occurrence3D(background_spe,
                                       reference_cell_type = "A",
                                       target_cell_types = c("A", "B"),
                                       radii = seq(10, 100, 10),
                                       feature_colname = "Cell.Type",
                                       plot_image = FALSE),
    iterations = 10,
    check = F
  )
  
})

plot_results_n_cells(results, n_cells)


# Testing calculate_cell_proportion_grid_metrics3D with increasing number of cells ----
n_cells <- c(1e3, 1e4, 1e5, 1e6, 1e7)

background_metadata <- spaSim3D::spe_metadata_background_template("random")
background_metadata$background$length <- 100
background_metadata$background$width <- 100
background_metadata$background$height <- 100
background_metadata$background$minimum_distance_between_cells <- 0
background_metadata$background$cell_types <- c("A", "B", "O")
background_metadata$background$cell_proportions <- c(0.4, 0.4, 0.2)

results <- lapply(n_cells, function(n) {
  background_metadata$background$n_cells <- n
  
  background_spe <- spaSim3D::simulate_spe_metadata3D(background_metadata,
                                                      plot_image = FALSE)
  
  mark(
    SPIAT3D::calculate_cell_proportion_grid_metrics3D(background_spe,
                                                      n_splits = 10,
                                                      reference_cell_types = "A",
                                                      target_cell_types = "B",
                                                      feature_colname = "Cell.Type",
                                                      plot_image = FALSE),
    iterations = 10,
    check = F
  )
  
})

plot_results_n_cells(results, n_cells)


# Testing calculate_prevalence_gradient3D with increasing number of cells ----
n_cells <- c(1e3, 1e4, 1e5, 1e6, 1e7)

background_metadata <- spaSim3D::spe_metadata_background_template("random")
background_metadata$background$length <- 100
background_metadata$background$width <- 100
background_metadata$background$height <- 100
background_metadata$background$minimum_distance_between_cells <- 0
background_metadata$background$cell_types <- c("A", "B", "O")
background_metadata$background$cell_proportions <- c(0.4, 0.4, 0.2)

results <- lapply(n_cells, function(n) {
  background_metadata$background$n_cells <- n
  
  background_spe <- spaSim3D::simulate_spe_metadata3D(background_metadata,
                                                      plot_image = FALSE)
  
  proportion_grid_metrics <- SPIAT3D::calculate_cell_proportion_grid_metrics3D(background_spe,
                                                                               n_splits = 10,
                                                                               reference_cell_types = "A",
                                                                               target_cell_types = "B",
                                                                               feature_colname = "Cell.Type",
                                                                               plot_image = FALSE)
  
  mark(
    SPIAT3D::calculate_prevalence_gradient3D(grid_metrics = proportion_grid_metrics,
                                             metric_colname = "proportion",
                                             show_AUC = FALSE,
                                             plot_image = FALSE),
    iterations = 10,
    check = F
  )
  
})

plot_results_n_cells(results, n_cells)



# Testing calculate_spatial_autocorrelation3D with increasing number of cells ----
n_cells <- c(1e3, 1e4, 1e5, 1e6, 1e7)

background_metadata <- spaSim3D::spe_metadata_background_template("random")
background_metadata$background$length <- 100
background_metadata$background$width <- 100
background_metadata$background$height <- 100
background_metadata$background$minimum_distance_between_cells <- 0
background_metadata$background$cell_types <- c("A", "B", "O")
background_metadata$background$cell_proportions <- c(0.4, 0.4, 0.2)

results <- lapply(n_cells, function(n) {
  background_metadata$background$n_cells <- n
  
  background_spe <- spaSim3D::simulate_spe_metadata3D(background_metadata,
                                                      plot_image = FALSE)
  
  proportion_grid_metrics <- SPIAT3D::calculate_cell_proportion_grid_metrics3D(background_spe,
                                                                               n_splits = 10,
                                                                               reference_cell_types = "A",
                                                                               target_cell_types = "B",
                                                                               feature_colname = "Cell.Type",
                                                                               plot_image = FALSE)
  
  mark(
    SPIAT3D::calculate_spatial_autocorrelation(grid_metrics = proportion_grid_metrics,
                                               metric_colname = "proportion",
                                               weight_method = "queen"),
    iterations = 10,
    check = F
  )
  
})

plot_results_n_cells(results, n_cells)



# Testing calculate_neighbourhood_counts_gradient3D with increasing number of radii values ----
background_metadata <- spaSim3D::spe_metadata_background_template("random")
background_metadata$background$n_cells <- 1e6
background_metadata$background$length <- 100
background_metadata$background$width <- 100
background_metadata$background$height <- 100
background_metadata$background$minimum_distance_between_cells <- 0
background_metadata$background$cell_types <- c("A", "B", "O")
background_metadata$background$cell_proportions <- c(0.4, 0.4, 0.2)

n_radii <- seq(5, 25, 5)

results <- lapply(n_radii, function(n) {
  
  background_spe <- spaSim3D::simulate_spe_metadata3D(background_metadata,
                                                      plot_image = FALSE)
  
  mark(
    SPIAT3D::calculate_neighbourhood_counts_gradient3D(background_spe,
                                                       reference_cell_type = "A",
                                                       target_cell_types = c("A", "B"),
                                                       radii = seq(10, 100, length.out = n),
                                                       feature_colname = "Cell.Type",
                                                       plot_image = FALSE),
    iterations = 10,
    check = F
  )
  
})

# Combine results
df <- map2_dfr(results, n_radii, ~{
  tibble(
    n_radii = .y,
    median_time = as.numeric(.x$median),      # seconds
    mem_alloc = as.numeric(.x$mem_alloc)      # bytes
  )
})

# Runtime plot
p_runtime <- ggplot(df, aes(n_radii, median_time)) +
  geom_line() +
  geom_point(size = 3) +
  labs(
    x = "Number of radii",
    y = "Median runtime (s)"
  ) +
  theme_minimal()

# Memory plot
p_memory <- ggplot(df, aes(n_radii, mem_alloc / 1024^2)) +
  geom_line() +
  geom_point(size = 3) +
  labs(
    x = "Number of radii",
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


# Testing calculate_cell_proportion_grid_metrics3D with increasing number of splits ----
background_metadata <- spaSim3D::spe_metadata_background_template("random")
background_metadata$background$n_cells <- 1e6
background_metadata$background$length <- 100
background_metadata$background$width <- 100
background_metadata$background$height <- 100
background_metadata$background$minimum_distance_between_cells <- 0
background_metadata$background$cell_types <- c("A", "B", "O")
background_metadata$background$cell_proportions <- c(0.4, 0.4, 0.2)

n_splits <- seq(5, 25, 5)

results <- lapply(n_splits, function(n) {
  
  background_spe <- spaSim3D::simulate_spe_metadata3D(background_metadata,
                                                      plot_image = FALSE)
  
  mark(
    SPIAT3D::calculate_cell_proportion_grid_metrics3D(background_spe,
                                                      n_splits = n,
                                                      reference_cell_types = "A",
                                                      target_cell_types = "B",
                                                      feature_colname = "Cell.Type",
                                                      plot_image = FALSE),
    iterations = 10,
    check = F
  )
  
})

# Combine results
df <- map2_dfr(results, n_splits, ~{
  tibble(
    n_splits = .y,
    median_time = as.numeric(.x$median),      # seconds
    mem_alloc = as.numeric(.x$mem_alloc)      # bytes
  )
})

# Runtime plot
p_runtime <- ggplot(df, aes(n_splits, median_time)) +
  geom_line() +
  geom_point(size = 3) +
  labs(
    x = "Number of splits",
    y = "Median runtime (s)"
  ) +
  theme_minimal()

# Memory plot
p_memory <- ggplot(df, aes(n_splits, mem_alloc / 1024^2)) +
  geom_line() +
  geom_point(size = 3) +
  labs(
    x = "Number of splits",
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
# Testing calculate_minimum_distances_between_cell_types3D with increasing number of cell type comparisons ----
background_metadata <- spaSim3D::spe_metadata_background_template("random")
background_metadata$background$n_cells <- 1e6
background_metadata$background$length <- 100
background_metadata$background$width <- 100
background_metadata$background$height <- 100
background_metadata$background$minimum_distance_between_cells <- 0

n_cell_types <- seq(2, 10, 2)

results <- lapply(n_cell_types, function(n) {
  
  background_metadata$background$cell_types <- LETTERS[(seq(n))]
  background_metadata$background$cell_proportions <- rexp(n)
  background_metadata$background$cell_proportions <- round(background_metadata$background$cell_proportions / sum(background_metadata$background$cell_proportions), 5)
  
  background_spe <- spaSim3D::simulate_spe_metadata3D(background_metadata,
                                                      plot_image = FALSE)
  
  mark(
    SPIAT3D::calculate_minimum_distances_between_cell_types3D(background_spe,
                                                              cell_types_of_interest = NULL,
                                                              feature_colname = "Cell.Type",
                                                              show_summary = FALSE,
                                                              plot_image = FALSE),
    iterations = 10,
    check = F
  )
  
})

# Combine results
df <- map2_dfr(results, n_cell_types, ~{
  tibble(
    n_cell_types = .y,
    median_time = as.numeric(.x$median),      # seconds
    mem_alloc = as.numeric(.x$mem_alloc)      # bytes
  )
})

# Runtime plot
p_runtime <- ggplot(df, aes(n_cell_types, median_time)) +
  geom_line() +
  geom_point(size = 3) +
  labs(
    x = "Number of cell types",
    y = "Median runtime (s)"
  ) +
  theme_minimal()

# Memory plot
p_memory <- ggplot(df, aes(n_cell_types, mem_alloc / 1024^2)) +
  geom_line() +
  geom_point(size = 3) +
  labs(
    x = "Number of cell types",
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

# Testing calculate_neighbourhood_counts_gradient3D with increasing number of cell type comparisons ----
background_metadata <- spaSim3D::spe_metadata_background_template("random")
background_metadata$background$n_cells <- 1e6
background_metadata$background$length <- 100
background_metadata$background$width <- 100
background_metadata$background$height <- 100
background_metadata$background$minimum_distance_between_cells <- 0

n_cell_types <- seq(2, 10, 2)

results <- lapply(n_cell_types, function(n) {
  
  background_metadata$background$cell_types <- LETTERS[(seq(n))]
  background_metadata$background$cell_proportions <- rexp(n)
  background_metadata$background$cell_proportions <- round(background_metadata$background$cell_proportions / sum(background_metadata$background$cell_proportions), 5)
  
  background_spe <- spaSim3D::simulate_spe_metadata3D(background_metadata,
                                                      plot_image = FALSE)
  
  mark(
    SPIAT3D::calculate_neighbourhood_counts_gradient3D(background_spe,
                                                       reference_cell_type = "A",
                                                       target_cell_types = LETTERS[(seq(n))],
                                                       radii = seq(10, 100, 10),
                                                       feature_colname = "Cell.Type",
                                                       plot_image = FALSE),
    iterations = 10,
    check = F
  )
  
})

# Combine results
df <- map2_dfr(results, n_cell_types, ~{
  tibble(
    n_cell_types = .y,
    median_time = as.numeric(.x$median),      # seconds
    mem_alloc = as.numeric(.x$mem_alloc)      # bytes
  )
})

# Runtime plot
p_runtime <- ggplot(df, aes(n_cell_types, median_time)) +
  geom_line() +
  geom_point(size = 3) +
  labs(
    x = "Number of cell types",
    y = "Median runtime (s)"
  ) +
  theme_minimal()

# Memory plot
p_memory <- ggplot(df, aes(n_cell_types, mem_alloc / 1024^2)) +
  geom_line() +
  geom_point(size = 3) +
  labs(
    x = "Number of cell types",
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

