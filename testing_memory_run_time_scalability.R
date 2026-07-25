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

# Libraries ----
library(bench)
library(dplyr)
library(purrr)
library(ggplot2)


# Testing simulate_random_background_cells function with increasing number of cells ----
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


# Combine results
df <- map2_dfr(results, n_cells, ~{
  tibble(
    size = .y,
    median_time = as.numeric(.x$median),      # seconds
    mem_alloc = as.numeric(.x$mem_alloc)      # bytes
  )
})

df

# Plot runtime
ggplot(df, aes(size, median_time)) +
  geom_line() +
  geom_point(size = 3) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Input size",
    y = "Median runtime (s)"
  ) +
  theme_minimal()


# Plot memory usage
ggplot(df, aes(size, mem_alloc / 1024^2)) +
  geom_line() +
  geom_point(size = 3) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Input size",
    y = "Memory allocated (MB)"
  ) +
  theme_minimal()

# Testing simulate_ordered_background_cells function with increasing number of cells ----
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


# Combine results
df <- map2_dfr(results, n_cells, ~{
  tibble(
    size = .y,
    median_time = as.numeric(.x$median),      # seconds
    mem_alloc = as.numeric(.x$mem_alloc)      # bytes
  )
})

df

# Plot runtime
ggplot(df, aes(size, median_time)) +
  geom_line() +
  geom_point(size = 3) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Input size",
    y = "Median runtime (s)"
  ) +
  theme_minimal()


# Plot memory usage
ggplot(df, aes(size, mem_alloc / 1024^2)) +
  geom_line() +
  geom_point(size = 3) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Input size",
    y = "Memory allocated (MB)"
  ) +
  theme_minimal()
