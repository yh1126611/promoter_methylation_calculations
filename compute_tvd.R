compute_tvd_pvalues <- function(wg_file, element_file, n_perm = 1000) {
  library(data.table)
  compute_tvd <- function(x, y, breaks = 99) { # Nested function
    if (length(y) == 0) { warning("Empty y vector - returning NA"); return(NA) }
    if (anyNA(y)) { warning("NA values in y - removing them"); y <- y[!is.na(y)] }
    y_range <- range(y); x_range <- range(x)
    combined_range <- if (y_range[1] < x_range[1] || y_range[2] > x_range[2]) {
      warning("Y values outside X range - adjusting breaks"); range(c(x, y))
    } else x_range
    breaks_vec <- seq(combined_range[1], combined_range[2], length.out = breaks + 1)
    hist_x <- hist(x, breaks = breaks_vec, plot = FALSE)
    hist_y <- hist(y, breaks = breaks_vec, plot = FALSE)
    p_x <- hist_x$counts / length(x); p_y <- hist_y$counts / length(y)
    0.5 * sum(abs(p_x - p_y))
  }
  wg_distrbtn <- fread(wg_file)[, V4]; mp <- fread(element_file) # Load population & sample
  p_values <- data.frame(spcdist = integer(), p_value = numeric()) # Create empty data frame “p_values”
  for (spcdist in -10000:10000) {
    spcdist_distrbtn <- mp[V2 == spcdist, V3]
    if (length(spcdist_distrbtn) == 0) { warning(paste("No data for spcdist", spcdist)); next }
    obs_tvd <- compute_tvd(wg_distrbtn, spcdist_distrbtn); n <- length(spcdist_distrbtn)
    perm_tvd <- replicate(n_perm, {
      shuffled <- sample(wg_distrbtn, size = 2 * n)
      compute_tvd(shuffled[1:n], shuffled[(n + 1):(2 * n)])
    })
    p_values <- rbind(p_values, data.frame(spcdist = spcdist, p_value = mean(perm_tvd >= obs_tvd)))
  }
  return(p_values)
}
