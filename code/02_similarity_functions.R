
# Szymkiewicz-Simpson (Overlap Coefficient) similarity --------------------


overlap_coef_bin <- function(x) {
  # Calculate the pairwise sums of non-zero matrix elements to find the number
  # of intersecting elements between each input column
  bin_intersect_mat <- t(x) %*% x
  
  # Calculate the input column sums to find the individual size of each set
  col_sum <- apply(x, 2, function(xx)
    sum(xx != 0))
  
  # Find the smaller of each pair of sets by taking the matrix outer minimum
  min_set_size_mat <- outer(col_sum, col_sum, FUN = pmin)
  
  # Szymkiewicz-Simpson overlap coefficient is the pairwise intersection of two
  # sets divided by the size of the smaller set
  res <- bin_intersect_mat / min_set_size_mat
  
  # Set diagonal to identity
  diag(res) <- 1L
  
  # Assign input column names to rows and columns of return matrix
  dimnames(res) <- list(colnames(x), colnames(x))
  
  return(res)
}

# Sorenson-Dice similarity ------------------------------------------------


soren_dice_sim_bin <- function(x) {
  # Calculate the pairwise sums of non-zero matrix elements to find the number
  # of intersecting elements between each input column
  bin_intersect_mat <- t(x) %*% x
  
  # Calculate the input column sums, to find the individual size of each set
  col_sum <- apply(x, 2, function(xx)
    sum(xx != 0))
  
  # Calculate the matrix outer sums for pairwise sum of set sizes
  set_size_sum_mat <- outer(col_sum, col_sum, FUN = "+")
  
  # Sorenson-Dice index is twice the size of the intersection divided by the
  # sum of the size for each set
  res <- (2 * bin_intersect_mat) / set_size_sum_mat
  
  # Set diagonal to identity
  diag(res) <- 1L
  
  # Assign input column names to rows and columns of return matrix
  dimnames(res) <- list(colnames(x), colnames(x))
  return(res)
  
}

# Jaccard similarity ------------------------------------------------------


jaccard_sim_bin <- function(x) {
  # Calculate the pairwise sums of non-zero matrix elements to find the number
  # of intersecting elements between each input column
  bin_intersect_mat <- t(x) %*% x
  
  # Calculate the input column sums to find the individual size of each set
  x_col_sum <- apply(x, 2, function(xx)
    sum(xx != 0))
  
  # Calculate the matrix outer sums for pairwise sum of set sizes
  set_size_sum_mat <- outer(x_col_sum, x_col_sum, FUN = "+")
  
  # Jaccard index is intersection of set sizes over the size of the union of
  # sets
  res <- bin_intersect_mat / (set_size_sum_mat - bin_intersect_mat)
  
  # Set diagonal to identity
  diag(res) <- 1L
  
  # Assign input column names to rows and columns of return matrix
  dimnames(res) <- list(colnames(x), colnames(x))
  
  return(res)
}

# TOM Adjacency -----------------------------------------------------------


topological_overlap <- function(adj_mat) {
  l_mat <- adj_mat %*% t(adj_mat)
  
  k_row <- rowSums(adj_mat)
  k_col <- colSums(adj_mat)
  k_min <- outer(k_row, k_col, FUN = pmin)
  
  tom <- (l_mat + adj_mat) / (k_min + 1 - adj_mat)
  
  diag(tom) <- 1
  
  dimnames(tom) <- list(colnames(adj_mat), colnames(adj_mat))
  
  return(tom)
}

# Hard threshold ----------------------------------------------------------

signum_adj <- function(x, tau = 0.5) {
  ifelse(x < tau, 0, 1)
}

# Soft thresholds ---------------------------------------------------------


## Sigmoid adjacency ------------------------------------------------------

sigmoid_adj <- function(x, alpha, mu) {
  1.0 / (1.0 + exp(-alpha * (x - mu)))
}

## Power adjacency --------------------------------------------------------

power_adj <- function(x, beta = 1) {
  abs(x) ^ beta
}

