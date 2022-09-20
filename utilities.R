library(pracma)

detach.all.packages <- function() {
  basic.packages <- c(
    "package:stats",
    "package:graphics",
    "package:grDevices",
    "package:utils",
    "package:datasets",
    "package:methods",
    "package:base"
  )
  package.list <- search()[ifelse(unlist(gregexpr("package:",search())) == 1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list) > 0)  for (package in package.list) detach(package, character.only=TRUE)
}

utils.log <- function(log.str) {
  
  print(paste(format(Sys.time(), "[%X]"), log.str))
}

is.valid.value <- function(value)
{
  return(!is.na(value) & value != "")
}

get.cross.manifold.similarity.matrix.score <- function(W)
{
  stopifnot(ncol(W) == nrow(W))
  # Measure the weights matrix score - same samples from both manifolds should be
  # close to each other, i.e. if X_i and Y_j represents the same sample, W[i,j] should be high
  X_rankings <- sapply(1:nrow(W), function(i) which(order(W[i,], decreasing = TRUE) == i))
  Y_rankings <- sapply(1:ncol(W), function(j) which(order(W[,j], decreasing = TRUE) == j))
  
  # First set the rankings to be in range [0,#samples-1]
  X_rankings <- X_rankings - 1 
  Y_rankings <- Y_rankings - 1
  
  # Then normalize to get a score between 0 to 1
  # For the sake of this score, we assume that we have the same number of samples from both datasets
  number_of_samples <- length(X_rankings) 
  score <- (sum(X_rankings) + sum(Y_rankings)) / (2 * number_of_samples * (number_of_samples - 1))
  return(score)
}

# returns the similarity matrix which describes the similarity between the columns of X and the columns of Y
get.inter.similarity.matrix <- function(X,Y)
{
  distances = distmat(t(X),t(Y))
  
  # An RBF-Kernel with sigma = mean(distances)/sqrt(2)
  W <- exp(-(distances^2) / (mean(distances)^2))
  return(W)
}

get.multi.subtypes.name <- function(subtypes)
{
  paste(subtypes, collapse = "_")
}

get.highly.variable.features <- function(data, number_of_features_to_keep) {
  
  var_data <- apply(data, 1, var)
  
  top_variance_features_to_consider <- min(number_of_features_to_keep, dim(data)[1])
  indices_of_features_with_top_variance <-
    order(var_data, decreasing = TRUE)[1:top_variance_features_to_consider]
  
  # Also remove zero variance genes if left by intersecting with non-zero variance rows
  indices_of_features_with_non_zero_variance <- which(var_data > 0)
  indices_of_features_with_top_variance <-
    intersect(
      indices_of_features_with_top_variance,
      indices_of_features_with_non_zero_variance
    )
  
  return(indices_of_features_with_top_variance)
}

select.highly.variable.features <- function(data, number_of_features_to_keep) {
  
  indices_of_features_with_top_variance <- get.highly.variable.features(data, number_of_features_to_keep)
  
  data_with_top_variance_features <- data[indices_of_features_with_top_variance, ]
  
  return(data_with_top_variance_features)
}

scale.data <- function(data) {
  # scale scales and centers the columns, so we transpose the data matrix
  data <- t(scale(t(data), center = TRUE, scale = TRUE))
  return(data.matrix(data))
}
