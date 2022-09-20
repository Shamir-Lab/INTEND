library(Matrix)
library(pracma)
library(igraph)
library(mstknnclust)
library(parallel)
library(Rfast)
library(dnet)

source('utilities.R')
source('join_manifolds_matrices_computation.R')

get_projection_cross_manifold_similarity_scores_for_d_range <- function(X_projected, Y_projected, d_range)
{
  scores <- mcmapply(function(d) {
    return(get_projection_cross_manifold_similarity_score(X_projected[1:d,], Y_projected[1:d,]))
  },
  d_range)

  names(scores) <- d_range
  
  return(scores)
}

get_projection_cross_manifold_similarity_score_per_sample <- function(X_projected, Y_projected)
{
  # Measure if the alignment succeeded in putting the same samples from both manifolds
  # close to each other in the projection space
  distances_in_projection_space <- distmat(t(X_projected), t(Y_projected))
  
  X_rankings <- sapply(1:nrow(distances_in_projection_space), function(i) which(order(distances_in_projection_space[i,]) == i))
  Y_rankings <- sapply(1:ncol(distances_in_projection_space), function(j) which(order(distances_in_projection_space[,j]) == j))
  
  # First set the rankings to be in range [0,#samples-1]
  X_rankings <- X_rankings - 1
  Y_rankings <- Y_rankings - 1
  
  # Then normalize to get a score between 0 to 1
  # For the sake of this score, we assume that we have the same number of samples from both datasets
  number_of_samples <- length(X_rankings) 
  
  score.per.sample <- c(X_rankings, Y_rankings)/(number_of_samples - 1)
  return(score.per.sample)
}

get_projection_cross_manifold_similarity_score <- function(X_projected, Y_projected)
{
  return(mean(get_projection_cross_manifold_similarity_score_per_sample(X_projected, Y_projected)))
}

get_projection_manifold_similarity_preservation_score <- function(data, projected_data)
{
  original_similarity_matrix <- get_similarity_matrix(data)
  projected_space_similarity_matrix <- get_similarity_matrix(projected_data)
  
  # A vector with entry for each sample which describes the preservation of the similarities of this sample to
  # other samples after the projection
  similarity_persevation_scores <- 
    sapply(1:ncol(original_similarity_matrix),
           function(i) cor(original_similarity_matrix[,i], projected_space_similarity_matrix[, i]))
  
  # Each entry in similarity_persevation_scores is in the range [-1,1], where 1 is the optimal score as the
  # new similarities in the projection are correlated well with the original similarities. -1 is the worse
  # score as the new similarities are anti-correlated with the original similarities. The expected score for a
  # random projection is 0, with no correlation between the original similarities and the projection ones.
  mean_score <- mean(similarity_persevation_scores)
  return(mean_score)
}

# returns the similarity matrix which describes the similarity between the columns of X
get_similarity_matrix <- function(X)
{
  return(get.inter.similarity.matrix(X,X))
}

my.generate.knn <- function(edges.complete.graph, suggested.k) 
{
  grafo_knn = list()
  arista_vecinos_unidas = list()
  grafo_knn_conectado = vector()
  nodos.object_i <- unique(edges.complete.graph$object_i)
  nodos.object_j <- unique(edges.complete.graph$object_j)
  nodos <- c(nodos.object_i, nodos.object_j)
  nodos <- unique(nodos)
  n <- length(nodos)
  k = 1
  while (k <= (n - 1)) {
    aristas.ordenadas <- edges.complete.graph[order(edges.complete.graph$d_ij), 
    ]
    aristas.ordenadas$object_i <- as.factor(aristas.ordenadas$object_i)
    aristas.ordenadas$object_j <- as.factor(aristas.ordenadas$object_j)
    vecinos.nodos.object_i <- do.call(rbind, lapply(split(aristas.ordenadas, 
                                                          aristas.ordenadas$object_i), function(x) {
                                                            return(x[1:k, ])
                                                          }))
    vecinos.nodos.object_i <- stats::na.omit(vecinos.nodos.object_i)
    vecinos.nodos.object_j <- do.call(rbind, lapply(split(aristas.ordenadas, 
                                                          aristas.ordenadas$object_j), function(x) {
                                                            return(x[1:k, ])
                                                          }))
    vecinos.nodos.object_j <- stats::na.omit(vecinos.nodos.object_j)
    vecinos.nodos.object_j <- vecinos.nodos.object_j[, c(2, 
                                                         1, 3)]
    colnames(vecinos.nodos.object_j) <- c("object_i", "object_j", 
                                          "d_ij")
    vecinos.nodos.ambos <- rbind(vecinos.nodos.object_i, 
                                 vecinos.nodos.object_j)
    ambos.ordenados = vecinos.nodos.ambos[order(vecinos.nodos.ambos$d_ij), 
    ]
    vecinos.final <- do.call(rbind, lapply(split(ambos.ordenados, 
                                                 ambos.ordenados$object_i), function(x) {
                                                   return(x[1:k, ])
                                                 }))
    vecinos.final <- stats::na.omit(vecinos.final)
    arista_vecinos_unidas[[k]] = vecinos.final[!duplicated(vecinos.final), 
    ]
    grafo_knn[[k]] = igraph::graph.data.frame(d = arista_vecinos_unidas[[k]][, 
                                                                             1:2], directed = FALSE)
    grafo_knn[[k]] = igraph::simplify(grafo_knn[[k]], remove.loops = TRUE, 
                                      remove.multiple = FALSE)
    grafo_knn_conectado[k] = igraph::is.connected(grafo_knn[[k]])
    if (grafo_knn_conectado[k] == TRUE && (missing(suggested.k) || k == suggested.k)) {
      k = n
    }
    else {
      k = k + 1
    }
  }
  evaluacion_k = which(grafo_knn_conectado == TRUE)
  if (length(evaluacion_k) > 0) {
    k_conectado = min(evaluacion_k)
  }
  else {
    cat("\n In any k the graph knn can be connected. It will use as k the log(n). \n")
    k_conectado = n
  }
  k_natural = floor(log(n))
  if (missing(suggested.k)) {
    valor_k = min(k_natural, k_conectado)
  }
  else {
    valor_k = suggested.k
  }
  if (valor_k == 0) {
    valor_k = 1
  }
  return(list(edges.knn.graph = arista_vecinos_unidas[[valor_k]], 
              k = valor_k, knn.graph = grafo_knn[[valor_k]]))
}

get.knn.graph.adjacency.matrix <- function(X,Y) {
  
  distances <- distmat(t(X), t(Y))
  complete.graph <- generate.complete.graph(nodes.list = 1:ncol(distances), distance.matrix = distances)
  knn <- my.generate.knn(edges.complete.graph = complete.graph)
  edges <- knn$edges.knn.graph
  knn.distances <- matrix(data = 0, nrow = nrow(distances), ncol = ncol(distances))
  
  for (edge in 1:nrow(edges)) {
    
    first_node <- as.integer(as.character(edges$object_i[edge]))
    second_node <- as.integer(as.character(edges$object_j[edge]))
    edge_weight <- edges$d_ij[edge]
    
    knn.distances[first_node, second_node] <- edge_weight
    
    # Our graph is not directed, and so we want the knn_distances matrix to by symmetric
    knn.distances[second_node, first_node] <- edge_weight
  }
  
  # floyd function expect NA for missing edge between 2 nodes (and not zero)
  knn.distances[knn.distances == 0] <- NA
  
  return(knn.distances)
}

get.inter.similarity.matrix.based.on.geodesic.ditance.in.knn.graph <- function(X, Y) {
  
  knn.distances <- get.knn.graph.adjacency.matrix(X, Y)
  
  shortest.paths <- floyd(knn.distances)
  
  W <- exp(-(shortest.paths^2) / (mean(shortest.paths)^2))
  return(W)
}

get.similarity.matrix.based.on.geodesic.ditance.in.knn.graph <- function(X) {
  
  return(get.inter.similarity.matrix.based.on.geodesic.ditance.in.knn.graph(X,X))
}

get.inter.similarity.matrix.based.on.RWR.in.knn.graph <- function(X, Y) {
  
  distances <- distmat(t(X), t(Y))
  complete.graph <- generate.complete.graph(nodes.list = 1:ncol(distances), distance.matrix = distances)
  knn <- my.generate.knn(edges.complete.graph = complete.graph)

  edges.data.frame <- knn$edges.knn.graph
  colnames(edges.data.frame) <- c("from", "to", "weight")
  knn.weighted.graph <- graph_from_data_frame(edges.data.frame, directed = FALSE, vertices = NULL)
  
  affinity.matrix <- dRWR(g = knn.weighted.graph, restart = 0.999)
  affinity.matrix <- 0.5*(affinity.matrix + t(affinity.matrix))
  
  affinity.matrix <- as.matrix(affinity.matrix)
  return(affinity.matrix)
}

get.similarity.matrix.based.on.RWR.in.knn.graph <- function(X) {
  
  return(get.inter.similarity.matrix.based.on.RWR.in.knn.graph(X,X))
}


generate_projection_mappings <- function(X, Y, W, d, mu)
{
  #
  # X is a P*M matrix (P=#feature, M=#samples)
  # Y is a Q*N matrix (Q=#feature, N=#samples)
  # W is a M*N Matrix, where W(i,j)= similarity of Input1(i) and Input2(j)
  # d is the number of dimensions in the projection space
  # mu is a parameter used to balance two goals: matching corresponding pairs and preserving manifold topology
  #
  
  P <- nrow(X)
  Q <- nrow(Y)
  
  # W_x <- get.similarity.matrix.based.on.geodesic.ditance.in.knn.graph(X)
  # W_x <- get.similarity.matrix.based.on.RWR.in.knn.graph(X)
  W_x <- get_similarity_matrix(X)
  diag(W_x) <- 0
  # W_y <- get.similarity.matrix.based.on.geodesic.ditance.in.knn.graph(Y)
  # W_y <- get.similarity.matrix.based.on.RWR.in.knn.graph(Y)
  W_y <- get_similarity_matrix(Y)
  diag(W_y) <- 0
  
  # normalize similarity matrices
  W_x <- W_x / sum(W_x)
  W_y <- W_y / sum(W_y)
  W <- W / sum(W)
  
  D_x <- compute_diagonal_similarity_matrix(W_x)
  D_y <- compute_diagonal_similarity_matrix(W_y)
  
  # normalize mu (probably redundant now, as we normalize the similarity matrices to be of sum=1)
  mu <- mu * (sum(W_x) + sum(W_y)) / (2 * sum(W))
  
  Z <- compute_Z(X, Y)                      # (P+Q)*(M+N) matrix
  D <- compute_D(D_x, D_y)                  # (M+N)*(M+N) matrix
  L <- compute_L(W_x, D_x, W_y, D_y, W, mu) # (M+N)*(M+N) matrix
  
  A <- Z %*% L %*% t(Z)
  B <- Z %*% D %*% t(Z)
  
  # Make A and B symmetric, to overcome the issue of floating point operations percision
  A <- 0.5 * (A + t(A))
  B <- 0.5 * (B + t(B))
  
  s <- svd(B)
  rankB <- rankMatrix(B)
  F <- s$u[, 1:rankB] %*% diag(sqrt(s$d[1:rankB]))
  Fplus <- pinv(F)
  
  T <- Fplus %*% A %*% t(Fplus)
  T <- 0.5 * (T + t(T))
  eigen_decom <- eigen(T, symmetric = TRUE)
  min_eigenvectors <- eigen_decom$vectors[, order(eigen_decom$values)[1:d]]
  min_eigenvectors <- t(Fplus) %*% min_eigenvectors
  
  # alpha is a mapping of p-dimensional vectors in space X to d-dimensional vectors in the new shared space Z
  alpha <- min_eigenvectors[1:P, ] # P*d matrix
  
  # beta is a mapping of q-dimensional vectors in space Y to d-dimensional vectors in the new shared space Z
  beta <- min_eigenvectors[(P + 1):(P + Q),] # Q*d matrix
  
  mappings <- (list(alpha = alpha, beta = beta, eigenvalues = eigen_decom$values[order(eigen_decom$values)])) 
  return(mappings)
}
