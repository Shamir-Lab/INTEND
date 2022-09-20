compute_Z <- function(X, Y)
{
  # X is a P*M matrix 
  P <- nrow(X)
  M <- ncol(X)
  
  # Y is a Q*N matrix 
  Q <- nrow(Y)
  N <- ncol(Y)
  
  # zeros_top_right is a P*N matrix
  zeros_top_right <- matrix(data = 0, nrow = P, ncol = N)
  
  # zeros_bottom_left is a Q*M matrix
  zeros_bottom_left <- matrix(data = 0, nrow = Q, ncol = M)
  
  # Z is a (P+Q)*(M+N) matrix
  Z <- rbind(cbind(X, zeros_top_right), cbind(zeros_bottom_left, Y))
  return(Z)  
}

compute_diagonal_similarity_matrix <- function(similarity_matrix)
{
  # similarity_matrix[i,j] contains the similarity between sample i and sample j diagonal_similarity_matrix is
  # a diagonal matrix, and diagonal_similarity_matrix[i,i] = sum(similarity_matrix[i,])
  diagonal_similarity_matrix <- diag(apply(similarity_matrix, 1, sum))
  return(diagonal_similarity_matrix)
}

compute_D <- function(D_x, D_y)
{
  # D_x is a diagonal M*M matrix 
  M <- ncol(D_x)
  
  # D_y is a diagonal N*N matrix 
  N <- ncol(D_y)
  
  # zeros_top_right is a M*N matrix
  zeros_top_right <- matrix(data = 0, nrow = M, ncol = N)
  
  # zeros_bottom_left is a N*M matrix
  zeros_bottom_left <- t(zeros_top_right)
  
  # D is a (M+N)*(M+N) matrix
  D <- rbind(cbind(D_x, zeros_top_right), cbind(zeros_bottom_left, D_y))
  return(D)
}

compute_L <- function(W_x, D_x, W_y, D_y, W, mu)
{
  # W: M*N Matrix, where W[i,j] = similarity of local geometries of samples X_i and Y_j
  # mu is a parameter used to balance two goals: matching corresponding pairs and preserving manifold topology   

  # M*M matrix (M=#samples in X)  
  L_x <- D_x - W_x
  
  # N*N matrix (N=#samples in Y)  
  L_y <- D_y - W_y
  
  # omega_1 is an M*M diagonal matrix where omega_1[i,i] = sum over j of W[i,j]
  omega_1 <- diag(apply(W, 1, sum))
  
  # omega_2 is an M*N matrix, and omega_2 = W
  omega_2 <- W
  
  # omega_3 is an N*M matrix, and omega_3 = t(W)
  omega_3 <- t(W)
  
  # omega_4 is an N*N diagonal matrix where omega_4[i,i] = sum over j of W[j,i]
  omega_4 <- diag(apply(W, 2, sum))
  
  # L is a (M+N)*(M+N) matrix
  L <- rbind(cbind((L_x + mu * omega_1), (-mu * omega_2)),
              cbind((-mu * omega_3),       (L_y + mu * omega_4)))
  
  return(L)
}