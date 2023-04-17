create_data = function(G, k, rho, gamma, chi_square){

  N <- 400 * G
  N_g <- rep(0, G)

  denom_tmp <- sum(exp(gamma * (1:G) / G))
  for(g in 1:G-1){
    N_g[g] <- floor(N * exp(gamma * g / G) / denom_tmp)
  }
  N_g[G] <- N - sum(N_g)

  clusters <- factor(rep(1:G, N_g))

  X <- matrix(NA, N, k)
  #X <- mvtnorm::rmvnorm(N, rep(0, k), diag(k))
  X[,1] <- 1
  for(i in 2:k){
    X[,i] = fabricatr::draw_normal_icc(clusters = clusters, ICC = rho)
  }

  u <- fabricatr::draw_normal_icc(clusters = clusters, ICC = 0.1)

  beta <- rep(1, k)
  beta[2] <- 0

  if(chi_square){
    X[,2] <- X[,2]^2
  }

  y <- X %*% beta + u


  list(
    X = X,
    y = y,
    cluster = clusters,
    beta = beta
  )
  #df = data.frame(cbind(y, X[,-1]), clusters)
  #names(df) <- c('y', paste0("X", 1:(k-1)), "cluster")



}


