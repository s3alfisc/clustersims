cluster_sim <- function(G, k, rho, gamma, chi_square){

  # dreamerr::check_arg(G, "numeric")
  # dreamerr::check_arg(k, "rho")
  # dreamerr::check_arg(rho, "numeric")
  # dreamerr::check_arg(gamma, "numeric")
  # dreamerr::check_arg(B, "integer")



  data <- create_data(
              G = G,
              k = k,
              rho = rho,
              gamma = gamma,
              chi_square = chi_square
          )

  y = data$y
  X = data$X
  cluster = data$cluster
  beta_true = data$beta

  data = data.frame(y = y, X = X, cluster = cluster)
  names(data) <- c("y", paste0("X", 1:k), "cluster")
  G <- length(unique(cluster))

  fit <- lm(y ~ X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9, data = data)


  if(G <= 9){
    type <- "webb"
  } else {
    type <- "rademacher"
  }
  boot_c <- boottest(
    fit,
    param = "X2",
    clustid = ~cluster,
    #bootstrap_type = "11",
    B = 399,
    conf_int= FALSE
  )

  p_val_C <- boot_c$p_val
  p_val_crv1 <- 2*pnorm(-abs(boot_c$t_stat))

  p_val_S <- boottest(
    fit,
    param = "X2",
    clustid = ~cluster,
    bootstrap_type = "31",
    B = 399,
    conf_int= FALSE
  )$p_val

  boot_v <- boottest(
    fit,
    param = "X2",
    clustid = ~cluster,
    bootstrap_type = "13",
    B = 399,
    conf_int= FALSE
  )

  p_val_V <- boot_v$p_val
  p_val_crv3 <- 2 * min(pt(boot_v$t_stat, G - 1), 1 - pt(boot_v$t_stat, G - 1))

  p_vals <- c(
    p_val_crv1,
    p_val_crv3,
    p_val_C,
    p_val_S,
    p_val_V
  )

  p_vals

}



get_rejection_frequency <- function(
    n_sims,
    G,
    k,
    rho,
    gamma,
    workers,
    chi_square
  ){

  #' compute rejection frequencies
  #' @importFrom future plan multisession
  #' @importFrom future.apply future_lapply
  #' @importFrom progressr progressor handlers
  #'
  #' @examples
  #' get_rejection_frequency(n_sims = 10, G = 10, k = 9, rho = 0.1, gamma = 2, workers = 1)

  future::plan(future::multisession, workers = workers)
  progressr::handlers("progress")

  with_progress({
    p <- progressr::progressor(along = 1:n_sims)
    res <-
      future.apply::future_lapply(X = 1:n_sims,
                                  future.seed = sample(1:10000000, 1),
                                  FUN = function(x){
                                    p(sprintf("x=%g", x))
                                    cluster_sim(
                                      G = G,
                                      k = k,
                                      rho = rho,
                                      gamma = gamma,
                                      chi_square = chi_square
                                      )
                                  })

  })

  df_res <- as.data.frame(Reduce(rbind, res))
  names(df_res) <- c("CRV1", "CRV3", "bootstrap-C", "bootstrap-S", "bootstrap-V")
  #df_res$n_sims = 1:n_sims
  rejection_frequency <- colMeans(df_res < 0.05)

  rejection_frequency


}

#
# get_rejection_frequency(
#   n_sims = 100,
#   G = 24,
#   k = 9,
#   rho = 0.51,
#   gamma = 1,
#   inference = "bootstrap-S",
#   workers = 1
# )
#
#
