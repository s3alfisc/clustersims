run_clustersim <- function(n_sims,
                           workers = 6,
                           G = seq(6, 84, 6),
                           k = 9,
                           rho = 0.5,
                           gamma = 2,
                           filename = NULL,
                           chi_square = FALSE
){

  #' @import ggplot2
  #' @export

  res <- matrix(NA, length(G), 5)

  for(g in seq_along(G)){

      cat("Currently, the number of clusters is", G[g])

     res[g, ] <-
       get_rejection_frequency(
        n_sims = n_sims,
        G = G[g],
        k = k,
        rho = rho,
        gamma = gamma,
        workers = workers,
        chi_square = chi_square
      )

  }

  inferences <- c("CRV1", "CRV3", "bootstrap-C", "bootstrap-S", "bootstrap-V")
  res_df <- data.frame(
    rejection_frequency = c(res),
    Inference = rep(inferences, length(G)),
    G = sort(rep(G, length(inferences)))
  )

  p <- ggplot(res_df, aes(x = G, y = rejection_frequency, group = Inference, color = Inference)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0.05, color = "red", linetype = "dotted") +
    theme_bw() +
    ylab("Rejection Frequency") +
    xlab("Number of Clusters G")

  if(!is.null(filename)){
    path = paste0("results/", filename, ".csv")
    path_plot = paste0("results/", filename, ".png")

    write.csv(res_df, path)
    ggsave(
      filename = path_plot,
      plot = p,
      width = 8,
      height = 5
    )
  }

  p

}


