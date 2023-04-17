#' ### SINGLE ITERATION OF THE SIM: CREATE A SAMPLE, ESTIMATE MODEL, STORE RESULTS
#'
#' # This function nests the data generation function.
#' # For each simulation iteration, simulate clustered data with the provided params.
#' # Then estimate the linear regression model, store the beta, and the standard errors + confidence intervals.
#' # Calculate the inference naively, with regular clustering, or fast and wild cluster bootstrap (good for few clusters).
#' cluster_sim <- function(param = c(.1, .5),
#'                         n = 1000,
#'                         n_cluster = 50,
#'                         rho = .5,
#'                         inference = "regular",
#'                         FE = FALSE,
#'                         balanced_cluster = TRUE){
#'
#'   # Required packages: multiwayvcov, fwildclusterboot
#'
#'   # Note we have to do some weirdness here because of the way the boottest function is coded.
#'   # We need df_it to be a global variable, because boottest looks for it implicitly in global memory.
#'   df_it <<-
#'     create_data(
#'       G,
#'       k,
#'       rho,
#'       gamma
#'     )
#'
#'   x_err_cor <- cor(df_it$x,df_it$error)
#'
#'   if (FE == TRUE) {
#'     fit <- lm(data = df_it, y ~ x + factor(cluster))
#'   } else {
#'     fit <- lm(data = df_it, y ~ x)
#'   }
#'
#'   b1 <- coef(fit)[2]
#'
#'   if (inference == "regular") {
#'     Sigma <- stats::vcov(fit)
#'     b1_ci95 <- stats::confint(fit)[2, ]
#'     pval = summary(fit)[4][[1]]["x", 4]
#'   } else if (inference == "CRV1") {
#'     Sigma <- suppressMessages(clubSandwich::vcovCR(fit,
#'                                                    cluster = df_it$cluster,
#'                                                    type = "CR1S"))
#'     conf_int <- suppressMessages(clubSandwich::conf_int(obj = fit,
#'                                                         vcov = Sigma ,
#'                                                         test = "z",
#'                                                         level = 0.95,
#'                                                         coefs = "x"))
#'     #pval <-
#'
#'   } else if (inference == "fast_n_wild_11") {
#'     # Per Roodman et al. 2019, employ Webb weights if cluster count <= 12.
#'     if (n_cluster <= 12) {
#'       weights = "webb"
#'     } else {
#'       weights = "rademacher"
#'     }
#'     suppressWarnings(boot <-
#'                        fwildclusterboot::boottest(
#'                          fit,
#'                          B = 9999,
#'                          param = "x",
#'                          clustid = "cluster",
#'                          type = weights,
#'                          engine = "R",
#'                          conf_int = FALSE
#'                        ))
#'
#'     pval = boot$pval
#'
#'   } else if (inference == "fast_n_wild_31") {
#'     # Per Roodman et al. 2019, employ Webb weights if cluster count <= 12.
#'     if (n_cluster <= 12) {
#'       weights = "webb"
#'     } else {
#'       weights = "rademacher"
#'     }
#'     suppressWarnings(boot <-
#'                        fwildclusterboot::boottest(
#'                          fit,
#'                          B = 9999,
#'                          param = "x",
#'                          clustid = "cluster",
#'                          type = weights,
#'                          bootstrap_type = "31",
#'                          engine = "R",
#'                          conf_int = FALSE
#'                        ))
#'
#'     pval = boot$pval
#'   } else if(inference == "CRV3"){
#'
#'     Sigma <- summclust::summclust(
#'       fit,
#'       params = "x",
#'       cluster = ~ cluster,
#'       type = "CRV3"
#'     )
#'     conf_int <- tidy(Sigma)[c('conf_int_l', 'conf_int_u')]
#'     pval = tidy(Sigma)['p_val']
#'
#'   }
#'
#'   return(pval)
#' }
#'
#' # Lastly, a function to iterate / repeat the simulation.
#' # A data frame is returned containing all the betas, se's, CIs, etc.
#' # We are also creating, for each sim, a binary indicator of whether the CIs contain the true parameter for b1.
#' run_cluster_sim <-
#'   function(n_sims = 1000,
#'            param = c(.1, 0),
#'            n = 1000,
#'            n_cluster = 50,
#'            rho = .5,
#'            inference = "regular",
#'            FE = FALSE,
#'            balanced_cluster = TRUE,
#'            workers = 8) {
#'     # Required packages: dplyr
#'
#'     cat(paste0("inference: ", inference), "\n")
#'     #if(inference != "fast_n_wild"){
#'     future::plan(multisession, workers = workers)
#'     handlers("progress")
#'
#'     with_progress({
#'       p <- progressor(along = 1:n_sims)
#'       df_res <-
#'         future.apply::future_lapply(X = 1:n_sims,
#'                                     future.seed = sample(1:10000000, 1),
#'                                     FUN = function(x){
#'                                       p(sprintf("x=%g", x))
#'                                       cluster_sim(param = param,
#'                                                   n = n,
#'                                                   rho = rho,
#'                                                   n_cluster = n_cluster,
#'                                                   inference = inference,
#'                                                   FE = FE)
#'
#'                                     })
#'
#'     })
#'
#'     df_res <- as.data.frame(Reduce(rbind, df_res))
#'     names(df_res) <- c('pval')
#'     df_res$nsim = 1:n_sims
#'     df_res$reject <- df_res$pval < 0.05
#'     df_res = df_res[, c("nsim", "reject")]
#'
#'     return(df_res)
#'   }
#'
#'
#' library(data.table)
#' # source https://en.wikipedia.org/wiki/List_of_U.S._states_and_territories_by_population
#'
#' get_state_freq <- function(){
#'   library(data.table)
#'   state_population <-
#'     c("California",	39538223,
#'       "Texas",	29145505,
#'       "Florida",	21538187,
#'       "New York",	20201249,
#'       "Pennsylvania",	13002700,
#'       "Illinois",	12812508,
#'       "Ohio",	11799448,
#'       "Georgia",	10711908,
#'       "North Carolina",	10439388,
#'       "Michigan",	10077331,
#'       "New Jersey",	9288994,
#'       "Virginia",	8631393,
#'       "Washington",	7705281,
#'       "Arizona",	7151502,
#'       "Massachusetts",	7029917,
#'       "Tennessee",	6910840,
#'       "Indiana",	6785528,
#'       "Maryland",	6177224,
#'       "Missouri",	6154913,
#'       "Wisconsin",	5893718,
#'       "Colorado",	5773714,
#'       "Minnesota",	5706494,
#'       "South Carolina",	5118425,
#'       "Alabama",	5024279,
#'       "Louisiana",	4657757,
#'       "Kentucky",	4505836,
#'       "Oregon",	4237256,
#'       "Oklahoma",	3959353,
#'       "Connecticut",	3605944,
#'       "Puerto Rico",	3285874,
#'       "Utah",	3271616,
#'       "Iowa",	3190369,
#'       "Nevada",	3104614,
#'       "Arkansas",	3011524,
#'       "Mississippi",	2961279,
#'       "Kansas",	2937880,
#'       "New Mexico",	2117522,
#'       "Nebraska",	1961504,
#'       "Idaho",	1839106,
#'       "West Virginia",	1793716,
#'       "Hawaii",	1455271,
#'       "New Hampshire",	1377529,
#'       "Maine",	1362359,
#'       "Rhode Island",	1097379,
#'       "Montana",	1084225,
#'       "Delaware",	989948,
#'       "South Dakota",	886667,
#'       "North Dakota",	779094,
#'       "Alaska",	733391,
#'       "DC",	689545,
#'       "Vermont", 643077,
#'       "Wyoming", 576851
#'     )
#'
#'   dt = data.table(state = state_population[seq(1, length(state_population), 2)],
#'                   population = as.numeric(state_population[seq(2, length(state_population), 2)]))
#'   dt <- dt[!(state %in% c("Puerto Rico", "DC"))]
#'   total_population <- sum(dt$population)
#'   dt[, population_share := population / total_population]
#'
#'   dt[, list(state, population_share)]
#' }
#'
#'
#'
#' compute_freject <- function(
#'     rho,
#'     inference,
#'     workers = 8,
#'     n_clusters = 5:10){
#'
#'   #' @export
#'
#'   mean_reject <- rep(0, 49)
#'   for (j in 2:50) {
#'     tryCatch(
#'       {
#'         sim_clustered_fix_tmp <-
#'           run_cluster_sim(rho = rho,
#'                           n_cluster = j,
#'                           inference = inference,
#'                           workers = workers)
#'         mean_reject[j-1] <- mean(sim_clustered_fix_tmp$pval)
#'       },
#'       error=function(cond) {
#'         message(paste("\n Inference failed.\n"))
#'         message("Here's the original error message:")
#'         message(cond)
#'         # Choose a return value in case of error
#'         return(NA)
#'       })
#'   }
#'
#'   res <- data.frame(coverage = mean_reject,
#'                     clusters = seq(2:50))
#'
#'   res
#' }
#'
#'
#' visualize_freject = function(){
#'   # Step 6: Collect all results
#'
#'   perfPalette_fw <- colorRampPalette(brewer.pal(9, "Greens"))(length(unique(cluster_SE_by_clustN_fw$coverage)))
#'   unique_coverages_fw <- cluster_SE_by_clustN_fw %>%
#'     group_by(coverage) %>%
#'     summarize(coverage = first(coverage)) %>%
#'     arrange(desc(coverage))
#'   unique_coverages_fw$heat <- perfPalette_fw
#'
#'   cluster_SE_by_clustN_fw2 <- cluster_SE_by_clustN_fw %>%
#'     merge(unique_coverages_fw, by="coverage")
#'
#'   # Now we can plot the SE performance side-by-side.
#'   ggplot(data = cluster_SE_by_clustN_fw2, aes(y=coverage*100,x=(clusters+1),alpha=0.7)) +
#'     geom_hline(yintercept=95,color="red",size=1,linetype="dashed")+
#'     #geom_vline(xintercept=30,color="black",size=1)+
#'     geom_line(aes(color="green"),size=2) +
#'     geom_line(data=cluster_SE_by_clustN, aes(color="blue"),size=2)+
#'     geom_line(data=regular_SE_by_clustN, aes(color="purple"),size=2)+
#'     #geom_line(data = cluster_SE_by_clustN_saddlepoint, aes(color = "violet"), size = 2) +
#'     geom_line(data = cluster_SE_by_clustN_satterthwaite, aes(color = "orange"), size = 2) +
#'     xlab("# of Clusters (1000 Sims Each)") +
#'     ylab("Avg 95% CI Coverage of True Beta") +
#'     theme_classic() +
#'     theme(
#'       text = element_text(size=18,family = "Economica"),
#'       axis.title = element_text(size = 16, margin = margin(
#'         t = 10,
#'         b = 0,
#'         r = 20,
#'         l = 0
#'       )),
#'       axis.text = element_text(size = 14)
#'     ) +
#'     ylim(-5,97) +
#'     scale_x_continuous(breaks=c(2,10,20,30,40,50))+
#'     scale_y_continuous(breaks=seq(5,95,by=10))+
#'     scale_alpha_continuous(guide=FALSE)+
#'     #ggtitle("") +
#'     scale_color_manual(name="SE Type",values=c("blue","green", "purple", "orange"),labels=c("Clustered","Wild Clustered", "Satterthwaite","Homoskedastic iid"))
#'   NULL
#'
#' }
#'
