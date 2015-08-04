measles <- function(input = c(1, 0.97, 30))
{

    myparams <- model_params
    seed <- input[1]
    if (length(input) > 1) { myparams$alpha = input[2] }
    if (length(input) > 2) { myparams$beta = input[3] }
    if (length(input) > 3) { myparams$theta = input[4] }
    if (length(input) > 4) { myparams$rho = input[5] }
    if (length(input) > 5) { myparams$tau1 = input[5] }
    if (length(input) > 6) { myparams$tau2 = input[6] }

    accept <- 1

    summaryStats <- c()

    result <-
        tsir(model_params = myparams, model_options = model_options,
             sim_params = sim_params, init = init, summary = summary,
             data = data, seed = seed)

    if (result$extinction)
        {
            accept <- 0
        }
    
    c(accept, result$summaryStats)
}

scaledmeasles <- function(input = c(1, 0.97, 30))
{

    myparams <- model_params
    seed <- input[1]
    scale <- 1
    if (length(input) > 1) { myparams$alpha = input[2] }
    if (length(input) > 2) { myparams$beta = input[3] }
    if (length(input) > 3) { init$S = input[4] }
    if (length(input) > 4) { scale = input[5] }

    accept <- 1

    summaryStats <- c()

    result <-
        tsir(model_params = myparams, model_options = model_options,
             sim_params = sim_params, init = init, summary = summary,
             data = data, seed = seed)

    if (result$extinction)
        {
            accept <- 0
        }
    
    c(accept, scale * result$summaryStats)
}

seasonalmeasles <- function(input = c(1, 0.97, rep(1, 26)))
{

    myparams <- model_params
    seed <- input[1]
    if (length(input) > 1) { myparams$alpha = input[2] }
    if (length(input) > 2) { data$seasonality = exp(input[3:28]) }
    if (length(input) > 28) { myparams$theta = input[29] }
    if (length(input) > 29) { myparams$rho = input[30] }
    if (length(input) > 30) { myparams$tau1 = input[31] }
    if (length(input) > 31) { myparams$tau2 = input[32] }

    accept <- 1

    summaryStats <- c()

    result <-
        tsir(model_params = myparams, model_options = model_options, sim_params = sim_params,
             init = init, summary = summary, data = data, seed = seed)

    if (result$extinction)
        {
            accept <- 0
        }
    
    c(accept, result$summaryStats)
}

seasonalscaledmeasles <- function(input = c(1, 0.97, rep(1, 26), 0.52, 250000))
{

    myparams <- model_params
    seed <- input[1]
    scale <- 1
    if (length(input) > 1) { myparams$alpha = input[2] }
    if (length(input) > 2) { data$seasonality = exp(input[3:28]) }
    if (length(input) > 28) { scale = input[29] }
    if (length(input) > 29) { init$S = input[30] }
#    if (length(input) > 28) { myparams$theta = input[29] }
#    if (length(input) > 29) { myparams$rho = input[30] }
#    if (length(input) > 30) { myparams$tau1 = input[31] }
#    if (length(input) > 31) { myparams$tau2 = input[32] }

    accept <- 1

    summaryStats <- c()

    result <-
        tsir(model_params = myparams, model_options = model_options, sim_params = sim_params,
             init = init, summary = summary, data = data, seed = seed)

    if (result$extinction)
        {
            accept <- 0
        }

    # fix 
    c(accept, scale * result$summaryStats)
}

iterate_mcmc <- function(init, index=0, saveRun=T, tab_normalization=stat, stat=stat, dist_max=0, n_rec = 10) {

    ABC_Marjoram <- ABC_mcmc(method = "Marjoram_original",
                             model = seasonalscaledmeasles,
                             prior = prior,
                             summary_stat_target = stat,
                             n_rec = n_rec,
                             n_between_sampling = 5,
                             use_seed = T,
                             verbose = T,
                             dist_max = dist_max,
                             tab_normalization = tab_normalization,
                             progress_bar = T,
                             rejection = T,
                             init_param = init,
                             inside_prior = F)

    ABC_Marjoram$param[,2:27] <- exp(ABC_Marjoram$param[,2:27])

    if (saveRun) {
        saveRDS(ABC_Marjoram, file = paste("marjoram", index, ".rds", sep=""))
    }

    index <- which(ABC_Marjoram$dist==min(ABC_Marjoram$dist))
    list(param = ABC_Marjoram$param[min(index),], dist= ABC_Marjoram$dist[min(index)])
}
