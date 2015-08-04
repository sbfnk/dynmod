

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

eval.theta <- function() {
    theta <-
        data.table(
            read.csv(
                "~/Research/Analysis/Measles/EnglandWales/theta_estimates_4494.csv", header=T, sep=','
                )
            )

    ggplot(theta, aes(x=theta, y=R))+
                                        #  geom_point()+
        stat_summary(fun.y="mean", geom="point")+
                                        #  geom_line(aes(x=theta, y=n/354*0.56))+
            facet_grid(beta~.)+
                scale_color_brewer(palette="Set1")+
                    scale_x_log10()
}

eval.extinctions <- function() {
    ex <- data.table(read.csv("~/Research/Analysis/Measles/EnglandWales/extinctions_4494_30.csv"))

    ex.histograms <- data.table(coverage=numeric(),
                                year=numeric(),
                                count=numeric(),
                                cumcount=numeric())

    for (current.coverage in unique(ex$coverage)) {
        h <- hist(ex[coverage==current.coverage]$year, breaks=seq(1950,2100,5), plot=F)
        ex.histograms <- rbind(ex.histograms,
                               data.table(coverage=rep(current.coverage, length(h$counts)),
                                          year=h$mids,
                                          count=h$counts,
                                          cumcount=cumsum(h$counts))
                               )
    }

    ex.histograms$uptake <- ex.histograms$coverage/100
    pdf("extinctions_4494.pdf")
    ggplot(ex.histograms, aes(x=year, y=count/100, fill=factor(coverage)))+
        theme_bw()+
            geom_bar(binwidth=5, color="black", stat="identity")+
                geom_vline(xintercept=2012, linetype="longdash")+
                    facet_grid(coverage~.)+
                        theme(legend.position="none", panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
                            scale_y_continuous("Fraction of extinctions", limits=c(0,1))+
                                geom_line(aes(x=year, y=cumcount/100))
    dev.off()
}

osd <- function() {
    hist.data <- data.table(size=integer(0),
                            rfrequency=integer(0),
                            frequency=integer(0),
                            uptake=numeric(0))

    for (i in 85:100) {
        cat(i,"\n")
        osdfile <- paste("~/Research/Analysis/Measles/osd_last10_4494_30_", i, ".csv", sep="")

        if (file.info(osdfile)$size > 0) {
            osd <- data.table(
                t(
                    read.table(
                        osdfile,
                        fill=TRUE,
                        sep=",",
                        row.names=1,
                        col.names=1:max(
                        count.fields(
                            osdfile,
                            sep=",")
                        )
                        )
                    )
                )

            dt <- data.table(size=osd[[1]])[size>0]

            hist.data <- rbind(
                hist.data,
                data.table(
                    size=rev(unique(dt$size)),
                    rfrequency=recdf(dt$size)(rev(unique(dt$size)))*length(dt$size),
                    frequency=recdf(dt$size)(rev(unique(dt$size))),
                    uptake=i
                    )
                )
        }
    }

    pdf("size_freq_4494.pdf", width=10)
    ggplot(hist.data, aes(x=size, y=frequency, color=factor(uptake)))+
        geom_point()+
            scale_x_log10("Outbreak size")+
                scale_y_log10("Frequency")+
                    theme_bw(20)+
                        theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
                            scale_color_discrete("Vaccine uptake")
    dev.off()
}

