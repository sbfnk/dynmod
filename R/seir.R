##' Takes a vector of initial conditions (in S, E(optionally), I, R)
##' and returns a vector of initial conditions to use with S(E)IR
##' models of arbitrary number of infectious stages and agegroups
##'
##' Multiple exposed or infectious stages might be introduced to
##' change the distribution of waiting times, or because there are
##' biologically distinct stages. Multiple age groups might be
##' introduced to follow infection across age groups.
##'
##' The number of age groups is determined by the number of S classes
##' in the state vector
##'
##' @title Initial conditions for the SIR model
##' @param state named state vector
##' @param parameters parameters of the model
##' @param additional additional classes
##' @return a named vector of initial conditions
##' @author Sebastian Funk
initSIR <- function(state, parameters,
                    additional = c("Z")) {

    nameVec <- c()
    initVec <- c()

    if (is.null(parameters[["nE"]]))
    {
        nE <- 0
    } else
    {
        nE <- parameters[["nE"]]
    }
    
    if (is.null(parameters[["nC"]]))
    {
        nC <- 0
    } else
    {
        nC <- parameters[["nC"]]
    }
    
    if (is.null(parameters[["nI"]]))
    {
        nI <- 1
    } else
    {
        nI <- parameters[["nI"]]
    }
    
    if ("B" %in% names(state)) {
        Bvec <- c("B")
        nameVec <- c(nameVec, "B")
        initVec <- c(initVec, state[Bvec])
    }

    # get number of age groups from the number of S classes
    if ("S" %in% names(state)) {
        nAgeGroups <- 1
        names(state) <- paste(names(state), "1", sep = "")
    } else {
        nAgeGroups <- length(grep("S[0-9]+$", names(state)))
    }

    # construct vector of names and initial conditions
    Svec <- c(paste(rep("S", nAgeGroups), seq(1, nAgeGroups), sep = ""))
    nameVec <- c(nameVec, Svec)
    initVec <- c(initVec, state[Svec])
    Evec <- paste(rep("E", nAgeGroups), seq(1, nAgeGroups), sep = "")
    if (nE > 0) {
        nameVec <- c(nameVec, paste(rep(Evec, each = nE),
                                    rep(seq(1, nE), nAgeGroups), sep = "."))
        initVec <- c(initVec,
                     c(t(matrix(c(state[Evec], rep(0, (nE - 1) * nAgeGroups)),
                                ncol = nE))))
    }
    Cvec <- paste(rep("C", nAgeGroups), seq(1, nAgeGroups), sep = "")
    if (nC > 0) {
        nameVec <- c(nameVec, paste(rep(Cvec, each = nC),
                                    rep(seq(1, nC), nAgeGroups), sep = "."))
        initVec <- c(initVec,
                     c(t(matrix(c(state[Cvec], rep(0, (nC - 1) * nAgeGroups)),
                                ncol = nC))))
    }
    Ivec <- paste(rep("I", nAgeGroups), seq(1, nAgeGroups), sep = "")
    nameVec <- c(nameVec, paste(rep(Ivec, each = nI),
                                rep(seq(1, nI), nAgeGroups), sep = "."))
    initVec <- c(initVec,
                 c(t(matrix(c(state[Ivec], rep(0, (nI - 1) * nAgeGroups)),
                            ncol = nI))))
    Rvec <- paste(rep("R", nAgeGroups), seq(1, nAgeGroups), sep = "")
    
    nameVec <- c(nameVec, Rvec)
    initVec <- c(initVec, state[Rvec])
        
    if (!is.null(parameters[["hospital.entry"]]))
    {
        nameVec <- c(nameVec,
                     paste(rep("H", nAgeGroups),
                           seq_len(nAgeGroups), sep = ""))
        initVec <- c(initVec, rep(0, nAgeGroups))
    }
    
    if (!is.null(parameters[["antibody.development"]]))
    {
        nameVec <- c(nameVec,
                     paste(rep("A", nAgeGroups),
                           seq_len(nAgeGroups), sep = ""))
        initVec <- c(initVec, rep(0, nAgeGroups))
    }

    if (!is.null(parameters[["vaccine.campaign"]]) |
        !is.null(parameters[["vaccine.uptake"]]))
    {
        nameVec <- c(nameVec,
                     paste(rep("V", nAgeGroups),
                           seq_len(nAgeGroups), sep = ""))
        initVec <- c(initVec, rep(0, nAgeGroups))
    }
    if (!is.null(parameters[["booster.campaign"]]))
    {
        nameVec <- c(nameVec,
                     paste(rep("W", nAgeGroups),
                           seq_len(nAgeGroups), sep = ""))
        initVec <- c(initVec, rep(0, nAgeGroups))
    }

    if (length(additional) > 0) {
        for (i in 1:length(additional)) {
            nameVec <- c(nameVec,
                         paste(rep(additional, nAgeGroups),
                               seq_len(nAgeGroups), sep = ""))
            initVec <- c(initVec, rep(0, nAgeGroups))
        }
    }

    names(initVec) <- nameVec
    return(initVec)
}

##' Runs an SIR model
##'
##' Runs an SIR model and returns the trajectory and initial state.
##' @param parameters SIR parmaeters
##' @param init initial state
##' @param times times at which to produce the trajectory
##' @param age.labels labels of age groups
##' @param thickening whether to thicken the data (use finer integration steps)
##' @param compiled whether to use the compiled version of the SIR model
##' @param jacobian whether to use the Jacobian
##' @return The initial state and trajectory
##' @export
##' @useDynLib dynmod
##' @author Sebastian Funk
##' @import deSolve
runSIR <- function(parameters, init, times, age.labels = NULL,
                   thickening = 1)
{

    parameters <- as.list(parameters)

    if (thickening > 1)
    {
        integration.times <- c(sapply(seq_len(length(times) - 1), function (x)
        {
            interval <- (times[x + 1] - times[x]) / thickening
            c(seq(times[x], times[x + 1] - interval,
                  by = interval))
        }))
        integration.times <- c(integration.times, times[length(times)])
    } else {
        integration.times <- times
    }

    if (is.null(parameters[["mixing"]]))
    {
        parameters[["mixing"]] <- 1
    }
    parameters[["mixing"]] <- matrix(parameters[["mixing"]],
                                     nrow = sqrt(length(parameters[["mixing"]])))

    if (!is.null(parameters[["mixing.dynamics"]]))
    {
        if (parameters[["mixing.dynamics"]] == "brownian")
        {
             # generate brownian motion
             parameters[["child.mixing"]] <-
                  brownian(parameters,
                           max(times) - min(times) + 1)

         } else if (parameters[["mixing.dynamics"]] == "linear") {

            mixing.values <-
                unname(unlist(parameters[sort(grep("^mixing\\.[0-9]+$",
                                                   names(parameters), value = T))]))
            change.years <- 
                unname(unlist(parameters[sort(grep("^change\\.[0-9]+$",
                                                   names(parameters), value = T))]))
            parameters[["child.mixing"]] <-
                piecewise.linear(mixing.values = mixing.values,
                                 change.times = change.years - min(times) + 1,
                                 n = max(times) - min(times) + 1)
        }
    }

    if (is.null(parameters[["nE"]]))
    {
        if (!is.null(parameters[["rho"]]))
        {
            parameters[["nE"]] <- as.integer(1)
        } else {
            parameters[["nE"]] <- as.integer(0)
        }
    }

    if (is.null(parameters[["nC"]]))
    {
        if (!is.null(parameters[["epsilon"]]))
        {
            parameters[["nC"]] <- as.integer(1)
        } else {
            parameters[["nC"]] <- as.integer(0)
        }
    }

    if (is.null(parameters[["nI"]]))
    {
         parameters[["nI"]] <- as.integer(1)
    }

    init.state <- initSIR(init, parameters)

    mixing.R0 <- max(eigen(parameters[["mixing"]])$values)

    parameters[["beta"]] <- parameters[["gamma"]] * parameters[["R0"]] / mixing.R0

    jacfunc <- NULL
    func <- "seir_derivs"
    dllname <- ""
    initfunc <- "seir_initmod"
    agegroups <- nrow(parameters[["mixing"]])
    outnames <- c(paste("N", seq(1, agegroups), sep = "."), "N")
    nout <- length(outnames)

    if ("runin.time" %in% names(parameters) && parameters[["runin.time"]] > 0)
    {
        first.diff <- times[2] - times[1]
        runin.step <- ceiling(first.diff)
        
        runin.times <-
            c(seq(times[1] - parameters[["runin.time"]] * runin.step,
                  times[1], by = runin.step / thickening))
        runin.parameters <- parameters

        for (vector in c("births", "deaths", "child.mixing"))
        {
            if (!is.null(parameters[[vector]]))
            {
                runin.parameters[[vector]] <- parameters[[vector]][1]
            }
        }

        for (matrix in c("vaccine.uptake"))
        {
            if (!is.null(parameters[[matrix]]))
            {
                parameters[[matrix]] <- parameters[[matrix]][1, ]
            }
        }

        runin.traj <- data.table(ode(y = init.state, times = runin.times,
                                     func = func, parms = runin.parameters,
                                     hini = 1, jacfunc = jacfunc,
                                     dllname = dllname, initfunc = initfunc,
                                     nout = nout, outnames = outnames,
                                     method = "ode45"))

        init.state <- unlist(runin.traj[nrow(runin.traj),
                                        seq(2, length(init.state) + 1), with = F])
        init.state[grep("^Z[0-9]+$", names(init.state))] <- 0
        # rescale populations
        
        scaling <- unlist(runin.traj[1,
                                     grep("^N\\.[0-9]+",
                                        colnames(runin.traj)), with = F]) /
                       unlist(runin.traj[nrow(runin.traj),
                                         grep("^N\\.[0-9]+",
                                            colnames(runin.traj)), with = F])

        init.scaling <- rep(scaling, length(init.state) %/% length(scaling))
        init.diff <- length(init) - length(init.scaling)
        if (init.diff > 0)
        {
            init.scaling <- c(rep(init.scaling[1], init.diff), init.scaling)
        }

        init.state <- init.state * init.scaling
    }

    traj <- data.table(ode(y = init.state, times = integration.times,
                           func = func, parms = parameters,
                           hini = 1, jacfunc = jacfunc,
                           dllname = dllname, initfunc = initfunc,
                           nout = nout, outnames = outnames,
                           method = "ode45"))

    traj <- traj[time %in% times]
        
    mtraj <- melt.trajectory(traj, age.labels = age.labels)
    mtraj <- mtraj[!is.na(abs.incidence)]

    return(mtraj)
}

##' calculates a stochastic trajectory of the S(E)IR model
##'
##' This is being used by \code{\link{runSIR}}.
##' @param parameters parameters of the S(E)IR mode
##' @param init initial state
##' @param seed random seed. If set to NA, will use R's current seed
##' @param transient transient time. If this is >0, the model will be
##' run deterministically for the time specified here before running
##' the stochastic model. Only the outcome of the stochastic model is
##' recorded
##' @param init.params whether the initial conditions are passed as parameters
##' @param tmax the maximum time to which to run the model (form 0)
##' @return a stochastic trajectory of the model.
##' @author Sebastian Funk.
seir.stochastic <- function(parameters, init, seed = NA,
                            transient = 0, tmax = 2128) {

    nAgeGroups <- nrow(parameters$mixing)

    statenames <- c(paste(rep("S", nAgeGroups), seq(1, nAgeGroups), sep = ""),
                    paste(paste(rep("E", nAgeGroups), seq(1, nAgeGroups),
                                sep = ""),
                          rep(seq(1, parameters$nI), nAgeGroups), sep = "."),
                    paste(paste(rep("I", nAgeGroups), seq(1, nAgeGroups),
                                sep = ""),
                          rep(seq(1, parameters$nI), nAgeGroups), sep = "."),
                    paste(rep("R", nAgeGroups), seq(1, nAgeGroups), sep = ""),
                    paste(rep("Z", nAgeGroups), seq(1, nAgeGroups), sep = ""))

    trans.list <- list(statenames)

    # count transitions, start at 2 because first parameter to
    # ssa.maketrans is the vector of state names
    counter <- 2

    if (parameters$nE > 0) {
        infections.into <- "E"
    } else {
        infections.into <- "I"
    }

    # births
    for (i in seq(1, nAgeGroups)) {
        trans.list[[counter]] <- rbind(paste("S", i, sep = ""), +1)
        counter <- counter + 1
    }

    # susceptible deaths
    for (i in seq(1, nAgeGroups)) {
        trans.list[[counter]] <- rbind(paste("S", i, sep = ""), -1)
        counter <- counter + 1
    }

    # infections
    for (i in seq(1, nAgeGroups)) {
        trans.list[[counter]] <-
            rbind(paste("S", i, sep = ""), -1,
                  paste(infections.into, i, ".1", sep = ""), +1,
                  paste("Z", i, sep = ""), +1)
        counter <- counter + 1
    }

    # progression through exposed states
    for (i in seq(1, nAgeGroups)) {
        if (parameters$nE > 1) {
            for (j in seq(1, parameters$nE - 1)) {
                trans.list[[counter]] <-
                    rbind(paste("E", i, ".", j, sep = ""), -1,
                          paste("E", i, ".", j + 1, sep = ""), +1)
                counter <- counter + 1
            }
        }
        # symptomatic infection
        if (parameters$nE > 0) {
            trans.list[[counter]] <-
                rbind(paste("E", i, ".", parameters$nE, sep = ""), -1,
                      paste("I", i, ".1", sep = ""), +1)
            counter <- counter + 1
        }
    }

    # exposed deaths
    for (i in seq(1, nAgeGroups)) {
        for (j in seq(1, parameters$nE)) {
            trans.list[[counter]] <- rbind(paste("E", i, ".", j, sep = ""), -1)
            counter <- counter + 1
        }
    }

    # progression through infected states
    for (i in seq(1, nAgeGroups)) {
        if (parameters$nI > 1) {
            for (j in seq(1, parameters$nI - 1)) {
                trans.list[[counter]] <-
                    rbind(paste("I", i, ".", j, sep = ""), -1,
                          paste("I", i, ".", j + 1, sep = ""), +1)
                counter <- counter + 1
            }
        }
        # recovery
        trans.list[[counter]] <-
            rbind(paste("I", i, ".", parameters$nI, sep = ""), -1,
                  paste("R", i, sep = ""), +1)
        counter <- counter + 1

    }

    # infected deaths
    for (i in seq(1, nAgeGroups)) {
        for (j in seq(1, parameters$nI)) {
            trans.list[[counter]] <- rbind(paste("I", i, ".", j, sep = ""), -1)
            counter <- counter + 1
        }
    }

    # recovered deaths
    for (i in seq(1, nAgeGroups)) {
        trans.list[[counter]] <- rbind(paste("R", i, sep = ""), -1)
        counter <- counter + 1
    }

    transitions <- do.call(ssa.maketrans, trans.list)

    rate.function <- function(x, parameters, time) {
        with(as.list(parameters), {

            rates <- c()

            if (!is.null(parameters$mixing)) {
                current.mixing <- mixing
            } else {
                current.mixing <- 1
            }

            if (!is.null(parameters$termtime.forcing)) {

                if (time %% 52 >= termtime.forcing$holiday.start &&
                        time %% 52 < termtime.forcint$holiday.end) {
                    current.mixing <- termtime.forcing$holiday.mixing
                }
            }

            if (!is.null(parameters$sinusoidal.forcing)) {

                current.mixing <- current.mixing *
                    sinusoidal.forcing$amplitude *
                        sin(2 * pi * (time - sinusoidal.forcing$t0)/
                                sinusoidal.forcing$period)
            }

            current.births <- births[min(time %/% 52 + 1, length(births))]
            current.mu <- mu[min(time %/% 52 + 1, length(mu))]

            N <- rep(0, nAgeGroups)
            for (i in 1:nAgeGroups) {
                N[i] <- x[paste("S", i, sep = "")] +
                    x[paste("R", i, sep = "")]
                for (j in 1:nE) {
                    N[i] <- N[i] + x[paste("E", i, ".", j, sep = "")]
                }
                for (j in 1:nI) {
                    N[i] <- N[i] + x[paste("I", i, ".", j, sep = "")]
                }
            }

            Isum <- rep(0, nAgeGroups)
            for (i in 1:nAgeGroups) {
                for (j in 1:nE) {
                    Isum[i] <- Isum[i] + x[paste("I", i, ".", j, sep = "")]
                }
            }

            # births
            for (i in seq(1, nAgeGroups)) {
                rates <- c(rates, current.births)
            }

            # susceptible deaths
            for (i in seq(1, nAgeGroups)) {
                rates <- c(rates, current.mu * x[paste("S", i, sep = "")])
            }

            # infections
            for (i in seq(1, nAgeGroups)) {
                foi <- 0
                for (j in 1:nAgeGroups) {
                    foi <- foi + beta * current.mixing[i,j] * Isum[j] / N[j]
                }

                rates <- c(rates, foi * x[paste("S", i, sep = "")])
            }

            # progression through exposed states
            for (i in seq(1, nAgeGroups)) {
                for (j in seq(1, parameters$nE)) {
                    rates <- c(rates, nE * rho *
                                   x[paste("E", i, ".", j, sep = "")])
                }
            }

            # exposed deaths
            for (i in seq(1, nAgeGroups)) {
                for (j in seq(1, parameters$nE)) {
                    rates <- c(rates, current.mu * x[paste("E", i, ".", j,
                                                           sep = "")])
                }
            }

            # progression through infected states
            for (i in seq(1, nAgeGroups)) {
                for (j in seq(1, parameters$nI)) {
                    rates <- c(rates, nI * gamma *
                                   x[paste("I", i, ".", j, sep = "")])
                }
            }

            # infected deaths
            for (i in seq(1, nAgeGroups)) {
                for (j in seq(1, parameters$nI)) {
                    rates <- c(rates, current.mu * x[paste("I", i, ".", j,
                                                           sep = "")])
                }
            }

            # recovered deaths
            for (i in seq(1, nAgeGroups)) {
                rates <- c(rates, current.mu * x[paste("R", i, sep = "")])
            }

            return(unname(rates))

        })
    }

    if (!is.na(seed)) {
        set.seed(seed)
    }

    true.init <- init

    if (transient > 0) {
        transient.parameters <- parameters
        transient.parameters$births <- parameters$births[1]
        transient.parameters$mu <- parameters$mu[1]
        transient <- (transient + 52 - (transient %% 52))
        transient.ode <- runSIR(parameters, seq(1, transient), init)
        true.init <- round(unlist(transient.ode[nrow(transient.ode),
                                                2:(ncol(transient.ode) - 1),
                                                with = F]))
        names(true.init) <- names(init)

        for (i in seq(1, nAgeGroups)) {
            true.init[paste("Z", i, sep = "")] <- 0
        }

        true.init[1] <- true.init[1] + N - sum(true.init)
    }


    ssa.adaptivetau(init = true.init, transitions, rate.function, parameters,
                    tf = tmax)
}

##' Run MCMC on the SEIR model
##'
##' @param target target function (takes one parameter, a numeric vector)
##' @param init.theta initial value of the target function
##' @param proposal.sd standard deviation of the proposal kernel
##' @param n.iterations number of iterations to run the MCMC for
##' @param print.info.every print info every n steps
##' @param acceptance.rate.window window to calculate acceptance rate
##' @param verbose if TRUE, print verbose output
##' @param max.scaling.sd maximal scaling of the standard dveiation in the adaptive MCMC
##' @param adapt.size.start after how many iterations to start adapting the size
##' @param adapt.size.cooling cooling schedule of size adaptation
##' @param adapt.shape.start after how many samples to start adapting the shape
##' @return the trace
##' @import mvtnorm
##' @export
##' @author Anton Camacho, Sebastian Funk
seir.mcmc <- function(target, init.theta, proposal.sd = NULL,
                     n.iterations,
                     print.info.every = n.iterations/100,
                     acceptance.rate.window = NULL,
                     verbose = FALSE, max.scaling.sd = 50, 
                     adapt.size.start = NULL, adapt.size.cooling = 0.99,
                     adapt.shape.start = NULL) {

    # initialise theta
    theta.current <- init.theta
    theta.propose <- init.theta

    # reorder vector and matrix by names, set to default if necessary
    theta.names <- names(init.theta)

    if (is.null(proposal.sd)) {
        proposal.sd <- init.theta/10
    }

    # covmat init
    covmat.proposal <-
        matrix(diag(proposal.sd[theta.names]^2, nrow = length(theta.names)),
               nrow = length(theta.names),
               dimnames = list(theta.names, theta.names))
    covmat.proposal.init <- covmat.proposal

    adapting.size <- FALSE # will be set to TRUE once we start
                           # adapting the size

    adapting.shape <- FALSE # will be set to TRUE once we start
                            # adapting the shape

    # find estimated theta
    theta.estimated.names <- names(which(proposal.sd > 0))

    # evaluate target at theta init
    target.theta.current <- target(theta.current)


    # if return value is a vector, set log.density and trace
    if (class(target.theta.current) == "numeric") {
        suppressWarnings(target.theta.current$log.density <-
            target.theta.current)
        suppressWarnings(target.theta.current$trace <-
            theta.current)
    }

    trace <- data.frame(t(target.theta.current[["trace"]]))

    # acceptance rate
    acceptance.rate <- 0
    if (!is.null(acceptance.rate.window))
    {
        acceptances <- c()
    }


    # scaling factor for covmat size
    scaling.sd  <- 1

    # scaling multiplier
    scaling.multiplier <- 1

    # empirical covariance matrix (0 everywhere initially)
    covmat.empirical <- covmat.proposal
    covmat.empirical[,] <- 0

    # empirical mean vector
    theta.mean <- theta.current

    # if print.info.every is null never print info
    if (is.null(print.info.every)) {
        print.info.every <- n.iterations + 1
    } else
    {
        message("Init: ", printNamedVector(theta.current[theta.estimated.names]),
                ", target: ", target.theta.current[["log.density"]], ", time: ",
                format.POSIXct(Sys.time()))
    }

    start_iteration_time <- Sys.time()

    for (i.iteration in seq_len(n.iterations)) {

        # adaptive step
        if (!is.null(adapt.size.start) && i.iteration >= adapt.size.start &&
            (is.null(adapt.shape.start) ||
                 acceptance.rate*i.iteration < adapt.shape.start)) {
            if (!adapting.size) {
                message("\n---> Start adapting size of covariance matrix")
                adapting.size <- TRUE
            }
            # adapt size of covmat until we get enough accepted jumps
            scaling.multiplier <-
                exp(adapt.size.cooling^(i.iteration-adapt.size.start) *
                        (acceptance.rate - 0.234))
            scaling.sd <- scaling.sd * scaling.multiplier
            scaling.sd <- min(c(scaling.sd,max.scaling.sd))
            # only scale if it doesn't reduce the covariance matrix to 0
            covmat.proposal.new <- scaling.sd^2*covmat.proposal.init
            if (!(any(diag(covmat.proposal.new)[theta.estimated.names] <
                .Machine$double.eps))) {
                covmat.proposal <- covmat.proposal.new
            }

        } else if (!is.null(adapt.shape.start) &&
           acceptance.rate*i.iteration >= adapt.shape.start) {
            if (!adapting.shape) {
                message("\n---> Start adapting shape of covariance matrix")
                # flush.console()
                adapting.shape <- TRUE
            }
            # adapt shape of covmat using optimal scaling factor for multivariate target distributions
            scaling.sd <- 2.38/sqrt(length(theta.estimated.names))

            covmat.proposal <- scaling.sd^2 * covmat.empirical
        }

        # print info
        if (i.iteration %% ceiling(print.info.every) == 0) {
            ## end_iteration_time <- Sys.time()
            state.mcmc <- trace[nrow(trace),]
            message("Iteration: ",i.iteration,"/", n.iterations,
                    ", acceptance rate: ",
                    sprintf("%.3f",acceptance.rate), appendLF=FALSE)
            if (!is.null(adapt.size.start) || !is.null(adapt.shape.start)) {
                message(", scaling.sd: ", sprintf("%.3f", scaling.sd),
                    ", scaling.multiplier: ", sprintf("%.3f", scaling.multiplier),
                    appendLF=FALSE)
            }
            message(", state: ",printNamedVector(state.mcmc),
                    ", target: ", target.theta.current[["log.density"]], 
                    ", time: ", format.POSIXct(Sys.time()))
            ## start_iteration_time <- end_iteration_time
        }

        # propose another parameter set
        if (any(proposal.sd[theta.estimated.names] <
                    .Machine$double.eps)) {
            print(proposal.sd[theta.estimated.names])
            stop("non-positive definite covmat",call.=FALSE)
        }
        if (length(theta.estimated.names) > 0) {
            theta.propose[theta.estimated.names] <-
                as.vector(rmvnorm(1,
                                   mean =
                                       theta.current[theta.estimated.names],
                                   sigma =
                                       covmat.proposal[theta.estimated.names,theta.estimated.names]))
        }

        # evaluate posterior of proposed parameter
        target.theta.propose <- target(theta.propose)
        # if return value is a vector, set log.density and trace
        if (class(target.theta.propose) == "numeric") {
            suppressWarnings(target.theta.propose$log.density <-
                target.theta.propose)
            suppressWarnings(target.theta.propose$trace <-
                theta.propose)
        }

        if (!is.finite(target.theta.propose$log.density)) {
            # if posterior is 0 then do not compute anything else and don't accept
            log.acceptance <- -Inf

        }else{

            # compute Metropolis-Hastings ratio (acceptance probability)
            log.acceptance <- target.theta.propose$log.density -
                target.theta.current$log.density
            ## log.acceptance <- log.acceptance +
            ##     sum(dlnorm(x = theta.current[theta.estimated.names],
            ##                meanlog =
            ##                    theta.propose[theta.estimated.names],
            ##                sdlog =
            ##                    proposal.sd[theta.estimated.names],
            ##                log= TRUE))
            ## log.acceptance <- log.acceptance -
            ##     sum(dlnorm(x = exp(theta.propose[theta.estimated.names]),
            ##                mean = theta.current[theta.estimated.names],
            ##                sdlog =
            ##                    proposal.sd[theta.estimated.names],
            ##                log = TRUE))
            log.acceptance <- log.acceptance +
                dmvnorm(x = theta.current[theta.estimated.names],
                         mean =
                             theta.propose[theta.estimated.names],
                         sigma =
                             covmat.proposal[theta.estimated.names,
                                             theta.estimated.names], 
                         log = TRUE)
            log.acceptance <- log.acceptance -
                dmvnorm(x = theta.propose[theta.estimated.names],
                         mean = theta.current[theta.estimated.names],
                         sigma =
                             covmat.proposal[theta.estimated.names,
                                             theta.estimated.names],
                         log = TRUE)
        }

        if (verbose) {
            message("Propose: ", theta.propose[theta.estimated.names],
                    ", target: ", target.theta.propose[["log.density"]],
                    ", acc prob: ", exp(log.acceptance), ", ",
                    appendLF = FALSE)
        }

        if (is.accepted <- (log(runif (1)) < log.acceptance)) {
            # accept proposed parameter set
            theta.current <- theta.propose
            target.theta.current <- target.theta.propose
            if (verbose) {
                message("accepted")
            }
        } else if (verbose) {
            message("rejected")
        }
        trace <- rbind(trace,c(target.theta.current$trace))

        if (!is.null(acceptance.rate.window))
        {
            acceptances <- c(is.accepted, acceptances)
            if (length(acceptances) > acceptance.rate.window) {
                acceptances <- acceptances[1:acceptance.rate.window]
            }
        }

        # update acceptance rate
        if (i.iteration == 1) {
            acceptance.rate <- is.accepted
        } else {
            if (!is.null(acceptance.rate.window)) {
                acceptance.rate <- mean(acceptances)
            } else {
                acceptance.rate <- acceptance.rate +
                    (is.accepted - acceptance.rate) / i.iteration
            }
        }

        # update empirical covariance matrix
        tmp <- updateCovmat(covmat.empirical, theta.mean,
            theta.current, i.iteration)
        covmat.empirical <- tmp$covmat
        theta.mean <- tmp$theta.mean            
        
    }

    return(list(trace = trace,
                acceptance.rate = acceptance.rate))

}
