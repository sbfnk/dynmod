##' Get parameter values from libbi results file
##'
##' Loop over all one-dimensional entries of the supplied ncdf file
##' (except the first variable, which is time) and store them in a
##' data table
##'
##' @title Get parameter values from libbi results file
##' @param filename name of the libbi results file (in ncdf format)
##' @return data table of parameter values
##' @author Sebastian Funk
##' @import ncdf
libbi.get.params <- function(filename) {

    try_require("data.table")

    nc <- open.ncdf(filename)

    res <- data.table(iteration = 1:nc$dim[[2]]$len)

    for (i in nc$varid2Rindex) {
        if (i > 1 && nc$var[[i]]$ndims == 1) {
            name <- nc$var[[i]]$name
            res[, temp := get.var.ncdf(nc, name)]
            setnames(res, "temp", name)
        }
    }

    close.ncdf(nc)
    res
}

##' Get a trajectory from libbi results file
##'
##' Gets the trajectory with the lowest largest log-likelihood or, if a row
##' is specified, a specific fit from a libbi results file
##' @title Get a fit from libbi results file
##' @param filename name of the libbi results file (in ncdf format)
##' @param row if not NA, the number of the fit one wants to get
##' @return a list of variables, with a trajectory for each
##' @author Sebastian Funk
##' @import ncdf4
##' @export
libbi.get.trajectory <- function(filename, row = NA_integer_) {

    nc <- nc_open(filename)

    varnames <- names(nc$var)

    if (is.na(row)) {
        if ("likelihood" %in% varnames)
        {
            ll <- ncvar_get(nc, "likelihood")
            row <- min(which(ll==max(ll[ll<0])))
        } else
        {
            row <- nc[["dim"]][["np"]][["len"]]
        }
    }

    res <- list()
    for (name in names(nc$var)) {
        if (name == "time") {
            res[[name]] <- ncvar_get(nc, name)
        } else {
            var <- ncvar_get(nc, name)
            if (class(var) == "matrix") {
                res[[name]] <- var[row,]
            } else {
                res[[name]] <- var[row]
            }
        }
    }

    nc_close(nc)
    res
}

##' Get likelihoods from a libbi results file
##'
##' Get likelihoods and associated parameter values from a libbi results file
##' @title Get likelihoods from a libbi results file
##' @param filename name of the libbi results file, in ncdf format
##' @param initial whether initial conditions were also fitted
##' @return a data frame of log-likelihoods and corresnponding parameters
##' @author Sebastian Funk
##' @import ncdf
libbi.get.likelihoods <- function(filename, initial = F) {

    nc <- open.ncdf(filename)

    res <- data.frame(loglikelihood=rep(0, nc$dim$np$len))

    for (i in nc$varid2Rindex) {
        name <- nc$var[[i]]$name
        if (name != "time") {
            if (nc$var[[i]]$ndims == 1) {
                # parameters
                var <- get.var.ncdf(nc, name)
                res[[name]] <- c(var)
            } else if (initial == T & nc$var[[i]]$ndims == 2) {
                # initial conditions
                var <- get.var.ncdf(nc, name)
                res[[name]] <- c(var[, 1])
            }
        }
    }

    res
}

##' Make ncdf file (for libbi) from data
##'
##' @param data dataset to convert
##' @param filename filename
##' @param date.var date variable
##' @param start time at which to start looking at the data
##' @param value.var value variable
##' @param data.var what to call the data
##' @param fill.date whether to fill the data with zeroes
##' @param t0 time 0
##' @import reshape2 ncdf
##' @export
##' @return the outcome of \code{close.ncdf}
##' @author Sebastian Funk
make.ncdf <- function(data, filename, date.var = "rash.date",
                      value.var = "cases", data.var = NA_character_,
                      fill.date = NA_integer_, start.date = NULL,
                      t0 = 0, delta = 1) {
    if (is.na(data.var)) {
        data.var <- value.var
    }

    extra.columns <- colnames(data)[!(colnames(data) %in% c(date.var, value.var))]

    if (length(extra.columns) > 0) {
        dt <-
            data.table(dcast(data,
                             as.formula(paste(date.var,
                                              paste(extra.columns, collapse = "+"),
                                              sep = "~")),
                             value.var = value.var, drop = F))
    } else {
        dt <- data.table(data)
    }

    if (!is.na(fill.date)) {
        time <- seq(dt[1, get(date.var)], dt[nrow(dt), get(date.var)], by = 1)

        if (length(time[!(time %in% dt[,get(date.var)])]) > 0) {
            fill.dt <- cbind(
                data.table(
                    date.var = time[!(time %in% dt[,get(date.var)])]
                ),
                apply(dt[,2:ncol(dt),with=F], 2, function(x) {
                    rep(0, length(time[!(time %in% dt[,get(date.var)])]))
                }
                      )
            )
            setnames(fill.dt, "date.var", date.var)

            dt <- rbind(dt,fill.dt)
        }
    }

    dt <- dt[order(get(date.var))]

    dim1 <- dim.def.ncdf("nr", '', seq(1, nrow(dt)), create_dimvar=F)

    varI <- list()
    vartimes <- list()

    if (length(extra.columns) > 0) {
        for (i in seq(1, ncol(dt) - 1)) {
            varI[[i]] <- var.def.ncdf(paste(data.var, colnames(dt)[i + 1], sep="_"), "",
                                      dim1, -1)
            vartimes[[i]] <- var.def.ncdf(paste("time", colnames(dt)[i + 1], sep="_"), "",
                                          dim1, -1)
        }
    } else {
        varI[[1]] <- var.def.ncdf(data.var, "", dim1, -1)
        vartimes[[1]] <- var.def.ncdf("time", "", dim1, -1)
    }

    nc <- create.ncdf(filename, c(vartimes, varI))
    for (i in seq(1, ncol(dt) - 1)) {
        ncdf.times <- c()
        min.date <- min(dt[, get(date.var)])
        subdata <- na.omit(dt[, c(1, i + 1), with = FALSE])
        if (is.null(start.date))
        {
            ncdf.times <- subdata[, as.integer(((get(date.var) - min.date) / delta) + t0)]
        } else
        {
            ncdf.times <- subdata[, as.integer(((get(date.var) - start.date) / delta) + t0)]
        }
        put.var.ncdf(nc, vartimes[[i]], ncdf.times)

        ncdf.cases <- c(dt[, i + 1, with=F])[[1]]

        put.var.ncdf(nc, varI[[i]], ncdf.cases)
    }

    close.ncdf(nc)
}

summarise.parameters <- function(data, breaks = 50) {

    ret <- data.frame(loglikelihood = data$loglikelihood, logprior = data$logprior)
    for (name in names(data)) {
        if (!(name %in% c("loglikelihood", "logprior"))) {
            ret[[name]] <- strapply(as.character(cut(data[[name]], breaks = breaks)),
                                    "[[:digit:].]+", as.numeric,
                                    simplify = ~colMeans(cbind(...)))
        }
    }
    ret
}

plot.model.file <- function(data, filename, row = 0, ...) {

    res <- get.fit(filename, row)
    plot.model.data.cpx(data, res, ...)

}

plot.model.data.measles <- function(data, sims, id.column = "rash.date",
                                    deterministic = F, N = 73640000, n = 3,
                                    smooth = F, all = F, model = T) {

    days <- seq(data[[id.column]][1], data[[id.column]][nrow(data)], by = 1)

    model.data <- data.table()
    if ("agegroup" %in% names(data)) {
        model.data <- dcast(
            data, ... ~ agegroup, value.var = "cases", drop = F
        )
    } else {
        model.data <- data.table(data)
    }

    model.data <- data.table(rbind.fill(model.data,
                                        data.frame(
                                            rash.date =
                                                days[!(days %in% data[[id.column]])]
                                        )
                                        ))
    model.data <- model.data[order(model.data[[id.column]])]

    ## model.data[days %in% dt[[column]],
    ##            cases := dt[,cases]]
    start <- length(days) - length(sims$time) + 1
    numberZ <- length(grep("^Z", names(sims)))
    if (numberZ > 1) {
        sim.data <- data.table(model.data[[id.column]])
        for (i in 1:numberZ) {
            sim.data[start:length(days),
                     paste("Z", i, sep="") := sims[[paste("Z", i, sep="")]],
                     with = F]
        }
        if ("agegroup" %in% names(data)) {
            setnames(sim.data, names(sim.data), names(model.data))
            model.data[, data := "cases"]
            sim.data[, data := "model"]
            model.data <- rbind(model.data, sim.data)
            model.data <- data.table(melt(model.data, id.vars=c("rash.date", "data")))
            setnames(model.data, "variable", "agegroup")
        } else {
            model.data[, model := rowSums(sim.data[, 2:ncol(sim.data), with = F])]
            model.data <- model.data[, list(rash.date, cases, model)]
            model.data <- data.table(melt(model.data, id.vars=c(id.column)))
            setnames(model.data, "variable", "data")
        }
    } else {
        model.data[start:length(days), model := sims$Z]
        model.data <- data.table(melt(model.data, id.vars=c(id.column)))
        setnames(model.data, "variable", "data")
    }

    if (all == F) {
        model.data <- model.data[rash.date >= days[start]]
    }

    model.data[data=="model" & !is.na(value),
               value := sapply(value, function(x) { rbinom(1, round(x), sims$rep) })]

    model.data[value == 0, value := NA]
    setnames(model.data, "value", "cases")

    if (model == F) {
        model.data <- model.data[data == "cases"]
    }

    p <- ggplot(model.data, aes(x=rash.date, y=cases, color=data))+
        geom_point()+
            ## scale_y_continuous(limits = c(0,
            ##                   max(model.data$cases, na.rm = T) * 1.1))+
            scale_x_date("Date of rash")+
                theme_bw(20)+
                    theme(legend.position="none")


    if ("agegroup" %in% names(data)) {
        p <- p + facet_grid(agegroup ~ ., scale = "free")
    }
    if (deterministic) {

        parameters <- c(beta = mean(res$b), gamma = res$g, N = res$n * N, n=n)
        times <- seq(start, nrow(model.data))
        init <- initSIR(c(S = res$R / res$R0 * res$n * N, I = 1,
                          R = (1 - res$prop_sus) * N), n=n)

        out <- data.table(ode(y = init, times = times, func = sir, parms = parameters))

        model.data[out$time, incidence := c(0, res$rep * (out$incidence[2:nrow(out)] -
                                                              out$incidence[1:(nrow(out)-1)]))]

        p <- p+ geom_line(data=model.data, aes(y=incidence,
                          color="deterministic"))
    }

    p <- p + scale_colour_brewer(palette="Set1", name="")
    if (smooth) {
        p <- p + geom_smooth()
    }
    p
}


plot.likelihoods <- function(data, breaks = 50, summarise = F,  smooth = F, density = F) {

    if (summarise) {
        library(gsubfn)
        df <- summarise.parameters(data, breaks = breaks)
    }

    df <- data.table(melt(data, id.vars=c("loglikelihood", "logprior")))

    if (density) {
        p <- ggplot(df[loglikelihood > -Inf],
                    aes(x = value, color = variable)) +
                        geom_density()
    } else {
        p <- ggplot(df[loglikelihood > -Inf],
                    aes(x = value, y = loglikelihood, color = variable))+
                        geom_point()
    }

    p <-  p + facet_wrap( ~ variable, scales="free_x")+
        theme(legend.position="none")
    if (smooth) {
        p <- p + geom_smooth()
    }

    p

}

cpx.ncdf <- function(transient) {

    make.ncdf(weekly.incidence, "cpx0.nc", date.var = "date", t0 = 0,
              value.var = "incidence", data.var = "Inc")
    make.ncdf(weekly.incidence, "cpx.nc", date.var = "date",
              t0 = (transient %/% 52) * 52 , value.var = "incidence",
              data.var = "Inc")
}

get_erlang_groups <- function(erlang_list)
{
    stage_list <- lapply(erlang_list, function(x) seq_len(x[["n"]]))
    return(apply(expand.grid(stage_list, stringsAsFactors = FALSE), 1,
                 function(x)
                 {
                     paste(x, collapse = ",")
                 }))
}

##' Write a bi file for an SEIR model
##'
##' @param model list of things -- see code
##' @param groups groups (if there are any)
##' @param delta delta for the transition block
##' @param sum.agegroups whether to sum up age groups for incidence
##' @param ... parameters to pass to generate_bi
##' @return nothing -- a libbi file is written
##' @export
##' @author Sebastian Funk
generate_bi_seir <- function(name = "SEIR", model, desc = NULL, print = FALSE)
{
    bi_params <- model[["params"]]
    names(bi_params) <- sub("\\.", "_", names(bi_params))
    model[["params"]] <- bi_params

    bi_options <- model[["options"]]

    groups <- NULL
    group_string <- ""
    if (!is.null(model[["groups"]]))
    {
        groups <- model[["groups"]]
        group.table <- data.table(expand.grid(model[["groups"]]))
        group_string <- paste(names(groups), collapse = ",")
    }

    ## checks
    for (param_name in names(bi_params))
    {
        param <- bi_params[[param_name]]
        if (!is.data.frame(param))
        {
            if (is.null(param[["string"]]))
            {
                stop("Parameter ", param_name, " needs a string")
            }
            if (!is.null(param[["proposal"]]))
            {
                proposal <- param[["proposal"]]
                if (is.null(proposal[["string"]]))
                {
                    stop("Parameter ", param_name, " proposal needs a string")
                }
            }
        }
    }

    # set up possible compartments,  and what can happen to them
    reported_compartments <- c("I", "I_H", "H", "R", "H_R")

    # container for erlang transitions
    transitions = list()

    # force of infection
    foi_string <- "foi"
    if (nchar(group_string) > 0)
    {
        foi_string <- paste(foi_string, "[", group_string, "]", sep = "")
    }

    transitions[["infection"]] <-
        list(from = "S", to = "I",
             rate = foi_string)
    symptom_transition <- "infection"

    # set up erlang transitions based on given parameters
    if ("gamma" %in% names(bi_params))
    {
        transitions[["recovery"]] <-
            list(from = "I", to = "R",
                 rate = "gamma")
    }
    if ("rho" %in% names(bi_params))
    {
        transitions[["symptoms"]] <-
            list(from = "E", to = "I",
                 rate = "rho")
        transitions[["infection"]][["to"]] <- "E"
        symptom_transition <- "symptoms"
    }
    if ("kappa" %in% names(bi_params))
    {
        transitions[["ready.hosp"]] <-
            list(from = "I", to = "I_H",
                 rate = "kappa")
        symptom_transition <- "ready.hosp"
    }
    if ("tau" %in% names(bi_params) && "K" %in% names(bi_params))
    {
        if ("kappa" %in% names(bi_params))
        {
            transitions[["hospitalisation"]] <-
                list(from = "I_H", to = "H",
                     rate = "tau * (1 - H/K)")
        } else
        {
            transitions[["hospitalisation"]] <-
                list(from = "I", to = "H",
                     rate = "tau * (1 - H/K)")
        }
        if ("gamma" %in% names(bi_params))
        {
            transitions[["hosp.recovery"]] <-
                list(from = "H", to = "H_R",
                     rate = "gamma")
        }
    } else if ("tau" %in% names(bi_params) || "K" %in% names(bi_params))
      {
          warning("tau and K both have to be given to include hospital transmission")
      }

    if ("alpha" %in% names(bi_params))
    {
        nalt <- length(transitions[[symptom_transition]][["alternative"]])
        if (nalt == 0)
        {
            transitions[[symptom_transition]][["alternative"]] <- list()
        }
        transitions[[symptom_transition]][["alternative"]][[nalt + 1]] <-
            list(to = "R", factor = "alpha")
    }

    # gather compartments
    compartments <- unique(c(unlist(sapply(transitions, function(x)
    {
        c(x[["from"]], x[["to"]])
    })), unlist(sapply(transitions, function(x)
    {
        comp <- c()
        for (alt_idx in seq_along(x[["alternative"]]))
        {
            comp <- c(comp, x[["alternative"]][[alt_idx]][["to"]])
        }
        comp
    }))))
    # filter compartments by what is actually there
    reported_compartments <- intersect(reported_compartments, compartments)

    if ("epsilon" %in% names(bi_params))
    {
        for (rep_comp in reported_compartments)
        {
            rep_name <- paste("unreported", rep_comp, sep = "_")
            transitions[[rep_name]] <-
                list(from = rep_comp, rate = "epsilon")
            nalt <- length(transitions[[rep_name]][["alternative"]])
            if (nalt == 0)
            {
                transitions[[rep_name]][["alternative"]] <- list()
            }
            to_inc <- "Z"
            if (length(grep("^H", rep_comp)) > 0)
            {
                to_inc <- "Z_H"
            }
            transitions[[rep_name]][["alternative"]][[nalt + 1]] <-
                list(to = to_inc, factor = transitions[[rep_name]][["factor"]])
        }
    }
    else
    {
        # incidence
        nalt <- length(transitions[[symptom_transition]][["alternative"]])
        if (nalt == 0)
        {
            transitions[[symptom_transition]][["alternative"]] <- list()
        }
        transitions[[symptom_transition]][["alternative"]][[nalt + 1]] <-
            list(to = "Z", factor = transitions[[symptom_transition]][["factor"]])
    }

    # update compartments
    compartments <- unique(c(unlist(sapply(transitions, function(x)
    {
        c(x[["from"]], x[["to"]])
    })), unlist(sapply(transitions, function(x)
    {
        comp <- c()
        for (alt_idx in seq_along(x[["alternative"]]))
        {
            comp <- c(comp, x[["alternative"]][[alt_idx]][["to"]])
        }
        comp
    }))))

    erlang_groups <- list()
    dimensions <- list()

    # fill erlang_groups
    for (transition.name in names(transitions))
    {
        transition <- transitions[[transition.name]]
        from <- transition[["from"]]
        rate <- transition[["rate"]]
        if (is.null(bi_params[[rate]]))
        {
            n <- 1
        } else
        {
            n <- bi_params[[rate]][["erlang"]]
            if (is.null(n))
            {
                n <- 1
            }
        }

        ngroups <- length(erlang_groups[[from]])

        to_class <- list(transition[["rate"]])
        names(to_class) <- transition[["to"]]
        if (!is.null(transition[["alternative"]]))
        {
            for (alt_idx in seq_along(transition[["alternative"]]))
            {
                alt.transition <- transition[["alternative"]][[alt_idx]]
                factor <- alt.transition[["factor"]]
                if (is.null(factor))
                {
                    alternative.rate <- to_class[[1]]

                } else
                {
                    alternative.rate <- paste(factor, "*", to_class[[1]])
                    other.factor <- paste("(1 - ", factor, ")", sep = "")
                    for (to_idx in seq_along(to_class))
                    {
                        to_class[[to_idx]] <-
                            paste(other.factor, "*", to_class[[to_idx]])
                    }
                }
                to_class[[length(to_class) + 1]] <- alternative.rate
                names(to_class)[length(to_class)] <- alt.transition[["to"]]
            }
        }
        erlang_groups[[from]][[ngroups + 1]] <-
            list(name = transition.name,
                 to = to_class,
                 rate = rate,
                 n = n)
        dimensions[[paste(rate, "erlang", sep = "_")]] <- n
    }

    for (group_name in names(groups))
    {
        dimensions[[group_name]] <- length(groups[[group_name]])
    }

    ## open file
    filename <- paste(name, "bi", sep = ".")
    f <- file(description = filename, open = "w")
    if (!is.null(desc))
    {
        bi_write(c("/**", paste(" *", desc), " */", ""),
                 con = f, print = print)
    }
    bi_write(c(paste("model", name, "{"), ""), con = f, print = print)

    if (length(dimensions) > 0)
    {
        for (dim_name in names(dimensions))
        {
            if (dimensions[[dim_name]] > 1)
            {
                bi_write(paste("  const",
                               paste("n_", dim_name, sep = ""),
                               "=", dimensions[[dim_name]]),
                         con = f, print = print)
            }
        }
        bi_write("", con = f, print = print)
        for (dim_name in names(dimensions))
        {
            if (dimensions[[dim_name]] > 1)
            {
                bi_write(paste("  dim",
                               paste(dim_name, "(n_", dim_name, ")", sep = "")),
                         con = f, print = print)
            }
        }
        bi_write("", con = f, print = print)
    }

    param_dims <- list()
    for (param_name in names(bi_params))
    {
        param <- bi_params[[param_name]]
        param_dims[[param_name]] <- c()
        if (is.data.frame(param))
        {
            for (group_name in sub("\\..*$", "", colnames(param)))
            {
                if (group_name %in% names(groups))
                {
                    param_dims[[param_name]] <-
                        c(param_dims[[param_name]], group_name)
                }
            }
        } else if (!is.null(param[["group"]]))
          {
              param_dims[[param_name]] <- param[["group"]]
          }

        param_string <- param_name
        if (length(param_dims[[param_name]]) > 0)
        {
            param_string <-
                paste(param_string,
                      "[", paste(param_dims[[param_name]], collapse = ","), "]",
                      sep = "")
        }

        bi_write(paste("  param", param_string))
    }
    bi_write("", con = f, print = print)

    bi_write(paste("  state", foi_string))

    state_dims <- list()
    for (state in compartments)
    {
        if (length(groups) > 0)
        {
            state_dims[[state]] <- names(groups)
        }
        for (group in erlang_groups[[state]])
        {
            if (group[["n"]] > 1)
            {
                state_dims[[state]] <- c(state_dims[[state]],
                                         paste(group[["rate"]], "erlang", sep = "_"))
            }
        }
        state_string <- state
        if (length(state_dims[[state]] > 0))
        {
            state_string <-
                paste(state_string,
                      "[", paste(state_dims[[state]], collapse = ","), "]",
                      sep = "")
        }
        bi_write(paste("  state", state_string), con = f, print = print)
    }
    bi_write("", con = f, print = print)

    inc_string <- paste("Inc")
    if (!is.null(model[["observations"]]))
    {
        inc_string <-
            paste(inc_string, "[",
                  paste(model[["observations"]], collapse = ","),
                  "]", sep = "")
    }
    bi_write(paste("  obs", inc_string), con = f, print = print)

    transition_string <- ""
    if (!is.null(bi_options[["transition"]]))
    {
        transition_opts <- bi_options[["transition"]]
        opts <- sapply(names(transition_opts), function(option)
        {
            paste(option, "=", transition_opts[[option]])
        })
        transition_string <- paste("(", paste(opts, collapse = ", "), ") ",
                                   sep = "")
    }
    bi_write(c("", paste("  sub transition ", transition_string, "{",
                         sep = "")), con = f, print = print)

    bi_write("", con = f, print = print)

    for (inc_comp in grep("^Z", compartments, value = TRUE))
    {
        state_str <- paste("   ", inc_comp)
        if (length(state_dims[[inc_comp]]) > 0)
        {
            state_str <- paste(state_str, "[",
                               paste(state_dims[[inc_comp]], collapse = ","),
                               "]", sep = "")
        }
        bi_write(paste(state_str, " <- ", 0, sep = ""),
                 con = f, print = print)
    }

    # force of infection
    foi_assign <- "R0"
    if (length(param_dims[["R0"]]) > 0)
    {
        foi_assign <- paste(foi_assign, "[",
                            paste(param_dims[["R0"]], collapse = ","),
                            "]", sep = "")
    }
    foi_assign <- paste(foi_assign, "*", "gamma", "*")
    if (!is.null(bi_params[["mixing"]]))
    {
        mixing_names <- colnames(bi_params[["mixing"]])
        mixing_dims <- list()

        for (group_name in mixing_names)
        {
            stripped_name <- sub("\\.[0-9]*$", "", group_name)
            if (stripped_name %in% names(groups))
            {
                if (length(grep("\\.2$", group_name)) > 0)
                {
                    mixing_dims[[length(mixing_dims) + 1]] <-
                        seq_along(groups[[stripped_name]]) - 1
                    names(mixing_dims)[length(mixing_dims)] <- group_name
                } else
                {
                    mixing_dims[[length(mixing_dims) + 1]] <-
                        stripped_name
                    names(mixing_dims)[length(mixing_dims)] <- group_name
                }
            }
        }
        sum_fields <- expand.grid(mixing_dims, stringsAsFactors = FALSE)
        sum_str <- c()
        if (nrow(sum_fields) > 0)
        {
            for (i in seq_len(nrow(sum_fields)))
            {
                I_comps <- list()
                for (I_dim in state_dims[["I"]])
                {
                    if (paste(I_dim, 2, sep = ".") %in% colnames(sum_fields))
                    {
                        I_comps <-
                            c(I_comps, sum_fields[i, paste(I_dim, 2, sep = ".")])
                    } else if (I_dim %in% colnames(sum_fields))
                      {
                          I_comps <-
                              c(I_comps, sum_fields[i, I_dim])
                      } else
                      {
                          I_comps <-
                              c(I_comps, seq_along(dimensions[[I_dim]]) - 1)
                      }
                }
                sum_str[i] <-
                    paste("mix[", paste(sum_fields[i, ], collapse = ","),
                          "] *", sep = "")
                I_fields <- expand.grid(I_comps)
                I_str <- "I"
                if (nrow(I_fields) > 0)
                {
                    I_str <- paste(I_str, "[",
                                   apply(I_fields, 1, function(x)
                                   {
                                       paste(x, collapse = ",")
                                   }), "]", sep = "")
                }
                if (nrow(I_fields) > 1)
                {
                    sum_str[i] <- paste(sum_str[i], "(",
                                        paste(I_str, collapse = " + "),
                                        ")", sep = "")
                } else
                {
                    sum_str[i] <- paste(sum_str[i], I_str)
                }
            }
        } else {
            sum_str <- paste("mix", "*")
            I_comps <- list()
            for (I_dim in state_dims[["I"]])
            {
                I_comps <- c(I_comps, seq_along(dimensions[[I_dim]]) - 1)
            }
            I_fields <- expand.grid(I_comps)
            I_str <- "I"
            if (nrow(I_fields) > 0)
            {
                I_str <- paste(I_str, "[",
                               apply(I_fields, 1, function(x)
                               {
                                   paste(x, collapse = ",")
                               }), "]", sep = "")
            }
            if (length(I_comps) > 1)
            {
                sum_str[i] <- paste(sum_str[i], "(",
                                    paste(I_str, collapse = " + "),
                                    ")", sep = "")
            } else
            {
                sum_str[i] <- paste(sum_str[i], I_str)
            }
        }
        if (length(sum_str) > 1)
        {
            foi_assign <- paste(foi_assign, " (",
                                paste(sum_str, collapse = " + "),
                                ")", sep = "")
        } else
        {
            foi_assign <- paste(foi_assign, sum_str)
        }
    } else
    {
        if (length(state_dims[["I"]]) > 0)
        {
            foi_elements <- c()
            for (group in get_erlang_groups(erlang_groups[["I"]]))
            {
                foi_elements <- c(foi_elements,
                                  paste("I[", group, "]", sep = ""))
            }
            foi_assign <- paste(foi_assign, " (",
                                paste(foi_elements, collapse = " + "),
                                ")", sep = "")
        } else
        {
            foi_assign <- paste(foi_assign, "I")
        }
    }
    bi_write(paste("   ", foi_string, "<-", foi_assign), con = f, print = print)

    bi_write(c("", "    ode {", ""), con = f, print = print)

    # gather transitions and assign to state derivatives
    derivs <- list()
    # froms
    for (state in intersect(compartments, names(erlang_groups)))
    {
        for (trans_idx in seq_along(erlang_groups[[state]]))
        {
            transition <- erlang_groups[[state]][[trans_idx]]
            erlang_id <- paste(transition[["rate"]], "erlang", sep = "_")
            prefix <- ""
            if (transition[["n"]] > 1)
            {
                prefix <- paste(transition[["n"]], "*")
            }
            state_name <- state
            if (length(state_dims[[state]]) > 0)
            {
                state_name <- paste(state_name, "[",
                                    paste(state_dims[[state]], collapse = ","),
                                    "]", sep = "")
            }
            if (nchar(prefix) > 0)
            {
                rate <- paste(prefix, transition[["rate"]])
            } else
            {
                rate <- transition[["rate"]]
            }
            derivs[[state_name]] <-
                paste(derivs[[state_name]], "-", rate, "*", state_name)

            if (length(state_dims[[state]]) > length(groups))
            {
                add_dims <- state_dims[[state]]
                for (dim_idx in seq_along(add_dims))
                {
                    if (add_dims[[dim_idx]] == erlang_id)
                    {
                        add_dims[[dim_idx]] <- paste(erlang_id, "-", 1)
                    }
                }
                state_name_from <- paste(state, "[",
                                         paste(add_dims, collapse = ","),
                                         "]", sep = "")
                derivs[[state_name]] <-
                    paste(derivs[[state_name]], "+", rate, "*", state_name_from)
            }
        }
    }

    for (state in names(erlang_groups))
    {
        for (trans_idx in seq_along(erlang_groups[[state]]))
        {
            transition <- erlang_groups[[state]][[trans_idx]]
            erlang_id <- paste(transition[["rate"]], "erlang", sep = "_")
            prefix <- ""
            if (transition[["n"]] > 1)
            {
                prefix <- paste(transition[["n"]], "*")
            }
            state_name <- state
            if (length(state_dims[[state]]) > 0)
            {
                state_name <- paste(state_name, "[",
                                    paste(state_dims[[state]], collapse = ","),
                                    "]", sep = "")
            }
            for (to_state in names(transition[["to"]]))
            {
                individual_flag <- FALSE
                to_dims <- state_dims[[to_state]]
                for (dim_idx in seq_along(to_dims))
                {
                    if (!(to_dims[dim_idx] %in% state_dims[[state]]))
                    {
                        to_dims[dim_idx] <- 0
                        individual_flag <- TRUE
                    }
                }
                to_name <- to_state
                if (length(to_dims) > 0)
                {
                    to_name <- paste(to_name, "[",
                                     paste(to_dims, collapse = ","),
                                     "]", sep = "")
                }
                rate <- paste(transition[["to"]][[to_state]], "*", rate)
                if (nchar(prefix) > 0)
                {
                    rate <- paste(prefix, transition[["rate"]])
                } else
                {
                    rate <- transition[["rate"]]
                }
                derivs[[to_name]] <-
                    paste(derivs[[to_name]], "+", rate, "*", state_name)
                if (individual_flag)
                {
                    for (trans_idx in seq_along(erlang_groups[[to_state]]))
                    {
                        transition <- erlang_groups[[to_state]][[trans_idx]]
                        erlang_id <- paste(transition[["rate"]], "erlang", sep = "_")
                        prefix <- ""
                        if (transition[["n"]] > 1)
                        {
                            prefix <- paste(transition[["n"]], "*")
                        }
                        if (nchar(prefix) > 0)
                        {
                            rate <- paste(prefix, transition[["rate"]])
                        } else
                        {
                            rate <- transition[["rate"]]
                        }
                        derivs[[to_name]] <-
                            paste(derivs[[to_name]], "-", rate, "*", to_name)
                    }
                }
            }
        }
    }

    for (deriv_idx in seq_along(derivs))
    {
        bi_write(paste("      d", names(derivs)[deriv_idx], "/dt = ",
                       derivs[[deriv_idx]], sep = ""),
                 con = f, print = print)
    }

    bi_write(c("    }", "  }", "", "  sub parameter {"), con = f, print = print)

    for (param in bi_params)
    {
    }
    bi_write(c("  }", "", "  sub initial {"), con = f, print = print)
    

    ##     # initial conditions
    ## for (grp_idx in seq_along(us_groups))
    ## {
    ##     group <- us_groups[grp_idx]
    ##     age_dist <- eigen(mixing[[grp_idx]])$vectors[, 1] /
    ##         sum(eigen(mixing[[grp_idx]])$vectors[, 1])
    ##     for (age_idx in seq_along(us_agegroups[[grp_idx]]))
    ##     {
    ##         agegroup <- us_agegroups[[grp_idx]][[age_idx]]
    ##         state_name_S <- paste("S", group, agegroup, sep = "")
    ##         init_S <- paste("N", group, agegroup, sep = "")
    ##         bi_states[[state_name_S]] <- list(string = "")

    ##         init_I_common <- ""
    ##         if ("initI" %in% names(bi_params))
    ##         {
    ##             init_I_common <- paste(init_I_common, "initI", "*",
    ##                                    age_dist[age_idx])
    ##         } else if (all(paste("initI", groups, sep = "_") %in%
    ##                        names(bi_params)))
    ##           {
    ##               init_I_common <- paste(init_I_common, "initI", group, " * ",
    ##                                      age_dist[age_idx], sep = "")
    ##           } else
    ##           {
    ##               init_I_common <-
    ##                   paste(init_I_common, "initI", group, agegroup, sep = "")
    ##           }

    ##         e_groups <- get_erlang_groups(erlang_groups[["E"]])
    ##         for (erlanggroup in e_groups)
    ##         {
    ##             init_S <- paste(init_S, " - E", erlanggroup, sep = "")
    ##             state_name_E <- paste("E", group, agegroup, erlanggroup, sep = "")
    ##             init_E <- paste("gamma", "/", "rho", "*")
    ##             if (length(n_compartments[["E"]]) > 1)
    ##             {
    ##                 init_E <- paste(1, "/",  length(n_compartments[["E"]]), "*", init_E)
    ##             }
    ##             i_comp <- c()
    ##             init_E <- paste(init_E, init_I_common)
    ##             bi_states[[state_name_E]] <- list(string = init_E)
    ##         }

    ##         i_groups <- get_erlang_groups(erlang_groups[["I"]])
    ##         for (erlanggroup in i_groups)
    ##         {
    ##             state_name_I <-
    ##                 paste("I", group, agegroup, erlanggroup, sep = "")
    ##             init_S <- paste(init_S, "-", state_name_I)
    ##             if (length(grep("_u(_|$)", erlanggroup)) > 0)
    ##             {
    ##                 init_I <- paste("gamma", "/", "(epsilon", "+", "gamma)", "*", init_I_common)
    ##                 if (n_compartments[["I"]] > 1)
    ##                 {
    ##                     init_I <- paste(1, "/",  n_compartments[["I"]], "*", init_I)
    ##                 }
    ##             }
    ##             else
    ##             {
    ##                 init_I <- init_I_common
    ##             }

    ##             bi_states[[state_name_I]] <- list(string = init_I)
    ##         }

    ##         ih_groups <- get_erlang_groups(erlang_groups[["I_H"]])
    ##         for (erlanggroup in ih_groups)
    ##         {
    ##             state_name_IH <-
    ##                 paste("I_H", group, agegroup, ih_groups, sep = "")
    ##             init_S <- paste(init_S, "-", state_name_IH)
    ##             bi_states[[state_name_IH]] <- list(string = 0)
    ##         }
    ##         h_groups <- get_erlang_groups(erlang_groups[["H"]])
    ##         for (erlanggroup in h_groups)
    ##         {
    ##             state_name_H <-
    ##                 paste("H", group, agegroup, erlanggroup, sep = "")
    ##             init_S <- paste(init_S, "-", state_name_H)
    ##             bi_states[[state_name_H]] <- list(string = 0)
    ##         }
    ##         r_groups <- get_erlang_groups(erlang_groups[["R"]])
    ##         for (erlanggroup in r_groups)
    ##         {
    ##             state_name_R <-
    ##                 paste("R", group, agegroup, erlanggroup, sep = "")
    ##             init_S <- paste(init_S, "-", state_name_R)
    ##             bi_states[[state_name_R]] <- list(string = 0)
    ##         }
    ##         hr_groups <- get_erlang_groups(erlang_groups[["H_R"]])
    ##         for (erlanggroup in hr_groups)
    ##         {
    ##             state_name_HR <-
    ##                 paste("H_R", group, agegroup,
    ##                       erlanggroup, sep = "")
    ##             init_S <- paste(init_S, "-", state_name_HR)
    ##             bi_states[[state_name_HR]] <- list(string = 0)
    ##         }
    ##         bi_states[[state_name_S]][["string"]] <- init_S
    ##         state_name_Z <- paste("Z", group, agegroup, sep = "")
    ##         bi_states[[state_name_Z]] <- list(string = 0)
    ##         bi_states[[state_name_Z]][["reset"]] <- list(string = 0)
    ##         if ("H" %in% compartments)
    ##         {
    ##             state_name_ZH <- paste("Z_H", group, agegroup, sep = "")
    ##             bi_states[[state_name_ZH]] <- list(string = 0)
    ##             bi_states[[state_name_ZH]][["reset"]] <- list(string = 0)
    ##         }
    ##     }
    ## }

    ## if (sum.agegroups)
    ## {
    ##     for (grp_idx in seq_along(us_groups))
    ##     {
    ##         group <- us_groups[grp_idx]
    ##         obs_name <- paste("Inc", group, sep = "")
    ##         inc_term <- c()
    ##         inc_h_term <- c()
    ##         for (agegroup in us_agegroups[[grp_idx]])
    ##         {
    ##             inc_term <- c(inc_term,
    ##                           paste("Z", group, agegroup, sep = ""))
    ##             inc_h_term <- c(inc_h_term,
    ##                             paste("Z_H", group, agegroup, sep = ""))
    ##         }
    ##         obs_args <-
    ##             c(paste(paste(inc_h_term, collapse = " + "), " + ",
    ##                     "rep * (", paste(inc_term, collapse = " + "), ")",
    ##                     sep = ""),
    ##               paste("sqrt(rep * (1 - rep) * (",
    ##                     paste(inc_term, collapse = " + "),
    ##                     "))", sep = ""),
    ##               "lower = 0")
    ##         bi_obs[[obs_name]] <- list(string = "truncated_normal",
    ##                                    argument = obs_args)
    ##     }
    ## } else
    ## {
    ##     for (grp_idx in seq_along(us_groups))
    ##     {
    ##         group <- us_groups[grp_idx]
    ##         for (agegroup in us_agegroups[[grp_idx]])
    ##         {
    ##             obs_name <- paste("Inc", group, agegroup, sep = "")
    ##             obs_args <-
    ##                 c(paste("Z_H", group, agegroup,
    ##                         " + rep * Z", group, agegroup,
    ##                         sep = ""),
    ##                   paste("sqrt(rep * (1 - rep) * Z",
    ##                         group, agegroup,
    ##                         ")", sep = ""),
    ##                   "lower = 0")
    ##             bi_obs[[obs_name]] <- list(string = "truncated_normal",
    ##                                        arguments = obs_args)
    ##         }
    ##     }
    ## }

    ## ## transitions
    ## for (grp_idx in seq_along(us_groups))
    ## {
    ##     group <- us_groups[grp_idx]
    ##     beta_name <- paste("beta", group, sep = "")
    ##     if (!(R0_name %in% names(bi_params)))
    ##     {
    ##         beta_name <- "beta"
    ##     }
    ##     for (age_idx in seq_along(us_agegroups[[grp_idx]]))
    ##     {
    ##         agegroup <- us_agegroups[[grp_idx]][[age_idx]]
    ##         inf_comp_all <- c()
    ##         i_groups <- get_erlang_groups(erlang_groups[["I"]])
    ##         for (age_idx2 in seq_along(us_agegroups[[grp_idx]]))
    ##         {
    ##             agegroup2 <- us_agegroups[[grp_idx]][[age_idx2]]
    ##             N_param_name <-
    ##                 paste("N", us_groups[grp_idx],
    ##                       us_agegroups[[grp_idx]][age_idx2], sep = "")
    ##             if (length(us_agegroups[[grp_idx]]) > 1)
    ##             {
    ##                 mix_term <- paste("mix", us_groups[grp_idx], "_",
    ##                                   age_idx, age_idx2, " * ",
    ##                                   sep = "")
    ##             } else
    ##             {
    ##                 mix_term <- ""
    ##             }
    ##             inf_sub_comp <- c()
    ##             for (erlanggroup in i_groups)
    ##             {
    ##                 inf_sub_comp <-
    ##                     c(inf_sub_comp,
    ##                       paste("I", group, agegroup2,
    ##                             erlanggroup, sep = ""))
    ##             }
    ##             inf_comp_all <-
    ##                 c(inf_comp_all,
    ##                   paste(mix_term, "(",
    ##                         paste(inf_sub_comp, collapse = " + "), ") / ",
    ##                         N_param_name, sep = ""))
    ##         }
    ##         if (length(i_groups) > 1)
    ##         {
    ##             inf_term <-
    ##                 paste(beta_name, " * ", "(",
    ##                       paste(inf_comp_all,
    ##                             collapse = " + "), ")",
    ##                       sep = "")
    ##         } else {
    ##             inf_term <-
    ##                 paste(beta_name, "*", paste("I", group, agegroup, sep = ""))
    ##         }
    ##         ## infection
    ##         state_name_S <- paste("S", group, agegroup, sep = "")
    ##         to_class <- list()
    ##         state_name_inf <- ""
    ##         inf_compartment <- infection_compartments[1]
    ##         inf_groups <- get_erlang_groups(erlang_groups[[inf_compartment]])
    ##         state_name_inf <-
    ##             paste(inf_compartment, group, agegroup,
    ##                   inf_groups[1], sep = "")
    ##         if (!is.null(erlang_groups[["E"]]))
    ##         {
    ##             if ("alpha" %in% names(bi_params))
    ##             {
    ##                 to_class[[state_name_inf]] <- "(1 - alpha)"
    ##                 E_groups <- sapply(erlang_groups[["E"]], function(x) { x[["short_name"]] })
    ##                 R_groups <- sapply(erlang_groups[["R"]], function(x) { x[["short_name"]] })
    ##                 matching_class_R <- get_erlang_groups(erlang_groups[["R"]])
    ##                 common_groups_R <- intersect(E_groups, R_groups)
    ##                 for (grp in common_groups_R)
    ##                 {
    ##                     from_grp <- sub(paste("^.*(_", grp, "_[0-9]+)(_.*$|$)", sep = ""), "\\1", from_group)
    ##                     matching_class_R <- grep(from_grp, matching_class_R, value = TRUE)
    ##                 }
    ##                 reported <- grep("_rep(_|$)", matching_class_R)
    ##                 if (length(reported) > 0)
    ##                 {
    ##                     ## asymptomatic,  won't be recorded
    ##                     matching_class_R <- grep("_rep(_|$)", matching_class_R, value = TRUE)
    ##                 }
    ##                 state_name_R <- paste("R", group, agegroup, matching_class_R[1], sep = "")
    ##                 to_class[[state_name_R]] <- "alpha"
    ##             } else {
    ##                 to_class[[state_name_inf]] <- 1
    ##             }
    ##         } else
    ##         {
    ##             to_class[[state_name_inf]] <- 1
    ##         }
    ##         transition_name <- paste("infection", group, agegroup, sep = "")
    ##         bi_transitions[[transition_name]] <-
    ##             list(from = state_name_S,
    ##                  to = to_class,
    ##                  rate = paste(inf_term, "*", state_name_S))

    ##         ## hospitalisation
    ##         ih_groups <- get_erlang_groups(erlang_groups[["I_H"]])
    ##         h_groups <- get_erlang_groups(erlang_groups[["H"]])
    ##         for (erlanggroup in ih_groups)
    ##         {
    ##             transition_name <- paste("hospitalisation", group,
    ##                                      agegroup, erlanggroup, sep = "")
    ##             state_name_I <-
    ##                 paste("I_H", group, agegroup, erlanggroup, sep = "")
    ##             state_name_H <-
    ##                 paste("H", group, agegroup, h_groups[1], sep = "")
    ##             rate <- paste("tau", "*", "(K > 0", "?",
    ##                           paste("(1 - ", state_name_H, "/K)", sep = ""),
    ##                           ": 0)", "*", state_name_I)
    ##             bi_transitions[[transition_name]] <-
    ##                 list(from = state_name_I,
    ##                      to = state_name_H,
    ##                      rate = rate)
    ##         }

    ##         ## moving through classes
    ##         for (state in names(erlang_groups))
    ##         {
    ##             erlang <- erlang_groups[[state]]
    ##             for (trans_idx in seq_along(erlang))
    ##             {
    ##                 transition <- erlang[[trans_idx]]
    ##                 if (length(transition[["stages"]]) > 0)
    ##                 {
    ##                     term <- transition[["rate"]]
    ##                     stages <- transition[["stages"]]
    ##                     n_stages <- transition[["n"]]
    ##                     if (n_stages > 1) {
    ##                         term <- paste(n_stages, "*", term)
    ##                     }

    ##                     all_groups <- get_erlang_groups(erlang_groups[[state]])
    ##                     for (grp in seq_len(length(stages) - 1))
    ##                     {
    ##                         from_groups <- list()
    ##                         for (from_group in grep(stages[grp], all_groups, value = T))
    ##                         {
    ##                             transition_name <-
    ##                                 paste(transition[["name"]], group, agegroup,
    ##                                       sub(paste("^", state, "_", sep = ""),
    ##                                           paste(transition[["name"]], "_", sep = ""),
    ##                                           from_group),
    ##                                       sep = "")
    ##                             to_class_main <- paste(state, group, agegroup,
    ##                                                    sub(stages[grp], stages[grp + 1], from_group),
    ##                                                    sep = "")
    ##                             to_class <- list()
    ##                             to_class[[to_class_main]] <- 1
    ##                             if (length(grep("_rep(_|$)", to_class_main)) > 0)
    ##                             {
    ##                                 # we're entering a report compartment
    ##                                 state_name_Z <- ""
    ##                                 if (length(grep("^H", to_class_main)) > 0)
    ##                                 {
    ##                                     state_name_Z <- paste("Z_H", group, agegroup, sep = "")
    ##                                 } else
    ##                                 {
    ##                                     state_name_Z <- paste("Z", group, agegroup, sep = "")
    ##                                 }
    ##                                 to_class[[state_name_Z]] <- 1
    ##                             }
    ##                             from_class <- paste(state, group, agegroup, from_group, sep = "")
    ##                             bi_transitions[[transition_name]] <-
    ##                                 list(from = from_class,
    ##                                      to = to_class,
    ##                                      rate = paste(term, "*", from_class))
    ##                         }
    ##                     }

    ##                     ## transitioning to the next compartment
    ##                     if (!is.null(transition[["final"]]))
    ##                     {
    ##                         final_transition <- transition[["final"]]

    ##                         other_transitions <- erlang
    ##                         other_transitions[[trans_idx]] <- NULL
    ##                         to_transitions <- erlang_groups[[final_transition[["to"]]]]

    ##                         other_groups <- sapply(other_transitions, function(x) { x[["short_name"]] })
    ##                         to_groups <- sapply(to_transitions, function(x) { x[["short_name"]] })

    ##                         common_groups <- intersect(other_groups, to_groups)

    ##                         for (from_group in grep(stages[length(stages)], all_groups, value = T))
    ##                         {
    ##                             matching_class <- get_erlang_groups(to_transitions)
    ##                             transition_name <- paste(final_transition[["name"]], group, agegroup, sep = "")
    ##                             for (grp in common_groups)
    ##                             {
    ##                                 from_group <-
    ##                                     sub(paste("^.*(_", grp, "_[0-9]*)([_a-zA-Z]|$)", sep = ""), "\\1", from_group)
    ##                                 matching_class <- grep(from_group, matching_class, value = TRUE)
    ##                                 transition_name <-
    ##                                     paste(transition_name, from_group, sep = "")
    ##                             }
    ##                             to_class_main <- paste(final_transition[["to"]], group, agegroup,
    ##                                                    matching_class[1], sep = "")
    ##                             to_class <- list()

    ##                             if (state == "E")
    ##                             {
    ##                                 if ("alpha" %in% names(bi_params))
    ##                                 {
    ##                                     to_class[[to_class_main]] <- "(1 - alpha)"
    ##                                     R_groups <- sapply(erlang_groups[["R"]], function(x) { x[["short_name"]] })
    ##                                     matching_class_R <- get_erlang_groups(erlang_groups[["R"]])
    ##                                     common_groups_R <- intersect(other_groups, R_groups)
    ##                                     for (grp in common_groups_R)
    ##                                     {
    ##                                         from_grp <- sub(paste("^.*(_", grp, "_[0-9]+)(_.*$|$)", sep = ""), "\\1", from_group)
    ##                                         matching_class_R <- grep(from_grp, matching_class_R, value = TRUE)
    ##                                     }
    ##                                     reported <- grep("_rep(_|$)", matching_class_R)
    ##                                     if (length(reported) > 0)
    ##                                     {
    ##                                         ## asymptomatic,  won't be recorded
    ##                                         matching_class_R <- grep("_rep(_|$)", matching_class_R, value = TRUE)
    ##                                     }
    ##                                     state_name_R <- paste("R", group, agegroup, matching_class_R[1], sep = "")
    ##                                     to_class[[state_name_R]] <- "alpha"
    ##                                 } else
    ##                                 {
    ##                                     to_class[[to_class_main]] <- 1
    ##                                 }
    ##                                 if (length(grep("_rep(_|$)", i_groups)) == 0)
    ##                                 {
    ##                                     state_name_Z <- paste("Z", group, agegroup, sep = "")
    ##                                     to_class[[state_name_Z]] <- to_class[[to_class_main]]
    ##                                 }
    ##                             } else
    ##                             {
    ##                                 to_class[[to_class_main]] <- 1
    ##                             }

    ##                             from_class <- paste(state, group, agegroup, from_group, sep = "")
    ##                             bi_transitions[[transition_name]] <-
    ##                                 list(from = from_class,
    ##                                      to = to_class,
    ##                                      rate = paste(term, "*", from_class))
    ##                         }
    ##                     }
    ##                 }
    ##             }
    ##         }
    ##     }
    ## }
    ## bi_model <- list(params = bi_params,
    ##                  inlines = bi_inlines,
    ##                  states = bi_states,
    ##                  transitions = bi_transitions,
    ##                  observations = bi_obs,
    ##                  options = bi_options)
    ## return(generate_bi(model = bi_model, ...))

    close(f)
}

##' Generates a bi string from a string and a vector
##'
##' Internal function used bi generate_bi
##' @param l list with a 'string' (character) and 'argument- (vector or list)
##' @return bi string
##' @author Sebastian Funk
bi_string <- function(l)
{
    str <- l[["string"]]
    vec <- l[["argument"]]
    if (is.null(vec))
    {
        return(paste("<-", str))
    } else
    {
        argument_str <-
            paste(str, "(", paste(vec, collapse = ", "), ")", sep = "")
        return(paste("~", argument_str))
    }
}

##' Writes a vector to a file connector (and to the screen, if descired)
##'
##' Internal function used bi generate_bi
##' @param text text to print
##' @param con file connector (stdout by default)
##' @param print whether to also print the text to screen (ignored if
##' \code{con} is set to stdout)
##' @return result of the write to the file connector
##' @author Sebastian Funk
bi_write <- function(text, con = stdout(), print = FALSE, ...)
{
    if (con == stdout())
    {
        print <- FALSE
    }
    if (print == TRUE)
    {
        writeLines(text, ...)
    }
    return(writeLines(text, con, ...))
}

##' Write a bi file for an SEIR model
##'
##' @param name name of the model (and file)
##' @param model list of things -- see code
##' @param desc description of the model
##' @return nothing -- a libbi file is written (unless \code{print} is
##' set to TRUE, in which case the output is also printed to screen
##' @export
##' @author Sebastian Funk
generate_bi <- function(name = "SEIR", model, desc = NULL, print = FALSE)
{
    bi_params <- model[["params"]]
    bi_inlines <- model[["inlines"]]
    bi_states <- model[["states"]]
    bi_transitions <- model[["transitions"]]
    bi_obs <- model[["observations"]]
    bi_options <- model[["options"]]

    names(bi_params) <- sub(" ", "_", names(bi_params))
    names(bi_states) <- sub(" ", "_", names(bi_states))
    names(bi_transitions) <- sub(" ", "_", names(bi_transitions))
    names(bi_obs) <- sub(" ", "_", names(bi_obs))
    ## checks
    for (param_name in names(bi_params))
    {
        param <- bi_params[[param_name]]
        if (is.null(param[["string"]]))
        {
            stop("Parameter ", param_name, " needs a string")
        }
        if (!is.null(param[["proposal"]]))
        {
            proposal <- param[["proposal"]]
            if (is.null(proposal[["string"]]))
            {
                stop("Parameter ", param_name, " proposal needs a string")
            }
        }
    }
    for (state_name in names(bi_states))
    {
        state <- bi_states[[state_name]]
        if (is.null(state[["string"]]))
        {
            stop("State ", state_name, " needs a string")
        }
        if (!is.null(state[["reset"]]))
        {
            reset <- state[["reset"]]
            if (is.null(reset[["string"]]))
            {
                stop("State ", state_name, " reset needs a string")
            }
        }
    }
    for (transition_name in names(bi_transitions))
    {
        transition <- bi_transitions[[transition_name]]
        if (is.null(transition[["from"]]) &&
            is.null(transition[["to"]]))
        {
            stop("Transition ", transition_name, " needs a from or to")
        }
        if (is.null(transition[["rate"]]))
        {
            stop("Transition ", transition_name, " needs a rate")
        }
    }
    for (obs_name in names(bi_obs))
    {
        obs <- bi_obs[[obs_name]]
        if (is.null(obs[["string"]]))
        {
            stop("Observation ", obs_name, " needs a string")
        }
    }

    ## open file
    filename <- paste(name, "bi", sep = ".")
    f <- file(description = filename, open = "w")
    if (!is.null(desc))
    {
        bi_write(c("/**", paste(" *", desc), " */", ""),
                 con = f, print = print)
    }
    bi_write(c(paste("model", name, "{"), ""), con = f, print = print)

    ## all parameters that are not lists are considered to be fixed
    fixed_params <- unlist(sapply(names(bi_params), function(param)
    {
        if (is.null(bi_params[[param]][["argument"]]))
        {
            param
        }
    }, USE.NAMES = FALSE))

    for (param_name in fixed_params)
    {
        param <- bi_params[[param_name]]
        bi_write(paste("  const", param_name, "=", param[["string"]],
                       sep = " "), con = f, print = print)
        bi_params[[param_name]] <- NULL
    }
    bi_write("", con = f, print = print)

    for (param_name in names(bi_params))
    {
        bi_write(paste("  param", param_name), con = f, print = print)
    }
    bi_write("", con = f, print = print)

    noise <- NULL
    if (!is.null(bi_options[["noise"]]))
    {
        noise <- bi_options[["noise"]]
        if (noise == "scaled")
        {
            for (transition in names(bi_transitions))
            {
                bi_write(paste("  noise",
                               paste("noise", transition, sep = "_")),
                         con = f, print = print)
            }
        }
        bi_write("", con = f, print = print)
    }

    for (state_name in names(bi_states))
    {
        bi_write(paste("  state", state_name), con = f, print = print)
    }
    bi_write("", con = f, print = print)

    for (obs_name in names(bi_obs))
    {
        bi_write(paste("  obs", obs_name), con = f, print = print)
    }
    bi_write("", con = f, print = print)

    for (inline_name in names(bi_inlines))
    {
        inline <- bi_inlines[[inline_name]]
        bi_write(paste("  inline", inline_name, "=", inline),
                 con = f, print = print)
    }

    transition_string <- ""
    if (!is.null(bi_options[["transition"]]))
    {
        transition_opts <- bi_options[["transition"]]
        opts <- sapply(names(transition_opts), function(option)
        {
            paste(option, "=", transition_opts[[option]])
        })
        transition_string <- paste("(", paste(opts, collapse = ", "), ") ",
                                   sep = "")
    }
    bi_write(c("", paste("  sub transition ", transition_string, "{",
                         sep = "")), con = f, print = print)
    bi_write("", con = f, print = print)

    if (!is.null(noise) && noise == "scaled")
    {
        for (transition_name in names(bi_transitions))
        {
            bi_write(paste("    ",
                           paste("noise", transition_name, sep = "_"),
                           " ~ wiener()", sep = ""),
                     con = f, print = print)
        }
        bi_write("", con = f, print = print)
    }

    reset_any <- FALSE
    for (state_name in names(bi_states))
    {
        state <- bi_states[[state_name]]
        if (!is.null(state[["reset"]]))
        {
            reset <- state[["reset"]]
            bi_write(paste("    ", state_name, " ", bi_string(reset), sep = ""),
                     con = f, print = print)
            reset_any <- TRUE
        }
    }
    if (reset_any)
    {
        bi_write("", con = f, print = print)
    }

    bi_write(c("    ode {", ""), con = f, print = print)

    ## transitions
    for (state_name in names(bi_states))
    {
        bi_states[[state_name]][["deriv"]] <- ""
    }

    for (transition_name in names(bi_transitions))
    {
        transition <- bi_transitions[[transition_name]]
        sign <- list(from = "-", to = "+")
        for (direction in c("from", "to"))
        {
            if (!is.null(transition[[direction]]))
            {
                dir_states <- transition[[direction]]
                if (!(class(dir_states) == "list"))
                {
                    name <- dir_states
                    dir_states <- list()
                    dir_states[[name]] <- 1
                }
                for (state in names(dir_states))
                {
                    if (is.null(bi_states[[state]]))
                    {
                        warning("Dir state ", state, " does not exit")
                    } else
                    {
                        rate_str <- ""
                        if (dir_states[[state]] != 1)
                        {
                            rate_str = paste(dir_states[[state]], " * ", sep = "")
                        }
                        rate_str <- paste(rate_str, transition[["rate"]], sep = "")
                        bi_states[[state]][["deriv"]] <-
                            paste(bi_states[[state]][["deriv"]], sign[[direction]], rate_str)
                        if (!is.null(noise) && noise == "scaled")
                        {
                            noise_term <- paste("sqrt(", rate_str, ") * ",
                                                "noise_", transition_name, sep = "")
                            bi_states[[state]][["deriv"]] <-
                                paste(bi_states[[state]][["deriv"]],
                                      "+", noise_term)
                        }
                    }
                }
            }
        }
    }

    for (state_name in names(bi_states))
    {
        state <- bi_states[[state_name]]
        eq <- paste("      d", state_name, "/dt =", state[["deriv"]], sep = "")
        bi_write(eq, con = f, print = print)
    }

    bi_write(c("    }", "  }", "", "  sub parameter {"), con = f, print = print)

    for (param_name in names(bi_params))
    {
        param <- bi_params[[param_name]]
        bi_write(paste("    ", param_name, " ", bi_string(param),
                       sep = ""), con = f, print = print)
    }
    bi_write(c("  }", "", "  sub proposal_parameter {"), con = f, print = print)

    for (param in names(bi_params))
    {
        proposal <- bi_params[[param]][["proposal"]]
        if (!is.null(proposal))
        {
            bi_write(paste("    ", param, " ", bi_string(proposal),
                           sep = ""), con = f, print = print)
        }
    }
    bi_write(c("  }", "", "  sub initial {"), con = f, print = print)

    for (state in names(bi_states))
    {
        bi_write(paste("    ", state, " ", bi_string(bi_states[[state]]),
                       sep = ""), con = f, print = print)
    }

    bi_write(c("  }", "", "  sub observation {"), con = f, print = print)

    for (obs_name in names(bi_obs))
    {
        obs <- bi_obs[[obs_name]]
        bi_write(paste("    ", obs_name, " ", bi_string(obs),
                       sep = ""), con = f, print = print)
    }
    bi_write(c("  }", "}"), con = f, print = print)

    close(f)
}
