##' Generate a trajectory of brownian motion
##'
##' Generate a trajectory of brownian motion, given an initial value,
##' as well as drift and scale parameters
##' @param parameters A named vector or list that includes the initial
##' value ("init"), drift ("drift") and scale ("scale") as named
##' parameters
##' @param n Number of steps to compute
##' @return Brownian trajectory
##' @author Sebastian Funk
brownian <- function(parameters, n) {

    res <- parameters[["brownian.init"]]
    steps <- floor(n) - 1

    if (steps > 0) {
        std.br <- cumsum(rnorm(steps))
        res <- c(res, sapply(seq_along(std.br), function(x) {
            parameters[["brownian.init"]] +
                parameters[["brownian.drift"]] * x +
                    parameters[["brownian.scale"]] * std.br[x]
        }))
    }

    return(res)
}

##' Generate a piecewise linear trajectory
##'
##' @param mixing.values the mixing values, including initial and final values (so this has to have length two more than \code{change.times}.
##' @param change.times the time points at which the linear trajectory changes
##' @param n the total number of time points to generate
##' @return the piecewise linear trajectory
##' @author Sebastian Funk
piecewise.linear <- function(mixing.values, change.times, n)
{
    if (length(mixing.values) == 0)
    {
        mixing <- c()
    } else if (length(mixing.values) == 1)
      {
          # if only one value is given,  we don't interpolate
          mixing <- rep(mixing.values,n)
      } else
      {
          # interpolate
          times <- c(1, change.times, n)
          mixing <- approx(times, mixing.values, seq_len(n))$y
      }

    return(mixing)
}

##' Estimate vaccine-induced immunity in the population from reported MMR coverage
##'
##' Estimate vaccine-induced immunity in the population from reported MMR coverage
##' @param vaccine.efficacy Efficacy of one dose of the vaccien
##' @param newborn.immunity Immunity of newborns for hte first year of
##' their lives. If not given, this is extracted from serology data
##' @param years the years for which to calculate immunity
##' @return A data table with vaccine-induced immunity in all age groups
##' @import Hmisc reshape2
##' @author Sebastian Funk
estimate.immunity.mmr <- function(vaccine.efficacy = 0.9,
                                  newborn.immunity = NULL,
                                  years = NULL,
                                  age.limits = NULL) {

    data(vaccine_ew)
    data(ms_sero)
    data(ms_ew_age)
    data(cpx_ew_age)

    if (is.null(newborn.immunity)) {
        newborn.immunity <-
            ms.sero[country == "UK" & exact.age == 1, mean(non.negative)]
    }

    if (is.null(age.limits))
    {
        age.limits <- unique(ms.ew.age[, lower.age.limit])
    }

    ## allocate vaccination data by birth cohort
    vaccine.ew.yob <- data.table(vaccine.ew)
    ## MCV1 coverage is measured at 2 years of age
    vaccine.ew.yob <- vaccine.ew.yob[, year := year - 2]
    setnames(vaccine.ew.yob, "year", "year.of.birth")

    max.age <- max(vaccine.ew[, year]) - min(vaccine.ew.yob[, year.of.birth])

    ## MCV1 coverage at 5 years of age
    vaccine.ew.yob <-
        vaccine.ew.yob[,
                       MCV15 := c(vaccine.ew.yob[4:nrow(vaccine.ew.yob), MCV15],
                             rep(NA, 3))]
    ## MCV2 coverage at 5 years of age
    vaccine.ew.yob <-
        vaccine.ew.yob[,
                       MCV2 := c(vaccine.ew.yob[4:nrow(vaccine.ew.yob), MCV2],
                            rep(NA, 3))]
    ## estimate increase in MCV1 coverage between 2 and 5 years of age
    vaccine.ew.yob <- vaccine.ew.yob[, annual.mcv1 := (MCV15 - MCV1) / 3]

    ## calculate vaccine coverage by birth cohort
    ## first year: maternal antibodies -- read from ESEN serology
    vaccine.ew.yob <-
        vaccine.ew.yob[, year.1 := newborn.immunity]

    ## fill in missing data with newborn immunity

    if (min(vaccine.ew.yob[, year.of.birth]) > min(cpx.ew.age[, year])) {
        for (year in seq(min(cpx.ew.age[, year]) - 1,
                         min(vaccine.ew.yob[, year.of.birth]) - 1)) {
            vaccine.ew.yob <-
                rbind(data.table(year.of.birth = year,
                                 MCV1 = 0,
                                 MCV15 = NA,
                                 MCV2 = NA,
                                 annual.mcv1 = NA,
                                 year.1 = newborn.immunity),
                      vaccine.ew.yob)
        }
    }
    ## second year: year 2 MCV1 coverage
    vaccine.ew.yob <-
        vaccine.ew.yob[, year.2 := MCV1 * vaccine.efficacy / 100]
    ## year 3 and 4: estimated from year 2 MCV1 coverage and year 2 MCV1 coverage
    vaccine.ew.yob <-
        vaccine.ew.yob[, year.3 := (MCV1 + annual.mcv1) * vaccine.efficacy / 100]
    vaccine.ew.yob <-
        vaccine.ew.yob[, year.4 := (MCV1 + 2 * annual.mcv1) * vaccine.efficacy / 100]
    ## year 5: measured MCV1 and MCV2 coverage (where avilable)
    vaccine.ew.yob <-
        vaccine.ew.yob[, year.5 := (MCV15 * vaccine.efficacy / 100 +
                                        MCV2 * (1 - vaccine.efficacy) / 100 *
                                            vaccine.efficacy * MCV15/100)]
    mvac.ew.yob <- melt(vaccine.ew.yob, id.vars = c("year.of.birth"),
                        measure.vars = grep("^year.[0-9]+", names(vaccine.ew.yob),
                        value = T))

    mvac.ew.yob <- mvac.ew.yob[, year := as.numeric(gsub("year\\.", "", variable))]
    mvac.ew.yob <- mvac.ew.yob[!is.na(value)]

    min.year <- min(vaccine.ew.yob[, year.of.birth]) + 1
    max.year <- max(vaccine.ew[, year])
    vaccination.cohorts <-
        c(sapply(seq(min.year, max.year), function (y) {
            sapply(seq_len(max.age), function(x) {
                ifelse(nrow(mvac.ew.yob[year.of.birth == y - x]) == 0, 0,
                       ifelse(length(mvac.ew.yob[year.of.birth == y - x & year == x, value]) == 0,
                              mvac.ew.yob[year.of.birth == y - x & year == max(mvac.ew.yob[year.of.birth == y - x, year]), value],
                              mvac.ew.yob[year.of.birth == y - x & year == x, value]))
            })
        }))

    immunity.estimate <-
        data.table(year = rep(seq(min.year, max.year), each = max.age),
                   age = rep(seq(1, max.age), times = max.year - min.year + 1),
                   vaccination = vaccination.cohorts)
    immunity.estimate <- immunity.estimate[, lower.age.limit := age - 1]
    immunity.estimate[, lower.age.limit :=
                            reduce.agegroups(lower.age.limit,
                                             limits = age.limits)]
    immunity.estimate <- immunity.estimate[, list(vaccination = mean(vaccination)),
                                           by = list(year, lower.age.limit)]

    ## bring back into wide format
    immunity.estimate.table <-
        data.table(dcast(immunity.estimate, year  ~ lower.age.limit,
                         value.var = "vaccination"))

    ## fill years before
    if (!is.null(years) && min(years) < min(immunity.estimate.table[, year])) {
        for (year in rev(years[years < min(immunity.estimate.table[, year])])) {
            new.year <- data.table(year,
                                   immunity.estimate.table[1, 2, with = F],
                                   t(rep(0, ncol(immunity.estimate.table) - 2)))
            setnames(new.year, names(new.year), names(immunity.estimate.table))
            immunity.estimate.table <- rbind(new.year, immunity.estimate.table)
        }
    }

    ## if no years were selected, select all
    if (is.null(years)) {
        years <- immunity.estimate.table[, year]
    }

    ## fill missing ages
    for (col in setdiff(age.limits, colnames(immunity.estimate.table)))
    {
        immunity.estimate.table[, paste(col) := as.numeric(0)]
    }

    return(immunity.estimate.table[year %in% years])

}

##' Estimate vaccine uptake by age group in the population from reported MMR coverage
##'
##' @param vaccine.efficacy Efficacy of one dose of the vaccien
##' @param years the years for which to calculate immunity
##' @return A data table with vaccine-induced immunity in all age groups
##' @author Sebastian Funk
##' @export
estimate.uptake.mmr <- function() {

    data(vaccine_ew)
    data(ms_sero)
    data(ms_ew_age)

    ## allocate vaccination data by birth cohort
    vaccine.ew.yob <- data.table(vaccine.ew)
    ## MCV1 coverage is measured at 2 years of age
    vaccine.ew.yob <- vaccine.ew.yob[, year := year - 2]
    setnames(vaccine.ew.yob, "year", "year.of.birth")

    max.age <- max(vaccine.ew[, year]) - min(vaccine.ew.yob[, year.of.birth])

    ## MCV1 coverage at 5 years of age
    vaccine.ew.yob <-
        vaccine.ew.yob[,
                       MCV15 := c(vaccine.ew.yob[4:nrow(vaccine.ew.yob), MCV15],
                             rep(NA, 3))]
    ## MCV2 coverage at 5 years of age
    vaccine.ew.yob <-
        vaccine.ew.yob[,
                       MCV2 := c(vaccine.ew.yob[4:nrow(vaccine.ew.yob), MCV2],
                            rep(NA, 3))]
    ## we take the maximum of the two (one approach) (or we could
    ## estimate the increase from the data)
    vaccine.ew.yob  <-
        vaccine.ew.yob[, list("0" = max(MCV1, MCV15)),
                       by = list(year.of.birth, MCV2, MCV15, MCV1)]

    increase.estimate <-
        vaccine.ew.yob[!is.na(MCV15), sum(MCV15) / sum(MCV1)]

    vaccine.ew.yob[is.na(MCV15), "0" := MCV1 * increase.estimate]
    vaccine.ew.yob <-
        vaccine.ew.yob[, list("0" = get("0"), "5" = 0, "10" = 0, "15" = 0, "25" = 0),
                       by = list(year.of.birth)]

    min.diff <- min(vaccine.ew.yob[, year.of.birth]) - min(ms.ew.age[, year])
    if (min.diff > 0)
    {
        for (i in seq_len(min.diff))
        {
            vaccine.ew.yob <- rbind(vaccine.ew.yob[1, ], vaccine.ew.yob)
            vaccine.ew.yob[1, 2:ncol(vaccine.ew.yob)] <- 0
            vaccine.ew.yob[1, 1] <- vaccine.ew.yob[1, 1] - 1
        }
    }

    vaccine.ew.yob[, year.of.birth := NULL]

    return(as.matrix(vaccine.ew.yob / 100))

}

##' Chickenpox age distribution of cases from the Wallinga (2004) method
##'
##' @param ... parameters to be passed to fit.polymod
##' @import reshape2
##' @return the results of fit.polymod
##' @author Sebastian Funk
cpx.age.wallinga <- function(...) {

    data(polymod)
    data(cpx_ew_age)
    data(pop_ew_age)

    ## only use british data (for now)
    participants <- polymod$participants[country == "GB"]
    contacts <- polymod$contacts[country == "GB"]

    ages <- pop.ew.age[year == 2006]
    ages[, lower.age.limit :=
               reduce.agegroups(lower.age.limit,
                                unique(cpx.ew.age[, lower.age.limit]))]
    ages <- ages[, list(population = sum(population)), by = lower.age.limit]

    cpx.annual <- cpx.ew.age[, list(number.cases = sum(abs.incidence),
                                    rel.incidence = sum(rel.incidence)),
                             by = list(year, lower.age.limit)]

    annual.cases <- cpx.annual[, list(annual.cases = sum(number.cases)),
                               by = list(year)]
    setkey(cpx.annual, year)
    setkey(annual.cases, year)

    cpx.annual <- merge(cpx.annual, annual.cases)
    cpx.annual <- cpx.annual[, proportion.cases := number.cases / annual.cases]

    cpx.annual.tab <- dcast(cpx.annual, year ~ lower.age.limit,
                            value.var="rel.incidence")

    fit.polymod(participants, contacts, ages, cpx.annual.tab[,-1],
                sample = T, var = "y", ...)
}

##' Measles age distribution of cases from the Wallinga (2004) method
##'
##' @param ... parameters to be passed to fit.polymod
##' @import reshape2
##' @return the results of fit.polymod
##' @author Sebastian Funk
ms.age.structure.wallinga <- function(...) {

    data(polymod)
    data(ms_ew_age)
    data(pop_ew_age)

    ## only use british data (for now)
    participants <- polymod$participants[country == "GB"]
    contacts <- polymod$contacts[country == "GB"]

    ages <- pop.ew.age[year == 2006]
    ages[, lower.age.limit :=
               reduce.agegroups(lower.age.limit,
                                unique(ms.ew.age[, lower.age.limit]))]
    ages <- ages[, list(population = sum(population)), by = lower.age.limit]

    ms.target <- ms.ew.age[, list(number.cases = sum(abs.incidence),
                                  rel.incidence = sum(rel.incidence)),
                           by = list(year, lower.age.limit)]

    annual.cases <- ms.target[, list(annual.cases = sum(number.cases)),
                              by = list(year)]
    setkey(ms.target, year)
    setkey(annual.cases, year)

    ms.target <- merge(ms.target, annual.cases)
    ms.target <- ms.target[, proportion.cases := number.cases / annual.cases]

    ms.target.tab <- dcast(ms.target, year ~ lower.age.limit,
                           value.var="proportion.cases")

    ## ms.target.tab <- dcast(ms.target, year ~ lower.age.limit,
    ##                         value.var="rel.incidence")

    yby <- fit.polymod(participants, contacts, ages, ms.target.tab[,-1],
                       sample = T, yearbyyear = T, var = "y",
                       ...)
    one <- fit.polymod(participants, contacts, ages, ms.target.tab[,-1],
                       sample = T, var = "y", ...)

    list(yby = yby, one = one)
}

##' Calculate previous states of a population, using a demographic model
##'
##' @param nyears number of years to backcalculate
##' @param lower.age.limits lower limits of age groups
##' @return named vector of states
##' @author Sebastian Funk
##' @import deSolve
backcalc.population <- function(nyears = 1, lower.age.limits = NULL) {

    data(pop_ew_age)
    data(births_ew)
    data(deaths_ew)

    if (is.null(lower.age.limits)) {
        lower.age.limits <- unique(pop.ew.age[, lower.age.limit])
    }

    ## reduce population to given age groups
    pop.ew.age[, lower.age.limit :=
                     reduce.agegroups(lower.age.limit, lower.age.limits)]
    pop.ew.age <- pop.ew.age[, list(population = sum(population)),
                             by = list(year, lower.age.limit)]

    min.year <- pop.ew.age[, min(year)]
    first.pop.state <- pop.ew.age[year == min.year, population]
    agewidths <- diff(lower.age.limits)

    pop.times <- seq(0, nyears)
    ew.pop.parameters <- list(births = rep(births.ew[year == 1971, births], nyears + 1),
                           deaths = rep(deaths.ew[year == 1971, deaths], nyears + 1),
                           agewidths = agewidths,
                           reverse = T)
    trajectory <- data.table(ode(y = first.pop.state, times = pop.times,
                                 func = "demo_derivs", parms = ew.pop.parameters,
                                 hini = 1, dllname = "",
                                 initfunc = "demo_initmod",
                                 nout = 1, outnames = "N"))
    setnames(trajectory, 2:ncol(trajectory), c(lower.age.limits, "total"))

    round(trajectory[-1, 2:(ncol(trajectory) - 1), with = F])
}

##' Get population parameters for England & Wales
##'
##' Fill parameters for births, deaths and ageing, with England &
##' Wales population data
##' @param age.limits Lower limits of the age groups
##' @param year.limits first and last year
##' @param interpolate whether to ineterpolate between years
##' @param df whether to return time-dependent variables as data frame
##' @return list of parameters
##' @export
##' @author Sebastian Funk
pop.parameters <- function(age.limits, births, deaths, year.limits = NULL, interpolate = FALSE, df = FALSE) {

    data(births_ew)
    data(deaths_ew)

    if (missing(births)) {
        births <- births.ew
    }
    if (missing(deaths)) {
        deaths <- deaths.ew
    }


    if (is.null(year.limits))
    {
        year.limits <- intersect(births[, year], deaths[, year])
    }

    if (interpolate)
    {
        year.limits[which.min(year.limits)] <- min(year.limits) - 1
        year.limits[which.max(year.limits)] <- max(year.limits) + 1
    }

    parameters <- list()
    parameters[["births"]] <-
        births[year >= min(year.limits) & year <= max(year.limits), births]
    parameters[["deaths"]] <-
        deaths[year >= min(year.limits) & year <= max(year.limits), deaths]

    ## ageing rates depend on the widths of age bands
    parameters[["ageing"]] <- 1 / diff(unique(age.limits))

    years <-
        unique(intersect(births[year >= min(year.limits) & year <= max(year.limits), year],
                         deaths[year >= min(year.limits) & year <= max(year.limits), year]))
    if (df)
    {
        parameters[["births"]] <-
            data.table(year = years,
                       value = births[year >= min(year.limits) & year <= max(year.limits), births])
        parameters[["deaths"]] <-
            data.table(year = years,
                       value = deaths[year >= min(year.limits) & year <= max(year.limits), deaths])
        parameters[["ageing"]] <-
          data.table(age = limits.to.agegroups(age.limits),
                     value = c(parameters[["ageing"]], 0))

        setkey(parameters[["births"]], year)
        setkey(parameters[["deaths"]], year)
        setkey(parameters[["ageing"]], age)
    } else
    {
        parameters[["years"]] <- years
        ## we're looking at annual data, births and deaths have an entry
        ## per time point
        parameters[["paramStep"]] <- 1
    }

    return(parameters)

}

## legacy function(deprecated)
ew.pop.parameters <- function(...) {
    pop.parameters(...)
}

##' Generate initial conditions for measles or chickenpox from given proportion
##' of susceptibles
##'
##' @param parameters Model parameters
##' @param infection "measles" or "chickenpox"
##' @return Named vector of initial conditions
##' @author Sebastian Funk
##' @export
gen.init <- function(parameters, infection = c("measles", "chickenpox")) {

    infection <- match.arg(infection)

    init <- c()

    if (infection == "measles")
    {
        data(ms_ew_age)

        infectious.period <- 1 / parameters[["gamma"]]
        latent.period <- 1 / parameters[["rho"]]

        lower.limits <- agegroups.to.limits(colnames(parameters[["mixing"]]))

        ms.ew.age[, lower.age.limit :=
                        reduce.agegroups(lower.age.limit, lower.limits)]
        ms.ew.age <- ms.ew.age[, list(abs.incidence = sum(abs.incidence),
                                      population = sum(population)),
                               by = list(year, lower.age.limit)]

        min.year <- ms.ew.age[, min(year)]
        ms.init.infected <-
            ms.ew.age[year == min.year,
                      list(abs.prevalence = sum(abs.incidence) /
                               parameters[["underreporting"]] *
                               (latent.period + infectious.period),
                           population = mean(population)),
                      by = list(lower.age.limit)]

        init.susc <- unname(unlist(parameters[grep("^init\\.susc\\.[0-9]+$",
                                                   names(parameters))]))

        init.I <- round(ms.init.infected[, abs.prevalence] * infectious.period /
                        (infectious.period + latent.period))
        init.E <- round(ms.init.infected[, abs.prevalence] * latent.period /
                        (infectious.period + latent.period))
        init.S <- pmin(round(init.susc * ms.init.infected[, population]),
                       ms.init.infected[, population] - init.I - init.E)
        init.R <- ms.init.infected[, population] - init.S - init.I - init.E
        init.Z <- rep(0, nrow(ms.init.infected))

        ## make sure all is >0
        init.S[which(init.R < 0)] <-
            init.S[which(init.R < 0)] + init.R[which(init.R < 0)]
        init.R[which(init.R < 0)] <-  0

        if ("maternal.immunity" %in% names(parameters) &&
            parameters[["maternal.immunity"]] > 0)
        {
            init.B <- init.S[1] * parameters[["maternal.immunity"]]
            init.S[1] <- init.S[1] * (1 - parameters[["maternal.immunity"]])
            init <- c(init.B, init.S, init.E, init.I, init.R, init.Z)
            names(init) <- c("B1",
                             paste("S", seq_len(nrow(ms.init.infected)), sep = ""),
                             paste("E", seq_len(nrow(ms.init.infected)), sep = ""),
                             paste("I", seq_len(nrow(ms.init.infected)), sep = ""),
                             paste("R", seq_len(nrow(ms.init.infected)), sep = ""),
                             paste("Z", seq_len(nrow(ms.init.infected)), sep = ""))
        } else
        {
            init <- c(init.S, init.E, init.I, init.R, init.Z)
            names(init) <- c(paste("S", seq_len(nrow(ms.init.infected)), sep = ""),
                             paste("E", seq_len(nrow(ms.init.infected)), sep = ""),
                             paste("I", seq_len(nrow(ms.init.infected)), sep = ""),
                             paste("R", seq_len(nrow(ms.init.infected)), sep = ""),
                             paste("Z", seq_len(nrow(ms.init.infected)), sep = ""))
        }
    } else if (infection == "chickenpox")
    {
        data(cpx_ew_age)

        infectious.period <- 1 / parameters[["gamma"]]
        latent.period <- 1 / parameters[["rho"]]

        lower.limits <- agegroups.to.limits(colnames(parameters[["mixing"]]))

        cpx.ew.age[, lower.age.limit := reduce.agegroups(lower.age.limit,
                                                         lower.limits)]
        min.year <- cpx.ew.age[, min(year)]
        cpx.init.infected <-
            cpx.ew.age[year == min.year,
                       list(abs.prevalence = sum(abs.incidence) /
                                parameters[["underreporting"]] *
                                (latent.period + infectious.period),
                            population = mean(population)),
                       by = list(lower.age.limit)]

        init.susc <- unname(unlist(parameters[grep("^init\\.susc\\.[0-9]+$",
                                                   names(parameters))]))
        init.S <-
            round(init.susc * cpx.init.infected[, population])

        if (any(grep("^init\\.inf\\.[0-9]+$", names(parameters))))
        {
            init.inf <- unname(unlist(parameters[grep("^init\\.inf\\.[0-9]+$",
                                                      names(parameters))]))
            init.I <-
                round(init.inf * cpx.init.infected[, population])
        } else
        {
            init.I <- round(cpx.init.infected[, abs.prevalence] * infectious.period /
                            (infectious.period + latent.period))
        }
        if (any(grep("^init\\.exp\\.[0-9]+$", names(parameters))))
        {
            init.exp <- unname(unlist(parameters[grep("^init\\.exp\\.[0-9]+$",
                                                      names(parameters))]))
            init.E <-
                round(init.exp * cpx.init.infected[, population])
        } else
        {
            init.E <- round(cpx.init.infected[, abs.prevalence] * latent.period /
                            (infectious.period + latent.period))
        }
        init.R <- cpx.init.infected[, population] - init.S - init.I - init.E
        init.Z <- rep(0, nrow(cpx.init.infected))

        ## make sure all is >0
        init.S[which(init.R < 0)] <-
            init.S[which(init.R < 0)] + init.R[which(init.R < 0)]
        init.R[which(init.R < 0)] <-  0

        if ("maternal.immunity" %in% names(parameters) &&
            parameters[["maternal.immunity"]] > 0)
        {
            init.B <- init.S[1] * parameters[["maternal.immunity"]]
            init.S[1] <- init.S[1] * (1 - parameters[["maternal.immunity"]])
            init <- c(init.B, init.S, init.E, init.I, init.R, init.Z)
            names(init) <- c("B1",
                             paste("S", seq_len(nrow(ms.init.infected)), sep = ""),
                             paste("E", seq_len(nrow(ms.init.infected)), sep = ""),
                             paste("I", seq_len(nrow(ms.init.infected)), sep = ""),
                             paste("R", seq_len(nrow(ms.init.infected)), sep = ""),
                             paste("Z", seq_len(nrow(ms.init.infected)), sep = ""))
        } else
        {
            init <- c(init.S, init.E, init.I, init.R, init.Z)
            names(init) <- c(paste("S", seq_len(nrow(cpx.init.infected)), sep = ""),
                             paste("E", seq_len(nrow(cpx.init.infected)), sep = ""),
                             paste("I", seq_len(nrow(cpx.init.infected)), sep = ""),
                             paste("R", seq_len(nrow(cpx.init.infected)), sep = ""),
                             paste("Z", seq_len(nrow(cpx.init.infected)), sep = ""))
        }

    }

    return(init)
}

##' Interpolate serology
##'
##' This splits age classes into other age classes by assuming that
##' results were randomly distributed between ages
##' @param age.limits The lower limits of the required age classes
##' @return interpolated serologies
##' @author Sebastian Funk
interpolate.serology <- function(age.limits) {

    data(cpx_sero)

    cpx.sero.age.limits <- c(unique(cpx.sero[, lower.age.limit]), 40)

    interpolating <-
        age.limits[!(age.limits %in% cpx.sero.age.limits) &
                       (age.limits < max(cpx.sero.age.limits))]

    ## create empty copy
    additional <- cpx.sero[lower.age.limit == interpolating[1]]

    for (age in interpolating) {
        ## for ones that sit between age limits, create a new one
        max.lower.age <- max(cpx.sero[lower.age.limit < age, lower.age.limit])
        this.additional <-
            cpx.sero[lower.age.limit == max.lower.age,
                     list(lower.age.limit = age,
                          antibodies = as.integer(round(antibodies / 2)),
                          samples = as.integer(round(samples / 2))),
                     by = list(year)]
        cpx.sero[lower.age.limit == max.lower.age,
                 antibodies := as.integer(round(antibodies / 2))]
        cpx.sero[lower.age.limit == max.lower.age,
                 samples := as.integer(round(samples / 2))]
        additional <- rbind(additional, this.additional)

    }

    cpx.sero <- rbind(cpx.sero, additional)

    setkey(cpx.sero, year, lower.age.limit)

    ## now, find in which interval of serum lower age limits the lower
    ## age limits are
    cpx.sero.age.limits <- c(unique(cpx.sero[, lower.age.limit]))

    intervals <- findInterval(cpx.sero[, lower.age.limit], age.limits)

    ## convert them to limits
    lower.limits <- age.limits[intervals]

    cpx.sero <- cpx.sero[, lower.age.limit := lower.limits]
    cpx.sero <- cpx.sero[, list(antibodies = sum(antibodies, na.rm = T),
                                samples = sum(samples, na.rm = T)),
                         by = list(year, lower.age.limit)]

    return(cpx.sero)
}

##' Chickenpox prior for initial conditions
##'
##' @param init initial conditions
##' @return (logged) prior
##' @author Sebastian Funk
cpx.init.prior <- function(init) {

    data(cpx_sero)
    data(cpx_ew_age)

    prior <- 0

    nagegroups <- max(as.numeric(gsub("[^0-9]", "", names(init))))

    serology <- interpolate.serology(unique(cpx.ew.age[, lower.age.limit]))
    first.serology <- serology[year < min(cpx.ew.age[, year])][year == max(year)]

    nrec <- unname(init[grep("^R", names(init))])
    N <- sapply(seq_len(nagegroups), function(x) {
        sum(init[x + nagegroups * seq(0, (length(init) / nagegroups) - 1)])
    })

    prior <- prior + sum(sapply(seq_len(nrow(first.serology)), function(x) {
        dhyper(first.serology[x, antibodies],
               nrec[x], N[x] - nrec[x], first.serology[x, samples],
               log = T)
    }))

    if (nrow(first.serology) < nagegroups) {
        prior <- prior +
            sum(sapply(seq(nrow(first.serology) + 1, nagegroups), function(x) {
                dunif(nrec[x], min = round(N[x] * nrec[x - 1] / N[x - 1]),
                      max = N[x], log = T)
            }))
    }

    prior

}

##' Evaluate parameters under a prior distribution for chickenpox
##'
##' @param parameters parameters
##' @return prior probability density
##' @author Sebastian Funk
##' @export
cpx.prior <- function(parameters)
{

    prior <- 0

    prior <- prior + dunif(parameters[["R0"]], 1, 20, log = T)

    if ("brownian.init" %in% names(parameters))
    {
        prior <- prior + dunif(parameters[["brownian.init"]], 0, 2, log = T)
    }
    if ("brownian.drift" %in% names(parameters))
    {
        prior <- prior + dunif(parameters[["brownian.drift"]], -0.1, 1, log = T)
    }
    if ("brownian.scale" %in% names(parameters))
    {
        prior <- prior + dunif(parameters[["brownian.scale"]], 0, 0.5, log = T)
    }
    if ("early.mixing" %in% names(parameters))
    {
        prior <- prior + dunif(parameters[["early.mixing"]], 0.5, 5, log = T)
    }
    if ("late.mixing" %in% names(parameters))
    {
        prior <- prior + dunif(parameters[["late.mixing"]], 0.5, 5, log = T)
    }
    if ("first.change" %in% names(parameters))
    {
        prior <- prior + dunif(parameters[["first.change"]], 1970, 1990,
                               log = T)
        if (is.finite(prior) &&
            "second.change" %in% names(parameters) &&
            parameters[["second.change"]] >
            parameters[["first.change"]] + 10)
        {
            prior <- prior + dunif(parameters[["second.change"]],
                                   parameters[["first.change"]] + 10, 2010,
                                   log = T)
        } else {
            prior <- -Inf
        }
    }

    prior <-
        prior + dunif(parameters[["gamma"]], 1/(11/365.25), 1/(10/365.25), log = T)
    prior <-
        prior + dunif(parameters[["rho"]], 1/(12/365.25), 1/(8/365.25), log = T)
    prior <-
        prior + dunif(parameters[["underreporting"]], 0.1, 0.9, log = T)

    if ("maternal.immunity" %in% names(parameters))
    {
        prior <-
            prior + dunif(parameters[["maternal.immunity"]], 0, 1, log = T)
    }

    for (init.param in grep("^init\\.", names(parameters), value = TRUE))
    {
        prior <-
            prior + dunif(parameters[[init.param]], 0, 1, log = T)
    }

    prior
}

##' Randomly sample initial conditions for measles or chickenpox on the basis of serology
##'
##' @param parameters given (maternal immunity and mixing) parameters
##' @param infection "measles" or "chickenpox"
##' @return proportion susceptible in each age group
##' @author Sebastian Funk
sample.init <- function(parameters, infection = c("measles", "chickenpox")) {

    infection <- match.arg(infection)

    prop.susc <- c()

    if (infection == "measles")
    {

        data(ms_sero_fine_1950)
        data(ms_ew_age)

        lower.limits <- agegroups.to.limits(colnames(parameters[["mixing"]]))

        min.year <- ms.ew.age[, min(year)]
        first.pop.state <- ms.ew.age[year == min.year,
                                     list(population = mean(population)),
                                     by = lower.age.limit]
        first.pop.state[, lower.age.limit :=
                              reduce.agegroups(lower.age.limit, lower.limits)]
        first.pop.state <- first.pop.state[, list(population = sum(population)),
                                           by = lower.age.limit]

        if ("maternal.immunity" %in% names(parameters))
        {
            maternal.immunity <- parameters[["maternal.immunity"]]
        } else
        {
            maternal.immunity <- 0
        }
        first.serology <-
            rbind(data.table(age = 0, immunity = maternal.immunity),
                  ms.sero.fine)
        setnames(first.serology, "age", "lower.age.limit")
        first.serology[, lower.age.limit :=
                             reduce.agegroups(lower.age.limit, lower.limits)]
        first.serology <-
            first.serology[, list(immunity = mean(immunity)), by = lower.age.limit]

        setkey(first.serology, lower.age.limit)
        setkey(first.pop.state, lower.age.limit)

        first.serology <- merge(first.serology, first.pop.state, all.y = T)

        if (max(first.serology[, lower.age.limit]) > 0)
        {
            max.immunity <- max(first.serology[lower.age.limit > 0, immunity], na.rm  = T)
        }
        first.serology <- first.serology[is.na(immunity), immunity := max.immunity]

        prop.susc <- first.serology[, 1 - immunity]
    } else if (infection == "chickenpox") {
        data(cpx_sero)
        data(cpx_ew_age)

        lower.limits <- agegroups.to.limits(colnames(parameters[["mixing"]]))

        serology <- interpolate.serology(lower.limits)
        before.serology <- serology[year < min(cpx.ew.age[, year])]
        max.before <- max(before.serology[, year])
        first.serology <- before.serology[year == max.before]
        min.year <- min(cpx.ew.age[, year])
        first.pop.state <- cpx.ew.age[year == min.year,
                                      list(population = mean(population)),
                                      by = lower.age.limit]
        first.pop.state[, lower.age.limit :=
                              reduce.agegroups(lower.age.limit, lower.limits)]
        first.pop.state <- first.pop.state[, list(population = sum(population)),
                                           by = lower.age.limit]

        setkey(first.serology, lower.age.limit)
        setkey(first.pop.state, lower.age.limit)

        first.serology <- merge(first.serology, first.pop.state)

        first.sampler <-
            first.serology[, list(mean = antibodies * population / samples,
                                  sd = sqrt((population -
                                             antibodies * population / samples) *
                                            antibodies * population / (samples^2)*
                                            (population - samples) /(population - 1)),
                                  population = population,
                                  prop.immune = -1),
                           by = lower.age.limit]

        while (any(first.sampler[, prop.immune > 1 | prop.immune < 0])) {
            first.sampler <-
                first.sampler[, list(prop.immune = ifelse(prop.immune > 1 |
                                                          prop.immune < 0,
                                                          rnorm(1, mean, sd) / population, prop.immune),
                                     mean = mean, sd = sd, population = population),
                              by = lower.age.limit]
        }

        first.sampler <- first.sampler[, list(lower.age.limit, prop.immune)]

        limits.diff <- setdiff(lower.limits,
                               unique(first.sampler[, lower.age.limit]))
        if (length(limits.diff) > 0) {
            missing.sampler <- data.table(lower.age.limit = limits.diff)
            missing.sampler <-
                missing.sampler[1, prop.immune := runif(1, first.sampler[nrow(first.sampler), prop.immune], 1)]

            if (nrow(missing.sampler > 1)) {
                for (i in 2:nrow(missing.sampler)) {
                    missing.sampler <-
                        missing.sampler[i, prop.immune := runif(1, missing.sampler[i - 1, prop.immune], 1), by = lower.age.limit]
                }
            }
            first.sampler <- rbind(first.sampler, missing.sampler)
        }

        prop.susc <- first.sampler[, 1 - prop.immune]

    }
    names(prop.susc) <- paste("init.susc", seq_along(prop.susc), sep = ".")

    return(prop.susc)
}

##' Evaluate parameters under a prior distribution for measles
##'
##' @param parameters parameters
##' @return prior probability density
##' @author Sebastian Funk
##' @export
ms.prior <- function(parameters) {

    prior <- 0

    prior <- prior + dunif(parameters[["R0"]], 1, 30, log = T)

    prior <-
        prior + dunif(parameters[["gamma"]], 1/(7/365.25), 6/(10/365.25), log = T)
    prior <-
        prior + dunif(parameters[["rho"]], 1/(7/365.25), 1/(6/365.25), log = T)
    prior <-
        prior + dunif(parameters[["underreporting"]], 0.1, 0.9, log = T)
    if ("child.mixing" %in% names(parameters))
    {
        prior <-
            prior + dunif(parameters[["child.mixing"]], 0.5, 5, log = T)
    }
    if ("maternal.immunity" %in% names(parameters))
    {
        prior <-
            prior + dunif(parameters[["maternal.immunity"]], 0, 1, log = T)
    }

    return(prior)
}

##' Sample parameters from a prior distribution for measles or chickenpox
##'
##' @param mixing given mixing matrix, if this is not to sample
##' @param infection "measles" or "chickenpox"
##' @return prior parameters
##' @author Sebastian Funk
##' @export
sample.prior <- function(fixed.parameters = list(), infection = c("measles", "chickenpox"))
{

    infection <- match.arg(infection)

    parameters <- fixed.parameters

    if (infection == "measles")
    {
        data(ms_ew_age)

        if (!("R0" %in% names(fixed.parameters)))
        {
            parameters[["R0"]] <- runif(1, 1, 30)
        }

        if (!("gamma" %in% names(fixed.parameters)))
        {
            parameters[["gamma"]] <- runif(1, 1/(7/365.25), 1/(6/365.25))
        }

        if (!("rho" %in% names(fixed.parameters)))
        {
            parameters[["rho"]] <- runif(1, 1/(7/365.25), 1/(6/365.25))
        }

        if (!("underreporting" %in% names(fixed.parameters)))
        {
            parameters[["underreporting"]] <- runif(1, 0.1, 0.9)
        }

        if ("child.mixing.group" %in% names(fixed.parameters) &
            !("child.mixing" %in% names(fixed.parameters)))
        {
            parameters[["child.mixing"]] <- runif(1, 0.5, 2)
        }

        if (!("maternal.immunity" %in% names(fixed.parameters)))
        {
            parameters[["maternal.immunity"]] <- runif(1, 0, 1)
        }

    } else if (infection == "chickenpox")
    {
        data(cpx_ew_age)

        if (!is.null(fixed.parameters[["mixing.dynamics"]]))
        {
            if (fixed.parameters[["mixing.dynamics"]] == "brownian")
            {
                if (!("brownian.init" %in% names(fixed.parameters)))
                {
                    parameters[["brownian.init"]] <- runif(1, 0, 2)
                }
                if (!("brownian.drift" %in% names(fixed.parameters)))
                {
                    parameters[["brownian.drift"]] <- runif(1, -0.1, 0.1)
                }
                if (!("brownian.scale" %in% names(fixed.parameters)))
                {
                    parameters[["brownian.scale"]] <- runif(1, 0, 0.5)
                }
            } else if (fixed.parameters[["mixing.dynamics"]] == "linear")
            {
                if (!is.null(fixed.parameters[["mixing.knots"]]))
                {
                    mixing.knots <- fixed.parameters[["mixing.knots"]]
                } else
                {
                    mixing.knots <- 1
                }

                for (i in seq_len(mixing.knots))
                {
                    param.name <- paste("mixing", i, sep = ".")

                    if (!(param.name %in% names(fixed.parameters)))
                    {
                        parameters[[param.name]] <- runif(1, 0.5, 5)
                    }
                }

                if (mixing.knots > 2)
                {
                    years <- unique(cpx.ew.age[, year])
                    intermediate.years <-
                        setdiff(years, c(min(years), max(years)))
                    change.points <-
                        sort(sample(intermediate.years, mixing.knots - 2))

                    for (i in seq_along(change.points))
                    {
                        param.name <- paste("change", i, sep = ".")

                        if (!(param.name %in% names(fixed.parameters)))
                        {
                            parameters[[param.name]] <- change.points[i]
                        }
                    }
                }
            }
        }

        if (!("R0" %in% names(fixed.parameters)))
        {
            parameters[["R0"]] <- runif(1, 1, 20)
        }

        if (!("gamma" %in% names(fixed.parameters)))
        {
            parameters[["gamma"]] <- runif(1, 1/(11/365.25), 1/(10/365.25))
        }
        if (!("rho" %in% names(fixed.parameters)))
        {
            parameters[["rho"]] <- runif(1, 1/(12/365.25), 1/(8/365.25))
        }
        if (!("underreporting" %in% names(fixed.parameters)))
        {
            parameters[["underreporting"]] <- runif(1, 0.1, 0.9)
        }
        if (!("maternal.immunity" %in% names(fixed.parameters)))
        {
            parameters[["maternal.immunity"]] <- runif(1, 0, 1)
        }
    }

    parameters <- c(parameters,
                    sample.init(fixed.parameters, infection))

    return(parameters)
}

##' Log-likelihood for a model run of an age-structured
##' childhood infection model
##'
##' @param trajectory The model trajectory, in a data.table with a
##' column "year", a column "lower.age.limit", and at least one of
##' "rel.incidence" or "abs.incidence"
##' @param data.cases The data trajectory, in the same format as the model trajectory
##' @param data.sero Serology, if any
##' @param rel Whether to fit relative incidences
##' @param sample Whether the data are the result of a sampling procedure
##' @param proportion.reported The average proportion of cases reported
##' @param normal.approximation Whether to use the normal
##' approximation to the binomial distribution
##' @return A number, the log-likelihood
##' @export
##' @author Sebastian Funk
age.likelihood <- function(trajectory, data.cases, data.sero = NULL, rel = F,
                           sample = NULL, proportion.reported = 1,
                           normal.approximation = F, time.column = "time") {

    if (rel)
    {
        trajectory[, abs.incidence := rel.incidence * population]
        data.cases[, abs.incidence := rel.incidence * population]
    }

    trajectory[, lower.age.limit :=
                     reduce.agegroups(lower.age.limit,
                                      unique(data.cases[, lower.age.limit]))]
    trajectory[, list(abs.incidence = sum(abs.incidence)),
               by = list(get(time.column), lower.age.limit)]

    setkeyv(trajectory, c(time.column, "lower.age.limit"))
    setkeyv(data.cases, c(time.column, "lower.age.limit"))

    ll <- NA

    model.data <- merge(trajectory, data.cases, suffixes = c(".model", ".data"))
    if (is.null(sample))
    {
        if (normal.approximation)
        {
            model.data <-
                model.data[, likelihood :=
                               dnorm(x = abs.incidence.data,
                                     mean = abs.incidence.model *
                                         proportion.reported,
                                     sd = abs.incidence.model *
                                         proportion.reported *
                                             (1 - proportion.reported),
                                     log = T)]
        }
        else
        {
            model.data <-
                model.data[, likelihood :=
                               dbinom(abs.incidence.data,
                                      abs.incidence.model,
                                      proportion.reported, log = T)]
        }
    }
    else
    {
        model.data <-
            model.data[, likelihood :=
                           dhyper(round(abs.incidence * sample),
                                  as.integer(abs.incidence.model),
                                  as.integer(population - abs.incidence.model),
                                  round(sample), log = T)]
    }

    if (nrow(model.data) > 0)
    {
        ll <- sum(model.data[, likelihood], na.rm = T)
    }

    if (!is.null(data.sero))
    {
        setkeyv(data.sero, c(time.column, "lower.age.limit"))
        model.sero <- merge(trajectory, data.sero, suffixes = c(".model", ".sero"))
        model.sero <-
            model.sero[, likelihood := dhyper(antibodies, as.integer(R),
                                    as.integer(N - R), samples, log = T),
                       by = list(time.column, "lower.age.limit")]
        if (nrow(model.sero) > 0)
        {
            ll <- sum(model.data[, likelihood], na.rm = T)
        }
    }

    return(ll)
}

##' Generate sample observations for the age-structured model
##'
##' @param parameters model parameters
##' @param incidence model incidence
##' @param normal.approximation whether to use the normal approximation
##' @return a sample trajcetory
##' @author Sebastian Funk
##' @export
sample.age.obs <- function(parameters, incidence, normal.approximation = F)
{
    if (normal.approximation)
    {
        obs <- rnorm(length(incidence),
                     mean = incidence * parameters[["underreporting"]],
                     sd = incidence * parameters[["underreporting"]] *
                         (1 - parameters[["underreporting"]]))
    }
    else
    {
        obs <- rbinom(length(incidence),
                      size = incidence,
                      prob = parameters[["underreporting"]])
    }

    return(obs)
}

##' Plot a model and data
##'
##' @param parameters model parameters
##' @param data data to plot
##' @param init initial conditions
##' @param time.column time column in the data
##' @param sample.obs if given, this well be used to sample observations
##' @return plot
##' @author Sebastian Funk
##' @export plot.model.data
##' @import ggplot2
plot.model.data <- function(parameters, data, init, time.column = "time",
                            sample.obs = NULL, traj = FALSE) {

    mtraj <- runSIR(parameters = parameters,
                    init = init,
                    times = unique(data[, get(time.column)]),
                    age.labels = unique(data[, lower.age.limit]))
    mtraj <- data.table(mtraj[, list(time,
                                      lower.age.limit = factor(lower.age.limit),
                                      abs.incidence)],
                         source = "model")

    model.data <-
        rbind(data.table(mtraj[, list(time,
                                      lower.age.limit = factor(lower.age.limit),
                                      abs.incidence)],
                         source = "model"),
              data.table(data[, list(time = get(time.column),
                                     lower.age.limit = factor(lower.age.limit),
                                     abs.incidence)],
                         source = "data"))
    model.data <-
        model.data[!(time %in% unique(model.data[is.na(abs.incidence), time]))]

    if (!is.null(sample.obs)) {
        model.data[source == "model",
                   abs.incidence := sample.obs(parameters, abs.incidence)]
    }

    setnames(model.data, "time", time.column)

    age.labels <- limits.to.agegroups(unique(data[, lower.age.limit]))
    p <- ggplot(model.data, aes_string(x = time.column,
                                       y = "abs.incidence",
                                       linetype = "source",
                                       color = "lower.age.limit",
                                       shape = "source"))
    if (length(unique(model.data[, get(time.column)])) == 1) {
        p <- p + geom_point(size = 4)
    } else {
        p <- p + geom_line(lwd = 1.5)
    }
    p <- p + scale_color_brewer("Age group", palette = "Set1",
                                labels = age.labels)
    p <- p + theme_bw(16)
    p <- p + theme(legend.position = "bottom")

    if (traj)
    {
        return(list(plot = p, traj = mtraj))
    } else
    {
        return(p)
    }

}

##' Plot a sampled run and data
##'
##' @param theta parameter vector
##' @param data The data to plot
##' @param init.sample function to sample initial conditions
##' @param time.column time column the data
##' @param fixed.parameters fixed parameters
##' @param sample.obs sample observation function
##' @param ... further parameters for \code{plot.model.data}
##' @return plot
##' @author Sebastian Funk
##' @export plot.sampled.run
plot.sampled.run <- function(theta, data, init.sample, time.column = "time",
                             fixed.parameters = NULL, sample.obs = NULL, ...) {

    ew.pop.parameters <-
        pop.parameters(year.limits = unique(data[, get(time.column)]),
                       age.limits = unique(data[, lower.age.limit]),
                       interpolate = !is.null(fixed.parameters[["interpolate"]]))

    parameters <-
        c(theta,
          fixed.parameters,
          ew.pop.parameters)

    plot.model.data(parameters, data, init.sample(parameters), time.column,
                    sample.obs, ...)
}

##' Fit the chickenpox data
##'
##' Fit a model with a free contact rate between <5 year olds to the
##' England & Wales chickenpox data, assuming everything is in
##' equilibrium every year (a method similar to Wallinga, Teunis and
##' Kretzschmar, 2006, "Using Data on Social Contacts to Estimate
##' Age-specific Transmission Parameters for Respiratory-spread
##' Infectious Agents", Am J Epidemiol, 164, 936-944. This estimates
##' q, the reporting rate, and the mixing between children every year.
##' @param years The years to fit (of not given, fit all years)
##' @param contact.matrix The contact matrix to use (if not given,
##' sample one from UK Polymod)
##' @return A list with the best fits (q, rep and child.mixing), the
##' used contact.matrix, and the return of fitting q and rep with any
##' messages and convergence information that may have been reported.
fit.chickenpox.wallinga <- function(years = NULL, contact.matrix = NULL) {

    data(cpx_ew_age)

    ## if contact matrix is not given, sample one from polymod
    if (is.null(contact.matrix)) {
        contact.matrix <- sample.uk.polymod(unique(cpx.ew.age[, lower.age.limit]))
    }

    ## if years is not given, fit all years
    if (is.null(years)) {
        years <- unique(cpx.ew.age[, year])
    }

    ## convert chickenpox data into annual data
    cpx.annual <- cpx.ew.age[, list(abs.incidence = sum(abs.incidence),
                                    rel.incidence = sum(rel.incidence),
                                    population = sum(population)),
                             by = list(year, lower.age.limit)]

    ## do the fit
    fit <- optim(par = c(-0.5, -0.5), function(x) {
        sum(sapply(years, function(y) {
            optimize(lower = -1, upper = 1, function(z) {
                age.dist <- endemic.age.dist(contact.matrix, 10^x[1], 10^z)
                year.fit <-
                    data.table(year = y, rel.incidence = age.dist[["y"]],
                               lower.age.limit = unique(cpx.ew.age[, lower.age.limit]))
                a <- -age.likelihood(year.fit, cpx.annual, reporting.rate = 10^x[2],
                                     rel = T)
                a
            })$objective
        }))
    })

    q <- 10^fit$par[1]
    rep <- 10^fit$par[2]

    child.mixing <- sapply(years, function(y) {
        10^optimize(lower = -1, upper = 1, function(z) {
            age.dist <-
                endemic.age.dist(contact.matrix, q, 10^z)
            year.fit <-
                data.table(year = y, rel.incidence = age.dist[["y"]],
                           lower.age.limit = unique(cpx.ew.age[, lower.age.limit]))
            -age.likelihood(year.fit, cpx.annual, reporting.rate = rep, rel = T)
        })$minimum
    })

    return(list(q = q, rep = rep, child.mixing = child.mixing,
                contact.matrix = contact.matrix, fit = fit))
}

##' Fit the measles data
##'
##' Fit a model with a free contact rate between <5 year olds to the
##' England & Wales measles data, assuming everything is in
##' equilibrium every year (a method similar to Wallinga, Teunis and
##' Kretzschmar, 2006, "Using Data on Social Contacts to Estimate
##' Age-specific Transmission Parameters for Respiratory-spread
##' Infectious Agents", Am J Epidemiol, 164, 936-944. This estimates
##' q, the reporting rate, and the mixing between children every year.
##' @param years The years to fit (of not given, fit all years)
##' @param contact.matrix The contact matrix to use (if not given,
##' sample one from UK Polymod)
##' @return A list with the best fits (q, rep and child.mixing), the
##' used contact.matrix, and the return of fitting q and rep with any
##' messages and convergence information that may have been reported.
fit.measles.wallinga <- function(years = NULL, contact.matrix = NULL) {

    data(ms_ew_age)

    age.groups <-
        as.integer(c(0, unique(ms.ew.age[lower.age.limit >= 5, lower.age.limit])))

    if (is.null(contact.matrix)) {
        contact.matrix <- sample.uk.polymod(age.groups)
    }
    if (is.null(years)) {
        years <- unique(ms.ew.age[, year])
    }

    ms.ew.age[, lower.age.limit :=
                    reduce.agegroups(lower.age.limit, age.groups)]
    ms.ew.age <- ms.ew.age[, list(abs.incidence = sum(abs.incidence),
                                  population = sum(population)),
                           by = list(year, lower.age.limit)]
    ms.ew.age <- ms.ew.age[, rel.incidence := abs.incidence / population]

    mmr <- estimate.immunity.mmr(years = years)

    fit <- optim(par = c(log10(0.5), log10(0.5)), function(x) {
        s <- sum(sapply(years, function(y) {
            optim(par = log10(contact.matrix[1,1]), fn = function(z) {
                cat(y, z, "\n")
                age.dist <-
                    endemic.age.dist(contact.matrix, 10^x[1], 10^z,
                                     vaccination = unname(unlist(mmr[year == y, -1,
                                     with = F])))
                year.fit <-
                    data.table(year = y, rel.incidence = age.dist[["y"]],
                               lower.age.limit = unique(ms.ew.age[, lower.age.limit]))
                a <- -age.likelihood(year.fit, ms.ew.age, sample = 10^x[2],
                                     rel = T)
                cat("a =", a, "\n")
                a
            })$value
        }))
        cat(x, s, "\n")
        s
    })

    q <- 10^fit$par[1]
    rep <- 10^fit$par[2]

    child.mixing <- sapply(years, function(y) {
        cat(y, "\n")
        10^optim(par = 0, function(z) {
            age.dist <-
                endemic.age.dist(contact.matrix, q, 10^z,
                                 vaccination = unname(unlist(mmr[year == y, -1, with = F])))
            year.fit <-
                data.table(year = y, rel.incidence = age.dist[["y"]],
                           lower.age.limit = unique(ms.ew.age[, lower.age.limit]))
            -age.likelihood(year.fit, ms.ew.age, sample = rep, rel = T)
        })$minimum
    })

    return(list(q = q, rep = rep, child.mixing = child.mixing))
}

##' Calculate age distribution in equilibrium
##'
##' Calculates the age distribution in equilibrium of an age-strutured
##' SIR-model, by year or as a function of R0
##' @param agegroups Vector of the lower limits of the age groups
##' @param R0 The value(s) of R0. If a vactor of length>1, will test
##' the different values
##' @param vaccination Vaccination levels, given as a data table with
##' a "year" column, and a column for each of the age groups
##' @param years Years for which to calculate age distribution. This
##' takes into account vaccination in the different years. If this is
##' given, R0 should be set to a single value
##' @return a data.table with the different relative incidences,
##' population sizes, absolute incidences and proportions of cases
##' @author Sebastian Funk
age.distribution <- function(contact.matrix, R0, vaccination = NULL, years = NULL) {

    data(pop_ew_age)

    ## check consistency of parameters
    if ((length(R0) > 1) && !is.null(years)) {
        warning("distribution by year requires a single R0 value")
        R0 <- R0[1]
    }
    if (!is.null(vaccination)) {
        if (is.null(years)) {
            if (length(vaccination == 1)) {
                vaccination <- rep(vaccination, ncol(contact.matrix))
            } else if (!(length(vaccination) == ncol(contact.matrix))) {
                stop("vaccination must have the as many elements as the contact matrix columns")
            }
        } else {
            if (!("year" %in% colnames(vaccination))) {
                stop("vaccination must have a 'year' column")
            }
            if (!(ncol(vaccination) == ncol(contact.matrix) + 1)) {
                stop("vaccination must have one more column than the contact matrix")
            }
        }
    }

    ## get the lower age limits from the headings of the contact matrix
    agegroups <-
        as.integer(gsub("^\\[([0-9]*),[0-9]*\\)", "\\1", colnames(contact.matrix)))

    ## get populations by age group -- if no years are given, take 2006 (like polymod)
    age.years <- years
    if (is.null(age.years)) {
        age.years <- 2006
    }
    ages <- pop.ew.age[year %in% age.years]
    ages[, lower.age.limit := reduce.agegroups(lower.age.limit, agegroups)]
    if (is.null(years)) {
        ages <- ages[, list(population = sum(population)), by = lower.age.limit]
    } else {
        ages <- ages[, list(population = sum(population)),
                     by = list(year, lower.age.limit)]
        ## fill older population structure with the oldest known
        if (min(ages[, year]) - min(years) > 0) {
            for (test.year in years[years < min(ages[, year])]) {
                new.year <- ages[year == min(year)]
                new.year <- new.year[, year := test.year]
                ages <- rbind(new.year, ages)
            }
        }
    }

    ## get q
    ## calculate age distributions
    if (is.null(years)) {
        ## to go from the contact matrix to a next-generation matrix, we
        ## need to multiply by N_i/N_j
        ni.matrix <- matrix(rep(ages[, population], each = ncol(contact.matrix)),
                            ncol = ncol(contact.matrix), byrow = T)
        nj.matrix <- matrix(rep(ages[, population], each = nrow(contact.matrix)),
                            nrow = nrow(contact.matrix))
        ngm.multiplier <- ni.matrix / nj.matrix
        ## get the matrix-R0 (R0/q)
        R0.matrix <- max(Re(eigen(contact.matrix * ngm.multiplier)$value))

        ## get q
        q <- R0 / R0.matrix

        if (!is.null(vaccination)) {
            current.vaccination <- vaccination
        } else {
            current.vaccination <- rep(0, length(agegroups))
        }

        age.dists <- data.table(t(sapply(q, function(x) {
            c(R0 = x * R0.matrix, endemic.age.dist(contact.matrix, x,
              vaccination = current.vaccination)$y,
              vaccination = current.vaccination)
        })))
    } else {
        age.dists <- data.table(t(sapply(years, function(y) {
            ## to go from the contact matrix to a next-generation matrix, we
            ## need to multiply by N_i/N_j
            ni.matrix <- matrix(rep(ages[year == y, population],
                                    each = ncol(contact.matrix)),
                                ncol = ncol(contact.matrix), byrow = T)
            nj.matrix <- matrix(rep(ages[year == y, population],
                                    each = nrow(contact.matrix)),
                                nrow = nrow(contact.matrix))
            ngm.multiplier <- ni.matrix / nj.matrix
            ## get the matrix-R0 (R0/q)
            R0.matrix <- max(Re(eigen(contact.matrix * ngm.multiplier)$value))

            ## get q
            q <- R0 / R0.matrix

            if (!is.null(vaccination)) {
                current.vaccination <-
                    unname(unlist(vaccination[year == y, !"year", with = F]))
            } else {
                current.vaccination <- rep(0, length(agegroups))
            }
            c(year = y, endemic.age.dist(contact.matrix, q,
              vaccination = current.vaccination)$y,
              vaccination = current.vaccination)
        })))
    }

    ## year or R0
    by <- names(age.dists)[1]

    ## construct results data table
    setnames(age.dists, seq(2, 1 + ncol(contact.matrix)), as.character(agegroups))
    if (is.null(vaccination)) {
        age.dists <- age.dists[, seq_len(1 + ncol(contact.matrix)), with = F]
    } else {
        setnames(age.dists, seq(2 + ncol(contact.matrix), 1 + 2 * ncol(contact.matrix)),
                 paste("vacc", as.character(agegroups), sep = "."))
    }
    m.age.dists <-
        melt(age.dists, id.vars = by, variable.name = "lower.age.limit")
    m.age.dists <- m.age.dists[, variable := "rel.incidence"]
    m.age.dists <- m.age.dists[grep("^vacc\\.[0-9]+$", lower.age.limit),
                               variable := "vaccination"]
    m.age.dists[, lower.age.limit := sub("^vacc\\.", "", lower.age.limit)]

    m.age.dists <-
        m.age.dists[, lower.age.limit := as.integer(as.character(lower.age.limit))]
    if (is.null(years))
    {
        m.age.dists <- data.table(dcast(m.age.dists, R0 + lower.age.limit ~ variable))
    } else {
        m.age.dists <- data.table(dcast(m.age.dists, year + lower.age.limit ~ variable))
    }

    if (is.null(years)) {
        setkey(m.age.dists, lower.age.limit)
        setkey(ages, lower.age.limit)
    } else {
        setkey(m.age.dists, year, lower.age.limit)
        setkey(ages, year, lower.age.limit)
    }

    m.age.dists <- merge(m.age.dists, ages)
    m.age.dists <- m.age.dists[, abs.incidence := round(rel.incidence * population)]

    m.age.dists.total <-
        m.age.dists[, list(abs.incidence = sum(abs.incidence)), by = list(get(by))]
    setnames(m.age.dists.total, "get", by)
    setkeyv(m.age.dists.total, by)
    setkeyv(m.age.dists, by)
    m.age.dists <- merge(m.age.dists, m.age.dists.total, suffixes = c("", ".total"))
    m.age.dists <-
        m.age.dists[, proportion.cases := abs.incidence / abs.incidence.total]

    m.age.dists[, abs.incidence.total := NULL]

    return(m.age.dists)
}

##' Take random samples from the posterior
##'
##' @return prior, likelihood, posterior, parameters, initial values
##' @author Sebastian Funk
##' @param sample.prior a function to sample from prior
##' @param eval.prior a function to evaluate the prior
##' @param gen.init a function to generate the initial conditions
##' @param data the data to compare the model to
##' @param test.init initial conditions (proportoin susceptible in
##' each each group)
##' @param fixed.parameters (list of) fixed parameters of the model
##' @param likelihood.arguments (list of) parameters to be passed to
##' age.likelihood
##' @param mcmc.arguments (list of) parameters to be passed to
##' mcmcMH
##' @param time.column the name of the time column in the data
##' @param tmax maximum time to run the simulations for
##' @param thickening thickening innterval
##' @param interpolate whether to intepolate births/deaths etc
##' @param nSamples number of samples (if this is 1, the sample will
##' be from the prior)
##' @param update.init whether to update the initial conditions
##' @param return.traj whether to return the whol trajectory
##' @export
sample.posterior <- function(sample.prior, eval.prior,
                             gen.init, data,
                             test.init = NULL,
                             fixed.parameters = list(),
                             likelihood.arguments = list(),
                             mcmc.arguments = list(),
                             time.column = "time", tmax = NULL,
                             thickening = 1, interpolate = T,
                             nSamples = 1, update.init = F)
{

    prior.parameters <-
        sample.prior(fixed.parameters = fixed.parameters)

    ew.pop.parameters <-
        pop.parameters(year.limits = unique(data[, get(time.column)]),
                       age.limits = unique(data[, lower.age.limit]),
                       interpolate = !is.null(fixed.parameters[["interpolate"]]))

    times <- unique(data[, get(time.column)])
    if (!is.null(tmax))
    {
        times <- c(times, tmax)
    }

    posterior.func <- function(theta)
    {
        parameters <-
            c(theta,
              fixed.parameters,
              ew.pop.parameters)

        init <- gen.init(parameters)

        lp <- eval.prior(parameters)
        if (is.finite(lp))
        {
            mtraj <- runSIR(parameters = parameters,
                            init = init,
                            times = times,
                            age.labels = unique(data[, lower.age.limit]),
                            thickening = thickening)
            mtraj <- mtraj[time %in% times]
            setnames(mtraj, "time", time.column)
            ll.args <- c(list(trajectory = mtraj, data.cases = data,
                              proportion.reported = parameters[["underreporting"]],
                              time.column = time.column), likelihood.arguments)
            ll <- do.call(what = age.likelihood, args = ll.args)
            log.posterior  <- lp + ll
        } else
        {
            ll <- NA
            log.posterior <- lp
            mtraj <- NULL
        }

        if (update.init)
        {
            init.susc <- mtraj[get(time.column) == times[1], S / (S + E + I + R)]
            theta[paste("init", "susc", seq_along(unique(data[, lower.age.limit])),
                        sep = ".")] <- init.susc
            init.exp <- mtraj[get(time.column) == times[1], E / (S + E + I + R)]
            theta[paste("init", "exp", seq_along(unique(data[, lower.age.limit])),
                        sep = ".")] <- init.exp
            init.inf <- mtraj[get(time.column) == times[1], I / (S + E + I + R)]
            theta[paste("init", "inf", seq_along(unique(data[, lower.age.limit])),
                        sep = ".")] <- init.inf
        }

        return(list(log.density = ll + lp,
                    trace = c(theta,
                    log.likelihood = ll, log.prior = lp,
                    log.posterior = log.posterior)))
    }


    res <-
        do.call(seir.mcmc, c(list(target = posterior.func,
                                 init.theta = prior.parameters,
                                 n.iterations = nSamples - 1),
                            mcmc.arguments))

    return(res)
}

##' Update initial conditions to the start of the actual run (not run-in)
##'
##' @param parameters model parameters
##' @param data data (needed for run times and age limits etc)
##' @param gen.init function to generate initial conditions
##' @param thickening whether to thicken the runs
##' @param time.column time column in the data
##' @return updated initial conditions
##' @author Sebastian Funk
##' @export
forward.initial.conditions <- function(parameters, data, gen.init, thickening = 0,
                                       time.column = "time")
{

    ew.pop.parameters <-
        pop.parameters(year.limits = unique(data[, get(time.column)]),
                       age.limits = unique(data[, lower.age.limit]),
                       interpolate = !is.null(parameters[["interpolate"]]))

    parameters <- c(parameters, ew.pop.parameters)

    init <- gen.init(parameters)

    times <- unique(data[, get(time.column)])[c(1, 2)]

    mtraj <- runSIR(parameters = parameters,
                    init = init,
                    times = times,
                    age.labels = unique(data[, lower.age.limit]),
                    thickening = thickening)

    prop.susc <-
        mtraj[, S] / rowSums(mtraj[, seq(4, ncol(mtraj) - 1), with = F], na.rm = T)
    names(prop.susc) <- paste("init.susc", seq_along(prop.susc), sep = ".")

    return(prop.susc)
}

