##' Samples polymod using weights following Baguelin et al. (2013), but
##' using a bootstrap; first, contacts and ages are sampled
##'
##' @param n number of matrices to sample
##' @param age.limits Lower limits of the age groups; if not given, will use 5 year age limits as in the population data
##' @param survey either a survey ("POLYMOD") or a list of 'participants' and 'contacts' (both data frames) to sample from
##' @param countries limit to one or more countries; if not given, will use all countries in the survey
##' @param survey.pop survey population -- either a data frame with columns lower.age.limit and population, or a character vector giving the name(s) to use with the 2013 WHO population; if not given, will use the country populations from the desired countries, or all countries in the survey if \code{countries} is not given
##' @param bootstrap whether to sample using a bootstrap; will be set to TRUE if n > 1
##' @param symmetric whether to make matrix symmetric
##' @param normalise whether to normalise to eigenvalue 1
##' @param split whether to split the number of contacts and assortativity
##' @param part.age.column column indicating age in participant data
##' @param contact.age.column column indicating age in contact data
##' @param id.column column to match participants with contacts
##' @param dayofweek.column column indicating the day of the week
##' @param country.column column indicating the country
##' @param add.weights additional weight columns (e.g., minutes etc
##' @param symmetry make contact matrix symmetric
##' @return a list of sampled contact matrices
##' @import wpp2015 
##' @importFrom reshape2 melt
##' @export
##' @author Sebastian Funk
contact_matrix <- function(n = 1, age.limits, survey = "polymod", countries, survey.pop, mixing.pop, bootstrap = FALSE,  symmetric = TRUE, normalise = FALSE, split = FALSE, add.weights = c(), part.age.column = "participant_age", contact.age.column = "cnt_age_mean", id.column = "global_id", dayofweek.column = "day_of_week", country.column = "country", year.column = "year")
{
    ## load population data if necessary
    if ((missing(survey.pop) || is.character(survey.pop)) &&
        (missing(mixing.pop) || is.character(mixing.pop)))
    {
        data(popM)
        data(popF)

        popM <- data.table(popM)
        popF <- data.table(popF)

        popM <- popM[, sex := "male"]
        popF <- popF[, sex := "female"]

        pop <- rbind(popM, popF)

        pop <- melt(pop, id.vars = c("country", "country_code", "age", "sex"), variable.name = "year")
        pop <- data.table(dcast(pop, country + country_code + age + year ~ sex, value.var = "value"))
        pop <- pop[, year := as.integer(as.character(year))]
        pop <- pop[, lower.age.limit := as.integer(sub("[-+].*$", "", age))]
        pop <- pop[, list(country, lower.age.limit, year, population = female + male)]
    }

    ## check if survey is given as character
    if (is.character(survey))
    {
        if (tolower(survey) == "polymod")
        {
            survey_data <- list(participants = polymod$participants,
                                contacts = polymod$contacts)
        } else
        {
            stop("Unknown survey, '", survey, "''")
        }
    } else
    {
        if (!is.list(survey) || is.null(names(survey)) || !(all(names(survey) %in% c("participants", "contacts"))))
        {
            stop("'survey' must be either a character string or a named list with elements named 'participants' and 'contacts'")
        }
        survey_data <- survey
    }
    survey_data <- lapply(survey_data, data.table)

    if (!missing(countries))
    {
        survey_data <- lapply(survey_data, function (x) {x[get(country.column %in% countries)]})
    }

    if (missing(survey.pop))
    {
        if (missing(countries))
        {
            if (country.column %in% names(survey_data[["participants"]]))
            {
                survey.countries <- unique(survey_data[["participants"]][[country.column]])
            } else
            {
                stop("No 'survey.pop' and or 'countries' given, and no country column found in the data. I don't know which population this is from.")
            }
        } else
        {
            survey.countries <- countries
        }

        if (year.column %in% names(survey_data[["participants"]]))
        {
            survey.year <- round(mean(survey_data[["participants"]][[year.column]], na.rm = TRUE) / 5) * 5
        } else if (missing(year.column))
        {
            survey.year <- pop[, max(year)]
            warning("No year column found in the data. Will use ", survey.year, " data.")
        }

        missing.countries <- setdiff(survey.countries, survey_data[["participants"]][[country.column]])
        if (length(missing.countries) > 0)
        {
            warning("Could not find population data for ", paste(missing.countries, collapse = ", "), ". ",
                    " Use countries() to get a list of country names.")
        }

        survey.pop <- pop[country %in% survey.countries & year == survey.year][, list(population = sum(population) * 1000), by = "lower.age.limit"]
    }

    ages <- survey.pop

    if (missing(age.limits)) {
        age.limits <- unique(survey.pop$lower.age.limits)
    } else {
        ages[, lower.age.limit := reduce.agegroups(lower.age.limit, age.limits)]
        ages <- ages[, list(population = sum(population)), by = lower.age.limit]
    }
    setkey(ages, lower.age.limit)

    ret <- list()

    ## clean participant data of age NAs
    participants <- copy(data.table(survey_data[["participants"]]))
    contacts <- copy(data.table(survey_data[["contacts"]]))
    if (nrow(participants[is.na(get(part.age.column))]) > 0) {
        warning("removing participants without age")
    }
    participants <- participants[!is.na(get(part.age.column))]

    ## check maximum age in the data
    max.age <- min(max(participants[, get(part.age.column)]),
                   max(contacts[, get(contact.age.column)], na.rm = TRUE)) + 1

    ## possibly adjust age groups according to maximum age (so as not to have empty age groups)
    ages[, lower.age.limit := reduce.agegroups(lower.age.limit, lower.age.limit[lower.age.limit < max.age])]
    ages <- ages[, list(population = sum(population)), by = lower.age.limit]

    ## assign age group to participants
    participants[, lower.age.limit := reduce.agegroups(get(part.age.column), ages$lower.age.limit)]
    present.lower.age.limits <-
        participants[, .N, by = lower.age.limit][N > 1]$lower.age.limit
    present.lower.age.limits <-
        present.lower.age.limits[order(present.lower.age.limits)]

    ## reduce to all lower limits that exist in the data
    ages[, lower.age.limit := reduce.agegroups(lower.age.limit, present.lower.age.limits)]
    ages <- ages[, list(population = sum(population)), by = lower.age.limit]

    ## set upper age limits
    ages[, upper.age.limit := c(ages$lower.age.limit[-1], max.age)]

    ## reduce age groups in contacts
    contacts[, lower.age.limit := reduce.agegroups(get(part.age.column), ages$lower.age.limit)]
    contacts <- merge(contacts, ages[, list(lower.age.limit, upper.age.limit)], by = "lower.age.limit")

    participants[, agegroup := cut(participants[, get(part.age.column)],
                                   breaks = union(present.lower.age.limits, max.age),
                                   right = FALSE)]

    ## assign weights to participants, to account for weekend/weekday
    ## variation
    if (dayofweek.column %in% colnames(participants))
    {
        participants[get(dayofweek.column) %in% 1:5, weight := 5 / nrow(participants[get(dayofweek.column) %in% 1:5])]
        participants[!(get(dayofweek.column) %in% 1:5),
                    weight := 2 / nrow(participants[!(get(dayofweek.column) %in% 1:5)])]
    } else
    {
        participants[, weight := 1]
    }

    ## get number of participants in each age group
    participants.age <- unname(table(participants[, lower.age.limit]))

    if (n > 1)
    {
        if (missing(bootstrap))
        {
            bootstrap <- TRUE
        } else if (!bootstrap)
        {
            warning("n > 1 does not make sense if not bootstrapping. Will return just one sample.")
            n <- 1
        }
    }

    for (i in seq_len(n))
    {
        if (bootstrap)
        {
            ## take a bootstrap sample from the participants
            part.sample <- participants[sample(nrow(participants), replace = T)]
        } else
        {
            ## just use all participants
            part.sample <- copy(participants)
        }

        ## gather contacts for sampled participants
        contacts.sample <- data.table(merge(contacts, part.sample, by = id.column,
                                            all = F, allow.cartesian = T))

        ## sample contacts
        for (this.agegroup in unique(contacts.sample[is.na(get(contact.age.column)), agegroup]))
        {
            if (nrow(contacts.sample[!is.na(get(contact.age.column)) & agegroup == this.agegroup]) > 0)
            {
                contacts.sample[is.na(get(contact.age.column)) & agegroup == this.agegroup,
                                paste(contact.age.column) :=
                                    sample(contacts.sample[!is.na(get(contact.age.column)) & agegroup == this.agegroup,
                                                           get(contact.age.column)], size = .N, replace = TRUE)]
            } else {
                contacts.sample[is.na(get(contact.age.column)) & agegroup == this.agegroup,
                                paste(contact.age.column) := runif(.N, min = lower.age.limit, max = upper.age.limit - 1)]
            }
        }
        ## age groups
        contacts.sample[, cnt.agegroup := cut(get(contact.age.column),
                                              breaks = union(present.lower.age.limits, max.age),
                                              right = FALSE)]

        ## further weigh contacts if columns are specified
        if (length(add.weights) > 0) {
            for (i in 1:length(add.weights)) {
                contacts.sample[, weight := weight * get(add.weights)]
            }
        }

        ## calculate weighted contact matrix
        weighted.matrix <- xtabs(data = contacts.sample,
                                 formula = weight ~ cnt.agegroup + agegroup)
        ## calculate normalisation vector
        norm.vector <- xtabs(data = part.sample, formula = weight ~ agegroup)

        ## normalise contact matrix
        weighted.matrix <- t(apply(weighted.matrix, 1, function(x) { x / norm.vector} ))
        ## get rid of name but preserve row and column names
        cols <- rownames(weighted.matrix)
        weighted.matrix <- unname(weighted.matrix)

        if (symmetric & prod(dim(as.matrix(weighted.matrix))) > 1) {
            ## set C_{ij} N_j and C_{ji} N_i (which should both be equal) to
            ## 0.5 * their sum; then C_{ij} is that sum / N_j
            normalised.weighted.matrix <- t(apply(weighted.matrix, 1,
                                                  function(x) { x * ages$population }))
            weighted.matrix <- t(apply(0.5 * (normalised.weighted.matrix +
                                              t(normalised.weighted.matrix)),
                                       1, function(x) { x / ages$population }))
        }

        rownames(weighted.matrix) <- cols
        colnames(weighted.matrix) <- cols

        ret[[i]] <- list()

        if (normalise)
        {
            if (!any(is.na(weighted.matrix)))
            {
                spectrum <- eigen(weighted.matrix, only.values = TRUE)$values[1]
                weighted.matrix <- weighted.matrix / spectrum
                ret[[i]][["normalisation"]] <- spectrum
            } else
            {
                ret[[i]][["normalisation"]] <- NA_real_
            }
        }

        if (split)
        {
            contacts <- apply(weighted.matrix, 2, sum)
            age_proportions <- ages$population / sum(ages$population)
            weighted.matrix <- t(t(weighted.matrix / contacts) / age_proportions)
            ret[[i]][["contacts"]] <- contacts
        }

        ret[[i]][["matrix"]] <- weighted.matrix
    }

    if (length(ret) > 1)
        return(list(matrices = ret, demography = ages))
    else if (length(ret) == 1)
        return(list(matrix = ret[[1]], demography = ages))
    else
        stop("No matrix.")
}

##' calculates the age distribution in an endemic setting using the iterative method
##' of Wallinga (2006) based on POLYMOD data;
##'
##' this is specific to the POLYMOD data format
##' @param participants list of participants
##' @param contacts list of contacts
##' @param ages age groups
##' @return endemic age distribution as calculated by endemic.age.dist
##' @author Sebastian Funk
polymod.endemic.age.dist <- function(participants, contacts, ages, ...) {
    ## if an age range is reported, we take a random sample from that range
    sample.contacts <- data.table(contacts)
    sample.contacts[!is.na(cnt_age_l) & !is.na(cnt_age_r),
                    contact_age:= sample(cnt_age_l:cnt_age_r, 1),
                    by = rownames(sample.contacts[!is.na(cnt_age_l) &
                                                      !is.na(cnt_age_r)])]

    contact.matrix <- sample.contact.matrix(participants,
                                            sample.contacts,
                                            ages,
                                            symmetry = T)

    endemic.age.dist(contact.matrix, ...)
}

##' calculates the age distribution in an epidemic setting using the iterative method
##' of Wallinga (2006) based on POLYMOD data;
##'
##' this is specific to the POLYMOD data format
##' @param participants list of participants
##' @param contacts (list) of contacts
##' @param agegroups age groups
##' @return epidemic age distribution as calculated by epidemic.age.dist
##' @author Sebastian Funk
polymod.epidemic.age.dist <- function(participants, contacts, ages, ...) {
    ## if an age range is reported, we take a random sample from that range
    sample.contacts <- data.table(contacts)
    sample.contacts[!is.na(cnt_age_l) & !is.na(cnt_age_r),
                    contact_age:= sample(cnt_age_l:cnt_age_r, 1),
                    by = rownames(sample.contacts[!is.na(cnt_age_l) &
                                                      !is.na(cnt_age_r)])]

    contact.matrix <- sample.contact.matrix(participants,
                                            sample.contacts,
                                            ages,
                                            symmetry = T)

    epidemic.age.dist(contact.matrix, ...)
}

##' Scan a range of q values
##'
##' Plots the results
##' @param participants participants
##' @param contacts contacts
##' @param ages ages
##' @param child.mixing mixing between under 5 year olds
##' @return a (ggplot) plot with the results
##' @author Sebastian Funk
scan.age.dist <- function(participants, contacts, ages,
                          child.mixing = NA) {
    sample.contacts <- data.table(contacts)
    sample.contacts[!is.na(cnt_age_l) & !is.na(cnt_age_r),
                    contact_age:= sample(cnt_age_l:cnt_age_r, 1),
                    by = rownames(sample.contacts[!is.na(cnt_age_l) &
                                                      !is.na(cnt_age_r)])]

    contact.matrix <- sample.contact.matrix(participants,
                                            sample.contacts,
                                            ages,
                                            symmetry = T)
    if (!is.na(child.mixing)) {
        contact.matrix[1,1] <- child.mixing
    }

    ## scan a range of q values
    age.dist.q <- data.table(q = 10^(seq(from = -5, to=-1, length.out = 100)),
                             t(sapply(10^(seq(from = -5, to=-1, length.out = 100)),
                                      FUN = function(x) {
                                          endemic.age.dist(contact.matrix,
                                                           x)$normalised.y
                                      })))
    setnames(age.dist.q, 2:(length(contact.matrix$ages)+1),
             names(contact.matrix$ages))

    age.dist.q.plot <- melt(age.dist.q, id.vars="q")
    p <- ggplot(age.dist.q.plot, aes(x = q, y = value, fill = variable,
                                     color = variable))
    p <- p + geom_bar(stat = "identity")
    p <- p + scale_fill_brewer(palette="Set1")
    p <- p + scale_color_brewer(palette="Set1")
    p <- p + scale_x_log10()
    p
}

##' fits the endemic age distribution from a list of reported contacts;
##'
##' this is specific to the POLYMOD data format
##' in:
##' @param contact.matrix contact matrix and age distribution
##' @param target target relative age distribution
##' @param q transmission probability
##' @param child.mixing mixing between children
##' @param var variable to fit
##' @param rel assess goodness of fit using relative instead of absolute differences
##' @return best fitting age distribution
fit.endemic.age.dist <- function(contact.matrix, target,
                                 q = NA, child.mixing = NA, reporting = 1,
                                 vaccination = NA,
                                 var = "y",
                                 rel = F, sep = ifelse(is.na(q), F, T)) {
    absdiff <- function(x, y) {
        if (rel == T) {
            sum(abs((x - y) / y))
        } else {
            sum(abs(x - y))
        }
    }
    if (is.na(q)) {
        if (sep == F) {
            nParams <- 1
        } else {
            nParams <- length(child.mixing)
        }

        if (is.na(reporting)) {
            best <- optim(par = c(0.2, 1), function(x) {
                res <- unlist(endemic.age.dist(contact.matrix,
                                               q = 10^x[1],
                                               child.mixing = child.mixing,
                                               reporting = x[2],
                                               vaccination = vaccination)[var])
                absdiff(res, target)
                ## }, lower = c(-5, 0.2), upper = c(-1, 1))
            })

            best.fit <- c(q = 10^best$par[1], reporting = best$par[2],
                          diff = best$value,
                          unlist(endemic.age.dist(contact.matrix,
                                                  q = 10^best$par[1],
                                                  child.mixing)[[var]] * 10^best$par[2]))
            names(best.fit)[3:7] <- colnames(contact.matrix)
        } else {
            best <- optimize(function(x) {
                res <- unlist(endemic.age.dist(contact.matrix,
                                               q = 10^x,
                                               child.mixing = child.mixing,
                                               reporting = reporting,
                                               vaccination = vaccination)[var])
                absdiff(res, target)
            }, interval = c(-5, 0))

            best.fit <- c(q = 10^best$minimum, diff = best$objective,
                          unlist(endemic.age.dist(contact.matrix, 10^best$minimum,
                                                  child.mixing)[var]))
            names(best.fit)[3:7] <- colnames(contact.matrix)
        }
    } else if (is.na(reporting)) {
        best.fit <- t(sapply(1:nrow(target), function(x) {
            best <- optim(par = c(0.2, 1), function(y) {
                res <- unlist(endemic.age.dist(contact.matrix, q, 10^y[1])[var])
                absdiff(res * y[2], target[x,])
                ## }, lower = c(-5, 0.2), upper = c(1, 1))
            })

            ind.best.fit <- c(child.mixing = 10^best$par[1], reporting = best$par[2],
                              diff = best$value,
                              unlist(endemic.age.dist(contact.matrix, q,
                                                      10^best$minimum)[[var]] *
                                                          10^best$par[2]))
            names(ind.best.fit)[3:7] <- colnames(contact.matrix)
            ind.best.fit
        }))
    } else {
        best.fit <- t(sapply(1:nrow(target), function(x) {
            best <- optimize(function(y) {
                res <- unlist(endemic.age.dist(contact.matrix, q, 10^y)[var])
                absdiff(res * reporting, target[x,])
            }, interval = c(-5, 1))

            ind.best.fit <- c(child.mixing = 10^best$minimum, diff = best$objective,
                              unlist(endemic.age.dist(contact.matrix, q,
                                                      10^best$minimum)[var]))
            names(ind.best.fit)[3:7] <- colnames(contact.matrix)
            ind.best.fit
        }))
    }

    best.fit
}

##' fits the epidemic age distribution from a list of reported contacts;
##'
##' this is specific to the POLYMOD data format
##' in:
##' @param contact.matrix contact matrix and age distribution
##' @param target target relative age distribution
##' @param var variable to fit
##' @param rel assess goodness of fit using relative instead of absolute differences
##' @return best fitting age distribution
fit.epidemic.age.dist <- function(contact.matrix, target,
                                  populations, rel = F) {
    absdiff <- function(x, y) {
        if (rel == T) {
            sum(abs((x - y) / y))
        } else {
            sum(abs(x - y))
        }
    }

    best <- optimize(interval = c(-2, 0), function(x) {
        res <- unlist(epidemic.age.dist(contact.matrix,
                                        q = 10^x[1])[["z"]])
        proportions <- res * populations / sum(res * populations)
        absdiff(proportions, target)
    })

    dist <- unlist(epidemic.age.dist(contact.matrix, q = 10^best$minimum)[["z"]])

    best.fit <- c(q = 10^best$minimum,
                  diff = best$objective,
                  dist * populations / sum(dist * populations))
    names(best.fit)[3:length(best.fit)] <- colnames(contact.matrix)

    return(best.fit)
}

##' Fit POLYMOD matrix and age population to target distribution of cases
##'
##' Takes a polymod matrix (participants and contacts) and a distribution of ages in
##' the population and fits it to a target distribution of cases among ages (fitting
##' the transmission parameter q and mixing between the first age group)
##'
##' This is specific to the POLYMOD data format
##' @param participants POLYMOD participants
##' @param contacts POLIYMOD contacts
##' @param ages age groups and number of people in each of these
##' groups. This is needed to normalise the contact matrix
##' @param target target relative age distribution
##' @param nFits nubmer of fits to produce
##' @param sample whether to take a sample from the age of contacts
##' @param yearbyyear whether to fit q year-by-year instead of fitting one parameter
##' @return a list containing
##'   q: best fitting q
##'   child.mixing: best fixing mixing in first age group
##'   R0: R0 for the best fitting parameters
##'   diff: difference between the fitted and target age distributions
fit.polymod <- function(participants, contacts, ages, target, nFits = 1,
                        sample = F, yearbyyear = F, verbose = F, ...) {
    for (i in 1:nFits) {
        if (verbose) {
            cat(i, " ", format(Sys.time(), "%H:%m:%S"), "\n")
        }
        if (sample == T) {
            contacts <- contacts[!is.na(cnt_age_l) & !is.na(cnt_age_r)]
            contacts <- contacts[, id := seq(1, nrow(contacts))]
            contacts <- contacts[, contact_age := sample(cnt_age_l:cnt_age_r, 1),
                                 by = id]
        }

        contact.matrix <- sample.contact.matrix(participants, contacts, ages,
                                                symmetry = T)

        if (yearbyyear) {
            fits <- apply(target, 1, function(x) {
                q.fit <- optimize(function(y) {
                    sum(fit.endemic.age.dist(contact.matrix, t(x), 10^y, ...)[,"diff"])
                }, interval = c(-5,0))
                child.fit <-
                    fit.endemic.age.dist(contact.matrix, t(x), 10^q.fit$minimum, ...)
                c(q = 10^q.fit$minimum,
                  child.mixing = child.fit[,"child.mixing"],
                  R0 = max(Re(eigen(contact.matrix * 10^q.fit$minimum *
                                        100)$values)),
                  diff = child.fit[,"diff"])
            })
            fits <- t(fits)
        } else {
            q.fit <- optimize(function(x) {
                sum(fit.endemic.age.dist(contact.matrix, target, 10^x, ...)[,"diff"])
            }, interval = c(-5,0))
            child.fit <-
                fit.endemic.age.dist(contact.matrix, target, 10^q.fit$minimum, ...)
            new.fit <- c(q = 10^q.fit$minimum,
                         child.mixing = child.fit[,"child.mixing"],
                         R0 = max(Re(eigen(contact.matrix * 10^q.fit$minimum)$values)),
                         diff = child.fit[,"diff"])
            if (i == 1) {
                fits <- new.fit
            } else {
                fits <- rbind(fits, new.fit)
            }
        }
    }
    rownames(fits) <- NULL
    fits
}

##' Generate a homogeneous mixing matrix from a social contact matrix.
##'
##' Generates a homogeneous mixing matrix with the same age groups as
##' a given contact matrix
##' @param contact.matrix the matrix to homogenise
##' @return homogeneous mixing matrix
##' @author Sebastian Funk
homogeneous.mixing.matrix <- function(contact.matrix) {

    ## get the lower age limits from the headings of the contact matrix
    agegroups <-
        as.integer(gsub("^\\[([0-9]*),[0-9]*\\)", "\\1", colnames(contact.matrix)))

    ## get populations by age group
    ages <- pop.ew.age[year == 2006]
    ages[, lower.age.limit := reduce.agegroups(lower.age.limit, agegroups)]
    ages <- ages[, list(population = sum(population)), by = lower.age.limit]

    hom.matrix <-
        matrix(rep(ages[, population] / ages[, sum(population)],
                   each = ncol(contact.matrix)),
               ncol = ncol(contact.matrix))

    colnames(hom.matrix) <- colnames(contact.matrix)
    rownames(hom.matrix) <- rownames(contact.matrix)

    return(hom.matrix)

}
