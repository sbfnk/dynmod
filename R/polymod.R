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
