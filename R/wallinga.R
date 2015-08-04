##' calculates the age distribution in an endemic setting using the iterative method
##' of Wallinga (2006)
##'
##' calculates the age distribution in an endemic setting using the iterative method
##' of Wallinga (2006)
##' @param contact.matrix: social contact matrix; this must have the
##' age groups as column names (with a comma between limits)
##' @param q: vector of probabilies of transmission per contact
##' @param child.mixing: vector of first entries of contact matrix
##' @param reporting: vector of reporting probabilities
##' @param vaccination: proportion vaccinated (between 0 and 1) in each age group
##' @param incidence.start: starting value for inidence
##' @param tol: tolerance for stopping the iteration
##' @return the endemic age distribution, a list composed of:
##'   y: incidence in each age group
##'   normalised.y: proportion of cases in each age group
##'   lambda: force of infection by age group
##'   prop.immune: proportion immune in each age group
##'   prop.immune.age: proportion immune at each age
endemic.age.dist <- function(contact.matrix, q,
                             child.mixing = NA,
                             reporting = NA,
                             vaccination = NA,
                             incidence.start = 0.01,
                             tol = 1e-5) {

  l <- strsplit(gsub("[^0-9,]", "", colnames(contact.matrix)), ",")
  ## get the sizes of the age compartments from the names of the age groups
  age.ranges <- unlist(lapply(l, function(x)
                          { as.numeric(x[2]) - as.numeric(x[1]) }))

  first.run <- T

  y <- rep(incidence.start, nrow(contact.matrix))

  ## set to greater than tol for first time the loop is run
  current.diff <- tol + 1
  nIter <- 0

  if (!is.na(child.mixing)) {
      contact.matrix[1,1] <- child.mixing
  }

  if (any(is.na(vaccination))) {
      vaccination <- rep(0, length(age.ranges))
  }

  if (any(is.na(reporting))) {
      reporting <- rep(1, length(age.ranges))
  }

  ## loop until difference between estimates is smaller than tolerance
  while (current.diff > tol) {
    nIter <- nIter + 1
    lambda <- y %*% (q * contact.matrix)
    prop.immune <- (1 - exp(-cumsum(lambda * age.ranges))) * (1 - vaccination)
    y <- (prop.immune - c(0, prop.immune[-length(prop.immune)])) /
        (age.ranges)

    ## correction for implausible situations because more people are
    ## vaccinated in an older age group
    y[y<0] <- 0
    if (first.run) {
      ## run loop at least two times
      current.diff <- tol + 1
      first.run <- F
    } else {
      current.diff <- sum(abs(last.y-y))
    }
    last.y <- y
  }
  lambda.age <- c()
  for (i in 1:length(lambda)) {
    lambda.age <- c(lambda.age, rep(lambda[i], age.ranges[i]))
  }
  list(y = y * reporting, true.y = y, lambda = lambda, prop.immune = prop.immune,
       iteration = nIter)

}

##' calculates the age distribution in an epidemic setting using the iterative method
##' of Wallinga (2006)
##'
##' calculates the age distribution in an epidemic setting using the iterative method
##' of Wallinga (2006)
##' @param contact.matrix: social contact matrix
##' @param age.sizes: sizes of (number of people in) the different age groups
##' @param q: probability of transmission per contact
##' @param child.mixing: first entry of contact matrix
##' @param immunity: proportion immune before the epidemic
##' @param incidence.start: starting value for inidence
##' @param tol: tolerance for stopping the iteration
##' @export
##' @return A list composed of
##'   z: final size in each age group
##'   normalised.z: proportion of cases in each age group
epidemic.age.dist <- function(contact.matrix, q,
                              child.mixing = NA,
                              immunity = 0,
                              final.size.start = 0.01,

                              tol = 1e-5) {
  ## initialise variables
  z <- rep(final.size.start, nrow(contact.matrix))
  last.z <- rep(0, nrow(contact.matrix))
  first.run <- T
  if (!is.na(child.mixing)) {
    contact.matrix[1,1] <- child.mixing
  }
  if (length(immunity) == 1) {
    immunity <- rep(immunity, nrow(contact.matrix))
  }
  contact.matrix <- t(apply(contact.matrix, 1, function(x) {
    x * (1 - immunity)
  }))

  ## set to greater than tol for first time the loop is run
  current.diff <- tol + 1

  ## loop until difference between estimates is smaller than tolerance
  while (current.diff > tol) {
    z <- 1 - exp(- z %*% (q * 100 * contact.matrix))
    last.z <- z
    if (first.run == T) {
      ## run loop at least two times
      current.diff <- tol + 1
      first.run <- F
    } else {
      current.diff <- sum(abs(last.z-z))
    }
  }
  list(z = z, normalised.z = z / sum(z))
}
