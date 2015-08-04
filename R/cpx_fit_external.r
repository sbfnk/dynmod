cpx.pomp <- function(parameters) {

    cpx.skel <- function(x, t, params, ...) {
        with(
            as.list(c(x, params)), {
                beta <- R0 / gamma
                current.mixing <- mixing
                if ((t %% 52) > holiday.start && (t %% 52) <= holiday.end) {
                    current.mixing <- sigma * mixing
                }
                c(S = births - beta * current.mixing * S * I/N - mu * S,
                  E = beta * current.mixing * S * I/N - rho * E - mu * E,
                  I = rho * E - gamma * I - mu * I,
                  R = gamma * I - mu * R,
                  Z = rho * E
                  )
            })
    }

    cpx.initializer <- function(params, t0, ...) {
        with(
            as.list(params), {
                cat(params, "\n")
                S.0 <- (N * (mu + gamma) * (mu + rho)) / (R0 * gamma * rho)
                I.0 <- (N * (births - mu * S.0)) / (gamma * R0 * S.0)
                E.0 <- (I.0 * (mu + gamma)) / rho
                R.0 <- (I.0 * gamma) / mu
                beta <- R0 / gamma
                holiday.mixing <- sigma
                return(c(S = S.0, I = I.0, E = E.0, R = R.0, Z = 0))
            })
    }

    cpx.measure <- function (y, x, t, params, log, ...) {
        ## state at time t:
        X <- x["X"]
        ## observation at time t:
        Y <- y["Y"]
        ## compute the likelihood of Y|X,tau
        f <- dnorm(x=Y, mean=X, sd=1 + X)
        return(f)
    }
    R0 <- parameters[1]
    sigma <- parameters[2]
    N <- 45e+6
    infper <- 10.5 / 7
    latper <- 10 / 7
    annual.births <- 796858 # from Bryan's data

    annual.mortality <- annual.births / N

    gamma <- 1 / infper
    rho <- 1 / latper
    mu <- annual.mortality / 52.5
    births <- annual.births / 52.5

    params <- c(N = N,
                gamma = gamma,
                rho = rho,
                mu = mu,
                births = births,
                holiday.start = 30,
                holiday.end = 36,
                mixing = 1,
                R0 = R0,
                sigma = sigma
                )

    target <- readRDS("cpx_target_1967.rds")
    target[, week := week + 20800]
    my.pomp <- pomp(data = target,
                    times = "week",
                    t0 = 0,
                    dmeasure = cpx.measure,
                    skeleton.type = "vectorfield",
                    skeleton = cpx.skel,
                    params = params,
                    initializer = cpx.initializer,
                    comp.names = c("S", "E", "I", "R", "Z"),
                    ic.names = c("S.0", "E.0", "I.0", "R.0", "Z.0"),
                    zeronames = "Z")
    return(my.pomp)
}

cpx.pomp.age <- function(parameters) {
    library(data.table)

    cpx.skel <- function(x, t, params, ...) {
        with(
            as.list(c(x, params)), {
                beta <- R0 / gamma
                current.mixing <- mixing
                if ((t %% 52) > holiday.start && (t %% 52) <= holiday.end) {
                    current.mixing <- sigma * mixing
                }
                c(S = births - beta * current.mixing * S * I/N - mu * S,
                  E = beta * current.mixing * S * I/N - rho * E - mu * E,
                  I = rho * E - gamma * I - mu * I,
                  R = gamma * I - mu * R,
                  Z = rho * E
                  )
            })
    }

    cpx.initializer <- function(params, t0, ...) {
        with(
            as.list(params), {
                cat(params, "\n")
                S.0 <- (N * (mu + gamma) * (mu + rho)) / (R0 * gamma * rho)
                I.0 <- (N * (births - mu * S.0)) / (gamma * R0 * S.0)
                E.0 <- (I.0 * (mu + gamma)) / rho
                R.0 <- (I.0 * gamma) / mu
                beta <- R0 / gamma
                holiday.mixing <- sigma
                return(c(S = S.0, I = I.0, E = E.0, R = R.0, Z = 0))
            })
    }

    cpx.measure <- function (y, x, t, params, log, ...) {
        ## state at time t:
        X <- x["X"]
        ## observation at time t:
        Y <- y["Y"]
        ## compute the likelihood of Y|X,tau
        f <- dnorm(x=Y, mean=X, sd=1 + X)
        return(f)
    }
    R0 <- parameters[1]
    sigma <- parameters[2]
    N <- 45e+6
    infper <- 10.5 / 7
    latper <- 10 / 7
    annual.births <- 796858 # from Bryan's data

    annual.mortality <- annual.births / N

    gamma <- 1 / infper
    rho <- 1 / latper
    mu <- annual.mortality / 52.5
    births <- annual.births / 52.5

    params <- c(N = N,
                gamma = gamma,
                rho = rho,
                mu = mu,
                births = births,
                holiday.start = 30,
                holiday.end = 36,
                mixing = 1,
                R0 = R0,
                sigma = sigma
                )

    target <- readRDS("cpx_target_1967.rds")
    target[, week := week + 20800]
    my.pomp <- pomp(data = target,
                    times = "week",
                    t0 = 0,
                    dmeasure = cpx.measure,
                    skeleton.type = "vectorfield",
                    skeleton = cpx.skel,
                    params = params,
                    initializer = cpx.initializer,
                    comp.names = c("S", "E", "I", "R", "Z"),
                    ic.names = c("S.0", "E.0", "I.0", "R.0", "Z.0"),
                    zeronames = "Z")
    return(my.pomp)
}

