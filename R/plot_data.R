##' Plot England & Wales measles data
##'
##' @param ext plot extensions
##' @author Sebastian Funk
##' @import ggplot2 Hmisc reshape2
plotMeaslesEWData <- function(ext = NA) {

    data(ms_ew_age)

    ms.ages <-
      ms.ew.age[, list(mean = wtd.mean(lower.age.limit, abs.incidence),
                       low = wtd.quantile(lower.age.limit, abs.incidence,
                         c(0.1)),
                       high = wtd.quantile(lower.age.limit, abs.incidence,
                         c(0.9)),
                       sd = sqrt(wtd.var(lower.age.limit, abs.incidence))),
                by = year]
    ms.ages.m <- melt(ms.ages, measure.vars = c("mean", "sd"))
    p <- ggplot(ms.ages.m, aes(x = year, y = value, color = variable))
    p <- p + expand_limits(y=0)
    p <- p + geom_line(lwd = 1.5)
    p <- p + scale_color_brewer("", palette="Set1",
                                labels = c("mean", "standard deviation"))
    p <- p + facet_grid(variable ~ ., scales = "free")
    p <- p + scale_y_continuous("Age")
    p <- p + theme_bw(20)
    p <- p + theme(legend.position = "bottom")
    p <- p + ggtitle("Measles")
    p <- p + geom_vline(xintercept=1997.5)
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("measles_annual_meanage.", ext, sep = ""), p)
    }

    p <- ggplot(ms.ew.age, aes(x = year, y = rel.incidence,
                               colour = factor(lower.age.limit)))
    p <- p + geom_line(lwd = 1.5)
    p <- p + theme_bw(20)
    p <- p + scale_color_brewer("age group", palette="Set1",
                  labels = limits.to.agegroups(ms.ew.age[, lower.age.limit]))
    p <- p + scale_y_continuous("Incidence per 100,000")
    p <- p + theme(legend.position="bottom")
    p <- p + ggtitle("Measles")
    if (is.na(ext)) {
        quartz(width = 15)
        p
    } else {
        ggsave(paste("measles_annual_incidence_age.", ext, sep = ""),
               p, width = 15)
    }

    p <- ggplot(ms.ew.age[year > 1997], aes(x = year, y = rel.incidence,
                                           colour = factor(lower.age.limit)))
    p <- p + geom_line(lwd = 1.5)
    p <- p + theme_bw(16)
    p <- p + scale_color_brewer("age group", palette="Set1",
                  labels = limits.to.agegroups(ms.ew.age[, lower.age.limit]))
    p <- p + scale_y_continuous("Incidence")
    p <- p + theme(legend.position="bottom")
    p <- p + ggtitle("Measles")
    if (is.na(ext)) {
        quartz()
        p
    } else {
        ggsave(paste("measles_annual_incidence_age_post1998.",
                     ext, sep = ""), p)
    }
}

##' Plot European measles maps
##'
##' @param ext plot extensions
##' @author Sebastian Funk
##' @import ggplot2 
plotMeaslesEuropeanMaps <- function(ext) {

    try_require("data.table")
    try_require("ggplot2")

    ## get map of europe
    data(europe)

    ## population
    data(pop_europe)

    ## put cases on the map
    data(ms_europe_cases)

    ms.europe.cases.map <- ms.europe.cases
    ms.europe.cases.map <-
        ms.europe.cases.map[, year.week := paste(year, sprintf("%02i", week),
                                        sep = "-")]
    ms.europe.cases.map <-
        ms.europe.cases.map[, year.month := paste(year, sprintf("%02i", month),
                                         sep = "-")]

    ## plot.scales <- c("year", "year.month", "year.week")
    plot.scales <- c("year.week")
    for (scale in plot.scales) {

        counts <- ms.europe.cases.map[, list(cases = length(date)),
                                      by = c(scale, "place")]
        setkey(counts, "place")
        map.data <- merge(counts, europe, allow.cartesian = TRUE)
        map.data <- map.data[, paste(scale) := as.character(get(scale))]

        setnames(map.data, "STAT_LEVL_", "stat.level")
        setkeyv(map.data, c("place", scale))

        ## set the population size for earlier years from what's
        ## available
        pop.europe <- pop.europe[year == min(year),
                                 year := min(ms.europe.cases[, year])]
        pop.europe <- pop.europe[, paste(scale) := as.character(year)]
        setkeyv(pop.europe, c("place", scale))
        map.data <- pop.europe[map.data, roll = T]

        for (i in unique(map.data[, stat.level])) {
            max.cases <- max(map.data[stat.level == i, cases])
            map.data <- map.data[stat.level == i,
                                 scaled.cases := cases/max.cases]
        }

        setkeyv(map.data, scale)
        dates <- unique(map.data[, scale, with = F])

        for (date in unname(unlist(dates))) {
            current.map <- map.data[get(scale) == date]
            map <- ggplot(current.map, aes(long, lat, group = group,
                                           fill = log10(cases / population)))
            map <- map + geom_polygon(data = europe, fill="grey")
            map <- map + geom_polygon()
            map <- map + scale_fill_continuous(low = "grey", high = "darkred",
                                               limits = c(-8,-2.9))
            map <- map + coord_equal()
            map <- map + scale_x_continuous(limits = c(-25,45))
            map <- map + scale_y_continuous(limits = c(35,70))
            map <- map + theme_bw()
            map <- map + theme(legend.position="none")
            map <- map + ggtitle(as.character(date))
            ggsave(filename = paste("europe_", gsub("\\.", "_", scale),
                   sprintf("_%s", date), ".", ext, sep=""),
                   plot = map, width = 6, height = 3)
        }
    }
}

##' Plot England & Wales chickenpox data
##'
##' @param ext plot extensions
##' @author Sebastian Funk
##' @import ggplot2 reshape2 Hmisc lubridate
plotChickenpoxEWData <- function(ext = NA) {

    data(cpx_ew)

    ## convert week/year to date
    cpx.ew <- cpx.ew[, date := year.week.to.date(year, week)]

    ## number of cases each year
    p <- ggplot(cpx.ew, aes(x = date, y = rel.incidence))
    p <- p + geom_line(lwd = 1)
    p <- p + theme_bw(20)
    p <- p + scale_y_continuous("Incidence / 100,000")
    p <- p + ggtitle("Chickenpox")
    if (is.na(ext)) {
        quartz()
        p
    } else {
        ggsave(paste("cpx_incidence.", ext, sep = ""), p, width = 15)
    }

    data(cpx_ew_age)

    ## convert week/year to date
    cpx.ew.age <- cpx.ew.age[, date := year.week.to.date(year, week)]

    p <- ggplot(cpx.ew.age, aes(x = date, y = rel.incidence,
                                colour = factor(lower.age.limit)))
    p <- p + geom_line()
    p <- p + theme_bw(20)
    p <- p + scale_color_brewer("age group", palette="Set1",
                  labels = limits.to.agegroups(cpx.ew.age[, lower.age.limit]))
    p <- p + scale_y_continuous("Incidence per 100,000")
    p <- p + theme(legend.position="bottom")
    p <- p + ggtitle("Chickenpox")
    if (is.na(ext)) {
        quartz(width = 15)
        p
    } else {
        ggsave(paste("cpx_incidence_age.", ext, sep = ""), p, width = 15)
    }

    cpx.ew.age.total <- cpx.ew.age[, list(abs.pop.incidence = sum(abs.incidence),
                                          population = sum(population)),
                                   by = list(year, week)]
    cpx.ew.age.total <-
        cpx.ew.age.total[, rel.pop.incidence := abs.incidence / population]
    setkey(cpx.ew.age.total, year, week)
    cpx.ew.age <- merge(cpx.ew.age, cpx.ew.age.total, by = "year")
    p <- ggplot(cpx.ew.age, aes(x = date, y = abs.pop.incidence / pop.incidence,
                                colour = factor(lower.age.limit)))
    p <- p + geom_line()
    p <- p + theme_bw(20)
    p <- p + scale_color_brewer("age group", palette="Set1",
                  labels = limits.to.agegroups(cpx.ew.age[, lower.age.limit]))
    p <- p + scale_y_continuous("Proportion of cases", limits = c(0, 1))
    p <- p + theme(legend.position="bottom")
    p <- p + ggtitle("Chickenpox")
    if (is.na(ext)) {
        quartz(width = 15)
        p
    } else {
        ggsave(paste("cpx_proportion_age.", ext, sep = ""), p, width = 15)
    }

    cpx.ew.age.annual <-
        cpx.ew.age[, list(rel.incidence = sum(rel.incidence),
                          abs.incidence = sum(abs.incidence)),
                   by = list(year, lower.age.limit)]

    p <- ggplot(cpx.ew.age.annual, aes(x = year, y = rel.incidence,
                                       color = factor(lower.age.limit)))
    p <- p + geom_line(lwd = 1.5)
    p <- p + facet_grid(lower.age.limit ~ ., scales = "free")
    p <- p + expand_limits(y=0)
    p <- p + scale_color_brewer("age group", palette="Set1",
                  labels = limits.to.agegroups(.ew.age[, lower.age.limit]))
    p <- p + scale_y_continuous("Incidence per 100,000")
    p <- p + theme_bw(20)
    p <- p + theme(legend.position = "bottom")
    p <- p + ggtitle("Chickenpox")
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("cpx_annual_age.", ext, sep = ""), p)
    }

    cpx.ew.age.annual.total <-
        cpx.ew.age.annual[, list(pop.incidence = sum(abs.incidence)), by = year]
    setkey(cpx.ew.age.annual.total, year)
    setkey(cpx.ew.age.annual, year)
    cpx.ew.age.annual <- merge(cpx.ew.age.annual, cpx.ew.age.annual.total, by = "year")
    p <- ggplot(cpx.ew.age.annual, aes(x = year, y = abs.incidence / pop.incidence,
                                       colour = factor(lower.age.limit)))
    p <- p + geom_line(lwd = 1.5)
    p <- p + theme_bw(20)
    p <- p + scale_color_brewer("age group", palette="Set1",
                  labels = limits.to.agegroups(cpx.ew.age[, lower.age.limit]))
    p <- p + scale_y_continuous("Proportion of reported cases", limits = c(0, 1))
    p <- p + theme(legend.position="bottom")
    p <- p + ggtitle("Chickenpox")
    if (is.na(ext)) {
        quartz(width = 15)
        p
    } else {
        ggsave(paste("cpx_annual_proportion_age.", ext, sep = ""), p)
    }

    cpx.ages <-
        cpx.ew.age[, list(mean = wtd.mean(lower.age.limit, abs.incidence),
                          low = wtd.quantile(lower.age.limit, abs.incidence,
                          c(0.1)),
                          high = wtd.quantile(lower.age.limit, abs.incidence,
                          c(0.9)),
                          sd = sqrt(wtd.var(lower.age.limit, rel.incidence))),
                   by = year]

    cpx.ages.m <- melt(cpx.ages, measure.vars = c("mean", "sd"))
    p <- ggplot(cpx.ages.m, aes(x = year, y = value, color = variable))
    p <- p + expand_limits(y=0)
    p <- p + geom_line(lwd = 1.5)
    p <- p + scale_color_brewer("", palette="Set1",
                                labels = c("mean", "standard deviation"))
    p <- p + facet_grid(variable ~ ., scales = "free")
    p <- p + scale_y_continuous("Age")
    p <- p + theme_bw(20)
    p <- p + theme(legend.position = "bottom")
    p <- p + ggtitle("Chickenpox")
    if (is.na(ext)) {
        quartz()
        p
    } else {
        ggsave(paste("cpx_annual_meanage.", ext, sep = ""), p)
    }

    data(cpx_ew)

    temp.cpx.ew.acf <- acf(cpx.ew$rel.incidence, lag.max=100, plot = F)
    cpx.ew.acf <- data.table(autocorrelation = temp.cpx.ew.acf$acf,
                             lag = temp.cpx.ew.acf$lag)
    temp.cpx.ew.spectrum <- spectrum(cpx.ew$rel.incidence, plot = F)
    cpx.ew.spectrum <- data.table(spectrum = temp.cpx.ew.spectrum$spec,
                                  frequency = temp.cpx.ew.spectrum$freq)

    p <- ggplot(data.frame(cpx.ew.acf), aes(x=lag, y=autocorrelation))
    p <- p + geom_line()
    p <- p + theme_bw(20)
    p <- p + scale_x_continuous("lag (in weeks)")
    p <- p + ggtitle("Chickenpox")
    if (is.na(ext)) {
        quartz()
        p
    } else {
        ggsave(paste("cpx_autocorrelations.", ext, sep = ""), p)
    }

    p <- ggplot(data.frame(cpx.ew.spectrum), aes(x = frequency, y = spectrum))
    p <- p + geom_line()
    p <- p + theme_bw(20)
    p <- p + ggtitle("Chickenpox")
    if (is.na(ext)) {
        quartz()
        p
    } else {
        ggsave(paste("cpx_spectrum.", ext, sep = ""), p)
    }

    p <- ggplot(data.frame(cpx.ew.spectrum), aes(x = 1 / frequency, y = spectrum))
    p <- p + geom_line()
    p <- p + theme_bw(20)
    p <- p + scale_x_continuous("period", limits=c(1, 150))
    p <- p + scale_y_continuous("amplitude")
    p <- p + ggtitle("Chickenpox")
    if (is.na(ext)) {
        quartz()
        p
    } else {
        ggsave(paste("cpx_period.", ext, sep = ""), p)
    }

    incidence.by.week <- cpx.ew[, list(rel.incidence = mean(rel.incidence),
                                       sd = sd(rel.incidence)),
                                by = list(week)]

    p <- ggplot(incidence.by.week, aes(x = week, y = rel.incidence,
                                       ymin = rel.incidence - sd,
                                       ymax = rel.incidence + sd))
    p <- p + geom_point()
    p <- p + geom_errorbar()
    p <- p + theme_bw(20)
    p <- p + coord_cartesian(ylim = c(0,
                             max(incidence.by.week[, rel.incidence + sd]) * 1.1))
    p <- p + scale_y_continuous("Mean incidence")
    if (is.na(ext)) {
        quartz()
        p
    } else {
        ggsave(paste("cpx_weeks.", ext, sep = ""), p)
    }

    incidence.by.week.shifted <- copy(incidence.by.week)
    incidence.by.week.shifted <-
        incidence.by.week.shifted[week >= 36, week := as.integer(week - 52)]
    p <- ggplot(incidence.by.week.shifted, aes(x = week, y = rel.incidence,
                                               ymin = rel.incidence - sd,
                                               ymax = rel.incidence + sd))
    p <- p + geom_point()
    p <- p + geom_errorbar()
    p <- p + theme_bw(20)
    p <- p + coord_cartesian(ylim = c(0,
                             max(incidence.by.week[, rel.incidence + sd]) * 1.1))
    p <- p + scale_y_continuous("Mean incidence")
    if (is.na(ext)) {
        quartz()
        p
    } else {
        ggsave(paste("cpx_weeks_shifted.", ext, sep = ""), p)
    }

    cpx.ew[, era := cut(cpx.ew[, year],
                 breaks = c(seq(1966, 1996, by = 10), 2014),
                 labels = c("pre-1977", "1977-86", "1987-96", "post-1996"))]

    incidence.by.week.era <- cpx.ew[, list(rel.incidence = mean(rel.incidence),
                                           sd = sd(rel.incidence)),
                                    by = list(week, era)]

    p <- ggplot(incidence.by.week.era, aes(x = week, y = rel.incidence,
                                           ymin = rel.incidence - sd,
                                           ymax = rel.incidence + sd,
                                           color = era))
    p <- p + geom_point()
    p <- p + geom_errorbar()
    p <- p + theme_bw(20)
    p <- p + coord_cartesian(ylim = c(0,
                             max(incidence.by.week.era[, rel.incidence + sd]) * 1.1))
    p <- p + scale_y_continuous("Mean incidence")
    p <- p + scale_color_brewer(palette = "Set1")
    p <- p + theme(legend.position = "bottom")
    p <- p + facet_grid(era ~ .)
    if (is.na(ext)) {
        quartz()
        p
    } else {
        ggsave(paste("cpx_weeks_era.", ext, sep = ""), p)
    }

    incidence.by.week.era.shifted <- copy(incidence.by.week.era)
    incidence.by.week.era.shifted <-
        incidence.by.week.era.shifted[week >= 36, week := as.integer(week - 52)]
    p <- ggplot(incidence.by.week.era.shifted, aes(x = week, y = rel.incidence,
                                                   ymin = rel.incidence - sd,
                                                   ymax = rel.incidence + sd,
                                                   color = era))
    p <- p + geom_point()
    p <- p + geom_errorbar()
    p <- p + theme_bw(20)
    p <- p + coord_cartesian(ylim = c(0,
                             max(incidence.by.week.era[, rel.incidence + sd]) * 1.1))
    p <- p + scale_y_continuous("Mean incidence")
    p <- p + scale_color_brewer(palette = "Set1")
    p <- p + theme(legend.position = "bottom")
    p <- p + facet_grid(era ~ .)
    if (is.na(ext)) {
        quartz()
        p
    } else {
        ggsave(paste("cpx_weeks_era_shifted.", ext, sep = ""), p)
    }

}

##' Plot chickenpox serology
##'
##' @param ext plot extensions
##' @author Sebastian Funk
##' @import ggplot2
plotChickenpoxSerology <- function(ext = NA) {
    data(cpx_sero)

    try_require("ggplot2")

    p <- ggplot(cpx.sero, aes(x = year, y = antibodies / samples,
                              fill = factor(lower.age.limit)))
    p <- p + geom_bar(stat = "identity", position = "dodge")
    p <- p + scale_fill_brewer("age group", palette="Set1",
                  labels = limits.to.agegroups(cpx.sero[, lower.age.limit]))
    p <- p + theme_bw(16)
    p <- p + theme(legend.position = "bottom")
    p <- p + scale_y_continuous("proportion with antibodies")
    p <- p + ggtitle("Chickenpox")
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("chickenpox_serology.", ext, sep = ""), p)
    }
}

##' Plot England and Wales measles and chicknepox data
##'
##' @param ext plot extensions
##' @author Sebastian Funk
##' @import ggplot2
plotMSCPXEWData <- function(ext = NA) {
    data(cpx_ew_age)
    data(ms_ew_age)

    ## find common lower age limits
    common.lower.limits <- intersect(unique(cpx.ew.age[, lower.age.limit]),
                                     unique(ms.ew.age[, lower.age.limit]))

    ## find in which interval of common lower age limits the chickenpox
    ## lower age limits are
    cpx.intervals <- findInterval(cpx.ew.age[, lower.age.limit],
                                  common.lower.limits)
    ## convert them to limits
    lower.limits <- common.lower.limits[cpx.intervals]
    ## assign them
    cpx.ew.age <- cpx.ew.age[, lower.age.limit := lower.limits]
    ## synthesise them
    cpx.ew.age <- cpx.ew.age[, list(abs.incidence = sum(abs.incidence)),
                             by = list(year, lower.age.limit)]

    ## find in which interval of common lower age limits the chickenpox
    ## lower age limits are
    ms.intervals <- findInterval(ms.ew.age[, lower.age.limit],
                                 common.lower.limits)
    ## convert them to limits
    lower.limits <- common.lower.limits[ms.intervals]
    ## assign them
    ms.ew.age <- ms.ew.age[, lower.age.limit := lower.limits]
    ## synthesise them
    ms.ew.age <- ms.ew.age[, list(abs.incidence = sum(abs.incidence)),
                           by = list(year, lower.age.limit)]

    ## merge the datasets
    cpx.ew.age <- cpx.ew.age[, disease := "chickenpox"]
    ms.ew.age <- ms.ew.age[, disease := "measles"]
    cpx.ms.ew.age <- rbind(cpx.ew.age, ms.ew.age)

    ## calculate proportion in different age groups
    annual.sums <-
        cpx.ms.ew.age[, list(total.abs.incidence = sum(abs.incidence)),
                      by = list(year, disease)]
    setkey(cpx.ms.ew.age, year, disease)
    setkey(annual.sums, year, disease)
    cpx.ms.ew.age <- merge(cpx.ms.ew.age, annual.sums)

    p <- ggplot(cpx.ms.ew.age,
                aes(x = year, y = abs.incidence / total.abs.incidence,
                    color = factor(lower.age.limit),
                    fill = factor(lower.age.limit)))
    p <- p + geom_bar(stat = "identity")
    p <- p + facet_grid(disease~., scales = "fixed")
    p <- p + expand_limits(y=0)
    p <- p + theme_bw(16)
    p <- p + scale_color_brewer(palette = "Set1", name = "age group",
                               labels = c("0-4", "5-14", "15+"))
    p <- p + scale_fill_brewer(palette = "Set1", name = "age group",
                               labels = c("0-4", "5-14", "15+"))
    p <- p + scale_y_continuous("proportion of cases in age group", limits=c(0,1))
    p <- p + theme(legend.position = "bottom")
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("measles_chickenpox_ages.", ext, sep = ""), p)
    }
}

##' Plot England & Wales measles confirmation data
##'
##' @param ext plot extensions
##' @author Sebastian Funk
##' @import ggplot2
plotMeaslesEWConfData <- function(ext = NA) {

    data(ms_ew_conf)
    ## convert year/quarter to date
    ms.ew.conf <-
        ms.ew.conf[, date := as.Date(paste(year, sprintf("%02i",
                          as.numeric(substr(quarter, 1, 1)) * 3), "01",
                          sep = "-"))]

    p <- ggplot(ms.ew.conf[problem != "x"], aes(x = uncorrected.notified.cases,
                   y = percent.confirmed))
    p <- p + geom_jitter()
    p <- p + theme_bw(16)
    p <- p + scale_x_continuous("Number of cases notified")
    p <- p + scale_y_continuous("Percent of cases tested that were confirmed")
    p <- p + ggtitle("Measles in England & Wales")
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("measles_ew_confirmations.", ext, sep = ""), p)
    }

    p <- ggplot(ms.ew.conf[problem != "x"], aes(x = date,
                   y = percent.confirmed))
    p <- p + geom_jitter()
    p <- p + theme_bw(16)
    p <- p + scale_x_date("Year")
    p <- p + scale_y_continuous("Percent of cases tested that were confirmed")
    p <- p + ggtitle("Measles in England & Wales")
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("measles_ew_confirmations.", ext, sep = ""), p)
    }
}

##' Plot WHO measles data
##'
##' @param ext plot extensions
##' @author Sebastian Funk
##' @import ggplot2
plotMeaslesWHOData <- function(ext = NA) {

    data(ms_who_age)

    ms.turkey.age <-
        ms.who.age[country == "Turkey" & week.date >= "2012-10-01"]
    setkey(ms.turkey.age, week.date, year, week, lower.age.limit)

    p <- ggplot(ms.turkey.age, aes(x = week.date, y = rel.incidence,
                                colour = factor(lower.age.limit)))
    p <- p + geom_line(lwd = 1.5)
    p <- p + theme_bw(20)
    p <- p + scale_color_brewer("age group", palette="Set1",
                  labels = limits.to.agegroups(ms.turkey.age[, lower.age.limit]))
    p <- p + scale_y_continuous("Incidence per 100,000")
    p <- p + theme(legend.position="bottom")
    p <- p + ggtitle("Measles (Turkey)")
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("ms_turkey_incidence_age.", ext, sep = ""), p)
    }

    ms.georgia.age <-
        ms.who.age[country == "Georgia" & week.date >= "2013-01-01"]
    setkey(ms.georgia.age, week.date, year, week, lower.age.limit)

    p <- ggplot(ms.georgia.age, aes(x = week.date, y = rel.incidence,
                                colour = factor(lower.age.limit)))
    p <- p + geom_line(lwd = 1.5)
    p <- p + theme_bw(20)
    p <- p + scale_color_brewer("age group", palette="Set1",
                  labels = limits.to.agegroups(ms.georgia.age[, lower.age.limit]))
    p <- p + scale_y_continuous("Incidence per 100,000")
    p <- p + theme(legend.position="bottom")
    p <- p + ggtitle("Measles (Georgia)")
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("ms_georgia_incidence_age.", ext, sep = ""), p)
    }

    ms.bulgaria.age <-
        ms.who.age[country == "Bulgaria" &
                   week.date >= "2009-01-01" & week.date <= "2011-06-01"]
    setkey(ms.bulgaria.age, week.date, year, week, lower.age.limit)

    p <- ggplot(ms.bulgaria.age, aes(x = week.date, y = rel.incidence,
                                colour = factor(lower.age.limit)))
    p <- p + geom_line(lwd = 1.5)
    p <- p + theme_bw(20)
    p <- p + scale_color_brewer("age group", palette="Set1",
                  labels = limits.to.agegroups(ms.bulgaria.age[, lower.age.limit]))
    p <- p + scale_y_continuous("Incidence per 100,000")
    p <- p + theme(legend.position="bottom")
    p <- p + ggtitle("Measles (Bulgaria)")
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("ms_bulgaria_incidence_age.", ext, sep = ""), p)
    }
}

##' Plot measles serology
##'
##' @param ext plot extensions
##' @author Sebastian Funk
##' @import ggplot2
plotMeaslesSerology <- function(ext = NA) {

    data(ms_sero)

    ms.sero.age <-
        ms.sero[, list(non.negative = mean(non.negative)),
                by = list(country, age1)]
}

##' Plot MMR data
##'
##' @param ext plot extensions
##' @author Sebastian Funk
##' @import ggplot2
plotMMRData <- function(ext = NA) {

    data(vaccine_ew)

    mvac <- melt(vaccine.ew, id.vars = "year")
    mvac12 <- mvac[variable %in% c("MCV1", "MCV2")]
    mvac12 <- mvac12[, variable := factor(variable)]
    p <- ggplot(mvac12, aes(x = year, y = value, color = variable))
    p <- p + geom_line(lwd = 1.5)
    p <- p + theme_bw(20)
    p <- p + scale_color_brewer("dose", palette="Set1",
                                labels = c("1st (1y)", "2nd (5y)"))
    p <- p + scale_y_continuous("Uptake")
    p <- p + theme(legend.position="bottom")
    p <- p + ggtitle("MMR")
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("mmr_uptake.", ext, sep = ""), p)
    }

}

##' Plot England & Wales chickenpox before/after ratios
##'
##' @param ext plot extensions
##' @author Sebastian Funk
##' @import ggplot2 zoo
plotChickenpoxEWBeforeAfterRatios <- function(ext = NA, lines = T) {

    data(cpx_ew)

    start.year <- 40

    cpx.ew <- cpx.ew[week >= start.year, year := as.integer(year + 1)]
    cpx.ew <- cpx.ew[week >= start.year, week := as.integer(week - 52)]
    cpx.ew <- cpx.ew[, week := week + 53 - start.year]
    cpx.ew <- cpx.ew[year > min(year) & year < max(year)]
    setkey(cpx.ew, year, week)

    data(cpx_ew_age)

    cpx.ew.age <- cpx.ew.age[week >= start.year, year := as.integer(year + 1)]
    cpx.ew.age <- cpx.ew.age[week >= start.year, week := as.integer(week - 52)]
    cpx.ew.age <- cpx.ew.age[, week := week + 53 - start.year]
    cpx.ew.age <- cpx.ew.age[year > min(year) & year < max(year)]
    setkey(cpx.ew.age, year, week)

    ## 3 week max
    maxes <-
        cpx.ew[, list(max = week[which.max(c(NA, rollapply(rel.incidence,
                      3, sum), NA))]), by = year]
    setnames(maxes, "year", "max.year")

    cpx.ew.age.before.after <-
        cpx.ew.age[, list(before = sum(abs.incidence[week < (maxes[max.year == year,
                          max] - 1)]),
                          after = sum(abs.incidence[week > (maxes[max.year == year,
                          max] + 1)])),
                   by = list(year, lower.age.limit)]
    setkey(cpx.ew.age.before.after, year, lower.age.limit)

    cpx.ew.total.before.after <-
        cpx.ew.age.before.after[, list(before.total = sum(before),
                                       after.total = sum(after)),
                                by = year]
    setkey(cpx.ew.total.before.after, year)

    cpx.ew.age.before.after <-
        merge(cpx.ew.age.before.after, cpx.ew.total.before.after)

    cpx.ew.age.before.after <-
        cpx.ew.age.before.after[, fraction.before := before / before.total]
    cpx.ew.age.before.after <-
        cpx.ew.age.before.after[, fraction.after := after / after.total]

    cpx.ew.age.before.after <-
        cpx.ew.age.before.after[, ratio := fraction.before / fraction.after]

    cpx.ew.age.before.after <-
        cpx.ew.age.before.after[ratio == Inf, ratio := NA]

    p <- ggplot(cpx.ew.age.before.after,
                aes(x = year, y = ratio, color = factor(lower.age.limit)))
    if (lines) {
        p <- p + geom_line(lwd = 1.5)
    } else {
        p <- p + geom_point(size = 2.5)
    }
    p <- p + scale_color_brewer(palette = "Set1", name = "age group",
            labels = limits.to.agegroups(cpx.ew.age.before.after[, lower.age.limit]))
    p <- p + theme_bw(16)
    p <- p + theme(legend.position = "bottom")
    p <- p + geom_hline(aes(yintercept = 1))
    p <- p + facet_grid(lower.age.limit ~ ., scale = "free_y")
    p <- p + ggtitle("Chickenpox")
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("cpx_before_after_3w.", ext, sep = ""), p)
    }

    ## 1 week max
    maxes <- cpx.ew[, list(max = week[which.max(rel.incidence)]), by = year]
    setnames(maxes, "year", "max.year")

    cpx.ew.age.before.after <-
        cpx.ew.age[, list(before = sum(abs.incidence[week < (maxes[max.year == year,
                          max] - 1)]),
                          after = sum(abs.incidence[week > (maxes[max.year == year,
                          max] + 1)])),
                   by = list(year, lower.age.limit)]
    setkey(cpx.ew.age.before.after, year, lower.age.limit)

    cpx.ew.total.before.after <-
        cpx.ew.age.before.after[, list(before.total = sum(before),
                                       after.total = sum(after)),
                                by = year]
    setkey(cpx.ew.total.before.after, year)

    cpx.ew.age.before.after <-
        merge(cpx.ew.age.before.after, cpx.ew.total.before.after)

    cpx.ew.age.before.after <-
        cpx.ew.age.before.after[, fraction.before := before / before.total]
    cpx.ew.age.before.after <-
        cpx.ew.age.before.after[, fraction.after := after / after.total]

    cpx.ew.age.before.after <-
        cpx.ew.age.before.after[, ratio := fraction.before / fraction.after]

    cpx.ew.age.before.after <-
        cpx.ew.age.before.after[ratio == Inf, ratio := NA]

    p <- ggplot(cpx.ew.age.before.after,
                aes(x = year, y = ratio, color = factor(lower.age.limit)))
    if (lines) {
        p <- p + geom_line(lwd = 1.5)
    } else {
        p <- p + geom_point(size = 2.5)
    }
    p <- p + scale_color_brewer(palette = "Set1", name = "age group",
            labels = limits.to.agegroups(cpx.ew.age.before.after[, lower.age.limit]))
    p <- p + theme_bw(16)
    p <- p + theme(legend.position = "bottom")
    p <- p + geom_hline(aes(yintercept = 1))
    p <- p + facet_grid(lower.age.limit ~ ., scale = "free_y")
    p <- p + ggtitle("Chickenpox")
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("cpx_before_after_1w.", ext, sep = ""), p)
    }

    ## single one week max
    maxes <- maxes[, max := round(mean(maxes[, max]))]

    cpx.ew.age.before.after <-
        cpx.ew.age[, list(before = sum(abs.incidence[week < (maxes[max.year == year,
                          max] - 1)]),
                          after = sum(abs.incidence[week > (maxes[max.year == year,
                          max] + 1)])),
                   by = list(year, lower.age.limit)]
    setkey(cpx.ew.age.before.after, year, lower.age.limit)

    cpx.ew.total.before.after <-
        cpx.ew.age.before.after[, list(before.total = sum(before),
                                       after.total = sum(after)),
                                by = year]
    setkey(cpx.ew.total.before.after, year)

    cpx.ew.age.before.after <-
        merge(cpx.ew.age.before.after, cpx.ew.total.before.after)

    cpx.ew.age.before.after <-
        cpx.ew.age.before.after[, fraction.before := before / before.total]
    cpx.ew.age.before.after <-
        cpx.ew.age.before.after[, fraction.after := after / after.total]

    cpx.ew.age.before.after <-
        cpx.ew.age.before.after[, ratio := fraction.before / fraction.after]

    cpx.ew.age.before.after <-
        cpx.ew.age.before.after[ratio == Inf, ratio := NA]

    p <- ggplot(cpx.ew.age.before.after,
                aes(x = year, y = ratio, color = factor(lower.age.limit)))
    if (lines) {
        p <- p + geom_line(lwd = 1.5)
    } else {
        p <- p + geom_point(size = 2.5)
    }
    p <- p + scale_color_brewer(palette = "Set1", name = "age group",
            labels = limits.to.agegroups(cpx.ew.age.before.after[, lower.age.limit]))
    p <- p + theme_bw(16)
    p <- p + theme(legend.position = "bottom")
    p <- p + geom_hline(aes(yintercept = 1))
    p <- p + facet_grid(lower.age.limit ~ ., scale = "free_y")
    p <- p + ggtitle("Chickenpox")
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("cpx_before_after_same.", ext, sep = ""), p)
    }

    ## max at the start of the holiday
    start.holiday <- 30 + 53 - start.year
    maxes <- maxes[, max := start.holiday]

    cpx.ew.age.before.after <-
        cpx.ew.age[, list(before = sum(abs.incidence[week < (maxes[max.year == year,
                          max] - 1)]),
                          after = sum(abs.incidence[week > (maxes[max.year == year,
                          max] + 1)])),
                   by = list(year, lower.age.limit)]
    setkey(cpx.ew.age.before.after, year, lower.age.limit)

    cpx.ew.total.before.after <-
        cpx.ew.age.before.after[, list(before.total = sum(before),
                                       after.total = sum(after)),
                                by = year]
    setkey(cpx.ew.total.before.after, year)

    cpx.ew.age.before.after <-
        merge(cpx.ew.age.before.after, cpx.ew.total.before.after)

    cpx.ew.age.before.after <-
        cpx.ew.age.before.after[, fraction.before := before / before.total]
    cpx.ew.age.before.after <-
        cpx.ew.age.before.after[, fraction.after := after / after.total]

    cpx.ew.age.before.after <-
        cpx.ew.age.before.after[, ratio := fraction.before / fraction.after]

    cpx.ew.age.before.after <-
        cpx.ew.age.before.after[ratio == Inf, ratio := NA]

    p <- ggplot(cpx.ew.age.before.after,
                aes(x = year, y = ratio, color = factor(lower.age.limit)))
    if (lines) {
        p <- p + geom_line(lwd = 1.5)
    } else {
        p <- p + geom_point(size = 2.5)
    }
    p <- p + scale_color_brewer(palette = "Set1", name = "age group",
            labels = limits.to.agegroups(cpx.ew.age.before.after[, lower.age.limit]))
    p <- p + theme_bw(16)
    p <- p + theme(legend.position = "bottom")
    p <- p + geom_hline(aes(yintercept = 1))
    p <- p + facet_grid(lower.age.limit ~ ., scale = "free_y")
    p <- p + ggtitle("Chickenpox")
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("cpx_before_after_holidays.", ext, sep = ""), p)
    }

    ## 3-week max, ignore holiday
    end.holiday <- start.holiday + 5
    maxes <-
        cpx.ew[, list(max = week[which.max(c(NA, rollapply(rel.incidence,
                      3, sum), NA))]), by = year]
    setnames(maxes, "year", "max.year")

    cpx.ew.age.before.after <-
        cpx.ew.age[, list(before = sum(abs.incidence[week < (maxes[max.year == year,
                          max] - 1) & !(week %in% seq(start.holiday, end.holiday))]),
                          after = sum(abs.incidence[week > (maxes[max.year == year,
                          max] + 1) & !(week %in% seq(start.holiday, end.holiday))])),
                   by = list(year, lower.age.limit)]
    setkey(cpx.ew.age.before.after, year, lower.age.limit)

    cpx.ew.total.before.after <-
        cpx.ew.age.before.after[, list(before.total = sum(before),
                                       after.total = sum(after)),
                                by = year]
    setkey(cpx.ew.total.before.after, year)

    cpx.ew.age.before.after <-
        merge(cpx.ew.age.before.after, cpx.ew.total.before.after)

    cpx.ew.age.before.after <-
        cpx.ew.age.before.after[, fraction.before := before / before.total]
    cpx.ew.age.before.after <-
        cpx.ew.age.before.after[, fraction.after := after / after.total]

    cpx.ew.age.before.after <-
        cpx.ew.age.before.after[, ratio := fraction.before / fraction.after]

    cpx.ew.age.before.after <-
        cpx.ew.age.before.after[ratio == Inf, ratio := NA]

    p <- ggplot(cpx.ew.age.before.after,
                aes(x = year, y = ratio, color = factor(lower.age.limit)))
    if (lines) {
        p <- p + geom_line(lwd = 1.5)
    } else {
        p <- p + geom_point(size = 2.5)
    }
    p <- p + scale_color_brewer(palette = "Set1", name = "age group",
            labels = limits.to.agegroups(cpx.ew.age.before.after[, lower.age.limit]))
    p <- p + theme_bw(16)
    p <- p + theme(legend.position = "bottom")
    p <- p + geom_hline(aes(yintercept = 1))
    p <- p + facet_grid(lower.age.limit ~ ., scale = "free_y")
    p <- p + ggtitle("Chickenpox")
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("cpx_before_after_3w_nohol.", ext, sep = ""), p)
    }

}

##' Plot England & Wales demographic data
##'
##' @param ext plot extensions
##' @author Sebastian Funk
##' @import ggplot2
plotDemographicEWData <- function(ext = NA) {

    data(pop_ew_age)
    data(cpx_ew_age)

    pop.ew.age <-
        reduce.agegroups(pop.ew.age, unique(cpx.ew.age[, lower.age.limit]))

    p <- ggplot(pop.ew.age, aes(x = year, fill = factor(lower.age.limit),
                                y = population))
    p <- p + scale_fill_brewer("Age group", palette = "Set1",
                               labels = limits.to.agegroups(unique(cpx.ew.age[,
                               lower.age.limit])))
    p <- p + geom_bar(stat = "identity", position = "stack")
    p <- p + theme_bw(16)
    p <- p + theme(legend.position = "bottom")
    p <- p + ggtitle("Demographics")
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("ew_population.", ext, sep = ""), p)
    }

    data(births_ew)

    p <- ggplot(births.ew, aes(x = year, y = total))
    p <- p + geom_line(lwd = 1.5)
    p <- p + theme_bw(16)
    p <- p + coord_cartesian(ylim = c(0, max(births.ew[, total]) * 1.1) )
    p <- p + ggtitle("Births")
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("ew_births.", ext, sep = ""), p)
    }

    p <- ggplot(births.ew, aes(x = year, y = cbr))
    p <- p + geom_line(lwd = 1.5)
    p <- p + theme_bw(16)
    p <- p + coord_cartesian(ylim = c(0, max(births.ew[, cbr]) * 1.1) )
    p <- p + ggtitle("Births")
    p <- p + scale_y_continuous("Crude birth rate")
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("ew_cbr.", ext, sep = ""), p)
    }

}

##' Plot age distributions for measles and chickenpox according to
##' different mixing assumptions
##'
##' @param ext plot extensions
##' @param classic.life.expectancy life expectancy in the "classic" model
##' @param kids.mult mixing multiplier for children
##' @param max.age maximum age
##' @author Sebastian Funk
##' @import ggplot2 reshape2 Hmisc
##' @export
plot_mixing_vacc_ages <- function(ext = NA, classic.life.expectancy = 80,
                                   kids.mult = 5, max.age = 99) {
    data(pop_ew_age)
    data(cpx_ew_age)
    data(ms_ew_age)
    data(vaccine_ew)
    data(births_ew)

    # set crude age limits
    lower.limits <- c(0, 5, 15)

    ## first, chickenpox
    contact.matrix <- sample.uk.polymod(lower.limits)
    hom.matrix <- homogeneous.mixing.matrix(contact.matrix)
    diagonal.matrix <- diag(diag(contact.matrix))
    colnames(diagonal.matrix) <- colnames(contact.matrix)
    rownames(diagonal.matrix) <- rownames(contact.matrix)

    R0 <- c(1.01, 1.1, seq(1.5, 19.5, by = 0.5))

    social.age.dist <- age.distribution(contact.matrix, R0 = R0)
    hom.age.dist <- age.distribution(hom.matrix, R0 = R0)
    diagonal.age.dist <- age.distribution(diagonal.matrix, R0 = R0)

    social.age.dist <-
        social.age.dist[, upper.age.limit :=
                            lower.to.upper.limits(lower.age.limit, lower.limits,
                                                  max.age)]
    hom.age.dist <-
        hom.age.dist[, upper.age.limit :=
                         lower.to.upper.limits(lower.age.limit, lower.limits, max.age)]

    diagonal.age.dist <-
        diagonal.age.dist[, upper.age.limit :=
                         lower.to.upper.limits(lower.age.limit, lower.limits, max.age)]

    cpx.ew.age <- reduce.agegroups(cpx.ew.age, lower.limits)
    cpx.ew.age <-
        cpx.ew.age[, upper.age.limit :=
                       lower.to.upper.limits(lower.age.limit, lower.limits,
                                             max.age)]
    cpx.ew.age <- cpx.ew.age[, list(abs.incidence = sum(abs.incidence),
                                    upper.age.limit = mean(upper.age.limit),
                                    population = sum(population)),
                             by = list(year, week, lower.age.limit)]
    cpx.ew.age.annual <- cpx.ew.age[, list(abs.incidence = sum(abs.incidence),
                                           upper.age.limit = mean(upper.age.limit),
                                           population = mean(population)),
                                    by = list(year, lower.age.limit)]
    cpx.ew.age.all <-
        cpx.ew.age.annual[, list(abs.incidence = sum(abs.incidence),
                                 upper.age.limit = mean(upper.age.limit)),
                          by = list(lower.age.limit)]

    ms.ew.age <- reduce.agegroups(ms.ew.age, lower.limits)
    ms.ew.age <-
        ms.ew.age[, upper.age.limit :=
                      lower.to.upper.limits(lower.age.limit, lower.limits,
                                             max.age)]
    ms.ew.age.annual <- ms.ew.age[, list(abs.incidence = sum(abs.incidence),
                                  upper.age.limit = mean(upper.age.limit),
                                  population = sum(population)),
                                  by = list(year, lower.age.limit)]
    ms.ew.age.all <-
        ms.ew.age.annual[, list(abs.incidence = sum(abs.incidence),
                                upper.age.limit = mean(upper.age.limit)),
                         by = list(lower.age.limit)]

    
    ## mean age at infection
    meanage.hom <-
        hom.age.dist[, list(low.mean.age = sum(proportion.cases * lower.age.limit),
                            high.mean.age = sum(proportion.cases * upper.age.limit)),
                     by = list(R0)]
    meanage.classic <-
        data.table(R0 = meanage.hom[, R0],
                   low.mean.age = classic.life.expectancy / (meanage.hom[, R0]),
                   high.mean.age = classic.life.expectancy / (meanage.hom[, R0]))
    meanage.social <-
        social.age.dist[, list(low.mean.age = sum(proportion.cases *
                                                         lower.age.limit),
                                  high.mean.age = sum(proportion.cases *
                                                          upper.age.limit)),
                        by = list(R0)]
    meanage.diagonal <- 
        diagonal.age.dist[, list(low.mean.age = sum(proportion.cases *
                                                         lower.age.limit),
                                  high.mean.age = sum(proportion.cases *
                                                          upper.age.limit)),
                        by = list(R0)]
    

    meanage.classic <- meanage.classic[, model := "classic"]
    meanage.hom <- meanage.hom[, model := "homogeneous"]
    meanage.social <- meanage.social[, model := "social"]
    meanage.diagonal <- meanage.diagonal[, model := "diagonal"]
    meanage <- rbind(meanage.classic, meanage.hom, meanage.social, meanage.diagonal)
    meanage <- meanage[, mean.mean.age := (low.mean.age + high.mean.age) / 2]

    meanage.data.cpx <-
        (cpx.ew.age.all[, wtd.mean(lower.age.limit, abs.incidence)] +
             cpx.ew.age.all[, wtd.mean(upper.age.limit, abs.incidence)]) / 2
    meanage.data.ms <-
        (ms.ew.age.all[, wtd.mean(lower.age.limit, abs.incidence)] +
             ms.ew.age.all[, wtd.mean(upper.age.limit, abs.incidence)]) / 2

    p <- ggplot(meanage, aes(x = R0, ymin = low.mean.age, y = mean.mean.age, 
                             ymax = high.mean.age, color = model, fill = model))
    p <- p + geom_line(lwd = 1.2)
    p <- p + geom_hline(yintercept = meanage.data.cpx, lwd = 1.2, color = "black")
    p <- p + geom_hline(yintercept = meanage.data.ms, lwd = 1.2, color = "black",
                        linetype = 2)
    p <- p + theme_bw(20)
    p <- p + theme(legend.position = "bottom")
    p <- p + scale_color_brewer(palette = "Set1")
    p <- p + scale_y_continuous("mean age at infection estimate")
    p <- p + scale_x_continuous(expression(R[0]))
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("meanage_eq_model", ext, sep = "."), p, width = 9)
    }
    p <- ggplot(meanage, aes(x = R0, ymin = low.mean.age, y = mean.mean.age, 
                             ymax = high.mean.age, color = model, fill = model))
    p <- p + geom_ribbon(alpha = 0.6, lwd = 1.2)
    p <- p + geom_hline(yintercept = meanage.data.cpx, lwd = 1.2, color = "black")
    p <- p + geom_hline(yintercept = meanage.data.ms, lwd = 1.2, color = "black",
                        linetype = 2)
    p <- p + theme_bw(20)
    p <- p + theme(legend.position = "bottom")
    p <- p + scale_color_brewer(palette = "Set1")
    p <- p + scale_fill_brewer(palette = "Set1")
    p <- p + scale_y_continuous("mean age at infection estimate")
    p <- p + scale_x_continuous(expression(R[0]))
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("meanage_eq_model_ribbons", ext, sep = "."), p, width = 9)
    }

    hom.age.dist <- hom.age.dist[, model := "homogeneous"]
    social.age.dist <- social.age.dist[, model := "social"]
    diagonal.age.dist <- diagonal.age.dist[, model := "diagonal"]

    age.dist <- rbind(hom.age.dist, social.age.dist, diagonal.age.dist)
    age.labels <- limits.to.agegroups(cpx.ew.age[, lower.age.limit])

    p <- ggplot(age.dist, aes(x = R0, y = proportion.cases,
                                 color = factor(lower.age.limit)))
    p <- p + geom_line(lwd = 1.2)
    p <- p + theme_bw(20)
    p <- p + theme(legend.position = "bottom")
    p <- p + scale_color_brewer("age group", palette = "Set1",
                                labels = age.labels)
    p <- p + facet_grid(model ~ ., scales = "free")
    p <- p + scale_y_continuous("proportion of cases", limits = c(0, 1))
    p <- p + scale_x_continuous(expression(R[0]))
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("prop_cases_eq_model", ext, sep = "."), p, width = 9)
    }

    cpx.ew.age.total <-
        cpx.ew.age.annual[, list(abs.pop.incidence = sum(abs.incidence),
                                abs.population = sum(population)),
                         by = list(year)]
    setkey(cpx.ew.age.annual, year)
    setkey(cpx.ew.age.total, year)
    cpx.ew.age.annual <- merge(cpx.ew.age.annual, cpx.ew.age.total, by = "year")
    cpx.ew.age.annual <-
        cpx.ew.age.annual[, proportion.cases := abs.incidence / abs.pop.incidence]
    cpx.ew.age.annual[, infection := "chickenpox"]

    ms.ew.age.total <-
        ms.ew.age.annual[, list(abs.pop.incidence = sum(abs.incidence),
                                abs.population = sum(population)),
                         by = list(year)]
    setkey(ms.ew.age.annual, year)
    setkey(ms.ew.age.total, year)
    ms.ew.age.annual <- merge(ms.ew.age.annual, ms.ew.age.total, by = "year")
    ms.ew.age.annual <-
        ms.ew.age.annual[, proportion.cases := abs.incidence / abs.pop.incidence]
    ms.ew.age.annual[, infection := "measles"]

    ch.ew.age.annual <- rbind(cpx.ew.age.annual, ms.ew.age.annual)

    p <- ggplot(ch.ew.age.annual, aes(x = year, y = proportion.cases,
                                      color = factor(lower.age.limit)))
    p <- p + geom_line(lwd = 1.2)
    p <- p + theme_bw(20)
    p <- p + theme(legend.position = "bottom")
    p <- p + scale_color_brewer("age group", palette = "Set1",
                                labels = age.labels)
    p <- p + facet_grid(infection ~ ., scales = "free")
    p <- p + scale_y_continuous("proportion of cases", limits = c(0, 1))
    p <- p + scale_x_continuous("")
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("prop_cases_data", ext, sep = "."), p, width = 9)
    }

    ## now, vaccination

    vacc <- seq(from = 0, to = 0.95, by = 0.05)

    ## 1) measles R0
    R0.ms.diagonal <- approx(meanage[model == "diagonal", mean.mean.age],
                           meanage[model == "diagonal", R0], meanage.data.ms)$y
    R0.ms.social <- approx(meanage[model == "social", mean.mean.age],
                           meanage[model == "social", R0], meanage.data.ms)$y
    R0.ms.hom <- approx(meanage[model == "homogeneous", mean.mean.age],
                        meanage[model == "homogeneous", R0], meanage.data.ms)$y
    R0.ms.classic <- classic.life.expectancy / meanage.data.ms

    ## 2) chickenpox R0
    R0.cpx.diagonal <- approx(meanage[model == "diagonal", mean.mean.age],
                           meanage[model == "diagonal", R0], meanage.data.cpx)$y
    R0.cpx.social <- approx(meanage[model == "social", mean.mean.age],
                           meanage[model == "social", R0], meanage.data.cpx)$y
    R0.cpx.hom <- approx(meanage[model == "homogeneous", mean.mean.age],
                        meanage[model == "homogeneous", R0], meanage.data.cpx)$y
    R0.cpx.classic <- classic.life.expectancy / meanage.data.cpx

    hom.ms.vacc <- rbindlist(lapply(vacc, function(x)
    {
        age.distribution(hom.matrix, R0 = R0.ms.hom,
                         vaccination = x)
    }))
    social.ms.vacc <- rbindlist(lapply(vacc, function(x)
    {
        age.distribution(contact.matrix, R0 = R0.ms.social,
                         vaccination = x)
    }))
    diagonal.ms.vacc <- rbindlist(lapply(vacc, function(x)
    {
        age.distribution(diagonal.matrix, R0 = R0.ms.diagonal,
                         vaccination = x)
    }))
    hom.cpx.vacc <- rbindlist(lapply(vacc, function(x)
    {
        age.distribution(hom.matrix, R0 = R0.cpx.hom,
                         vaccination = x)
    }))
    social.cpx.vacc <- rbindlist(lapply(vacc, function(x)
    {
        age.distribution(contact.matrix, R0 = R0.cpx.social,
                         vaccination = x)
    }))
    diagonal.cpx.vacc <- rbindlist(lapply(vacc, function(x)
    {
        age.distribution(diagonal.matrix, R0 = R0.cpx.diagonal,
                         vaccination = x)
    }))

    hom.ms.vacc <- hom.ms.vacc[, mixing := "homogeneous"]
    social.ms.vacc <- social.ms.vacc[, mixing := "social"]
    diagonal.ms.vacc <- diagonal.ms.vacc[, mixing := "diagonal"]

    hom.cpx.vacc <- hom.cpx.vacc[, mixing := "homogeneous"]
    social.cpx.vacc <- social.cpx.vacc[, mixing := "social"]
    diagonal.cpx.vacc <- diagonal.cpx.vacc[, mixing := "diagonal"]

    ms.vacc <- rbind(hom.ms.vacc, social.ms.vacc)
    cpx.vacc <- rbind(hom.cpx.vacc, social.cpx.vacc)

    ms.vacc[, infection := "measles"]
    cpx.vacc[, infection := "chickenpox"]

    vacc <- rbind(ms.vacc, cpx.vacc)

    p <- ggplot(vacc, aes(x = vaccination, y = proportion.cases,
                          color = factor(lower.age.limit)))
    p <- p + geom_line(lwd = 1.2)
    p <- p + theme_bw(20)
    p <- p + theme(legend.position = "bottom")
    p <- p + scale_color_brewer("age group", palette = "Set1",
                                labels = age.labels)
    p <- p + scale_y_continuous("proportion of cases", limits = c(0, 1))
    p <- p + scale_x_continuous("vaccination coverage")
    p <- p + facet_grid(infection ~ mixing)

    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("prop_cases_vacc", ext, sep = "."), p, width = 9)
    }

    dynamic.uptake <- c(0.1, 0.5, 0.8, 0.9, 0.95)
    
    fixed.parameters <- list(mixing = contact.matrix,
                             interpolate = F,
                             runin.time = 100,
                             mu = births.ew[year == 2010, total] /
                                 sum(pop.ew.age[year == 2010, population]),
                             maternal.immunity = 0,
                             ageing = 1/diff(lower.limits))

    l <- lapply(dynamic.uptake, function(x)
    {
        fixed.parameters[["vaccine.uptake"]] <- matrix(c(rep(c(0, 0, 0), 10), c(x, 0, 0)),
                                                       ncol = length(lower.limits), byrow = T)

        # 1) Measles
        # first: social mixing
        theta.ms <- sample.ms.prior(fixed.parameters)
        theta.ms[["R0"]] <- R0.ms.social

        parameters.ms <- c(theta.ms, fixed.parameters)
    
        init <- ms.init(parameters.ms)
        ms.social.traj <- runSIR(parameters = parameters.ms,
                                 init = init,
                                 times = seq_len(50),
                                 age.labels = lower.limits,
                                 thickening = 10)
        ms.social.traj[, mixing := "social"]

        # second: homogeneous mixing
        theta.ms[["R0"]] <- R0.ms.hom

        parameters.ms <- c(theta.ms, fixed.parameters)
    
        init <- ms.init(parameters.ms)
        ms.hom.traj <- runSIR(parameters = parameters.ms,
                                 init = init,
                                 times = seq_len(50),
                                 age.labels = lower.limits,
                                 thickening = 10)
        ms.hom.traj[, mixing := "homogeneous"]

        # third: diagonal mixing
        theta.ms[["R0"]] <- R0.ms.diag

        parameters.ms <- c(theta.ms, fixed.parameters)
    
        init <- ms.init(parameters.ms)
        ms.diag.traj <- runSIR(parameters = parameters.ms,
                                 init = init,
                                 times = seq_len(50),
                                 age.labels = lower.limits,
                                 thickening = 10)
        ms.diag.traj[, mixing := "diagonal"]


        # 1) Chickenpox
        # first: social mixing
        theta.cpx <- sample.cpx.prior(fixed.parameters)

        theta.cpx[["R0"]] <- R0.cpx.social

        parameters.cpx <- c(theta.cpx, fixed.parameters)
    
        init <- cpx.init(parameters.cpx)
        cpx.social.traj <- runSIR(parameters = parameters.cpx,
                                 init = init,
                                 times = seq_len(50),
                                 age.labels = lower.limits,
                                 thickening = 10)
        cpx.social.traj[, mixing := "social"]

        # second: homogeneous mixing
        theta.cpx[["R0"]] <- R0.cpx.hom

        parameters.cpx <- c(theta.cpx, fixed.parameters)
    
        init <- cpx.init(parameters.cpx)
        cpx.hom.traj <- runSIR(parameters = parameters.cpx,
                                 init = init,
                                 times = seq_len(50),
                                 age.labels = lower.limits,
                                 thickening = 10)
        cpx.hom.traj[, mixing := "homogeneous"]

        ms.traj <- rbind(ms.social.traj, ms.hom.traj)
        ms.traj[, infection := "measles"]

        cpx.traj <- rbind(cpx.social.traj, cpx.hom.traj)
        cpx.traj[, infection := "chickenpox"]
        
        dynamic.vacc <- rbind(ms.traj, cpx.traj)
        dynamic.vacc <- dynamic.vacc[, uptake := x]
    })

    dynamic.vacc.all <- rbindlist(l)
    dynamic.vacc.all <- dynamic.vacc.all[!is.na(abs.incidence)]
    
    dynamic.vacc.all.total <-
        dynamic.vacc.all[, list(abs.pop.incidence = sum(abs.incidence, na.rm = T)),
                         by = list(time, mixing, infection, uptake)]
    dynamic.vacc.all.total <- dynamic.vacc.all.total[abs.pop.incidence > 0]
                                
    setkey(dynamic.vacc.all, time, mixing, infection, uptake)
    setkey(dynamic.vacc.all.total, time, mixing, infection, uptake)

    dynamic.vacc.all <- merge(dynamic.vacc.all, dynamic.vacc.all.total, by = c("time", "mixing", "infection", "uptake"))
    dynamic.vacc.all <-
        dynamic.vacc.all[, proportion.cases := abs.incidence / abs.pop.incidence]
    
    p <- ggplot(dynamic.vacc.all[infection == "measles"],
                aes(x = time, y = proportion.cases,
                    color = factor(lower.age.limit)))
    p <- p + geom_line(lwd = 1.2)
    p <- p + theme_bw(16)
    p <- p + theme(legend.position = "bottom")
    p <- p + scale_color_brewer("age group", palette = "Set1",
                                labels = age.labels)
    p <- p + scale_y_continuous("proportion of cases", limits = c(0, 1))
    p <- p + scale_x_continuous("vaccination coverage")
    p <- p + facet_grid(uptake ~ mixing)

    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("prop_cases_ms_vacc_dynamic", ext, sep = "."), p, width = 9)
    }

    p <- ggplot(dynamic.vacc.all[infection == "chickenpox"],
                aes(x = time, y = proportion.cases,
                    color = factor(lower.age.limit)))
    p <- p + geom_line(lwd = 1.2)
    p <- p + theme_bw(16)
    p <- p + theme(legend.position = "bottom")
    p <- p + scale_color_brewer("age group", palette = "Set1",
                                labels = age.labels)
    p <- p + scale_y_continuous("proportion of cases", limits = c(0, 1))
    p <- p + scale_x_continuous("years")
    p <- p + facet_grid(uptake ~ mixing)

    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("prop_cases_cpx_vacc_dynamic", ext, sep = "."), p, width = 9)
    }

    ## run dynamic model from equilibrium, with vaccination at a few
    ## different levels introduced at birth (e.g., 50, 80, 90, 95)
    ## observe change in age distribution
    
    ## compare with measles

    years <- seq(1944, 2010)

    ## polymod: approximate R0 from mean age of infection
    mmr <- estimate.immunity.mmr(years = years, age.limits = lower.limits)

    social.age.dist.ms.vacc <-
        age.distribution(contact.matrix, R0 = R0.ms.social,
                         vaccination = mmr, years = years)
    hom.age.dist.ms.vacc <- age.distribution(hom.matrix, R0 = R0.ms.hom,
                                             vaccination = mmr, years = years)

    ## for the classic age of infection, we test 1) the overall
    ## vaccination level and 2) vaccination levels at birth
    min.diff <- min(pop.ew.age[, year]) - min(vaccine.ew[, year])
    if (min.diff > 0) {
        previous.pop <-
            backcalc.population(min.diff, unique(ms.ew.age[, lower.age.limit]))
        previous.pop <-
            previous.pop[, year := vaccine.ew[year < min(pop.ew.age[, year]), year]]
        m.previous.pop <- melt(previous.pop, id.vars = "year",
                               variable.name = "lower.age.limit")
        setnames(m.previous.pop, "value", "population")
        m.previous.pop <-
            m.previous.pop[, lower.age.limit :=
                               as.numeric(as.character(lower.age.limit))]
        pop.ew.age <- rbind(m.previous.pop, pop.ew.age)
    }

    pop.ew.age <- reduce.agegroups(pop.ew.age, unique(ms.ew.age[, lower.age.limit]))
    pop.ew.age <- pop.ew.age[, list(population = sum(population)),
                             by = list(year, lower.age.limit)]
    setkey(pop.ew.age, year, lower.age.limit)

    ## calculate overall vaccination levels
    m.mmr <- melt(mmr, id.vars = "year", variable.name = "lower.age.limit")
    setnames(m.mmr, "value", "rel.vaccinated")
    m.mmr <- m.mmr[, lower.age.limit := as.numeric(as.character(lower.age.limit))]
    setkey(m.mmr, year, lower.age.limit)

    m.mmr <- merge(m.mmr, pop.ew.age, all.x = T)

    m.mmr <- m.mmr[, list(rel.vaccinated = wtd.mean(rel.vaccinated, population)),
                   by = list(year)]
    m.mmr <- m.mmr[!is.finite(rel.vaccinated), rel.vaccinated := 0]

    ## mean age at infection
    meanage.hom.ms.vacc <-
        hom.age.dist.ms.vacc[,
                             list(mean.age = sum(proportion.cases * lower.age.limit)),
                             by = list(year)]
    meanage.classic.ms.vacc.all <-
        data.table(year = years,
                   mean.age = classic.life.expectancy /
                       (R0.ms.classic * (1 - m.mmr[, rel.vaccinated])))
    ## meanage.classic.ms.vacc.newborn <-
    ##     data.table(year = meanage.hom.ms.vacc[, year],
    ##                mean.age = classic.life.expectancy / (R0.ms.classic *
    ##                    (1 - c(rep(0, min(vaccine.ew[, year]) - min(years)),
    ##                           vaccine.ew[year %in% years, MCV1 / 100]))))
    meanage.social.ms.vacc <-
        social.age.dist.ms.vacc[, list(mean.age =
                                           sum(proportion.cases * lower.age.limit)),
                                by = list(year)]
    meanage.data.ms.vacc <-
        data.table(year = years,
                   mean.age = sapply(years, function(x) {
                       ms.ew.age[year == x, wtd.mean(lower.age.limit, abs.incidence)]
                   }))

    meanage.classic.ms.vacc.all <-
        meanage.classic.ms.vacc.all[, model := "classic"]
    ## meanage.classic.ms.vacc.newborn <-
    ##     meanage.classic.ms.vacc.newborn[, model := "classic (birth vaccinations)"]
    meanage.hom.ms.vacc <- meanage.hom.ms.vacc[, model := "homogenous"]
    meanage.social.ms.vacc <- meanage.social.ms.vacc[, model := "social"]
    meanage.data.ms.vacc <- meanage.data.ms.vacc[, model := "data"]

    meanage.ms.vacc <-
        rbind(meanage.classic.ms.vacc.all, meanage.hom.ms.vacc,
              meanage.social.ms.vacc, meanage.data.ms.vacc)

    p <- ggplot(meanage.ms.vacc, aes(x = year, y = mean.age, color = model))
    p <- p + geom_line(lwd = 1.2)
    p <- p + theme_bw(20)
    p <- p + theme(legend.position = "bottom")
    p <- p + scale_color_brewer(palette = "Set1")
    p <- p + scale_y_continuous("mean age at infection estimate")
    p <- p + scale_x_continuous(expression(R[0]))
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("meanage_ms_vacc_eq_model", ext, sep = "."), p)
    }

    hom.age.dist.ms.vacc <- hom.age.dist.ms.vacc[, model := "homogeneous"]
    social.age.dist.ms.vacc <- social.age.dist.ms.vacc[, model := "social"]

    ## get proportion of cases in the real world
    ms.ew.age.total <- ms.ew.age[, list(abs.pop.incidence = sum(abs.incidence),
                                        abs.population = sum(population)),
                                 by = list(year)]
    setkey(ms.ew.age.total, year)
    ms.ew.age <- merge(ms.ew.age, ms.ew.age.total, by = "year")
    ms.ew.age <- ms.ew.age[, proportion.cases := abs.incidence / abs.pop.incidence]
    ms.ew.age <- ms.ew.age[, abs.pop.incidence := NULL]
    ms.ew.age <- ms.ew.age[, abs.population := NULL]
    ms.ew.age <- ms.ew.age[, model := "data"]

    

    age.dist.ms.vacc <- rbind(hom.age.dist.ms.vacc, social.age.dist.ms.vacc,
                              ms.ew.age[year < max(years)], use.names = T,
                              fill = TRUE)

    p <- ggplot(age.dist.ms.vacc, aes(x = year, y = proportion.cases,
                                 color = factor(lower.age.limit)))
    p <- p + geom_line(lwd = 1.2)
    p <- p + theme_bw(20)
    p <- p + theme(legend.position = "bottom")
    p <- p +
        scale_color_brewer("age group", palette = "Dark2",
                           labels = limits.to.agegroups(ms.ew.age[, lower.age.limit]))
    p <- p + facet_grid(model ~ .)
    p <- p + scale_y_continuous("proportion of cases")
    p <- p + scale_x_continuous("year")
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("prop_cases_ms_vacc_eq_model", ext, sep = "."), p)
    }

    R0.ms <- c(1.01, 1.1, seq(1.5, 39.5, by = 0.5))

    ## strengthen diagonal
    contact.matrix.diag <- diag(diag(contact.matrix))
    colnames(contact.matrix.diag) <- colnames(contact.matrix)
    rownames(contact.matrix.diag) <- rownames(contact.matrix)

    hom.matrix.diag <- diag(diag(hom.matrix))
    colnames(hom.matrix.diag) <- colnames(hom.matrix)
    rownames(hom.matrix.diag) <- rownames(hom.matrix)

    social.age.dist.ms.diag <- age.distribution(contact.matrix.diag, R0 = R0.ms)
    hom.age.dist.ms.diag <- age.distribution(hom.matrix.diag, R0 = R0.ms)

    ## mean age at infection
    meanage.hom.ms.diag <-
        hom.age.dist.ms.diag[, list(mean.age = sum(proportion.cases *
                                                       lower.age.limit)),
                        by = list(R0)]
    meanage.classic.ms.diag <-
        data.table(R0 = meanage.hom.ms.diag[, R0],
                   mean.age = classic.life.expectancy / (meanage.hom.ms.diag[, R0]))
    meanage.social.ms.diag <-
        social.age.dist.ms.diag[, list(mean.age =
                                           sum(proportion.cases * lower.age.limit)),
                           by = list(R0)]

    meanage.classic.ms.diag <- meanage.classic.ms.diag[, model := "classic"]
    meanage.hom.ms.diag <- meanage.hom.ms.diag[, model := "homogeneous"]
    meanage.social.ms.diag <- meanage.social.ms.diag[, model := "social"]
    meanage.ms.diag <-
        rbind(meanage.classic.ms.diag, meanage.hom.ms.diag, meanage.social.ms.diag)

    p <- ggplot(meanage.ms.diag, aes(x = R0, y = mean.age, color = model))
    p <- p + geom_line(lwd = 1.2)
    p <- p + geom_hline(yintercept = meanage.data.ms, lwd = 1.2, color = "black")
    p <- p + theme_bw(20)
    p <- p + theme(legend.position = "bottom")
    p <- p + scale_color_brewer(palette = "Set1")
    p <- p + scale_y_continuous("mean age at infection estimate")
    p <- p + scale_x_continuous(expression(R[0]))
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("meanage_ms_eq_diag_model", ext, sep = "."), p)
    }

    hom.age.dist.ms.diag <- hom.age.dist.ms.diag[, model := "homogeneous"]
    social.age.dist.ms.diag <- social.age.dist.ms.diag[, model := "social"]

    age.dist.ms.diag <- rbind(hom.age.dist.ms.diag, social.age.dist.ms.diag)

    p <- ggplot(age.dist.ms.diag, aes(x = R0, y = proportion.cases,
                                 color = factor(lower.age.limit)))
    p <- p + geom_line(lwd = 1.2)
    p <- p + theme_bw(20)
    p <- p + theme(legend.position = "bottom")
    p <- p +
        scale_color_brewer("age group", palette = "Dark2",
                           labels = limits.to.agegroups(ms.ew.age[, lower.age.limit]))
    p <- p + facet_grid(model ~ ., scales = "free")
    p <- p + scale_y_continuous("proportion of cases")
    p <- p + scale_x_continuous(expression(R[0]))
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("prop_cases_ms_diag_eq_model", ext, sep = "."), p)
    }

    ## now, vaccination

    years <- seq(1944, 2010)

    ## polymod: approximate R0 from mean age of infection
    R0.ms.diag.social <- approx(meanage.social.ms.diag[, mean.age],
                                meanage.social.ms.diag[, R0], meanage.data.ms)$y
    R0.ms.diag.hom <- approx(meanage.hom.ms.diag[, mean.age],
                             meanage.hom.ms.diag[, R0], meanage.data.ms)$y
    R0.ms.diag.classic <- classic.life.expectancy / meanage.data.ms

    mmr <- estimate.immunity.mmr(years = years, age.limits = lower.limits)

    social.age.dist.ms.diag.vacc <-
        age.distribution(contact.matrix.diag, R0 = R0.ms.diag.social,
                         vaccination = mmr, years = years)
    hom.age.dist.ms.diag.vacc <-
        age.distribution(hom.matrix.diag, R0 = R0.ms.diag.hom,
                         vaccination = mmr, years = years)

    social.age.dist.ms.diag.vacc <-
        age.distribution(contact.matrix.diag, R0 = R0.ms.diag.social,
                         vaccination = mmr, years = years)
    hom.age.dist.ms.diag.vacc <- age.distribution(hom.matrix.diag, R0 = R0.ms.diag.hom,
                                             vaccination = mmr, years = years)

    ## mean age at infection
    meanage.hom.ms.diag.vacc <-
        hom.age.dist.ms.diag.vacc[,
                             list(mean.age = sum(proportion.cases * lower.age.limit)),
                             by = list(year)]
    meanage.classic.ms.diag.vacc.all <-
        data.table(year = years,
                   mean.age = classic.life.expectancy /
                       (R0.ms.diag.classic * (1 - m.mmr[, rel.vaccinated])))
    ## meanage.classic.ms.diag.vacc.newborn <-
    ##     data.table(year = meanage.hom.ms.diag.vacc[, year],
    ##                mean.age = classic.life.expectancy / (R0.ms.diag.classic *
    ##                    (1 - c(rep(0, min(vaccine.ew[, year]) - min(years)),
    ##                           vaccine.ew[year %in% years, MCV1 / 100]))))
    meanage.social.ms.diag.vacc <-
        social.age.dist.ms.diag.vacc[, list(mean.age = sum(proportion.cases *
                                                               lower.age.limit)),
                                     by = list(year)]
    meanage.data.ms.diag.vacc <-
        data.table(year = years,
                   mean.age = sapply(years, function(x) {
                       ms.ew.age[year == x, wtd.mean(lower.age.limit, abs.incidence)]
                   }))

    meanage.classic.ms.diag.vacc.all <-
        meanage.classic.ms.diag.vacc.all[, model := "classic"]
    ## meanage.classic.ms.diag.vacc.newborn <-
    ##     meanage.classic.ms.diag.vacc.newborn[, model := "classic (birth vaccinations)"]
    meanage.hom.ms.diag.vacc <- meanage.hom.ms.diag.vacc[, model := "homogenous"]
    meanage.social.ms.diag.vacc <- meanage.social.ms.diag.vacc[, model := "social"]
    meanage.data.ms.diag.vacc <- meanage.data.ms.diag.vacc[, model := "data"]

    meanage.ms.diag.vacc <-
        rbind(meanage.classic.ms.diag.vacc.all, meanage.hom.ms.diag.vacc,
              meanage.social.ms.diag.vacc, meanage.data.ms.diag.vacc)

    p <- ggplot(meanage.ms.diag.vacc, aes(x = year, y = mean.age, color = model))
    p <- p + geom_line(lwd = 1.2)
    p <- p + theme_bw(20)
    p <- p + theme(legend.position = "bottom")
    p <- p + scale_color_brewer(palette = "Set1")
    p <- p + scale_y_continuous("mean age at infection estimate")
    p <- p + scale_x_continuous(expression(R[0]))
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("meanage_ms_diag_vacc_eq_model", ext, sep = "."), p)
    }

    hom.age.dist.ms.diag.vacc <- hom.age.dist.ms.diag.vacc[, model := "homogeneous"]
    social.age.dist.ms.diag.vacc <- social.age.dist.ms.diag.vacc[, model := "social"]

    ## get proportion of cases in the real world
    ms.ew.age.total <- ms.ew.age[, list(abs.pop.incidence = sum(abs.incidence),
                                        abs.population = sum(population)),
                                 by = list(year)]
    setkey(ms.ew.age.total, year)
    ms.ew.age <- merge(ms.ew.age, ms.ew.age.total, by = "year")
    ms.ew.age <- ms.ew.age[, proportion.cases := abs.incidence / abs.pop.incidence]
    ms.ew.age <- ms.ew.age[, abs.pop.incidence := NULL]
    ms.ew.age <- ms.ew.age[, abs.population := NULL]
    ms.ew.age <- ms.ew.age[, model := "data"]

    age.dist.ms.diag.vacc <-
        rbind(hom.age.dist.ms.diag.vacc, social.age.dist.ms.diag.vacc,
              ms.ew.age[year < max(years)], use.names = T, fill = TRUE)

    p <- ggplot(age.dist.ms.diag.vacc, aes(x = year, y = proportion.cases,
                                 color = factor(lower.age.limit)))
    p <- p + geom_line(lwd = 1.2)
    p <- p + theme_bw(20)
    p <- p + theme(legend.position = "bottom")
    p <-
        p + scale_color_brewer("age group", palette = "Dark2",
                               labels = limits.to.agegroups(ms.ew.age[,
                               lower.age.limit]))
    p <- p + facet_grid(model ~ .)
    p <- p + scale_y_continuous("proportion of cases")
    p <- p + scale_x_continuous("year")
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("prop_cases_ms_diag_vacc_eq_model", ext, sep = "."), p)
    }

    contact.matrix.kids <- copy(contact.matrix)

    contact.matrix.kids[1,1] <- contact.matrix.kids[1, 1] * kids.mult
    hom.matrix.kids <- copy(hom.matrix)

    social.age.dist.ms.kids <- age.distribution(contact.matrix.kids, R0 = R0.ms)
    hom.age.dist.ms.kids <- age.distribution(hom.matrix.kids, R0 = R0.ms)

    ## mean age at infection
    meanage.hom.ms.kids <-
        hom.age.dist.ms.kids[, list(mean.age = sum(proportion.cases *
                                                       lower.age.limit)),
                        by = list(R0)]
    meanage.classic.ms.kids <-
        data.table(R0 = meanage.hom.ms.kids[, R0],
                   mean.age = classic.life.expectancy / (meanage.hom.ms.kids[, R0]))
    meanage.social.ms.kids <-
        social.age.dist.ms.kids[, list(mean.age = sum(proportion.cases *
                                                          lower.age.limit)),
                           by = list(R0)]

    meanage.classic.ms.kids <- meanage.classic.ms.kids[, model := "classic"]
    meanage.hom.ms.kids <- meanage.hom.ms.kids[, model := "homogeneous"]
    meanage.social.ms.kids <- meanage.social.ms.kids[, model := "social"]
    meanage.ms.kids <-
        rbind(meanage.classic.ms.kids, meanage.hom.ms.kids, meanage.social.ms.kids)

    p <- ggplot(meanage.ms.kids, aes(x = R0, y = mean.age, color = model))
    p <- p + geom_line(lwd = 1.2)
    p <- p + geom_hline(yintercept = meanage.data.ms, lwd = 1.2, color = "black")
    p <- p + theme_bw(20)
    p <- p + theme(legend.position = "bottom")
    p <- p + scale_color_brewer(palette = "Set1")
    p <- p + scale_y_continuous("mean age at infection estimate")
    p <- p + scale_x_continuous(expression(R[0]))
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("meanage_ms_eq_kids_model", ext, sep = "."), p)
    }

    hom.age.dist.ms.kids <- hom.age.dist.ms.kids[, model := "homogeneous"]
    social.age.dist.ms.kids <- social.age.dist.ms.kids[, model := "social"]

    age.dist.ms.kids <- rbind(hom.age.dist.ms.kids, social.age.dist.ms.kids)

    p <- ggplot(age.dist.ms.kids, aes(x = R0, y = proportion.cases,
                                 color = factor(lower.age.limit)))
    p <- p + geom_line(lwd = 1.2)
    p <- p + theme_bw(20)
    p <- p + theme(legend.position = "bottom")
    p <- p +
        scale_color_brewer("age group", palette = "Dark2",
                           labels = limits.to.agegroups(ms.ew.age[, lower.age.limit]))
    p <- p + facet_grid(model ~ ., scales = "free")
    p <- p + scale_y_continuous("proportion of cases")
    p <- p + scale_x_continuous(expression(R[0]))
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("prop_cases_ms_kids_eq_model", ext, sep = "."), p)
    }

    ## now, vaccination

    years <- seq(1944, 2010)

    ## polymod: approximate R0 from mean age of infection
    R0.ms.kids.social <- approx(meanage.social.ms.kids[, mean.age],
                                meanage.social.ms.kids[, R0], meanage.data.ms)$y
    R0.ms.kids.hom <- approx(meanage.hom.ms.kids[, mean.age],
                             meanage.hom.ms.kids[, R0], meanage.data.ms)$y
    R0.ms.kids.classic <- classic.life.expectancy / meanage.data.ms

    mmr <- estimate.immunity.mmr(years = years, age.limits = lower.limits)

    social.age.dist.ms.kids.vacc <-
        age.distribution(contact.matrix.kids, R0 = R0.ms.kids.social,
                         vaccination = mmr, years = years)
    hom.age.dist.ms.kids.vacc <-
        age.distribution(hom.matrix.kids, R0 = R0.ms.kids.hom,
                         vaccination = mmr, years = years)

    social.age.dist.ms.kids.vacc <-
        age.distribution(contact.matrix.kids, R0 = R0.ms.kids.social,
                         vaccination = mmr, years = years)
    hom.age.dist.ms.kids.vacc <-
        age.distribution(hom.matrix.kids, R0 = R0.ms.kids.hom,
                         vaccination = mmr, years = years)

    ## mean age at infection
    meanage.hom.ms.kids.vacc <-
        hom.age.dist.ms.kids.vacc[,
                             list(mean.age = sum(proportion.cases * lower.age.limit)),
                             by = list(year)]
    meanage.classic.ms.kids.vacc.all <-
        data.table(year = years,
                   mean.age = classic.life.expectancy /
                       (R0.ms.kids.classic * (1 - m.mmr[, rel.vaccinated])))
    ## meanage.classic.ms.kids.vacc.newborn <-
    ##     data.table(year = meanage.hom.ms.kids.vacc[, year],
    ##                mean.age = classic.life.expectancy / (R0.ms.kids.classic *
    ##                    (1 - c(rep(0, min(vaccine.ew[, year]) - min(years)),
    ##                           vaccine.ew[year %in% years, MCV1 / 100]))))
    meanage.social.ms.kids.vacc <-
        social.age.dist.ms.kids.vacc[, list(mean.age = sum(proportion.cases *
                                                               lower.age.limit)),
                                     by = list(year)]
    meanage.data.ms.kids.vacc <-
        data.table(year = years,
                   mean.age = sapply(years, function(x) {
                       ms.ew.age[year == x, wtd.mean(lower.age.limit, abs.incidence)]
                   }))

    meanage.classic.ms.kids.vacc.all <-
        meanage.classic.ms.kids.vacc.all[, model := "classic"]
    ## meanage.classic.ms.kids.vacc.newborn <-
    ##     meanage.classic.ms.kids.vacc.newborn[, model := "classic (birth vaccinations)"]
    meanage.hom.ms.kids.vacc <- meanage.hom.ms.kids.vacc[, model := "homogenous"]
    meanage.social.ms.kids.vacc <- meanage.social.ms.kids.vacc[, model := "social"]
    meanage.data.ms.kids.vacc <- meanage.data.ms.kids.vacc[, model := "data"]

    meanage.ms.kids.vacc <-
        rbind(meanage.classic.ms.kids.vacc.all, meanage.hom.ms.kids.vacc,
              meanage.social.ms.kids.vacc, meanage.data.ms.kids.vacc)

    p <- ggplot(meanage.ms.kids.vacc, aes(x = year, y = mean.age, color = model))
    p <- p + geom_line(lwd = 1.2)
    p <- p + theme_bw(20)
    p <- p + theme(legend.position = "bottom")
    p <- p + scale_color_brewer(palette = "Set1")
    p <- p + scale_y_continuous("mean age at infection estimate")
    p <- p + scale_x_continuous(expression(R[0]))
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("meanage_ms_kids_vacc_eq_model", ext, sep = "."), p)
    }

    hom.age.dist.ms.kids.vacc <- hom.age.dist.ms.kids.vacc[, model := "homogeneous"]
    social.age.dist.ms.kids.vacc <- social.age.dist.ms.kids.vacc[, model := "social"]

    ## get proportion of cases in the real world
    ms.ew.age.total <- ms.ew.age[, list(abs.pop.incidence = sum(abs.incidence),
                                        abs.population = sum(population)),
                                 by = list(year)]
    setkey(ms.ew.age.total, year)
    ms.ew.age <- merge(ms.ew.age, ms.ew.age.total, by = "year")
    ms.ew.age <- ms.ew.age[, proportion.cases := abs.incidence / abs.pop.incidence]
    ms.ew.age <- ms.ew.age[, abs.pop.incidence := NULL]
    ms.ew.age <- ms.ew.age[, abs.population := NULL]
    ms.ew.age <- ms.ew.age[, model := "data"]

    age.dist.ms.kids.vacc <-
        rbind(hom.age.dist.ms.kids.vacc, social.age.dist.ms.kids.vacc,
              ms.ew.age[year < max(years)], use.names = T, fill = TRUE)

    p <- ggplot(age.dist.ms.kids.vacc, aes(x = year, y = proportion.cases,
                                 color = factor(lower.age.limit)))
    p <- p + geom_line(lwd = 1.2)
    p <- p + theme_bw(20)
    p <- p + theme(legend.position = "bottom")
    p <- p +
        scale_color_brewer("age group", palette = "Dark2",
                           labels = limits.to.agegroups(ms.ew.age[, lower.age.limit]))
    p <- p + facet_grid(model ~ .)
    p <- p + scale_y_continuous("proportion of cases")
    p <- p + scale_x_continuous("year")
    if (is.na(ext)) {
        p
    } else {
        ggsave(paste("prop_cases_ms_kids_vacc_eq_model", ext, sep = "."), p)
    }
}

