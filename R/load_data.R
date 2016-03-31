##' Load measles England & Wales data
##'
##' @param path Path in which to find the raw data
##' @author Sebastian Funk
##' @import reshape2
loadBryanData <- function(path) {

    load("ewm")
    
    ms.ew.4494 <- ewM4494
    times <- as.numeric(rownames(ms.ew.4494)) + 1900
    years <- ceiling(times)
    weeks <- sapply(unique(years), function(x) sum(times > (x - 1) & times <= x))
    weeks <- unlist(sapply(weeks, function(x) rep(x, x)))
    
    ms.ew.4494$time <- parse_date_time(paste((years - 1), round((times - years + 1) * weeks), 1, sep = "-"), "%Y %W %w")
    setnames(ms.ew.4494, "rn", "time")
    ms.ew.4494[, time := as.numeric(as.character(time)) + 1900]

    births.ew.4494 <- data.table(ewB4494, keep.rownames = TRUE)
    setnames(births.ew.4494, "rn", "time")
    births.ew.4494[, time := as.numeric(as.character(time)) + 1900]

    deaths.ew.4494 <- data.table(ewD4494, keep.rownames = TRUE)
    setnames(deaths.ew.4494, "rn", "time")
    deaths.ew.4494[, time := as.numeric(as.character(time)) + 1900]

    pop.ew.4494 <- data.table(ewP4494, keep.rownames = TRUE)
    setnames(pop.ew.4494, "rn", "time")
    pop.ew.4494[, time := as.numeric(as.character(time)) + 1900]

    coordinates.ew.4494 <- data.table(ewXY4494)

    save(ms.ew.4494, births.ew.4494, deaths.ew.4494, pop.ew.4494,
         coordinates.ew.4494, file = "ew_4494.rdata")
}

##' Load age-specific measles England & Wales data
##'
##' @param path Path in which to find the raw data
##' @author Sebastian Funk
##' @import reshape2
loadMeaslesEWAge <- function(path) {

    try_require("reshape2")

    ## bands into which we have the data split
    yearbands <- c("4455", "5675", "7685", "8686", "8792", "9394", "9597")
    bandnames <- sapply(yearbands, function(x) {
        paste("19", substr(x, 1, 2), "-", substr(x, 3, 4), sep = "")
    })
    ## read age-specific case data, 1944-1997 England & Wales
    ms.ew.l <- list()
    for (i in 1:length(yearbands)) {
        ms.ew.l[[bandnames[i]]] <-
            data.table(read.csv(paste(path, "/ages_", yearbands[i], ".csv",
                                      sep = "")))
    }

    ## time variables that we might be interested in
    time.variables <- c("year", "quarter", "week")
    ## time variables that they all have in common (for merging later)
    common.time.columns <- time.variables
    ## variable to hold a list of the different datasets
    ms.ew.age.l <- list()
    ## extract year, quarter and age group data (ignoring ones of unknown age)
    for (i in 1:length(ms.ew.l)) {
        ## extract age columns (they begin with an "x")
        age.columns <- names(ms.ew.l[[i]])[grep("^X", names(ms.ew.l[[i]]))]
        ## extract time columns
        time.columns <- intersect(time.variables, tolower(names(ms.ew.l[[i]])))
        ## convert colum names of time columns to lower case
        time.column.numbers <- match(time.columns, tolower(names(ms.ew.l[[i]])))
        setnames(ms.ew.l[[i]], time.column.numbers,
                 tolower(names(ms.ew.l[[i]])[time.column.numbers]))
        ## update common time columns
        common.time.columns <- intersect(common.time.columns, time.columns)
        ms.ew.l[[i]] <- ms.ew.l[[i]][, c(time.columns, age.columns), with = F]
        ms.ew.age.l[[i]] <- data.table(melt(ms.ew.l[[i]], id.vars = time.columns))
        setnames(ms.ew.age.l[[i]], "value", "abs.incidence")
        ## convert all columns to lower case
        ## only select variables that we are interested in

        ms.ew.age.l[[i]] <- ms.ew.age.l[[i]][, lower.age.limit := NA_real_]
        ## extract number from original age data if there is a lower
        ## case m in the variable name, it's talking about month so we
        ## get the age by dividing by 12
        ms.ew.age.l[[i]] <-
            ms.ew.age.l[[i]][grep("X.*m", variable),
                             lower.age.limit := as.numeric(sub("\\..*", "",
                                             sub("X", "", variable)))/12]
        ## extract number from original age data if there is no lower
        ## case m in the variable name, it's talking about years so we
        ## get the age directly
        ms.ew.age.l[[i]] <-
            ms.ew.age.l[[i]][grep("^X[^m]*$", variable),
                             lower.age.limit := as.numeric(sub("\\..*", "",
                                             sub("X", "", variable)))]

        ## sum incidences for each age group
        ms.ew.age.l[[i]] <-
            ms.ew.age.l[[i]][, list(abs.incidence = sum(abs.incidence)),
                             by = c(time.columns, "lower.age.limit")]

        ## if year is given in [0, 100] add 1900
        ms.ew.age.l[[i]] <-
            ms.ew.age.l[[i]][year <= 100, year := as.integer(year + 1900)]

        ## update common age bands
        if (i == 1) {
            ## first iteration, we start with the age bands in that group
            common.lower.age.limits <- unique(ms.ew.age.l[[i]][, lower.age.limit])
        } else {
            common.lower.age.limits <-
                intersect(common.lower.age.limits,
                          unique(ms.ew.age.l[[i]][, lower.age.limit]))
        }
        setkeyv(ms.ew.age.l[[i]], c(time.columns, "lower.age.limit"))
    }
    names(ms.ew.age.l) <- names(ms.ew.l)
    ms.ew.list.pre1998 <- ms.ew.age.l
    save(ms.ew.list.pre1998, file = "ms_ew_age_pre1998.RData")

    ## construct aggregate data frame
    for (i in 1:length(ms.ew.age.l)) {
        temp.ms.ew.age <-
            ms.ew.age.l[[i]][, c(common.time.columns, "abs.incidence",
                                 "lower.age.limit"), with = F]

        ## find in which interval of common lower age limits the lower
        ## age limits are
        intervals <- findInterval(temp.ms.ew.age[, lower.age.limit],
                                  common.lower.age.limits)
        ## convert them to limits
        lower.limits <- common.lower.age.limits[intervals]
        ## assign them
        temp.ms.ew.age <- temp.ms.ew.age[, lower.age.limit := lower.limits]
        temp.ms.ew.age <-
            temp.ms.ew.age[, list(abs.incidence = sum(abs.incidence)),
                           by = c(common.time.columns, "lower.age.limit")]

        if (i == 1) {
            ms.ew.age <- temp.ms.ew.age
        } else {
            ms.ew.age <- rbind(ms.ew.age, temp.ms.ew.age)
        }
    }

    ## get ECDC data (for recent E&W)
    data(ms_europe_cases)
    ## reduce to England and Wales
    ms.ew.cases.post2000 <-
        ms.europe.cases[ReportingCountry == "UK" &
                            place != "UKM" & place != "UKN"]
    ## convert age to numeric
    ms.ew.cases.post2000 <- ms.ew.cases.post2000[Age != "NULL"]
    ms.ew.cases.post2000 <-
        ms.ew.cases.post2000[, numeric.age := as.integer(as.character(Age))]

    ## find in which interval of common lower age limits the lower
    ## age limits are
    intervals <- findInterval(ms.ew.cases.post2000[, numeric.age],
                              common.lower.age.limits)
    ## convert them to limits
    lower.limits <- common.lower.age.limits[intervals]
    ## assign them
    ms.ew.cases.post2000 <-
        ms.ew.cases.post2000[, lower.age.limit := lower.limits]
    ms.ew.post2000 <-
        ms.ew.cases.post2000[, list(abs.incidence = length(numeric.age)),
                             by = c(common.time.columns, "lower.age.limit")]

    ms.ew.age <- ms.ew.age[, source := "Registrar"]
    ms.ew.post2000 <- ms.ew.post2000[, source := "ECDC"]
    ms.ew.age <- rbind(ms.ew.age, ms.ew.post2000)

    ## get population data for relative incidence
    data(pop_ew_age)

    ## find common lower age limits
    common.lower.age.limits <- intersect(unique(pop.ew.age[, lower.age.limit]),
                                         unique(ms.ew.age[, lower.age.limit]))

    ## find in which interval of common lower age limits the population
    ## lower age limits are
    pop.intervals <- findInterval(pop.ew.age[, lower.age.limit],
                                  common.lower.age.limits)
    ## convert them to limits
    lower.limits <- common.lower.age.limits[pop.intervals]
    ## assign them
    pop.ew.age <- pop.ew.age[, lower.age.limit := lower.limits]
    ## synthesise them
    pop.ew.age <- pop.ew.age[, list(population = sum(population)),
                             by = list(year, lower.age.limit)]

    ## find in which interval of common lower age limits the chickenpox
    ## lower age limits are
    ms.intervals <- findInterval(ms.ew.age[, lower.age.limit],
                                 common.lower.age.limits)
    ## convert them to limits
    lower.limits <- common.lower.age.limits[ms.intervals]
    ## assign them
    ms.ew.age <- ms.ew.age[, lower.age.limit := lower.limits]
    ## synthesise them
    ms.ew.age <- ms.ew.age[, list(abs.incidence = sum(abs.incidence)),
                           by = list(year, lower.age.limit)]

    ## fill missing population years with later/earlier reported data
    missing.years <- setdiff(unique(ms.ew.age[, year]),
                             unique(pop.ew.age[, year]))
    ## forward loop (fill in future missing data from past years)
    early.missing.years <-
        missing.years[missing.years < min(unique(pop.ew.age[, year]))]
    earliest.present <- pop.ew.age[year == min(year)]
    for (early.year in early.missing.years) {
        pop.ew.age <-
            rbind(earliest.present[, year := early.year], pop.ew.age)
    }
    missing.years <- setdiff(missing.years, early.missing.years)
    for (missing.year in missing.years) {
        closest.year <- max(pop.ew.age[year < missing.year, year])
        closest.present <- pop.ew.age[year == closest.year]
        pop.ew.age <-
            rbind( pop.ew.age, closest.present[, year := missing.year])
    }
    setkey(pop.ew.age, year, lower.age.limit)
    setkey(ms.ew.age, year, lower.age.limit)

    ## now, merge them
    ms.ew.age <- merge(ms.ew.age, pop.ew.age)
    ms.ew.age <- ms.ew.age[, rel.incidence := abs.incidence / population]

    setkeyv(ms.ew.age, c(common.time.columns, "lower.age.limit"))
    save(ms.ew.age, file = "ms_ew_age.RData")

}

##' Load measles Europe data
##'
##' @param path Path in which to find the raw data
##' @author Sebastian Funk
##' @import ISOweek
loadMeaslesEurope <- function(path) {

    ## ECDC case data
    ms.europe.cases <-
        data.table(read.csv(paste(path, "/ecdc_measles.csv", sep = "")))

    ## get best possible spatial resolution
    ms.europe.cases <- ms.europe.cases[, place := PlaceOfNotification]
    ms.europe.cases <- ms.europe.cases[place == "UNK", place := PlaceOfResidence]
    ms.europe.cases <- ms.europe.cases[place == "NULL", place := PlaceOfResidence]
    ms.europe.cases <- ms.europe.cases[place == "UNK", place := ReportingCountry]
    ms.europe.cases <- ms.europe.cases[place == "NULL", place := ReportingCountry]

    ## set date
    ms.europe.cases <- ms.europe.cases[, date := as.Date(DateUsedForStatisticsISO)]
    ms.europe.cases <-
        ms.europe.cases[grep("W", DateUsedForStatisticsISO),
                        date := ISOweek2date(paste(DateUsedForStatisticsISO,
                             1, sep = "-"))]
    ## set date parts
    ms.europe.cases <-
        ms.europe.cases[, week := as.integer(substr(DateUsedForStatisticsWeek, 6, 7))]
    ms.europe.cases <-
        ms.europe.cases[, month := as.integer(DateUsedForStatisticsMonth)]
    ms.europe.cases <-
        ms.europe.cases[, quarter := as.integer(DateUsedForStatisticsQuarter)]
    ms.europe.cases <-
        ms.europe.cases[, year := as.integer(DateUsedForStatisticsYear)]

    setkey(ms.europe.cases, date)
    save(ms.europe.cases, file = "measles_europe_cases.RData")

    ## spatial data
}

##' Load measles WHO data
##'
##' @param path Path in which to find the raw data
##' @author Sebastian Funk
##' @import lubridate
loadMeaslesWHO <- function(path) {

    try_require("lubridate")

    ## denominator for the relative incidence data
    denominator <- 1e+5

    ## load data and define useful variables
    ms.who.cases <-
        data.table(read.csv(paste(path, "/who_measles.csv", sep = "")))
    ms.who.cases <-
        ms.who.cases[, rash.date := as.Date(DRash, "%d/%m/%Y %H:%M")]
    ms.who.cases <- ms.who.cases[, week := isoweek(rash.date)]
    ms.who.cases <- ms.who.cases[, year := year(rash.date)]
    ms.who.cases <-
        ms.who.cases[, ISOweek := paste(year, "-W",
                               sprintf("%02i", week), "-1", sep="")]
    ms.who.cases <- ms.who.cases[, week.date := ISOweek2date(ISOweek)]

    setnames(ms.who.cases, c("ISO3", "Region", "AgeAtRashOnset", "NumOfVaccines"),
             c("country", "region", "age", "nvaccines"))

    save(ms.who.cases, file = "ms_who_cases.RData")

    ## store only those that have a sensible data (100 seems to
    ## indicate unknown)
    ms.who.age.cases <- ms.who.cases[age >= 0 & age < 100]

    ## get population data (for relative incidences)
    data(pop_world_age)
    common.lower.age.limits <-
        intersect(unique(pop.world.age[, lower.age.limit]),
                  unique(ms.who.age.cases[, lower.age.limit]))

    ## find the age groups that cases are in
    ms.intervals <- findInterval(ms.who.age.cases[, age],
                                 common.lower.age.limits)
    ## convert them to limits
    ms.lower.limits <- common.lower.age.limits[intervals]

    ## assign
    ms.who.age.cases <- ms.who.age.cases[, lower.age.limit := ms.lower.limits]
    save(ms.who.age.cases, file = "ms_who_age_cases.RData")

    ## find the age groups that cases are in
    pop.intervals <- findInterval(pop.world.age[, lower.age.limit],
                                  common.lower.age.limits)
    ## convert them to limits
    pop.lower.limits <- common.lower.age.limits[pop.intervals]

    ## assign
    pop.world.age <- pop.world.age[, lower.age.limit := pop.lower.limits]

    ## reduce
    pop.world.age <-
        pop.world.age[, list(population = sum(both.sexes.population)),
                      by = list(country, lower.age.limit)]

    ms.who.age <-
        ms.who.age.cases[, list(abs.incidence = length(DRYEAR)),
                         by = list(country, week.date, year, week,
                         lower.age.limit)]
    setkey(ms.who.age, country, lower.age.limit)
    setkey(pop.world.age, country, lower.age.limit)

    ms.who.age <- merge(ms.who.age, pop.world.age, allow.cartesian = T)
    ms.who.age <- ms.who.age[, rel.incidence := abs.incidence * denominator /
                                 population]

    ms.who.age[, population := NULL]
    setkey(ms.who.age, country, year, week)

    save(ms.who.age, file = "ms_who_age.RData")
}

##' Load measles serology data
##'
##' @param path Path in which to find the raw data
##' @author Sebastian Funk
loadMeaslesSerology <- function(path) {

    ## get 2000 serology data
    ms.sero <-
        data.table(read.csv(paste(path, "/measles_serology.csv", sep = "")))
    setnames(ms.sero, names(ms.sero), tolower(names(ms.sero)))

    ms.sero <- ms.sero[, exact.age := as.numeric(as.character(exact.age))]
    ms.sero <- ms.sero[, non.negative := 1]
    ms.sero <- ms.sero[stdres == "NEG", non.negative := 0]

    setkey(ms.sero, country, age1)
    save(ms.sero, file = "ms_sero.RData")
}

##' Load measles serology data from Fine(1950)
##'
##' @param path Path in which to find the raw data
##' @author Sebastian Funk
loadMeaslesSerologyFine1950 <- function(path) {

    ## get 2000 serology data
    ms.sero.temp <-
        data.table(read.csv(paste(path, "/measles_serology_1950_fine.csv", sep = "")))
    ms.sero.temp <- rbind(data.table(age = 0, immunity = 1), ms.sero)

    ms.approx <- approx(ms.sero.temp[, age], ms.sero.temp[, immunity], 1:20)
    ms.sero.fine <- data.table(matrix(unlist(ms.approx), ncol = 2))
    setnames(ms.sero.fine, 1:2, names(ms.sero.temp))

    ms.sero.fine <- ms.sero.fine[, immunity := immunity / max(immunity)]

    setkey(ms.sero.fine, age)

    save(ms.sero.fine, file = "ms_sero_fine_1950.RData")
}

##' Load MMR world data
##'
##' @param path Path in which to find the raw data
##' @author Sebastian Funk
##' @import reshape2 XLConnect
loadMMRWorld <- function(path) {

    clean_mcv_data <- function(workbook, sheet)
    {
        mcv.workbook <- loadWorkbook(workbook)
        mcv.table <- data.table(readWorksheet(mcv.workbook, sheet))
        setnames(mcv.table, names(mcv.table), tolower(names(mcv.table)))
        setnames(mcv.table, "cname", "country")
        mcv.empty.columns <- apply(mcv.table, 2, function(x) {all(is.na(x))})
        if (any(mcv.empty.columns))
        {
            mcv.table <-
                mcv.table[, !mcv.empty.columns, with = F]
        }
        mcv.region <-
            grep("region", colnames(mcv.table), value = TRUE)
        setnames(mcv.table, mcv.region, "who_region")
        mcv.long <-
            data.table(melt(mcv.table,
                            id.vars = c("who_region", "iso_code", "country", "vaccine")))
        setnames(mcv.long, "variable", "year")
        setnames(mcv.long, "value", "uptake")
        mcv.long <- mcv.long[, year := as.integer(gsub("x", "", year))]
        mcv.long <- mcv.long[!is.na(uptake)]
        mcv.long <- clean_countries(mcv.long)
        return(mcv.long)
    }

    mcv1.world <- clean_mcv_data(paste(path, "/coverage_series.xls", sep = ""), "MCV")
    mcv2.world <- clean_mcv_data(paste(path, "/coverage_series.xls", sep = ""), "MCV2")

    mcv.world <- rbind(mcv1.world, mcv2.world)
    save(mcv.world, file = "mcv_world.RData")

    mcv1.world.estimate <- clean_mcv_data(paste(path, "/coverage_estimates_series.xls", sep = ""), "MCV")
    mcv2.world.estimate <- clean_mcv_data(paste(path, "/coverage_estimates_series.xls", sep = ""), "MCV2")

    mcv.world.estimate <- rbind(mcv1.world.estimate, mcv2.world.estimate)
    save(mcv.world.estimate, file = "mcv_world_estimate.RData")

}

##' Clean country names
##'
##' @param table data table which contains country names
##' @param column the column with country names
##' @return cleaned data table
##' @author Sebastian Funk
clean_countries <- function(table, column = "country")
{
    dt <- data.table(table)
    
    dt[get(column) == "Czech Republic (the)", paste(column) := "Czech Republic"]
    dt[get(column) == "Netherlands (the)", paste(column) := "Netherlands"]
    dt[get(column) == "Moldova", paste(column) := "Rep of Moldova"]
    dt[get(column) == "Republic of Moldova (the)", paste(column) := "Rep of Moldova"]
    dt[get(column) == "Russia", paste(column) := "Russian Federation"]
    dt[get(column) == "Russian Federation (the)",
       paste(column) := "Russian Federation"]
    dt[get(column) == "Macedonia", paste(column) := "FYR Macedonia"]
    dt[get(column) == "The former Yugoslav Republic of Macedonia",
       paste(column) := "FYR Macedonia"]
    dt[get(column) == "United Kingdom of Great Britain and Northern Ireland (the)",
       paste(column) := "United Kingdom"]
    
    return(dt)
}

##' Load MMR England & Wales data
##'
##' @param path Path in which to find the raw data
##' @author Sebastian Funk
loadMMREW <- function(path){

    ## join vaccination data
    vaccine.ew <-
        data.table(read.csv(paste(path, "/vaccine_ew.csv", sep = "")))
    setnames(vaccine.ew, "X", "year")

    save(vaccine.ew, file = "vaccine_ew.RData")
}

##' Load Chickenpox England & Wales age data
##'
##' @param path Path in which to find the raw data
##' @author Sebastian Funk
loadChickenpoxEWAge <- function(path) {

    ## denominator for the absolute incidence data
    denominator <- 1e+5
    ## age groups
    agegroups <- c("0_4", "5_14", "15_44", "45_64", "65_")
    cpx.ew.age <- data.table(week = integer(0), year = integer(0),
                             agegroup = character(0),
                             rel.incidence = numeric(0))

    ## read in files
    for (age in agegroups) {
        temp.cpx.ew.age <-
            data.table(read.csv(paste(path, "/cpx_data_", age, ".csv", sep = ""),
                                header = T, row.names = 1))

        for (i in 1:ncol(temp.cpx.ew.age)) {
            year <- as.integer(sub("X", "", names(temp.cpx.ew.age)[i]))
            incidence <- unlist(unname(temp.cpx.ew.age[, i, with = F]))
            cpx.ew.age <-
                rbind(cpx.ew.age,
                      data.frame(week = as.integer(rownames(temp.cpx.ew.age)),
                                 year = year,
                                 agegroup = age,
                                 rel.incidence = incidence / denominator))
        }
    }

    ## clean up data set and define useful variables
    cpx.ew.age <- cpx.ew.age[!is.na(rel.incidence)]

    cpx.ew.age <-
        cpx.ew.age[, lower.age.limit := as.numeric(sub("_.*", "", agegroup))]

    ## get population data for absolute incidence
    data(pop_ew_age)

    ## determine lower age limits of chickenpox data
    lower.age.limits <- unique(cpx.ew.age[, lower.age.limit])
    ## find in which interval of lower age limits the population lower age limits
    pop.intervals <- findInterval(pop.ew.age[, lower.age.limit],
                                  lower.age.limits)
    ## convert them to limits
    lower.limits <- lower.age.limits[pop.intervals]
    ## assign them
    pop.ew.age <-
        pop.ew.age[, lower.age.limit := lower.limits]
    pop.ew.age <-
        pop.ew.age[, list(population = sum(population)),
                   by = list(year, lower.age.limit)]

    setkey(cpx.ew.age, year, week, lower.age.limit)

    ## fill missing population years with later/earlier reported data
    missing.years <- setdiff(unique(cpx.ew.age[, year]),
                             unique(pop.ew.age[, year]))
    ## forward loop (fill in future missing data from past years)
    early.missing.years <-
        missing.years[missing.years < min(unique(pop.ew.age[, year]))]
    earliest.present <- pop.ew.age[year == min(year)]
    for (early.year in early.missing.years) {
        pop.ew.age <-
            rbind(earliest.present[, year := early.year], pop.ew.age)
    }
    missing.years <- setdiff(missing.years, early.missing.years)
    for (missing.year in missing.years) {
        closest.year <- max(pop.ew.age[year < missing.year, year])
        closest.present <- pop.ew.age[year == closest.year]
        pop.ew.age <-
            rbind( pop.ew.age, closest.present[, year := missing.year])
    }
    setkey(pop.ew.age, year, lower.age.limit)

    ## merge chickenpox data with population data
    cpx.ew.age <- merge(cpx.ew.age, pop.ew.age)
    cpx.ew.age <- cpx.ew.age[, abs.incidence :=
                                 round(rel.incidence * population)]

    cpx.ew.age <- cpx.ew.age[, agegroup := NULL]
    setkey(cpx.ew.age, year, week, lower.age.limit)

    save(cpx.ew.age, file = "cpx_ew_age.RData")

    ## overall data
    cpx.ew <- data.table(week = integer(0), year = integer(0),
                         rel.incidence = numeric(0))
    ## read in files
    temp.cpx.ew <-
        data.table(read.csv(paste(path, "/cpx_data_all.csv", sep = ""),
                            header = T, row.names = 1))

    for (i in 1:ncol(temp.cpx.ew)) {
        year <- as.integer(sub("X", "", names(temp.cpx.ew)[i]))
        incidence <- unlist(unname(temp.cpx.ew[, i, with = F]))
        cpx.ew <-
            rbind(cpx.ew, data.frame(week = as.integer(rownames(temp.cpx.ew)),
                                     year = year,
                                     rel.incidence = incidence / denominator))
    }

    ## clean up data set and define useful variables
    cpx.ew <- cpx.ew[!is.na(rel.incidence)]

    setkey(cpx.ew, year, week)

    ## get population data for absolute incidence
    data(pop_ew)

    ## fill missing population years with later/earlier reported data
    missing.years <- setdiff(unique(cpx.ew[, year]),
                             unique(pop.ew[, year]))
    ## forward loop (fill in future missing data from past years)
    early.missing.years <-
        missing.years[missing.years < min(unique(pop.ew[, year]))]
    earliest.present <- pop.ew[year == min(year)]
    for (early.year in early.missing.years) {
        pop.ew <-
            rbind(earliest.present[, year := early.year], pop.ew)
    }
    missing.years <- setdiff(missing.years, early.missing.years)
    for (missing.year in missing.years) {
        closest.year <- max(pop.ew[year < missing.year, year])
        closest.present <- pop.ew[year == closest.year]
        pop.ew <-
            rbind(pop.ew, closest.present[, year := missing.year])
    }
    setkey(pop.ew, year)

    ## merge chickenpox data with population data
    cpx.ew <- merge(cpx.ew, pop.ew)
    cpx.ew <- cpx.ew[, abs.incidence :=
                         round(rel.incidence * population)]
    setkey(cpx.ew, year, week)

    save(cpx.ew, file = "cpx_ew.RData")
}

##' Load Chickenpox serology data
##'
##' @param path Path in which to find the raw data
##' @author Sebastian Funk
loadChickenpoxSerology <- function(path) {

    cpx.sero <-
        data.table(read.csv(paste(path, "/cpx_serology.csv", sep = "")))

    cpx.sero <- melt(cpx.sero, id.vars = "year")
    cpx.sero <- cpx.sero[, lower.age.limit := as.numeric(gsub("X", "", variable))]
    cpx.sero <- cpx.sero[, antibodies := as.integer(gsub("/.*$", "", value))]
    cpx.sero <- cpx.sero[, samples := as.integer(gsub("^.*/", "", value))]

    cpx.sero <- cpx.sero[, variable := NULL]
    cpx.sero <- cpx.sero[, value := NULL]

    setkey(cpx.sero, year, lower.age.limit)

    save(cpx.sero, file = "cpx_sero.RData")

}

##' Load Copenhagen disease data
##'
##' @param path Path in which to find the raw data
##' @author Sebastian Funk
##' @import reshape2
loadCopenhagenDisease <- function(path) {

    ## load csv file
    cph <- data.table(read.csv(paste(path, "/4other_cph_diseases.csv",
                                     sep = "")))

    ## get columns in csv file that correspond to diseases
    disease.columns <- which(grepl("^[^(X|tot|date)]", names(cph)))
    disease.names <- tolower(names(cph)[disease.columns])
    ## for the loop later
    disease.columns <- c(disease.columns, ncol(cph) + 1)
    date.column <- which(grepl("^date", names(cph)))

    ## get date
    setnames(cph, names(cph), as.character(unlist(unname(cph[1]))))
    setnames(cph, date.column, "date")
    cph <- cph[-1]

    cph.disease.total.l <- list()
    for (i in 1:length(disease.names)) {
        cph.disease.total.l[[disease.names[i]]] <-
            cph[, c("date", "tot"), with = F]
        cph.disease.total.l[[i]] <-
            cph.disease.total.l[[i]][, cases := as.integer(as.character(tot))]
        cph.disease.total.l[[i]] <-
            cph.disease.total.l[[i]][!is.na(cases)]
        cph.disease.total.l[[i]] <-
            cph.disease.total.l[[i]][, name := disease.names[i]]
        cph.disease.total.l[[i]] <-
            cph.disease.total.l[[i]][, list(date, cases, name)]
    }

    if (length(cph.disease.total.l) > 0) {
        cph.disease.total <- cph.disease.total.l[[1]]
        for (disease in names(cph.disease.total.l)[-1]) {
            cph.disease.total <-
                rbind(cph.disease.total, cph.disease.total.l[[disease]])
        }
    }

    setkey(cph.disease.total, date)
    save(cph.disease.total, file = "cph_disease.RData")

    cph.disease.total <- cph.disease.total.l[[1]]
    for (i in names(cph.disease.total.l)[-1]) {
        cph.disease.total <- rbind(cph.disease.total, cph.disease.total.l[[i]])
    }

    ## get different diseases
    cph.disease.age.l <- list()
    for (i in seq(1, length(disease.columns) - 1)) {
        cph.disease.age.l[[disease.names[[i]]]] <-
            cph[, c(date.column,
                    seq(disease.columns[i], disease.columns[i + 1] - 1)),
                with = F]
        ## get age groups
        cph.disease.age.l[[i]] <-
            cph.disease.age.l[[i]][, name := disease.names[i]]
        cph.disease.age.l[[i]] <-
            data.table(melt(cph.disease.age.l[[i]], id.vars=c("date", "name")))
        cph.disease.age.l[[i]] <-
            cph.disease.age.l[[i]][value != ""]
        cph.disease.age.l[[i]] <-
            cph.disease.age.l[[i]][, cases := as.integer(as.character(value))]
        cph.disease.age.l[[i]] <-
            cph.disease.age.l[[i]][, date := as.Date(date, format="%m/%d/%Y")]
        setnames(cph.disease.age.l[[i]], "variable", "agegroup")
        cph.disease.age.l[[i]] <-
            cph.disease.age.l[[i]][grepl("[[:digit:]]", agegroup)]
        cph.disease.age.l[[i]] <-
            cph.disease.age.l[[i]][, lower.age.limit := 0]
        agegroups <- as.character(unique(cph.disease.age.l[[i]][, agegroup]))
        for (current.agegroup in agegroups) {
            if (!grepl("^<", current.agegroup)) {
                digits <- regexpr("[[:digit:]]+", current.agegroup)
                if (digits > 0) {
                    limit <- regmatches(current.agegroup, digits)
                    cph.disease.age.l[[i]] <-
                        cph.disease.age.l[[i]][agegroup == current.agegroup,
                                               lower.age.limit := as.numeric(limit)]
                }
            }
        }
        ## sum within age groups
        cph.disease.age.l[[i]] <-
            cph.disease.age.l[[i]][, list(cases = sum(cases)),
                                   by = list(date, name,
                                   lower.age.limit)]
        ## clean
        cph.disease.age.l[[i]] <- cph.disease.age.l[[i]][!is.na(cases)]
    }


    if (length(cph.disease.age.l) > 0) {
        cph.disease.age <- cph.disease.age.l[[1]]
        for (disease in names(cph.disease.age.l)[-1]) {
            cph.disease.age <-
                rbind(cph.disease.age, cph.disease.age.l[[disease]])
        }
    }

    setkey(cph.disease.age, date, name, lower.age.limit)
    save(cph.disease.age, file = "cph_disease_age.RData")

}

##' Load Measles England & Wales confirmation data
##'
##' @param path Path in which to find the raw data
##' @author Sebastian Funk
loadMeaslesEWConfirmation <- function(path) {

    ms.ew.conf <-
        data.table(read.csv(paste(path, "/hpa_confirmations.csv", sep = "")))
    setnames(ms.ew.conf, names(ms.ew.conf), tolower(names(ms.ew.conf)))
    ms.ew.conf <-
        ms.ew.conf[, percent.positive := number.confirmed / number.tested]
    ms.ew.conf <- ms.ew.conf[, quarter := as.integer(substr(quarter, 1, 1))]

    setkey(ms.ew.conf, year, quarter)

    save(ms.ew.conf, file = "ms_ew_conf.RData")
}

##' Load population England & Wales age data
##'
##' @param path Path in which to find the raw data
##' @author Sebastian Funk
##' @import reshape2
loadPopulationEWAge <- function(path, mult.factor = 1000) {

    try_require("reshape2")

    countries <- c("england", "wales")

    ## first dataset

    for (i in 1:length(countries)) {
        ## read file
        temp.pop.ew.age <-
            data.table(read.csv(paste(path, "/pop_age_", countries[i],
                                      ".csv", sep = ""),
                                header = T), country = countries[i])

        ## melt
        temp.pop.ew.age <- melt(temp.pop.ew.age, id.vars = c("year", "country"))
        ## start new data frame if first country, otherwise append
        if (i == 1) {
            pop.ew.age <- temp.pop.ew.age
        } else {
            pop.ew.age <- rbind(pop.ew.age, temp.pop.ew.age)
        }
    }
    ## convert age classses to lower age limits
    pop.ew.age <-
        pop.ew.age[grep("^X[^m]*$", variable),
                   lower.age.limit := as.numeric(sub("\\..*", "",
                                   sub("X", "", variable)))]

    ## sum by age year and age group, multiply by factor according to
    ## how the data is stored
    pop.ew.age <-
        pop.ew.age[, list(population = sum(value) * mult.factor),
                   by = list(year, lower.age.limit)]

    ## clean up
    pop.ew.age <- pop.ew.age[!is.na(lower.age.limit)]
    pop.ew.age <- pop.ew.age[!is.na(population)]

    ## second dataset

    pop.ew.age.2 <-
        data.table(read.csv(paste(path, "/pop_old.csv", sep = "")))

    pop.ew.age.2[, 2:nrow(pop.ew.age.2)] <-
        as.data.table(lapply(pop.ew.age.2[, 2:nrow(pop.ew.age.2), with = F],
                             as.numeric))
    setnames(pop.ew.age.2, "Year", "lower.age.limit")

    pop.ew.age.2 <- melt(pop.ew.age.2[c(4:8, 10:nrow(pop.ew.age.2))],
                         id.vars = "lower.age.limit", variable.name = "year",
                         value.name = "population")
    pop.ew.age.2 <- pop.ew.age.2[, year := as.integer(gsub("^X", "", year))]
    pop.ew.age.2 <-
        pop.ew.age.2[, lower.age.limit := as.integer(as.character(lower.age.limit))]
    pop.ew.age.2 <- pop.ew.age.2[, population := as.integer(mult.factor * population)]

    pop.ew.age.2[, lower.age.limit := 
                     reduce.agegroups(lower.age.limit,
                                      unique(pop.ew.age[, lower.age.limit]))]
    pop.ew.age.2 <- pop.ew.age.2[, list(population = sum(population)),
                                 by = list(year, lower.age.limit)]

    pop.ew.age <- rbind(pop.ew.age.2, pop.ew.age)

    setkey(pop.ew.age, year, lower.age.limit)

    save(pop.ew.age, file = "pop_ew_age.RData")

    pop.ew <- pop.ew.age[, list(population = sum(population)), by = list(year)]
    save(pop.ew, file = "pop_ew.RData")

}

##' Load population Europe data
##'
##' @param path Path in which to find the raw data
##' @author Sebastian Funk
loadPopulationEurope <- function(path) {

    ## load European population CSV
    pop.europe <-
        data.table(read.csv(paste(path, "/nuts_populations.csv", sep = "")))
    setnames(pop.europe, "TIME", "year")

    ## clean data
    pop.europe <- pop.europe[, place := as.character(NUTS)]
    pop.europe <- pop.europe[!(place == "")]

    ## fill missing data in both directions
    ## mark filled data in the data
    pop.europe <- pop.europe[, missing := F]
    years <- unique(pop.europe[, year])

    ## forward loop (fill in future missing data from past years)
    for (i in 2:(length(years))) {

        missing.places <- pop.europe[Value == ":" & year == years[i], place]
        pop.europe[place %in% missing.places, missing := T]
        ## if places have missing data, look up the previous year in data
        lookup.previous.year <-
            pop.europe[place %in% missing.places & year == years[i - 1], Value]

        pop.europe <- pop.europe[place %in% missing.places & year == years[i],
                                 Value := lookup.previous.year]
    }
    ## reverse loop (fill in past missing data from future years)
    for (i in rev(1:(length(years) - 1))) {

        missing.places <- pop.europe[Value == ":" & year == years[i], place]
        ## if places have missing data, look up the next year in data
        pop.europe[place %in% missing.places, missing := T]
        lookup.next.year <-
            pop.europe[place %in% missing.places & year == years[i + 1], Value]

        pop.europe <- pop.europe[place %in% missing.places & year == years[i],
                                 Value := lookup.next.year]
    }

    pop.europe <- pop.europe[, population := as.numeric(gsub(" ", "", Value))]
    pop.europe <- pop.europe[, list(year, place, population)]
    setkey(pop.europe, year, place)

    save(pop.europe, file = "pop_europe.RData")

}

##' Load population world data 
##'
##' @param path Path in which to find the raw data
##' @author Sebastian Funk
loadPopulationWorldAge <- function(path) {

    ## load European population CSV
    pop.world.age <-
        data.table(read.csv(paste(path, "/country_age_structure.csv",
                                  sep = "")))
    setnames(pop.world.age, names(pop.world.age), tolower(names(pop.world.age)))

    ## get lower age limits
    lower.limits <- sub("[-+].*$", "", pop.world.age[, age])
    pop.world.age <- pop.world.age[, lower.age.limit := as.numeric(lower.limits)]
    pop.world.age <- pop.world.age[!is.na(lower.age.limit)]

    pop.world.age <- clean_countries(pop.world.age)
    
    save(pop.world.age, file = "pop_world_age.RData")

}

##' Load demographic England & Wales data 
##'
##' @param path Path in which to find the raw data
##' @author Sebastian Funk
loadDemographicsEW <- function(path) {

    births.ew <- data.table(read.csv(paste(path, "births_ew.csv", sep = "")))
    setnames(births.ew, names(births.ew), tolower(names(births.ew)))
    births.ew <- births.ew[, cbr := cbr1000 / 1000]
    births.ew <- births.ew[, cbr1000 := NULL]
    setkey(births.ew, year)

    save(births.ew, file = "births_ew.RData")

    deaths.ew <- data.table(read.csv(paste(path, "deaths_ew.csv", sep = "")))
    setnames(deaths.ew, names(deaths.ew), tolower(names(deaths.ew)))
    deaths.ew <- deaths.ew[, cdr := cdr1000 / 1000]
    deaths.ew <- deaths.ew[, cdr1000 := NULL]
    deaths.ew <- deaths.ew[, infdr := infdr1000 / 1000]
    deaths.ew <- deaths.ew[, infdr1000 := NULL]
    deaths.ew <- deaths.ew[, neodr := neodr1000 / 1000]
    deaths.ew <- deaths.ew[, neodr1000 := NULL]
    setkey(deaths.ew, year)

    save(deaths.ew, file = "deaths_ew.RData")
}

##' Load POLYMOD data 
##'
##' @param path Path in which to find the raw data
##' @author Sebastian Funk
loadPolymod <- function(path) {

    ## read in participant data
    participants <-
        data.table(read.csv(paste(path, "/participants_final_v3.csv", sep = "")))

    ## read in contact data
    contacts <-
        data.table(read.csv(paste(path, "/contacts_final_v2_with_weights.csv",
                                  sep = "")))

    participants[country == "BE", country := "Belgium"]
    participants[country == "DE", country := "Germany"]
    participants[country == "FI", country := "Finland"]
    participants[country == "GB", country := "United Kingdom"]
    participants[country == "IT", country := "Italy"]
    participants[country == "LU", country := "Luxembourg"]
    participants[country == "NL", country := "Netherlands"]
    participants[country == "PL", country := "Poland"]
    participants[, country := factor(country)]

    contacts[country == "BE", country := "Belgium"]
    contacts[country == "DE", country := "Germany"]
    contacts[country == "FI", country := "Finland"]
    contacts[country == "GB", country := "United Kingdom"]
    contacts[country == "IT", country := "Italy"]
    contacts[country == "LU", country := "Luxembourg"]
    contacts[country == "NL", country := "Netherlands"]
    contacts[country == "PL", country := "Poland"]
    contacts[, country := factor(country)]

    polymod <- list(participants = participants, contacts = contacts)
    save(polymod, file = "polymod.RData")

}

##' Load European maps 
##'
##' @param path Path in which to find the raw data
##' @author Sebastian Funk
##' @import maptools ggplot2
loadEuropeanMaps <- function(shapePath) {

    try_require("maptools")
    try_require("ggplot2")

    ## prepare map of Europe
    europe.shape <-
        readShapePoly(paste(shapePath, "/NUTS_RG_60M_2010.shp", sep = ""))
    europe.data <- data.table(europe.shape@data)
    setnames(europe.data, "NUTS_ID", "place")
    setkey(europe.data, "place")
    gpclibPermit()
    europe <- data.table(fortify(europe.shape, region = "NUTS_ID"))
    setnames(europe, "id", "place")
    setkey(europe, "place")
    europe <- merge(europe, europe.data, by="place")

    save(europe, file = "europe.RData")
}

##' Load MMR European timing data 
##'
##' @param path Path in which to find the raw data
##' @author Sebastian Funk
loadMMREuropeTimings <- function(path) {
    ## read in vaccination timings

    mmr.timings <-
        data.table(read.csv(paste(path, "/europe_mmr_timings.csv", sep = "")))

    for (coverage in c("mcv1", "mcv2")) {
        temp.mmr.timings <-
            unlist(strsplit(as.character(mmr.timings[, get(coverage)]), " "))
        dt.mmr.timings <- data.table(matrix(temp.mmr.timings, ncol = 2, byrow = T))

        setnames(dt.mmr.timings, "V1", "number")
        setnames(dt.mmr.timings, "V2", "unit")
        dt.mmr.timings <-
            dt.mmr.timings[, min := as.numeric(gsub("-.*$", "",
                                 as.character(dt.mmr.timings[, number])))]
        dt.mmr.timings <-
            dt.mmr.timings[, max := as.numeric(gsub("^.*-", "",
                                 as.character(dt.mmr.timings[, number])))]
        dt.mmr.timings[grep("month", unit), min := min / 12]
        dt.mmr.timings[grep("month", unit), max := max / 12]
        setnames(dt.mmr.timings, "min", paste(coverage, "min", sep = "."))
        setnames(dt.mmr.timings, "max", paste(coverage, "max", sep = "."))
        dt.mmr.timings <- dt.mmr.timings[, number := NULL]
        dt.mmr.timings <- dt.mmr.timings[, unit := NULL]
        mmr.timings <- cbind(mmr.timings, dt.mmr.timings)
    }

    mmr.timings <- clean_countries(mmr.timings)

    save(mmr.timings, file = "europe_mmr_timings.RData")
}

##' Load measles world data 
##'
##' @param path Path in which to find the raw data
##' @import XLConnect
##' @author Sebastian Funk
loadWorldMeaslesData <- function(path) {
    ## read in world measles data

    try_require("XLConnect")

    incidence.workbook <- loadWorkbook(paste(path, "/who_incidence.xls", sep = ""),
                                       create = FALSE)
    ms.world <- data.table(readWorksheet(incidence.workbook, "Measles"))
    setnames(ms.world, names(ms.world), tolower(names(ms.world)))
    ms.world[, disease := NULL]
    ms.world <-
        data.table(melt(ms.world, id.vars = c("who_region", "iso_code", "cname")))
    ms.world[, variable := as.integer(sub("^x", "", variable))]
    setnames(ms.world, "cname", "country")
    setnames(ms.world, "variable", "year")
    setnames(ms.world, "value", "cases")

    ms.world <- clean_countries(ms.world)

    save(ms.world, file = "ms_world.RData")
}

##' Load SIA data
##'
##' @param path Path in which to find the raw data
##' @import XLConnect
##' @author Sebastian Funk
loadSIAData <- function(path)
{

    try_require("XLConnect")

    sia.workbook <-
        loadWorkbook(paste(path, "/RVC-vaccination-info-from-JRF-WHO-HQ.xlsx",
                           sep = ""))
    sia <- data.table(readWorksheet(sia.workbook, "SIA 2000-2014"))
    setnames(sia, colnames(sia), tolower(colnames(sia)))

    sia <- clean_countries(sia)

    for (col in colnames(sia))
    {
        if ("character" %in% class(sia[, get(col)]))
        {
            sia[, paste(col) := factor(get(col))]
        } else if ("POSIXct" %in% class(sia[, get(col)]))
          {
              sia[, paste(col) := as.Date(get(col))]
          }
    }

    sia[, age.lower := sub("-.*$", "", age.group)]
    sia[, age.upper := sub("^.*-", "", age.group)]
    
    sia <- sia[grep("^<", age.lower), age.lower := "0"]
    sia <- sia[, age.lower := gsub("\\+", "", age.lower)]
    sia <- sia[grep(" M", age.lower), months := TRUE]
    sia <- sia[, age.lower := gsub(" M", "", age.lower)]
    sia <- sia[, age.lower := gsub(" Y", "", age.lower)]
    sia <- sia[, age.lower := as.integer(age.lower)]
    sia <- sia[months == TRUE, age.lower := as.integer(age.lower / 12)]

    sia <- sia[, age.upper := gsub("<", "", age.upper)]
    sia <- sia[grep("\\+$", age.upper), age.upper := "123"]
    sia <- sia[grep(" Y", age.upper), months := FALSE]
    sia <- sia[grep(" M", age.upper), months := TRUE]
    sia <- sia[, age.upper := gsub(" M", "", age.upper)]
    sia <- sia[, age.upper := gsub(" Y", "", age.upper)]
    sia <- sia[, age.upper := as.integer(age.upper)]
    sia <- sia[months == TRUE, age.upper := as.integer(age.upper / 12)]

    sia[, months := NULL]

    sia[is.na(age.lower) & is.na(age.upper) & age.group == "School-age",
        c("age.lower", "age.upper") := list(6, 18)]

    sia <- clean_countries(sia)

    save(sia, file = "sia.RData")
}

##' Load SIA data
##'
##' @param path Path in which to find the raw data
##' @import XLConnect reshape2
##' @author Sebastian Funk
loadRVCVaccinationData <- function(path)
{
    try_require("XLConnect")

    rvc.workbook <-
        loadWorkbook(paste(path, "/Historical-coverages-reported-to-RVC_28_10.xlsx",
                           sep = ""))
    rvc.measles <- data.table(readWorksheet(rvc.workbook, "Measles coverage to RVC"))

    cols1.nb <- grep("^Col[0-9]+", colnames(rvc.measles), invert = TRUE)
    cols1.names <- colnames(rvc.measles)[cols1.nb]
    if (cols1.nb[length(cols1.nb)] < ncol(rvc.measles))
    {
        cols1.nb <- c(cols1.nb, ncol(rvc.measles))
    }
    cols1 <- c()
    if (cols1.nb[1] > 1)
    {
        cols1 <- rep("", cols1.nb[1] -  1)
    }

    for (i in seq_len(length(cols1.nb) - 1))
    {
        cols1 <- c(cols1, rep(cols1.names[i], diff(cols1.nb)[i]))
    }
    cols2 <- which(!is.na(rvc.measles[1, ]))

    setnames(rvc.measles, cols2,
             paste(cols1[cols2],
                   unname(unlist(rvc.measles[1, cols2, with = FALSE])), sep = "."))

    
    setnames(rvc.measles, colnames(rvc.measles), tolower(colnames(rvc.measles)))
    setnames(rvc.measles, colnames(rvc.measles), gsub(" ", ".", colnames(rvc.measles)))
    rvc.measles <- rvc.measles[!is.na(country)]
    rvc.measles[, year := sub("\\/[0-9]+$", "", year)]
    rvc.measles[, year := gsub("(\\*|\\.)", "", year)]
    rvc.measles[, year := gsub("…", "", year)]
    rvc.measles[, year := gsub("г", "", year)]
    rvc.measles[, year := as.integer(year)]

    rvc.measles <- melt(rvc.measles, id.vars = c("country", "year"))

    rvc.measles[grep("1st", variable), dose := 1]
    rvc.measles[grep("2nd", variable), dose := 2]
    rvc.measles <- rvc.measles[!is.na(dose)]
    rvc.measles[, variable := sub("^measles.vaccine..*dose\\.", "", variable)]

    rvc.measles <- data.table(dcast(rvc.measles, country + year + dose ~ variable,
                                    value.var = "value"))

    rvc.mmr <- data.table(readWorksheet(rvc.workbook, "MMR to RVC"))

    

    save(rvc.measles, file = "rvc.RData")
}

##' Load worldwide serology data
##'
##' @param path Path in which to find the raw data
##' @import XLConnect reshape2
##' @author Sebastian Funk
loadWorldPopulationProspectData <- function(path)
{
    try_require("XLConnect")

    sheets <- getSheets(wpp.workbook)
    sheets <- setdiff(sheets, "NOTES")

    wpp.workbook <-
        loadWorkbook(paste0(path, "/WPP2012_POP_F07_1_POPULATION_BY_AGE_BOTH_SEXES.XLS"))

    wpp.table <- NULL

    for (sheet in sheets)
    {
        temp.table <- data.table(readWorksheet(wpp.workbook, sheet))
        temp.table <- temp.table[!is.na(Col1)]
        setnames(temp.table, colnames(temp.table), unlist(temp.table[1]))
        temp.table <- temp.table[-1]

        temp.table[, Index := NULL]
        if (is.null(wpp.table))
        {
            wpp.table <- temp.table
        } else
        {
            common.columns <- intersect(colnames(wpp.table),
                                        colnames(temp.table))
            wpp.table <- rbind(wpp.table[, common.columns, with = F],
                               temp.table[, common.columns, with = F])
        }
    }

    setnames(wpp.table, "Major area, region, country or area *", "area")
    setnames(wpp.table, "Reference date (as of 1 July)", "year")
    wpp.table[, Notes := NULL]
    for (col in 3:ncol(wpp.table))
    {
        colname <- colnames(wpp.table)[col]
        wpp.table[, paste(colname) := as.integer(get(colname))]
    }
    wpp <- melt(wpp.table, id.vars = c("Variant", "area", "Country code", "year"),
                variable.name = "age_group", value.name = "population")
    setnames(wpp, colnames(wpp), gsub(" ", "_", tolower(colnames(wpp))))
    
    save(wpp, file = "wpp.RData")
}

##' Load Childcare data
##'
##' @param path Path in which to find the raw data
##' @author Sebastian Funk
loadChickenpoxSerology <- function(path) {

    childcare_raw <-
        data.table(read.csv(paste(path, "/childcare.csv", sep = "")))
    years <- seq(ceiling(min(childcare$year)), floor(max(childcare$year)))

    approx_childcare <- approx(childcare$year, childcare$childcare, years)
    childcare <- data.table(year = approx_childcare$x,
                            childcare = approx_childcare$y)

    setkey(childcare, year)

    save(childcare, file = "childcare.RData")

}

