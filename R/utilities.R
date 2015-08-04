.datatable.aware <- TRUE

##' Reduce the number of age groups given a broader set of limtis
##'
##' Operates on lower limits
##' @param x data table in which to reduce age groups
##' @param limits new limits
##' @param column column in x that denotes age
##' @return data table with updated age groups
##' @export
reduce.agegroups <- function(x, limits, column = "lower.age.limit") {

    x <- x[, new.agegroup := limits[findInterval(x[, get(column)], limits)]]
    x <- x[, paste(column) := NULL]
    setnames(x, "new.agegroup", column)

    x
}

##' Convert lower age limits to age groups.
##'
##' Mostly used for plot labelling
##' @param limits Lower age limits
##' @return Age groups (limits separated by dashes)
##' @export
limits.to.agegroups <- function(limits) {
    unique.limits <- unique(limits)
    if (length(unique.limits) > 1)
    {
        agegroups <- c(sapply(seq(1, length(unique.limits) - 1), function(x) {
            if ((unique.limits[x+1] - 1) > unique.limits[x]) {
                paste(unique.limits[x], unique.limits[x+1] - 1, sep = "-")
            } else {
                unique.limits[x]
            }
        }), paste(unique.limits[length(unique.limits)], "+", sep = ""))
    } else
    {
        agegroups <- c("all")
    }
    return(agegroups)
}

##' Convert age groups to lower age limits.
##'
##' Mostly used to go from a contact matrix to lower limits
##' @param groups age groups
##' @return lower age limits
agegroups.to.limits <- function(groups) {
    lower.limits <- sub("^\\[", "", groups)
    lower.limits <- sub(",[0-9]+\\)$", "", lower.limits)
    lower.limits <- as.integer(lower.limits)
    return(lower.limits)
}

##' Convert year and week to date.
##'
##' Mostly used to create a date field for ggplot
##' @param year Year(s)
##' @param week Week(s)
##' @return Date(s)
##' @import ISOweek
year.week.to.date <- function(year, week) {

    dates <- data.table(year = year, week = week)

    dates <-
        dates[, ISOweek := paste(year, "-W",
                           sprintf("%02i", week), "-1", sep="")]
    dates <- dates[, date := ISOweek2date(ISOweek)]
    return(dates[, date])
}

##' Convert lower to upper age limits
##'
##' @param x the age group to convert
##' @param limits overall lower age limits
##' @param max maximum age
##' @return upper age limits
##' @author Sebastian Funk
lower.to.upper.limits <- function(x, limits, max = 100) {
    upper.limits <- c(limits[-1], max)
    return(upper.limits[match(x, limits)])
}

##' Combine compartments that have been split up to get an Erlang distribution
##'
##' @param traj trajectory data frame
##' @return trajectory data frame with Erlang compartments summed up
##' @author Sebastian Funk
##' @import reshape2
combine.compartments <- function(traj) {

    mtraj <- melt(traj, "time")
    mtraj <- mtraj[, variable := sub("\\..+$", "", variable)]
    mtraj <- mtraj[, list(value = sum(value)), by = list(time, variable)]
    traj <- data.table(dcast(mtraj, time ~ variable))

    return(traj)
}

##' Melt a trajectory in terms of time and age groups
##'
##' @param traj Trajectory
##' @param age.labels labels of the age groups (optional)
##' @return a nice, melted trajectory
##' @author Sebastian Funk
##' @import reshape2
melt.trajectory <- function(traj, age.labels = NULL, time.label = NULL) {

    traj <- combine.compartments(traj)

    for (Z.group in grep("Z", names(traj), value = T)) {
        traj <-
            traj[, gsub("Z", "abs.incidence", Z.group) :=
                     c(round(diff(get(Z.group))), NA)]
    }

    mtraj <- melt(traj, "time")
    mtraj <- mtraj[, agegroup := as.integer(sub("[^0-9]+", "", variable))]
    mtraj <- mtraj[!is.na(agegroup) & !(agegroup == 0)]
    mtraj <- mtraj[, variable := sub("[0-9]", "", variable)]
    traj <- data.table(dcast(mtraj, time + agegroup ~ variable))

    if (is.null(time.label)) {
        time.key <- "time"
    } else {
        setnames(traj, "time", time.label)
        time.key  <- time.label
    }

    if (is.null(age.labels)) {
        age.key <- "agegroup"
    } else {
        traj <-
            traj[, agegroup := as.integer(as.character(factor(traj$agegroup,
                            labels = age.labels)))]
        setnames(traj, "agegroup", "lower.age.limit")
        age.key <- "lower.age.limit"
    }

    setkeyv(traj, c(time.key, age.key))
    if (length(unique(traj[, get(age.key)])) == 1)
    {
        traj[, paste(age.key) := NULL]
    }

    return(traj)
}

##' Format a number padded with zeroes
##'
##' @param no the number to print
##' @param max the maximal number (which determines the numbef of zeroes
##' @return padded character
##' @author Sebastian Funk
##' @export
pad.zeroes <- function(no, max)
{
    len <- nchar(as.character(max))
    return(sprintf(paste("%0", len, "d", sep = ""), no))
}
