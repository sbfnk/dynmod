##' Convert a \code{data.frame} into a \code{ts} time series object
##'
##' @param data data to convert
##' @param interval data interval (day, week, month, quarter, year, or an arbitrary number)
##' @param date.column date column; if not given will be guessed as the first column of class \code{Date}; if not given and no column is of class \code{Date}, will quit with an error
##' @param value.column value column; if not given, will assume that dates are individual data per date (which will be added up)
##' @param start starting date
##' @param end end date
##' @import data.table
##' @export
##' @return A time series (\code{ts}) object
##' @author Sebastian Funk
df_to_ts <- function(data, interval = "day", date.column = NULL, value.column = NULL, start = NULL, end = NULL)
{

    data <- data.table(data)

    if (is.null(date.column))
    {
        date.column <- colnames(data)[min(which(sapply(colnames(dt), function(x) { class(dt[[x]])}) == "Date"))]
    }

    
}
