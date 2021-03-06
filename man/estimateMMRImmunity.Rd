% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/susceptibility.R
\name{estimateMMRImmunity}
\alias{estimateMMRImmunity}
\title{Estimate immunity from MMR and case data}
\usage{
estimateMMRImmunity(countries = NULL, vaccine.efficacy = 0.95,
  estimation.year = 2014, age.dist.cases = NULL, reporting = 1,
  average.last = 0, coverage.estimate = FALSE, use.sia = FALSE)
}
\arguments{
\item{countries}{the countries to investigate. If set to NULL
(default), will estimate immunity for all countries}

\item{vaccine.efficacy}{vaccine efficacy}

\item{estimation.year}{year for which to estimate the profile}

\item{age.dist.cases}{age distribution of cases}

\item{reporting}{proportion of cases that get reported}

\item{average.last}{average over the last X years to get further
years of cases (0 if no averaging is to be done)}
}
\value{
A data table with immunity by age in 2014
}
\description{
Estimate immunity from MMR and case data
}
\author{
Sebastian Funk
}

