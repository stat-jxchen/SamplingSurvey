% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unequal.R
\name{PPS_HH_total}
\alias{PPS_HH_total}
\title{Hansen-Hurwitz Estimator for Population Total}
\usage{
PPS_HH_total(y, z, alpha = 0.05)
}
\arguments{
\item{y}{a vector, the sample data.}

\item{z}{a vector, the probability that the corresponding sample is selected
in a single sampling.}

\item{alpha}{a number, the confidence level.}
}
\description{
This function estimates the population total with Hansen-Hurwitz estimator
when PPS sampling is used.
}
