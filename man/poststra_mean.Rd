% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subpop_and_poststra.R
\name{poststra_mean}
\alias{poststra_mean}
\title{Estimate Population Mean with Poststratification}
\usage{
poststra_mean(Wh, nh, ybarh, s2h, N = Inf, alpha = 0.05)
}
\arguments{
\item{Wh}{a vector, the weight of each strata.}

\item{nh}{a vector, the sample size of each strata.}

\item{ybarh}{a vector, the sample mean of each strata.}

\item{s2h}{a vector, the sample variance of each strata.}

\item{N}{a number, the population size.}

\item{alpha}{a number, the confidence level.}
}
\description{
This function is used to estimate the population mean in simple random sampling with poststratification method.
}
