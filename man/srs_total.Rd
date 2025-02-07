% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/srs.R
\name{srs_total}
\alias{srs_total}
\title{Estimate Population Total}
\usage{
srs_total(sample_data, N, alpha = 0.05)
}
\arguments{
\item{sample_data}{a vector, the sample data.}

\item{N}{a number, the size of population.}

\item{alpha}{a number, the confidence level.}
}
\value{
A list containing:\tabular{ll}{
\code{ytot_result} \tab a list containing estimation of the population total, estimation of the variance of the estimator and the estimation of the standard error of the estimator. \cr
\tab \cr
\code{ci_est} \tab a list contain absolute/relative error limit and the estimation of lower/upper confidence bound of the estimator.\cr
\tab \cr
}
}
\description{
This function is used to estimate the population total in simple random sampling.
}
