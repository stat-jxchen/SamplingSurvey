% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/two_stage.R
\name{two_stage_eqM_mean}
\alias{two_stage_eqM_mean}
\title{Estimate Population Mean Under Two-Stage Sampling}
\usage{
two_stage_eqM_mean(M, m, ybari, s2i, N = Inf, alpha = 0.05)
}
\arguments{
\item{M}{a number, the number of secondary units in each primary unit in population.}

\item{m}{a number, the number of secondary units in each primary unit in sample.}

\item{ybari}{a vector, the sample mean of each primary unit in sample.}

\item{s2i}{a vector, the sample variance of each primary unit in sample.}

\item{N}{a number, the number of primary units in population.}

\item{alpha}{a number, the confidence level.}
}
\description{
This function estimates the population mean when two-stage sampling is used and
the cluster sizes are equal.
}
