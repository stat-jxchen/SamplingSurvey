% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stra.R
\name{stra_sample_extract}
\alias{stra_sample_extract}
\title{Sample Information Extraction}
\usage{
stra_sample_extract(sample_data, type = "mean")
}
\arguments{
\item{sample_data}{a list, each element denotes a strata sample.}

\item{type}{string, "mean" for \code{stra_srs_mean}, "prop" for \code{stra_srs_prop}.}
}
\description{
This function is an auxiliary function of \code{stra_srs_mean} and
\code{stra_srs_prop}, which is used to extract some information from the sample.
}
