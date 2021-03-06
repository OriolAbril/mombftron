% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_data.R
\name{get_grouped_corr_data}
\alias{get_grouped_corr_data}
\title{Generate data correlated by groups}
\usage{
get_grouped_corr_data(n, m_indep, m_corr, corr_mat, seed)
}
\arguments{
\item{n}{Number of observations in generated data}

\item{m_indep}{Number of independent covariates}

\item{m_corr}{Number of covariates correlated in group}

\item{corr_mat}{Correlation matix, has to be of dim m_corr. If missing
it uses a random correlation matrix.}

\item{seed}{Integer to use as random seed. If missing, no seed is set.}
}
\description{
Generates `m_indep * m_corr` covariates, with a grouped
correlation matrix with group size `m_corr`.
}
\examples{
get_grouped_corr_data(10, 3, 2, seed=7)
}
