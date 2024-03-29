% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Phi_prime_fun.R
\name{Phi_prime_fun}
\alias{Phi_prime_fun}
\title{Computes the Deflated Relative Occurrences}
\usage{
Phi_prime_fun(
  p,
  T_stop,
  num_dummies,
  phi_T_mat,
  Phi,
  eps = .Machine$double.eps
)
}
\arguments{
\item{p}{Number of candidate variables.}

\item{T_stop}{Number of included dummies after which the random experiments (i.e., forward selection processes) are stopped.}

\item{num_dummies}{Number of dummies}

\item{phi_T_mat}{Matrix of relative occurrences for all variables (i.e., j = 1,..., p) and for T = 1, ..., T_stop.}

\item{Phi}{Vector of relative occurrences for all variables (i.e., j = 1,..., p) at T = T_stop.}

\item{eps}{Numerical zero.}
}
\value{
Vector of deflated relative occurrences for all variables (i.e., j = 1,..., p) and T = T_stop.
}
\description{
Computes the vector of deflated relative occurrences for all variables (i.e., j = 1,..., p) and T = T_stop.
}
