% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{kernel_gauss_dC}
\alias{kernel_gauss_dC}
\title{Derivative of Gaussian kernel covariance matrix in C}
\usage{
kernel_gauss_dC(x, theta, C_nonug, s2_est, beta_est, lenparams_D, s2_nug)
}
\arguments{
\item{x}{Matrix x}

\item{theta}{Theta vector}

\item{C_nonug}{cov mat without nugget}

\item{s2_est}{whether s2 is being estimated}

\item{beta_est}{Whether theta/beta is being estimated}

\item{lenparams_D}{Number of parameters the derivative is being calculated for}

\item{s2_nug}{s2 times the nug}
}
\value{
Correlation matrix
}
\description{
Derivative of Gaussian kernel covariance matrix in C
}
