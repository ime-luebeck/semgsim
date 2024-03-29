% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IAP_Rosenfalck.r
\name{psi_Rosenfalck_transformed}
\alias{psi_Rosenfalck_transformed}
\title{Implementation of the unitary Fourier transform of \code{psi.Rosenfalck}}
\usage{
psi_Rosenfalck_transformed(fz, A = 96 * 0.001, B = -90 * 0.001)
}
\arguments{
\item{fz}{Ordinary frequency [Hz], scalar or vector.}
}
\value{
Scalar or vector, depending on the input size.
}
\description{
Implements a Fourier Transform of Rosenfalck's psi function, scaled to take
its input variable in SI unit m, instead of in mm.
}
