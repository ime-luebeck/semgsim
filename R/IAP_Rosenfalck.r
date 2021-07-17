#' Rosenfalck's IAP model function
#'
#' An implementation of a variant of Rosenfalck's model of the Intracellular
#' Action Potential (IAP) wave shape.  The model is simply
#'
#'    Vm(z[mm]) = Az^3 exp(-z) + B, if z > 0 and
#'    Vm(z[mm]) = B, if z <= 0.
#'
#' Note that the front of the IAP is towards z = 0, whereas the tail of the wave
#' is around z = 15. The wave is running in negative z direction.
#'
#' @param z The spatial variable. Can be a scalar or a vector. Should be in
#'   \code{mm} when employing the default parameters for \code{A} and \code{B}.
#' @param A Model parameter. Must be a scalar value in [V/mm^3]. Default value
#'   is 96 mv/mm^3.
#' @param B Model parameter. Must be a scalar value in SI units. Default value
#'   is -90mV.
#' @return A numerical vector containing the value(s) of Rosenfalcks function at
#'   position(s) z.
#' @export
IAP_Rosenfalck <- function(z, A = 96 * 1e-3, B = -90 * 1e-3) {
    
    ## Determine positive z values
    z_pos_mask <- z >= 0
    z_pos <- z[z_pos_mask]

    ## Generally, calling this with negative z arguments will probably not be
    ## intended behavior. Give out a warning.
    if (length(z_pos) < length(z))
        warning("Rosenfalck's function called with negative z argument.")

    ## This is the actual implementation of the formula.
    res <- vector("numeric", length(z))
    res[z_pos_mask] <- A * z_pos^3 * exp(-z_pos) + B
    res[!z_pos_mask] <- B

    ## Return the result.
    res
}


#' An implementation of Rosenfalck's model of the first derivative of the
#' Intracellular Action Potential (IAP) wave shape, evaluated in the negative
#' z direction. This quantity is called psi.
#'
#'    psi(z[mm]) = (Vm(-z))'
#'           = | -(z+3) A z^2 exp(z)                                if z < 0
#'             | 0                                                  if z > 0
#'             | -(z+3) A z^2 exp(z) + A z^3 exp(z) delta(z), if z = 0.
#'
#' Here, we implement psi(0) = 0.
#' Note that the front of the IAP is towards z = 0, whereas the tail of the wave
#' is around z = -15. It is running in positive z direction.
#'
#' @param z The spatial variable. Can be a scalar or a vector. Should be in
#'   \code{mm} when employing the default parameters for \code{A} and \code{B}.
#' @param A Model parameter. Must be a scalar value in [V/mm^3]. Default value
#'   is 96 mv/mm^3.
#' @param B Model parameter. Must be a scalar value in SI units. Default value
#'   is -90mV.
#' @return A numerical vector containing the value(s) of Rosenfalcks psi
#'   function at position(s) z.
#' @export
psi_Rosenfalck <- function(z, A = 96 * 1e-3, B = -90 * 1e-3) {

    z_pos_mask <- z >= 0
    z_neg <- z[!z_pos_mask]

    ## This is the actual implementation of the formula.
    ## Note that at z = 0 (where a Dirac delta would occur), 0 is returned.
    res <- vector("numeric", length(z))
    res[z_pos_mask] <- 0
    res[!z_pos_mask] <- -(z_neg + 3) * A * z_neg^2 * exp(z_neg)

    ## Return the result.
    res
}

##' Implementation of the unitary Fourier transform of \code{psi.Rosenfalck}
##'
##' Implements a Fourier Transform of Rosenfalck's psi function, scaled to take
##' its input variable in SI unit m, instead of in mm.
##'
##' @param fz Ordinary frequency [Hz], scalar or vector.
##' @return Scalar or vector, depending on the input size.
##' @export
psi_Rosenfalck_transformed <- function(fz, A = 96 * 1e-3, B = -90 * 1e-3)
    750i * 1e6 * A * pi * fz / (500i + pi * fz)^4
