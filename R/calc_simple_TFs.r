#' Calculation of simple MU - electrode transfer functions (wrapper)
#' 
#' Calculate transfer functions from each motor unit (MU) to each electrode in
#' a simple, distance-based, normalized, instantaneous fashion.
#' 
#' The 'transfer functions' from a single muscle fiber to an electrode are given
#' by
#'
#'    tf{MF, Electrode}(kt) = (dist_min / dist(MF, Electrode))^2 * j * kt,
#'
#' where dist_min is a normalization term which is defined as the minimum over
#' all MUs of
#'
#'    dist_MU = sqrt(sum(dist(MF, Electrode)^-2)),
#'
#' with the sum being taken over all muscle fibers in that MU. This formula
#' results from assuming a point-source solution
#'
#'    phi(xe, ye, ze) = (dist_min / dist(MF, Electrode))^2
#'
#' and actually considering point sources with time course
#'
#'    i(t) = d/dz psi(-vt).
#'
#' Fourier transformation then yields the above transfer function. The TF of a
#' MU now results from the summation of the TFs of all its muscle fibers,
#' yielding
#'
#'    tf{MU, Electrode}(kt) = dist_min^2 * j * kt * sum(dist(MF, Electrode)^-2).
#'
#' Note that this simple summation is only valid if we assume the same
#' conduction velocity (cv) in all muscle fibers located in the same MU.
#'
#' @param .muscles A data.frame containing data on the muscle geometry,
#'   especially a column 'MU.obj' containing MU objects.
#' @param .electrodes A data.frame containing columns 'electrode' with the
#'   electrode ID and 'electrode.obj' containing the actual object.
#' @param .volume.conductor A volume conductor object. In fact, this argument is
#'   not required for the simple TFs computed here, and is íncluded solely for
#'   the purpose of interface consistency.
#' @param .freqs A numerical vector containing the frequencies at which the
#'   transfer function is to be evaluated. Results will be returned in order.
#'
#' @return A data.frame containing the sampled transfer functions in column
#'   'TF', represented as vectors of complex values.
#'
calculate.simple.TFs <- function(.muscles, .electrodes, .vol.conductor,
                                 .freqs) {

    stopifnot(is.vector(.freqs),
              is.numeric(.freqs),
              all(Im(.freqs) == 0),
              is.data.frame(.muscles),
              is.data.frame(.electrodes))
    
    geom <- merge(.muscles, .electrodes)
    
    geom$TF.raw <- mapply(simple.TF,
                          geom$MU.obj,
                          geom$electrode.obj,
                          MoreArgs = list(2 * pi * .freqs),
                          SIMPLIFY = FALSE)

    amp.max <- max(sapply(geom$TF.raw, function(TF.raw) TF.raw$amp))

    geom$TF <- lapply(geom$TF.raw,
                      function(TF.raw) TF.raw$sampled.data / amp.max)

    geom[, c("muscle", "MU", "electrode", "TF")]
}


#' Calculation of simple MU - electrode transfer functions (inner wrapper)
#'
#' Compute discrete, unnormalized values of the transfer function from the motor
#' unit .MU to the electrode .electrode. For the transfer function model, see
#' comments above.
#'
#' @param .MU The motor unit which is the source of the surface potential.
#' @param .electrode The electrode at which the potential is measured.
#' @param .kt A vector containing the discrete angular frequencies at which the
#'   transfer function is to be evaluated.
#' 
#' @return A list object providing two properties: First, a property
#'   '$sampled.data' containing the discrete values of the transfer function
#'   from motor unit .MU to electrode .electrode obtained by evaluation at the
#'   discrete angular frequencies kt. Values are _not_ normalized according to
#'   .dist.min. And second, a property '$amp' that contains the sum over all
#'   muscle fibers of the squared inverse of the distance from that fiber's
#'   innvervation zone to the electrode position. This will be used for
#'   normalization later on.
#'
simple.TF <- function(.MU, .electrode, .kt) {

    TF <- list()

    TF$amp <- sum(
        sapply(.MU$fibers,
               function(fiber) dist.emg(fiber$innervation.zone,
                                        c(.electrode$position[1], 0,
                                          .electrode$position[2]))^-2))
    
    TF$sampled.data <- .kt * 1i * TF$amp

    TF
}
