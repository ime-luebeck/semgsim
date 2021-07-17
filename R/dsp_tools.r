#' sinc function of frequency f.
sinc <- function(x, f)
    if (x == 0) {
        2 * pi * f
    } else {
        sin(2 * pi * f * x) / x
    }

#' Blackman window from 0 to m. (Center is at m/2+1.)
blackman <- function(m)
    0.42 - 0.5 * cos(2 * pi * (0:m) / m) + 0.08 * cos(4 * pi * (0:m) / m)

#' Hamming window from 0 to m. (Center is at m/2+1.)
hamming <- function(m)
    0.54 - 0.46 * cos(2 * pi * (0:m) / m)

#' Boxcar window from 0 to m. (Boxcar is equal to one everywhere.)
boxcar <- function(m)
    rep(1, m + 1)


#' Fast Fourier Transform (FFT) of a function
#'
#' Provides a utility wrapper for the builtin FFT function that is applied
#' directly to function objects and that allows for the specification of desired
#' sampling intervals, number of sampling points and spatial and amplitude
#' scaling.
#'
#' @param func The function to be fft-ed.
#' @param NFFT The number of sample points to be taken, thus the number of
#'   points in the fft.
#' @param range The (spatial/temporal) range from which the samples are to be
#'   selected. Expected to be a two-element numerical vector. First element
#'   should be the smaller of the two boundaries. Note that the lower boundary
#'   will not be sampled since lower_bound + stepsize is taken as the starting
#'   point.
#' @param scale_spatial (Optional) Can be given to take the fft of
#'   func(scale_spatial * x) instead of func(x). Yields exactly the same result
#'   as 'calc_fft(function(x) func(2*x), NFFT, range)'.
#' @param scale_amp (Optional) Scale the amplitude of func by taking the fft of
#'   scale_amp * func(x) instead of func(x).
#'
#' @return An object of class "fft" that has two properties: $fft, which
#'   contains the actual result of the fft, and $vars, which contains a variety
#'   of parameters involved with the fft, such as the sampled frequencies, etc.
#'
calc_fft <- function(func, NFFT, range, scale_spatial = 1.0, scale_amp = 1.0) {

    ## Sanity checks
    stopifnot(is.function(func),
              is.nice.scalar(NFFT),
              is.vector(range),
              is.numeric(range),
              length(range) == 2,
              is.nice.scalar(scale_spatial),
              is.nice.scalar(scale_amp))

    samples <- scale_spatial * calc_fft_samples(NFFT, range)

    ft <- calc_fft_freqs(NFFT, range)

    kt <- 2 * pi * ft

    vars <- environment()

    fft_obj <- list()
    class(fft_obj) <- "fft"

    fft_obj$fft <- fft(func(samples)) / NFFT  * scale_amp

    fft_obj$vars <- vars

    fft_obj
}


calc_fft_freqs <- function(NFFT, range = NULL, Fs = NULL) {

    stopifnot(is.not.null(range) || is.not.null(Fs))

    if (is.not.null(range)) {

        step <- (range[2] - range[1]) / NFFT
        fs <- 1 / step

    } else {
        
        fs <- Fs
    }
    
    if (is.even(NFFT)) {
        ft <- c(seq(0, fs/2, fs/NFFT),
                seq(-fs/2 + fs/NFFT, -fs/NFFT, fs/NFFT))
    } else {
        ft <- c(seq(0, (NFFT-1)*fs/2/NFFT, fs/NFFT),
                seq(-(NFFT-1)*fs/2/NFFT, -fs/NFFT, fs/NFFT))
    }

    stopifnot(length(ft) == NFFT)
    
    ft
}


calc_fft_samples <- function(NFFT, range) {

   step <- (range[2] - range[1]) / NFFT
   samples <- seq(range[1] + step, range[2], step)
   samples
}


append_fft <- function(obj) {

    NFFT <- obj$nSamples
    
    ## Calculate the sample time
    obj$fft.timestep <- obj$max.position / NFFT

    ## Calculate the sampling frequency
    obj$fft.fs <- 1 / obj$fft.timestep
    Fs <- obj$fft.fs

    ## Calculate the sampling instants
    obj$fft.samples <- seq(0 + obj$fft.timestep,
                            obj$max.position, obj$fft.timestep)

    ## Calculate the fft at these instants
    obj$fft <- fft(obj$wave.shape(obj$fft.samples)) / NFFT

    ## Derive the sampling points in the frequency domain
    if (is.even(NFFT)) {
        obj$fft.ft <- c(seq(0, Fs/2, Fs/NFFT),
                         seq(-Fs/2 + Fs/NFFT, -Fs/NFFT, Fs/NFFT))
    } else {
        obj$fft.ft <- c(seq(0, (NFFT-1)*Fs/2/NFFT, Fs/NFFT),
                         seq(-(NFFT-1)*Fs/2/NFFT, -Fs/NFFT, Fs/NFFT))
    }
              
    ## Convert to angular frequency
    obj$fft.kt <- 2 * pi * obj$fft.ft

    ## Return the original obj with all these new properties attached.
    obj
}
    
#' Perform a shift in time domain on fft values
#'
#' Takes already calculated values of an fft and applies the relation between a
#' shift in time domain and a multiplication by a complex factor in the
#' frequency domain to obtain a shifted time signal.
#'
#' @note This function will work but most likely produce undesired results when
#'   applied to signals of even length that have a non-zero coefficient at Fs/2.
#'   This coefficient is required to have a purely real coefficient for the
#'   inverse transformed signal to be purely real. However, the time-shift
#'   operation will possibly multiply this coefficient by a complex factor, thus
#'   resulting in a complex (as opposed to purely real) inverse FFT.
#' 
#' @param fft_vals Values of the fourier transform of a function calculated at
#'   the frequencies given by sampling$fft_freqs.
#' @param sampling An object providing information on the sampling paramters
#'   chosen. Usually constructed by the function create.sampling.
#' @param time_shift The desired time shift in SI units, i.e. seconds.
#' 
#' @return FFT values of the same frequencies as in fft_vals, but s.t. they
#'   produce a time-shifted signal when inverse transformed.
#' 
#' @examples
#' t <- seq(0.1, 1, 0.1)
#' y <- sin(2*pi*t)
#' fft_shifted <- semg.simulatr:::time_shift_fft(
#'     fft_vals = fft(y) / length(y),
#'     sampling = semg.simulatr:::create_sampling(Fs = 10, NFFT = 10),
#'     time_shift = 0.2)
#' dev.new()
#' plot(t,y)
#' points(t,Re(fft(fft_shifted, inverse = TRUE)), col = "red")
#'
time_shift_fft <- function(fft_vals, sampling, time_shift) {

    stopifnot(is.nice.scalar(time_shift))
    
    rotation <- exp(-2 * pi * 1i * time_shift * sampling$fft_freqs)
    rel_err <- calc_abs_rel_err(1, abs(rotation))
    stopifnot(all(rel_err < 1e-14))
    
    fft_vals_shifted <- fft_vals * rotation
    fft_vals_shifted
}

#' Mirror a set of FFT values by assuming the fourier transform to be hermitian.
#'
#' The FT of a purely real function is hermitian. Exploit this to obtain a full
#' set of FFT values from the values at positive frequencies only (and f=0).
#'
#' @note Should be clear, but once again: This function produces false results
#'   but **does not throw an error** if the underlying function is not purely
#'   real!
#' 
#' @param fft_vals A complex vector containing the FFT values of a function for
#'   the frequencies given by sampling$freqs_to_calc (i.e. f=0 as well as the
#'   positive ones).
#' @param sampling A sampling object, as usually created by
#'   \code{create.sampling}.
#' @return A complex vector containing the full set of FFT values at all
#'   frequencies (and in the same order as) given by sampling$fft_freqs.
#' 
mirror_hermitian_fft <- function(fft_vals, sampling) {

    ## Exploit hermitian symmetry due to realness of the resulting signal.
    N <- length(fft_vals)
    
    if (is.even(sampling$NFFT))
        fft_vals_mirrored <- c(fft_vals, Conj(rev(fft_vals[-c(1,N)])))
    else
        fft_vals_mirrored <- c(fft_vals, Conj(rev(fft_vals[-1])))

    stopifnot(length(fft_vals_mirrored) == sampling$NFFT)

    fft_vals_mirrored
}



is_fft_of_real_signal <- function(fft_vals, sampling, abstol = 1e-14) {

    pos_freqs <- which(sampling$fft_freqs > 0)

    if (is.even(sampling$NFFT)) {

        if (abs(Im(fft_vals[1])) > abstol ||
            abs(Im(fft_vals[sampling$NFFT/2 + 1])) > abstol)
            return(FALSE)
        
        pos_freqs <- which(sampling$fft_freqs > 0)
        pos_freqs <- pos_freqs[1:(length(pos_freqs)-1)]
        neg_freqs <- which(sampling$fft_freqs < 0)

    } else {

        if (abs(Im(fft_vals[1])) > abstol)
            return(FALSE)
        
        pos_freqs <- which(sampling$fft_freqs > 0)
        neg_freqs <- which(sampling$fft_freqs < 0)
    }

    stopifnot(length(pos_freqs) == length(neg_freqs))

    abs_err <- abs(fft_vals[pos_freqs] - Conj(rev(fft_vals[neg_freqs])))
    
    all(abs_err < abstol)
}


#' univariate Akima interpolation
#'
#' Creates an interpolating function of irregular, univariate data.
#'
#' Employs Akima's algorithm for univariate interpolation. This is a Hermite
#' interpolation method, prescribing slopes at data points depending on the
#' neighbouring points. It is a local method.
#' 
#' Adapted from the "akimaInterp" function in the pracma package, written by Hans W. Borchers.
#' (Which is also licensed under GPLv3.)
#'
#' @param x Independent variable data
#' @param y Dependent variable data
#' @references Akima, Hiroshi, "A New Method of Interpolation and Smooth Curve
#' Fitting Based on Local Procedures", 1970.
#' @return An interpolating function, in the spirit of \link{stats::splinefun}.
#' @export
akimafun <- function(x,y) {
    stopifnot(is.numeric(x), is.numeric(y), is.numeric(xi), is.vector(x), 
        is.vector(y), is.vector(xi))
    n <- length(x)
    if (length(y) != n) 
        stop("Vectors 'x' and 'y' must be of the same length.")
    dx <- diff(x)
    if (any(dx <= 0)) 
        stop("Argument 'x' must be an in strictly ascending order.")
    if (any(xi < x[1]) || any(xi > x[n])) 
        stop("All points in 'xi' must lie between x[1] and x[n].")
    m <- diff(y)/dx
    mm <- 2 * m[1] - m[2]
    mmm <- 2 * mm - m[1]
    mp <- 2 * m[n - 1] - m[n - 2]
    mpp <- 2 * mp - m[n - 1]
    m1 <- c(mmm, mm, m, mp, mpp)
    dm <- abs(diff(m1))
    f1 <- dm[3:(n + 2)]
    f2 <- dm[1:n]
    f12 <- f1 + f2
    id <- which(f12 > 1e-08 * max(f12))
    b <- m1[2:(n + 1)]
    b[id] <- (f1[id] * m1[id + 1] + f2[id] * m1[id + 2])/f12[id]
    e <- (3 * m - 2 * b[1:(n - 1)] - b[2:n])/dx
    d <- (b[1:(n - 1)] + b[2:n] - 2 * m)/dx^2

    function(xi) {
        bin <- findInterval(xi, x)
        bin <- pmin(bin, n - 1)
        bb <- bin[1:length(xi)]
        wj <- xi - x[bb]
        yi <- ((wj * d[bb] + e[bb]) * wj + b[bb]) * wj + y[bb]
        return(yi)
    }
}
