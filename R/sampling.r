#' Set up tranfer function sampling
#'
#' Determine the frequencies at which the transfer functions are to be sampled.
#'
#' Uses information on the motor units to obtain an estimate for the maximum
#' length of a MU firing response. This in turn defines - for given sampling
#' frequency - the number of samples that needs to be computed for each firing
#' response. And thus the number of points to be computed for the inverse FFT is
#' defined.
#'
#' @note This function will always choose an odd NFFT since even NFFTs lead to
#'   spectra that - when inverse transformed - yield non-zero imaginary parts in
#'   the resulting time signal. This is due to the presence of the single
#'   frequency coefficient at Fs/2. This coefficient is required to be purely
#'   real for purely real time domain signals. This can not be guaranteed,
#'   however, when the coefficient stems from the evaluation of analytical
#'   expressions for the fourier transform, as it is the case here. By choosing
#'   NFFT to be odd, we circumvent this problem as then the frequency bin at
#'   Fs/2 is simply not present.
#'
#' @param muscles A data.frame containing muscles as typically created by
#'   \code{\link{convert_muscle_list_to_df}}.
#' @param f_target The temporal sampling frequency in Hz that is desired for the
#'   final output.
#' @param psi_length The length of the relevant (non-zero) section of the
#'   function psi, that represents the shape of the IAP. This is relevant since
#'   it influences the length of the firing response of a MU.
#' @param max_freq_dom_sample_dist If given (unit is Hz), specifies an upper
#'   bound on the allowed distance between two sampled points in the frequency
#'   domain. This is achieved by adjusting NFFT accordingly.
#' @return A "sampling" object with members "Fs", "NFFT", "fft_freqs" and
#'   "freqs_to_calc". All frequencies are in Hz.
#' 
setup_sampling <- function(muscles, f_target, psi_length,
                           max_freq_dom_sample_dist = Inf,
						   unsampled_firings = TRUE) {

    min_cond_vel <- muscles$muscle.obj %>%
        lapply(function(muscle.obj) muscle.obj$MU_size_params$cv$lower_bound) %>%
        simplify2array %>%
        min
    
    max_fiber_half_length <- muscles$muscle.obj %>%
        plyr::laply(function(muscle) muscle$fiber_length_max *
                    muscle$fiber_rel_half_length_max) %>%
                        max

    max_response_length <- (max_fiber_half_length + psi_length) / min_cond_vel

    ## Number of samples necessary to represent the longest possible MU firing
    ## response. Add one extra sample just to be sure.
    max_response_length_samples <- ceiling(max_response_length * f_target) + 1

    ## Number of samples necessary to comply with the user-defined maximum
    ## inter-point distance in the frequency domain
    min_num_samples <- ceiling(f_target / max_freq_dom_sample_dist)

    ## Note how we only allow odd NFFTs
    num_samples <- nextn(max(max_response_length_samples,
                                       min_num_samples),
                                   factors = c(3, 5))

    sampling <- create_sampling(Fs = f_target, NFFT = num_samples)
	sampling$unsampled_firings <- unsampled_firings
    sampling
}


create_sampling <- function(Fs, NFFT) {

    sampling <- list()
    sampling$Fs <- Fs
    sampling$NFFT <- NFFT

    sampling$fft_freqs <- calc_fft_freqs(sampling$NFFT, Fs = Fs)
    
    if (is.even(sampling$NFFT)) {

        sampling$freqs_to_calc <-
            sampling$fft_freqs[1:(round(sampling$NFFT/2) + 1)]
        
    } else {
        
        sampling$freqs_to_calc <-
            sampling$fft_freqs[1:ceiling(sampling$NFFT/2)]
    }

    sampling
}
