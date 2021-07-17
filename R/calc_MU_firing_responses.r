
#' Simple worker function wrapper to facilitate progression display
calc_MU_firing_response_fftvals_wrapper = function(x, p, ...) {
  #p(message=sprintf("x=%g", x))
  setTxtProgressBar(p, x)
  #message(paste0("    * thread <",x,"> mem <",sprintf("%.1f MB",pryr::mem_used()/1e6),">"))
  calc_MU_firing_response_fftvals(...)
}


##' Calculate the values of the FFT of the firing responses of all MUs at all
##' electrodes
##'
##' This is just a wrapper that calls \link{calc_MU_firing_response_fftvals}
##' for all MU-electrode pairs.
##' 
##' @param MUs A data.frame containing a column 'MU.obj' containing MU objects
##'   and index columns 'MU' and 'muscle'.
##' @param TFs_MU_input_to_surf_potential A data.frame containing the values of
##'   the transfer functions at frequencies \code{sampling$freqs_to_calc} in
##'   column 'TF', represented as vectors of complex values.
##' @param psi_trafo A function returning the values of the unitary Fourier
##'   transform of psi, given a vector of ordinary input frequencies.
##' @param num_cores The number of cores that shall be used for computations. If
##'   equal to one (default), calculations are executed sequentially as usual.
##'   Otherwise, the parallel package is used for multi-process computations.
##' @param sampling A sampling object as created by \link{setup_sampling}.
##' @return A data.frame containing index columns 'muscle', 'MU' and 'electrode'
##'   as well as a data column 'firing_response_fftvals' that contains vectors
##'   that store the FFT values of the firing response of the particular MU at
##'   the particular electrode.
calc_MU_firing_responses_fftvals <- function(MUs,
                                             TFs_MU_input_to_surf_potential,
                                             psi_trafo,  num_cores = 1, sampling) {
    
    data <- merge(MUs, TFs_MU_input_to_surf_potential)
    
    message(paste0("    ",whoami(2)," : Running <",length(data$MU.obj),"> threads on <",num_cores,"> cores ..."))
    
    #p <- progressr::progressor(steps=length(data$MU.obj))
    p <- txtProgressBar(max=length(data$MU.obj),style=3)
    
    if (num_cores > 1) {
        
        cl <- parallel::makeCluster(num_cores, outfile="")
        future::plan(future::cluster, workers=cl, gc=TRUE)
        
        data$firing_response_fftvals <-
            with(data,
                 #parallel::clusterMap(cl,
                 future.apply::future_mapply(
                     calc_MU_firing_response_fftvals_wrapper,
                     seq_along(MU.obj),
                     MU_electrode_TF = TF,
                     MU_obj = MU.obj,
                     MoreArgs = list(
                         psi_trafo = psi_trafo,
                         sampling = sampling, 
                         p = p),
                     SIMPLIFY = FALSE,
                     future.seed = TRUE))
        
        parallel::stopCluster(cl)
        
    } else {
        
        data$firing_response_fftvals <-
            with(data,
                 mapply(
                     calc_MU_firing_response_fftvals_wrapper,
                     seq_along(MU.obj),
                     MU_electrode_TF = TF,
                     MU_obj = MU.obj,
                     MoreArgs = list(
                         psi_trafo = psi_trafo,
                         sampling = sampling,
                         p = p),
                     SIMPLIFY = FALSE))
    }
    
    close(p)
    data[, c("muscle", "MU", "electrode", "firing_response_fftvals")]
}



##' Calculate the values of the FFT of the firing response of a MU at an
##' electrode
##'
##' Take information on the transfer function MU-electrode and information on
##' the excitatory input signal fed into the MU and combine this to attain the
##' firing response of the given MU at the given electrode.
##'
##' This does not involve calculating an FFT! Actually, it is just
##' \itemize{
##'  \item{1.}{Evaluate the Fourier transform of the excitatory input signal
##'          (given in analytical form) at the same frequencies at which the
##'          MU-electrode TF has been computed.}
##'  \item{2.}{Multiply the complex conjugates of these values with the
##'            corresponding TF values.}
##'  \item{3.}{Obtain the second half of the FFT by mirroring the result,
##'            exploiting the hermitianness of the FFT of a real signal (which
##'            we expect the resulting firing response to be).}
##' }
##'
##' @param MU_electrode_TF A complex vector containing the value of the transfer
##'   function from MU to electrode at frequencies \code{sampling$freqs_to_calc}
##' @param MU_obj A MU object.
##' @param psi_trafo A function returning the values of the unitary Fourier
##'   transform of psi, given a vector of ordinary input frequencies.
##' @param sampling A sampling object as created by \link{setup_sampling}.
##' @return A complex vector that stores the FFT values of the firing response
##'   of MU at electrode.
calc_MU_firing_response_fftvals <- function(MU_electrode_TF, MU_obj, psi_trafo,
                                            sampling) {

    freqs <- sampling$freqs_to_calc / MU_obj$cv

    N <- length(freqs)

    psi_vals <- psi_trafo(freqs) * MU_obj$IAP_scale_amp

    stopifnot(length(MU_electrode_TF) == length(psi_vals))

    response_fftvals_half <- MU_electrode_TF * Conj(psi_vals) / MU_obj$cv

    response_fftvals <- mirror_hermitian_fft(response_fftvals_half, sampling)
    
    response_fftvals
}    



## -- Data restructuring --
## Currently, firing_responses_fftvals has columns `muscle`, `MU`,
## `electrode` and `firing_response_fftvals`, where the last column
## contains a list of the fftvals in each row of the top-level dataframe.
## What we want is the same thing, but with the fftvals for all electrodes
## coerced into a single list item, i.e.
##   $muscle num
##   $MU num
##   $firing_response_fftvals list of num_electrodes elements, each a vector
##      of the corresponding fftvals.
## This structure is required so as to allow passing a single row of the
## top-level df to subroutines simulating a single motor unit. This is
## helpful, since firing event handling and force simulation needs only be
## done once for each MU; hence electrode-wise handling appears suboptimal.
restructure_firing_responses_fftvals <- function(MU_firing_responses_fftvals) {
  fftvals_dfs <- list()
  i <- 1
  for (muscle_id in unique(MU_firing_responses_fftvals$muscle)) {
    muscle_reduced_df <-
      MU_firing_responses_fftvals %>% dplyr::filter(muscle==muscle_id)
    for (MU_id in unique(muscle_reduced_df$MU)) {
      MU_reduced_df <- muscle_reduced_df %>% dplyr::filter(MU==MU_id)
      fftvals <- list()
      fftvals_dfs[[i]] <- data.frame(muscle=muscle_id, MU=MU_id)
      for (electrode_id in unique(MU_reduced_df$electrode)) {
		electrode_reduced_df <- MU_reduced_df %>% dplyr::filter(electrode == electrode_id)
        fftvals[[electrode_id]] <-
        electrode_reduced_df$firing_response_fftvals[[1]]
	  }
      fftvals_dfs[[i]]$firing_response_fftvals <- list(fftvals)
      i <- i+1
    }
  }
  
  plyr::rbind.fill(fftvals_dfs)
}
