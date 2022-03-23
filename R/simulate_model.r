
#' Simple worker function wrapper to facilitate progression display
simulate_MU_wrapper = function(x, p, ...) {
  #p(message=sprintf("x=%g", x))
  setTxtProgressBar(p, x)
  #message(paste0("    * thread <",x,"> mem <",sprintf("%.1f MB",pryr::mem_used()/1e6),">"))
  simulate_MU(...)
}


simulate_MUs = function(data, sampling, num_cores) {
  
  message(paste0("    ",whoami(2)," : Running <",length(data$MU.obj),"> threads on <",num_cores,"> cores ..."))
  
  #p <- progressr::progressor(steps=length(data$MU.obj))
  p <- txtProgressBar(max=length(data$MU.obj),style=3)

  if (num_cores > 1) {
    cl <- parallel::makeCluster(num_cores, outfile="") 
    future::plan(future::cluster, workers=cl, gc=TRUE)
    
    #res <- parallel::clusterMap(cl, 
    res <- future.apply::future_mapply(
      simulate_MU_wrapper, 
      x = seq_along(data$MU.obj),
      MU_obj = data$MU.obj,
      sim_state = data$sim_state,
      MU_contribs = data$contribs,
      MU_impulse_train = data$impulse_train,
      force = data$force,
      excitation = data$excitation,
      MU_rate_function = data$MU_rate_function,
      firing_response_fftvals = data$firing_response_fftvals,
      MoreArgs = list(sampling = sampling, 
                      p = p),
      SIMPLIFY = FALSE,
      future.seed = TRUE)
    
    parallel::stopCluster(cl)
    
  } else {
    
    res <- mapply(
      simulate_MU_wrapper,
      seq_along(data$MU.obj),
      MU_obj = data$MU.obj,
      sim_state = data$sim_state,
      MU_contribs = data$contribs,
      MU_impulse_train = data$impulse_train,
      force = data$force,
      excitation = data$excitation,
      MU_rate_function = data$MU_rate_function,
      firing_response_fftvals = data$firing_response_fftvals,
      MoreArgs = list(sampling = sampling, 
                      p = p),
      SIMPLIFY = FALSE)
  } 
  
  close(p)
  res
}


simulate_MU <- function(MU_obj, sim_state, MU_contribs, MU_impulse_train,
                        force, excitation, MU_rate_function, firing_response_fftvals, 
                        sampling) {
  
    for (t_idx in seq_along(force)) {
        sim_state = sample_MU_contribs(MU_obj,
                                       sim_state,
                                       force[t_idx],
                                       excitation[t_idx],
                                       MU_rate_function,
                                       firing_response_fftvals,
                                       sampling)
        
        ## Assign emg values
        for (electrode_id in seq_along(sim_state$contribs$emg)) {
            MU_contribs$emg[[electrode_id]][t_idx] <- sim_state$contribs$emg[[electrode_id]]
            ## Assign force value
            MU_contribs$force[t_idx] <- sim_state$contribs$force
        }
        
        # Update MU impulse trains.
        # Note that these must only be used for broad analyses of simulation behavior,
        # as they do not capture the exact course of the simulation due to being restrained
        # to firing events occuring at discrete sampling moments. In the actual simulation,
        # firing events may also occur between sampling instants.
        if (sim_state$time_since_last_firing <= 0) {
            # There is a firing event between the current and the next sample
            # (including the former, excluding the latter)
            fired = TRUE
        } else
            fired = FALSE
        MU_impulse_train[t_idx] = fired
        
    }
    
    list(impulse_train = MU_impulse_train,
         contribs = MU_contribs)
}


simulate_sample <- function(MUs, sim_state, load_vals, MUs_rate_functions,
                            firing_responses_fftvals_restructured, sampling, cluster=NULL) {
    
    data <- MUs %>% merge(sim_state, all=TRUE) %>% merge(load_vals) %>%
        merge(MUs_rate_functions) %>% merge(firing_responses_fftvals_restructured)

    ## calculate the contributions of all MUs to the different electrodes

    if (is.not.null(cluster)) {
      
      future::plan(future::cluster, workers=cluster)
      
      data$sim_state <- with(data,
                             future.apply::future_mapply(
                                    sample_MU_contribs,
                                    MU_obj = MU.obj,
                                    sim_state = sim_state,
                                    force = force,
                                    excitation = excitation,
                                    rate_function = MU_rate_function,
                                    firing_response_fftvals =
                                      firing_response_fftvals,
                                    MoreArgs = list(sampling = sampling),
                                    SIMPLIFY = FALSE))
      
    } else {
        
        data$sim_state <- with(data,
                             mapply(sample_MU_contribs,
                                    MU_obj = MU.obj,
                                    sim_state = sim_state,
                                    force = force,
                                    excitation = excitation,
                                    rate_function = MU_rate_function,
                                    firing_response_fftvals =
                                        firing_response_fftvals,
                                    MoreArgs = list(sampling = sampling),
                                    SIMPLIFY = FALSE))
    }

    data[, c("muscle", "MU", "sim_state")]
}


setup_sim_cluster <- function(num_cores) {

    if (file.exists("parlog_contrib.txt"))
        file.remove("parlog_contrib.txt")
        
    cl <- parallel::makeCluster(num_cores, outfile="") #outfile="parlog_contrib.txt")

    parallel::clusterExport(
        cl,
        varlist = list("sample_MU_contribs",
                       "calc_shifted_emg_response", "rnorm_bounded",
                       "is_fft_of_real_signal", "is.even", "time_shift_fft",
                       "is.nice.scalar", "calc_abs_rel_err",
                       "sum_unequal_lengths"),
        envir = environment())

    cl
}

#' Calculate the EMG and force contributions of a MU at the next sample.
#'
#' Take the previous simulation history into account; determine whether the MU
#' has fired between the last and this sample, and determine the EMG and force
#' contributions that this MU does in this measurement sample.
#'
#' Examines the time since the last firing event, the current load value, the
#' firing rate relationship of this MU and inter-spike variability to determine
#' whether the MU has fired since the last measurement sample.
#' Then, takes previous firing events and - if applicable - the current firing
#' event into account to calculate the contributions of this MU to overall EMG
#' and force measurements in this measurement sample.
#'
#' @param MU_obj A motor unit object containing information on force twitch
#'   characteristics.
#' @param sim_state A list containing elements `time_since_last_firing`,
#'   `Z_score`, `contribs` and `future_contribs`. Describes the current state
#'   of the simulation of this MU, i.e., summarizes the simulation history.
#' @param load_val Current value of the muscle load function; describes muscle
#'   activation.
#' @param rate_function The rate coding function of this MU.
#' @param firing_response_fftvals A list of num_electrodes vectors containing
#'   firing response fft values. The index of each vector indicates the
#'   corresponding electrode ID.
#' @param sampling An object providing information on the sampling parameters
#'   chosen. Usually constructed by the function create_sampling.
#' @return An updated `sim_state` of this MU.
#'
sample_MU_contribs <- function(MU_obj, sim_state, force, excitation, rate_function,
                              firing_response_fftvals, sampling) {

    ## -- Calculate current EMG and force contributions from past firing events
    for (electrode in seq_along(firing_response_fftvals)) {
        sim_state$contribs$emg[[electrode]] <-
            sim_state$future_contribs$emg[[electrode]][1]
        if (length(sim_state$future_contribs$emg[[electrode]]) > 1)
            sim_state$future_contribs$emg[[electrode]] <-
                tail(sim_state$future_contribs$emg[[electrode]], -1)
        else
            sim_state$future_contribs$emg[[electrode]] <- 0
    }

    sim_state$contribs$force <- sim_state$future_contribs$force[1]    
    if (length(sim_state$future_contribs$force) > 1)
        sim_state$future_contribs$force <-
            tail(sim_state$future_contribs$force, -1)
    else
        sim_state$future_contribs$force <- 0


    ## -- Detect and handle new firing events
    
    ## Nominal inter-spike interval
    isi_nom <- 1 / rate_function(excitation)
    
    # Special behavior if the MU is newly recruited
    if (is.infinite(isi_nom)) # MU currently not recruited
      sim_state$is_recruited <- FALSE
    else { # MU currently recruited
      if (!sim_state$is_recruited) {# MU is recruited now but wasn't before. This is what we want to deal with.
        sim_state$is_recruited <- TRUE
        # To prevent all MUs from firing simultaneously when rapidly increasing muscle activation, reset 
        # time_since_last_firing to a random time based on the current isi_nom. 
        # The "min(time_since_last_firing, ...)" part is to prevent wrong behavior when MUs are recruited, then
        # un-recruited for a very short time, and then re-recruited again. (Probably a very rare scenario, but still.)
        sim_state$time_since_last_firing <- min(sim_state$time_since_last_firing,
                                                runif(n=1, min=0, max=isi_nom)) 
      }
    }

    ## ISI varibility is dependent on current level of activation, see
    ## Moritz, Barry, Pascoe, Enoka (2005). Discharge Rate Variability...
    coeffvar <- 0.1 + 0.2 * exp(-(force - MU_obj$rec_thresh_force) * 100 / 2.5)
    
    isi_var <- isi_nom + isi_nom * coeffvar * sim_state$Z_score

    tdiff <- sim_state$time_since_last_firing + 2 / sampling$Fs - isi_var

    if (is.nan(tdiff))
        tdiff <- -1
    
    if (tdiff >= 0) {
        ## Between the current and the next sample, a firing event will occur!
        ##
        ## Due to changing muscle activation it is possible that tdiff > 1/Fs,
        ## i.e., it is detected that "there should have been a firing before the
        ## current sample". This is a) complicated to implement, and b) probably
        ## not necessary. Hence, we restrict tdiff to one sample.
        ## * Note that this is a result of _not_ using an iterative solution
        ## approach for equation (4) in Fuglevand et al. (1993), as proposed
        ## there. We do this for efficiency reasons.
        tdiff <- min(tdiff, 1 / sampling$Fs)
		
		## Do we want to allow firing instants at times that are not sampling times?
		## If no, fix tdiff to one sampling period.
		if (!sampling$unsampled_firings)
			tdiff <- 1 / sampling$Fs
        
        sim_state$time_since_last_firing <- tdiff - 1 / sampling$Fs
        sim_state$Z_score <- rnorm_bounded(min=-3.9, max=3.9, n_sigma=3.9)
        
        ## Calculate EMG contributions
        for (electrode in seq_along(firing_response_fftvals)) {
            fftvals <- firing_response_fftvals[[electrode]]
            emg_contrib <-
                calc_shifted_emg_response(
                    firing_response_fft = fftvals,
                    time_shift = tdiff,
                    sampling = sampling)
            emg_contrib_prev <- sim_state$future_contribs$emg[[electrode]]

            sim_state$future_contribs$emg[[electrode]] <-
                sum_unequal_lengths(emg_contrib_prev, emg_contrib)
        }
        
        ## Calculate force contribution
        t <- seq(tdiff, MU_obj$twitch_t_lead + 2 * MU_obj$twitch_t_half_rel,
                 1 / sampling$Fs)
        force_contrib <- MU_obj$force_twitch(t)
        gainfac <- MU_obj$gain_fac_fun(isi_var)
        force_contrib <- gainfac * force_contrib
        sim_state$future_contribs$force <-
            sum_unequal_lengths(sim_state$future_contribs$force, force_contrib)
        
    } else
        sim_state$time_since_last_firing <-
            sim_state$time_since_last_firing + 1 / sampling$Fs

    sim_state
}


initialize_sim_state <- function(num_electrodes) {
    sim_state <- list()
    sim_state$Z_score <- rnorm_bounded(min=-3.9, max=3.9, n_sigma=3.9)
    sim_state$time_since_last_firing <- Inf
    sim_state$is_recruited <- FALSE
    sim_state$future_contribs <- list()
    sim_state$future_contribs$force <- 0
    sim_state$future_contribs$emg <- list()
    for (i in seq(1, num_electrodes))
        sim_state$future_contribs$emg[[i]] <- 0

    sim_state
}


#' Time-shift a signal given by its FFT and compute the IFFT.
#'
#' Take FFT values of a signal, perform a small time-shift on it and
#' compute the inverse FFT.
#'
#' @note Requires the underlying signal to be purely real and thus the FFT to be
#'   hermitian.
#' 
#' @param firing_response_fft FFT values of a signal, computed at the
#'   frequencies given by sampling$fft_freqs.
#' @param time_shift Amount of time by which the signal should be shifted in
#'   time. Must be zero or positive and smaller than one sampling interval.
#' @param sampling An object providing information on the sampling paramters
#'   chosen. Usually constructed by the function create_sampling.
#' @return A real, numerical vector of length sampling$NFFT.
#' 
#' @examples
#' set.seed(1)
#' t <- seq(0.1, 1, 0.1)
#' y <- sin(2*pi*t)
#' resp <- semg.simulatr:::calc_shifted_emg_response(
#'     firing_response_fft = fft(y) / length(y), 1/20,
#'     sampling = semg.simulatr:::create_sampling(Fs = 10, NFFT = 10))
#' dev.new()
#' plot(t, resp)
#' points(t, y, col="red")
#' 
calc_shifted_emg_response <- function(firing_response_fft, time_shift,
                                      sampling) {

    stopifnot(is_fft_of_real_signal(firing_response_fft, sampling))
    stopifnot(0 <= time_shift, time_shift <= 1/sampling$Fs)

    abstol <- 1e-12
    reltol <- 1e-6
    
    if (is.even(sampling$NFFT))
        if (abs(firing_response_fft[sampling$NFFT/2 + 1]) > abstol)
            stop(paste("Non-zero coefficient for frequency Fs/2 found. This",
                       "leads to problems with the time-shift operation on",
                       "real signals. Using an odd NFFT solves this problem."))
        
    NFFT <- sampling$NFFT
    
    firing_response_fft_shifted <-
        time_shift_fft(fft_vals = firing_response_fft,
                       sampling = sampling,
                       time_shift = time_shift)

    firing_response_cmplx <- fft(firing_response_fft_shifted, inverse = TRUE) * sampling$Fs / NFFT
        
    stopifnot(length(firing_response_cmplx) == NFFT)
    
    ## Numerical inaccuracies lead to non-zero imaginary parts. Cut these
    ## off, as they are purely computational artifacts.
    rel_err_max <-
        (Im(firing_response_cmplx) / Re(firing_response_cmplx)) %>% abs %>% max
    
    if(rel_err_max > reltol)
        warning(paste("Relative error greater than expected, was", rel_err_max, "."))
    
    firing_response <- Re(firing_response_cmplx)
    
    firing_response
}


append_MU_sample <- function(MU_contribs, sim_state, t_idx) {
    for (i in as.numeric(rownames(MU_contribs))) {
        sim_state_row <-
            which(sim_state$muscle == MU_contribs$muscle[[i]] &
                  sim_state$MU == MU_contribs$MU[[i]])
        ## Assign emg values
        for (electrode_id in
             seq_along(sim_state$sim_state[[sim_state_row]]$contribs$emg))
            MU_contribs$contribs[[i]]$emg[[electrode_id]][t_idx] <-
                sim_state$sim_state[[sim_state_row]]$contribs$emg[[electrode_id]]
        ## Assign force value
        MU_contribs$contribs[[i]]$force[t_idx] <-
            sim_state$sim_state[[sim_state_row]]$contribs$force
    }
    
    MU_contribs
}


aggregate_sim_results <- function(muscles, muscle_contribs, surface_potentials) {

    ## Extract time column and initialize df
    df <- muscle_contribs$emg %>% dplyr::filter(muscle==1, electrode==1) %>% `[`('time')

    ## Add force columns
    for (muscle_id in muscles$muscle) {
        colname <- paste('Fmus_', muscles$muscle.obj[[muscle_id]]$identifier,
                         sep="")
        df[colname] <- muscle_contribs$force %>%
            dplyr::filter(muscle==muscle_id) %>% `[[`('force')
    }

    ## Add semg columns
    for (electrode_id in unique(surface_potentials$electrode)) {
        colname <- paste('semg_', sprintf("%02d", electrode_id),
                         sep="")
        df[colname] <- surface_potentials %>%
            dplyr::filter(electrode==electrode_id) %>% `[[`('potential')
    }

    df
}