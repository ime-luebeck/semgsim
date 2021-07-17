##' Apply electrode configurations
##'
##' Calculate combinations of different, monopolar electrode signals.
##'
##' This allows for the implementation of differential measurements, etc.
##'
##' @note This implementation is probably highly inefficient in every regard.
##' @param muscle_potentials A \code{data.frame} as usually created by a call to
##'   \code{\link{sum_muscle_contribs}}, containing the contributions of
##'   different muscles to the monopolar signals measured at each electrode.
##' @param electrode_configs A list of vectors, where each vector should have as
##'   many elements as there are electrodes. The elements of the vector
##'   determine the weights given to the different monopolar signals for
##'   calculating the combined signal.
##' @return A \code{data.frame} containing the contributions of the different
##'   muscles to the signals obtained using the different electrode
##'   configurations. These are still specified via the categorical variable
##'   \code{electrode} (instead of using "electrode_config", or something
##'   similar), in order to allow for easy usage of the same analysis functions.
apply_electrode_configs <- function(muscle_potentials, electrode_configs) {

    potentials_split <- split(muscle_potentials, muscle_potentials$electrode)

    for (i in seq_along(potentials_split)) {
        potentials_split[[i]]$electrode <- NULL
        potential_col_idx <- which(colnames(potentials_split[[i]]) == "potential")
        colnames(potentials_split[[i]])[potential_col_idx] <- paste("potential", i, sep = "")
    }

    potentials_comb <- Reduce(merge, potentials_split)

    potentials_lst <- list()

    for (i in seq_along(electrode_configs)) {

        potential <- rep(0, nrow(potentials_comb))
        
        for (j in seq_along(electrode_configs[[i]])) {
            potential_name <- paste("potential", j, sep = "")
            potential <- potential +
              potentials_comb[[potential_name]] * electrode_configs[[i]][j]
        }

        potentials_lst[[i]] <- data.frame(electrode = rep(i, nrow(potentials_comb)),
                                          time = potentials_comb$time,
                                          potential = potential)

        # Re-add "muscle", "MU" columns to be able to handle MU/muscle contrib dfs as well.
        for (colname in colnames(potentials_comb))
            if (colname %in% c("muscle", "MU"))
                potentials_lst[[i]][colname] <- potentials_comb[colname]
        }

    res_final <- potentials_lst %>% (dplyr::bind_rows)

    res_final
}


post_process_wrap <- function(potentials, post_process_fun, ...) {

    if ("muscle" %in% colnames(potentials))
    
        potentials_split <- split(potentials,
                                list(potentials$muscle,
                                     potentials$electrode))
    else
        potentials_split <- split(potentials, potentials$electrode)        

    potentials_processed_lst <-
        potentials_split %>%
          plyr::llply(function(potential) post_process_fun(potential, ...))

    result <- list()
    
    potentials_processed_lst %>% dplyr::bind_rows(.)
}


post_process_rms <- function(muscle_potentials, win_length = 0.25) {

    muscle_potentials <- muscle_potentials %>% plyr::arrange(time)

    dT <- muscle_potentials$time[2] - muscle_potentials$time[1]
    win_samples <- ceiling(win_length / dT)

    win_samples <- win_samples + 1

    muscle_potentials$rms <- muscle_potentials$potential %>%
      magrittr::raise_to_power(2) %>%
      stats::filter(filter = rep(1 / win_samples, win_samples),
             method = "convolution",
             sides = 2,
             circular = FALSE) %>% sqrt %>% replace_na(0)

	if ("muscle" %in% colnames(muscle_potentials))
		muscle_potentials[, c("muscle", "electrode", "time", "rms")]
	else
		muscle_potentials[, c("electrode", "time", "rms")]
}

get_post_process_config_rms <- function(for_latex = FALSE) {
    if (for_latex)
        ylabel <- "Root Mean Square of Electric Potential (V)"
    else
        ylabel <- "Root Mean Square of Electric Potential (V)"
    list(analysis_label = "RMS",
         x_name = "time",
         x_label = "Time (s)",
         y_name = "rms",
         y_label = ylabel,
         geom_fun = ggplot2::geom_line)
}

post_process_fourier <- function(muscle_potentials) {

    n <- nrow(muscle_potentials)
    
    ## Detrend, demean
	fit <- lm(potential ~ time, muscle_potentials)
	potentials_cleaned <- muscle_potentials$potential - fit$coefficients[1] -
		fit$coefficients[2] * muscle_potentials$time

    ## trim tails before first and beyond last zero crossing
    pos_mask <- potentials_cleaned > 0
    zero_crossing_bool <- pos_mask[1:(n-1)] != pos_mask[2:n]
    start_idx <- min(which(zero_crossing_bool == TRUE)) + 1
    end_idx <- max(which(zero_crossing_bool == TRUE))
    potentials_trimmed <- potentials_cleaned[start_idx:end_idx]    
    
    ## Pad with zeros
    NFFT <- nextn(length(potentials_trimmed))
    diff <- NFFT - length(potentials_trimmed)
    potentials_padded <- c(potentials_trimmed, rep(0, diff))

    ## Perform actual FFT calculations
    N <- ceiling((NFFT + 1) / 2)
    dT <- muscle_potentials$time[2] - muscle_potentials$time[1]
    Fs <- 1 / dT
    fft_vals <- fft(potentials_padded) / NFFT
    fft_freqs <- calc_fft_freqs(NFFT, Fs = Fs)

    if ("muscle" %in% colnames(muscle_potentials)) {
        data.frame(muscle = rep(muscle_potentials$muscle[1], N),
                   electrode = rep(muscle_potentials$electrode[1], N),
                   frequency = fft_freqs[fft_freqs >= 0],
                   amplitude = abs(fft_vals[fft_freqs >= 0]))
    } else {
        data.frame(electrode = rep(muscle_potentials$electrode[1], N),                   
                   frequency = fft_freqs[fft_freqs >= 0],
                   amplitude = abs(fft_vals[fft_freqs >= 0]))
    }
}

get_post_process_config_fourier <- function(for_latex = FALSE) {
    if (for_latex)
        ylabel <- "Amplitude (V)"
    else
        ylabel <- "Amplitude (V)"
    
    list(analysis_label = "Fourier Analysis",
         x_name = "frequency",
         x_label = "Frequency (Hz)",
         y_name = "amplitude",
         y_label = ylabel,
         geom_fun = function() ggplot2::geom_freqpoly(stat = "identity"))
}

post_process_density <- function(muscle_potentials) {
    muscle_potentials
}

get_post_process_config_density <- function(for_latex = FALSE) {
    if (for_latex)
        xlabel <- "EMG Amplitude (V)"
    else
        xlabel <- "EMG Amplitude (V)"
    
    list(analysis_label = "Amplitude PDF",
         x_name = "potential",
         x_label = xlabel,
         y_label = "Probability Density",
         geom_fun = ggplot2::geom_density,
         addon_fun = function(plot_obj, df) {

             vals <- with(df, seq(min(potential), max(potential), length = 100))
			 
             if ("muscle" %in% colnames(df))
                 vars <- c("muscle", "electrode")
             else
                 vars <- c("electrode")
             
             ref_dens <-
                 plyr::ddply(.data = df,
                             .variables = vars,
                             .fun = function(df) {
                                 data.frame(potential = vals,
                                            norm =
                                              dnorm(vals,
                                                    mean(df$potential),
                                                    sd(df$potential)),
                                            lap =
                                              dlaplace(vals,
                                                       mean(df$potential),
                                                       sd(df$potential) / sqrt(2))
                                            )
                             })

             plot_obj +
               ggplot2::geom_line(aes(y = norm),
                                  data = ref_dens, linetype = "dashed") +
               ggplot2::geom_line(aes(y = lap),
                                  data = ref_dens, linetype = "dotted")
         })
}


calc_emg_force_rel <- function(surface_potentials, stairs_start, step_length,
                               num_steps, transient_length, electrodes=NA,
                               electrode_confs=NULL) {
    
    stopifnot(is.nice.scalar(stairs_start),
              is.nice.scalar(step_length),
              is.nice.scalar(num_steps),
              num_steps > 0,
              is.nice.scalar(transient_length),
              transient_length < step_length,
              is.data.frame(surface_potentials))
    
    if (!is.numeric(electrodes))
      electrodes <- unique(surface_potentials$electrode)
    else
      electrodes <- unique(c(0, electrodes)) # force "electrode" should be kept in any case
  
    surface_potentials <- surface_potentials %>%
      dplyr::filter(electrode %in% electrodes)
  
    if (!is.null(electrode_confs)) {
      force_df <- surface_potentials %>% dplyr::filter(electrode == 0)
      surface_potentials <- apply_electrode_configs(surface_potentials,
                                                    electrode_confs)
      ## Re-add the force "electrode"
      surface_potentials <- dplyr::bind_rows(surface_potentials, force_df)
      
      # Re-determine electrodes to be considered
      electrodes <- unique(surface_potentials$electrode)
    }
    
    dfs <- list()
    i <- 1

    ## Calculate mean force and emg in each staircase step and for each different electrode (setup)
    for (step in seq(1, num_steps)) {
        force_mean <- 0
        for (ele in sort(electrodes)) {
            subset <- surface_potentials %>% dplyr::filter(electrode == ele) %>%
                dplyr::filter(time >= stairs_start + (step-1)*step_length + transient_length & time < stairs_start + step*step_length)
            meanval <- mean(abs(subset$potential))
            if (ele == 0)
                ## electrode 0 is the force signal
                force_mean <- meanval
            else {
                dfs[[i]] <- data.frame(step = step,
                                       force_mean = force_mean,
                                       emg_mav = meanval,
                                       electrode = ele)
                i <- i+1
            }
        }
    }

    emg_force_rels <- plyr::rbind.fill(dfs)

    ## Normalize each column individually
    emg_force_rels$force_mean <- (emg_force_rels$force_mean - min(emg_force_rels$force_mean)) / (max(emg_force_rels$force_mean) - min(emg_force_rels$force_mean))
    for (ele in electrodes) {
        if (ele == 0)
            next
        subset <- emg_force_rels %>% dplyr::filter(electrode == ele)
        emg_force_rels[emg_force_rels$electrode == ele,]$emg_mav <-
            (emg_force_rels[emg_force_rels$electrode == ele,]$emg_mav - min(subset$emg_mav)) / (max(subset$emg_mav) - min(subset$emg_mav))
    }
    
    emg_force_rels <- emg_force_rels %>%
      dplyr::mutate(electrode = as.factor(electrode))

    if (!all(0 <= emg_force_rels$emg_mav)) {
      print("Here I am")
    }
    stopifnot(all(0 <= emg_force_rels$force_mean),
              all(emg_force_rels$force_mean <= 1),
              all(0 <= emg_force_rels$emg_mav),
              all(emg_force_rels$emg_mav <= 1))
    
    emg_force_rels
}


calc_emg_force_rel_summary <- function(list_of_surface_potentials, stairs_start, step_length,
                                       num_steps, transient_length, electrodes=NA,
                                       electrode_confs=NULL) {
  
  stopifnot(is.data.frame(list_of_surface_potentials[[1]]))
  emg_force_rels_lst = list()
  for (ii in seq_along(list_of_surface_potentials))
    emg_force_rels_lst[[ii]] <- calc_emg_force_rel(
      surface_potentials=list_of_surface_potentials[[ii]],
      stairs_start = stairs_start, 
      step_length = step_length,
      num_steps = num_steps, 
      transient_length = transient_length, 
      electrodes = electrodes,
      electrode_confs = electrode_confs)
  emg_force_rels <- dplyr::bind_rows(emg_force_rels_lst)
  emg_force_rels <- emg_force_rels %>% dplyr::group_by(step, electrode) %>%
    dplyr::summarise_all(list("mean", "sd"))
  emg_force_rels <- emg_force_rels %>% dplyr::rename(force_mean = force_mean_mean)
  emg_force_rels <- emg_force_rels %>% dplyr::rename(emg_mav = emg_mav_mean)
  emg_force_rels <- emg_force_rels %>%
    dplyr::mutate(electrode = as.factor(electrode))
  emg_force_rels
}


calc_force_variability <- function(force_contribs, stairs_start, step_length,
                                   num_steps, transient_length) {
  
  stopifnot(is.nice.scalar(stairs_start),
            is.nice.scalar(step_length),
            is.nice.scalar(num_steps),
            num_steps > 0,
            is.nice.scalar(transient_length),
            transient_length < step_length,
            is.data.frame(force_contribs))
  
  muscles <- unique(force_contribs$muscle)
  
  dfs <- list()
  i <- 1
  
  ## Calculate force variability in each staircase step for each muscles
  for (step in seq(1, num_steps)) {
    for (muscle_id in sort(muscles)) {
      subset <- force_contribs %>% dplyr::filter(muscle == muscle_id) %>%
        dplyr::filter(time >= stairs_start + (step-1)*step_length + transient_length & time < stairs_start + step*step_length)
      force_mean = mean(subset$force)
      force_sd = sd(subset$force)
      dfs[[i]] <- data.frame(step = step,
                             force_mean = force_mean,
                             force_cv = force_sd/force_mean,
                             force_sd = force_sd,
                             muscle = muscle_id)
      i <- i+1
    }
  }
  
  force_variability <- plyr::rbind.fill(dfs)
  
  ## Normalize force mean and sd, _don't_ normalize cv as this is already a relative measure
  for (muscle_id in muscles) {
    subset <- force_variability %>% dplyr::filter(muscle == muscle_id)
    force_variability[force_variability$muscle == muscle_id,]$force_mean <-
      (force_variability[force_variability$muscle == muscle_id,]$force_mean - min(subset$force_mean)) / (max(subset$force_mean) - min(subset$force_mean))
    force_variability[force_variability$muscle == muscle_id,]$force_sd <-
      (force_variability[force_variability$muscle == muscle_id,]$force_sd - min(subset$force_sd)) / (max(subset$force_sd) - min(subset$force_sd))
  }
  
  stopifnot(all(0 <= force_variability$force_cv),
            all(0 <= force_variability$force_mean),
            all(force_variability$force_mean <= 1),
            all(0 <= force_variability$force_sd))
  
  force_variability
}