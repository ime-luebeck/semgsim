#' Simulate a surface electromyogram and force resulting from muscle acticity
#'
#' Given a virtual experimental set-up, simulate the sEMG and force data
#' resulting at a set of electrodes from activity in several muscles.
#'
#' Implements a physiologically sound and detailed mathematical model of motor
#' unit activity and the volume conductors represented by skin and fat tissues
#' between muscle and electrodes. Many (in fact, practically all except the
#' underlying base model) physiological parameters can be configured by
#' providing an appropriate config file.
#'
#' For more details regarding the simulation model, refer to the paper
#' Petersen and Rostalski (2019): "A Comprehensive Mathematical Model of Motor
#' Unit Pool Organization, Surface Electromyography and Force Generation",
#' doi: 10.3389/fphys.2019.00176.
#' preint
#' @note Currently, effects resulting from fatigue are not implemented.
#' 
#' @param config_file Specifies the simulation config file to be read. Must
#'   specify a number of physiological parameters. See section Details for more
#'   information.
#' @param out_file Optional. If given, specifies a file to save the simulation
#'   \code{environment()}. Results can be loaded via \code{readRDS(out_file)}
#'   afterwards.
#' @param MoreArgs A list of parameter name - value pairs. Can be used to
#'   provide necessary model parameters that are missing in the
#'   \code{config_file}.
#' @param num_cores Number of cores to use for various parts of the simulatino
#'   that can be performed in parallel. Optional, default=1.
#' @author Eike Petersen
#' @export
#' 
simulate_muscle <- function(config_file, deterministic = FALSE, out_file = NULL,
                            num_cores = parallelly::availableCores(omit=1)) {

  message(paste0(whoami(2)," : Simulation started with <",num_cores,"> cores ... "))
  
  #progressr::handlers(global = TRUE)
  #progressr::handlers("progress")
  
  if (deterministic)
        set.seed(1)

    model <- setup_muscles(config_file = config_file,
                           deterministic = deterministic,
                           out_file = out_file)

    model <- calculate_semg_TFs(model = model,
                                deterministic = deterministic,
                                out_file = out_file,
                                num_cores = num_cores)

    model <- simulate_model(model = model,
                            muscle_force_functions = model$force_functions,
                            time_limits = model$time_limits,
                            deterministic = deterministic,
                            out_file = out_file,
                            num_cores = num_cores,
                            parallel = "chunk")

    if (is.not.null(out_file))
        saveRDS(model, out_file)

    invisible(model)
}


#' @export
setup_muscles <- function(config_file, deterministic = FALSE, out_file = NULL) {

    message(paste0(whoami(2)," : Running ...")) 
  
    if (deterministic)
        set.seed(1)
    
    source(config_file, local = TRUE, chdir = T)

    stopifnot(length(muscles_list) == length(force_functions))

    muscles <- convert_muscle_list_to_df(muscles_list)

    MUs <- generate_MUs(muscles = muscles)

    MUs_rate_functions <-
        calc_MU_firing_rate_funcs(muscles = muscles, MUs = MUs)
    
    act_funcs_df = calc_activation_transformations(MUs = MUs, MUs_rate_functions = MUs_rate_functions)

    MUs <- update_MU_rec_threshs(MUs = MUs, act_funcs_df = act_funcs_df)
    
    all_vars <- ls()
    model <- mget(all_vars[all_vars != "model"])

    if (is.not.null(out_file))
        saveRDS(model, out_file)

    invisible(model)
}


#' @export
calculate_semg_TFs <- function(model, deterministic = FALSE, out_file = NULL,
                               num_cores = 1) {

    message(paste0(whoami(2)," : Running ...")) 
  
    if (deterministic)
        set.seed(1)

    for (name in names(model))
        if (!(name %in% c("deterministic", "out_file", "num_cores")))
            assign(name, model[[name]])


	if (!exists('unsampled_firings'))
		sampling <- setup_sampling(muscles = muscles,
								   f_target = Fs,
								   psi_length = psi$length,
								   max_freq_dom_sample_dist =
									   max_freq_dom_sample_dist)
	else
		sampling <- setup_sampling(muscles = muscles,
								   f_target = Fs,
								   psi_length = psi$length,
								   max_freq_dom_sample_dist =
									   max_freq_dom_sample_dist,
								   unsampled_firings = unsampled_firings)
    
    TFs_MU_input_to_surf_potential <-
        calc_MU_electrode_time_TFs(MUs = MUs,
                                   electrodes = electrodes,
                                   vol_conductor = volume_conductor,
                                   freqs = sampling$freqs_to_calc,
                                   num_cores = num_cores,
                                   B_sampling = B_sampling)

    MU_firing_responses_fftvals <-
        calc_MU_firing_responses_fftvals(MUs = MUs,
                                         TFs_MU_input_to_surf_potential =
                                             TFs_MU_input_to_surf_potential,
                                         psi_trafo = psi$transform,
                                         num_cores = num_cores,
                                         sampling = sampling)

    all_vars <- ls()
    model <- mget(all_vars[all_vars != "model"])

    if (is.not.null(out_file))
        saveRDS(model, out_file)

    invisible(model)
}



#' @export
simulate_model <- function(model, muscle_force_functions, time_limits,
                           deterministic = FALSE, out_file = NULL,
                           num_cores = 1, parallel = c("samplewise","chunk")) {

    message(paste0(whoami(2)," : Running ...")) 
  
    if (deterministic)
        set.seed(1)

    for (name in names(model))
        if (!(name %in% c("deterministic", "out_file", "num_cores")))
            assign(name, model[[name]])

    t <- seq(time_limits[1], time_limits[2], 1 / Fs)

    ## Initialize simulation state
    sim_state <- MUs[, c("muscle", "MU")]
    MU_contribs <- MUs[, c("muscle", "MU")]
    MU_impulse_trains <- MUs[, c("muscle", "MU")]
    for (i in seq(1, nrow(sim_state))) {
         sim_state$sim_state[[i]] <- initialize_sim_state(num_electrodes =
                                                              nrow(electrodes))
         MU_contribs$contribs[[i]] <- list(force=rep(NA, length(t)), emg=list())
         for(j in seq(1, nrow(electrodes)))
             MU_contribs$contribs[[i]]$emg[[j]] <- rep(NA, length(t))
         MU_impulse_trains$impulse_train[[i]] <- rep(NA, length(t))
    }
    
    # Restructure firing response fftvals for the following analyses
    MU_firing_responses_fftvals_restructured <-
      restructure_firing_responses_fftvals(MU_firing_responses_fftvals)

    if (parallel[1] == "samplewise") {
        
        if (num_cores > 1)
            ## seems to be inefficient anyway
            ##cl <- setup_sim_cluster(num_cores)
            cl <- NULL
        else
            cl <- NULL
        
        #p <- progressr::progressor(steps=length(t))
        p <- txtProgressBar(max=length(t),style=3)
        
        for (t_idx in seq_along(t)) {
            
            #p(message = sprintf("x=%g", t_idx))
            setTxtProgressBar(p, x)
            
            force_vals <- lapply(seq_along(muscle_force_functions), 
                                 function(forcefun_idx) muscle_force_functions[[forcefun_idx]](t[t_idx])) %>% 
                simplify2array
            
            excitation_vals <-
                lapply(seq_along(muscle_force_functions), 
                       function(forcefun_idx) act_funcs_df[[forcefun_idx, 'force_to_act_fun']](eval(muscle_force_functions[[forcefun_idx]](t[t_idx])))) %>% 
                simplify2array
            
            load_vals <- data.frame(muscle = seq_along(muscle_force_functions),
                                    force = force_vals,
                                    excitation = excitation_vals)
            
            sim_state <- simulate_sample(MUs, sim_state, load_vals,
                                         MUs_rate_functions,
                                         firing_responses_fftvals_restructured =
                                             MU_firing_responses_fftvals_restructured,
                                         sampling = sampling,
                                         cluster = cl)
            
            MU_contribs <- append_MU_sample(MU_contribs, sim_state, t_idx)
            
            # Update MU impulse trains.
            # Note that these must only be used for broad analyses of simulation behavior,
            # as they do not capture the exact course of the simulation due to being restrained
            # to firing events occuring at discrete sampling moments. In the actual simulation,
            # firing events may also occur between sampling instants.
            for (i in seq(1, nrow(MU_impulse_trains))) {
                if (sim_state$sim_state[[i]]$time_since_last_firing <= 0) {
                    # There is a firing event between the current and the next sample
                    # (including the former, excluding the latter)
                    fired = TRUE
                } else
                    fired = FALSE
                MU_impulse_trains$impulse_train[[i]][t_idx] = fired
            }
            
            ##print(t_idx)
        }
        
        close(p)
        ##if (num_cores > 1)
        ##    parallel::stopCluster(cl)
        
    } else if (parallel[1] == "chunk") {
        
        load_vals <- data.frame(muscle = seq_along(muscle_force_functions))
        load_vals$force <- lapply(seq_along(muscle_force_functions), 
                                  function(forcefun_idx) muscle_force_functions[[forcefun_idx]](t))
        load_vals$excitation <- lapply(seq_along(muscle_force_functions), 
                                       function(forcefun_idx) act_funcs_df[[forcefun_idx, 'force_to_act_fun']](eval(muscle_force_functions[[forcefun_idx]](t))))
        
        data <- MUs %>% 
            merge(sim_state) %>% 
            merge(MU_contribs) %>% 
            merge(MU_impulse_trains) %>% 
            merge(load_vals) %>% 
            merge(MUs_rate_functions) %>% 
            merge(MU_firing_responses_fftvals_restructured)
        
        aux = simulate_MUs(data, sampling, num_cores)
        
        MU_contribs$contribs <- lapply(aux, function(x)x$contribs)
        MU_impulse_trains$impulse_train <- lapply(aux, function(x)x$impulse_train)
        
    } else {
        stop(paste0("Wrong value <",parallel,"> for parameter <parallel>."))
    }
    
	  ## contains $force and $emg 
    muscle_contribs <- sum_MU_contribs(MU_firing_contribs = MU_contribs,
                                       time = t,
                                       sampling = sampling)

    surface_potentials <- sum_muscle_contribs(muscle_firing_contribs = muscle_contribs$emg)

    sim_results <- aggregate_sim_results(muscles, muscle_contribs, surface_potentials)
    
    all_vars <-  setdiff(ls(), c("model", "aux", "sim_state", "data"))   # alternativ: diese Variablen vor ls() aus dem Environment löschen
    res <- mget(all_vars)

    if (is.not.null(out_file))
        saveRDS(res, out_file)
    
    invisible(res)
}
