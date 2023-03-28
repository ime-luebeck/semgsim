#--------------------------------------------------------------------
#' Simulate a surface electromyogram and force resulting from muscle activity
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
#' 
#' @note Currently, effects resulting from fatigue are not implemented.
#' 
#' @param config_file Specifies the simulation config file to be read.
#'   Must specify a number of physiological parameters including the
#'   force function for each muscle. See section Details for more information.
#' @param deterministic Option sets seed for reproducible random 
#'   numbers.
#' @param out_file Optional. If given, specifies a file to save 
#'   the simulation \code{environment()}. Results can be loaded via 
#'   \code{readRDS(out_file)} afterwards.
#' @param num_cores Optional. Number of cores to use for various 
#'   parts of the simulation that can be performed in parallel. 
#'   Default = all available cores minus 1.
#'   
#' @author Eike Petersen
#' @export
#--------------------------------------------------------------------
#' @examples
#' model <- simulate_muscle("config/abdominals.r", deterministic=TRUE, num_cores=3, out_file = "simout/abdominals.RDS")
#' model <- simulate_muscle("config/integration.r", deterministic = TRUE) 
#--------------------------------------------------------------------
simulate_muscle <- function(config_file, deterministic = FALSE, out_file = NULL,
                            num_cores = future::availableCores(omit=1)) {
  
  message(paste0(whoami(2)," : Simulation started with <",num_cores,"> cores ... "))
  
  #progressr::handlers(global = TRUE)
  #progressr::handlers("progress")
  
  model <- setup_muscles(config_file = config_file,
                         deterministic = deterministic)
  
  model <- calculate_semg_TFs(model = model,
                              deterministic = deterministic,
                              num_cores = num_cores)
  
  model <- simulate_model(model = model,
                          muscle_force = model$force_functions,
                          time_limits = model$time_limits,
                          deterministic = deterministic,
                          num_cores = num_cores,
                          parallel = "chunk")
  
  if (is.not.null(out_file))
    saveRDS(model, out_file)
  
  invisible(model)
}



#--------------------------------------------------------------------
#' Muscle setup
#' 
#' @export
#--------------------------------------------------------------------
#' @examples
#' model <- setup_muscles("config/integration.r", deterministic = TRUE) 
#--------------------------------------------------------------------
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


#--------------------------------------------------------------------
#' Calculate transfer functions
#' 
#' @note  The config file may "link" [calculate_TFs_Far_Mer()] to 
#'   `calc_MU_electrode_time_TFs()` which is called from below. 
#' @export
#--------------------------------------------------------------------
#' @examples
#' model <- calculate_semg_TFs(setup_muscles("config/integration.r", deterministic = TRUE), deterministic = TRUE) 
#--------------------------------------------------------------------
calculate_semg_TFs <- function(model, deterministic = FALSE, out_file = NULL,
                               num_cores = future::availableCores(omit=1)) {
  
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


#--------------------------------------------------------------------
#' Run simulation based on existing `model`
#' 
#' The simulation is based on already existing setup and transfer 
#' functions available via parameter `model`. 
#' 
#' @details
#' This simulation is equivalent to [init.simulate_model_stepwise()]
#' below except that execution is performed in one step. Moreover, 
#' motor unit (MU) impulse trains and contribs are saved in `out_file`
#' or are returned. 
#' 
#' @param model Produced by [calculate_semg_TFs()].
#' @param muscle_force Either a list of muscle force functions 
#'   (one function per muscle) or a matrix with dimension
#'   \[time x muscles\] or a vector (only one muscle).
#' @param time_limits Optional. Two-elements vector limiting the 
#'   time axis. If `NULL` then `muscle_force` is assumed to be a 
#'   matrix or data frame which determines the time axis.
#' @param deterministic Option sets seed for reproducible random 
#'   numbers. Note that the results are not identical if
#'   [simulate_model()] is called with `parallel="sample"` or with 
#'   `parallel="chunk"`. The reason is that the sequence of calling 
#'   the random number generator is different.
#' @param out_file Optional. If given, specifies a file to save 
#'   the simulation \code{environment()}. Results can be loaded via 
#'   \code{readRDS(out_file)} afterwards.
#' @param num_cores Optional. Number of cores to use for various 
#'   parts of the simulation that can be performed in parallel. 
#'   Default = all available cores minus 1.
#' @param parallel Specifies whether calculation is performed
#'   `sample` by sample or as one `chunk`. The `future` framework 
#'   is used for parallel computation. `sample` by sample performs
#'   simulation samplewise internally. However, it should not 
#'   be confounded with [init.simulate_model_stepwise()] which 
#'   provides an interface to explicitly allow stepwise (e.g.,
#'   samplewise) calculatione, e.g., for realtime processing.
#' 
#' @export
#--------------------------------------------------------------------
#' @examples
#' # This may take long because transfer functions have to be calculated, initially:
#' model <- calculate_semg_TFs(setup_muscles("config/integration.r", deterministic = TRUE), deterministic = TRUE) 
#' saveRDS(model, file="integration.rds")
#' 
#' # Now run the simulation:
#' model <- simulate_model(model, model$force_functions, model$time_limits, num_cores=1, parallel="chunk", deterministic=TRUE)
#' 
#' # Run simulation in "chunk" mode on input array of force data:
#' force_vals = as.data.frame(model$load_vals$force)[1:20,]
#' model <- simulate_model(model, as.matrix(force_vals), num_cores=1, parallel="chunk", deterministic=TRUE)
#' 
#' # Run simulation in "sample" mode:
#' model <- simulate_model(model, as.matrix(force_vals), num_cores=1, parallel="sample", deterministic=TRUE)
#'
#' # return muscle force
#' model$muscle_contribs$force
#--------------------------------------------------------------------
simulate_model <- function(model, muscle_force = NULL, time_limits = NULL,
                           deterministic = FALSE, out_file = NULL,
                           num_cores = future::availableCores(omit=1), 
                           parallel = c("sample","chunk")) {
  
  message(paste0(whoami(2)," : Running ...")) 
  if (deterministic)
    set.seed(1)
  
  for (name in names(model))      # add all components of model to current environment ...
    if (!(name %in% c("muscle_force", "time_limits", "deterministic", "out_file", "num_cores", "parallel")))  # ... except these
      assign(name, model[[name]])

  if (is.null(muscle_force)) {    # use default muscle force functions provided in the model configuration
    muscle_force = model$force_functions  
  }
  
  if (is.list(muscle_force) &&    # muscle_force is a list of functions (as many as muscles), the force values are generated based on the time axis
      is.function(muscle_force[[1]]) && 
      !is.null(time_limits)) 
  {
    stopifnot(length(model$muscles$muscle) == length(muscle_force))
    t <- seq(time_limits[1], time_limits[2], 1 / Fs)
    force_vals <- sapply(model$muscles$muscle, function(forcefun_idx) muscle_force[[forcefun_idx]](t))
  } else if (is.numeric(muscle_force)) {  # muscle force is a vector or matrix [time x muscles] containing the force values which defines the time axis
    if (is.vector(muscle_force)) {        # if muscle force is a vector convert to matrix with 1 column
      muscle_force <- matrix(muscle_force, ncol=1)
    }
    stopifnot(length(model$muscles$muscle) == ncol(muscle_force))
    force_vals <- muscle_force
    time_limits <- c(0, (nrow(force_vals)-1) / Fs)    # define time axis range based on number of rows (samples) and sampling frequency
    t <- seq(time_limits[1], time_limits[2], 1 / Fs)  # construct the time axis
  } else {
    stop("muscle_force is neither a list of force functions nor a matrix of values")
  }
  
  excitation_vals <- sapply(model$muscles$muscle, 
                            function(forcefun_idx) act_funcs_df[[forcefun_idx, 'force_to_act_fun']](force_vals[,forcefun_idx]))
  
  ## Initialize simulation state
  sim_state <- MUs[, c("muscle", "MU")]
  MU_contribs <- MUs[, c("muscle", "MU")]
  MU_impulse_trains <- MUs[, c("muscle", "MU")]
  for (i in seq(1, nrow(sim_state))) {
    sim_state$sim_state[[i]] <- initialize_sim_state(num_electrodes = nrow(electrodes))
    MU_contribs$contribs[[i]] <- list(force=rep(NA, length(t)), emg=list())
    for(j in seq(1, nrow(electrodes)))
      MU_contribs$contribs[[i]]$emg[[j]] <- rep(NA, length(t))
    MU_impulse_trains$impulse_train[[i]] <- rep(NA, length(t))
  }
  
  # Restructure firing response fftvals for the following analyses
  MU_firing_responses_fftvals_restructured <-
    restructure_firing_responses_fftvals(MU_firing_responses_fftvals)
  
  if (parallel[1] == "sample") {
    
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
      setTxtProgressBar(p, t_idx)
      
      load_vals <- data.frame(muscle = model$muscles$muscle,
                              force = force_vals[t_idx,],
                              excitation = excitation_vals[t_idx,])
      
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
      
    }
    
    close(p)
    ##if (num_cores > 1)
    ##    parallel::stopCluster(cl)
    
  } else if (parallel[1] == "chunk") {
    
    load_vals <- data.frame(muscle = model$muscles$muscle)
    load_vals$force <- as.list(as.data.frame(force_vals))
    load_vals$excitation <- as.list(as.data.frame(excitation_vals))
    
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


#--------------------------------------------------------------------
#' Initialize stepwise simulation
#' 
#' It is based on already existing setup and transfer functions 
#' available via parameter `model`. 
#' 
#' @details 
#' This simulation is equivalent to `simulate_model()` above except
#' that execution is performed stepwise (sample by sample) and motor
#' unit (MU) impulse trains and contribs are not saved. The examples
#' below show how a simulation `instance` is created and simulation is
#' run.
#' 
#' @param model Produced by `calculate_semg_TFs()`.
#' @param deterministic Option sets seed for reproducible random 
#'   numbers. Note that the the results are not identical if
#'   `instance$run()` is called repeatedly (samplewise) or once 
#'   providing a chunk. The reason is that the sequence of calling 
#'   the random number generator is different
#' @param out_file Optional. If given, specifies a file to save 
#'   the simulation \code{environment()}. Results can be loaded via 
#'   \code{readRDS(out_file)} afterwards.
#' @export
#--------------------------------------------------------------------
#' @examples
#' # This may take long because transfer functions have to be calculated, initially:
#' model <- calculate_semg_TFs(setup_muscles("config/integration.r", deterministic = TRUE), deterministic = TRUE) 
#' saveRDS(model, file="integration.rds")
#'
#' # Run simulation stepwise by repeatedly calling instance$run():
#' instance <- init.simulate_model_stepwise(model, deterministic=TRUE)
#' res <- apply(as.data.frame(model$load_vals$force), 1, instance$run)
#' 
#' # As an alternative return also the force and w/o summing (using sapply just for fun):
#' res <- sapply(as.data.frame(t(as.data.frame(model$load_vals$force))),
#'               function(x) {
#'                 instance$run(x)
#'                 return(c(emg=unlist(instance$getEmg()),force=instance$getForce()))   # note: emg and force are NOT summed across muscles
#'               })
#'         
#' # Run simulation on input array of force data by calling instance$run() once:
#' force_vals = as.data.frame(model$load_vals$force)[1:20,]
#' instance <- init.simulate_model_stepwise(model, deterministic=TRUE)
#' instance$run(force_vals)
#' 
#' # Performance tests:
#' microbenchmark::microbenchmark(instance$run(rnorm(2)))
#--------------------------------------------------------------------
init.simulate_model_stepwise = function(model, deterministic = FALSE, out_file = NULL) {
  message(paste0(whoami(2)," : Initializing ..."))
  if (deterministic)  # must be before initialize_sim_state()
    set.seed(1)
  
  emg <- NULL
  result <- NULL
  sim_state <- lapply(seq(1, nrow(model$MUs)),function(i) initialize_sim_state(num_electrodes = nrow(model$electrodes)))
  MU_firing_responses_fftvals_restructured <- restructure_firing_responses_fftvals(model$MU_firing_responses_fftvals)
  
  
  instance <- list(
    run = function(force_vals) {
      if (is.null(ncol(force_vals)))   # here, only one time step/sample. Generally, force has dimension [time x muscles]
        force_vals <- matrix(force_vals, ncol=length(model$muscles$muscle)) 
      excitation_vals <- matrix(sapply(1:ncol(force_vals), function(i) model$act_funcs_df[[i,'force_to_act_fun']](force_vals[,i])), nrow=nrow(force_vals))
      force_vals <- as.data.frame(matrix(sapply(model$MUs$muscle,function(x)force_vals[,x]),nrow=nrow(force_vals)))            # replicate muscle-specific force_vals for each motor unit to feed mapply below
      excitation_vals <- as.data.frame(matrix(sapply(model$MUs$muscle,function(x)excitation_vals[,x]),nrow=nrow(force_vals)))  # replicate muscle-specific excitation_vals for each motor unit to feed mapply below
      
      result <<- mapply(chunk_MU_contribs,   # apply simulation for a chunk of time steps across all motor units (MUs)
                        model$MUs$MU.obj,
                        sim_state,
                        force_vals,
                        excitation_vals,
                        model$MUs_rate_functions$MU_rate_function,
                        MU_firing_responses_fftvals_restructured$firing_response_fftvals,
                        MoreArgs = list(sampling = model$sampling),
                        SIMPLIFY = FALSE)
      
      sim_state <<- sapply(result, function(x)x$sim_state, simplify=FALSE)   # must be global to ensure that updated sim_state can be used in next call
      
      emg <<- sapply(model$muscles$muscle, function(m) {
        aux = sapply(result[which(model$MUs$muscle==m)], function(x)x$emg,simplify=FALSE)
        Reduce("+",aux)/length(aux)
      }, simplify=FALSE)
      
      return(Reduce("+",emg)/length(emg))   # return emg - averaged across muscles
    },
    
    getEmg = function() {   
      return(emg)          # return emg - *not* averaged across muscles
    },
    
    getForce = function() {
      force <- sapply(model$muscles$muscle, function(m) rowSums(matrix(sapply(result[which(model$MUs$muscle==m)], function(x)x$force), nrow=nrow(emg[[1]]))))
      return(force)        # return force - not averaged across muscle
    },
    
    getModel = function() {
      return(model)
    }
  )
  
  return(instance)
}



#--------------------------------------------------------------------
#' Run stepwise simulation
#' 
#' It is based on already existing `instance` for stepwise 
#' simulation generated by [init.simulate_model_stepwise()]. 
#' The function may not be useful per se because stepwise 
#' simulation usually avoids being called as in chunk based processing.
#' However, it allows direct comparison to the conventional method
#' [simulate_model()].
#' 
#' @details 
#' This simulation is equivalent to `simulate_model()` above.
#' Motor unit (MU) impulse trains and contribs are not saved. 
#' The examples below show how a simulation `instance` is created 
#' and simulation is run.
#' 
#' @param instance Produced by `run.simulate_model_stepwise()`.
#' @param muscle_force Either a list of muscle force functions 
#'   (one function per muscle) or a matrix with dimension
#'   \[time x muscles\] or a vector (only one muscle).
#' @param time_limits Optional. Two-elements vector limiting the 
#'   time axis. If `NULL` then `muscle_force` is assumed to be a 
#'   matrix or data frame which determines the time axis.
#' @export
#--------------------------------------------------------------------
#' @examples
#' # This may take long because transfer functions have to be calculated:
#' model <- calculate_semg_TFs(setup_muscles("config/integration.r", deterministic = TRUE), deterministic = TRUE) 
#' saveRDS(model, file="integration.rds")
#' 
#' # Now run the simulation conventionally:
#' model <- simulate_model(model, model$force_functions, model$time_limits, num_cores=1, parallel="chunk", deterministic=TRUE)
#' 
#' # Re-run stepwise simulation (as a chunk, though):
#' instance <- init.simulate_model_stepwise(model, deterministic=TRUE)
#' model <- run.simulate_model_stepwise(instance, model$force_functions, model$time_limits)
#'         
#' # Thus, the available semgsim plot routines can be used for displaying the results 
#' # (except MU specific plots because MU specific results were not saved):
#' plot_surface_potentials(model$surface_potentials, electrodes=c(1,2,3,4), electrode_confs = list(c(1, -1, 0, 0), c(0, 0, 1, -1)))
#' plot_muscle_contribs(model$muscle_contribs)
#'         
#' # Performance tests:
#' system.time(simulate_model(model, model$force_functions, model$time_limits, num_cores=1, parallel="sample", deterministic=TRUE))
#' system.time(run.simulate_model_stepwise(instance, model$force_functions, model$time_limits))
#--------------------------------------------------------------------
run.simulate_model_stepwise = function(instance, muscle_force = NULL, time_limits = NULL) {
  message(paste0(whoami(2)," : Running ...")) 
  
  model <- instance$getModel()

  if (is.null(muscle_force)) {    # use default muscle force functions provided in the model configuration
    muscle_force = model$force_functions  
  }
  
  if (is.list(muscle_force) &&    # muscle_force is a list of functions (as many as muscles), the force values are generated based on the time axis
      is.function(muscle_force[[1]]) && 
      !is.null(time_limits)) 
  {
    stopifnot(length(model$muscles$muscle) == length(muscle_force))
    t <- seq(time_limits[1], time_limits[2], 1 / model$Fs)
    force_vals <- sapply(model$muscles$muscle, function(forcefun_idx) muscle_force[[forcefun_idx]](t))
  } else if (is.numeric(muscle_force)) {  # muscle force is a vector or matrix [time x muscles] containing the force values which defines the time axis
    if (is.vector(muscle_force)) {        # if muscle force is a vector convert to matrix with 1 column
      muscle_force <- matrix(muscle_force, ncol=1)
    }
    stopifnot(length(model$muscles$muscle) == ncol(muscle_force))
    force_vals <- muscle_force
    time_limits <- c(0, (nrow(force_vals)-1) / model$Fs)    # define time axis range based on number of rows (samples) and sampling frequency
    t <- seq(time_limits[1], time_limits[2], 1 / model$Fs)  # construct the time axis
  } else {
    stop("muscle_force is neither a list of force functions nor a matrix of values")
  }
  
  res <- sapply(as.data.frame(t(force_vals)),
                function(x) {
                  instance$run(x)
                  return(c(emg=unlist(instance$getEmg()),force=instance$getForce()))   # note: emg and force are NOT summed across muscles
                })
  
  # Now, manually save the results into model (equivalent to what simulate_model() does, above):
  model$time_limits <- time_limits
  model$t <- seq(model$time_limits[1], model$time_limits[2], 1 / model$Fs)      # generate time axis
  emg <- as.data.frame(t(rbind(res[grep("^emg", row.names(res)),], t=model$t))) # transpose result and add time axis
  NM <- length(model$muscles$muscle); NE <- length(model$electrodes$electrode)  # shortcuts
  
  emgAux <- reshape(emg, varying=as.list(as.data.frame(matrix(1:(NE*NM),ncol=NM))), direction="long", timevar="electrode")   # reshape data.frame into long format - first melting electrodes
  emgAux <- reshape(emgAux, varying=list(grep("^emg",names(emgAux))), drop="id", direction="long", timevar="muscle")         # reshape a 2nd time melting muscles
  model$muscle_contribs$emg <- data.frame(muscle=emgAux$muscle, electrode=emgAux$electrode, time=emgAux$t, potential=emgAux[,grep("^emg", names(emgAux))])   # write back into model results
  
  force <- as.data.frame(t(rbind(res[grep("^force", row.names(res)),], t=model$t)))
  forceAux <- reshape(force, varying=list(grep("^force",names(force))), direction="long", timevar="muscle")   # reshape data.frame into long format melting muscles
  model$muscle_contribs$force <- data.frame(muscle=forceAux$muscle, time=forceAux$t, force=forceAux[,grep("^force", names(forceAux))])   # write back into model results
  
  model$surface_potentials <- sum_muscle_contribs(muscle_firing_contribs = model$muscle_contribs$emg)
  model$sim_results <- aggregate_sim_results(model$muscles, model$muscle_contribs, model$surface_potentials)
  
  return(model)  
}



#--------------------------------------------------------------------
#' Calibrate force values
#' 
#' The amplitude of the simulated force is not necessarily similar
#' to the amplitude of the force target provided by the parameter
#' muscle force. While the target is in the range between 0 and 1
#' the amplitude of the resulting simulation is not known beforehand.
#' Thus, a calibration is needed to determine the gain.
#' 
#' @param model Produced by `calculate_semg_TFs()`.
#' @param deterministic Option sets seed for reproducible random 
#'   numbers.
#' @param muscle_force Either a list of muscle force functions 
#'   (one function per muscle) or a matrix with dimension
#'   \[time x muscles\] or a vector (only one muscle).
#' @param time_limits Optional. Two-elements vector limiting the 
#'   time axis. If `NULL` then `muscle_force` is assumed to be a 
#'   matrix or data frame which determines the time axis.
#' @param delay Rough estimate of delay between force target and
#'   simulated force. 
#' @export
#--------------------------------------------------------------------
#' @examples
#' model = readRDS("integration.rds")
#' model = calibrate_model_force(model)
#' cat(model$muscles$force_gain)
#--------------------------------------------------------------------
calibrate_model_force = function(model, time_limits=c(0,20),
                                 muscle_force = sapply(row.names(model$muscles),function(x)   # for each muscle define the same default force function
                                 function(t) {                   # default force function used to perform force calibration
                                   dT = 1/model$Fs
                                   T = 4                         # period: 4 seconds
                                   aux = function(t) {
                                     t = t %% T
                                     if (t < 0.25*T) {           # 1. second: zeros
                                       return(0)
                                     } else if (t < 0.5*T) {     # 2. second: linear rise from zero to one
                                       return(4/T*t - 1)
                                     } else if (t < 0.75*T) {    # 3. second: ones
                                       return(1)
                                     } else {                    # 4. second: linear decay from one to zero
                                       return(-4/T*t + 4)
                                     }
                                   }
                                   sapply(t,aux)
                                 }),
                                 deterministic = FALSE,
                                 delay = 0.06                    # rough estimate of delay between force target and resulting force
                                 ) 
{
  instance <- init.simulate_model_stepwise(model, deterministic=deterministic)

  if (is.null(muscle_force)) {    # use default muscle force functions provided in the model configuration
    muscle_force = model$force_functions  
  }
  
  if (is.list(muscle_force) &&    # muscle_force is a list of functions (as many as muscles), the force values are generated based on the time axis
      is.function(muscle_force[[1]]) && 
      !is.null(time_limits)) 
  {
    stopifnot(length(model$muscles$muscle) == length(muscle_force))
    t <- seq(time_limits[1], time_limits[2], 1 / model$Fs)
    force_vals <- sapply(model$muscles$muscle, function(forcefun_idx) muscle_force[[forcefun_idx]](t))
  } else if (is.numeric(muscle_force)) {  # muscle force is a vector or matrix [time x muscles] containing the force values which defines the time axis
    if (is.vector(muscle_force)) {        # if muscle force is a vector convert to matrix with 1 column
      muscle_force <- matrix(muscle_force, ncol=1)
    }
    stopifnot(length(model$muscles$muscle) == ncol(muscle_force))
    force_vals <- muscle_force
  } else {
    stop("muscle_force is neither a list of force functions nor a matrix of values")
  }
  
  message(paste0(whoami(2)," : Running ...")) 
  res <- sapply(as.data.frame(t(force_vals)),
                function(x) {
                  instance$run(x)
                  return(instance$getForce())   
                })
  
  n = floor(delay*model$Fs)    # determine number of samples to delay
  
  model$muscles$force_gain = c()
  for (i in model$muscles$muscle) {
    force_vals[,i] = c(rep(force_vals[1,i],n), head(force_vals[,i], -n))   # apply delay to force_vals and do zero order extrapolation at the beginning
    x = force_vals[,i]
    y = res[i,]
    fit = lm(y ~ 0 + x)
    model$muscles$force_gain[i] = coefficients(fit)[[1]]   # store gain factor in model
  }
     
  return(model)
}

