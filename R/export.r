#' @export
export_to_mat <- function(model, out_file, export_individual_impulse_trains=FALSE,
                          export_summed_impulse_trains=FALSE, export_emg_contribs=FALSE) {

	export_data <- setup_export_data(model, export_individual_impulse_trains, export_summed_impulse_trains, export_emg_contribs)
    
	# This ugly if-else-branching should be done nicer.
	if (export_emg_contribs) {
		if (export_individual_impulse_trains) {
		  if (export_summed_impulse_trains)
			R.matlab::writeMat(out_file,
							   emg_measurements = export_data$emg_measurements,
							   force_contributions = export_data$force_contributions,
							   force_profiles = export_data$force_profiles,
							   individual_impulse_trains = export_data$individual_impulse_trains,
							   summed_impulse_trains = export_data$summed_impulse_trains,
							   emg_contribs = export_data$emg_contribs,
							   time = export_data$time)
		  else
			R.matlab::writeMat(out_file,
							   emg_measurements = export_data$emg_measurements,
							   force_contributions = export_data$force_contributions,
							   force_profiles = export_data$force_profiles,
							   individual_impulse_trains = export_data$individual_impulse_trains,
							   emg_contribs = export_data$emg_contribs,
							   time = export_data$time)
		} else if (export_summed_impulse_trains) {
		  R.matlab::writeMat(out_file,
							 emg_measurements = export_data$emg_measurements,
							 force_contributions = export_data$force_contributions,
							 force_profiles = export_data$force_profiles,
							 summed_impulse_trains = export_data$summed_impulse_trains,
							 emg_contribs = export_data$emg_contribs,
							 time = export_data$time)
		
	  } else
		R.matlab::writeMat(out_file,
						   emg_measurements = export_data$emg_measurements,
						   force_contributions = export_data$force_contributions,
						   force_profiles = export_data$force_profiles,
						   emg_contribs = export_data$emg_contribs,
						   time = export_data$time)
	} else {
		if (export_individual_impulse_trains) {
		  if (export_summed_impulse_trains)
			R.matlab::writeMat(out_file,
							   emg_measurements = export_data$emg_measurements,
							   force_contributions = export_data$force_contributions,
							   force_profiles = export_data$force_profiles,
							   individual_impulse_trains = export_data$individual_impulse_trains,
							   summed_impulse_trains = export_data$summed_impulse_trains,
							   time = export_data$time)
		  else
			R.matlab::writeMat(out_file,
							   emg_measurements = export_data$emg_measurements,
							   force_contributions = export_data$force_contributions,
							   force_profiles = export_data$force_profiles,
							   individual_impulse_trains = export_data$individual_impulse_trains,
							   time = export_data$time)
		} else if (export_summed_impulse_trains) {
		  R.matlab::writeMat(out_file,
							 emg_measurements = export_data$emg_measurements,
							 force_contributions = export_data$force_contributions,
							 force_profiles = export_data$force_profiles,
							 summed_impulse_trains = export_data$summed_impulse_trains,
							 time = export_data$time)
		
	  } else
		R.matlab::writeMat(out_file,
						   emg_measurements = export_data$emg_measurements,
						   force_contributions = export_data$force_contributions,
						   force_profiles = export_data$force_profiles,
						   time = export_data$time)
	}						   
}


#' @export
export_to_csv <- function(model, out_file, export_individual_impulse_trains=FALSE, export_summed_impulse_trains=FALSE, export_emg_contribs=FALSE) {

	if (export_individual_impulse_trains || export_summed_impulse_trains || export_emg_contribs)
		error("These export options are currently only supported for .mat files, not for .csv. This could be easily changed, though.")

	export_data <- setup_export_data(model)
	
	emg_measurements_df <- data.frame(t(export_data$emg_measurements))
	for (electrode_id in 1:nrow(export_data$emg_measurements))
		colnames(emg_measurements_df)[electrode_id] = sprintf("EMG %d", electrode_id)
		
	force_contributions_df <- data.frame(t(export_data$force_contributions))
	for (muscle_id in 1:nrow(export_data$force_contributions))
		colnames(force_contributions_df)[muscle_id] = sprintf("Muscle force %d", muscle_id)
	
	force_profiles_df <- data.frame(t(export_data$force_profiles))
	for (muscle_id in 1:nrow(export_data$force_profiles))
		colnames(force_profiles_df)[muscle_id] = sprintf("Muscle force profile %d", muscle_id)
		
	time_df <- data.frame(export_data$time)
	colnames(time_df)[1] <- "time"
	
	export_df <- cbind(time_df, emg_measurements_df, force_contributions_df, force_profiles_df)
	
	readr::write_csv(export_df, out_file)
}


setup_export_data <- function(results, export_individual_impulse_trains=FALSE, export_summed_impulse_trains=FALSE, export_emg_contribs=FALSE) {
												  
  muscle_ids <- 
    results$muscle_contribs$emg$muscle %>% as.factor %>% levels %>% as.integer
  electrode_ids <- 
    results$surface_potentials$electrode %>% as.factor %>% levels %>% as.integer
  time <- dplyr::filter(results$muscle_contribs$emg,
                        muscle == muscle_ids[1], electrode == electrode_ids[1])$time

  ## Muscle activity profiles
  force_profiles <- array(dim=c(length(muscle_ids), length(time)),
                               dimnames=list(Muscle=NULL, Time=NULL))
  for (muscle_id in muscle_ids)
    force_profiles[muscle_id, ] <- results$force_functions[[muscle_id]](time)
  
  ## EMG signals
  emg_measurements <-
    array(dim=c(length(electrode_ids), length(time)),
          dimnames=list(Electrode=NULL, Time=NULL))  
  
  for (electrode_id in electrode_ids)
      emg_measurements[electrode_id, ] <- dplyr::filter(results$surface_potentials,
                                                        electrode == electrode_id)$potential
  
  stopifnot(all(dplyr::filter(results$surface_potentials,
                       electrode == electrode_ids[length(electrode_ids)])$potential[1:10]
                == emg_measurements[length(electrode_ids), 1:10]))
  
  ## Force signals
  force_contributions <- array(dim=c(length(muscle_ids), length(time)),
                               dimnames=list(Muscle=NULL, Time=NULL))

  for (muscle_id in muscle_ids)
    force_contributions[muscle_id, ] <- dplyr::filter(results$muscle_contribs$force,
                                                      muscle == muscle_id)$force
  
  export_data <- list()
  export_data$time <- time
  export_data$emg_measurements <- emg_measurements
  export_data$force_contributions <- force_contributions
  export_data$force_profiles <- force_profiles
  
  
  ## Impulse trains
  if (export_summed_impulse_trains || export_individual_impulse_trains) {
    
    if (export_summed_impulse_trains) {
      summed_impulse_trains <- array(dim=c(length(muscle_ids), length(time)),
                                     dimnames=list(Muscle=NULL, Time=NULL))
      summed_impulse_trains[,] <- 0
    }
    
    if (export_individual_impulse_trains) {
      
      max_num_MUs <- 0
      for (muscle in results$muscles$muscle.obj)
        max_num_MUs <- max(max_num_MUs, muscle$num_MUs)
      
      individual_impulse_trains <- array(dim=c(length(muscle_ids), 
                                               max_num_MUs,
                                               length(time)))
      individual_impulse_trains[,,] <- 0
    }
    
    for (muscle_id in muscle_ids) {
      muscle_impulse_trains <- dplyr::filter(results$MU_impulse_trains, muscle == muscle_id)
      kk <- 1
      for (MU_impulse_train in muscle_impulse_trains$impulse_train) {
        if (export_individual_impulse_trains) {
          individual_impulse_trains[muscle_id, kk,] <- MU_impulse_train
          kk <- kk + 1
        }
        if (export_summed_impulse_trains)
          summed_impulse_trains[muscle_id,] <- 
            summed_impulse_trains[muscle_id,] + MU_impulse_train
      }
    }
	
    if (export_individual_impulse_trains)
	  export_data$individual_impulse_trains <- individual_impulse_trains
    if (export_summed_impulse_trains)
	  export_data$summed_impulse_trains <- summed_impulse_trains
  }
  
  ## EMG_contribs
  if (export_emg_contribs) {
		emg_contribs <- array(dim=c(length(muscle_ids), length(electrode_ids), length(time)),
								   dimnames=list(Muscle=NULL, Time=NULL))			   
		for (muscle_id in muscle_ids) {
			for (electrode_id in electrode_ids) {
				emg_contribs[muscle_id, electrode_id,] <- dplyr::filter(results$muscle_contribs$emg,
                                                      muscle == muscle_id & electrode == electrode_id)$potential
			}
		}
		export_data$emg_contribs <- emg_contribs
  }
  
  return(export_data)
}