#' @importFrom plyr "."
sum_MU_contribs <- function(MU_firing_contribs, time, sampling) {

    emg_dfs <- list()
	force_dfs <- list()
    i <- 1
	j <- 1
	
    for (muscle_id in unique(MU_firing_contribs$muscle)) {
        muscle_reduced_df <- MU_firing_contribs %>% dplyr::filter(muscle == muscle_id)
		
		## FORCE
		force_dfs[[j]] <- data.frame(muscle = muscle_id,
							   time = time,
							   force = rep(0, length(t)))
        ## loop over MUs and add their respective force contributions
		for (row_id in as.numeric(rownames(muscle_reduced_df))) {
				force_dfs[[j]]$force <- force_dfs[[j]]$force +
					muscle_reduced_df$contribs[[row_id]]$force
		}
		
		j <- j + 1
		
		## EMG
        for (electrode_id in seq(1, length(muscle_reduced_df$contribs[[1]]$emg))) {
            emg_dfs[[i]] <- data.frame(muscle = muscle_id,
                                   electrode = electrode_id,
                                   time = time,
                                   potential = rep(0, length(t)))
            ## Now loop over MUs and add their respective contributions
            for (row_id in as.numeric(rownames(muscle_reduced_df)))
				emg_dfs[[i]]$potential <- emg_dfs[[i]]$potential +
					muscle_reduced_df$contribs[[row_id]]$emg[[electrode_id]]
            
            i <- i + 1
        }
    }

	muscle_contribs <- list()

    muscle_contribs$emg <- plyr::rbind.fill(emg_dfs)

	muscle_contribs$force <- plyr::rbind.fill(force_dfs)

    muscle_contribs
}

#' Combine the contributions of all individual MUs to each electrode to obtain the final electrode signals.
#'
sum_muscle_contribs <- function(muscle_firing_contribs) {

    electrode_potentials <-
        muscle_firing_contribs %>%
            dplyr::group_by(electrode, time) %>%
                dplyr::summarise(potential = sum(potential))
        
    electrode_potentials[, c("electrode", "time", "potential")]
}
