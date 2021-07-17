resp_insp_load <- function(t, breath_length, start_load, max_load, I_E_ratio, transient_length,
                           rel_length_insp_activity_in_exp = 1/5) {

    insp_length <- I_E_ratio * breath_length / (I_E_ratio + 1)
    exp_length <- breath_length - insp_length
    exp_activity_length <- rel_length_insp_activity_in_exp * exp_length
    
    single_breath_load <- function(t) {
        transient_mask <- t < transient_length
        insp_mask <- (t >= transient_length) & (t < insp_length)
        exp_start_mask <-
            (t >= insp_length) & (t < insp_length + exp_activity_length)
        load <- rep(0, length(t))
        load[transient_mask] <- t[transient_mask] / transient_length * start_load
        load[insp_mask] <-
            start_load + t[insp_mask] / insp_length * (max_load - start_load)
        load[exp_start_mask] <- max_load *
          (1 - (t[exp_start_mask] - insp_length) / exp_activity_length)
        load
    }

    single_breath_load(t %% breath_length)
}


resp_exp_load <- function(t, breath_length, start_load, max_load, I_E_ratio,
                          rel_length_to_max_load) {

    insp_length <- I_E_ratio * breath_length / (I_E_ratio + 1)
    exp_length <- breath_length - insp_length
    exp_rising_length <- exp_length * rel_length_to_max_load
    exp_falling_length <- exp_length - exp_rising_length
    turning_point <- insp_length + exp_rising_length
    
    single_breath_load <- function(t) {
        exp_rising_mask <- (t >= insp_length) & (t < turning_point)
        exp_falling_mask <- t >= turning_point
        load <- rep(0, length(t))
        load[exp_rising_mask] <- start_load + (max_load - start_load) * 
          (t[exp_rising_mask] - insp_length) / exp_rising_length
        load[exp_falling_mask] <- max_load *
          (1 - (t[exp_falling_mask] - turning_point) / exp_falling_length)
        load
    }

    single_breath_load(t %% breath_length)
}


sinusoid_activity_load <- function(t, t_start, length, start_load, max_load) {

    load <- rep(0, length(t))
    activity_mask <- (t >= t_start) & (t <= t_start + length)
    load[activity_mask] <- start_load +
      sin((t[activity_mask] - t_start) / length * pi)^2 * (max_load - start_load)
    load
}


cough_load <- function(t, t_start, cough_length, n_coughs) {

    load <- rep(0, length(t))
    t_start <- t_start + (0:(n_coughs - 1)) * 2 * cough_length
    for (i in 1:n_coughs)
        load <- load + 1 / i * sinusoid_activity_load(t,
                                                      t_start = t_start[i],
                                                      length = cough_length,
                                                      start_load = 0.6,
                                                      max_load = 0.9)
    load
}


stair_load <- function(t, t_start, start_load, max_load, load_step, step_length, transient_length, gap_length=0) {

	num_steps <- round((max_load - start_load) / load_step) + 1

  load <- rep(0, length(t))
	
	for (i in 1:num_steps) {
	
		step_mask <- (t >= t_start + (i-1) * step_length + transient_length) & (t <= t_start + i * step_length - gap_length)
		transient_mask <- (t >= t_start + (i-1) * step_length) & (t < t_start + (i-1) * step_length + transient_length)
		if (gap_length == 0)
		  load[transient_mask] <- start_load + (i-2) * load_step + 
		    (t[transient_mask] - (t_start + (i-1) * step_length)) / transient_length * load_step
		else
		  load[transient_mask] <- (t[transient_mask] - (t_start + (i-1) * step_length)) / transient_length * (
		    start_load + (i-1) * load_step)
		
		load[step_mask] <- start_load + (i-1) * load_step
	}
	
    load
}