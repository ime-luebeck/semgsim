#' Calculate nonlinear transformations of activations signals
#'
#' Calculate a nonlinear input transformation for each muscle that makes the 
#' muscle generated the desired relative output force.
#'
#' Depending on the firing rate, MU pool, and twitch force model, a nonlinear
#' transformation of the muscle activation input is calculated, such that the
#' generated relative output muscle force matches the input activation signal in 
#' the mean over time.
#' 
#' @param MUs A data.frame of MUs belonging to (potentially) multiple muscles.
#' @param MUs_rate_functions A data frame containing a column $rate.function that gives the
#'   firing rate as a function of muscle activation for each MU.
#' @return A list of functions that gives the muscle input as a function of the desired
#'   relative output force for each muscle.
#'
calc_activation_transformations <- function(MUs, MUs_rate_functions) {
  
  # Determine interpolating grid
  acts = seq(0, 1, 0.01)
  
  calc_mean_muscle_force <- function(muscle_id, act) {
    force <- 0
    # Loop over MUs and accumulate the generated for contributions
    muscle_MUs <- dplyr::filter(MUs, muscle == muscle_id)
    for (row_idx in seq(1, nrow(muscle_MUs))) {
      MU_row <- muscle_MUs[row_idx, ]
      fr_func <- dplyr::filter(MUs_rate_functions, muscle == muscle_id, MU == MU_row$MU)$MU_rate_function[[1]]
      fr <- fr_func(act)
      MU_obj <- MU_row$MU.obj[[1]]
      force <- force + MU_obj$total_twitch_power * fr * MU_obj$gain_fac_fun(1/fr)
    }
    
    stopifnot(force >= 0)
    force
  }

  # Calculate the desired nonlinear transformations for each muscle
  act_funcs_lst <- lapply(seq_along(unique(MUs$muscle)), function(muscle_id) {
    
    # Loop over activation levels and calculated the mean relative generated force
    forces <- lapply(acts, function(act) calc_mean_muscle_force(muscle_id, act)) %>% simplify2array
    
    # Normalize
    forces_rel <- forces / max(forces)
  
    # Generate and return interpolated forward and inverse relationships
    forward_interp <- suppressWarnings(splinefun(x = acts, y = forces_rel, method="monoH.FC"))
    inverse_interp <- suppressWarnings(splinefun(x = forces_rel, y = acts, method = "monoH.FC"))
    
    if (inverse_interp(0) < 0 || forward_interp(0) < 0) {
      stop()
    }
    
    list(forward_interp, inverse_interp)
  })
  
  
  act_funcs_df <- data.frame(muscle = seq_along(unique(MUs$muscle)))
  act_funcs_df$act_to_force_fun <- lapply(act_funcs_lst, function(lst) lst[[1]])
  act_funcs_df$force_to_act_fun <- lapply(act_funcs_lst, function(lst) lst[[2]])
  
  act_funcs_df
}


update_MU_rec_threshs <- function(MUs, act_funcs_df) {

  MUs$MU.obj <- plyr::alply(MUs, .margins = 1, .fun = function(df_row) {
    act_to_force_fun <- dplyr::filter(act_funcs_df, muscle == df_row$muscle)$act_to_force_fun[[1]]
    MU_obj <- df_row$MU.obj[[1]]
    MU_obj$rec_thresh_force <- act_to_force_fun(MU_obj$rec_thresh)
    MU_obj
  })
  
  MUs
}