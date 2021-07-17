#' Calculate the firing-rate-force relationship of all MUs
#'
#' For each MU calculate a function that gives its firing rate in dependence of
#' the current muscle work load.
#'
#' Note that this includes the different MU's activation threshold, as before
#' being activated the firing rate of a MU is just zero. The muscle work load
#' should be given in the ratio of the maximum voluntary contraction (MVC) of
#' the muscle, i.e., a value between 0 and 1.
#'
#' @param muscles A data frame containing information on the muscles.
#' @param MUs A data frame containing MU objects.
#'
#' @return A data frame containing a column $rate.function that gives the
#'   function described above for each MU.
#' 
calc_MU_firing_rate_funcs <- function(muscles, MUs) { 

    rate_function_creators <- create_rate_function_creators(muscles)

    MU_rate_functions <- merge(MUs, rate_function_creators)

    ## (Unnecessarily) return list of lists since do.call (see below) requires 
    ## a list of arguments.
    rec_threshs <- lapply(MU_rate_functions$MU.obj,
                          function(MU) list(MU$rec_thresh))
    
    MU_rate_functions$MU_rate_function <-
        mapply(do.call,
               MU_rate_functions$rate_function_creator,
               rec_threshs,
               SIMPLIFY = FALSE)

    MU_rate_functions[, c("muscle", "MU", "MU_rate_function")]
}

#' Set up rate coding for a set of muscles
#'
#' This wraps \link{create_rate_function_creator} by applying it to a data frame
#' of muscles and coercing the output into a data frame again.
#'
#' @param muscles A data frame containing the muscle objects for which the rate
#'   coding shall be set up.
#' 
#' @return A data frame containing for each muscle a function that takes the
#'   recruitment threshold of a MU in that muscle and returns a function giving
#'   the firing rate of that MU in dependence of the current muscle work load.
#' 
create_rate_function_creators <- function(muscles) {

    unique_idces <- get_unique_idces(muscles$muscle)

    rate_function_creators <- unique_idces %>%
        lapply(function(idx) {
            muscle_obj <- muscles$muscle.obj[[idx]]
            rate_function_creator <-
                create_rate_function_creator(muscle_obj$rate_coding_params)
            rate_function_creator
        })

    data.frame(muscle = muscles[unlist(unique_idces), "muscle"],
               rate_function_creator = I(rate_function_creators))
}

#' Set up the firing-rate-force relationships, i.e. the rate coding, for a
#' muscle
#'
#' Create a factory for firing-rate-force relationships of the MUs in the
#' muscle.
#'
#' @param params A list of values that parameterize the rate coding
#'   characteristics of a muscle. Must contain parameters "A" to "E" or "A" to "G" 
#'   if params$alternative_rate_coding is TRUE.
#' @return A function that takes the recruitment threshold of a MU and returns a
#'   function that gives that MU's firing rate in dependence of the current
#'   muscle work load.
#' 
create_rate_function_creator <- function(params) {

	if (!params$alternative_rate_coding) {
		#' Model of De Luca and Hostage (2010)
		with(params, 
			 stopifnot(A >= 0, B >= 0, C <= 0, D >= 0, E >= 0,
					   all(simplify2array(
						   lapply(list(A, B, C, D, E), is.nice.scalar))))) 
	
		rate_function_creator <- function(recruitment_threshold, params) {
			params$recruitment_threshold <- recruitment_threshold
			functional::Curry(rate_function, params = params)
		}
    } else {
		#' New model proposed in our frontiers in physiology modelling paper
		with(params, 
			 stopifnot(A >= 0, B >= 0, C >= 0, D >= 0, E >= 0, F >= 0, G >= 0,
					   all(simplify2array(
						   lapply(list(A, B, C, D, E, F, G), is.nice.scalar))))) 
	
		rate_function_creator <- function(recruitment_threshold, params) {
			params$recruitment_threshold <- recruitment_threshold
			functional::Curry(alternative_rate_function, params = params)
		}
	}
    functional::Curry(rate_function_creator, params = params)
}


#' This implements an exponential model of the firing-rate-force relationship
#' relationship, as proposed in "Relationship between firing rate and
#' recruitment threshold of motoneurons in voluntary isometric contractions"
#' [De Luca, Hostage] (2010). The model is
#'
#' fr_i(x) = Dx + (C - A exp(-x/B)) x_i + E
#'
#' where fr_i(x) denotes the firing rate of motoneuron i at excitation level x,
#' and x_i denotes the recruitment threshold of that motoneuron. The remaining
#' parameters are constant parameters of the rate coding model that describe the
#' rate coding behavior of a particular muscle.
#'
#' The model respects the following principles: \itemize{
#' \item The 'onion-skin' phenomenon: at each point in time, earlier recruited
#'   MUs maintain higher firing rates than later recruited ones.
#' \item Firing-rate curves converge to similar (here: the same) firing rates at
#'   high force levels.
#' \item At a given value of force, higher-threshold MUs display higher slopes
#'   in their firing-rate-force curves.
#' \item Initial firing rates of MUs increase with the recruitment threshold.
#' }
rate_function <- function(muscle_load, params) {

    with(params, {
        recruited_mask <- muscle_load >= recruitment_threshold
        
        rate <- rep(0, length(muscle_load))
        rate[recruited_mask] <- D * muscle_load[recruited_mask] +
            (C - A * exp(-muscle_load[recruited_mask] / B)) * recruitment_threshold + E

        rate
    })
}


#' This implements the rate coding model proposed in 
#' Petersen and Rostalski (2019): A Comprehensive Mathematical Model of Motor Unit Pool Organization, Surface
#' Electromyography, and Force Generation. frontiers in Physiology. https://doi.org/10.3389/fphys.2019.00176
#'
alternative_rate_function <- function(muscle_load, params) {

    with(params, {
        recruited_mask <- muscle_load >= recruitment_threshold
        
        rate <- rep(0, length(muscle_load))
    		rate[recruited_mask] <- -A * (B - muscle_load[recruited_mask]) * recruitment_threshold + 
    			C * muscle_load[recruited_mask] + D - 
    		  (E - F*recruitment_threshold) * exp(-(muscle_load[recruited_mask] - recruitment_threshold)/G)
        
    		rate
    })
}