##### Rectus Abdominis Superior Definitions #####

generate_rectus_part <- function(width, thickness, length, base_pos, emg_force_rel) {
    ## Create the muscle object.
    ## Set the base shape of the muscle endplate.
    ## Currently supported shapes are:
    ##    - "rect". Will be interpreted as the unit square [0,1]x[0,1].

    rectus_part <- create_rect_muscle(width = width,
                                      thickness = thickness)
    
    
    ## Define the geometrical transformation by which the endplate base shape will
    ## be turned into its final shape. Must be some function that takes a 2D vector
    ## and returns a 3D vector, as this is also used to position the endplate in 3D
    ## space.

    rectus_part$shape_transformation <-
        function(x) c(x[1] * width + base_pos[1],
                      -x[2] * thickness,
                      base_pos[2])


    ## Define the muscle fiber shape. This must be an object providing a method
    ## $pos(x) that implements a parametrization of the fiber shape curve. It takes
    ## a numerical argument x between 0 and 1 denoting the relative position on the
    ## fiber and returns a 3D vector containing the position of the curve in 3D
    ## space at that length. The curve is assumed to start at (0,0) and to be 
    ## normalized to length 1, i.e. the length of the curve given by pos(x) from x=0
    ## to x=1 is equal to one.
    ## Here, we set all fibers to be straight lines.

    rectus_part$fiber_shape <- create_straight_line()


    ## Define the rotation that will be applied to the muscle fiber's base shape,
    ## depending on their relative position on the muscle base shape.
    ## The absolute position in 3D space of a point on the fiber is then given in
    ## terms of its relative position on the fiber (which is between 0 and 1) as
    ##
    ##     x(x_rel) = fiber$innervation_zone + fiber$length * translation,
    ##
    ## where
    ##
    ##     translation =
    ##        rotate3D_lh(fiber_shape$pos(x.rel) -
    ##                       fiber_shape$pos(innervation_zone_rel),
    ##                    fiber_rotations(fiber_base_pos_rel_to_muscle))
    ##
    ## and rotate3D_lh(x, angles) rotates the vector x by the three angles specified
    ## in angles around the x, y and z axis, respectively. Positive angles result in
    ## counter-clockwise rotations.

    rectus_part$fiber_rotations <- function(pos_rel_to_muscle)
        c(0, 0, 0)

    ## Define how the muscle fibers will be placed inside a given MU.

    rectus_part$rel_fiber_distributor <- distribute_fibers_unif


    ## Define the base length of the muscles fibers as a function of their position
    ## on the muscle endplate.

    rectus_part$fiber_base_length_calculator <- function(base_pos_rel_to_muscle) length

    ## Define a random relative variation of the fiber length.
    ## Don't forget to keep the relation between this and fiber_length_max (see
    ## below) up to date.

    rectus_part$variation <- 0.025

    rectus_part$fiber_length_randomizer <-
        function() runif(n = 1,
                         min = 1 - rectus_part$variation,
                         max = 1 + rectus_part$variation)


    ## Give an upper bound for the value
    ##
    ##    fiber_base_length_calculator() * fiber_length_randomizer().
    ##
    ## This will be used for estimating an upper bound on the time scale length of
    ## an IAP impulse response as measured on the skin.

    rectus_part$fiber_length_max <-
        length * (1 + rectus_part$variation)


    ## Define, how the relative locations of the innervation zones of the muscle
    ## fibers will be determined.

    rectus_part$innervation_zone_placer <-
        function(fiber) rnorm_bounded(min = 0.4, max = 0.6, n_sigma = 3)


    ## Give an upper bound on
    ##
    ##     max(innervation_zone_placer(fiber), 1 - innervation_zone_placer(fiber)).
    ##
    ## This will be used for estimating an upper bound on the time scale length of
    ## an IAP impulse response as measured on the skin.

    rectus_part$fiber_rel_half_length_max <- 0.6


    ## Define whether the conduction velocity of the muscle fibers can change with
    ## time.
    ## Here, we take it to be constant.

    rectus_part$cv_variability <- "const"

    rectus_part$MU_class <- "MU_circ"


    ## MU Pool Organization ####

    rectus_part$num_MUs <- 50	
	
    rectus_part$random_MU_params <- list(
                      fiber_density = list(shape = 8, scale = 9 / mm_^2,
                                           min = 1 / mm_^2, max = 30 / mm_^2))

    ## Rate Coding and Recruitment ##
    ## Define the parameters that specify the characteristics of the rate coding
    ## strategy employed in the muscle.
    rectus_part$rate_coding_params <- list()
	if(exists('alternative_rate_coding') && alternative_rate_coding) {
		# Use our newly proposed rate coding model
		rectus_part$rate_coding_params$alternative_rate_coding <- TRUE
		rectus_part$rate_coding_params$A <- 20
		rectus_part$rate_coding_params$B <- 1.5
		rectus_part$rate_coding_params$C <- 30
		rectus_part$rate_coding_params$D <- 13
		rectus_part$rate_coding_params$E <- 8
		rectus_part$rate_coding_params$F <- 8
		rectus_part$rate_coding_params$G <- 0.05
	} else {
		# Use the model of De Luca and Hostage (2010)
		rectus_part$rate_coding_params$alternative_rate_coding <- FALSE
		
		## Define the minimum and maximum firing rate (in pps, pulses per second) of all
		## MUs in the muscle.
		rectus_part$rate_coding_params$A <- 100
		rectus_part$rate_coding_params$B <- 0.2
		rectus_part$rate_coding_params$C <- -25
		rectus_part$rate_coding_params$D <- 7
		rectus_part$rate_coding_params$E <- 20
	}
	
	## Define the ratio of maximum muscle force at which MU recruitment is complete.
	## Must be a value between 0 and 1.
	rec_stop <- 0.6	
	
	rectus_part$rate_coding_params$rec_model <- 2
	
	## Shape parameter of the recruitment model
	rectus_part$rate_coding_params$shape_b <- 0.2
	
	a1 <- log(rec_stop) / rectus_part$num_MUs
	a2 <- log(rec_stop/ rectus_part$rate_coding_params$shape_b) / rectus_part$num_MUs

	if (rectus_part$rate_coding_params$rec_model == 1) {
	  rec_start <- exp(a1 * 1)
	} else 
	  rec_start <- rectus_part$rate_coding_params$shape_b * 1 * exp(a2 * 1) / rectus_part$num_MUs

	## Size principle parameters
	LINEAR <- 1   # Variable is linearly related to the recruitment threshold
	SQUARE <- 2   # Variable is related to the square of the recruitment threshold
	
	# Intracellular conductivity, cf. Dimitrova (1974) http://www.ncbi.nlm.nih.gov/pubmed/4457323
	cond_ic <- 1

	# Range of fiber radius, cf. Hwang et al. (2013) [for index finger!]
	# https://doi.org/10.3109/2000656X.2012.755988
	min_fib_rad <- 12e-6
	max_fib_rad <- 20e-6
	
	rectus_part$MU_size_params <- list(rec_thresh = list(min=rec_start, max=rec_stop, lower_bound=0.005,
														 relation=LINEAR, round=FALSE),
									   peak_force = list(min=1, max=100, lower_bound=0.1,
														 relation=LINEAR, round=FALSE),
									   cv = list(min=2.5, max=5.5, lower_bound=2,
												 relation=LINEAR, round=FALSE),
									   twitch_t_rise = list(shape = 7, min = 50 * ms_, max = 30 * ms_, 
															lower_bound = 20 * ms_, relation=LINEAR, round=FALSE),
									   twitch_t_half_rel = list(shape = 7, min = 125 * ms_, max = 25 * ms_,
																lower_bound = 20 * ms_, relation=LINEAR, round=FALSE),
									   twitch_t_lead = list(shape = 10, min = 10 * ms_, max = 3.5 * ms_, lower_bound = 3 * ms_,
															relation=LINEAR, round=FALSE),
										IAP_scale_amp = list(shape = 12,
														min = cond_ic*pi*min_fib_rad^2,
														max = cond_ic*pi*max_fib_rad^2,
														lower_bound = 0.8*cond_ic*pi*min_fib_rad^2,
														relation=SQUARE,
														round=FALSE),
									   peak_fiber_force = list(shape = 8, min = 0.033, max = 0.1, lower_bound = 0.025,
															   relation=SQUARE, round=FALSE))
	
	
	
    ## Nonlinear force gain factor parameters
    rectus_part$twitch_r <- 0.85
    rectus_part$twitch_c <- 2.5

	
	## Define the way in which the motor unit (MU) regions are determined from the
	## muscle endplate region.	
    rectus_part$MU_distributor <-
        function(MUs) distribute_MUs_unif(MUs = MUs,
                                          muscle = rectus_part,
                                          tiling = c(3, 2))
    
    rectus_part
}

## For ref, refer to Rankin, Stokes, Newham / Abdominal Muscle Size and Symmetry in normal Subjects, 2006, Muscle & Nerve
## Right RA thickness 1.25, cross/sectional area 8.27
## Left RA thickness 1.24, cross/sectional area 8.15
## both for mean male subjects
## Teyhen, Chils et al. 2012 report mean over both left and right /male again/ to be 1.41, 8.00 which is in good agreement.
## Delp et al. 2001 report cumulative recti lengths of mean 34.3cm.
left_thickness = 1.24 * cm_
right_thickness = 1.25 * cm_
left_width = 8.15 / 1.24 * cm_
right_width = 8.27 / 1.25 * cm_
lower_length = 10 * cm_
mid_length = 12 * cm_
upper_length = 10 * cm_

rectus_la <- generate_rectus_part(width = left_width, thickness = left_thickness, length = lower_length, base_pos = c(-left_width - 1 * mm_, 0), emg_force_rel = emg_force_rel)
rectus_lb <- generate_rectus_part(width = left_width, thickness = left_thickness, length = mid_length, base_pos = c(-left_width - 1 * mm_, lower_length + 1 * cm_), emg_force_rel = emg_force_rel)
rectus_lc <- generate_rectus_part(width = left_width, thickness = left_thickness, length = upper_length, base_pos = c(-left_width - 1 * mm_, lower_length + mid_length + 2 * cm_), emg_force_rel = emg_force_rel)
rectus_ra <- generate_rectus_part(width = right_width, thickness = right_thickness, length = lower_length, base_pos = c(1 * mm_, 0), emg_force_rel = emg_force_rel)
rectus_rb <- generate_rectus_part(width = right_width, thickness = right_thickness, length = mid_length, base_pos = c(1 * mm_, lower_length + 1 * cm_), emg_force_rel = emg_force_rel)
rectus_rc <- generate_rectus_part(width = right_width, thickness = right_thickness, length = upper_length, base_pos = c(1 * mm_, lower_length + mid_length + 2 * cm_), emg_force_rel = emg_force_rel)
recti <- list(rectus_la, rectus_lb, rectus_lc, rectus_ra, rectus_rb, rectus_rc)
