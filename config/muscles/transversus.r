##### Transversus Abdominis #####
## Note that there's two of them; left and right.

geom$transversus <- list()
geom$transversus$angle <- 20 * deg_
geom$transversus$base_length <- 15 * cm_
geom$transversus$max_fiber_length <- 15 * cm_
geom$transversus$thickness <- 2 * mm_
geom$transversus$depth <- 0.5 * cm_
geom$transversus$mid_offset <- 3 * cm_
geom$transversus$base_pos <- c(0, 10 * cm_)


## Create the muscle object.
## Set the base shape of the muscle endplate.
## Currently supported shapes are:
##    - "rect". Will be interpreted as the unit square [0,1]x[0,1].

transversus_r <- create_rect_muscle(width = geom$transversus$base_length,
                                    thickness = geom$transversus$thickness)
transversus_l <- create_rect_muscle(width = geom$transversus$base_length,
                                    thickness = geom$transversus$thickness)

## Define the geometrical transformation by which the endplate base shape will
## be turned into its final shape. Must be some function that takes a 2D vector
## with values in [0,1] and returns a 3D vector, as this is also used to
## position the endplate in 3D space.

transversus_r$shape_transformation <-
    function(x) with(geom$transversus,
                     c(base_pos[1] + mid_offset,
                       - depth - x[2] * thickness,
                       base_pos[2] + x[1] * base_length))
transversus_l$shape_transformation <-
    function(x) with(geom$transversus,
                     c(base_pos[1] - mid_offset,
                       - depth - x[2] * thickness,
                       base_pos[2] + x[1] * base_length))

## Define the muscle fiber shape. This must be an object providing a method
## $pos(x) that implements a parametrization of the fiber shape curve. It takes
## a numerical argument x between 0 and 1 denoting the relative position on the
## fiber and returns a 3D vector containing the position of the curve in 3D
## space at that length. The curve is assumed to start at (0,0) and to be 
## normalized to length 1, i.e. the length of the curve given by pos(x) from x=0
## to x=1 is equal to one.
## Here, we set all fibers to be straight lines.

transversus_r$fiber_shape <- create_straight_line()
transversus_l$fiber_shape <- create_straight_line()

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

transversus_r$fiber_rotations <- function(pos_rel_to_muscle)
    c(0, -pi / 2 + geom$transversus$angle, 0)
transversus_l$fiber_rotations <- function(pos_rel_to_muscle)
    c(0, pi / 2 - geom$transversus$angle, 0)


transversus_r$MU_class <- "MU_circ"
transversus_l$MU_class <- "MU_circ"


## Define whether the conduction velocity of the muscle fibers can change with
## time.
## Here, we take it to be constant.

transversus_r$cv_variability <- "const"
transversus_l$cv_variability <- "const"

## Define the way in which the motor unit (MU) regions are determined from the
## muscle endplate region.

transversus_r$MU_distributor <-
    function(MUs) distribute_MUs_unif(MUs = MUs,
                                      muscle = transversus_r,
                                      tiling = c(5, 2))
transversus_l$MU_distributor <-
    function(MUs) distribute_MUs_unif(MUs = MUs,
                                      muscle = transversus_l,
                                      tiling = c(5, 2))
                                                             

## Define how the muscle fibers will be placed inside a given MU.

transversus_r$rel_fiber_distributor <- distribute_fibers_unif
transversus_l$rel_fiber_distributor <- distribute_fibers_unif


transversus_r$random_MU_params <- list(
					  fiber_density = list(shape = 8, scale = 9 / mm_^2,
										   min = 1 / mm_^2, max = 30 / mm_^2))
transversus_l$random_MU_params <- transversus_r$random_MU_params


## Define the base length of the muscles fibers as a function of their position
## on the muscle endplate.

transversus_r$fiber_base_length_calculator <- function(base_pos_rel_to_muscle) {
    start_point_abs <- transversus_r$shape_transformation(base_pos_rel_to_muscle)
    start_point_abs <- start_point_abs[c(1,3)] ## only x,z coordinates matter
    with(geom,
         min((start_point_abs[2] - 7.5 * cm_ -
                start_point_abs[1] * 0.5) /
                (0.5 - tan(transversus$angle)) *
                sqrt(1 + tan(transversus$angle)^2),
             transversus$max_fiber_length)
         )
}
transversus_l$fiber_base_length_calculator <-
    transversus_r$fiber_base_length_calculator

## Define a random relative variation of the fiber length.
## Don't forget to keep the relation between this and fiber_length_max (see
## below) up to date.

geom$transversus$variation <- 0.025

transversus_r$fiber_length_randomizer <-
    function() runif(n = 1,
                     min = 1 - geom$transversus$variation,
                     max = 1 + geom$transversus$variation)
transversus_l$fiber_length_randomizer <- transversus_r$fiber_length_randomizer

## Give an upper bound for the value
##
##    fiber_base_length_calculator() * fiber_length_randomizer().
##
## This will be used for estimating an upper bound on the time scale length of
## an IAP impulse response as measured on the skin.

transversus_r$fiber_length_max <-
    transversus_r$fiber_base_length_calculator(c(1,1)) *
      (1 + geom$transversus$variation)
transversus_l$fiber_length_max <- transversus_r$fiber_length_max


## Define, how the relative locations of the innervation zones of the muscle
## fibers will be determined.

transversus_r$innervation_zone_placer <-
    function(fiber) rnorm_bounded(min = 0.4, max = 0.6, n_sigma = 3)
transversus_l$innervation_zone_placer <- transversus_r$innervation_zone_placer


## Give an upper bound on
##
##     max(innervation_zone_placer(fiber), 1 - innervation_zone_placer(fiber)).
##
## This will be used for estimating an upper bound on the time scale length of
## an IAP impulse response as measured on the skin.

transversus_r$fiber_rel_half_length_max <- 0.6
transversus_l$fiber_rel_half_length_max <-
    transversus_r$fiber_rel_half_length_max
	
## MU Pool Organization

transversus_r$num_MUs <- 100
transversus_l$num_MUs <- 100


## Rate Coding and Recruitment ##
## Define the parameters that specify the characteristics of the rate coding
## strategy employed in the muscle.
transversus_rate_coding_params <- list()
if(exists('alternative_rate_coding') && alternative_rate_coding) {
	# Use our newly proposed rate coding model
	transversus_rate_coding_params$alternative_rate_coding <- TRUE
	transversus_rate_coding_params$A <- 20
	transversus_rate_coding_params$B <- 1.5
	transversus_rate_coding_params$C <- 30
	transversus_rate_coding_params$D <- 13
	transversus_rate_coding_params$E <- 8
	transversus_rate_coding_params$F <- 8
	transversus_rate_coding_params$G <- 0.05
} else {
	# Use the model of De Luca and Hostage (2010)
	transversus_rate_coding_params$alternative_rate_coding <- FALSE
	
	## Define the minimum and maximum firing rate (in pps, pulses per second) of all
	## MUs in the muscle.
	transversus_rate_coding_params$A <- 100
	transversus_rate_coding_params$B <- 0.2
	transversus_rate_coding_params$C <- -25
	transversus_rate_coding_params$D <- 7
	transversus_rate_coding_params$E <- 20
}


## Define the ratio of maximum muscle force at which MU recruitment is complete.
## Must be a value between 0 and 1.
rec_stop <- 0.6	

transversus_rate_coding_params$rec_model <- 2

## Shape parameter of the recruitment model
transversus_rate_coding_params$shape_b <- 0.2

a1 <- log(rec_stop) / transversus_l$num_MUs
a2 <- log(rec_stop/ transversus_rate_coding_params$shape_b) / transversus_l$num_MUs

if (transversus_rate_coding_params$rec_model == 1) {
  rec_start <- exp(a1 * 1)
} else 
  rec_start <- transversus_rate_coding_params$shape_b * 1 * exp(a2 * 1) / transversus_l$num_MUs

transversus_r$rate_coding_params <- transversus_rate_coding_params
transversus_l$rate_coding_params <- transversus_rate_coding_params

## Size principle parameters
LINEAR <- 1   # Variable is linearly related to the recruitment threshold
SQUARE <- 2   # Variable is related to the square of the recruitment threshold

# Intracellular conductivity, cf. Dimitrova (1974) http://www.ncbi.nlm.nih.gov/pubmed/4457323
cond_ic <- 1

# Range of fiber radius, cf. Hwang et al. (2013) [for index finger!]
# https://doi.org/10.3109/2000656X.2012.755988
min_fib_rad <- 12e-6
max_fib_rad <- 20e-6

transversus_MU_size_params <- list(rec_thresh = list(min=rec_start, max=rec_stop, lower_bound=0.005,
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
	
transversus_r$MU_size_params <- transversus_MU_size_params
transversus_l$MU_size_params <- transversus_MU_size_params
	
## Nonlinear force gain factor parameters
transversus_r$twitch_r <- 0.85
transversus_r$twitch_c <- 2.5
transversus_l$twitch_r <- 0.85
transversus_l$twitch_c <- 2.5