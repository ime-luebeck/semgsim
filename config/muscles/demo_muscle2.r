##### Muscle Definitions #####
# Exactly the same as demo_muscle.r, just at a slightly different position (for testing multi-muscle setups).

geom$demo2 <- list()
geom$demo2$width <- 5 * cm_
geom$demo2$length <- 15 * cm_
geom$demo2$thickness <- 4 * mm_
geom$demo2$depth <- 0 * cm_

## Create the muscle object.
## Set the base shape of the muscle endplate.
## Currently supported shapes are:
##    - "rect". Will be interpreted as the unit square [0,1]x[0,1].

demo2 <- create_rect_muscle(width = geom$demo2$width,
                           thickness = geom$demo2$thickness)

demo2$identifier <- 'demo2'


## Define the geometrical transformation by which the endplate base shape will
## be turned into its final shape. Must be some function that takes a 2D vector
## with values in [0,1] and returns a 3D vector, as this is also used to
## position the endplate in 3D space.

demo2$shape_transformation <-
  function(x) c(8*cm_ + x[1] * demo2$width,
                - geom$demo2$depth - x[2] * demo2$thickness,
                0)


## Define the muscle fiber shape. This must be an object providing a method
## $pos(x) that implements a parametrization of the fiber shape curve. It takes
## a numerical argument x between 0 and 1 denoting the relative position on the
## fiber and returns a 3D vector containing the position of the curve in 3D
## space at that length. The curve is assumed to start at (0,0) and to be 
## normalized to length 1, i.e. the length of the curve given by pos(x) from x=0
## to x=1 is equal to one.
## Here, we set all fibers to be straight lines.

demo2$fiber_shape <- create_straight_line()


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

demo2$fiber_rotations <- function(pos_rel_to_muscle)
  c(0, 0, 0)



demo2$MU_class <- "MU_circ"

## Define how the muscle fibers will be placed inside a given MU.

demo2$rel_fiber_distributor <- distribute_fibers_unif


## Define the base length of the muscles fibers, e.g.  as a function of their
## position on the muscle endplate.

demo2$fiber_base_length_calculator <- function(base_pos_rel_to_muscle)
  geom$demo2$length


## Define a random relative variation of the fiber length.
## Don't forget to keep the relation between this and fiber_length_max (see
## below) up to date.

geom$demo2$variation <- 0.025

demo2$fiber_length_randomizer <-
  function() runif(n = 1,
                   min = 1 - geom$demo2$variation,
                   max = 1 + geom$demo2$variation)


## Give an upper bound for the value
##
##    fiber_base_length_calculator(pos) * fiber_length_randomizer().
##
## This will be used for estimating an upper bound on the time scale length of
## an IAP impulse response as measured on the skin.

demo2$fiber_length_max <-
  geom$demo2$length * (1 + geom$demo2$variation)


## Define, how the relative locations of the innervation zones of the muscle
## fibers will be determined.

demo2$innervation_zone_placer <-
  function(fiber) rnorm_bounded(min = 0.4, max = 0.6, n_sigma = 3)


## Give an upper bound on
##
##     max(innervation_zone_placer(fiber), 1 - innervation_zone_placer(fiber)).
##
## This will be used for estimating an upper bound on the time scale length of
## an IAP impulse response as measured on the skin.

demo2$fiber_rel_half_length_max <- 0.6


## Define whether the conduction velocity of the muscle fibers can change with
## time.
## Here, we take it to be constant.

demo2$cv_variability <- "const"


#### MU Pool Organization ####

demo2$num_MUs <- 3

## Randomly sampled parameters ##
demo2$random_MU_params <- list(fiber_density = list(shape = 8, scale = 9 / mm_^2,
                                                   min = 1 / mm_^2, max = 30 / mm_^2))

### Rate Coding and Recruitment ###
## Define the parameters that specify the characteristics of the rate coding
## strategy employed in the muscle.
demo2$rate_coding_params <- list()
demo2$rate_coding_params$alternative_rate_coding <- FALSE
demo2$rate_coding_params$rec_model <- 2

## Shape parameter of the recruitment model
demo2$rate_coding_params$shape_b <- 0.2

## Define the minimum and maximum firing rate (in pps, pulses per second) of all
## MUs in the muscle.
demo2$rate_coding_params$A <- 100
demo2$rate_coding_params$B <- 0.2
demo2$rate_coding_params$C <- -25
demo2$rate_coding_params$D <- 7
demo2$rate_coding_params$E <- 20

## Define the ratio of maximum muscle force at which MU recruitment is complete.
## Must be a value between 0 and 1.
rec_stop <- 0.6

a1 <- log(rec_stop) / demo2$num_MUs
a2 <- log(rec_stop/ demo2$rate_coding_params$shape_b) / demo2$num_MUs

if (demo2$rate_coding_params$rec_model == 1) {
  rec_start <- exp(a1 * 1)
} else 
  rec_start <- demo2$rate_coding_params$shape_b * 1 * exp(a2 * 1) / demo2$num_MUs

## Size principle parameters
LINEAR <- 1   # Variable is linearly related to the recruitment threshold
SQUARE <- 2   # Variable is related to the square of the recruitment threshold

# Intracellular conductivity, cf. Dimitrova (1974) http://www.ncbi.nlm.nih.gov/pubmed/4457323
cond_ic <- 1

# Range of fiber radius, cf. Hwang et al. (2013) [for index finger!]
# https://doi.org/10.3109/2000656X.2012.755988
min_fib_rad <- 12e-6
max_fib_rad <- 20e-6

demo2$MU_size_params <- list(rec_thresh = list(min=rec_start, max=rec_stop, lower_bound=0.005,
                                              relation=LINEAR, round=FALSE),
                            peak_force = list(min=1, max=35, lower_bound=0.1,
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
demo2$twitch_r <- 0.85
demo2$twitch_c <- 2.5


## Define the way in which the motor unit (MU) regions are determined from the
## muscle endplate region.

demo2$MU_distributor <-
  function(MUs) distribute_MUs_unif(MUs = MUs,
                                    muscle = demo,
                                    tiling = c(5, 2))

