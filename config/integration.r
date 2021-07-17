##### A simple demo test case to check whether everything works as intended #####

source("load_curves.r", local = T)
geom <- list()
source("muscles/demo_muscle.r", local = T)
source("muscles/demo_muscle2.r", local = T)

electrodes <- create_point_electrodes(c(0.03, 0.10), c(0.03, 0.12), c(0.11, 0.10), c(0.11, 0.12))


##### Volume Conductor Definitions #####


fat_thickness <- 0.003
skin_thickness <- 0.001
volume_conductor <-
    create_planar_3layer_volume_conductor(fat_thickness = fat_thickness,
                                          skin_thickness = skin_thickness,
                                          fat_cond = 0.05,
                                          skin_cond = 1,
                                          muscle_cond_fiber = 0.5,
                                          muscle_cond_nofiber = 0.1)
                                                           
class(volume_conductor) <- "volume_conductor"

## Define the base geometrical model we are using.
## Currently implemented: layered.

volume_conductor$model <- "layered"


## Define the function that will be used to calculate the MU -> electrode
## transfer functions.

calc_MU_electrode_time_TFs <- calculate_TFs_Far_Mer



##### Intracellular Action Potential Wave Shape #####

psi <- list()

## Define a function describing the first derivative of an IAP running in a
## negative spatial direction, evaluated in reverse direction. Must be given as
## a function of space. Must accept vectors. We call this psi.
##
##    psi(z) = (Vm(-z))'
##
## (If Vm(z) is the IAP wave shape running from +inf to -inf.)
##
## Here, we use Rosenfalck's approximation.

psi$transform <- function(fz) psi_Rosenfalck_transformed(fz)

## Indicate the rough length of the relevant portion of psi, i.e. the length [m]
## after which the signal can be considered close to zero.
## This information is used to calculate the length of a firing response, by
## considering this in conjuction with the fibre length and the conduction
## velocity.

psi$length <- 0.015



##### Technical Parameters #####

Fs <- 1024
max_freq_dom_sample_dist <- Inf

B_sampling <- list()
## The following settings led to oscillating firing responses in some MU-
## electrode combinations, hence they were refined.
## B_sampling$min <- -2
## B_sampling$max <- 4
## B_sampling$dist <- 0.2
## B_sampling$method <- "natural"
B_sampling$min <- -3
B_sampling$max <- 5
B_sampling$dist <- 0.1
B_sampling$method <- "natural"



##### Muscle load profiles #####
# Define muscle load profiles: provide a function for each muscle that specifies the normalized target muscle force as a function of (continuous) time (in seconds).
# Values must be between 0 (inactive) and 1 (maximally active) at all times.

muscles_list <- list(demo, demo2)

time_limits <- c(0, 10)

# Force functions are defined in the same order the muscles appear in muscles_list.
force_functions <- list(function(t) sin(t)^2, function(t) cos(t)^2)
