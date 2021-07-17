##### Abdominal muscle simulation: rectus abdominis superior + transversus abdominis, two single-differential measurement channels #####

source("load_curves.r", local = T)

##### Muscle definitions #####

geom <- list()

emg_force_rel = 2

source("muscles/recti.r", local = T)
source("muscles/transversus.r", local = T)


##### Electrode Definitions #####

## Define the type and location of the electrodes.

electrodes <- create_point_electrodes(c(0, 7*cm_), c(0, 9*cm_), c(12*cm_, 20*cm_), c(14*cm_, 20.7 * cm_))


##### Volume Conductor Definitions #####


fat_thickness <- 0.01
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
B_sampling$min <- -2
B_sampling$max <- 4
B_sampling$dist <- 0.2
B_sampling$method <- "natural"

##### Muscle load profiles #####
# Define muscle load profiles: provide a function for each muscle that specifies the normalized target muscle force as a function of (continuous) time (in seconds).
# Values must be between 0 (inactive) and 1 (maximally active) at all times.

# We have six recti and two transversi
muscles_list <- c(recti, list(transversus_r, transversus_l))

time_limits = c(0, 30) * s_

recti_activation <- function(t) pmin(resp_exp_load(t,
                                  breath_length = 5 * s_,
                                  start_load = 0.1,
                                  max_load = 1/3,
                                  rel_length_to_max_load = 1/4,
                                  I_E_ratio = 1/2) +
                    cough_load(t,
                               t_start = 11.3 * s_,
                               cough_length = 1/3 * s_,
                               n_coughs = 3),
                    1)

# Force functions are defined in the same order the muscles appear in muscles_list.
force_functions <- c(recti_activation, recti_activation, recti_activation, 
							 recti_activation, recti_activation, recti_activation,
    function(t) pmin(0.05 + resp_exp_load(t,  # right transversus
                                  breath_length = 5 * s_,
                                  start_load = 0.1,
                                  max_load = 1/3,
                                  rel_length_to_max_load = 1/6,
                                  I_E_ratio = 1/2) +
                    cough_load(t,
                               t_start = 11.4 * s_,
                               cough_length = 1/3 * s_,
                               n_coughs = 3),
                    1),
    function(t) pmin(0.05 + resp_exp_load(t,  # left transversus
                                  breath_length = 5 * s_,
                                  start_load = 0.1,
                                  max_load = 1/3,
                                  rel_length_to_max_load = 1/6,
                                  I_E_ratio = 1/2) +
                    0.7 * cough_load(t,
                                     t_start = 11.6 * s_,
                                     cough_length = 1/3 * s_,
                                     n_coughs = 3),
                1))