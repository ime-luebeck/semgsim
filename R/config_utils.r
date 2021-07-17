create_rect_muscle <- function(width, thickness) {

    muscle <- list()
    class(muscle) <- "muscle"
    muscle$shape <- "rect"
    muscle$width <- width
    muscle$thickness <- thickness

    muscle
}


create_straight_line <- function() {

    line <- list()
    class(line) <- "line_straight"

    line$pos <- function(x)
        c(0, 0, x)

    line
}


#' Create a set of point electrodes
#'
#' Take an arbitrary number of positions on the 2D skin plane and return a data
#' frame containing point electrodes positioned there.
#'
#' @param ... Arbitrarily many 2D numerical vectors giving the desired electrode
#'   positions in 3D space with the y coordinate being automatically defined by
#'   the fact that the electrode is placed on the patient's skin. Positions are
#'   given as (x, z).
#'
#' @return A \code{data.frame} containing one row for each electrode. Columns
#'   are 'electrode' for the id of the electrode and 'electrode.obj' containing
#'   the actual object of type 'electrode.point'.
#' 
create_point_electrodes <- function(...) {

    positions <- list(...)
    
    for (position in positions)
        stopifnot(is.numeric(position),
                  is.vector(position),
                  length(position) == 2)

    electrode_objs <- lapply(positions, create_point_electrode)
    
    data.frame(electrode = seq_along(electrode_objs),
               electrode.obj = I(electrode_objs))
}


create_point_electrode <- function(position) {

    electrode <- list()
    class(electrode) <- "electrode_point"
    electrode$position <- position
    electrode$TF <- function(kx, kz, theta) 1
    electrode
}


create_load_curves <- function(time_limits, fs, ...) {

    stopifnot(is.vector(time_limits),
              is.numeric(time_limits),
              length(time_limits) == 2,
              time_limits[2] > time_limits[1],
              is.nice.scalar(fs))

    funcs <- list(...)

    for (func in funcs)
        stopifnot(is.function(func))

    samples <- seq(time_limits[1], time_limits[2], 1 / fs)
    
    load_curves <- lapply(funcs, function(func) {
        load <- func(samples)
        stopifnot(length(load) == length(samples))
        load
    })

    data.frame(muscle = seq_along(load_curves), load_curve = I(load_curves))
}   

#' Implements the planar three-layer volume conductor model of Farina and Rainoldi (1999),
#' https://doi.org/10.1016/S1350-4533(99)00075-2
#'
create_planar_3layer_volume_conductor <- function(fat_thickness, skin_thickness,
                                                  fat_cond, skin_cond,
                                                  muscle_cond_fiber,
                                                  muscle_cond_nofiber) {

    args <- as.list(environment())
    for (arg in args)
        stopifnot(is.nice.scalar(arg))
    
    vc <- list()
    class(vc) <- "vc.layered.planar"

    Ra <- muscle_cond_fiber / muscle_cond_nofiber
    Rc <- skin_cond / fat_cond
    Rm <- fat_cond / muscle_cond_nofiber

    alpha <- function(s, kya) kya + s * Rm * tanh(s)

    
    vc$TF <- function(kx, kz, y0) {

        ky <- sqrt(kx^2 + rep(kz^2, length(kx)))
        kya <- sqrt(kx^2 + rep(Ra * kz^2, length(kx)))

        res <- 2 / muscle_cond_nofiber * exp(-kya * abs(y0)) /
            ((1 + Rc) * cosh(ky * (fat_thickness + skin_thickness)) *
                 alpha(ky * (fat_thickness + skin_thickness), kya) +
                     (1 - Rc) * cosh(ky * (fat_thickness - skin_thickness)) *
                         alpha(ky * (fat_thickness - skin_thickness), kya))

        # For large values of kx, NaNs can occur due to multiple +/- Infs
        # cancelling each other. Symbolic analysis shows however that the
        # function converges monotonically towards 0 for kx -> +Inf.  Therefore
        # return 0 if this occurs.
        
        ifelse(is.nan(res), 0, res)
    }
    
    vc
}
                                                   
