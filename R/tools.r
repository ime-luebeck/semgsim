mode_env <- new.env()
mode_vars <- list("debug" = ".DEBUG",
                  "interactive" = ".INTERACTIVE")

mode_defaults <- list("debug" = FALSE,
                      "interactive" = TRUE)

#' @export
is_mode <- function(mode) {
    stopifnot(mode %in% names(mode_vars))
    if (exists(mode_vars[[mode]], envir = mode_env)) {
        if (get(mode_vars[[mode]], envir = mode_env) == TRUE)
            TRUE
        else
            FALSE
    } else {
        mode_defaults[[mode]]
    }

}

#' @export
set_mode <- function(mode, val) {
    stopifnot(mode %in% names(mode_vars),
              is.logical(val),
              length(val) == 1)
    assign(mode_vars[[mode]], val, mode_env)
}

#' @export
toggle_mode <- function(mode) {
    stopifnot(mode %in% names(mode_vars))
    if (is_mode(mode))
        set_mode(mode, FALSE)
    else
        set_mode(mode, TRUE)
}

is.even <- function(x)
    (x %% 2 == 0)

is.odd <- function(x)
    !is.even(x)

is.not.null <- function(x)
    ! is.null(x)

is.not.nan <- function(x)
    ! is.nan(x)

is.not.inf <- function(x) {
    if (length(x) > 1)
        sapply(x, function(element) is.not.inf(element))
    else
        (!identical(x, Inf)) && (!identical(x, -Inf))
}

is.not.na <- function(x)
    ! is.na(x)

is.nice.scalar <- function(val)
    is.numeric(val) &
      is.not.null(val) &
      is.vector(val) &
      is.not.nan(val) &
      is.not.na(val) &
      is.not.inf(val) &
      (length(val) == 1)

is.nice.vector <- function(vec)
    is.vector(vec) &
      is.numeric(vec) &
      is.not.null(vec) &
      all(is.not.nan(vec) &
            is.not.na(vec) &
            is.not.inf(vec))

is.equal <- function(...)
    isTRUE(all.equal(...))

is_outside_unit_box <- function(pos) {

    if (dim(pos) %>% is.not.null)
        pos[,1] < 0 | pos[,1] > 1 |  pos[,2] < 0 | pos[,2] > 1
    else
        pos[1] < 0 || pos[1] > 1 || pos[2] < 0 || pos[2] > 1
}

none <- function(...)
    all(not(...))

expand_grid_vec <- function(seq1, seq2)
    cbind(Var1 = rep.int(seq1, length(seq2)),
          Var2 = rep.int(seq2, rep.int(length(seq1), length(seq2))))

#' Returns the indices of the first occurences of unique elements in a vector.
#'
#' @note This is _highly_ optimizable.
#' 
get_unique_idces <- function(vec) {

    unique_elements <- unique(vec)
    
    unique_idces <- unique_elements %>%
      lapply(function(element) which(vec == element)[1])

    unique_idces
}



#' Rotate a 3D vector around the three axes in a left-handed system
#'
#' Rotates a three-dimensional vector by the three angles specified around the
#' three axes of a left-handed coordinate system.
#'
#' If more than one angle is non-zero, the rotations will be applied in the
#' order x-y-z.
#'
#' @note Rotation is not as easily invertible as might seem apparent:
#'   \code{rotate3D_lh(rotate3D_lh(vec, angles), -angles) != vec}
#'   at least if there is more than one non-zero component in angles.
#' 
#' @param x The three-component vector(s) that is (are) to be rotated. In case
#'   of multiple vectors, they should be in the form of a 3xn array or matrix.
#' @param angles A three-component vector that contains the angles by which to
#'   rotate around the axes in radians and in the order x, y, z. Positive angles
#'   rotate counter-clockwise in left-handed coordinate systems.
#' 
#' @return A three-component vector that contains the rotated x, or a matrix
#'   containing the rotated columns of x.
#' 
rotate3D_lh <- function(x, angles) {

    stopifnot(is.vector(angles),
              length(angles) == 3,
              is.vector(x) || is.matrix(x) || is.array(x))

    if (is.matrix(x) || is.array(x)) {
        stopifnot(dim(x)[1] == 3)
    } else {
        stopifnot(length(x) == 3)
    }

    ## Since we're in a left-handed coordinate system...
    angles <- -angles
    
    Rx <- matrix(c(1, 0, 0,
                   0, cos(angles[1]), -sin(angles[1]),
                   0, sin(angles[1]), cos(angles[1])),
                 nrow = 3,
                 byrow = TRUE)

    Ry <- matrix(c(cos(angles[2]), 0, sin(angles[2]),
                   0, 1, 0,
                   -sin(angles[2]), 0, cos(angles[2])),
                 nrow = 3,
                 byrow = TRUE)

    Rz <- matrix(c(cos(angles[3]), -sin(angles[3]), 0,
                   sin(angles[3]), cos(angles[3]), 0,
                   0, 0, 1),
                 nrow = 3,
                 byrow = TRUE)

    res <- Rz %*% Ry %*% Rx %*% x

    if (dim(res)[2] == 1) {
        as.vector(res)
    } else {
        res
    }
}

#' Rotate a 2D vector by an angle
#'
#' Rotates a two-dimensional vector about the origin.
#'
#' @param x The two-component vector(s) that is (are) to be rotated. In case
#'   of multiple vectors, they should be in the form of a 2xn array or matrix.
#' @param angle A single numeric specifying the angle by which to rotate about
#'   the origin, in radians.
#' 
#' @return A two-component vector that contains the rotated x, or a matrix
#'   containing the rotated columns of x.
#' 
##' @author Eike Petersen
##' 
rotate2D <- function(x, angle) {

    stopifnot(is.nice.scalar(angle),
              is.vector(x) || is.matrix(x) || is.array(x))

    if (is.matrix(x) || is.array(x)) {
        stopifnot(dim(x)[1] == 2)
    } else {
        stopifnot(length(x) == 2)
    }
    
    R <- matrix(c(cos(angle), -sin(angle),
                  sin(angle), cos(angle)),
                nrow = 2,
                byrow = TRUE)

    res <- R %*% x

    if (dim(res)[2] == 1) {
        as.vector(res)
    } else {
        res
    }
}


calc_abs_rel_err <- function(exact, estimate)
    abs(exact - estimate) / abs(exact)

sample_points_in_circ_unif <- function(nsamples) {

    circ_midpoint <- c(0.5, 0.5)
    
    ## Note the adjusted distribution of the radii to ensure uniform
    ## sampling over the full area of the circle.
    radii <- 0.5 * sqrt(runif(n = nsamples, min = 0, max = 1))
    angles <- runif(n = nsamples, min = 0, max = 2*pi)

    positions <- cbind(circ_midpoint[1] + radii * cos(angles),
                       circ_midpoint[2] + radii * sin(angles))

    if (nsamples == 1)
        as.vector(positions)
    else
        positions
}



#' Helper function for random number generation
#'
#' Takes an intended range of values for the property and returns a randomly
#' instantiated value.
#'
#' This function considers the given range \code{[min, max]} to be the
#' \code{n}-sigma range of a Gaussian probability distribution, with the mean
#' lying in the center of the range. Samples outside of this interval
#' (i.e. outside of the \code{n}-sigma range, hence e.g. less than 5\% of all
#' samples for \code{n = 2}) are rejected and resampled.
#'
#' @param min The lower boundary of the interval to be sampled from.
#' @param max The upper boundary of the interval to be sampled from.
#' @param n_sigma An integer number indicating the range of the Gaussian
#'   distribution be accepted without resampling.
#' @return A single numerical value lying in the given interval.
#' @export
rnorm_bounded <- function(min = 0, max = 1, n_sigma = 2) {
    
    stopifnot(is.nice.scalar(min),
              is.nice.scalar(max),
              min <= max,
              is.nice.scalar(n_sigma))
    
    sampled <- FALSE

    while (!sampled) {

        property_val <- rnorm(n = 1,
                              mean = (min + max) / 2,
                              sd = (max - min) / (2 * n_sigma))
        
        if (min > property_val ||
            max < property_val)
            
            sampled <- FALSE
        else
            sampled <- TRUE
    }
    
    property_val
}

sample_bounded <- function(dist, dist_params, bounds) {
    stopifnot(is.nice.scalar(bounds[[1]]),
              is.nice.scalar(bounds[[2]]),
              length(bounds) == 2,
              is.function(dist),
              dist_params %>% lapply(is.nice.scalar) %>% simplify2array %>% all)
    
    sampled <- FALSE
    
	# min and max could potentially be passed in reverse order; we want to accept that.
    min <- min(bounds[[1]], bounds[[2]])
    max <- max(bounds[[1]], bounds[[2]])

    while (!sampled) {

        property_val <- do.call(dist, dist_params)
        
        if (min > property_val ||
            max < property_val)
            
            sampled <- FALSE
        else
            sampled <- TRUE
    }
    
    property_val    
}

replace_na <- function(x, by = 0) {
    x[is.na(x)] <- by
    x
}

dlaplace <- function(x, mean = 0, b = 1) {
    1 / 2 / b * exp(-abs(x - mean) / b)
}

sum_unequal_lengths <- function(lst_a, lst_b) {
    N <- max(length(lst_a), length(lst_b))
    if (length(lst_a) < N)
        lst_res <- c(lst_a, rep(0, N - length(lst_a))) + lst_b
    else
        lst_res <- lst_a + c(lst_b, rep(0, N - length(lst_b)))

    lst_res
}


#--------------------------------------------------------------------
#' whoami
#'
#' Returns name of calling function.
#' @export
#--------------------------------------------------------------------
#' @examples
#' whoami()
#--------------------------------------------------------------------
whoami = function (level=0, skip.cat=TRUE) {
  a <- sys.calls()
  if (length(a) - level <= 1) 
    name <- "main"
  else {
    if(skip.cat) {
      cnt = 1
      repeat {
        b <- as.character(a[[length(a) - cnt - level]])[1]
        if (b != "cat")
          break
        cnt = cnt + 1
      }
    } else {
      b <- as.character(a[[length(a) - 1 - level]])[1]
    }
    name <- gsub("[\\(\\)]","",b)
  }
  name
}

