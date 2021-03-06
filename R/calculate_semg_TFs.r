
#' Simple worker function wrapper to facilitate progression display
calculate_TF_Far_Mer_MU_wrapper <- function(x, p, ...) {
  #p(message=sprintf("x=%g", x))
  setTxtProgressBar(p, x)
  #message(paste0("    * thread <",x,"> mem <",sprintf("%.1f MB",pryr::mem_used()/1e6),">"))
  calculate_TF_Far_Mer_MU(...)
}


#' Calculate MU - electrode time transfer functions by a precise model
#'
#' Compute the transfer functions between MUs and electrodes according to the
#' model published by Farina and Merletti, https://doi.org/10.1109/10.923782
#'
#' This calculates the values of the MU-electrode transfer functions for a set
#' of Motor Units, sampled at a given set of frequencies. Note that the
#' potentially differing conduction velocities of the MUs will be taken into
#' account here, s.t. they do neither need to be taken care of in the
#' specification of .freqs, nor later on in the processing stage.
#'
#' @param MUs A data.frame containing a column 'MU.obj' containing MU objects
#'   and index columns 'MU' and 'muscle'.
#' @param electrodes A data.frame containing columns 'electrode' with the
#'   electrode ID and 'electrode.obj' containing the actual object.
#' @param vol_conductor A volume conductor object. Specifies the properties
#'   (i.e. thickness, conductivity) of the tissue layers separating muscle
#'   fibers and electrodes.
#' @param freqs A numerical vector containing the (non-angular) frequencies at
#'   which the transfer function is to be evaluated. Results will be returned in
#'   order.
#'   Usually should be generated by a call to the \code{setup.sampling} function
#'   and passing \code{sampling$freqs.to.calc}.
#' @param num_cores The number of cores that shall be used for computations. If
#'   equal to one (default), calculations are executed sequentially as usual.
#'   Otherwise, the parallel package is used for multi-process computations.
#' @param B_sampling 
#' @return A data.frame containing the values of the transfer functions at
#'   frequencies \code{.freqs} in column 'TF', represented as vectors of complex
#'   values.
#'
calculate_TFs_Far_Mer <- function(MUs, electrodes, vol_conductor,
                                  freqs, num_cores = 1, B_sampling = FALSE) {
  
    stopifnot(is.vector(freqs),
              is.numeric(freqs),
              all(Im(freqs) == 0),
              is.data.frame(MUs),
              is.data.frame(electrodes),
              is.nice.scalar(num_cores))
  
    data <- merge(MUs, electrodes)

    message(paste0("    ",whoami(2)," : Running <",length(data$MU.obj),"> threads on <",num_cores,"> cores ..."))
    
    #p <- progressr::progressor(steps=length(data$MU.obj))
    p <- txtProgressBar(max=length(data$MU.obj),style=3)
    
    if (num_cores > 1) {

        ## TODO: employ tools::pskill and tools::psnice to improve child process
        ## handling
        
        #if (file.exists("parlog_TF.txt"))
        #   file.remove("parlog_TF.txt")
      
        cl <- parallel::makeCluster(num_cores, outfile="")  #outfile = "parlog_TF.txt")
        future::plan(future::cluster, workers=cl, gc=TRUE)

        #data$TF <- parallel::clusterMap(cl,
        data$TF <- future.apply::future_mapply(
            calculate_TF_Far_Mer_MU_wrapper,
            seq_along(data$MU.obj),
            data$MU.obj,
            data$electrode.obj,
            MoreArgs = list(
                freqs = freqs,
                vol_conductor = vol_conductor,
                B_sampling = B_sampling,
                p = p),
            SIMPLIFY = FALSE,
            future.seed = TRUE)
        
        parallel::stopCluster(cl)
        
    } else {
        
        data$TF <- mapply(
            calculate_TF_Far_Mer_MU_wrapper,
            seq_along(data$MU.obj),
            data$MU.obj,
            data$electrode.obj,
            MoreArgs = list(
                freqs = freqs,
                vol_conductor = vol_conductor,
                B_sampling = B_sampling,
                p = p),
            SIMPLIFY = FALSE)
    }
    
    close(p)
    data[, c("muscle", "MU", "electrode", "TF")]
}


calculate_TF_Far_Mer_MU <- function(MU_obj, electrode_obj, freqs,
                                    vol_conductor, B_sampling) {

    ## TODO optimize via implementing in C(++)?
    
    calculate_TF_Far_Mer <- Vectorize(calculate_TF_point_Far_Mer, 'ang_freq')

    #print(".")

    calc_TF_for_fiber <- function(fiber) {

        Bfun <- create_Bfun(electrode = electrode_obj,
                            fiber = fiber,
                            vol_conductor = vol_conductor,
                            B_sampling = B_sampling)
        
        calculate_TF_Far_Mer(ang_freq = freqs / MU_obj$cv * 2 * pi,
                             electrode = electrode_obj,
                             fiber = fiber,
                             vol_conductor = vol_conductor,
                             Bfun = Bfun)
    }

    fiber_TFs <- MU_obj$fibers %>% lapply(calc_TF_for_fiber)

    stopifnot(is.list(fiber_TFs))
    
    MU_TF <- fiber_TFs %>% do.call(cbind, .) %>% rowSums
    ## MU.TF <- fiber.TFs %>% do.call(rbind, .) %>% colSums

    stopifnot(is.vector(MU_TF), length(MU_TF) == length(freqs))

    MU_TF
}


create_Bfun <- function(electrode, fiber, vol_conductor, B_sampling,
                        reltol = 1e-5, abstol = 1e-20, abstol_inc = 1e-5) {
                        ## approx_regions = c(T, T, T),
 
    stopifnot(class(fiber$shape) == "line_straight",
              fiber$rotation[1] == 0,
              fiber$rotation[3] == 0,
              fiber$innervation_zone[2] < 0)

    theta <- fiber$rotation[2]

    electrode_pos <- calc_rotated_electrode_pos(electrode, theta)
    fiber_innervation_zone_pos <-
        calc_rotated_innervation_zone_pos(fiber, theta)
    
    y0 <- fiber_innervation_zone_pos[2]
    x0 <- electrode_pos[1] - fiber_innervation_zone_pos[1]
    
    H_global <- function(kx, kz)
        vol_conductor$TF(kx = kx, kz = kz, y0 = y0) *
          electrode$TF(kx = kx, kz = kz, theta = theta) 

    ## Check exemplarily whether H.global is an even (at least in kx), real
    ## function. If it is not, the below calculation of B(kz) by just
    ## considering the positive part of the integral is not correct.
    stopifnot(!is.complex(H_global(1, 1)),
              H_global(1, 1) == H_global(-1, 1))

    B_kernel <- function(kx, kz) H_global(kx, kz) * cos(kx * x0)

    Bfun <- function(kz) {
        val <- try_to_integrate(
            params = list("inner", kz),
            fun = function(kx) B_kernel(kx, kz),
            lower = 0,
            upper = Inf,
            rel.tol = reltol,
            abs.tol = abstol) / pi
        ## print(paste("kz:", kz, "val:", val))
        val
    }

    Bfun_vec <- Vectorize(Bfun, 'kz')
    
    if (!is.list(B_sampling) || B_sampling$method == "none") {

        return(Bfun_vec)
        
    } else {

        stopifnot(is.nice.scalar(B_sampling$min),
                  is.nice.scalar(B_sampling$max),
                  is.nice.scalar(B_sampling$dist),
                  is.character(B_sampling$method),
                  Bfun(1) == Bfun(-1))
        
        B_sampled_freqs_log <- seq(from = B_sampling$min,
                                   to = B_sampling$max,
                                   by = B_sampling$dist)

        B_sampled_freqs <- 10^B_sampled_freqs_log
        
        Bvals <- Bfun_vec(B_sampled_freqs)
        n <- length(Bvals)

        for (i in which(abs(Bvals) < 10 * abstol)) {

            safety <- 0
            abstol_local <- abstol
            num_refines <- 0 
            
            ## Reduce tolerance until the result can be trusted
            while(safety < 10 && num_refines < 5) {
                abstol_local <- abstol_local * abstol_inc
                num_refines <- num_refines + 1
                Bvals[i] <- try_to_integrate(
                    params = list("inner", kz),
                    fun = function(kx) B_kernel(kx, B_sampled_freqs[i]),
                    lower = 0,
                    upper = Inf,
                    rel.tol = reltol,
                    abs.tol = abstol_local) / pi
                safety <- abs(Bvals[i]) / abstol_local
            }
        }

        nice_mask <- Bvals %>% sapply(is.nice.scalar)
        Bvals_nice <- Bvals[nice_mask]
        B_sampled_freqs_nice <- B_sampled_freqs[nice_mask]
        B_sampled_freqs_log_nice <- B_sampled_freqs_log[nice_mask]
        n_nice <- length(Bvals_nice)
        
        ## Create an interpolating function in the semi-log scale
        if (B_sampling$method == "monoH.FC") {
            Bfun_semilog_interp <- splinefun(x = B_sampled_freqs_log_nice,
                                             y = Bvals_nice,
                                             method = "monoH.FC")
        } else if (B_sampling$method == "natural") {
            Bfun_semilog_interp <- splinefun(x = B_sampled_freqs_log_nice,
                                             y = Bvals_nice,
                                             method = "natural")
        } else if (B_sampling$method == "akima") {
            Bfun_semilog_interp <- akimafun(x = B_sampled_freqs_log_nice,
                                            y = Bvals_nice)
        } else
            stop("Invalid 'method' argument supplied.")

        ## Calculate extrapolation slopes (towards 0 and towards Inf)
        expo_slope_inf_loglog <-
            calc_extrapolation_slope(Bvals_nice %>% abs %>% log10,
                                     B_sampled_freqs_log_nice, "right")
        expo_slope_zero_semilog <-
            calc_extrapolation_slope(Bvals_nice, 
                                     B_sampled_freqs_log_nice, "left")

        B_zero <- Bfun(0)
        
        ## Bfun_interp <- function(kz) {
        ##     if (log10(abs(kz)) < B_sampling$min) {
        ##         if (approx_regions[1]) {
        ##             ## extrapolate towards 0
        ##             increase <- expo_slope_zero_semi_log * 
        ##               (log10(abs(kz)) - B_sampled_freqs_log[1])
        ##             val <- Bvals[1] + increase
        ##         } else
        ##             val <- Bfun(kz)

        ##     } else if (log10(abs(kz)) > B_sampling$max) {
        ##         if (approx_regions[3]) {
        ##             ## extrapolate towards Inf
        ##             decrease_log <- expo_slope_Inf_double_log *
        ##               (log10(abs(kz)) - B_sampled_freqs_log[n_freqs])
        ##             val <- Bvals[n_freqs] * 10^decrease_log
        ##         } else
        ##             val <- Bfun(kz)
                
        ##     } else {
        ##         if (approx_regions[2]) {
        ##             ## interpolate the calculated values
        ##             if (kz > 0)
        ##                 val <- Bfun_semilog_interp(log10(kz))
        ##             else if (kz < 0)
        ##                 val <- Bfun_semilog_interp(log10(-kz))
        ##             else
        ##                 val <- Bfun(kz)
        ##         } else
        ##             val <- Bfun(kz)
        ##     }

        ##     val
        ## }
        
        return(function(kzs) Bfun_interp(kzs,
                                         B_zero = B_zero,
                                         B_left = Bvals_nice[1],
                                         B_right = Bvals_nice[n_nice],
                                         B_semilog_interp_fun =
                                           Bfun_semilog_interp,
                                         kz_left_log =
                                           B_sampled_freqs_log_nice[1],
                                         kz_right_log =
                                           B_sampled_freqs_log_nice[n_nice],
                                         expo_slope_zero_semilog =
                                           expo_slope_zero_semilog,
                                         expo_slope_inf_loglog =
                                           expo_slope_inf_loglog))
    }
}


Bfun_interp <- function(kzs, B_zero, B_left, B_right, B_semilog_interp_fun,
                        kz_left_log, kz_right_log, expo_slope_zero_semilog,
                        expo_slope_inf_loglog) {
    kzs <- abs(kzs)
    kzs_log <- log10(kzs)
    expo_zero_mask <- (kzs_log < kz_left_log)
    expo_inf_mask <- (kzs_log > kz_right_log)
    zero_mask <- (kzs == 0)
    interp_mask <- !expo_zero_mask & !expo_inf_mask & !zero_mask
    vals <- rep(0, length(kzs))

    ## extrapolate towards 0
    vals[expo_zero_mask] <-
        B_left + expo_slope_zero_semilog *
          (kzs_log[expo_zero_mask] - kz_left_log)
    
    ## extrapolate towards Inf
    vals[expo_inf_mask] <-
        B_right * 10^(expo_slope_inf_loglog *
                        (kzs_log[expo_inf_mask] - kz_right_log))
    
    ## insert B(0)
    vals[zero_mask] <- B_zero
    
    ## interpolate the calculated values
    vals[interp_mask] <- B_semilog_interp_fun(kzs_log[interp_mask])

    vals
}


calc_extrapolation_slope <- function(Bvals_abs_log, freqs_log, side) {

    stopifnot(length(Bvals_abs_log) == length(freqs_log),
              is.nice.vector(freqs_log),
              is.vector(Bvals_abs_log),
              is.numeric(Bvals_abs_log))

    nice_mask <- Bvals_abs_log %>% sapply(is.nice.scalar)
    Bvals_log_cleaned <- Bvals_abs_log[nice_mask]
    freqs_log_cleaned <- freqs_log[nice_mask]
    
    if (side == "left") {
        left <- T
        m <- 1
        n <- 2
        
    } else if (side == "right") {
        left <- F
        n <- length(Bvals_log_cleaned)
        m <- n - 1
        
    } else
        stop("Invalid parameter passed for 'side'. Must be 'left' or 'right'.")

    slope <- (Bvals_log_cleaned[n] - Bvals_log_cleaned[m]) /
      (freqs_log_cleaned[n] - freqs_log_cleaned[m])

    while (slope >= 0) {
        ## Successively take more points into account until regression yields
        ## negative slope.

        if (left)
            if (n == length(Bvals_log_cleaned))
                stop("No negative regression slope could be achieved.")
            else
                n <- n + 1
        else
            if (m == 1)
                stop("No negative regression slope could be achieved.")
            else
                m <- m - 1
        
        slope <-
            lm(Bvals_log_cleaned[m:n] ~ freqs_log_cleaned[m:n])$coefficients[2]
    }

    slope
}



#' Evaluate a muscle fiber - electrode TF at a single angular (!) frequency
#'
#' Note that this implementation assumes an even, real electrode transfer
#' function!
#' 
calculate_TF_point_Far_Mer <- function(ang_freq, electrode, fiber,
                                       vol_conductor, Bfun) {

    stopifnot(is.nice.scalar(ang_freq),
              is.function(Bfun))

    theta <- fiber$rotation[2]
    electrode_pos <- calc_rotated_electrode_pos(electrode, theta)
    fiber_innervation_zone_pos <-
        calc_rotated_innervation_zone_pos(fiber, theta)
     
    L1 <- (1 - fiber$innervation_zone_rel) * fiber$length
    L2 <- fiber$innervation_zone_rel * fiber$length
    zi <- fiber_innervation_zone_pos[3]
    z0 <- electrode_pos[2]
    
    kernel <- function(kzs) {

        k_eps <- kzs + ang_freq
        k_beta <- kzs - ang_freq

        C <- 2 * exp(-1i * k_eps * L1 / 2) * sin(k_eps * L1 / 2) / k_eps -
            2 * exp(1i * k_beta * L2 / 2) * sin(k_beta * L2 / 2) / k_beta

        val <- 1i * kzs * exp(-1i * kzs * (zi - z0)) * Bfun(kzs) * C

        val
    }

    integrate_kernel(kernel, ang_freq)
}


calc_rotated_electrode_pos <- function(electrode, theta) {

   if (theta == 0) {
        
        electrode_pos <- electrode$position
        
    } else {

        ## The transfer function model used below assumes the fiber to run in z
        ## direction. Therefore rotate the whole system (fiber + electrode)
        ## accordingly, if this is not the case.
        electrode_pos <- electrode$position %>% rotate2D(-theta)
    }

   electrode_pos
}


calc_rotated_innervation_zone_pos <- function(fiber, theta) {

    if (theta == 0) {
        
        fiber_innervation_zone_pos <- fiber$innervation_zone
        
    } else {

        ## The transfer function model used below assumes the fiber to run in z
        ## direction. Therefore rotate the whole system (fiber + electrode)
        ## accordingly, if this is not the case.
        fiber_innervation_zone_pos <-
            fiber$innervation_zone %>% rotate3D_lh(c(0, -theta, 0))
    }

   fiber_innervation_zone_pos
}


integrate_kernel <- function(kernel_fun, freq) {

    ## The integrand has three removable singularities: One at 0 and one each at
    ## -freq and +freq, respectively, i.e at k_eps = 0 and k_beta = 0.
    ## Therefore split the integration such that the singularities are on
    ## interval boundaries, as integrate can handle integrable singularities on
    ## boundary points.

    .try_to_integrate <- function(fun, lower, upper, ...)
        try_to_integrate(params = list("outer", freq),
                         fun = fun,
                         lower = lower,
                         upper = upper,
                         ...)

    kernel_re <- function(s) Re(kernel_fun(s))
    kernel_im <- function(s) Im(kernel_fun(s))
    
    integral_re_1 <- .try_to_integrate(kernel_re, -Inf, -abs(freq))
    integral_im_1 <- .try_to_integrate(kernel_im, -Inf, -abs(freq))
    
    if (is.equal(freq, 0)) {

        integral_re_2 <- 0
        integral_re_3 <- 0
        integral_im_2 <- 0
        integral_im_3 <- 0
        
    } else {

        integral_re_2 <- .try_to_integrate(kernel_re, -abs(freq), 0)
        integral_re_3 <- .try_to_integrate(kernel_re, 0, abs(freq))
        integral_im_2 <- .try_to_integrate(kernel_im, -abs(freq), 0)
        integral_im_3 <- .try_to_integrate(kernel_im, 0, abs(freq))        
    }
    
    integral_re_4 <- .try_to_integrate(kernel_re,
                               abs(freq), Inf)
    integral_im_4 <- .try_to_integrate(kernel_im,
                               abs(freq), Inf)
    
    integral_re <- integral_re_1 + integral_re_2 + integral_re_3 + integral_re_4
    integral_im <- integral_im_1 + integral_im_2 + integral_im_3 + integral_im_4

    (integral_re + 1i * integral_im) / 2 / pi
}



try_to_integrate <- function(params, fun, lower, upper, ...) {

    int <- try(integrate(fun, lower, upper, ...), silent = !(is_mode("debug")))

    if (class(int) == 'try-error') {

        ## Probably, numerical problems due to very small numbers have led to
        ## the error. Just set the integral value to 0.
        res <- 0

        if (is_mode("debug")) {
            print("Integration problem occurred. Setting integral value to 0.")
            print(params)
        }
        
    } else {

        res <- int$value
    }

    res
}
