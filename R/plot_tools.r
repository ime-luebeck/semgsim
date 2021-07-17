#' @export
plot_geometry <- function(sim_result, bare = FALSE, merge_plots = TRUE,
                          dim_pairs = list(c(1, 2), c(3, 1), c(3, 2)),
                          electrode_labels = c(0, 4, 0)) {
    
    stopifnot(is.list(sim_result),
              is.logical(bare),
              is.logical(merge_plots),
              is.list(dim_pairs),
              length(dim_pairs[[1]]) == 2,
              length(electrode_labels) == length(dim_pairs),
              all(electrode_labels >= 0 & electrode_labels <= 4))

    
    ## Definitions
    
    skin_col <- "#FDB462"
    fat_col <- "#FFED6F"
    muscle_col <- "#FB8072"

    plot_margin <- 0.0075

    dim_names <- c("x", "y", "z")

    n_dim_pairs <- length(dim_pairs)
    n_muscles <- length(sim_result$muscles$muscle.obj)
    n_electrodes <- length(sim_result$electrodes$electrode.obj)

    
    ##### Calculate Electrode positions #####
    
    ## all electrodes are placed at the same y coord: on the skin
    electrode_ycoord <- sim_result$skin_thickness + sim_result$fat_thickness

    electrodes <- sim_result$electrodes
    electrodes_sorted <- electrodes[order(electrodes$electrode),]
    electrode_coords <- electrodes_sorted$electrode.obj %>%
        lapply(function(el) c(el$pos[1], electrode_ycoord, el$pos[2])) %>%
            unlist %>%
                array(dim = c(3, n_electrodes)) %>%
                    t

    
    ##### Calculate the positions of the 8 corners of each muscle.

    corners <- list()

    for (i in seq(1, n_muscles)) {
        corners[[i]] <- array(dim=c(8,3))

        muscle <- sim_result$muscles$muscle.obj[[i]]
        shape_trafo <- muscle$shape_transformation
        fiber_shape <- muscle$fiber_shape
        length_calculator <- muscle$fiber_base_length_calculator
        rotator <- muscle$fiber_rotations
        
        corners[[i]][1,] <- shape_trafo(c(0,0))
        corners[[i]][2,] <- shape_trafo(c(0,1))
        corners[[i]][3,] <- shape_trafo(c(1,0))
        corners[[i]][4,] <- shape_trafo(c(1,1))
        
        corners[[i]][5,] <- corners[[i]][1,] + length_calculator(c(0,0)) *
            rotate3D_lh(fiber_shape$pos(1), rotator(c(0, 0)))
        corners[[i]][6,] <- corners[[i]][2,] + length_calculator(c(0,1)) *
            rotate3D_lh(fiber_shape$pos(1), rotator(c(0, 1)))
        corners[[i]][7,] <- corners[[i]][3,] + length_calculator(c(1,0)) *
            rotate3D_lh(fiber_shape$pos(1), rotator(c(1, 0)))
        corners[[i]][8,] <- corners[[i]][4,] + length_calculator(c(1,1)) *
            rotate3D_lh(fiber_shape$pos(1), rotator(c(1, 1)))

    }

    
    ##### Calculate min and max x,y,z values to plot

    ## A bounding box (the convex hull) around the corner x,y,z values of the
    ## fibers is constructed. 
    
    xmax <- plyr::laply(corners, function(arr) arr[,1]) %>% max
    xmin <- plyr::laply(corners, function(arr) arr[,1]) %>% min
    ymax <- plyr::laply(corners, function(arr) arr[,2]) %>% max
    ymin <- plyr::laply(corners, function(arr) arr[,2]) %>% min
    zmax <- plyr::laply(corners, function(arr) arr[,3]) %>% max
    zmin <- plyr::laply(corners, function(arr) arr[,3]) %>% min

    ## Now calculate the necessary plot range from the distinct boundaries of
    ## fibers and electrodes. Note that this assumes that the fiber shapes are
    ## such that the convex hull constructed around the corner nodes is never
    ## left.

    plot_range <- list()

    ## xmin, xmax
    plot_range[[1]] <- c(min(c(xmin, electrode_coords[,1])) - plot_margin,
                         max(c(xmax, electrode_coords[,1])) + 2 * plot_margin)
    ## ymin, ymax
    plot_range[[2]] <- c(min(c(ymin, electrode_coords[,2])),
                         max(c(ymax, electrode_coords[,2])) + plot_margin)
    ## zmin, zmax
    plot_range[[3]] <- c(min(c(zmin, electrode_coords[,3])) - plot_margin,
                         max(c(zmax, electrode_coords[,3])) + 3 * plot_margin)


    ##### Perform the actual plot

    if (merge_plots) {
        if (is_mode("interactive")) dev.new()

        ## Tile the plot window into n_dim_pairs sections
        par(mfcol = c(n_dim_pairs, 1))

        if (bare)
            par(mar = c(0, 0, 0, 0), bty = "n", ## Remove margins, bounding box
                xaxs = "i", yaxs = "i", ## Disable autmatic axis extension
                xaxt = "n", yaxt = "n") ## Disable plotting the axis
    }

    for (i in seq_along(dim_pairs)) {

        dim_pair <- dim_pairs[[i]]
        
        if (!merge_plots) {
            
            if (is_mode("interactive")) dev.new()

            if (bare)
                par(mar = c(0, 0, 0, 0), bty = "n", 
                    xaxs = "i", yaxs = "i",
                    xaxt = "n", yaxt = "n")
        }
        
        plot(plot_range[[dim_pair[1]]], plot_range[[dim_pair[2]]],
             xlab = paste(dim_names[[dim_pair[1]]], "axis"),
             ylab = paste(dim_names[[dim_pair[2]]], "axis"),
             type = "n")

        if (!bare)
            title(paste(dim_names[[dim_pair[1]]], "-", dim_names[[dim_pair[2]]],
                        " plane", sep = ""))

        if (2 %in% dim_pair) {

            ## plot all layers
            depth_dim <- which(dim_pair == 2)
            add_layer(positions = c(0, sim_result$fat_thickness),
                      dim = depth_dim,
                      color = fat_col)
            add_layer(positions =
                        c(sim_result$fat_thickness,
                          sim_result$fat_thickness + sim_result$skin_thickness),
                      dim = depth_dim,
                      color = skin_col)
            add_layer(positions = c(-Inf, 0),
                      dim = depth_dim,
                      color = muscle_col)
        } else {

            ## only plot skin layer since no axis represents depth
            add_layer(positions = c(-Inf, +Inf),
                      dim = 1, ## doesn't matter
                      col = skin_col)
        }
        
        add_fibers(sim_result$MUs, dim_pair[1], dim_pair[2])
        add_electrode_markers(electrode_coords, dim_pair[1], dim_pair[2],
                              electrode_labels[i])
    }
}


add_layer <- function(positions, dim, color) {
    
    plot_lims <- par("usr")

    if (dim == 1) {
        
        if (positions[1] == -Inf)
            positions[1] <- plot_lims[1]
        if (positions[2] == Inf)
            positions[2] <- plot_lims[2]

        polygon(c(positions[1], positions[2], positions[2], positions[1]),
                c(plot_lims[3], plot_lims[3], plot_lims[4], plot_lims[4]),
                col = color,
                border = NA)

    } else if (dim == 2) {
        
        if (positions[1] == -Inf)
            positions[1] <- plot_lims[3]
        if (positions[2] == Inf)
            positions[2] <- plot_lims[4]

        polygon(c(plot_lims[1], plot_lims[1], plot_lims[2], plot_lims[2]),
                c(positions[1], positions[2], positions[2], positions[1]),
                col = color,
                border = NA)
    } else
        stop("Incorrect parameter passed for dim argument")
}


#' @importFrom plyr llply
add_fibers <- function(MUs, dim_1, dim_2) {
	# Do not plot ALL fibers because that takes ages and leads to memory overflow
	target_num_plot_fibers <- 5000
	MU_num_fibers = unlist(llply(MUs$MU.obj, function(MU.obj) MU.obj$num_fibers))
	sum(MU_num_fibers)
	total_num_fibers <- sum(MU_num_fibers)
	ratio <- target_num_plot_fibers / total_num_fibers
    for (MU in MUs$MU.obj)
        for (fiber in MU$fibers)
			if (ratio >= 1 || runif(1) < ratio)
				add_fiber(fiber, dim_1, dim_2)
}


add_fiber <- function(fiber, dim_1, dim_2) {

    xcoords <- c(fiber$base_pos_abs[dim_1],
                 fiber$innervation_zone[dim_1],
                 fiber$end_point[dim_1])
    
    ycoords <- c(fiber$base_pos_abs[dim_2],
                 fiber$innervation_zone[dim_2],
                 fiber$end_point[dim_2])
    
    lines(xcoords, ycoords, type = "l", lwd = 0.5)
    lines(xcoords, ycoords, type = "p", col = "#B3DE69", pch = 22)
}


add_electrode_markers <- function(electrode_coords, dim_1, dim_2, label_pos) {
    electrode_coords %>%
        cbind(seq(1, dim(electrode_coords)[1])) %>%
            plyr::a_ply(.margins = 1,
                        .fun = add_electrode_marker,
                        dim_1 = dim_1,
                        dim_2 = dim_2,
                        label_pos = label_pos)
}


add_electrode_marker <- function(data, dim_1, dim_2, label_pos) {

    points(x = data[dim_1], y = data[dim_2], type = "p", col = "#80B1D3",
           lwd = 5)

    if (label_pos != 0)
        text(x = data[dim_1],
             y = data[dim_2],
             labels = data[4],
             cex = 1.0,
             pos = label_pos)
}


#' Old & currently broken. Keeping this only for possible future use.
#'
#' @import ggplot2
#' @importFrom plyr "."
plot_TFs <- function(sim_result, MUs = 1:3, print_plot = TRUE) {

    freqs <- sim_result$sampling$freqs.to.calc
    
    ## Function that transforms sim results into a nicer format that is
    ## natively supported by ggplot2, plyr, dplyr, ...
    transformer <- function(df_row) {
        n <- length(df_row$TF[[1]])
        with(df_row,
             data.frame(muscle = muscle[[1]] %>% rep(n),
                         MU = MU[[1]] %>% rep(n),
                         electrode = electrode[[1]] %>% rep(n),
                         freq = freqs,
                         TF = TF[[1]]))
    }

    TFs <- plyr::ddply(.data = sim_result$TFs_MU_input_to_surf_potential %>%
                        dplyr::filter(MU %in% MUs),
                       .variables = .(muscle, MU, electrode),
                       .fun = transformer)

    TFs <- TFs %>%
        plyr::mutate(muscle = as.factor(muscle)) %>%
            plyr::mutate(MU = as.factor(MU)) %>%
                plyr::mutate(electrode = as.factor(electrode))
    
    TF_mag_plot <- ggplot(data = TFs,
                          aes(x = freq, y = abs(TF), colour = MU)) + 
                            geom_line(aes(group = MU)) +
                            facet_grid(electrode ~ muscle) +
                            ggtitle("Transfer Function Magnitudes") +
                            xlab("Frequency (Hz)") +
                            ylab("Magnitude (V)")

    TF_arg_plot <- ggplot(data = TFs,
                          aes(x = freq, y = Arg(TF), colour = MU)) + 
                            geom_line(aes(group = MU)) +
                            face_grid(electrode ~ muscle) +
                            ggtitle("Transfer Function Arguments") +
                            xlab("Frequency (Hz)") +
                            ylab("Argument (rad)")

    if (print_plot) {
        if (is_mode("interactive")) dev.new()
        print(TF_mag_plot)
        if (is_mode("interactive")) dev.new()
        print(TF_arg_plot)
    }

    invisible(list(TF_mag_plot = TF_mag_plot, TF_arg_plot = TF_arg_plot))
}

#' @import ggplot2
#' @importFrom plyr "."
#' @export
plot_firing_response_ffts <- function(sim_result, fmax = Inf, MUs = 1:3,
                                      print_plot = TRUE) {

    NFFT <- sim_result$sampling$NFFT

    transformer <- function(df_row) {
        with(df_row,
             data.frame(
                 muscle = rep(muscle[[1]], NFFT),
                 MU = MU[[1]] %>% rep(NFFT),
                 electrode = electrode[[1]] %>% rep(NFFT),
                 freq = sim_result$sampling$fft_freqs,
                 response =
                     abs(firing_response_fftvals[[1]])))
    }

    responses <- plyr::ddply(.data = sim_result$MU_firing_responses_fftvals %>%
                              dplyr::filter(MU %in% MUs),
                             .variables = .(muscle, MU, electrode),
                             .fun = transformer)

    responses <- responses %>%
        plyr::mutate(muscle = as.factor(muscle)) %>%
            plyr::mutate(MU = as.factor(MU)) %>%
                plyr::mutate(electrode = as.factor(electrode))
    
    plot_obj <-
        ggplot(data = responses %>%
                   dplyr::filter(freq > -fmax & freq < fmax),
               aes(x = freq, y = response, colour = MU)) + 
                 geom_line(aes(group = MU)) +
                 facet_grid(electrode ~ muscle, scales = "free_y") +
                 ggtitle("FFT Magnitude of MUAPs") +
                 xlab("Frequency (Hz)") +
                 ylab("MUAP Magnitude (V)")

    if (print_plot) {
        if (is_mode("interactive")) dev.new()
        print(plot_obj)
    }
    
    invisible(plot_obj)
}

#' @import ggplot2
#' @importFrom plyr "."
#' @export
plot_firing_responses <- function(sim_result, MUs = 1:3, nrow = 1, dir = "h",
                                  scales = "free_y", tlims = c(0, Inf), muscles = NA,
                                  electrode_confs = NA,
                                  electrodes = NA, print_plot = TRUE,
                                  bw = FALSE, facet_labels = TRUE, for_latex = FALSE) {
    ## TODO Implement free switching of variables between colour and facet dims
    
    NFFT <- sim_result$sampling$NFFT
    times <- seq(0, (NFFT - 1) / sim_result$Fs, 1 / sim_result$Fs)

    if (!is.numeric(electrodes))
        electrodes <- unique(sim_result$electrodes$electrode)

    if (!is.numeric(muscles))
        muscles <- unique(sim_result$muscles$muscle)
    
    transformer <- function(df_row) {
        with(df_row,
             data.frame(
                 muscle = rep(muscle[[1]], NFFT),
                 MU = MU[[1]] %>% rep(NFFT),
                 electrode = electrode[[1]] %>% rep(NFFT),
                 time = times,
                 potential =
                     firing_response_fftvals[[1]] %>% fft(inverse = TRUE) %>% Re))
    }

    responses <- plyr::ddply(.data = sim_result$MU_firing_responses_fftvals %>%
                              dplyr::filter(MU %in% MUs &
                                            electrode %in% electrodes &
                                            muscle %in% muscles),
                             .variables = .(muscle, MU, electrode),
                             .fun = transformer)

    if (!is.na(electrode_confs))
        responses <- apply_electrode_configs(responses,
                                             electrode_confs)

    responses <- responses %>%
        dplyr::mutate(muscle = as.factor(muscle)) %>%
            dplyr::mutate(MU = as.factor(MU)) %>%
        dplyr::mutate(electrode = as.factor(electrode))

    if (!bw)
        plot_obj <- 
        ggplot(data = responses %>%
                   dplyr::filter(time >= tlims[1] & time <= tlims[2]),
               aes(x = time / ms_, y = potential, colour = MU),
               environment = environment())
    else
        plot_obj <- ggplot(data = responses %>%
                               dplyr::filter(time >= tlims[1] & time <= tlims[2]),
                           aes(x = time / ms_, y = potential, group=electrode),
                           environment = environment())

    if (for_latex) {
      ylabel <- "$\\mathrm{SFAP}$ (V)"
      xlabel <- "Time $t$ (ms)"
    } else {
      ylabel <- "MUAP (V)"
      xlabel <- "Time (ms)"
    }
    
    if (length(MUs) == 1) {
      plot_obj <- plot_obj + 
        geom_line(aes(group=electrode, linetype = electrode)) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
        xlab(xlabel) + ylab(ylabel)
      
      if (nlevels(responses$muscle) == 1) {
        
        plot_obj <- plot_obj +
          ggtitle("Motor Unit Action Potentials")
        
      } else {
        
        plot_obj <- plot_obj + 
          facet_wrap(~ muscle, nrow = nrow, scales = scales) +
          ggtitle("MUAPs of Each Muscle")
        
      }
      
    } else {
      plot_obj <- plot_obj + 
          geom_line(aes(group = MU)) +
          scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
          xlab("Time (ms)") +
          ylab(ylabel)
      
      if (nlevels(responses$muscle) == 1 &&
          nlevels(responses$electrode) == 1) {
  
          plot_obj <- plot_obj +
            ggtitle("Motor Unit Action Potentials")
  
      } else if (nlevels(responses$muscle) > 1 &&
                 nlevels(responses$electrode) == 1) {
  
          ## plot_obj <- plot_obj + 
          ##   facet_wrap(~ muscle, nrow = nrow, scales = scales, dir = dir) +
          ##   ggtitle("MUAPs of Each Muscle")
          plot_obj <- plot_obj + 
            facet_wrap(~ muscle, nrow = nrow, scales = scales) +
            ggtitle("MUAPs of Each Muscle")
  
      } else if (nlevels(responses$muscle) == 1 &&
                 nlevels(responses$electrode) > 1) {
  
          ## plot_obj <- plot_obj + 
          ##   facet_wrap(~ electrode, nrow = nrow, scales = scales, dir = dir) +
          ##   ggtitle("MUAPs at Each Electrode")
          plot_obj <- plot_obj + 
            facet_wrap(~ electrode, nrow = nrow, scales = scales) +
            ggtitle("MUAPs at Each Electrode")
  
      } else {
  
          plot_obj <- plot_obj + 
            facet_grid(electrode ~ muscle, scales = scales) +
            ggtitle("MUAPs of Each Muscle at Each Electrode")
      }
    }
    
    if (!facet_labels)
        plot_obj <- plot_obj + theme(
          strip.background = element_blank(),
          strip.text.x = element_blank()
        )
    
    if (print_plot) {
        if (is_mode("interactive")) dev.new()
        print(plot_obj)
    }
    
    invisible(plot_obj)
}


#' @import ggplot2
#' @importFrom plyr "."
#' @export
plot_impulse_responses <- function(sim_result, MUs = 1:3, nrow = 1, dir = "h",
                                   scales = "free_y", print_plot = TRUE) {

    NFFT <- sim_result$sampling$NFFT
    times <- seq(0, (NFFT - 1) / sim_result$Fs, 1 / sim_result$Fs)

    transformer <- function(df_row) {
        fftvals <- df_row$TF[[1]] %>% mirror_hermitian_fft(sim_result$sampling)
        with(df_row,
             data.frame(muscle = muscle[[1]] %>% rep(NFFT),
                        MU = MU[[1]] %>% rep(NFFT),
                        electrode = electrode[[1]] %>% rep(NFFT),
                        time = times,
                        impulse = fftvals %>% fft(inverse = TRUE) %>% Re))
    }

    impulses <- plyr::ddply(.data =
                              sim_result$TFs_MU_input_to_surf_potential %>%
                              dplyr::filter(MU %in% MUs),
                            .variables = .(muscle, MU, electrode),
                            .fun = transformer)

    impulses <- impulses %>%
        dplyr::mutate(muscle = as.factor(muscle)) %>%
            dplyr::mutate(MU = as.factor(MU)) %>%
                dplyr::mutate(electrode = as.factor(electrode))
    
    plot_obj <-
        ggplot(data = impulses,
               aes(x = time / ms_, y = impulse, colour = MU),
               environment = environment()) + 
                 geom_line(aes(group = MU)) +
                 scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
                 xlab("Time (ms)") +
                 ylab("MU Impulse Response (V)")

    if (nlevels(sim_result$muscle_contribs$emg$muscle) == 1 &&
        nlevels(sim_result$muscle_contribs$emg$electrode) == 1) {

        plot_obj <- plot_obj +
          ggtitle("MU Input to Electrode Impulse Responses")

    } else if (nlevels(sim_result$muscle_contribs$emg$muscle) > 1 &&
               nlevels(sim_result$muscle_contribs$emg$electrode) == 1) {

        plot_obj <- plot_obj +
          facet_wrap(~ muscle, nrow = nrow, scales = scales, dir = dir) + 
          ggtitle("MU Input to Electrode Impulse Responses For Each Muscle")

    } else if (nlevels(sim_result$muscle_contribs$emg$muscle) == 1 &&
               nlevels(sim_result$muscle_contribs$emg$electrode) > 1) {

        plot_obj <- plot_obj +
          facet_wrap(~ electrode, nrow = nrow, scales = scales, dir = dir) + 
          ggtitle("MU Input to Each Electrode Impulse Responses")
                
    } else {

        plot_obj <- plot_obj +
          facet_grid(electrode ~ muscle, scales = scales) + 
          ggtitle("MU Input to Each Electrode Impulse Responses For Each Muscle")
    }

    if (print_plot) {
        if (is_mode("interactive")) dev.new()
        print(plot_obj)
    }
    
    invisible(plot_obj)
}

#' @import ggplot2
#' @importFrom plyr "."
#' @export
plot_firing_rate_funcs <- function(model, MUs = 1:3, use_transformation = TRUE, print_plot = TRUE,
                                   muscles = NA, bw = FALSE, for_latex = FALSE) {

    load <- c(seq(0, 0.05, 0.0001), seq(0.05, 0.1, 0.001), seq(0.1, 1, 0.01))
    n  <- length(load)
    
    if (!is.numeric(muscles))
      muscles <- unique(model$muscles$muscle)
    
    if (use_transformation)
        if (for_latex) {
          x_label <- "Normalized Muscle Force $\\tilde F$"
        } else {
          x_label <- "Normalized Muscle Force"
        }
    else
      if (for_latex) {
        x_label <- "Common Drive $\\mathrm{CD}$"
      } else {
        x_label <- "Common Drive"
      }
    
    if (for_latex)
      y_label <- "Firing Rate $\\lambda_i$ (imp/s)"
    else
      y_label <- "Firing Rate (imp/s)"
    
    transformer <- function(df_row) {
        if (use_transformation)
          transformation <- model$act_funcs_df$force_to_act_fun[[df_row$muscle[[1]]]]
        else
          transformation <- function(x) x
        
        with(df_row,
             data.frame(muscle = muscle[[1]] %>% rep(n),
                         MU = MU[[1]] %>% rep(n),
                         load = load,
                         rate = MU_rate_function[[1]](load %>% transformation)))
    }

    rates <- plyr::ddply(.data = model$MUs_rate_functions %>%
                          dplyr::filter(MU %in% MUs & muscle %in% muscles),
                         .variables = .(muscle, MU),
                         .fun = transformer)

    rates <- rates %>%
        dplyr::mutate(muscle = as.factor(muscle)) %>%
            dplyr::mutate(MU = as.factor(MU))

    if (!bw)
      plot_obj <-
          ggplot(data = rates %>% dplyr::filter(rate > 0),
                 aes(x = load, y = rate, colour = MU)) + 
                   geom_line(aes(group = MU)) + ylim(c(0, max(rates$rate))) + 
                   ggtitle("MU Firing Rate Functions") +
                   xlab(x_label) + ylab(y_label)
    else
      plot_obj <-
        ggplot(data = rates %>% dplyr::filter(rate > 0),
               aes(x = load, y = rate)) + 
        geom_line(aes(group = MU)) + ylim(c(0, max(rates$rate))) + 
        ggtitle("MU Firing Rate Functions") +
        xlab(x_label) + ylab(y_label)
    
    if (nlevels(rates$muscle) > 1)
        plot_obj <- plot_obj + facet_grid(muscle ~ .)

    if (print_plot) {
        if (is_mode("interactive")) dev.new()
        print(plot_obj)
    }
    
    invisible(plot_obj)
}


#' @import ggplot2
#' @importFrom ggpubr ggarrange
#' @export
plot_muscle_contribs <- function(muscle_contribs, tlims = c(0, Inf),
                                    xlims = NULL, ylims = NULL, electrodes = NA,
                                    electrode_confs = NULL,
                                    analysis = c("plain", "rms", "fourier",
                                        "density"),
                                    nrow = 1, scales = "free_y", dir = "h",
                                    print_plot = TRUE, 
                                    for_latex = FALSE, ...) {
    stopifnot(is.character(analysis),
              is.nice.scalar(nrow),
              is.character(scales),
              is.numeric(tlims),
              length(tlims) == 2,
              is.data.frame(muscle_contribs$emg),
			  is.data.frame(muscle_contribs$force))

    analysis <- analysis[1]
    
    if (!is.numeric(electrodes))
        electrodes <- unique(muscle_contribs$emg$electrode)
    
    emg_contribs <- muscle_contribs$emg %>%
      dplyr::filter(electrode %in% electrodes) %>%
      dplyr::filter(time >= tlims[1] & time <= tlims[2])
    
	force_contribs <- muscle_contribs$force %>%
      dplyr::filter(time >= tlims[1] & time <= tlims[2])
	
    if (!is.null(electrode_confs))
        emg_contribs <- apply_electrode_configs(emg_contribs, electrode_confs)

    emg_contribs <- emg_contribs %>%
      dplyr::mutate(muscle = as.factor(muscle)) %>%
      dplyr::mutate(electrode = as.factor(electrode))
    
    force_contribs <- force_contribs %>% dplyr::mutate(muscle = as.factor(muscle))
	
    if (analysis == "plain") {
        
        x_name <- "time"
        x_label <- "Time (s)"
        y_name <- "potential"
        if (for_latex)
          y_label <- "Electric Potential (V)"
        else
          y_label <- "Electric Potential (V)"
        geom_fun <- ggplot2::geom_line
        addon_fun <- NULL
        plot_data <- emg_contribs
        analysis_label <- "Raw data"

    } else {

        post_process_fun <- get(paste("post_process_", analysis, sep = ""))
        plot_data <- emg_contribs %>%
          post_process_wrap(post_process_fun, ...)
        post_process_config <- get(paste("get_post_process_config_", analysis,
                                         sep = ""))(for_latex)
        x_name <- post_process_config$x_name
        x_label <- post_process_config$x_label
        y_name <- post_process_config$y_name
        y_label <- post_process_config$y_label
        geom_fun <- post_process_config$geom_fun
        addon_fun <- post_process_config$addon_fun
        analysis_label <- post_process_config$analysis_label
    }

    if (is.null(addon_fun))
        addon_fun <- function(plot_obj, df) plot_obj

    emg_plot <-
        (ggplot(data = plot_data,
                aes_string(x = x_name, y = y_name)) +
                  geom_fun() +
                  coord_cartesian(xlim = xlims, ylim = ylims) +
                  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
                  xlab(x_label) +
                  ylab(y_label)) %>%
                  addon_fun(plot_data)
    
	force_plot <- ggplot(data = force_contribs, aes(x = time, y = force)) +
									geom_fun() + coord_cartesian(xlim = xlims, ylim = ylims) +
									scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
									xlab('Time (s)') +
									ylab('Force (k*N)')
	
    if (nlevels(emg_contribs$muscle) == 1 &&
        nlevels(emg_contribs$electrode) == 1) {

        emg_plot <- emg_plot +
          ggtitle(paste("Contribution of the muscle to the measured electrode",
                        "potential:", analysis_label))
						
		force_plot <- force_plot +
          ggtitle("Force generated by the muscle")
		  
		 plot_obj <- ggarrange(force_plot, emg_plot, ncol=2)

    } else if (nlevels(emg_contribs$muscle) > 1 &&
               nlevels(emg_contribs$electrode) == 1) {

        emg_plot <- emg_plot +
          facet_wrap(~ muscle, nrow = nrow, scales = scales, dir = dir) + 
          ggtitle(paste("Contribution of each muscle to the measured electrode",
                        "potential:", analysis_label))
						
		force_plot <- force_plot + facet_wrap(~ muscle, nrow = nrow, scales = scales, dir = dir) +
												ggtitle("Force generated by each muscle")
												
		plot_obj <- ggarrange(force_plot, emg_plot, nrow=2)

    } else if (nlevels(emg_contribs$muscle) == 1 &&
               nlevels(emg_contribs$electrode) > 1) {

        emg_plot <- emg_plot +
          facet_wrap(~ electrode, nrow = nrow, scales = scales, dir = dir) + 
          ggtitle(paste("Contributions of the muscle to each measured",
                        "electrode potential:", analysis_label))
        
        force_plot <- force_plot + ggtitle("Force generated by the muscle")
						
		plot_obj <- ggarrange(force_plot, emg_plot, ncol=2)
		
    } else {

        emg_plot <- emg_plot +
          facet_grid(electrode ~ muscle, scales = scales) + 
          ggtitle(paste("Contributions of each muscle to each measured", 
                        "electrode potential:", analysis_label))
						
        force_plot <- force_plot +
          facet_wrap( ~ muscle, scales = scales, nrow=1) + 
          ggtitle("Force generated by each muscle")					

		plot_obj <- ggarrange(force_plot, emg_plot, nrow=2, heights=c(2, nlevels(emg_contribs$electrode)))
    }

    if (print_plot) {
        if (is_mode("interactive")) dev.new()
        print(plot_obj)
    }
    
    invisible(plot_obj)
}


#' @import ggplot2
#' @export
plot_surface_potentials <- function(surface_potentials, tlims = c(0, Inf),
                                    xlims = NULL, ylims = NULL, electrodes = NA,
									                  electrode_confs = NULL,
                                    analysis = c("plain", "rms", "fourier",
                                        "density"),
                                    nrow = 1, scales = "free_y", dir = "h",
									                  downsampling_factor = 1,
                                    print_plot = TRUE, facet_labels = TRUE,
									                  for_latex = FALSE, ...) {
    stopifnot(is.character(analysis),
				  is.nice.scalar(nrow),
				  is.character(scales),
				  is.numeric(tlims),
				  length(tlims) == 2,
				  is.data.frame(surface_potentials))

    analysis <- analysis[1]

    if (!is.numeric(electrodes))
        electrodes <- unique(surface_potentials$electrode)
    
    Tsampling <- surface_potentials$time[2] - surface_potentials$time[1]
    stopifnot(Tsampling > 0)
    surface_potentials <- surface_potentials %>%
        dplyr::filter(electrode %in% electrodes) %>%
        dplyr::filter(time >= tlims[1] & 
                        time <= tlims[2] & 
                        time %% (Tsampling * downsampling_factor) == 0)
		
    if (!is.null(electrode_confs))
        surface_potentials <- apply_electrode_configs(surface_potentials,
                                                      electrode_confs)

    if (analysis == "plain") {
        
        x_name <- "time"
        x_label <- "Time (s)"
        y_name <- "potential"
        if (!for_latex)
          y_label <- "Electric Potential (V)"
        else
          y_label <- "Electric Potential (V)"
        geom_fun <- ggplot2::geom_line
        addon_fun <- NULL
        plot_data <- surface_potentials
        analysis_label <- "Raw data"

    } else {

        post_process_fun <- get(paste("post_process_", analysis, sep = ""))
        plot_data <- surface_potentials %>%
          post_process_wrap(post_process_fun, ...)
        post_process_config <- get(paste("get_post_process_config_", analysis,
                                         sep = ""))(for_latex)
        x_name <- post_process_config$x_name
        x_label <- post_process_config$x_label
        y_name <- post_process_config$y_name
        y_label <- post_process_config$y_label
        geom_fun <- post_process_config$geom_fun
        addon_fun <- post_process_config$addon_fun
        analysis_label <- post_process_config$analysis_label
    }

    if (is.null(addon_fun))
        addon_fun <- function(plot_obj, df) plot_obj

    if (is.not.null(y_name))
        obj_aes <- aes_string(x = x_name, y = y_name)
    else
        obj_aes <- aes_string(x = x_name)

    #plot_data <- plot_data %>%
    #    dplyr::mutate(electrode = as.factor(electrode))
    
    plot_obj <-
        (ggplot(data = plot_data,
                obj_aes) +
                  geom_fun() +
                  coord_cartesian(xlim = xlims, ylim = ylims) +
                  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
                  xlab(x_label) +
                  ylab(y_label)) %>%
                  addon_fun(plot_data)
    
    if (nlevels(plot_data$electrode) == 1) {

        plot_obj <- plot_obj +
          ggtitle(paste("Measured Potential:", analysis_label))

    } else {

        plot_obj <- plot_obj +
          facet_wrap(~ electrode, nrow = nrow, scales = scales, dir = dir) + 
          ggtitle(paste("Measured Potentials:", analysis_label))
                
    }
    
    if (!facet_labels)
      plot_obj <- plot_obj + theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
      )
    

    if (print_plot) {
        if (is_mode("interactive")) dev.new()
        print(plot_obj)
    }
    
    invisible(plot_obj)
}


#' @import ggplot2
#' @importFrom plyr "."
#' @export
plot_impulse_trains <- function(sim_result, MUs = 1:3, muscles = NA,
                                tlims = c(0, Inf), print_plot = TRUE) {

    if (!is.numeric(muscles))
        muscles <- 1:nrow(sim_result$muscles)

    N <- length(sim_result$MU_contribs$contribs[[1]]$force)
    times <- ((1:N) - 1) / sim_result$sampling$Fs
    times_mask <- (times >= tlims[1]) & (times <= tlims[2])
    n <- sum(times_mask)
    
    transformer <- function(df_row) {
        with(df_row,
             data.frame(muscle = muscle[[1]] %>% rep(n),
                        MU = MU[[1]] %>% rep(n),
                        time = times[times_mask],
                        impulse = impulse_train[[1]][times_mask]))
    }

    instants <- plyr::ddply(.data = sim_result$MU_impulse_trains %>%
                              dplyr::filter(MU %in% MUs & muscle %in% muscles),
                            .variables = .(muscle, MU),
                            .fun = transformer)

    instants <- instants %>%
      dplyr::mutate(muscle = as.factor(muscle)) %>%
      dplyr::mutate(MU = as.factor(MU)) %>%
      dplyr::mutate(impulse = as.numeric(impulse))

    plot_obj <- ggplot(data = instants,
                       aes(x = time, y = impulse, colour = MU),
                       environment = environment()) +
                         geom_col() +
                         geom_point() +
                         xlab("Time [s]") +
                         ylab("Firing Impulse")

    if (nlevels(instants$muscle) > 1)
        plot_obj <- plot_obj + facet_grid(muscle ~ .) +
          ggtitle("MU Impulse Trains For Each Muscle")
    else
        plot_obj <- plot_obj + ggtitle("MU Impulse Trains")

    if (print_plot) {
        if (is_mode("interactive")) dev.new()
        print(plot_obj)
    }
    
    invisible(plot_obj)
}


#' @import ggplot2
#' @importFrom plyr "."
#' @export
plot_force_functions <- function(sim_result, tlims = c(0, Inf),
                                muscles = NA, print_plot = TRUE) {

    if (!is.numeric(muscles))
        muscles <- 1:nrow(sim_result$muscles)

    N <- length(sim_result$MU_contribs$contribs[[1]]$force)
    times <- ((1:N) - 1) / sim_result$sampling$Fs
    times_mask <- (times >= tlims[1]) & (times <= tlims[2])
    n <- sum(times_mask)
    
    transformer <- function(list_elem) {
      func <- list_elem[[1]]
      id <- list_elem[[2]]
      data.frame(muscle = id %>% rep(n), time = times[times_mask], force=func(times[times_mask]))
    }

    # zip force functions with corresponding muscle ids to a single list
    funcs_and_ids <- mapply(list, sim_result$force_functions, muscles, SIMPLIFY = FALSE)
    
    loads <- plyr::ldply(funcs_and_ids, transformer)
    
    loads <- loads %>% dplyr::mutate(muscle = as.factor(muscle))
    
    plot_obj <- ggplot(data = loads,
                       aes(x = time, y = force, colour = muscle)) +
                         geom_line() +
                         ggtitle("Muscle Force Profiles") +
                         xlab("Time [s]") +
                         ylab("Relative muscle force")

    if (print_plot) {
        if (is_mode("interactive")) dev.new()
        print(plot_obj)
    }
    
    invisible(plot_obj)
}


#' @import ggplot2
#' @export
plot_emg_force_rel <- function(surface_potentials_or_list_of_those, stairs_start, step_length, 
                               num_steps, transient_length,
                               electrodes = NA, electrode_confs = NULL,
                               nrow = 1, surface_potentials_or_list_of_those_2 = NULL, errorbars = FALSE,
                               print_plot = TRUE, for_latex=FALSE, ...) {
  
  stopifnot(is.nice.scalar(nrow),
            is.list(surface_potentials_or_list_of_those),
            is.data.frame(surface_potentials_or_list_of_those) || 
              is.data.frame(surface_potentials_or_list_of_those[[1]]))
  
  if (!is.data.frame(surface_potentials_or_list_of_those)) {
    grouped_mode <- TRUE
    emg_force_fun <- calc_emg_force_rel_summary
  } else {
    grouped_mode <- FALSE
    emg_force_fun <- calc_emg_force_rel
  }
  
  emg_force_rels <- emg_force_fun(surface_potentials_or_list_of_those,
                                               stairs_start = stairs_start,
                                               step_length = step_length,
                                               num_steps = num_steps,
                                               transient_length = transient_length,
                                               electrodes = electrodes,
                                               electrode_confs = electrode_confs)

  
  if (!is.null(surface_potentials_or_list_of_those_2)) {
    emg_force_rels$setup <- 1
    emg_force_rels_2 <- 
      emg_force_fun(surface_potentials_or_list_of_those_2,
                    stairs_start = stairs_start,
                    step_length = step_length,
                    num_steps = num_steps,
                    transient_length = transient_length,
                    electrodes = electrodes,
                    electrode_confs = electrode_confs)
    emg_force_rels_2$setup <- 2
    emg_force_rels <- rbind(emg_force_rels, emg_force_rels_2)
    emg_force_rels <- emg_force_rels %>% dplyr::mutate(setup = as.factor(setup))
  }
  
  if (!for_latex) {
    y_label <- 'Normalized MAV[EMG]'
    x_label <- 'Normalized Muscle Force'
  } else {
    y_label <- 'Normalized MAV[EMG]'
    x_label <- 'Normalized Muscle Force $\\tilde\\tau$'
  }
  
  if (is.null(surface_potentials_or_list_of_those_2))
    plot_obj <- ggplot(data = emg_force_rels,
                       aes(x = force_mean, y = emg_mav)) + geom_line()
  else {
    plot_obj <- ggplot(data = emg_force_rels,
                       aes(x = force_mean, y = emg_mav, group=setup))
    plot_obj <- plot_obj + geom_line(aes(linetype = setup, group = setup))
  }
  
  if (grouped_mode && errorbars)
      plot_obj <- plot_obj + geom_errorbar(aes(ymin=emg_mav-emg_mav_sd, ymax=emg_mav+emg_mav_sd), width=0.05)
  
  plot_obj <- plot_obj + coord_cartesian(xlim = c(0, 1), ylim=c(0, 1)) + 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    ylab(y_label) + xlab(x_label)
  
  if (nlevels(emg_force_rels$electrode) == 1) {
    
    plot_obj <- plot_obj +
      ggtitle("Measured EMG-Force Relationship")
    
  } else {
    
    plot_obj <- plot_obj +
      facet_wrap(~ electrode, nrow = nrow) + 
      ggtitle("Measured EMG-Force Relationships")
    
  }
  
  if (print_plot) {
    if (is_mode("interactive")) dev.new()
    print(plot_obj)
  }
  
  invisible(plot_obj)
}

#' @import ggplot2
#' @importFrom plyr "."
#' @importFrom ggpubr ggarrange
#' @export
plot_electromechanical_properties <- function(MUs, muscles = NA, print_plot = TRUE) {
  
  if (!is.numeric(muscles))
    muscles <- unique(MUs$muscle)

  MUs <- MUs %>% dplyr::filter(muscle %in% muscles)  
  
  extract_property <- function(object_lst, property) {
    sapply(object_lst, function(obj) obj[[property]])
  }
  
  MU_props <- MUs %>% dplyr::mutate(num_fibers=extract_property(MU.obj, "num_fibers"),
                                    peak_force=extract_property(MU.obj, "peak_force"),
                                    rec_thresh=extract_property(MU.obj, "rec_thresh"),
                                    muscle=as.factor(muscle))
  
  plot_obj_1 <- ggplot(data = MU_props,
                       aes(x = peak_force, y = num_fibers)) +
    geom_point() +
    xlab("Peak MU Twitch Force") +
    ylab("Number of Muscle Fibers in MU")
  
  plot_obj_2 <- ggplot(data = MU_props,
                       aes(x = peak_force, y = rec_thresh)) +
    geom_point() +
    xlab("Peak MU Twitch Force") +
    ylab("MU Recruitment Threshold")
  
  if (nlevels(MU_props$muscle) > 1) {
    plot_obj_1 <- plot_obj_1 + facet_grid(muscle ~ .) +
      ggtitle("Electromechanical MU Properties of each Muscle")
    plot_obj_2 <- plot_obj_2 + facet_grid(muscle ~ .) +
      ggtitle("Force and Recruitment Properties of each Muscle")
    
  } else {
    plot_obj_1 <- plot_obj_1 + ggtitle("Electromechanical MU Properties")
    plot_obj_2 <- plot_obj_2 + ggtitle("Force and Recruitment Properties")
  }
  
  plots_combined = ggarrange(plot_obj_1, plot_obj_2, ncol=2)

  if (print_plot) {
    if (is_mode("interactive")) dev.new()
    print(plots_combined)
  }
  
  invisible(plots_combined)  
}


#' @import ggplot2
#' @importFrom ggpubr ggarrange
#' @export
plot_force_variability <- function(muscle_contribs_or_list_of_those, stairs_start, step_length, num_steps,
                                   transient_length, muscles = NA, nrow = 1, plot_mode=3,
                                   print_plot = TRUE, for_latex=FALSE, ...) {
  
  stopifnot(is.nice.scalar(nrow),
            is.list(muscle_contribs_or_list_of_those),
            plot_mode %in% c(1,2,3))
  
  process_force_contribs <- function(force_contribs) {
    if (!is.numeric(muscles))
      muscles <- unique(force_contribs$muscle)
    
    force_contribs <- force_contribs %>%
      dplyr::filter(muscle %in% muscles)
    
    force_variability <- calc_force_variability(force_contribs, stairs_start, step_length,
                                                num_steps, transient_length)
    force_variability
  }
  
  if (!is.data.frame(muscle_contribs_or_list_of_those$force)) {
    stopifnot(is.data.frame(muscle_contribs_or_list_of_those[[1]]$force))
    grouped_mode <- TRUE
    force_variability_lst <- list()
    for (ii in 1:length(muscle_contribs_or_list_of_those))
      force_variability_lst[[ii]] <- 
      process_force_contribs(muscle_contribs_or_list_of_those[[ii]]$force)
    force_variability <- dplyr::bind_rows(force_variability_lst)
    force_variability <- force_variability %>% dplyr::group_by(step) %>%
      dplyr::summarise_all(list("mean", "sd"))
    force_variability <- force_variability %>% dplyr::rename(force_mean = force_mean_mean)
    force_variability <- force_variability %>% dplyr::rename(force_cv = force_cv_mean)
    force_variability <- force_variability %>% dplyr::rename(force_sd = force_sd_mean)
    
  } else {
    
    grouped_mode <- FALSE
    force_variability <- process_force_contribs(muscle_contribs_or_list_of_those$force)  
  }
  
  if (!for_latex) {
    y_label_left <- 'CV of Muscle Force'
    y_label_right <- 'Normalized SD of Muscle Force '
    x_label <- 'Normalized Muscle Force'
  } else {
    y_label_left <- 'CV$\\,[\\tau]$'
    y_label_right <- 'Normalized SD$\\,[\\tau]$'
    x_label <- 'Normalized Muscle Force $\\tilde\\tau$'
  }
  
  plot_obj_left <-
    ggplot(data = force_variability,
           aes(x = force_mean, y = force_cv)) +
    geom_line() +
    coord_cartesian(xlim = c(0, 1)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    ylab(y_label_left) +
    xlab(x_label)
  
  plot_obj_right <-
    ggplot(data = force_variability,
           aes(x = force_mean, y = force_sd)) +
    geom_line() +
    coord_cartesian(xlim = c(0, 1)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    ylab(y_label_right) +
    xlab(x_label)
  
  if (grouped_mode) {
    plot_obj_left <- plot_obj_left + geom_errorbar(aes(ymin=force_cv-force_cv_sd,
                                                       ymax=force_cv+force_cv_sd))
    plot_obj_right <- plot_obj_right + geom_errorbar(aes(ymin=force_sd-force_sd_sd,
                                                       ymax=force_sd+force_sd_sd))
  }
  
  if (length(muscles) == 1) {
    plot_obj_left <- plot_obj_left + ggtitle("Coefficient of Force Variation")
    plot_obj_right <- plot_obj_right + ggtitle("Standard Deviation of Force")
    
  } else {
    
    plot_obj_left <- plot_obj_left +
      facet_grid(muscle ~ .) + 
      ggtitle("Coefficient of Force Variation")
    plot_obj_right <- plot_obj_right +
      facet_grid(muscle ~ .) + 
      ggtitle("Standard Deviation of Force")
  }
  
  plots_combined = ggarrange(plot_obj_left, plot_obj_right, ncol=2)
  
  if (plot_mode == 1) {
    # Only generate CoV plot
    if (print_plot) {
      if (is_mode("interactive")) dev.new()
      print(plot_obj_left)
    }
    
    invisible(plot_obj_left)  
    
  } else if (plot_mode == 2) {
    # Only generate SD plot
    if (print_plot) {
      if (is_mode("interactive")) dev.new()
      print(plot_obj_right)
    }
  
    invisible(plot_obj_right)
  } else {
    if (print_plot) {
      if (is_mode("interactive")) dev.new()
      print(plots_combined)
    }
    
    invisible(plots_combined) 
  }
}


#' @import ggplot2
#' @export
plot_force_twitches <- function(MUs, MUs_to_plot, muscle=1, print_plot=TRUE, for_latex=FALSE) {
  
  muscle_to_plot <- muscle
  MUs <- dplyr::filter(MUs, muscle == muscle_to_plot)
  
  force_twitches_lst <- list()
  t <- seq(0, 1, 0.002)
  for (ii in seq(1, length(MUs_to_plot))) {
    MU_id <- MUs_to_plot[ii]
    force_twitch_vals <- MUs$MU.obj[[MU_id]]$force_twitch(t)
    force_twitches_lst[[ii]] <- data.frame(time=t, twitch=force_twitch_vals, MU=MU_id)
  }
  force_twitches <- dplyr::bind_rows(force_twitches_lst) %>% dplyr::mutate(MU = as.factor(MU))
  
  if (!for_latex) {
    y_label <- 'MU Twitch Force'
    x_label <- 'Time (ms)'
  } else {
    y_label <- '$f_i$ ($k_f$ $\\cdot$ N)'
    x_label <- 'Time $t$ (ms)'
  }
  
  plot_obj <-
    ggplot(data = force_twitches,
           aes(x = time*1000, y = twitch, colour = MU, group=MU)) +
    geom_line(aes(group=MU)) +
    coord_cartesian(xlim = c(0, 1000)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    ylab(y_label) + xlab(x_label) + ggtitle("MU force twitches")

  if (print_plot) {
    if (is_mode("interactive")) dev.new()
    print(plot_obj)
  }
  
  invisible(plot_obj)  
}


#' @import ggplot2
#' @export
plot_exp_rel <- function(MUs, muscle=1, print_plot=TRUE, for_latex=FALSE) {
  
  muscle_to_plot <- muscle
  MUs <- dplyr::filter(MUs, muscle == muscle_to_plot)
  
  max_twitch_forces_lst <- list()
  t <- seq(0, 1, 0.002)
  for (ii in seq(1, nrow(MUs))) {
    max_twitch_forces_lst[[ii]] <- data.frame(id=MUs$MU[[ii]], max_force=MUs$MU.obj[[ii]]$peak_force)
  }
  max_twitch_forces <- dplyr::bind_rows(max_twitch_forces_lst)
  
  plot_obj <-
    ggplot(data = max_twitch_forces,
           aes(x = id, y = max_force)) +
    geom_point() +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    ylab('Peak twitch force') + xlab('MU ID') + ggtitle("Peak MU forces")
  
  if (print_plot) {
    if (is_mode("interactive")) dev.new()
    print(plot_obj)
  }
  
  invisible(plot_obj)  
}