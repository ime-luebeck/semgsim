#' @import ggplot2
#' @importFrom plyr "."
#' @export
compare_sim_runs <- function(sim_results, MU = 1, t_max = 0.03) {

    stopifnot(is.list(sim_results))
    
    transformer <- function(model_num) {

        NFFT <- sim_results[[model_num]]$sampling$NFFT
        times <- seq(from = 0,
                     to = (NFFT - 1) / sim_results[[model_num]]$Fs,
                     by = 1 / sim_results[[model_num]]$Fs)

        transform_model <- function(df_row) {
            with(df_row,
                 data.frame(
                     model = model_num %>% rep(NFFT) %>% as.character,
                     electrode = electrode[[1]] %>% rep(NFFT),
                     time = times,
                     response =
                       Re(fft(firing_response_fftvals[[1]] *
                                sim_results[[model_num]]$Fs / (NFFT-1),
                              inverse = TRUE)),
                     stringsAsFactors = FALSE))
        }
    }
    
    plot_data <- list()

    for (i in seq_along(sim_results))
        plot_data[[i]] <-
            plyr::ddply(.data = sim_results[[i]]$MU_firing_responses_fftvals %>%
                          subset(MU == MU),
                        .variables = .(electrode),
                        .fun = transformer(i))

    plot_data_TF <- dplyr::bind_rows(plot_data) %>%
      plyr::mutate(model = factor(model)) %>%
      subset(time < t_max)

    plot_obj <-
        ggplot(plot_data_TF, aes(x = time, y = response, colour = model)) +
          geom_line(aes(group = model), size = 1) +
          facet_grid(electrode ~ ., scales = "free") +
          ggtitle(paste("Firing response of MU", MU, "at different electrodes"))
    
    dev.new()
    print(plot_obj)
}
