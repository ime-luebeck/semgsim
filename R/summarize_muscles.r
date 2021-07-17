#' Summarize key properties of generated muscles to command line
#'
#' @importFrom plyr llply
#' @export
summarize_muscles <- function(sim_result, muscles=NA) {

    if (!is.numeric(muscles))
		muscles = unique(sim_result$MUs$muscle)

	for (i in muscles) {
		print(sprintf("--------- MUSCLE %.0f------------", i))
		MUs <- dplyr::filter(sim_result$MUs, muscle == i)
		print(sprintf("Number of MUs: %.0f", nrow(MUs)))
		MU_num_fibers = unlist(llply(MUs$MU.obj, function(MU.obj) MU.obj$num_fibers))
		print(sprintf("Number of fibers (total): %.0f", sum(MU_num_fibers)))
		print(sprintf("Number of fibers per MU (avg, min-max range): %.1f [%.0f - %.0f]", mean(MU_num_fibers), min(MU_num_fibers), max(MU_num_fibers)))
		
		MU_rec_threshs = unlist(llply(MUs$MU.obj, function(MU.obj) MU.obj$rec_thresh))
		print(sprintf("MU recruitment thresholds (avg, min-max range): %.2f [%.2f - %.2f]", mean(MU_rec_threshs), min(MU_rec_threshs), max(MU_rec_threshs)))
		
		MU_cutoff_ratios = unlist(llply(MUs$MU.obj, function(MU.obj) MU.obj$cutoff_ratio))
		print(sprintf("MU cut-off ratios (avg, min-max range): %.2f [%.2f - %.2f]", mean(MU_cutoff_ratios),min(MU_cutoff_ratios), max(MU_cutoff_ratios)))
	}
}