#' @export
split_semg_geometry <- function(model_geom, num_parts, file_name_prefix = "") {

    MUs_parts <- split(model_geom$MUs, (1:nrow(model_geom$MUs)) %% num_parts)

    if (file_name_prefix != "")
        for (i in seq_along(MUs_parts)) {
            model_part <- model_geom
            model_part$MUs <- MUs_parts[[i]]
            saveRDS(object = model_part,
                    file =
                      paste(file_name_prefix, as.character(i), ".RDS",
                            sep = ""))
        }

    MUs_parts
}

#' @export
merge_semg_geometry <- function(file_name_prefix) {


}
