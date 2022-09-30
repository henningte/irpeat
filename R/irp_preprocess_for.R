#' Preprocesses spectra as for the prediction of a peat property
#'
#' This function lets you extract the spectra after automated preprocessing for
#' the prediction of a specific peat property. This allows to inspect the
#' preprocessed spectra.
#'
#' @inheritParams irp_predict
#'
#' @return `x` where the spectra now a are the preprocessed spectra which would
#' be used in the next step to predict the peat property indicated by
#' `variable`.
#'
#' @examples
#' library(ir)
#'
#' irp_preprocess_for(
#'    x = irpeat_sample_data
#'    variable = "carbon_content_1"
#' ) %>%
#'   plot()
#'
#' @export
irp_preprocess_for <- function(x, variable) {

  stopifnot(inherits(x, "ir"))
  stopifnot(is.character(variable) && length(variable) == 1L)

  switch(variable,
         "klason_lignin_content_1" = ,
         "holocellulose_content_1" =
           x,
         "klason_lignin_content_2" =
           irp_klason_lignin_content_2(x = x, ...),
         "holocellulose_content_2" =
           irp_holocellulose_content_2(x = x, ...),
         "eac_1" =
           irp_eac_1(x = x, ...),
         "edc_1" =
           irp_edc_1(x = x, ...),
         "carbon_content_1" =,
         "nitrogen_content_1" =,
         "hydrogen_content_1" =,
         "oxygen_content_1" =,
         "phosphorous_content_1" =,
         "potassium_content_1" =,
         "sulfur_content_1" =,
         "titanium_content_1" =,
         "d13C_1" =,
         "d15N_1" =,
         "nosc_1" =,
         "dgf0_1" =,
         "bulk_density_1" =,
         "H_to_C_1" =,
         "O_to_C_1" =,
         "C_to_N_1" =,
         "volume_fraction_solids_1" =,
         "non_macroporosity_1" =,
         "macroporosity_1" =,
         "saturated_hydraulic_conductivity_1" = {
           check_irpeatmodels(version = "0.0.0")
           data(list = paste0("model_", variable, "_config"), package = "irpeatmodels", envir = environment())
           config <- get(x = paste0("model_", variable, "_config"), pos = -1)
           irp_preprocess_eb1079(x = x, config = config)
         },
         stop(paste0("Unknown value for `variable`: ", variable[[i]]))
  )

}

