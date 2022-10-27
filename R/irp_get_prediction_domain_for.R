#' Preprocesses spectra as for the prediction of a peat property
#'
#' This function lets you extract the spectra after automated preprocessing for
#' the prediction of a specific peat property. This allows to inspect the
#' preprocessed spectra.
#'
#' @inheritParams irp_carbon_content_1
#' @inheritParams irp_predict
#'
#' @return An object of class `irp_prediction_domain` with the specified
#' prediction domain for the selected variable.
#'
#' @examples
#' library(ir)
#'
#' irp_get_prediction_domain_for(
#'    variable = "carbon_content_1"
#' ) %>%
#'   plot()
#'
#' irp_get_prediction_domain_for(
#'    variable = "klason_lignin_content_1"
#' ) %>%
#'   plot()
#'
#' irp_get_prediction_domain_for(
#'    variable = "klason_lignin_content_2"
#' ) %>%
#'   plot()
#'
#' irp_get_prediction_domain_for(
#'    variable = "eac_1"
#' ) %>%
#'   plot()
#'
#' @export
irp_get_prediction_domain_for <- function(variable, check_prediction_domain = "train") {

  stopifnot(is.character(variable) && length(variable) == 1L)

  # placeholder: ---todo: remove later
  empty_prediction_domain <-
    tibble::tibble(
      x = numeric(),
      ymin = numeric(),
      ymax = numeric()
    ) %>%
    irp_as_irp_prediction_domain()

  # point to parent variable where necessary
  variable <-
    switch(
      variable,
      "volume_fraction_solids_1" =,
      "non_macroporosity_1" =,
      "macroporosity_1" =,
      "saturated_hydraulic_conductivity_1" =
        "bulk_density_1",
      variable
    )

  switch(
    variable,
    "klason_lignin_content_1" =,
    "holocellulose_content_1" =,
    "klason_lignin_content_2" =,
    "holocellulose_content_2" =,
    "eac_1" =,
    "edc_1" =,
    "carbon_content_1" =,
    "nitrogen_content_1" =,
    "hydrogen_content_1" =,
    "oxygen_content_1" =,
    "phosphorus_content_1" =,
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
      utils::data(list = paste0("model_", variable, "_prediction_domain"), package = "irpeatmodels", envir = environment())
      prediction_domain <- get(x = paste0("model_", variable, "_prediction_domain"), pos = -1)
    },
    stop(paste0("Unknown value for `variable`: ", variable))
  )

  switch(
    check_prediction_domain,
    "train" = prediction_domain$train,
    "test" = prediction_domain$test,
    "none" = empty_prediction_domain
  )

}
