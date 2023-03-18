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
#'    x = irpeat_sample_data,
#'    variable = "carbon_content_1"
#' ) %>%
#'   plot()
#'
#' @export
irp_preprocess_for <- function(x, y = NULL, variable) {

  stopifnot(inherits(x, "ir"))
  stopifnot(is.character(variable) && length(variable) == 1L)
  stopifnot(inherits(y, "ir")|| (is.null(y) && ! variable %in% c("microbial_nitrogen_content_1")))

  switch(variable,
         "klason_lignin_content_1" = ,
         "holocellulose_content_1" =
           x %>%
           irp_content_klh_hodgkins_preprocess(),
         "klason_lignin_content_2" = {
           check_irpeatmodels(version = "0.0.0")
           irp_preprocess_eb1014(x = x, config = irpeatmodels::model_klason_lignin_content_2_config)
         },
         "holocellulose_content_2" = {
           check_irpeatmodels(version = "0.0.0")
           irp_preprocess_eb1014(x = x, config = irpeatmodels::model_holocellulose_content_2_config)
         },
         "eac_1" = {
           check_irpeatmodels(version = "0.0.0")
           irp_preprocess_eb1014(x = x, config = irpeatmodels::model_eac_1_config)
         },
         "edc_1" = {
           check_irpeatmodels(version = "0.0.0")
           irp_preprocess_eb1014(x = x, config = irpeatmodels::model_edc_1_config)
         },
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
         "saturated_hydraulic_conductivity_1" =,
         "specific_heat_capacity_1" =,
         "dry_thermal_conductivity_1" = {
           check_irpeatmodels(version = "0.0.0")
           utils::data(list = paste0("model_", variable, "_config"), package = "irpeatmodels", envir = environment())
           config <- get(x = paste0("model_", variable, "_config"), pos = -1)
           irp_preprocess_eb1079(x = x, config = config)
         },
         "microbial_nitrogen_content_1" = {
           check_irpeatmodels(version = "0.0.0")
           utils::data(list = paste0("model_", variable, "_config"), package = "irpeatmodels", envir = environment())
           config <- get(x = paste0("model_", variable, "_config"), pos = -1)

           # preprocess both
           x <- irp_preprocess_unpack_config(x = x, config = config)
           y <- irp_preprocess_unpack_config(x = y, config = config)

           # make sure x matches wavenumber range in y
           x <- ir::ir_clip(x, range = data.frame(start = min(y$spectra[[1]]$x), end = max(y$spectra[[1]]$x)))

           # compute second derivative difference spectra
           ir::ir_subtract(x, y) %>%
             magrittr::multiply_by(-100) %>%
             ir::ir_smooth(method = "sg", p = 3, n = 65, ts = 1, m = 2) %>%
             ir::ir_get_intensity(wavenumber = 1220) %>%
             dplyr::mutate(
               spectra =
                 purrr::map(.data$intensity, function(.y) {
                   .y %>%
                     dplyr::mutate(
                       y =
                         y %>%
                         magrittr::multiply_by(config$data_scale$x_scale) %>%
                         magrittr::subtract(config$data_scale$x_center)
                     )
                 })
             )

         },
         stop(paste0("Unknown value for `variable`: ", variable))
  )

}

