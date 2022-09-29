#' Predicts peat properties from transmission mid infrared spectra
#'
#' Functions to predict peat properties from transmission mid infrared spectra.
#' All functions below have been computed using peat samples. For detailed
#' information on the underlying prediction models, see the details section.
#'
#' @include utils.R
#'
#' @name irp-predict-transmission-mir
#'
#' @details
#' The models use the models of the same name in the 'irpeatmodels' package. The
#' 'irpeatmodels' package provides information on the models.
#'
#' @param x An object of class [`ir`][ir::ir_new_ir] with transmission mid
#' infrared spectra. Some tests are applied to check if the supplied spectra
#' match the spectra used to fit the models (the spectral range is checked). The
#' spectral resolution of the original spectral data should not be smaller than
#' 4 cm\eqn{^{-1}} and it is not checked if this assumption is met.
#'
#' @param ... Additional arguments passed to
#' [rstanarm::posterior_predict.stanreg()] (`irp_eac_1()`,`irp_eac_2()`).
#'
#' @param do_summary A logical value indicating if the predicted values should
#' be returned in a summarized version (`TRUE`) or not (`FALSE`).
#' \itemize{
#'   \item If `do_summary = FALSE`, a list column is returned and each
#'   element of the list column is a numeric vector with draws from the
#'   posterior distribution, including the residual variance of the model.
#'   \item If `do_summary = TRUE`, each element is a
#'   [quantities::quantities()] object with the
#'   `error` attribute being the standard deviation of the unsummarized
#'   values.
#' }
#'
#' @param summary_function_mean A function used to summarize the predicted
#' values (average).
#'
#' @param summary_function_sd A function used to summarize the predicted
#' values (spread).
#'
#' @param check_prediction_domain A character value indicating if and how it
#' should be checked whether the spectra in `x` are within the prediction domain
#' of the model. One of:
#' \describe{
#'   \item{`"train"`}{It is checked whether the spectra in `x` are within the
#'   prediction domain formed by the training data for the model.}
#'   \item{`"test"`}{It is checked whether the spectra in `x` are within the
#'   prediction domain formed by the testing data for the model.}
#'   \item{`"none"`}{It is not checked whether the spectra in `x` are within the
#'   prediction domain for the model.}
#' }
#' Note: Currently ignored for `irp_eac_1()`, `irp_edc_1()`.
#'
#' @return `x` with a new column with the predicted peat property and a new
#' column with value `TRUE` if the respective spectrum is within the prediction
#' domain for the model and `FALSE` if not. If
#' `check_prediction_domain = "none"`, all values in this column are `NA`.
#'
#' @note
#' \describe{
#'   \item{`irp_eac_1()`, `irp_edc_1()`}{
#'     The model still has a relatively large uncertainty because it is fitted
#'     with few samples \insertCite{Teickner.2022}{irpeat}. For further
#'     limitations, see \insertCite{Teickner.2022;textual}{irpeat}.
#'   }
#' }
#'
#' @source
#' \describe{
#'   \item{`irp_eac_1()`, `irp_edc_1()`}{
#'     \insertCite{Teickner.2022;textual}{irpeat}.
#'   }
#' }
#'
#' @references
#'   \insertAllCited{}
#'
NULL


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' library(ir)
#'
#' x <- ir::ir_sample_data[1, ]
#'
#' ## make predictions
#'
#' # electron accepting capacity
#' x <- irpeat::irp_eac_1(
#'   x,
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_eac_1 <- function(x, ..., do_summary = FALSE, summary_function_mean = mean, summary_function_sd = stats::sd, check_prediction_domain = "train") {

  check_irpeatmodels(version = "0.0.0")
  if(! requireNamespace("rstanarm", quietly = TRUE)) {
    rlang::abort("You have to install the 'rstanarm' package to use this function.")
  }
  stopifnot(inherits(x, "ir"))
  stopifnot(is.logical(do_summary) && length(do_summary) == 1)

  x_or <- x

  # get data
  m <- irpeatmodels::model_eac_1
  config <- irpeatmodels::model_eac_1_config

  # check spectra
  x_flat <- ir::ir_flatten(x)
  if(x_flat$x[[1]] > config$irp_preprocess$clip_range$start) {
    rlang::warn(paste0("The minimum wavenumber value in `x` is ", x_flat$x[[1]], " , but should be ", config$irp_preprocess$clip_range$start, " or smaller."))
  }
  if(x_flat$x[[nrow(x_flat)]] < config$irp_preprocess$clip_range$end) {
    rlang::warn(paste0("The maximum wavenumber value in `x` is ", x_flat$x[[nrow(x)]], " , but should be ", config$irp_preprocess$clip_range$end, " or larger."))
  }

  # preprocess the spectra
  x <-
    irp_preprocess(
      x,
      do_interpolate = config$irp_preprocess$do_interpolate,
      interpolate_start = config$irp_preprocess$interpolate_start,
      interpolate_dw = config$irp_preprocess$interpolate_dw,
      do_clip = config$irp_preprocess$do_clip,
      clip_range = config$irp_preprocess$clip_range,
      do_interpolate_region = config$irp_preprocess$do_interpolate_region,
      interpolate_region_range = config$irp_preprocess$interpolate_region_range,
      do_bc = config$irp_preprocess$do_bc,
      bc_method = config$irp_preprocess$bc_method,
      bc_cutoff = config$irp_preprocess$bc_cutoff,
      bc_do_impute = config$irp_preprocess$bc_do_impute,
      do_smooth = config$irp_preprocess$do_smooth,
      do_normalise = config$irp_preprocess$do_normalise,
      normalise_method = config$irp_preprocess$normalise_method,
      do_bin = config$irp_preprocess$do_bin,
      bin_width = config$irp_preprocess$bin_width,
      bin_new_x_type = config$irp_preprocess$bin_new_x_type,
      do_scale = config$irp_preprocess$do_scale,
      scale_center = config$data_scale$x_center,
      scale_scale = config$data_scale$x_scale,
      do_return_as_ir = FALSE
    )

  # predict
  res <- as.data.frame(rstanarm::posterior_predict(m, newdata = data.frame(x = I(as.matrix(x)), stringsAsFactors = FALSE), ...))
  res <- res * config$data_scale$y_scale + config$data_scale$y_center

  # summarize and add unit
  res <-
    irp_summarize_predictions(
      x = res,
      x_unit = "umol/g",
      do_summary = do_summary,
      summary_function_mean = mean,
      summary_function_sd = stats::sd
    )

  x_or$eac <- res
  x_or

}


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # electron donating capacity
#' x <- irpeat::irp_edc_1(
#'   x,
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_edc_1 <- function(x, ..., do_summary = FALSE, summary_function_mean = mean, summary_function_sd = stats::sd, check_prediction_domain = "train") {

  check_irpeatmodels(version = "0.0.0")
  if(! requireNamespace("rstanarm", quietly = TRUE)) {
    rlang::abort("You have to install the 'rstanarm' package to use this function.")
  }
  stopifnot(inherits(x, "ir"))
  stopifnot(is.logical(do_summary) && length(do_summary) == 1)

  x_or <- x

  # get data
  m <- irpeatmodels::model_edc_1
  config <- irpeatmodels::model_edc_1_config

  # check spectra
  x_flat <- ir::ir_flatten(x)
  if(x_flat$x[[1]] > config$irp_preprocess$clip_range$start) {
    rlang::warn(paste0("The minimum wavenumber value in `x` is ", x_flat$x[[1]], " , but should be ", config$irp_preprocess$clip_range$start, " or smaller."))
  }
  if(x_flat$x[[nrow(x_flat)]] < config$irp_preprocess$clip_range$end) {
    rlang::warn(paste0("The maximum wavenumber value in `x` is ", x_flat$x[[nrow(x)]], " , but should be ", config$irp_preprocess$clip_range$end, " or larger."))
  }

  # preprocess the spectra
  x <-
    irp_preprocess(
      x,
      do_interpolate = config$irp_preprocess$do_interpolate,
      interpolate_start = config$irp_preprocess$interpolate_start,
      interpolate_dw = config$irp_preprocess$interpolate_dw,
      do_clip = config$irp_preprocess$do_clip,
      clip_range = config$irp_preprocess$clip_range,
      do_interpolate_region = config$irp_preprocess$do_interpolate_region,
      interpolate_region_range = config$irp_preprocess$interpolate_region_range,
      do_bc = config$irp_preprocess$do_bc,
      bc_method = config$irp_preprocess$bc_method,
      bc_cutoff = config$irp_preprocess$bc_cutoff,
      bc_do_impute = config$irp_preprocess$bc_do_impute,
      do_smooth = config$irp_preprocess$do_smooth,
      do_normalise = config$irp_preprocess$do_normalise,
      normalise_method = config$irp_preprocess$normalise_method,
      do_bin = config$irp_preprocess$do_bin,
      bin_width = config$irp_preprocess$bin_width,
      bin_new_x_type = config$irp_preprocess$bin_new_x_type,
      do_scale = config$irp_preprocess$do_scale,
      scale_center = config$data_scale$x_center,
      scale_scale = config$data_scale$x_scale,
      do_return_as_ir = FALSE
    )

  # predict
  res <- as.data.frame(rstanarm::posterior_predict(m, newdata = data.frame(x = I(as.matrix(x)), stringsAsFactors = FALSE), ...))
  res <- res * config$data_scale$y_scale + config$data_scale$y_center

  # summarize and add unit
  res <-
    irp_summarize_predictions(
      x = res,
      x_unit = "umol/g",
      do_summary = do_summary,
      summary_function_mean = mean,
      summary_function_sd = stats::sd
    )

  x_or$edc <- res
  x_or

}




#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # carbon content
#' x <- irpeat::irp_carbon_content_1(
#'   x,
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_carbon_content_1 <-
  irp_function_factory_eb1079(
    m =
      irpeatmodels::model_carbon_content_1_draws,
    m_pls =
      irpeatmodels::model_carbon_content_1_pls,
    config =
      irpeatmodels::model_carbon_content_1_config,
    prediction_domain =
      irpeatmodels::model_carbon_content_1_prediction_domain,
    target_variable_name =
      "carbon_content_1",
    x_unit =
      "g/g",
    irpeatmodels_required_version =
      "0.0.0",
    .f_check_packages =
      check_irpeatmodels_and_pls,
    .f_dimreduce =
      irp_make_predictions_plsr,
    .f_predict =
      irp_mcmc_predictions_beta_logit
  )

#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # nitrogen content
#' irpeat::irp_nitrogen_content_1(
#'   ir::ir_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_nitrogen_content_1 <-
  irp_function_factory_eb1079(
    m =
      irpeatmodels::model_nitrogen_content_1_draws,
    m_pls =
      irpeatmodels::model_nitrogen_content_1_pls,
    config =
      irpeatmodels::model_nitrogen_content_1_config,
    prediction_domain =
      irpeatmodels::model_nitrogen_content_1_prediction_domain,
    target_variable_name =
      "nitrogen_content_1",
    x_unit =
      "g/g",
    irpeatmodels_required_version =
      "0.0.0",
    .f_check_packages =
      check_irpeatmodels_and_pls,
    .f_dimreduce =
      irp_make_predictions_plsr,
    .f_predict =
      irp_mcmc_predictions_beta_logit
  )

#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # hydrogen content
#' irpeat::irp_hydrogen_content_1(
#'   ir::ir_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_hydrogen_content_1 <-
  irp_function_factory_eb1079(
    m =
      irpeatmodels::model_hydrogen_content_1_draws,
    m_pls =
      irpeatmodels::model_hydrogen_content_1_pls,
    config =
      irpeatmodels::model_hydrogen_content_1_config,
    prediction_domain =
      irpeatmodels::model_hydrogen_content_1_prediction_domain,
    target_variable_name =
      "hydrogen_content_1",
    x_unit =
      "g/g",
    irpeatmodels_required_version =
      "0.0.0",
    .f_check_packages =
      check_irpeatmodels_and_pls,
    .f_dimreduce =
      irp_make_predictions_plsr,
    .f_predict =
      irp_mcmc_predictions_beta_logit
  )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # oxygen content
#' irpeat::irp_oxygen_content_1(
#'   ir::ir_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_oxygen_content_1 <-
  irp_function_factory_eb1079(
    m =
      irpeatmodels::model_oxygen_content_1_draws,
    m_pls =
      irpeatmodels::model_oxygen_content_1_pls,
    config =
      irpeatmodels::model_oxygen_content_1_config,
    prediction_domain =
      irpeatmodels::model_oxygen_content_1_prediction_domain,
    target_variable_name =
      "oxygen_content_1",
    x_unit =
      "g/g",
    irpeatmodels_required_version =
      "0.0.0",
    .f_check_packages =
      check_irpeatmodels_and_pls,
    .f_dimreduce =
      irp_make_predictions_plsr,
    .f_predict =
      irp_mcmc_predictions_beta_logit
  )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # phosphorous content
#' irpeat::irp_phosphorous_content_1(
#'   ir::ir_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_phosphorous_content_1 <-
  irp_function_factory_eb1079(
    m =
      irpeatmodels::model_phosphorous_content_1_draws,
    m_pls =
      irpeatmodels::model_phosphorous_content_1_pls,
    config =
      irpeatmodels::model_phosphorous_content_1_config,
    prediction_domain =
      irpeatmodels::model_phosphorous_content_1_prediction_domain,
    target_variable_name =
      "phosphorous_content_1",
    x_unit =
      "g/g",
    irpeatmodels_required_version =
      "0.0.0",
    .f_check_packages =
      check_irpeatmodels_and_pls,
    .f_dimreduce =
      irp_make_predictions_plsr,
    .f_predict =
      irp_mcmc_predictions_beta_logit
  )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # potassium content
#' irpeat::irp_potassium_content_1(
#'   ir::ir_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_potassium_content_1 <-
  irp_function_factory_eb1079(
    m =
      irpeatmodels::model_potassium_content_1_draws,
    m_pls =
      irpeatmodels::model_potassium_content_1_pls,
    config =
      irpeatmodels::model_potassium_content_1_config,
    prediction_domain =
      irpeatmodels::model_potassium_content_1_prediction_domain,
    target_variable_name =
      "potassium_content_1",
    x_unit =
      "g/g",
    irpeatmodels_required_version =
      "0.0.0",
    .f_check_packages =
      check_irpeatmodels_and_pls,
    .f_dimreduce =
      irp_make_predictions_plsr,
    .f_predict =
      irp_mcmc_predictions_beta_logit
  )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # sulfur content
#' irpeat::irp_sulfur_content_1(
#'   ir::ir_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_sulfur_content_1 <-
  irp_function_factory_eb1079(
    m =
      irpeatmodels::model_sulfur_content_1_draws,
    m_pls =
      irpeatmodels::model_sulfur_content_1_pls,
    config =
      irpeatmodels::model_sulfur_content_1_config,
    prediction_domain =
      irpeatmodels::model_sulfur_content_1_prediction_domain,
    target_variable_name =
      "sulfur_content_1",
    x_unit =
      "g/g",
    irpeatmodels_required_version =
      "0.0.0",
    .f_check_packages =
      check_irpeatmodels_and_pls,
    .f_dimreduce =
      irp_make_predictions_plsr,
    .f_predict =
      irp_mcmc_predictions_beta_logit
  )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # titanium content
#' irpeat::irp_titanium_content_1(
#'   ir::ir_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_titanium_content_1 <-
  irp_function_factory_eb1079(
    m =
      irpeatmodels::model_titanium_content_1_draws,
    m_pls =
      irpeatmodels::model_titanium_content_1_pls,
    config =
      irpeatmodels::model_titanium_content_1_config,
    prediction_domain =
      irpeatmodels::model_titanium_content_1_prediction_domain,
    target_variable_name =
      "titanium_content_1",
    x_unit =
      "g/g",
    irpeatmodels_required_version =
      "0.0.0",
    .f_check_packages =
      check_irpeatmodels_and_pls,
    .f_dimreduce =
      irp_make_predictions_plsr,
    .f_predict =
      irp_mcmc_predictions_beta_logit
  )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # d13C values
#' irpeat::irp_d13C_1(
#'   ir::ir_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_d13C_1 <-
  irp_function_factory_eb1079(
    m =
      irpeatmodels::model_d13C_1_draws,
    m_pls =
      irpeatmodels::model_d13C_1_pls,
    config =
      irpeatmodels::model_d13C_1_config,
    prediction_domain =
      irpeatmodels::model_d13C_1_prediction_domain,
    target_variable_name =
      "d13C_1",
    x_unit =
      "1",
    irpeatmodels_required_version =
      "0.0.0",
    .f_check_packages =
      check_irpeatmodels_and_pls,
    .f_dimreduce =
      irp_make_predictions_plsr,
    .f_predict =
      irp_mcmc_predictions_beta_logit
  )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # d15N values
#' irpeat::irp_d15N_1(
#'   ir::ir_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_d15N_1 <-
  irp_function_factory_eb1079(
    m =
      irpeatmodels::model_d15N_1_draws,
    m_pls =
      irpeatmodels::model_d15N_1_pls,
    config =
      irpeatmodels::model_d15N_1_config,
    prediction_domain =
      irpeatmodels::model_d15N_1_prediction_domain,
    target_variable_name =
      "d15N_1",
    x_unit =
      "1",
    irpeatmodels_required_version =
      "0.0.0",
    .f_check_packages =
      check_irpeatmodels_and_pls,
    .f_dimreduce =
      irp_make_predictions_plsr,
    .f_predict =
      irp_mcmc_predictions_beta_logit
  )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # nominal oxidation state of carbon (NOSC)
#' irpeat::irp_nosc_1(
#'   ir::ir_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_nosc_1 <-
  irp_function_factory_eb1079(
    m =
      irpeatmodels::model_nosc_1_draws,
    m_pls =
      irpeatmodels::model_nosc_1_pls,
    config =
      irpeatmodels::model_nosc_1_config,
    prediction_domain =
      irpeatmodels::model_nosc_1_prediction_domain,
    target_variable_name =
      "nosc_1",
    x_unit =
      "mol/mol",
    irpeatmodels_required_version =
      "0.0.0",
    .f_check_packages =
      check_irpeatmodels_and_pls,
    .f_dimreduce =
      irp_make_predictions_plsr,
    .f_predict =
      irp_mcmc_predictions_beta_logit
  )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # Gibbs free energy of formation
#' irpeat::irp_dgf0_1(
#'   ir::ir_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_dgf0_1 <-
  irp_function_factory_eb1079(
    m =
      irpeatmodels::model_dgf0_1_draws,
    m_pls =
      irpeatmodels::model_dgf0_1_pls,
    config =
      irpeatmodels::model_dgf0_1_config,
    prediction_domain =
      irpeatmodels::model_dgf0_1_prediction_domain,
    target_variable_name =
      "dgf0_1",
    x_unit =
      "kJ/mol", #---todo: check
    irpeatmodels_required_version =
      "0.0.0",
    .f_check_packages =
      check_irpeatmodels_and_pls,
    .f_dimreduce =
      irp_make_predictions_plsr,
    .f_predict =
      irp_mcmc_predictions_beta_logit
  )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # dry bulk density
#' irpeat::irp_bulk_density_1(
#'   ir::ir_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_bulk_density_1 <-
  irp_function_factory_eb1079(
    m =
      irpeatmodels::model_bulk_density_1_draws,
    m_pls =
      irpeatmodels::model_bulk_density_1_pls,
    config =
      irpeatmodels::model_bulk_density_1_config,
    prediction_domain =
      irpeatmodels::model_bulk_density_1_prediction_domain,
    target_variable_name =
      "bulk_density_1",
    x_unit =
      "g/(cm^3)",
    irpeatmodels_required_version =
      "0.0.0",
    .f_check_packages =
      check_irpeatmodels_and_pls,
    .f_dimreduce =
      irp_make_predictions_plsr,
    .f_predict =
      irp_mcmc_predictions_beta_logit
  )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # O/C
#' irpeat::irp_O_to_C_1(
#'   ir::ir_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_O_to_C_1 <-
  irp_function_factory_eb1079(
    m =
      irpeatmodels::model_O_to_C_1_draws,
    m_pls =
      irpeatmodels::model_O_to_C_1_pls,
    config =
      irpeatmodels::model_O_to_C_1_config,
    prediction_domain =
      irpeatmodels::model_O_to_C_1_prediction_domain,
    target_variable_name =
      "O_to_C_1",
    x_unit =
      "g/g",
    irpeatmodels_required_version =
      "0.0.0",
    .f_check_packages =
      check_irpeatmodels_and_pls,
    .f_dimreduce =
      irp_make_predictions_plsr,
    .f_predict =
      irp_mcmc_predictions_beta_logit
  )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # H/C
#' irpeat::irp_H_to_C_1(
#'   ir::ir_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_H_to_C_1 <-
  irp_function_factory_eb1079(
    m =
      irpeatmodels::model_H_to_C_1_draws,
    m_pls =
      irpeatmodels::model_H_to_C_1_pls,
    config =
      irpeatmodels::model_H_to_C_1_config,
    prediction_domain =
      irpeatmodels::model_H_to_C_1_prediction_domain,
    target_variable_name =
      "H_to_C_1",
    x_unit =
      "g/g",
    irpeatmodels_required_version =
      "0.0.0",
    .f_check_packages =
      check_irpeatmodels_and_pls,
    .f_dimreduce =
      irp_make_predictions_plsr,
    .f_predict =
      irp_mcmc_predictions_beta_logit
  )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # C/N
#' irpeat::irp_C_to_N_1(
#'   ir::ir_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_C_to_N_1 <-
  irp_function_factory_eb1079(
    m =
      irpeatmodels::model_C_to_N_1_draws,
    m_pls =
      irpeatmodels::model_C_to_N_1_pls,
    config =
      irpeatmodels::model_C_to_N_1_config,
    prediction_domain =
      irpeatmodels::model_C_to_N_1_prediction_domain,
    target_variable_name =
      "C_to_N_1",
    x_unit =
      "g/g",
    irpeatmodels_required_version =
      "0.0.0",
    .f_check_packages =
      check_irpeatmodels_and_pls,
    .f_dimreduce =
      irp_make_predictions_plsr,
    .f_predict =
      irp_mcmc_predictions_beta_logit
  )


