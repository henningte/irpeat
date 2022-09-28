#' Predicts the carbon content from mid infrared spectra
#'
#' @inheritParams irp_eac_1
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
#'
#' @return `x` with a new column `carbon_content_1` with the predicted carbon
#' contents \[g/g\] and a new column `carbon_content_1_in_pd` with value `TRUE`
#' if the respective spectrum is within the prediction domain for the model and
#' `FALSE` if not. If `check_prediction_domain = "none"`, all values in
#' `carbon_content_1_in_pd` are `NA`.
#'
#' @source ---todo
#'
#' @examples
#' library(ir)
#'
#' # make predictions
#' irpeat::irp_carbon_content_1(
#'   ir::ir_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_carbon_content_1 <- function(x, do_summary = FALSE, summary_function_mean = mean, summary_function_sd = stats::sd, check_prediction_domain = "train") {

  check_irpeatmodels(version = "0.0.0")
  if(! requireNamespace("pls", quietly = TRUE)) {
    rlang::abort("You have to install the 'pls' package to use this function.")
  }
  stopifnot(inherits(x, "ir"))
  stopifnot(is.logical(do_summary) && length(do_summary) == 1)

  x_or <- x

  # get data
  m <- irpeatmodels::model_carbon_content_1_draws
  m_pls <- irpeatmodels::model_carbon_content_1_pls
  config <- irpeatmodels::model_carbon_content_1_config
  prediction_domain <- irpeatmodels::model_carbon_content_1_prediction_domain

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
      # bc_cutoff = config$irp_preprocess$bc_cutoff,
      bc_cutoff = 0,
      # bc_do_impute = config$irp_preprocess$bc_do_impute,
      bc_do_impute = TRUE,
      do_smooth = config$irp_preprocess$do_smooth,
      do_normalise = config$irp_preprocess$do_normalise,
      normalise_method = config$irp_preprocess$normalise_method,
      do_bin = config$irp_preprocess$do_bin,
      bin_width = config$irp_preprocess$bin_width,
      bin_new_x_type = config$irp_preprocess$bin_new_x_type,
      do_scale = config$irp_preprocess$do_scale,
      scale_center = config$data_scale$x_center,
      scale_scale = config$data_scale$x_scale,
      do_return_as_ir = TRUE
    )

  # check prediction domain
  prediction_domain <-
    switch(
      check_prediction_domain,
      "train" = prediction_domain$train,
      "test" = prediction_domain$test,
      "none" = NULL
    )

  if(check_prediction_domain != "none") {
    x_or$carbon_content_1_in_pd <-
      x %>%
      irp_is_in_prediction_domain(prediction_domain = prediction_domain) %>%
      dplyr::pull(is_in_prediction_domain)
  }

  # get plsr scores
  res <-
    irp_make_predictions_plsr(
      x = x,
      m_pls = m_pls,
      config = config
    )

  # predict
  res <-
    irp_mcmc_predictions_beta_logit(
      x = res,
      draws = m,
      config = config
    )

  # summarize and add unit
  res <-
    irp_summarize_predictions(
      x = res,
      x_unit = "g/g",
      do_summary = do_summary,
      summary_function_mean = mean,
      summary_function_sd = stats::sd
    )

  x_or$carbon_content_1 <- res
  x_or

}
