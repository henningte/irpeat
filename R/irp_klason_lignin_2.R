#' Predicts the Klason lignin content from mid infrared spectra
#'
#' `irp_klason_lignin_content_2` predicts the Klason lignin content from mid
#' infrared spectra for peat samples using the model
#' [irpeatmodels::model_klason_lignin_content_2]. This function may also work
#' for organic matter in general. Note that this is a preliminary model only
#' which has not been fully validated for peat samples yet and which has known
#' limitations in predicting contents for peat samples
#' \insertCite{Teickner.2022a}{irpeat}.
#'
#' @inheritParams irp_eac_1
#'
#' @param x An object of class [`ir`][ir::ir_new_ir] with FTIR-ATR mid
#' infrared spectra. Some tests are applied to check if the supplied spectra
#' match the spectra used to fit the models (the spectral range is checked). The
#' spectral resolution of the original spectral data should not be smaller than
#' 4 cm\eqn{^{-1}} and it is not checked if this assumption is met.
#'
#' @note Note that this is a preliminary model only which has not been fully
#' validated for peat samples yet and which has known limitations in predicting
#' contents for peat samples \insertCite{Teickner.2022a}{irpeat}.
#'
#' @return `x` with a new column "klason_lignin_content_2" with the predicted
#' Klason lignin contents \[g/g\].
#'
#' @examples
#' library(ir)
#'
#' irp_klason_lignin_content_2(ir::ir_sample_data[1, ], do_summary = TRUE)
#'
#' @references
#'   \insertAllCited{}
#'
#' @export
irp_klason_lignin_content_2 <- function(x, ..., do_summary = FALSE, summary_function_mean = mean, summary_function_sd = stats::sd) {

  check_irpeatmodels(version = "0.0.0")
  if(! requireNamespace("brms", quietly = TRUE)) {
    rlang::abort("You have to install the 'brms' package to use this function.")
  }
  stopifnot(inherits(x, "ir"))
  stopifnot(is.logical(do_summary) && length(do_summary) == 1)

  x_or <- x

  # get data
  m <- irpeatmodels::model_klason_lignin_content_2
  config <- irpeatmodels::model_klason_lignin_content_2_config

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
  res <- as.data.frame(brms::posterior_predict(m, newdata = data.frame(x = I(as.matrix(x)), stringsAsFactors = FALSE), ...))
  res <- res * config$data_scale$y_scale + config$data_scale$y_center

  # summarize and add unit
  res <-
    irp_summarize_predictions(
      x = res,
      x_unit = "g/g",
      do_summary = do_summary,
      summary_function_mean = mean,
      summary_function_sd = stats::sd
    )

  x_or$klason_lignin_content_2 <- res
  x_or

}

