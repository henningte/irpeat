#' Predicts the holocellulose content from mid infrared spectra
#'
#' `irp_holocellulose_content_2` predicts the holocellulose content from mid
#' infrared spectra for peat samples using the model
#' [irpeatmodels::model_holocellulose_content_2]. This function may also work
#' for organic matter in general. Note that this is a preliminary model only
#' which has not been fully validated for peat samples yet and which has known
#' limitations in predicting contents for peat samples
#' \insertCite{Teickner.2022a}{irpeat}.
#'
#' @inheritParams irp_klason_lignin_content_2
#'
#' @note Note that this is a preliminary model only which has not been fully
#' validated for peat samples yet and which has known limitations in predicting
#' contents for peat samples \insertCite{Teickner.2022a}{irpeat}.
#'
#' @return `x` with a new column "holocellulose_content_2" with the predicted
#' holocellulose contents \[g/g\].
#'
#' @examples
#' library(ir)
#'
#' irp_holocellulose_content_2(ir::ir_sample_data[1, ], do_summary = TRUE)
#'
#' @references
#'   \insertAllCited{}
#'
#' @export
irp_holocellulose_content_2 <- function(x, ..., do_summary = FALSE, summary_function_mean = mean, summary_function_sd = stats::sd, check_prediction_domain = "train") {

  check_irpeatmodels(version = "0.0.0")
  if(! requireNamespace("brms", quietly = TRUE)) {
    rlang::abort("You have to install the 'brms' package to use this function.")
  }
  stopifnot(inherits(x, "ir"))
  stopifnot(is.logical(do_summary) && length(do_summary) == 1)

  x_or <- x

  # get data
  m <- irpeatmodels::model_holocellulose_content_2
  config <- irpeatmodels::model_holocellulose_content_2_config

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
      do_return_as_ir = TRUE
    )

  # check prediction domain
  prediction_domain <-
    switch(
      check_prediction_domain,
      "train" = irpeatmodels::model_holocellulose_content_2_prediction_domain$train,
      "test" = irpeatmodels::model_holocellulose_content_2_prediction_domain$train,
      "none" = NULL
    )

  res_pd <-
    if(check_prediction_domain != "none") {
      tibble::tibble(
        y =
          x %>%
          irp_is_in_prediction_domain(prediction_domain = prediction_domain) %>%
          dplyr::pull(.data$is_in_prediction_domain)
      )
    } else {
      tibble::tibble(
        y = rep(NA, nrow(x_or))
      )
    }

  # reformat for predictions
  res <- ir::ir_flatten(x)
  res_colnames <- paste0("V", res[, 1, drop= TRUE])
  res <- as.data.frame(t(res[, -1, drop = FALSE]))
  attr(res, "scaled:center") <- attr(x, "scaled:center")
  attr(res, "scaled:scale") <- attr(x, "scaled:scale")
  colnames(res) <- res_colnames

  # predict
  res <- as.data.frame(brms::posterior_predict(m, newdata = data.frame(x = I(as.matrix(res)), stringsAsFactors = FALSE), ...))
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

  x_or$holocellulose_content_2 <- res
  cbind(x_or, res_pd %>% stats::setNames(nm = "holocellulose_content_2_in_pd"))

}

