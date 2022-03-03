#' Predicts the holocellulose content from mid infrared spectra.
#'
#' \code{irp_holocellulose_2} predicts the holocellulose content from mid
#' infrared spectra for peat samples using the model
#' \code{\link{model_holocellulose_2}}. This function may also work for organic matter
#' in general. Note that this is a preliminary model only which has not been
#' fully validated for peat samples yet and which has known limitations in
#' predicting contents for peat samples [--- todo: add reference].
#'
#' @inheritParams irp_klason_lignin_2
#'
#' @note Note that this is a preliminary model only which has not been fully
#' validated for peat samples yet and which has known limitations in predicting
#' contents for peat samples [--- todo: add reference].
#'
#' @return \code{x} with a new column "holocellulose_2" with the predicted
#' holocellulose contents [g/g].
#' @examples
#' \dontrun{
#' # get sample data
#' x <- ir::ir_sample_data
#'
#' # make predictions
#' x <- irpeat::irp_holocellulose_2(x, do_summary = FALSE)
#' }
#' @references
#'   \insertAllCited{}
#' @export
irp_holocellulose_2 <- function(x,
                    ...,
                    do_summary = FALSE)
{

  ir::ir_check_ir(x)
  stopifnot(is.logical(do_summary) && length(do_summary) == 1)

  x_or <- x

  # get data
  model_holocellulose_2 <- NULL
  model_holocellulose_2_config <- NULL
  utils::data("model_holocellulose_2", envir = environment())
  utils::data("model_holocellulose_2_config", envir = environment())
  m <- model_holocellulose_2
  config <- model_holocellulose_2_config

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
      do_smooth = config$irp_preprocess$do_smooth,
      do_normalise = config$irp_preprocess$do_normalise,
      normalise_method = config$irp_preprocess$normalise_method,
      do_bin = config$irp_preprocess$do_bin,
      bin_width = config$irp_preprocess$bin_width,
      do_scale = config$irp_preprocess$do_scale,
      scale_center = config$data_scale$x_center,
      scale_scale = config$data_scale$x_scale
    )

  # predict
  res <- as.data.frame(brms::posterior_predict(m, newdata = data.frame(x = I(as.matrix(x)), stringsAsFactors = FALSE), ...))
  res <- res * config$data_scale$y_scale + config$data_scale$y_center

  if(do_summary) {
    res <- quantities::set_quantities(purrr::map_dbl(res, mean),
                                      unit = "g/g",
                                      errors = purrr::map_dbl(res, stats::sd))
  } else {
    res <- as.list(res)
  }

  x_or$holocellulose_2 <- res
  x_or

}

