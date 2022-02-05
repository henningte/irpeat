#' Predicts the electron donating capacity from mid infrared spectra.
#'
#' \code{irp_edc_1} predicts the electron accepting capacity (EDC) from mid
#' infrared spectra of the peat samples. This function may also work for
#' organic matter in general \insertCite{Teickner.2022}{irpeat}.
#'
#' @inheritParams irp_eac_1
#' @return \code{x} with a new column "edc" with the predicted EDC values
#' [\eqn{\mu}mol g\eqn{_\text{C}^{-1}}].
#' @note The model still has a relatively large uncertainty because it is fitted
#' with few samples. Moreover, the model is known to produce biased predictions
#' \insertCite{Teickner.2022}{irpeat}. For further limitations, see
#' \insertCite{Teickner.2022;textual}{irpeat}.
#' @source \insertCite{Teickner.2022;textual}{irpeat}
#' @seealso \code{\link{model_edc_1}}.
#' @examples
#' \dontrun{
#' # get sample data
#' x <- ir::ir_sample_data
#'
#' # make predictions
#' x <- irpeat::irp_edc_1(x, do_summary = TRUE)
#' }
#' @references
#'   \insertAllCited{}
#' @export
irp_edc_1 <- function(x,
                    ...,
                    do_summary = FALSE)
{

  ir::ir_check_ir(x)
  stopifnot(is.logical(do_summary) && length(do_summary) == 1)

  x_or <- x

  # get data
  model_edc_1 <- NULL
  model_edc_1_config <- NULL
  utils::data("model_edc_1", envir = environment())
  utils::data("model_edc_1_config", envir = environment())
  m <- model_edc_1
  config <- model_edc_1_config

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
  res <- as.data.frame(rstanarm::posterior_predict(m, newdata = data.frame(x = I(as.matrix(x)), stringsAsFactors = FALSE), ...))
  res <- res * config$data_scale$y_scale + config$data_scale$y_center

  if(do_summary) {
    res <- quantities::set_quantities(purrr::map_dbl(res, mean),
                                      unit = "umol/g",
                                      errors = purrr::map_dbl(res, stats::sd))
  } else {
    res <- as.list(res)
  }

  x_or$edc <- res
  x_or

}

