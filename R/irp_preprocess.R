#' General purpose preprocessing of spectra for predictions
#'
#' `irp_preprocess` is a function that provides a general-purpose
#' preprocessing workflow for spectra. The workflow comprises interpolation,
#' clipping, baseline correction, smoothing (including optionally
#' derivatization), normalization, binning, and scaling. All these steps are
#' optionally, but occur on a fixed order that cannot be changed.
#'
#' @param x An object of class [`ir`][ir::ir_new_ir] (it is assumed
#' that `x` is not yet preprocessed).
#'
#' @param do_interpolate A logical value indicating if spectra should be
#' interpolated using [ir::ir_interpolate()].
#'
#' @param interpolate_start See `ir_interpolate` (parameter `start`).
#'
#' @param interpolate_dw See `ir_interpolate` (parameter `dw`).
#'
#' @param do_clip A logical value indicating if spectra should be
#' clipped using [ir::ir_clip()].
#'
#' @param clip_range See `ir_clip` (parameter `range`).
#'
#' @param do_interpolate_region A logical value indicating if spectra should
#' be linearly interpolated in selected regions using
#' [ir::ir_interpolate_region()].
#'
#' @param interpolate_region_range See `ir_interpolate_region` (parameter
#' `range`).
#'
#' @param do_bc A logical value indicating if spectra should be
#' baseline corrected using [ir::ir_bc()].
#'
#' @param bc_method See `ir_bc` (parameter `method`).
#'
#' @param bc_degree See `ir_bc` (parameter `degree`).
#'
#' @param bc_cutoff A numeric value representing the wavenumber
#' units to remove at the start and end of each spectrum in `x` during
#' baseline correction. This may be done to remove artifacts due to baseline
#' correction.
#'
#' @param do_smooth A logical value indicating if spectra should be
#' smoothed using [ir::ir_smooth()].
#'
#' @param smooth_method See `ir_smooth` (parameter `method`).
#'
#' @param smooth_p See `ir_smooth` (parameter `p`).
#'
#' @param smooth_n See `ir_smooth` (parameter `n`).
#'
#' @param smooth_ts See `ir_smooth` (parameter `ts`).
#'
#' @param smooth_m See `ir_smooth` (parameter `m`).
#'
#' @param smooth_k See `ir_smooth` (parameter `k`).
#'
#' @param do_normalise A logical value indicating if spectra should be
#' normalized using [ir::ir_normalize()].
#'
#' @param normalise_method See `ir_normalize` (parameter `method`).
#'
#' @param do_bin A logical value indicating if spectra should be
#' binned using [ir::ir_bin()].
#'
#' @param bin_width See `ir_bin` (parameter `width`).
#'
#' @param bin_new_x_type See `ir_bin` (parameter `new_x_type`).
#'
#' @param do_scale A logical value indicating if spectral variables should be
#' scaled using [base::scale()].
#'
#' @param scale_center See `scale` (parameter `center`). To scale each
#' spectral variable independently, provide a vector with length equal to the
#' number of spectral variables returned after preprocessing.
#'
#' @param scale_scale See `scale` (parameter `scale`). To scale each
#' spectral variable independently, provide a vector with length equal to the
#' number of spectral variables returned after preprocessing.
#'
#' @return A data frame with spectra in rows and a column for each spectral
#' variable after preprocessing.
#'
#' @examples
#' # get sample data
#' x <- ir::ir_sample_data[1, ]
#'
#' # example preprocessing
#' res <-
#'   irpeat::irp_preprocess(
#'     x,
#'     do_interpolate = TRUE,
#'     interpolate_start = NULL,
#'     interpolate_dw = 1,
#'     do_bc = TRUE,
#'     do_clip = FALSE,
#'     do_interpolate_region = FALSE,
#'     bc_method = "rubberband",
#'     bc_cutoff = 10,
#'     do_smooth = FALSE,
#'     do_normalise = TRUE,
#'     normalise_method = "area",
#'     do_bin = TRUE,
#'     bin_width = 10,
#'     bin_new_x_type = "start",
#'     do_scale = TRUE,
#'     scale_center = TRUE,
#'     scale_scale = TRUE
#'   )
#'
#' @export
irp_preprocess <- function(
  x,
  do_interpolate = TRUE,
  interpolate_start = NULL,
  interpolate_dw = 1,
  do_clip = TRUE,
  clip_range,
  do_interpolate_region = FALSE,
  interpolate_region_range,
  do_bc = TRUE,
  bc_method = "rubberband",
  bc_degree = 2,
  bc_cutoff = 0,
  do_smooth = TRUE,
  smooth_method = "sg",
  smooth_p = 3,
  smooth_n = smooth_p + 3 - smooth_p%%2,
  smooth_m = 1,
  smooth_ts = 0,
  smooth_k = 111,
  do_normalise = TRUE,
  normalise_method = "area",
  do_bin = TRUE,
  bin_width = 10,
  bin_new_x_type = "start",
  do_scale = TRUE,
  scale_center = TRUE,
  scale_scale = TRUE
) {

  # checks
  stopifnot(is.logical(do_interpolate)        && length(do_interpolate) == 1)
  stopifnot(is.logical(do_clip)               && length(do_clip) == 1)
  stopifnot(is.logical(do_interpolate_region) && length(do_interpolate_region) == 1)
  stopifnot(is.logical(do_bc)                 && length(do_bc) == 1)
  stopifnot(is.logical(do_smooth)             && length(do_smooth) == 1)
  stopifnot(is.logical(do_normalise)          && length(do_normalise) == 1)
  stopifnot(is.logical(do_bin)                && length(do_bin) == 1)
  stopifnot(is.logical(do_scale)              && length(do_scale) == 1)

  if(!is.numeric(bc_cutoff) && length(bc_cutoff) != 1 && do_bc) {
    rlang::abort("`bc_cutoff` must be a numeric value.")
  }


  # interpolation
  if(do_interpolate) {
    x <- ir::ir_interpolate(x,
                            start = interpolate_start,
                            dw = interpolate_dw)
  }

  # clipping
  if(do_clip) {
    x <- ir::ir_clip(x,
                     range = clip_range)
  }

  # regional interpolation
  if(do_interpolate_region) {
    x <- ir::ir_interpolate_region(x,
                                   range = interpolate_region_range)
  }

  # baseline correction
  if(do_bc) {
    x_spectra_x_range <-
      data.frame(
        start =
          min(purrr::map_dbl(x$spectra, function(.x) min(.x$x))),
        end =
          max(purrr::map_dbl(x$spectra, function(.x) max(.x$x)))
      )

    bc_clip_range = data.frame(start = x_spectra_x_range$start + bc_cutoff,
                               end = x_spectra_x_range$end - bc_cutoff,
                               stringsAsFactors = FALSE)
    x <- ir::ir_bc(x,
                   method = bc_method,
                   degree = bc_degree,
                   return_bl = FALSE)
    x <- ir::ir_clip(x,
                     range = bc_clip_range)
    x <- ir::ir_bc(x,
                   method = bc_method,
                   degree = bc_degree,
                   return_bl = FALSE)
  }

  # smoothing
  if(do_smooth) {
    x <- ir::ir_smooth(x,
                       method = smooth_method,
                       p = smooth_p,
                       n = smooth_n,
                       ts = smooth_ts,
                       m = smooth_m,
                       k = smooth_k)
  }

  # normalise
  if(do_normalise) {
    x <- ir::ir_normalize(x,
                          method = normalise_method)
  }

  # binning
  if(do_bin) {
    x <- ir::ir_bin(x, width = bin_width, new_x_type = bin_new_x_type)
  }

  res <- ir::ir_flatten(x)
  res_colnames <- paste0("V", res[, 1, drop= TRUE])
  res <- t(res[, -1, drop = FALSE])

  # scale
  if(do_scale) {
    res <- scale(res, center = scale_center, scale = scale_scale)
  }

  if(do_scale) {
    scale_center = attr(res, "scaled:center")
    scale_scale = attr(res, "scaled:scale")
  }
  res <- as.data.frame(res)
  if(do_scale) {
    attr(res, "scaled:center") <- scale_center
    attr(res, "scaled:scale") <- scale_scale
  }
  colnames(res) <- res_colnames
  res

}
