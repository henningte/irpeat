#' General purpose preprocessing of spectra for predictions.
#'
#' \code{irp_preprocess} is a function that provides a general-purpose
#' preprocessing workflow for spectra. The workflow comprises interpolation,
#' clipping, baseline correction, smoothing (including optionally
#' derivatization), normalization, binning, and scaling. All these steps are
#' optionally, but occur on a fixed order that cannot be changed.
#'
#' @param x An object of class \code{\link[ir:ir_new_ir]{ir}} (it is assumed
#' that \code{x} is not yet preprocessed).
#' @param do_interpolate A logical value indicating if spectra should be
#' interpolated using \code{\link[ir:ir_interpolate]{ir_interpolate}}.
#' @param interpolate_start See \code{ir_interpolate} (parameter \code{start}).
#' @param interpolate_dw See \code{ir_interpolate} (parameter \code{dw}).
#' @param do_clip A logical value indicating if spectra should be
#' clipped using \code{\link[ir:ir_clip]{ir_clip}}.
#' @param clip_range See \code{ir_clip} (parameter \code{range}).
#' @param do_interpolate_region A logical value indicating if spectra should
#' be linearly interpolated in selected regions using
#' \code{\link[ir:ir_interpolate_region]{ir_interpolate_region}}.
#' @param interpolate_region_range See \code{ir_interpolate_region} (parameter
#' \code{range}).
#' @param do_bc A logical value indicating if spectra should be
#' baseline corrected using \code{\link[ir:ir_bc]{ir_bc}}.
#' @param bc_method See \code{ir_bc} (parameter \code{method}).
#' @param bc_degree See \code{ir_bc} (parameter \code{degree}).
#' @param bc_cutoff A numeric value representing the wavenumber
#' units to remove at the start and end of each spectrum in \code{x} during
#' baseline correction. This may be done to remove artifacts due to baseline
#' correction.
#' @param do_smooth A logical value indicating if spectra should be
#' smoothed using \code{\link[ir:ir_smooth]{ir_smooth}}.
#' @param smooth_method See \code{ir_smooth} (parameter \code{method}).
#' @param smooth_p See \code{ir_smooth} (parameter \code{p}).
#' @param smooth_n See \code{ir_smooth} (parameter \code{n}).
#' @param smooth_ts See \code{ir_smooth} (parameter \code{ts}).
#' @param smooth_m See \code{ir_smooth} (parameter \code{m}).
#' @param smooth_k See \code{ir_smooth} (parameter \code{k}).
#' @param do_normalise A logical value indicating if spectra should be
#' normalised using \code{\link[ir:ir_normalise]{ir_normalise}}.
#' @param normalise_method See \code{ir_normalise} (parameter \code{method}).
#' @param do_bin A logical value indicating if spectra should be
#' binned using \code{\link[ir:ir_bin]{ir_bin}}.
#' @param bin_width See \code{ir_bin} (parameter \code{width}).
#' @param do_scale A logical value indicating if spectral variables should be
#' scaled using \code{\link[base:scale]{scale}}.
#' @param scale_center See \code{scale} (parameter \code{center}). To scale each
#' spectral variable independently, provide a vector with length equal to the
#' number of spectral variables returned after preprocessing.
#' @param scale_scale See \code{scale} (parameter \code{scale}). To scale each
#' spectral variable independently, provide a vector with length equal to the
#' number of spectral variables returned after preprocessing.
#' @return A data frame with spectra in rows and a column for each spectral
#' variable after preprocessing.
#' @examples
#' # get sample data
#' x <- ir::ir_sample_data
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
#'     do_scale = TRUE,
#'     scale_center = TRUE,
#'     scale_scale = TRUE
#'   )
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
    x_flat <- ir::ir_flatten(x)
    bc_clip_range = data.frame(start = x_flat$x[[1]] + bc_cutoff,
                               end = x_flat$x[[nrow(x_flat)]] - bc_cutoff,
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
    x <- ir::ir_normalise(x,
                          method = normalise_method)
  }

  # binning
  if(do_bin) {
    x <- ir::ir_bin(x,
                    width = bin_width)
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
