#' Compute humification indices from mid-infrared spectra
#'
#' `irp_hi` computes either custom humification indices or a predefined set of
#' humification indices as reported by \insertCite{Broder.2012;textual}{irpeat}
#' for mid-infrared spectra with the x values representing wavenumber values
#' (no checks are performed) (see the details section for the humification indices
#' computed by default). A humification index is the ratio of the intensity
#' values at two different x axis values (e.g. wavenumbers) and defined as:
#' \deqn{
#'   HI = \frac{y_x1}{y_x2}
#' },
#' where \eqn{y_x1} is the intensity of the spectrum at x axis value \eqn{x1} and
#' \eqn{y_x2} is the intensity of the spectrum at x axis value \eqn{x2}.
#'
#' @details
#' The following humification indices are computed by
#' default (if `x1 = NULL` and `x2 = NULL`) (values
#' represent wavenumbers \[cm\eqn{^{-1}}\]) \insertCite{Broder.2012}{irpeat} :
#' \describe{
#'   \item{hi1: 1420/1090}{OH and CO of phenols or CH of CH{_1} and CH{_3}
#'   groups (phenolic and aliphatic structures)/polysaccharides}
#'   \item{hi2: 1510/1090}{Aromatic C=C or C=O of amides/polysaccharides}
#'   \item{hi3: 1630/1090}{Aromatic C=C and COO\eqn{^{-1}} (aromatics and aromatic
#'   or aliphatic carboxylates)/polysaccharides}
#'   \item{hi4: 1720/1090}{carbonylic and carboxylic C=O (carboxylic acids and
#'   aromatic esters)/polysaccharides}
#' }
#'
#' @param x An object of class [`ir`][ir::ir_new_ir].
#'
#' @param x1 A numeric vector with values representing x axis values in the spectra
#' of `x`. This is \eqn{x1} in the equation displayed above. If
#' multiple humification indices for the same `x2` should be computed, `x1`
#' can contain multiple values. If `x1 = NULL`,
#' the default humification indices will be computed (see details section).
#'
#' @param x2 A numeric value representing an x axis value in the spectra
#' of `x`. This is \eqn{x2} in the equation displayed above. If `x2 = NULL`,
#' the default humification indices will be computed (see details section).
#'
#' @return An object of class [`ir`][ir::ir_new_ir] with additional
#' columns for additional humification indices. If `x1 = NULL` and `x2 = NULL`,
#' these are four new columns (hi1, hi2, hi3 and hi4) that correspond to the humification
#' indices defined in the details section. If `x1` and `x2` are not `NULL`,
#' the columns have names `"hi_x1_x2"` where `x1` and `x2` are replaced by the
#' respective values.
#'
#' @examples
#' # get sample data
#' library(ir)
#'
#' # compute default humification indices
#' d <- ir_sample_data[1:5, ]
#' d <- irp_hi(d)
#'
#' # compute custom humification index
#' d <- ir_sample_data
#' d <- irp_hi(d, x1 = 2900, x2 = 1090)
#'
#' # compute custom humification indices
#' d <- ir_sample_data
#' d <- irp_hi(d, x1 = c(2900, 1630), x2 = 1090)
#'
#' @references
#'   \insertAllCited{}
#'
#' @export
irp_hi <- function(x, x1 = NULL, x2 = NULL) {

  default <- TRUE
  stopifnot(inherits(x, "ir"))
  if((is.null(x1) && !is.null(x2)) || (!is.null(x1) && is.null(x2))) {
    rlang::abort("Only one of `x1` and `x2` is `NULL`. Both of `x1` and `x2` must be `NULL` if you want to compute the default humification indices. If you want to compute custom humification indices, both `x1` and `x2` must be numeric vectors.")
  }
  if(!is.null(x1) && !is.null(x2)) {
    if(!is.numeric(x1)) {
      rlang::abort("`x1` must be numeric or `NULL`.")
    }
    if(!is.numeric(x2)) {
      rlang::abort("`x2` must be numeric or `NULL`.")
    }
    cond <- length(x2)
    if(cond != 1) {
      rlang::abort(paste0("x2 must be a numeric value, but is of length ", cond, "."))
    }
    default <- FALSE
  } else {
   x1 <- c(1420, 1510, 1630, 1720)
   x2 <- 1090
  }

  hi_wn <- data.frame(numerator = x1,
                      denominator = rep(x2, length(x1)),
                      stringsAsFactors = FALSE)
  x_flat <- ir::ir_flatten(x)
  hi_index <- data.frame(
    numerator = ir::ir_get_wavenumberindex(x_flat, wavenumber = hi_wn$numerator, warn = TRUE),
    denominator = ir::ir_get_wavenumberindex(x_flat, wavenumber = hi_wn$denominator, warn = TRUE),
    stringsAsFactors = FALSE)
  hi <- as.data.frame(t(x_flat[hi_index$numerator, -1]/x_flat[hi_index$denominator, -1]))
  colnames(hi) <- paste0("hi_", x1, "_", x2)
  if(default) {
    colnames(hi) <- paste0("hi", seq_len(ncol(hi)))
  }
  dplyr::bind_cols(x, hi)

}
