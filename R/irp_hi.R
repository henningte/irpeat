#' Humification indices from mid infrared spectra.
#'
#' \code{irp_hi} computes various humification indices as reported by
#' \insertCite{Broder.2012;textual}{irpeat} for mid infrared
#' spectra with the x values representing wavenumber values (no checks are performed):
#' The following humification indices are computed (values
#' represent wavenumbers [cm\eqn{^{-1}}]) \insertCite{Broder.2012}{irpeat}:
#' \describe{
#'   \item{1420/1090}{OH and CO of phenols or CH of CH{_1} and CH{_3}
#'   groups (phenolic and aliphatic structures)/polysaccharides}
#'   \item{1510/1090}{Aromatic C=C or C=O of amides/polysaccharides}
#'   \item{1630/1090}{Aromatic C=C and COO\eqn{^{-1}} (aromatics and aromatic
#'   or aliphatic carboxylates)/polysaccharides}
#'   \item{1720/1090}{carbonylic and carboxylic C=O (carboxylic acids and
#'   aromatic esters)/polysaccharides}
#' }
#'
#' @param x An object of class \code{\link[ir:ir_new_ir]{ir}}.
#' @return An object of class \code{\link[ir:ir_new_ir]{ir}} with three new
#' columns (hi1, hi2, hi3 and hi4) that correspond to the humification indices
#' defined in the description.
#' @references
#'   \insertAllCited{}
#' @export
irp_hi <- function(x) {

  ir::ir_check_ir(x)
  x_flat <- ir::ir_flatten(x)
  hi_wn <- data.frame(numerator = c(1420, 1510, 1630, 1720),
                   denominator = c(1090, 1090, 1090, 1090),
                   stringsAsFactors = FALSE)
  hi_index <- data.frame(
    numerator = ir::ir_get_wavenumberindex(x_flat, wavenumber = hi_wn$numerator, warn = TRUE),
    denominator = ir::ir_get_wavenumberindex(x_flat, wavenumber = hi_wn$denominator, warn = TRUE),
    stringsAsFactors = FALSE)
  hi <- as.data.frame(t(x_flat[hi_index$numerator, -1]/x_flat[hi_index$denominator, -1]))
  colnames(hi) <- paste0("hi", seq_len(ncol(hi)))
  dplyr::bind_cols(x, hi)

}
