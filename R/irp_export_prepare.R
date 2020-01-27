#' Stores errors and units of error or units columns of objects of class \code{ir} in separate columns.
#'
#' \code{irp_export_prepare} stores errors (\code{\link[errors:errors]{errors}}) and
#' units (\code{\link[units:set_units]{units}}) of error or
#' units columns in an object of class \code{\link[ir:ir_new_ir]{ir}} in separate columns.
#'
#' @param x An object of class \code{\link[ir:ir_new_ir]{ir}}.
#' @return An object of class \code{\link[ir:ir_new_ir]{ir}} with a new column for
#' each column that is of class \code{"units"} and each column that is of class
#' \code{"errors"}. New columns are named old_colmn_name_errors and old_colmn_name_units,
#' respectively.
#' @export
irp_export_prepare <- function(x) {

  ir::ir_check_ir(x)

  d <- data.frame(errors = purrr::map_lgl(x, inherits, "errors"),
                  units = purrr::map_lgl(x, inherits, "units"),
                  stringsAsFactors = FALSE)

  d_errors <- purrr::map2(d$errors, seq_along(d$errors), function(y, z) {
    if(y) {
      d <- data.frame(x = quantities::quantities(x[, z, drop = TRUE])$errors,
                 stringsAsFactors = FALSE)
      colnames(d) <- paste0(colnames(x)[[z]], "_errors")
      d
    }
  })

  d_units <- purrr::map2(d$units, seq_along(d$units), function(y, z) {
    if(y) {
      d <- data.frame(x = rep(as.character(quantities::quantities(x[, z, drop = TRUE])$units), nrow(x)),
                      stringsAsFactors = FALSE)
      colnames(d) <- paste0(colnames(x)[[z]], "_units")
      d
    }
  })

  dplyr::bind_cols(x, d_errors, d_units)

}
