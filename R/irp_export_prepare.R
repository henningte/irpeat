#' Stores errors and units of error or units columns of objects of class `ir` in separate columns
#'
#' `irp_export_prepare` stores errors ([errors::errors]) and
#' units ([`units`][units::set_units]) of error or
#' units columns in an object of class [`ir`][ir::ir_new_ir] in separate columns.
#'
#' @param x An object of class [`ir`][ir::ir_new_ir].
#'
#' @return An object of class [`ir`][ir::ir_new_ir] with a new column for
#' each column that is of class `"units"` and each column that is of class
#' `"errors"`. New columns are named old_column_name_errors and old_column_name_units,
#' respectively.
#'
#' @examples
#' library(ir)
#'
#' irp_export_prepare(ir::ir_sample_data[1:5, ])
#'
#' @export
irp_export_prepare <- function(x) {

  stopifnot(inherits(x, "ir"))

  d <- data.frame(errors = purrr::map_lgl(x, inherits, "errors"),
                  units = purrr::map_lgl(x, inherits, "units"),
                  quantities = purrr::map_lgl(x, inherits, "quantities"),
                  stringsAsFactors = FALSE)

  d_errors <- purrr::map2(d$errors, seq_along(d$errors), function(y, z) {
    if(y) {
      d <- data.frame(x = quantities::quantities(x[, z, drop = TRUE])$errors,
                 stringsAsFactors = FALSE)
      colnames(d) <- paste0(colnames(x)[[z]], "_errors")
      d
    }
  })

  d_units_units <- purrr::map2(d$units  & !d$quantities, seq_along(d$quantities), function(y, z) {
    if(y) {
      d <- data.frame(x = rep(as.character(attr(x[, z, drop = TRUE], "units")), nrow(x)),
                      stringsAsFactors = FALSE)
      colnames(d) <- paste0(colnames(x)[[z]], "_units")
      d
    }
  })

  d_quantities_units <- purrr::map2(d$quantities, seq_along(d$quantities), function(y, z) {
    if(y) {
      d <- data.frame(x = rep(as.character(quantities::quantities(x[, z, drop = TRUE])$units), nrow(x)),
                      stringsAsFactors = FALSE)
      colnames(d) <- paste0(colnames(x)[[z]], "_units")
      d
    }
  })

  dplyr::bind_cols(x, d_errors, d_units_units, d_quantities_units)

}
