#' Converts atom percentages to corresponding delta values
#'
#' @param x A numeric vector with values representing the atom percentage of a
#' specific isotope.
#'
#' @param r A numeric value representing the delta value for the corresponding
#' isotope standard.
#'
#' @return `x` with values converted to delta values.
#'
#' @keywords Internal
#' @noRd
#'
#' @examples
#' # convert a 13C atom percentage to delta values
#' irp_atompercent_to_delta(x = 1.2, r = 0.0112372)
#'
#' @export
irp_atompercent_to_delta <- function(x, r) {

  stopifnot(is.numeric(x))
  stopifnot(is.numeric(r) && length(r) == 1L)

  -1000*(1-x/((100-x)*r))
}
