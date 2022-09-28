#' Checks whether the 'irpeatmodels' package is loaded
#'
#' @keywords internal
#' @noRd
check_irpeatmodels <- function(version = "0.1.0") {
  if(! requireNamespace("irpeatmodels", versionCheck = list(op = ">=", version = version), quietly = TRUE)) {
    rlang::abort(paste0("You have to install the 'irpeatmodels' package (>=", version,") to use this function."))
  }
}

#' Summarizes predictions of models where predictions are given as MCMC draws from the posterior predictive distribution
#'
#' @param x A data frame with the predicted values (rows represent draws and
#' columns samples).
#'
#' @param do_summary A logical value indicating whether values in `x` should be
#' summarized (`TRUE`) or not (`FALSE`). If `TRUE`, `summary_funtion_mean` and
#' `summary_function_sd` are used to summarize the values in `x` for each sample,
#' and a [`quantities::quantities`] object is returned. If `FALSE`, the values in
#' `x` are converted to [`units::units`] objects, and `x` is returned as list
#' (with elements representing the draws for each sample).
#'
#' @param summary_function_mean A function used to summarize the values in `x`
#' (average).
#'
#' @param summary_function_sd A function used to summarize the values in `x`
#' (spread).
#'
#' @param x_unit A character value giving the unit which will be assigned to
#' values in `x`.
#'
#' @keywords internal
#' @noRd
irp_summarize_predictions <- function(x, x_unit, do_summary, summary_function_mean = mean, summary_function_sd = stats::sd) {

  stopifnot(is.data.frame(x))
  stopifnot(is.character(x_unit) && length(x_unit) == 1L)
  stopifnot(is.logical(do_summary) && length(do_summary) == 1L)

  if(do_summary) {
    purrr::map_dbl(x, summary_function_mean) %>%
      quantities::set_quantities(
        unit = x_unit,
        errors = purrr::map_dbl(x, summary_function_sd),
        mode = "standard"
      )
  } else {
    purrr::map(x, units::set_units(x_unit, mode = "standard"))
  }

}

