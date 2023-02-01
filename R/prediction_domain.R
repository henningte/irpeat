#' Creates an object of class `irp_prediction_domain`
#'
#' @name irp_prediction_domain
#'
#' @param x A data frame with a row for each wavenumber value in the training
#' data and the following columns:
#' \describe{
#'   \item{`x`}{A numeric value representing the wavenumber value
#'   \[cm$^{-1}$\].}
#'   \item{`ymin`}{A numeric value representing the minimum predictor variable
#'   value in the training data at this wavenumber value.}
#'   \item{`ymax`}{A numeric value representing the maximum predictor variable
#'   value in the training data at this wavenumber value.}
#' }
#' `x` may contain additional columns.
#'
#' @return An object of class `irp_prediction_domain`. This is the same as `x`,
#' but with an additional class label.
#'
#' @export
new_irp_prediction_domain <- function(x) {

  stopifnot(is.data.frame(x))
  stopifnot(c("x", "ymin", "ymax") %in% colnames(x))
  stopifnot(is.numeric(x$x))
  stopifnot(is.numeric(x$ymin))
  stopifnot(is.numeric(x$ymax))
  stopifnot(all(x$x > 0))
  stopifnot(all(purrr::map2_lgl(x$ymax, x$ymin, magrittr::is_greater_than) + purrr::map2_lgl(x$ymax, x$ymin, magrittr::equals)))

  structure(x, class = c("irp_prediction_domain", class(x)))

}

#' Converts an object to class `irp_prediction_domain`
#'
#' @name irp_prediction_domain-conversion
#'
#' @param x An object.
#'
#' @param ... Further arguments passed to methods.
#'
#' @return An object of class `irp_prediction_domain`.
#'
#' @export
irp_as_irp_prediction_domain <- function(x, ...) {
  UseMethod("irp_as_irp_prediction_domain")
}

#' @rdname irp_prediction_domain-conversion
#' @export
irp_as_irp_prediction_domain.ir_flat <- function(x, ...) {

  x_wavenumbers <- x$x

  .x <-
    x %>%
    dplyr::select(-x) %>%
    t() %>%
    tibble::as_tibble(.name_repair = "minimal")

  tibble::tibble(
    x = x_wavenumbers,
    ymin = purrr::map_dbl(.x, min, na.rm = TRUE),
    ymax = purrr::map_dbl(.x, max, na.rm = TRUE)
  ) %>%
    new_irp_prediction_domain()

}

#' @rdname irp_prediction_domain-conversion
#' @export
irp_as_irp_prediction_domain.ir <- function(x, ...) {
  irp_as_irp_prediction_domain(ir::ir_flatten(x))
}

#' @rdname irp_prediction_domain-conversion
#' @export
irp_as_irp_prediction_domain.data.frame <- function(x, ...) {
  new_irp_prediction_domain(x)
}

#' Plotting method for `irp_prediction_domain` objects
#'
#' @param x An object of class `irp_prediction_domain`.
#' @param ... Further parameters; passed to `ggplot2::geom_ribbon`.
#'
#' @return A `ggplot` object.
#'
#' @export
plot.irp_prediction_domain <- function(x, ...) {

  ggplot2::ggplot(x, ggplot2::aes(x = .data$x)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax), ...)

}

#' Checks whether spectra are within a prediction domain
#'
#' @param x An object of class `ir`. The wavenumber values in the spectra of `x`
#' need to be identical.
#'
#' @param prediction_domain An object of class `irp_prediction_domain`. The
#' wavenumber values in `prediction_domain` need to be identical to those in `x`.
#'
#' @return `x` with an additional column `is_in_prediction_domain` with value
#' `TRUE` if a spectrum in `x` has intensity values which are within the
#' prediction domain and `FALSE` if not.
#'
#' @export
irp_is_in_prediction_domain <- function(x, prediction_domain) {

  stopifnot(inherits(x, "ir"))
  stopifnot(inherits(prediction_domain, "irp_prediction_domain"))

  # check matching wavenumbers
  index_nonempty_spectra <- purrr::map_lgl(x$spectra, function(.x) ! (nrow(.x) == 0 || all(is.na(.x$y))))
  if(! all(purrr::map_lgl(x$spectra[index_nonempty_spectra], function(.x) identical(.x$x, x$spectra[index_nonempty_spectra][[1]]$x)))) {
    rlang::abort("Not all non-empty spectra in `x` have the same x axis values. Make sure that all spectra have the same x axis values (except for empty spectra).")
  }
  stopifnot(identical(prediction_domain$x, x$spectra[index_nonempty_spectra][[1]]$x))

  x_flat <-
    x %>%
    ir::ir_flatten() %>%
    dplyr::select(-x)

  x_check_lower <-
    purrr::map_lgl(x_flat, function(.x) {
      all(prediction_domain$ymin <= .x)
    })

  x_check_upper <-
    purrr::map_lgl(x_flat, function(.x) {
      all(prediction_domain$ymax >= .x)
    })

  x %>%
    dplyr::mutate(
      is_in_prediction_domain = (x_check_lower + x_check_upper) == 2,
      is_in_prediction_domain = unname(.data$is_in_prediction_domain)
    )

}
