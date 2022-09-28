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

#' Predicts PLSR scores for new data (data format as in irpetpaper)
#'
#' @keywords internal
#' @noRd
irp_make_predictions_plsr <- function(x, m_pls, config) {

  # format data
  x <-
    tibble::tibble(
      x =
        x %>%
        ir::ir_flatten() %>%
        dplyr::select(-1) %>%
        t()  %>%
        tibble::as_tibble() %>%
        setNames(nm = paste0("V", x$spectra[[1]]$x %>% as.character())) %>% # ---note: this works only because all spectra were clipped to the same range
        as.matrix()
    )

  # get plsr scores
  res <-
    x %>%
    dplyr::mutate(
      x = {
        # get original scores to compute standard deviation of first compound (see below)
        res_or <-
          stats::predict(
            m_pls,
            ncomp = seq_len(config$pls$ncomp),
            type = "scores"
          ) %>%
          tibble::as_tibble()

        stats::predict(
          m_pls,
          ncomp = seq_len(config$pls$ncomp),
          type = "scores",
          newdata = x
        ) %>%
          tibble::as_tibble() %>%
          setNames(nm = stringr::str_remove_all(colnames(.), pattern = " ")) %>%
          purrr::map_dfc( function(.x) {
            #---note: scale the extracted components by dividing them by the standard deviation of the component with the largest variance. This is done to keep priors for slope coefficients roughly on the same scale (see Piironen.2020).
            .x/sd(res_or[, 1, drop = TRUE])
          }) %>%
          as.matrix()
      }
    )
}

#' Computes predictions using draws for model coefficients for a Beta regression model with constant dispersion parameter
#'
#' @keywords internal
#' @noRd
irp_mcmc_predictions_beta_logit <- function(x, draws, config) {

  # linear predictor
  mu <- draws$Intercept + as.matrix(draws %>% dplyr::select(dplyr::starts_with("b["))) %*% t(x$x)
  mu <- as.data.frame(config$likelihood$linkinv(mu), stringsAsFactors = FALSE)

  # predictions
  yhat <-
    purrr::map_dfc(mu, function(.x) {
      stats::rbeta(n = length(.x), shape1 = .x * draws$phi, shape2 = (1 - .x) * draws$phi)
    })

  # scale
  yhat * config$data_scale$y_scale + config$data_scale$y_center

}
