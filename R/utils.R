#' Checks whether the 'irpeatmodels' package is loaded
#'
#' @keywords internal
#' @noRd
check_irpeatmodels <- function(version = "0.1.0") {
  if(! requireNamespace("irpeatmodels", versionCheck = list(op = ">=", version = version), quietly = TRUE)) {
    rlang::abort(paste0("You have to install the 'irpeatmodels' package (>=", version,") to use this function."))
  }
}

#' `check_irpeatmodels` combined with a check for the 'pls' package
#'
#' @inheritParams irp_function_factory_eb1079
#'
#' @keywords internal
#' @noRd
check_irpeatmodels_and_pls <- function(irpeatmodels_required_version) {

  check_irpeatmodels(version = irpeatmodels_required_version)
  if(! requireNamespace("pls", quietly = TRUE)) {
    rlang::abort("You have to install the 'pls' package to use this function.")
  }


}



#' Function factory to generate prediction functions for models from project eb1079
#'
#' @param m A data frame with the model parameters (each column is a
#' parameter and each row a draw from the posterior distribution).
#'
#' @param m_pls A `pls::mvr` or `dimreduce::dimreduce` model.
#'
#' @param config A list with configuration parameters for the prediction.
#'
#' @param prediction_domain A list with two elements:
#' \describe{
#'   \item{`train`}{An
#'   [`irp_prediction_domain`](irpeat::new_irp_prediction_domain) object
#'   representing the prediction domain for the training data.}
#'   \item{`test`}{An
#'   [`irp_prediction_domain`](irpeat::new_irp_prediction_domain) object
#'   representing the prediction domain for the testing data.}
#' }
#'
#' @param target_variable_name A character value representing the column name
#' for the column with predicted values.
#'
#' @param x_unit A character value giving the unit which will be assigned to
#' predicted values.
#'
#' @param irpeatmodels_required_version A character value with format "x.y.z"
#' representing the minimum version of the 'irpeatmodels' package required
#' to make predictions with the function which is generated.
#'
#' @param .f_check_packages A function which checks whether the correct version
#' of 'irpeatmodels' is installed and of other packages which may be required.
#' Must take the argument `irpeatmodels_required_version` as input.
#'
#' @param .f_dimreduce A function which returns a data
#' frame with score values of the dimension reduction method predicted for new
#' data. See `irp_make_predictions_plsr` for the structure and value of the
#' function.
#'
#' @param .f_predict A function which returns a data frame with predicted values
#' (samples in columns, rows represent draws from the posterior predictive
#' distribution). See `irp_mcmc_predictions_beta_logit` for the structure and
#' value of the function.
#'
#' @return A function which allows to make predictions with a model computed in
#' project eb1079.
#'
#' @keywords internal
#' @noRd
irp_function_factory_eb1079 <- function(m, m_pls, config, prediction_domain, target_variable_name, x_unit, irpeatmodels_required_version, .f_check_packages, .f_dimreduce, .f_predict) {

  stopifnot(is.character(target_variable_name) && length(target_variable_name) == 1)
  stopifnot(is.character(x_unit) && length(x_unit) == 1)

  function(x, do_summary = FALSE, summary_function_mean = mean, summary_function_sd = stats::sd, check_prediction_domain = "train") {

    .f_check_packages(irpeatmodels_required_version = irpeatmodels_required_version)
    stopifnot(inherits(x, "ir"))
    stopifnot(is.logical(do_summary) && length(do_summary) == 1)

    x_or <- x

    # check spectra
    x_flat <- ir::ir_flatten(x)
    if(x_flat$x[[1]] > config$irp_preprocess$clip_range$start) {
      rlang::warn(paste0("The minimum wavenumber value in `x` is ", x_flat$x[[1]], " , but should be ", config$irp_preprocess$clip_range$start, " or smaller."))
    }
    if(x_flat$x[[nrow(x_flat)]] < config$irp_preprocess$clip_range$end) {
      rlang::warn(paste0("The maximum wavenumber value in `x` is ", x_flat$x[[nrow(x)]], " , but should be ", config$irp_preprocess$clip_range$end, " or larger."))
    }

    # preprocess the spectra
    x <-
      irp_preprocess(
        x,
        do_interpolate = config$irp_preprocess$do_interpolate,
        interpolate_start = config$irp_preprocess$interpolate_start,
        interpolate_dw = config$irp_preprocess$interpolate_dw,
        do_clip = config$irp_preprocess$do_clip,
        clip_range = config$irp_preprocess$clip_range,
        do_interpolate_region = config$irp_preprocess$do_interpolate_region,
        interpolate_region_range = config$irp_preprocess$interpolate_region_range,
        do_bc = config$irp_preprocess$do_bc,
        bc_method = config$irp_preprocess$bc_method,
        # bc_cutoff = config$irp_preprocess$bc_cutoff,
        bc_cutoff = 0,
        # bc_do_impute = config$irp_preprocess$bc_do_impute,
        bc_do_impute = TRUE,
        do_smooth = config$irp_preprocess$do_smooth,
        do_normalise = config$irp_preprocess$do_normalise,
        normalise_method = config$irp_preprocess$normalise_method,
        do_bin = config$irp_preprocess$do_bin,
        bin_width = config$irp_preprocess$bin_width,
        bin_new_x_type = config$irp_preprocess$bin_new_x_type,
        do_scale = config$irp_preprocess$do_scale,
        scale_center = config$data_scale$x_center,
        scale_scale = config$data_scale$x_scale,
        do_return_as_ir = TRUE
      )

    # check prediction domain
    prediction_domain <-
      switch(
        check_prediction_domain,
        "train" = prediction_domain$train,
        "test" = prediction_domain$test,
        "none" = NULL
      )

    res_pd <-
      if(check_prediction_domain != "none") {
        tibble::tibble(
          y =
            x %>%
            irp_is_in_prediction_domain(prediction_domain = prediction_domain) %>%
            dplyr::pull(is_in_prediction_domain)
        )
    } else {
      tibble::tibble(
        y = rep(NA, nrow(x_or))
      )
    }

    # get plsr scores
    res <-
      .f_dimreduce(
        x = x,
        m_pls = m_pls,
        config = config
      )

    # predict
    res <-
      .f_predict(
        x = res,
        draws = m,
        config = config
      )

    # summarize and add unit
    res <-
      irp_summarize_predictions(
        x = res,
        x_unit = x_unit,
        do_summary = do_summary,
        summary_function_mean = mean,
        summary_function_sd = stats::sd
      )

    res <-
      tibble::tibble(y = res) %>%
      setNames(nm = target_variable_name)

    cbind(x_or, res, res_pd %>% setNames(nm = paste0(target_variable_name, "_in_pd")))

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
    purrr::map(x, units::set_units, value = x_unit, mode = "standard")
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
        tibble::as_tibble(.name_repair = "minimal") %>%
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

  phi <-  draws$phi/config$parameter_scale$phi_scale

  # predictions
  yhat <-
    purrr::map_dfc(mu, function(.x) {
      stats::rbeta(n = length(.x), shape1 = .x * phi, shape2 = (1 - .x) * phi)
    })

  # scale
  yhat * config$data_scale$y_scale + config$data_scale$y_center

}


#' Computes predictions using draws for model coefficients for a Normal regression model with constant standard deviation
#'
#' @keywords internal
#' @noRd
irp_mcmc_predictions_normal_identity <- function(x, draws, config) {

  # linear predictor
  mu <- draws$Intercept + as.matrix(draws %>% dplyr::select(dplyr::starts_with("b["))) %*% t(x$x)
  mu <- as.data.frame(config$likelihood$linkinv(mu), stringsAsFactors = FALSE)

  # predictions
  yhat <-
    purrr::map_dfc(mu, function(.x) {
      stats::rnorm(n = length(.x), mean = .x, sd = draws$sigma/config$parameter_scale$sigma_scale)
    })

  # scale
  yhat * config$data_scale$y_scale + config$data_scale$y_center

}
