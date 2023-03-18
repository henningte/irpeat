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


#' `check_irpeatmodels` combined with a check for the 'dimreduce' package
#'
#' @inheritParams irp_function_factory_eb1079
#'
#' @keywords internal
#' @noRd
check_irpeatmodels_and_dimreduce <- function(irpeatmodels_required_version) {

  check_irpeatmodels(version = irpeatmodels_required_version)
  if(! requireNamespace("dimreduce", quietly = TRUE)) {
    rlang::abort("You have to install the 'dimreduce' package to use this function. The 'dimreduce' package is available from 'https://github.com/jpiironen/dimreduce'.")
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

    # preprocessing
    x <- irp_preprocess_for(x = x, variable = target_variable_name)

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
            dplyr::pull(.data$is_in_prediction_domain)
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

    # additional transformations
    res <-
      switch(
        target_variable_name,
        "d13C_1" = purrr::map_dfc(res, function(.x) irp_atompercent_to_delta(.x * 100, r = 0.0112372)), # standard atom percent for VPDB
        "d15N_1" = purrr::map_dfc(res, function(.x) irp_atompercent_to_delta(.x * 100, r = 0.0036765)), # standard atom percent for AIR,
        res
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
      stats::setNames(nm = target_variable_name)

    cbind(x_or, res, res_pd %>% stats::setNames(nm = paste0(target_variable_name, "_in_pd")))

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
        stats::setNames(nm = paste0("V", x$spectra[[1]]$x %>% as.character())) %>% # ---note: this works only because all spectra were clipped to the same range
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
          stats::setNames(nm = stringr::str_remove_all(colnames(.), pattern = " ")) %>%
          purrr::map_dfc( function(.x) {
            #---note: scale the extracted components by dividing them by the standard deviation of the component with the largest variance. This is done to keep priors for slope coefficients roughly on the same scale (see Piironen.2020).
            .x/stats::sd(res_or[, 1, drop = TRUE])
          }) %>%
          as.matrix()
      }
    )
}


#' Predicts ISPCA scores for new data (data format as in irpetpaper)
#'
#' @keywords internal
#' @noRd
irp_make_predictions_ispca <- function(x, m_pls, config) {

  #---todo (code copied from irp_make_predictions_plsr)

  # format data
  x <-
    tibble::tibble(
      x =
        x %>%
        ir::ir_flatten() %>%
        dplyr::select(-1) %>%
        t()  %>%
        tibble::as_tibble(.name_repair = "minimal") %>%
        stats::setNames(nm = paste0("V", x$spectra[[1]]$x %>% as.character())) %>% # ---note: this works only because all spectra were clipped to the same range
        as.matrix()
    )

  # get ispca scores
  res <-
    x %>%
    dplyr::mutate(
      x = {
        # get original scores to compute standard deviation of first compound (see below)
        # get original scores to compute standard deviation of first compound (see below)
        res_or <-
          m_pls$z %>%
          tibble::as_tibble(.name_repair = "minimal")

        stats::predict(
          m_pls,
          xnew = x
        ) %>%
          tibble::as_tibble() %>%
          purrr::map_dfc( function(.x) {
            #---note: scale the extracted components by dividing them by the standard deviation of the component with the largest variance. This is done to keep priors for slope coefficients roughly on the same scale (see Piironen.2020).
            .x/stats::sd(res_or[, 1, drop = TRUE])
          }) %>%
          dplyr::select(seq_len(config$pls$ncomp)) %>%
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
      if(all(is.na(.x))) {
        rep(NA_real_, length(.x))
      } else {
        stats::rbeta(n = length(.x), shape1 = .x * phi, shape2 = (1 - .x) * phi)
      }
    })

  # scale
  yhat * config$data_scale$y_scale + config$data_scale$y_center

}


#' Computes predictions using draws for model coefficients for a Normal regression model with constant standard deviation
#'
#' @keywords internal
#' @noRd
irp_mcmc_predictions_normal_identity_non_centered <- function(x, draws, config) {

  # scale intercept and slopes
  intercept <- draws$Intercept * config$parameter_scale$Intercept_sigma + config$parameter_scale$Intercept_mu
  b <- draws %>%
    dplyr::select(dplyr::starts_with("b[")) %>%
    magrittr::multiply_by(config$parameter_scale$b_sigma) %>%
    magrittr::add(config$parameter_scale$b_mu)

  # linear predictor
  mu <- intercept + as.matrix(b) %*% t(x$x)
  mu <- as.data.frame(config$likelihood$linkinv(mu), stringsAsFactors = FALSE)

  # predictions
  yhat <-
    purrr::map_dfc(mu, function(.x) {
      if(all(is.na(.x))) {
        rep(NA_real_, length(.x))
      } else {
        stats::rnorm(n = length(.x), mean = .x, sd = draws$sigma/config$parameter_scale$sigma_scale)
      }
    })

  # scale
  yhat * config$data_scale$y_scale + config$data_scale$y_center

}


#### preprocessing helper functions ####

#' Helper function to take care of the assignment of arguments in `irp_preprocess()` from a config list
#'
#' @keywords internal
#' @noRd
irp_preprocess_unpack_config <- function(x, config) {

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
    bc_cutoff = config$irp_preprocess$bc_cutoff,
    bc_do_impute = config$irp_preprocess$bc_do_impute,
    do_smooth = config$irp_preprocess$do_smooth,
    do_normalise = config$irp_preprocess$do_normalise,
    normalise_method = config$irp_preprocess$normalise_method,
    do_bin = config$irp_preprocess$do_bin,
    bin_width = config$irp_preprocess$bin_width,
    bin_new_x_type = config$irp_preprocess$bin_new_x_type,
    do_scale = config$irp_preprocess$do_scale,
    scale_center = config$data_scale$x_center,
    scale_scale = config$data_scale$x_scale,
    do_return_as_ir = config$irp_preprocess$do_return_as_ir
  )

}

#' Helper function to preprocess spectra for prediction with the models from project eb1079
#'
#' @keywords internal
#' @noRd
irp_preprocess_eb1079 <- function(x, config) {

  ## custom preprocessing special to the models for eb1079

  # First: get a baseline which can be subtracted from all spectra even with negative CO2 peaks. I have to do a SG smoothing and regional interpolation here to avoid negative CO2 peaks elsewhere to corrupt the baseline. I also have to interpolate linearly the CO2 peak around 670 cm$^{-1}$ since this peak corrupts the baseline and does not get smoothed out completely by the SG smoothing.
  x_bl <-
    x %>%
    ir::ir_interpolate(start = NULL, dw = 1) %>%
    ir::ir_smooth(method = "sg", n = 91) %>% #---note: new
    ir::ir_clip(range = config$irp_preprocess$clip_range) %>%
    ir::ir_interpolate_region(range = tibble::tibble(start = c(650, 2230), end = c(695, 2410))) %>%
    ir::ir_bc(method = "rubberband", return_bl = TRUE, do_impute = TRUE)

  # clipping and baseline correction
  x <-
    x %>%
    ir::ir_interpolate(start = NULL, dw = 1) %>%
    ir::ir_clip(range = config$irp_preprocess$clip_range) %>%
    ir::ir_subtract(x_bl) %>%
    ir::ir_bc(method = "rubberband", do_impute = TRUE)

  # preprocess the spectra
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

}


#' Helper function to preprocess spectra for prediction with the models from project eb1014 (eac_1, edc_1)
#'
#' @keywords internal
#' @noRd
irp_preprocess_eb1014 <- function(x, config) {

  # preprocess the spectra
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
    bc_cutoff = config$irp_preprocess$bc_cutoff,
    bc_do_impute = config$irp_preprocess$bc_do_impute,
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

}
