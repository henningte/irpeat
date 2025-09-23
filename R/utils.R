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
irp_summarize_predictions <- function(x, x_unit, do_summary, return_as_list, summary_function_mean = mean, summary_function_sd = stats::sd) {

  stopifnot(is.data.frame(x))
  stopifnot(is.character(x_unit) && length(x_unit) == 1L)
  stopifnot(is.logical(do_summary) && length(do_summary) == 1L)

  if(do_summary) {
    purrr::map_dbl(x, summary_function_mean) %>%
      quantities::set_quantities(
        unit = x_unit,
        errors = purrr::map_dbl(x, summary_function_sd),
        mode = "standard"
      ) %>%
      unname()
  } else if (return_as_list) {
    purrr::map(x, units::set_units, value = x_unit, mode = "standard") %>%
      unname()
  } else if (! return_as_list) {
    posterior::rvar(as.matrix(x)) |>
      irp_set_units_rvar(value = x_unit, mode = "standard")
  }

}

#' Predicts PLSR scores for new data (data format as in irpeatpaper)
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



#### eb1079 ####


#' Function factory for prediction functions from project eb1149 (degree of decomposition)
#'
#' @param model A [`brmsfit`](brms::brm) object.
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
#' @param irpeatmodels_required_version A character value with format "x.y.z"
#' representing the minimum version of the 'irpeatmodels' package required
#' to make predictions with the function which is generated.
#'
#' @param .f_check_packages A function which checks whether the correct version
#' of 'irpeatmodels' is installed and of other packages which may be required.
#' Must take the argument `irpeatmodels_required_version` as input.
#'
#' @return A function that makes predictions with a model computed in
#' project eb1079.
#'
#' @keywords Internal
#' @noRd
irp_function_factory_eb1079 <- function(target_variable, model, config, prediction_domain, irpeatmodels_required_version = "0.0.0") {

  .f_check_packages <- function() {

    check_irpeatmodels(version = irpeatmodels_required_version)
    rlang::is_installed("brms")
    if(! requireNamespace("posterior", versionCheck = list(op = ">=", version = "1.5.0"), quietly = TRUE)) {
      rlang::abort(paste0("You have to install the 'posterior' package (>=", "1.5.0",") to use this function."))
    }

  }

  function(x, do_summary = FALSE, summary_function_mean = mean, summary_function_sd = stats::sd, check_prediction_domain = "train", return_as_list = FALSE) {

    .f_check_packages()
    stopifnot(inherits(x, "ir"))
    stopifnot(is.logical(do_summary) && length(do_summary) == 1)
    stopifnot(is.logical(return_as_list) && length(return_as_list) == 1)
    if(do_summary && return_as_list) {
      stop("Both `do_summary` and `return_as_list` are set to `TRUE`, but only one of both must be `TRUE`.")
    }


    x_or <- x
    x <- irp_preprocess_eb1079(x = x, config = config)

    # check prediction domain
    prediction_domain <-
      switch(
        check_prediction_domain,
        "train" = prediction_domain$train,
        "test" = prediction_domain$test,
        "none" = NULL
      )

    x_in_pd <-
      if(check_prediction_domain != "none") {
        tibble::tibble(
          is_in_prediction_domain =
            x |>
            irp_is_in_prediction_domain(prediction_domain = prediction_domain) |>
            dplyr::pull(.data$is_in_prediction_domain)
        )
      } else {
        tibble::tibble(
          is_in_prediction_domain = rep(NA, nrow(x_or))
        )
      }

    newdata <-
      tibble::tibble(
        x =
          x |>
          ir::ir_flatten() |>
          dplyr::select(-1) |>
          t() |>
          tibble::as_tibble(.name_repair = "minimal") |>
          stats::setNames(nm = paste0("V", seq_along(x$spectra[[1]]$x))) |>
          as.matrix(),
        y_err =
          if(target_variable %in% c("dgf0_1")) {
            0.00001
          } else {
            NULL
          }
      )


    yhat <-
      brms::posterior_predict(object = model, newdata = newdata) |>
      as.data.frame()

    yhat <-
      irp_summarize_predictions(
        x = yhat * config$model_scale$y_scale + config$model_scale$y_center,
        x_unit = config$unit,
        do_summary = do_summary,
        return_as_list = return_as_list,
        summary_function_mean = summary_function_mean,
        summary_function_sd = summary_function_sd
      )

    x_or[[target_variable]] <- yhat
    x_or[[paste0(target_variable, "_in_pd")]] <- x_in_pd$is_in_prediction_domain

    x_or

  }

}


#### eb1149 ####

#' Function factory for prediction functions from project eb1149 (degree of decomposition)
#'
#' @param model A [`brmsfit`](brms::brm) object.
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
#' @param irpeatmodels_required_version A character value with format "x.y.z"
#' representing the minimum version of the 'irpeatmodels' package required
#' to make predictions with the function which is generated.
#'
#' @param .f_check_packages A function which checks whether the correct version
#' of 'irpeatmodels' is installed and of other packages which may be required.
#' Must take the argument `irpeatmodels_required_version` as input.
#'
#' @return A function that makes predictions with a model computed in
#' project eb1079.
#'
#' @keywords Internal
#' @noRd
irp_function_factory_eb1149 <- function(model, config, prediction_domain, target_variable_name, irpeatmodels_required_version = "0.0.0") {

  .f_check_packages <- function() {

    check_irpeatmodels(version = irpeatmodels_required_version)
    rlang::is_installed("brms")
    rlang::is_installed("posterior", version = "1.5.0")

  }

  function(x, do_summary = FALSE, summary_function_mean = mean, summary_function_sd = stats::sd, check_prediction_domain = "train", return_as_list = FALSE) {

    .f_check_packages()
    stopifnot(inherits(x, "ir"))
    stopifnot(is.logical(do_summary) && length(do_summary) == 1)
    stopifnot(is.logical(return_as_list) && length(return_as_list) == 1)
    if(do_summary && return_as_list) {
      stop("Both `do_summary` and `return_as_list` are set to `TRUE`, but only one of both must be `TRUE`.")
    }

    x_or <- x
    x <- irp_preprocess_eb1149(x = x, config = config)
    x_in_pd <- irpeat::irp_is_in_prediction_domain(x = x, prediction_domain = prediction_domain)

    newdata <-
      tibble::tibble(
        x =
          x |>
          ir::ir_flatten() |>
          dplyr::select(-1) |>
          t()
      )

    yhat <-
      brms::posterior_predict(object = model, newdata = newdata) |>
      posterior::rvar()

    posterior::draws_of(yhat) <- units::set_units(posterior::draws_of(yhat), value = "g/g", mode = "standard")

    if(return_as_list) {
      yhat <-
        posterior::draws_of(yhat) |>
        as.data.frame() |>
        purrr::map(function(.x) .x)
    } else if(do_summary) {
      yhat <-
        purrr::map(yhat, function(.x) {
          quantities::set_quantities(
            x = summary_function_mean(.x),
            unit = "g/g",
            errors = summary_function_sd(.x),
            mode = "standard"
          )
        })
      yhat <- do.call("c", yhat)
    }

    x_or[[target_variable_name]] <- yhat
    x_or[[paste0(target_variable_name, "_in_pd")]] <- x_in_pd$is_in_prediction_domain

    x_or

  }

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



#' Preprocesses spectra according to a preprocessing config object
#'
#' @param x An ir object to be preprocessed.
#'
#' @param config A list with arguments compatible with `irpeat::irp_preprocess()`.
#'
#' @keywords Internal
#' @noRd
irp_preprocess_eb1079 <- function (x, config) {

  res <- x

  # First: get a baseline which can be subtracted from all spectra even with negative CO2 peaks. I have to do a SG smoothing and regional interpolation here to avoid negative CO2 peaks elsewhere to corrupt the baseline. I also have to interpolate linearly the CO2 peak around 670 cm$^{-1}$ since this peak corrupts the baseline and does not get smoothed out completely by the SG smoothing.
  res_bl <-
    res |>
    ir::ir_interpolate(start = NULL, dw = 1) |>
    ir::ir_smooth(method = "sg", n = 91) |>
    ir::ir_clip(range = config$irp_preprocess$clip_range[[1]]) |>
    ir::ir_interpolate_region(range = tibble::tibble(start = c(645, 2230), end = c(695, 2410))) |>
    ir::ir_bc(method = "rubberband", do_impute = TRUE, return_bl = TRUE)

  ## Second: perform the actual correction

  # clipping and baseline correction
  res <-
    res |>
    ir::ir_interpolate(start = NULL, dw = 1) |>
    ir::ir_clip(range = config$irp_preprocess$clip_range[[1]]) |>
    ir::ir_subtract(res_bl)

  res |>
    #ir::ir_interpolate_region(range = tibble::tibble(start = c(645, 2250), end = c(695, 2400))) |>
    ir::ir_interpolate_region(range = tibble::tibble(start = c(650), end = c(695))) |>
    #ir::ir_bc(method = "rubberband", do_impute = TRUE) |>
    irpeat::irp_preprocess(
      do_interpolate = config$irp_preprocess$do_interpolate,
      interpolate_start = config$irp_preprocess$interpolate_start[[1]],
      interpolate_dw = config$irp_preprocess$interpolate_dw,
      do_clip = config$irp_preprocess$do_clip[[1]],
      clip_range = config$irp_preprocess$clip_range[[1]],
      do_interpolate_region = FALSE,
      interpolate_region_range = config$irp_preprocess$interpolate_region_range[[1]],
      do_bc = config$irp_preprocess$do_bc,
      bc_method = config$irp_preprocess$bc_method,
      bc_cutoff = config$irp_preprocess$bc_cutoff,
      bc_do_impute = config$irp_preprocess$bc_do_impute,
      do_smooth = config$irp_preprocess$do_smooth,
      smooth_method = config$irp_preprocess$smooth_method,
      smooth_p = config$irp_preprocess$smooth_p,
      smooth_n = config$irp_preprocess$smooth_n,
      smooth_m = config$irp_preprocess$smooth_m,
      smooth_ts = config$irp_preprocess$smooth_ts,
      smooth_k = config$irp_preprocess$smooth_k,
      do_normalise = config$irp_preprocess$do_normalise,
      normalise_method = config$irp_preprocess$normalise_method,
      do_bin = config$irp_preprocess$do_bin,
      bin_width = config$irp_preprocess$bin_width,
      bin_new_x_type = config$irp_preprocess$bin_new_x_type,
      do_scale = FALSE, #config$irp_preprocess$do_scale,
      scale_center = config$irp_preprocess$scale_center,
      scale_scale = config$irp_preprocess$scale_scale,
      do_return_as_ir = TRUE
    ) |>
    ir::ir_clip(range = data.frame(start = c(650, 2400), end = c(2250, 4000))) |>
    ir::ir_scale(center = config$irp_preprocess$scale_center, scale = config$irp_preprocess$scale_scale)

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


#' Preprocesses spectra according to a preprocessing config object
#'
#' @param x An ir object to be preprocessed.
#'
#' @param config A list with arguments compatible with `irpeat::irp_preprocess()`.
#'
#' @keywords Internal
#' @noRd
irp_preprocess_eb1149 <- function (x, config) {

  res <- x

  irpeat::irp_preprocess(
    x = res,
    do_interpolate = config$irp_preprocess$do_interpolate,
    interpolate_start = config$irp_preprocess$interpolate_start[[1]],
    interpolate_dw = config$irp_preprocess$interpolate_dw,
    do_clip = config$irp_preprocess$do_clip[[1]],
    clip_range = config$irp_preprocess$clip_range[[1]],
    do_interpolate_region = FALSE,
    interpolate_region_range = config$irp_preprocess$interpolate_region_range[[1]],
    do_bc = config$irp_preprocess$do_bc,
    bc_method = config$irp_preprocess$bc_method,
    bc_cutoff = config$irp_preprocess$bc_cutoff,
    bc_do_impute = config$irp_preprocess$bc_do_impute,
    do_smooth = config$irp_preprocess$do_smooth,
    smooth_method = config$irp_preprocess$smooth_method,
    smooth_p = config$irp_preprocess$smooth_p,
    smooth_n = config$irp_preprocess$smooth_n,
    smooth_m = config$irp_preprocess$smooth_m,
    smooth_ts = config$irp_preprocess$smooth_ts,
    smooth_k = config$irp_preprocess$smooth_k,
    do_normalise = config$irp_preprocess$do_normalise,
    normalise_method = config$irp_preprocess$normalise_method,
    do_bin = config$irp_preprocess$do_bin,
    bin_width = config$irp_preprocess$bin_width,
    bin_new_x_type = config$irp_preprocess$bin_new_x_type,
    do_scale = FALSE,
    scale_center = config$irp_preprocess$scale_center,
    scale_scale = config$irp_preprocess$scale_scale,
    do_return_as_ir = TRUE
  ) |>
    ir::ir_clip(range = data.frame(start = c(650, 2400), end = c(2250, 4000))) |>
    ir::ir_scale(center = config$irp_preprocess$scale_center, scale = config$irp_preprocess$scale_scale)

}


#### rvar ####

#' Sets units of the array of an rvar object
#'
#' @param ... Additional arguments passed to `units::set_units()`.
#'
#' @keywords Internal
#' @noRd
irp_set_units_rvar <- function(x, ...) {
  posterior::draws_of(x) <- units::set_units(posterior::draws_of(x), ...)
  x
}

#' Drops units of the array of an rvar object
#'
#' @keywords Internal
#' @noRd
irp_drop_units_rvar <- function(x) {
  posterior::draws_of(x) <- units::drop_units(posterior::draws_of(x))
  x
}
