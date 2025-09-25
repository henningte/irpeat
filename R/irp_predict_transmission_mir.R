#' Predicts peat properties from transmission mid-infrared spectra
#'
#' Functions to predict peat properties from transmission mid-infrared spectra.
#' All functions below have been computed using peat samples. For detailed
#' information on the underlying prediction models, see the details section.
#'
#' @include utils.R
#'
#' @name irp-predict-transmission-mir
#'
#' @details
#' The models use the models of the same name in the 'irpeatmodels' package. The
#' 'irpeatmodels' package provides information on the models.
#'
#' @param x An object of class [`ir`][ir::ir_new_ir] with transmission mid
#' infrared spectra. Some tests are applied to check if the supplied spectra
#' match the spectra used to fit the models (the spectral range is checked). The
#' spectral resolution of the original spectral data should not be smaller than
#' 4 cm\eqn{^{-1}} and it is not checked if this assumption is met. For the
#' following models, `x` has special meaning:
#' \describe{
#'   \item{`irp_microbial_nitrogen_content_1()`}{Here, `x` is a set of litter
#'   spectra after decomposition (and `y` is a set of litter spectra before
#'   decomposition). See the Details section.}
#' }
#'
#' @param y An object of class [`ir`][ir::ir_new_ir] with transmission mid
#' infrared spectra. This argument is required for the following models which
#' need more than one set of spectra to make predictions:
#' \describe{
#'   \item{`irp_microbial_nitrogen_content_1()`}{Here, `y` is a set of litter
#'   spectra before decomposition.}
#' }
#'
#' @param temperature For `irp_specific_heat_capacity_1()`: The temperature in K
#' for which to predict the specific heat capacity.
#'
#' @param bulk_density For `irp_porosity_1()`, `irp_non_macroporosity_1()`,
#' `irp_macroporosity_1()`, `irp_volume_fraction_of_solids_1()`,
#' `irp_saturated_hydraulic_conductivity_1()`,
#' `irp_dry_thermal_conductivity_1()`: One of:
#' \describe{
#'   \item{1.}{A numeric vector with the same number of elements as spectra in
#'   `x` with values for the dry bulk density in g cm\eqn{^{-3}}. These values
#'   will be used to predict the peat property.}
#'   \item{2.}{A list with the same number of elements as spectra in `x`. Each
#'   element must be a numeric vector with the same number of elements as there
#'   are MCMC draws in the corresponding model to use. Each numeric value is the
#'   the dry bulk density in g cm\eqn{^{-3}}. These values will be used to
#'   predict the peat property.}
#'   \item{3.}{`NULL`: Dry bulk density will be estimated from the spectra in `x`
#'   and these estimates will be used to predict the peat property.}
#' }
#'
#' @param nitrogen_content For `irp_specific_heat_capacity_1()`: One of:
#' \describe{
#'   \item{1.}{A numeric vector with the same number of elements as spectra in
#'   `x` with values for the nitrogen content in g g\eqn{^{-1}}. These values
#'   will be used to predict the peat property.}
#'   \item{2.}{A list with the same number of elements as spectra in `x`. Each
#'   element must be a numeric vector with the same number of elements as there
#'   are MCMC draws in the corresponding model to use. Each numeric value is the
#'   the nitrogen content in g g\eqn{^{-1}}. These values will be used to
#'   predict the peat property.}
#'   \item{3.}{`NULL`: Nitrogen content will be estimated from the spectra in `x`
#'   and these estimates will be used to predict the peat property.}
#' }
#'
#' @param ... Additional arguments passed to
#' [rstanarm::posterior_predict.stanreg()] (`irp_eac_1()`,`irp_eac_2()`).
#'
#' @param do_summary A logical value indicating if the predicted values should
#' be returned in a summarized version (`TRUE`) or not (`FALSE`).
#' \itemize{
#'   \item If `do_summary = FALSE`, a list column is returned and each
#'   element of the list column is a numeric vector, or an `rvar` object is
#'   returned as column in `x`, depending on the value of `return_as_list`. In
#'   both cases, the column contains draws from the posterior predictive
#'   distribution.
#'   \item If `do_summary = TRUE`, each element is a
#'   [quantities::quantities()] object with value and error summarized from
#'   posterior draws via `summary_function_mean` and `summary_function_sd`
#' }
#'
#' @param summary_function_mean A function used to summarize the predicted
#' values (average).
#'
#' @param summary_function_sd A function used to summarize the predicted
#' values (spread).
#'
#' @param return_as_list Logical value. If set to `TRUE`, the result will be
#' returned as list of draws, otherwise the result will be returned as `rvar`
#' object. This is a new argument currently only implemented for models
#' predicting the degree of decomposition.
#'
#' @param check_prediction_domain A character value indicating if and how it
#' should be checked whether the spectra in `x` are within the prediction domain
#' of the model. One of:
#' \describe{
#'   \item{`"train"`}{It is checked whether the spectra in `x` are within the
#'   prediction domain formed by the training data for the model.}
#'   \item{`"test"`}{It is checked whether the spectra in `x` are within the
#'   prediction domain formed by the testing data for the model.}
#'   \item{`"none"`}{It is not checked whether the spectra in `x` are within the
#'   prediction domain for the model.}
#' }
#'
#' @return `x` with a new column with the predicted peat property and a new
#' column with value `TRUE` if the respective spectrum is within the prediction
#' domain for the model and `FALSE` if not. If
#' `check_prediction_domain = "none"`, all values in this column are `NA`.
#'
#' @note
#' \describe{
#'   \item{`irp_eac_1()`, `irp_edc_1()`}{
#'     The model still has a relatively large uncertainty because it is fitted
#'     with few samples \insertCite{Teickner.2022}{irpeat}. For further
#'     limitations, see \insertCite{Teickner.2022;textual}{irpeat}.
#'   }
#'   \item{`irp_microbial_nitrogen_content_1()`}{
#'    \insertCite{Reuter.2020;textual}{irpeat} describes limitations and uncertainties:
#'    “Small method modifications should be considered for the applicability of
#'    the method in aerobic decomposition studies. These modifications include
#'    an optimization of the calibration curve, either through the addition of
#'    very low N litters to a decomposition study as calibration samples or
#'    through the artificial mixing of undecomposed litter with microbial
#'    biomass. Furthermore, the contribution of fungi must be considered, which
#'    we assumed to be negligible in anoxic soils. Differences in the amount of
#'    DNA per biomass units and in the C/N ratio should be considered for the
#'    decomposer biomass in aerobic systems. Finally, the applicability of the
#'    same calibration curve for decomposed litters of different plant species
#'    still has to be investigated.”
#'   }
#' }
#'
#' @source
#' \describe{
#'   \item{`irp_holocellulose_2()`, `irp_klason_lignin_2()`}{
#'     \insertCite{Teickner.2022f;textual}{irpeat}.
#'   }
#'   \item{`irp_eac_1()`, `irp_edc_1()`}{
#'     \insertCite{Teickner.2022;textual}{irpeat}.
#'   }
#'   \item{`irp_microbial_nitrogen_content_1()`}{
#'     \insertCite{Reuter.2020;textual}{irpeat}
#'   }
#'   \item{`irp_degree_of_decomposition_1()`, `irp_degree_of_decomposition_2()`, `irp_degree_of_decomposition_3()`}{
#'     \insertCite{Teickner.2025h;textual}{irpeat}
#'   }
#'   \item{All other models}{
#'     \insertCite{Teickner.2023;textual}{irpeat}
#'   }
#' }
#'
#' @seealso `irp_predict()`
#'
#' @references
#'   \insertAllCited{}
#'
NULL



#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' library(ir)
#'
#' x <- ir::ir_sample_data[1, ]
#'
#' ## make predictions
#'
#' # holocellulose content
#' x <- irpeat::irp_holocellulose_content_2(
#'   x,
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_holocellulose_content_2 <- function(x, ..., do_summary = FALSE, summary_function_mean = mean, summary_function_sd = stats::sd, return_as_list = FALSE, check_prediction_domain = "train") {

  check_irpeatmodels(version = "0.0.0")
  if(! requireNamespace("brms", quietly = TRUE)) {
    rlang::abort("You have to install the 'brms' package to use this function.")
  }
  stopifnot(inherits(x, "ir"))
  stopifnot(is.logical(do_summary) && length(do_summary) == 1)

  x_or <- x

  # get data
  m <- irpeatmodels::model_holocellulose_content_2
  config <- irpeatmodels::model_holocellulose_content_2_config

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

  # check prediction domain
  prediction_domain <-
    switch(
      check_prediction_domain,
      "train" = irpeatmodels::model_holocellulose_content_2_prediction_domain$train,
      "test" = irpeatmodels::model_holocellulose_content_2_prediction_domain$train,
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

  # reformat for predictions
  res <- ir::ir_flatten(x)
  res_colnames <- paste0("V", res[, 1, drop= TRUE])
  res <- as.data.frame(t(res[, -1, drop = FALSE]))
  attr(res, "scaled:center") <- attr(x, "scaled:center")
  attr(res, "scaled:scale") <- attr(x, "scaled:scale")
  colnames(res) <- res_colnames

  # predict
  res <- as.data.frame(brms::posterior_predict(m, newdata = data.frame(x = I(as.matrix(res)), stringsAsFactors = FALSE), ...))
  res <- res * config$data_scale$y_scale + config$data_scale$y_center

  # summarize and add unit
  res <-
    irp_summarize_predictions(
      x = res,
      x_unit = "g/g",
      do_summary = do_summary,
      return_as_list = return_as_list,
      summary_function_mean = mean,
      summary_function_sd = stats::sd
    )

  x_or$holocellulose_content_2 <- res
  cbind(x_or, res_pd %>% stats::setNames(nm = "holocellulose_content_2_in_pd"))

}

#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # Klason lignin content
#' x <- irpeat::irp_klason_lignin_content_2(
#'   x,
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_klason_lignin_content_2 <- function(x, ..., do_summary = FALSE, summary_function_mean = mean, summary_function_sd = stats::sd, return_as_list = FALSE, check_prediction_domain = "train") {

  check_irpeatmodels(version = "0.0.0")
  if(! requireNamespace("brms", quietly = TRUE)) {
    rlang::abort("You have to install the 'brms' package to use this function.")
  }
  stopifnot(inherits(x, "ir"))
  stopifnot(is.logical(do_summary) && length(do_summary) == 1)

  x_or <- x

  # get data
  m <- irpeatmodels::model_klason_lignin_content_2
  config <- irpeatmodels::model_klason_lignin_content_2_config

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

  # check prediction domain
  prediction_domain <-
    switch(
      check_prediction_domain,
      "train" = irpeatmodels::model_klason_lignin_content_2_prediction_domain$train,
      "test" = irpeatmodels::model_klason_lignin_content_2_prediction_domain$train,
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

  # reformat for predictions
  res <- ir::ir_flatten(x)
  res_colnames <- paste0("V", res[, 1, drop= TRUE])
  res <- as.data.frame(t(res[, -1, drop = FALSE]))
  attr(res, "scaled:center") <- attr(x, "scaled:center")
  attr(res, "scaled:scale") <- attr(x, "scaled:scale")
  colnames(res) <- res_colnames

  # predict
  res <- as.data.frame(brms::posterior_predict(m, newdata = data.frame(x = I(as.matrix(res)), stringsAsFactors = FALSE), ...))
  res <- res * config$data_scale$y_scale + config$data_scale$y_center

  # summarize and add unit
  res <-
    irp_summarize_predictions(
      x = res,
      x_unit = "g/g",
      do_summary = do_summary,
      return_as_list = return_as_list,
      summary_function_mean = mean,
      summary_function_sd = stats::sd
    )

  x_or$klason_lignin_content_2 <- res
  cbind(x_or, res_pd %>% stats::setNames(nm = "klason_lignin_content_2_in_pd"))

}



#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # electron accepting capacity
#' x <- irpeat::irp_eac_1(
#'   x,
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_eac_1 <- function(x, ..., do_summary = FALSE, summary_function_mean = mean, summary_function_sd = stats::sd, return_as_list = FALSE, check_prediction_domain = "train") {

  check_irpeatmodels(version = "0.0.0")
  rlang::is_installed("rstantools")
  if(! requireNamespace("rstanarm", quietly = TRUE)) {
    rlang::abort("You have to install the 'rstanarm' package to use this function.")
  }
  stopifnot(inherits(x, "ir"))
  stopifnot(is.logical(do_summary) && length(do_summary) == 1)

  x <-
    x %>%
    dplyr::mutate(
      `__id_sample__` = seq_len(nrow(.))
    )

  x_or <- x

  # get data
  m <- irpeatmodels::model_eac_1
  config <- irpeatmodels::model_eac_1_config

  # check spectra
  x_flat <- ir::ir_flatten(x)
  if(x_flat$x[[1]] > config$irp_preprocess$clip_range$start) {
    rlang::warn(paste0("The minimum wavenumber value in `x` is ", x_flat$x[[1]], " , but should be ", config$irp_preprocess$clip_range$start, " or smaller."))
  }
  if(x_flat$x[[nrow(x_flat)]] < config$irp_preprocess$clip_range$end) {
    rlang::warn(paste0("The maximum wavenumber value in `x` is ", x_flat$x[[nrow(x)]], " , but should be ", config$irp_preprocess$clip_range$end, " or larger."))
  }

  # identify and discard empty spectra
  x <-
    x %>%
    dplyr::filter(! ir::ir_identify_empty_spectra(x))

  # preprocessing
  x <-
    irp_preprocess_eb1014(x, config = config)

  # check prediction domain
  prediction_domain <-
    switch(
      check_prediction_domain,
      "train" = irpeatmodels::model_eac_1_prediction_domain$train,
      "test" = irpeatmodels::model_eac_1_prediction_domain$train,
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

  # reformat for predictions
  res <- ir::ir_flatten(x)
  res_colnames <- paste0("V", res[, 1, drop= TRUE])
  res <- as.data.frame(t(res[, -1, drop = FALSE]))
  attr(res, "scaled:center") <- attr(x, "scaled:center")
  attr(res, "scaled:scale") <- attr(x, "scaled:scale")
  colnames(res) <- res_colnames

  # predict
  res <- as.data.frame(rstantools::posterior_predict(m, newdata = data.frame(x = I(as.matrix(res)), stringsAsFactors = FALSE), ...))
  res <- res * config$data_scale$y_scale + config$data_scale$y_center

  # summarize and add unit
  res <-
    irp_summarize_predictions(
      x = res,
      x_unit = "umol/g",
      do_summary = do_summary,
      return_as_list = return_as_list,
      summary_function_mean = mean,
      summary_function_sd = stats::sd
    )

  # combine results
  dplyr::left_join(
    x_or,
    x %>%
      dplyr::select("__id_sample__") %>%
      dplyr::mutate(
        eac_1 = res
      ),
    by = "__id_sample__"
  ) %>%
    dplyr::left_join(
      x %>%
        dplyr::select("__id_sample__") %>%
        dplyr::mutate(
          eac_1_in_pd = res_pd$y
        ),
      by = "__id_sample__"
    ) %>%
    dplyr::select(-"__id_sample__")

}


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # electron donating capacity
#' x <- irpeat::irp_edc_1(
#'   x,
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_edc_1 <- function(x, ..., do_summary = FALSE, summary_function_mean = mean, summary_function_sd = stats::sd, return_as_list = FALSE, check_prediction_domain = "train") {

  check_irpeatmodels(version = "0.0.0")
  if(! requireNamespace("rstanarm", quietly = TRUE)) {
    rlang::abort("You have to install the 'rstanarm' package to use this function.")
  }
  stopifnot(inherits(x, "ir"))
  stopifnot(is.logical(do_summary) && length(do_summary) == 1)

  x <-
    x %>%
    dplyr::mutate(
      `__id_sample__` = seq_len(nrow(.))
    )

  x_or <- x

  # get data
  m <- irpeatmodels::model_edc_1
  config <- irpeatmodels::model_edc_1_config

  # check spectra
  x_flat <- ir::ir_flatten(x)
  if(x_flat$x[[1]] > config$irp_preprocess$clip_range$start) {
    rlang::warn(paste0("The minimum wavenumber value in `x` is ", x_flat$x[[1]], " , but should be ", config$irp_preprocess$clip_range$start, " or smaller."))
  }
  if(x_flat$x[[nrow(x_flat)]] < config$irp_preprocess$clip_range$end) {
    rlang::warn(paste0("The maximum wavenumber value in `x` is ", x_flat$x[[nrow(x)]], " , but should be ", config$irp_preprocess$clip_range$end, " or larger."))
  }

  # identify and discard empty spectra
  x <-
    x %>%
    dplyr::filter(! ir::ir_identify_empty_spectra(x))

  # preprocessing
  x <-
    irp_preprocess_eb1014(x, config = config)

  # check prediction domain
  prediction_domain <-
    switch(
      check_prediction_domain,
      "train" = irpeatmodels::model_edc_1_prediction_domain$train,
      "test" = irpeatmodels::model_edc_1_prediction_domain$train,
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

  # reformat for predictions
  res <- ir::ir_flatten(x)
  res_colnames <- paste0("V", res[, 1, drop= TRUE])
  res <- as.data.frame(t(res[, -1, drop = FALSE]))
  attr(res, "scaled:center") <- attr(x, "scaled:center")
  attr(res, "scaled:scale") <- attr(x, "scaled:scale")
  colnames(res) <- res_colnames

  # predict
  res <- as.data.frame(rstantools::posterior_predict(m, newdata = data.frame(x = I(as.matrix(res)), stringsAsFactors = FALSE), ...))
  res <- res * config$data_scale$y_scale + config$data_scale$y_center

  # summarize and add unit
  res <-
    irp_summarize_predictions(
      x = res,
      x_unit = "umol/g",
      do_summary = do_summary,
      return_as_list = return_as_list,
      summary_function_mean = mean,
      summary_function_sd = stats::sd
    )

  # combine results
  dplyr::left_join(
    x_or,
    x %>%
      dplyr::select("__id_sample__") %>%
      dplyr::mutate(
        edc_1 = res
      ),
    by = "__id_sample__"
  ) %>%
    dplyr::left_join(
      x %>%
        dplyr::select("__id_sample__") %>%
        dplyr::mutate(
          edc_1_in_pd = res_pd$y
        ),
      by = "__id_sample__"
    ) %>%
    dplyr::select(-"__id_sample__")

}




#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # carbon content
#' x <-
#'   irp_carbon_content_1(
#'     x,
#'     do_summary = TRUE,
#'     check_prediction_domain = "train"
#'   )
#'
#' @export
irp_carbon_content_1 <-
  irp_function_factory_eb1079(
    target_variable = "carbon_content_1",
    model = readRDS(system.file("extdata", paste0("model_", "carbon_content_1", ".rds"), package = "irpeatmodels")),
    config = readRDS(system.file("extdata", paste0("model_", "carbon_content_1", "_config.rds"), package = "irpeatmodels")),
    prediction_domain = readRDS(system.file("extdata", paste0("model_", "carbon_content_1", "_prediction_domain.rds"), package = "irpeatmodels")),
    irpeatmodels_required_version = "0.0.0"
  )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # nitrogen content
#' irpeat::irp_nitrogen_content_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_nitrogen_content_1 <-
  irp_function_factory_eb1079(
    target_variable = "nitrogen_content_1",
    model = readRDS(system.file("extdata", paste0("model_", "nitrogen_content_1", ".rds"), package = "irpeatmodels")),
    config = readRDS(system.file("extdata", paste0("model_", "nitrogen_content_1", "_config.rds"), package = "irpeatmodels")),
    prediction_domain = readRDS(system.file("extdata", paste0("model_", "nitrogen_content_1", "_prediction_domain.rds"), package = "irpeatmodels")),
    irpeatmodels_required_version = "0.0.0"
  )

#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # hydrogen content
#' irpeat::irp_hydrogen_content_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_hydrogen_content_1 <-
  irp_function_factory_eb1079(
    target_variable = "hydrogen_content_1",
    model = readRDS(system.file("extdata", paste0("model_", "hydrogen_content_1", ".rds"), package = "irpeatmodels")),
    config = readRDS(system.file("extdata", paste0("model_", "hydrogen_content_1", "_config.rds"), package = "irpeatmodels")),
    prediction_domain = readRDS(system.file("extdata", paste0("model_", "hydrogen_content_1", "_prediction_domain.rds"), package = "irpeatmodels")),
    irpeatmodels_required_version = "0.0.0"
  )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # oxygen content
#' irpeat::irp_oxygen_content_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_oxygen_content_1 <-
  irp_function_factory_eb1079(
    target_variable = "oxygen_content_1",
    model = readRDS(system.file("extdata", paste0("model_", "oxygen_content_1", ".rds"), package = "irpeatmodels")),
    config = readRDS(system.file("extdata", paste0("model_", "oxygen_content_1", "_config.rds"), package = "irpeatmodels")),
    prediction_domain = readRDS(system.file("extdata", paste0("model_", "oxygen_content_1", "_prediction_domain.rds"), package = "irpeatmodels")),
    irpeatmodels_required_version = "0.0.0"
  )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # phosphorus content
#' irpeat::irp_phosphorus_content_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_phosphorus_content_1 <-
  irp_function_factory_eb1079(
    target_variable = "phosphorus_content_1",
    model = readRDS(system.file("extdata", paste0("model_", "phosphorus_content_1", ".rds"), package = "irpeatmodels")),
    config = readRDS(system.file("extdata", paste0("model_", "phosphorus_content_1", "_config.rds"), package = "irpeatmodels")),
    prediction_domain = readRDS(system.file("extdata", paste0("model_", "phosphorus_content_1", "_prediction_domain.rds"), package = "irpeatmodels")),
    irpeatmodels_required_version = "0.0.0"
  )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # potassium content
#' irpeat::irp_potassium_content_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_potassium_content_1 <-
  irp_function_factory_eb1079(
    target_variable = "potassium_content_1",
    model = readRDS(system.file("extdata", paste0("model_", "potassium_content_1", ".rds"), package = "irpeatmodels")),
    config = readRDS(system.file("extdata", paste0("model_", "potassium_content_1", "_config.rds"), package = "irpeatmodels")),
    prediction_domain = readRDS(system.file("extdata", paste0("model_", "potassium_content_1", "_prediction_domain.rds"), package = "irpeatmodels")),
    irpeatmodels_required_version = "0.0.0"
  )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # sulfur content
#' irpeat::irp_sulfur_content_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_sulfur_content_1 <-
  irp_function_factory_eb1079(
    target_variable = "sulfur_content_1",
    model = readRDS(system.file("extdata", paste0("model_", "sulfur_content_1", ".rds"), package = "irpeatmodels")),
    config = readRDS(system.file("extdata", paste0("model_", "sulfur_content_1", "_config.rds"), package = "irpeatmodels")),
    prediction_domain = readRDS(system.file("extdata", paste0("model_", "sulfur_content_1", "_prediction_domain.rds"), package = "irpeatmodels")),
    irpeatmodels_required_version = "0.0.0"
  )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # titanium content
#' irpeat::irp_titanium_content_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_titanium_content_1 <-
  irp_function_factory_eb1079(
    target_variable = "titanium_content_1",
    model = readRDS(system.file("extdata", paste0("model_", "titanium_content_1", ".rds"), package = "irpeatmodels")),
    config = readRDS(system.file("extdata", paste0("model_", "titanium_content_1", "_config.rds"), package = "irpeatmodels")),
    prediction_domain = readRDS(system.file("extdata", paste0("model_", "titanium_content_1", "_prediction_domain.rds"), package = "irpeatmodels")),
    irpeatmodels_required_version = "0.0.0"
  )

#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # silicon content
#' irpeat::irp_silicon_content_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_silicon_content_1 <-
  irp_function_factory_eb1079(
    target_variable = "silicon_content_1",
    model = readRDS(system.file("extdata", paste0("model_", "silicon_content_1", ".rds"), package = "irpeatmodels")),
    config = readRDS(system.file("extdata", paste0("model_", "silicon_content_1", "_config.rds"), package = "irpeatmodels")),
    prediction_domain = readRDS(system.file("extdata", paste0("model_", "silicon_content_1", "_prediction_domain.rds"), package = "irpeatmodels")),
    irpeatmodels_required_version = "0.0.0"
  )

#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # calcium content
#' irpeat::irp_calcium_content_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_calcium_content_1 <-
  irp_function_factory_eb1079(
    target_variable = "calcium_content_1",
    model = readRDS(system.file("extdata", paste0("model_", "calcium_content_1", ".rds"), package = "irpeatmodels")),
    config = readRDS(system.file("extdata", paste0("model_", "calcium_content_1", "_config.rds"), package = "irpeatmodels")),
    prediction_domain = readRDS(system.file("extdata", paste0("model_", "calcium_content_1", "_prediction_domain.rds"), package = "irpeatmodels")),
    irpeatmodels_required_version = "0.0.0"
  )

#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # d13C values
#' irpeat::irp_d13C_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_d13C_1 <-
  irp_function_factory_eb1079(
    target_variable = "d13C_1",
    model = readRDS(system.file("extdata", paste0("model_", "d13C_1", ".rds"), package = "irpeatmodels")),
    config = readRDS(system.file("extdata", paste0("model_", "d13C_1", "_config.rds"), package = "irpeatmodels")),
    prediction_domain = readRDS(system.file("extdata", paste0("model_", "d13C_1", "_prediction_domain.rds"), package = "irpeatmodels")),
    irpeatmodels_required_version = "0.0.0"
  )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # d15N values
#' irpeat::irp_d15N_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_d15N_1 <-
  irp_function_factory_eb1079(
    target_variable = "d15N_1",
    model = readRDS(system.file("extdata", paste0("model_", "d15N_1", ".rds"), package = "irpeatmodels")),
    config = readRDS(system.file("extdata", paste0("model_", "d15N_1", "_config.rds"), package = "irpeatmodels")),
    prediction_domain = readRDS(system.file("extdata", paste0("model_", "d15N_1", "_prediction_domain.rds"), package = "irpeatmodels")),
    irpeatmodels_required_version = "0.0.0"
  )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # nominal oxidation state of carbon (NOSC)
#' irpeat::irp_nosc_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_nosc_1 <-
  irp_function_factory_eb1079(
    target_variable = "nosc_1",
    model = readRDS(system.file("extdata", paste0("model_", "nosc_1", ".rds"), package = "irpeatmodels")),
    config = readRDS(system.file("extdata", paste0("model_", "nosc_1", "_config.rds"), package = "irpeatmodels")),
    prediction_domain = readRDS(system.file("extdata", paste0("model_", "nosc_1", "_prediction_domain.rds"), package = "irpeatmodels")),
    irpeatmodels_required_version = "0.0.0"
  )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # Gibbs free energy of formation
#' irpeat::irp_dgf0_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_dgf0_1 <-
  target_variable <-
    irp_function_factory_eb1079(
      target_variable = "dgf0_1",
      model = readRDS(system.file("extdata", paste0("model_", "dgf0_1", ".rds"), package = "irpeatmodels")),
      config = readRDS(system.file("extdata", paste0("model_", "dgf0_1", "_config.rds"), package = "irpeatmodels")),
      prediction_domain = readRDS(system.file("extdata", paste0("model_", "dgf0_1", "_prediction_domain.rds"), package = "irpeatmodels")),
      irpeatmodels_required_version = "0.0.0"
    )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # dry bulk density
#' irpeat::irp_bulk_density_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_bulk_density_1 <-
  irp_function_factory_eb1079(
    target_variable = "bulk_density_1",
    model = readRDS(system.file("extdata", paste0("model_", "bulk_density_1", ".rds"), package = "irpeatmodels")),
    config = readRDS(system.file("extdata", paste0("model_", "bulk_density_1", "_config.rds"), package = "irpeatmodels")),
    prediction_domain = readRDS(system.file("extdata", paste0("model_", "bulk_density_1", "_prediction_domain.rds"), package = "irpeatmodels")),
    irpeatmodels_required_version = "0.0.0"
  )

#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # loss on ignition
#' irp_loss_on_ignition_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_loss_on_ignition_1 <-
  irp_function_factory_eb1079(
    target_variable = "loss_on_ignition_1",
    model = readRDS(system.file("extdata", paste0("model_", "loss_on_ignition_1", ".rds"), package = "irpeatmodels")),
    config = readRDS(system.file("extdata", paste0("model_", "loss_on_ignition_1", "_config.rds"), package = "irpeatmodels")),
    prediction_domain = readRDS(system.file("extdata", paste0("model_", "loss_on_ignition_1", "_prediction_domain.rds"), package = "irpeatmodels")),
    irpeatmodels_required_version = "0.0.0"
  )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # O/C
#' irpeat::irp_O_to_C_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_O_to_C_1 <-
  irp_function_factory_eb1079(
    target_variable = "O_to_C_1",
    model = readRDS(system.file("extdata", paste0("model_", "O_to_C_1", ".rds"), package = "irpeatmodels")),
    config = readRDS(system.file("extdata", paste0("model_", "O_to_C_1", "_config.rds"), package = "irpeatmodels")),
    prediction_domain = readRDS(system.file("extdata", paste0("model_", "O_to_C_1", "_prediction_domain.rds"), package = "irpeatmodels")),
    irpeatmodels_required_version = "0.0.0"
  )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # C/N
#' irpeat::irp_C_to_N_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_C_to_N_1 <-
  irp_function_factory_eb1079(
    target_variable = "C_to_N_1",
    model = readRDS(system.file("extdata", paste0("model_", "C_to_N_1", ".rds"), package = "irpeatmodels")),
    config = readRDS(system.file("extdata", paste0("model_", "C_to_N_1", "_config.rds"), package = "irpeatmodels")),
    prediction_domain = readRDS(system.file("extdata", paste0("model_", "C_to_N_1", "_prediction_domain.rds"), package = "irpeatmodels")),
    irpeatmodels_required_version = "0.0.0"
  )


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # H/C
#' irpeat::irp_H_to_C_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_H_to_C_1 <-
  irp_function_factory_eb1079(
    target_variable = "H_to_C_1",
    model = readRDS(system.file("extdata", paste0("model_", "H_to_C_1", ".rds"), package = "irpeatmodels")),
    config = readRDS(system.file("extdata", paste0("model_", "H_to_C_1", "_config.rds"), package = "irpeatmodels")),
    prediction_domain = readRDS(system.file("extdata", paste0("model_", "H_to_C_1", "_prediction_domain.rds"), package = "irpeatmodels")),
    irpeatmodels_required_version = "0.0.0"
  )



#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # volume fraction of solids, macroporosity, non-macroporosity
#' irpeat::irp_porosity_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' irpeat::irp_porosity_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train",
#'   bulk_density = 0.04
#' )
#'
#' @noRd
#' @keywords internal
irp_porosity_1 <- function(x, do_summary = FALSE, summary_function_mean = mean, summary_function_sd = stats::sd, check_prediction_domain = "train", return_as_list = FALSE, bulk_density = NULL) {

  bulk_density_1_in_pd <- NULL

  # check additional packages
  if(! requireNamespace("brms", quietly = TRUE)) {
    rlang::abort("You have to install the 'brms' package to use this function.")
  }
  if(!(is.null(bulk_density) || length(bulk_density) == nrow(x))) {
    rlang::abort("`bulk_density` must either be `NULL` or a numeric vector or list with an element for each row in `x`.")
  }

  x_or <- x

  # import data
  m_draws_porosity_1 <- irpeatmodels::model_porosity_1_draws
  config_porosity_1 <-  irpeatmodels::model_porosity_1_config

  # predict bulk density
  if(! is.null(bulk_density)) {
    if(requireNamespace("posterior", quietly = TRUE) && posterior::is_rvar(bulk_density)) {
      cond <- nrow(posterior::draws_of(bulk_density)) != nrow(m_draws_porosity_1)
      if(cond) {
        rlang::abort(paste0("If `bulk_density` is provided as `rvar`, it must contain ", nrow(m_draws_porosity_1), " draws."))
      }
      x[["bulk_density_1"]] <-
        posterior::draws_of(bulk_density) |>
        as.data.frame() |>
        purrr::map(function(.x) .x)
    } else if(is.list(bulk_density) && ! inherits(bulk_density, "rvar")) {
      cond <-
        purrr::map_lgl(bulk_density, function(.x) {
          length(.x) != nrow(m_draws_porosity_1)
        })
      if(any(cond)) {
        rlang::abort(paste0("If `bulk_density` is provided as list, it must contain ", nrow(m_draws_porosity_1), " values."))
      }
      x[["bulk_density_1"]] <- bulk_density
    } else {
      x[["bulk_density_1"]] <-
        purrr::map(bulk_density, function(.x) {
          rep(.x, nrow(irpeatmodels::model_porosity_1_draws))
        })
    }
    x[["non_macroporosity_1_in_pd"]] <- x[["macroporosity_1_in_pd"]] <- x[["volume_fraction_solids_1_in_pd"]] <- NA
  } else {
    x <-
      x %>%
      dplyr::select(! dplyr::any_of(c("bulk_density_1", "bulk_density_1_in_pd"))) %>%
      irp_bulk_density_1(check_prediction_domain = check_prediction_domain, do_summary = FALSE, return_as_list = TRUE) %>%
      dplyr::mutate(
        non_macroporosity_1_in_pd = bulk_density_1_in_pd,
        macroporosity_1_in_pd = .data$non_macroporosity_1_in_pd,
        volume_fraction_solids_1_in_pd = .data$non_macroporosity_1_in_pd
      )
  }


  ## predict porosity

  # get predictor matrix
  X <-
    as.data.frame(x$bulk_density_1) %>%
    as.matrix()
  logX <- log(X)
  X <- X - config_porosity_1$data_scale$x_center["b_bulkdensity"]
  logX <- logX - config_porosity_1$data_scale$x_center["b_logbulkdensity"]

  Intercept_non_macroporosity <-
    m_draws_porosity_1$b_munonmacroporosity_Intercept +
    m_draws_porosity_1$b_munonmacroporosity_bulk_density * config_porosity_1$data_scale$x_center[["b_bulkdensity"]] +
    m_draws_porosity_1$b_munonmacroporosity_logbulk_density * config_porosity_1$data_scale$x_center[["b_logbulkdensity"]]

  Intercept_macroporosity <-
    m_draws_porosity_1$b_mumacroporosity_Intercept +
    m_draws_porosity_1$b_mumacroporosity_bulk_density * config_porosity_1$data_scale$x_center[["b_bulkdensity"]] +
    m_draws_porosity_1$b_mumacroporosity_logbulk_density * config_porosity_1$data_scale$x_center[["b_logbulkdensity"]]


  # linear predictor
  mu <-
    list(
      non_macroporosity =
        Intercept_non_macroporosity +
        sweep(X, 1, (m_draws_porosity_1 %>% dplyr::pull(.data$b_munonmacroporosity_bulk_density)), FUN = "*") +
        sweep(logX, 1, (m_draws_porosity_1 %>% dplyr::pull(.data$b_munonmacroporosity_logbulk_density)), FUN = "*"),
      macroporosity =
        Intercept_macroporosity +
        sweep(X, 1, (m_draws_porosity_1 %>% dplyr::pull(.data$b_mumacroporosity_bulk_density)), FUN = "*") +
        sweep(logX, 1, (m_draws_porosity_1 %>% dplyr::pull(.data$b_mumacroporosity_logbulk_density)), FUN = "*")
    )
  mu <- simplify2array(mu)
  mu_inv <- brms:::inv_link_categorical(mu)
  mu_inv[mu_inv <= .Machine$double.eps] <- .Machine$double.eps

  # predictions
  res <-
    purrr::map(seq_len(dim(mu_inv)[[2]]), function(i) {
      brms::rdirichlet(n = dim(mu_inv)[[1]], alpha = mu_inv[, i, ] * m_draws_porosity_1$phi) %>%
        as.data.frame() %>%
        stats::setNames(nm = c("volume_fraction_solids_1", "non_macroporosity_1", "macroporosity_1"))
    }) %>%
    purrr::transpose() %>%
    purrr::map2_dfc(names(.), function(.x, .y) {
      tibble::tibble(
        y =
          irp_summarize_predictions(
            x = as.data.frame(.x),
            x_unit = "L/L",
            do_summary = do_summary,
            return_as_list = return_as_list,
            summary_function_mean = summary_function_mean,
            summary_function_sd = summary_function_sd
          ) %>%
          unname()
      ) %>%
        stats::setNames(nm = .y)
    })

  x_or[, c("volume_fraction_solids_1", "non_macroporosity_1", "macroporosity_1")] <- res
  x_or[, paste0(c("volume_fraction_solids_1", "non_macroporosity_1", "macroporosity_1"), "_in_pd")] <- x %>% dplyr::select(.data$volume_fraction_solids_1_in_pd, .data$non_macroporosity_1_in_pd, .data$macroporosity_1_in_pd)
  x_or

}


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # volume fraction of solids
#' irpeat::irp_volume_fraction_solids_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' irpeat::irp_volume_fraction_solids_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train",
#'   bulk_density = 0.04
#' )
#'
#' @export
irp_volume_fraction_solids_1 <- function(x, do_summary = FALSE, summary_function_mean = mean, summary_function_sd = stats::sd, check_prediction_domain = "train", return_as_list = FALSE, bulk_density = NULL) {

  dplyr::bind_cols( #---note: to avoid conflicts when some of the porosity-related variables have already been compted for x
    x |>
      dplyr::select(! dplyr::any_of(c("volume_fraction_solids_1", "volume_fraction_solids_1_in_pd"))),
    irp_porosity_1(
      x = x %>% dplyr::select(.data$spectra),
      do_summary = do_summary,
      summary_function_mean = summary_function_mean,
      summary_function_sd = summary_function_sd,
      check_prediction_domain = check_prediction_domain,
      return_as_list = return_as_list,
      bulk_density = bulk_density
    ) %>%
      dplyr::select(.data$volume_fraction_solids_1, .data$volume_fraction_solids_1_in_pd)
  )

}

#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # non-macroporosity
#' irpeat::irp_non_macroporosity_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' irpeat::irp_non_macroporosity_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train",
#'   bulk_density = 0.04
#' )
#'
#' @export
irp_non_macroporosity_1 <- function(x, do_summary = FALSE, summary_function_mean = mean, summary_function_sd = stats::sd, check_prediction_domain = "train", return_as_list = FALSE, bulk_density = NULL) {

  dplyr::bind_cols( #---note: to avoid conflicts when some of the porosity-related variables have already been compted for x
    x |>
      dplyr::select(! dplyr::any_of(c("non_macroporosity_1", "non_macroporosity_1_in_pd"))),
    irp_porosity_1(
      x = x %>% dplyr::select(.data$spectra),
      do_summary = do_summary,
      summary_function_mean = summary_function_mean,
      summary_function_sd = summary_function_sd,
      check_prediction_domain = check_prediction_domain,
      return_as_list = return_as_list,
      bulk_density = bulk_density
    ) %>%
      dplyr::select(.data$non_macroporosity_1, .data$non_macroporosity_1_in_pd)
  )

}


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # macroporosity
#' irpeat::irp_macroporosity_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' irpeat::irp_macroporosity_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train",
#'   bulk_density = 0.04
#' )
#'
#' @export
irp_macroporosity_1 <- function(x, do_summary = FALSE, summary_function_mean = mean, summary_function_sd = stats::sd, check_prediction_domain = "train", return_as_list = FALSE, bulk_density = NULL) {

  dplyr::bind_cols( #---note: to avoid conflicts when some of the porosity-related variables have already been compted for x
    x |>
      dplyr::select(! dplyr::any_of(c("macroporosity_1", "macroporosity_1_in_pd"))),
    irp_porosity_1(
      x = x %>% dplyr::select(.data$spectra),
      do_summary = do_summary,
      summary_function_mean = summary_function_mean,
      summary_function_sd = summary_function_sd,
      check_prediction_domain = check_prediction_domain,
      return_as_list = return_as_list,
      bulk_density = bulk_density
    ) %>%
      dplyr::select(.data$macroporosity_1, .data$macroporosity_1_in_pd)
  )

}



#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # saturated hydraulic conductivity
#' irpeat::irp_saturated_hydraulic_conductivity_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' irpeat::irp_saturated_hydraulic_conductivity_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train",
#'   bulk_density = 0.04
#' )
#'
#' @export
irp_saturated_hydraulic_conductivity_1 <- function(x, do_summary = FALSE, summary_function_mean = mean, summary_function_sd = stats::sd, check_prediction_domain = "train", return_as_list = FALSE, bulk_density = NULL) {

  # check additional packages
  if(! requireNamespace("brms", quietly = TRUE)) {
    rlang::abort("You have to install the 'brms' package to use this function.")
  }
  if(!(is.null(bulk_density) || length(bulk_density) == nrow(x))) {
    rlang::abort("`bulk_density` must either be `NULL` or a numeric vector with an element for each row in `x`.")
  }

  x_or <- x
  x_has_bulk_density <- any(colnames(x) == "bulk_density_1")

  # import data
  m_draws_ks_1 <- irpeatmodels::model_saturated_hydraulic_conductivity_1_draws
  config_ks_1 <-  irpeatmodels::model_saturated_hydraulic_conductivity_1_config

  # predict bulk density
  if(!is.null(bulk_density)) {
    if (requireNamespace("posterior", quietly = TRUE) && posterior::is_rvar(bulk_density)) {
      cond <- nrow(posterior::draws_of(bulk_density)) != nrow(m_draws_ks_1)
      if(cond) {
        rlang::abort(paste0("If `bulk_density` is provided as `rvar`, it must contain ", nrow(m_draws_ks_1), " draws."))
      }
      x[["bulk_density_1"]] <-
        posterior::draws_of(bulk_density) |>
        as.data.frame() |>
        purrr::map(function(.x) .x)
    } else if(is.list(bulk_density)) {
      cond <-
        purrr::map_lgl(bulk_density, function(.x) {
          length(.x) != nrow(m_draws_ks_1)
        })
      if(any(cond)) {
        rlang::abort(paste0("If `bulk_density` is provided as list, it must contain ", nrow(m_draws_ks_1), " values."))
      }
      x[["bulk_density_1"]] <- bulk_density
    } else {
      x[["bulk_density_1"]] <-
        purrr::map(bulk_density, function(.x) {
          rep(.x, nrow(m_draws_ks_1))
        })
    }
    x[["saturated_hydraulic_conductivity_1_in_pd"]] <- NA
  } else {
    x <-
      x %>%
      dplyr::select(! dplyr::any_of(c("bulk_density_1", "bulk_density_1_in_pd"))) %>%
      irp_bulk_density_1(do_summary = FALSE, check_prediction_domain = check_prediction_domain, return_as_list = TRUE) %>%
      dplyr::rename(saturated_hydraulic_conductivity_1_in_pd = "bulk_density_1_in_pd")
  }

  ## predict saturated hydraulic conductivity

  # get predictor matrix
  X <-
    as.data.frame(x$bulk_density_1) %>%
    as.matrix()
  logX <- log(X)
  X <- X - config_ks_1$data_scale$x_center["b_bulkdensity"]
  logX <- logX - config_ks_1$data_scale$x_center["b_logbulkdensity"]

  # correct intercept for scaled predictor variables (brms returns the intercepts for the unscaled predictor variables)
  Intercept_mu <-
    m_draws_ks_1$b_Intercept +
    m_draws_ks_1$b_bulk_density * config_ks_1$data_scale$x_center[["b_bulkdensity"]] +
    m_draws_ks_1$b_logbulk_density * config_ks_1$data_scale$x_center[["b_logbulkdensity"]]

  Intercept_phi <-
    m_draws_ks_1$b_phi_Intercept +
    m_draws_ks_1$b_phi_bulk_density * config_ks_1$data_scale$x_center[["b_bulkdensity"]] +
    m_draws_ks_1$b_phi_logbulk_density * config_ks_1$data_scale$x_center[["b_logbulkdensity"]]

  # linear predictor
  mu <-
    Intercept_mu +
    sweep(X, 1, (m_draws_ks_1 %>% dplyr::pull(.data$b_bulk_density)), FUN = "*") +
    sweep(logX, 1, (m_draws_ks_1 %>% dplyr::pull(.data$b_logbulk_density)), FUN = "*")

  phi <-
    Intercept_phi +
    sweep(X, 1, (m_draws_ks_1 %>% dplyr::pull(.data$b_phi_bulk_density)), FUN = "*") +
    sweep(logX, 1, (m_draws_ks_1 %>% dplyr::pull(.data$b_phi_logbulk_density)), FUN = "*")

  mu <- config_ks_1$likelihood$linkinv(mu)
  phi <- exp(phi) #---note: in the model, phi has a log link

  # predictions
  res <-
    purrr::map_dfc(seq_len(ncol(mu)), function(i) {

      res <- stats::rbeta(n = nrow(mu), shape1 = mu[, i] * phi[, i], shape2 = (1 - mu[, i]) * phi[, i])

      # scale
      tibble::tibble(x = res * config_ks_1$data_scale$y_scale + config_ks_1$data_scale$y_center) %>%
        stats::setNames(nm = paste0("V", i))

    })

  # summarize and add unit
  res <-
    irp_summarize_predictions(
      x = res,
      x_unit = "cm/h",
      do_summary = do_summary,
      return_as_list = return_as_list,
      summary_function_mean = summary_function_mean,
      summary_function_sd = summary_function_sd
    )

  res <-
    tibble::tibble(y = res) %>%
    stats::setNames(nm = "saturated_hydraulic_conductivity_1")

  cbind(x_or, res, x %>% dplyr::select(.data$saturated_hydraulic_conductivity_1_in_pd))

}




#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # specific heat capacity
#' irpeat::irp_specific_heat_capacity_1(
#'   irpeat_sample_data[1, ],
#'   temperature = 290,
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' irpeat::irp_specific_heat_capacity_1(
#'   irpeat_sample_data[1, ],
#'   temperature = 290,
#'   do_summary = TRUE,
#'   check_prediction_domain = "train",
#'   nitrogen_content = irpeat_sample_data$N[1]
#' )
#'
#' @export
irp_specific_heat_capacity_1 <- function(x, temperature = 273.15, do_summary = FALSE, summary_function_mean = mean, summary_function_sd = stats::sd, check_prediction_domain = "train", return_as_list = FALSE, nitrogen_content = NULL) {

  stopifnot(temperature >= 0)

  # check additional packages
  if(! requireNamespace("brms", quietly = TRUE)) {
    rlang::abort("You have to install the 'brms' package to use this function.")
  }
  if(!(is.null(nitrogen_content) || length(nitrogen_content) == nrow(x))) {
    rlang::abort("`nitrogen_content` must either be `NULL` or a numeric vector with an element for each row in `x`.")
  }

  x_or <- x
  x_has_nitrogen_content <- any(colnames(x) == "nitrogen_content_1")

  # import data
  m_draws_cp_1 <- irpeatmodels::model_specific_heat_capacity_1_draws
  config_cp_1 <-  irpeatmodels::model_specific_heat_capacity_1_config

  # predict N
  if(!is.null(nitrogen_content)) {
    if (requireNamespace("posterior", quietly = TRUE) && posterior::is_rvar(nitrogen_content)) {
      cond <- nrow(posterior::draws_of(nitrogen_content)) != nrow(m_draws_cp_1)
      if(cond) {
        rlang::abort(paste0("If `nitrogen_content` is provided as `rvar`, it must contain ", nrow(m_draws_cp_1), " draws."))
      }
      x[["nitrogen_content_1"]] <-
        posterior::draws_of(nitrogen_content) |>
        as.data.frame() |>
        purrr::map(function(.x) .x)
    } else if(is.list(nitrogen_content) && ! inherits(nitrogen_content, "rvar")) {
      cond <-
        purrr::map_lgl(nitrogen_content, function(.x) {
          length(.x) != nrow(m_draws_cp_1)
        })
      if(any(cond)) {
        rlang::abort(paste0("If `nitrogen_content` is provided as list, it must contain ", nrow(m_draws_cp_1), " values."))
      }
      x[["nitrogen_content_1"]] <- nitrogen_content
    } else {
      x[["nitrogen_content_1"]] <-
        purrr::map(nitrogen_content, function(.x) {
          rep(.x, nrow(m_draws_cp_1))
        })
    }
    x[["specific_heat_capacity_1_in_pd"]] <- NA
  } else {
    x <-
      x %>%
      dplyr::select(! dplyr::any_of(c("nitrogen_content_1", "nitrogen_content_1_in_pd"))) %>%
      irp_nitrogen_content_1(do_summary = FALSE, check_prediction_domain = check_prediction_domain, return_as_list = TRUE) %>%
      dplyr::rename(specific_heat_capacity_1_in_pd = "nitrogen_content_1_in_pd")
  }

  ## predict specific heat capacity

  # get predictor matrix
  X_nitrogen_content_1 <-
    as.data.frame(x$nitrogen_content_1) %>%
    as.matrix() %>%
    magrittr::subtract(config_cp_1$data_scale$x_center) %>%
    magrittr::divide_by(config_cp_1$data_scale$x_scale["b_N"])
  X_temperature <-
    matrix(data = temperature - 273.15, nrow = nrow(X_nitrogen_content_1), ncol = ncol(X_nitrogen_content_1), byrow = FALSE) %>%
    magrittr::subtract(config_cp_1$data_scale$x_center) %>%
    magrittr::divide_by(config_cp_1$data_scale$x_scale["b_temperature"])
  X_nitrogen_content_1_by_temperature <-
    X_nitrogen_content_1 * X_temperature

  Intercept_mu <-
    m_draws_cp_1$b_Intercept +
    m_draws_cp_1$b_N * config_cp_1$model_scale$x_center[["b_N"]] +
    m_draws_cp_1$b_temperature * config_cp_1$model_scale$x_center[["b_temperature"]] +
    m_draws_cp_1$`b_N:temperature` * config_cp_1$model_scale$x_center[["b_N"]] * config_cp_1$model_scale$x_center[["b_temperature"]]

  # linear predictor
  mu <-
    Intercept_mu +
    sweep(X_nitrogen_content_1, 1, (m_draws_cp_1 %>% dplyr::pull(.data$b_N)), FUN = "*") +
    sweep(X_temperature, 1, (m_draws_cp_1 %>% dplyr::pull(.data$b_temperature)), FUN = "*") +
    sweep(X_nitrogen_content_1_by_temperature, 1, (m_draws_cp_1 %>% dplyr::pull(.data$`b_N:temperature`)), FUN = "*")

  mu <- as.data.frame(config_cp_1$likelihood$linkinv(mu))
  shape <- m_draws_cp_1$shape

  # predictions
  res <-
    purrr::map_dfc(seq_len(ncol(mu)), function(i) {

      res <- stats::rgamma(n = nrow(mu), shape = shape, rate = shape/mu[, i])

      # scale
      tibble::tibble(x = res * config_cp_1$data_scale$y_scale + config_cp_1$data_scale$y_center) %>%
        stats::setNames(nm = paste0("V", i))

    })

  # summarize and add unit
  res <-
    irp_summarize_predictions(
      x = res,
      x_unit = "J/(g * K)",
      do_summary = do_summary,
      return_as_list = return_as_list,
      summary_function_mean = summary_function_mean,
      summary_function_sd = summary_function_sd
    )

  res <-
    tibble::tibble(y = res) %>%
    stats::setNames(nm = "specific_heat_capacity_1")

  cbind(x_or, res, x %>% dplyr::select(.data$specific_heat_capacity_1_in_pd))

}


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # dry thermal conductivity
#' irpeat::irp_dry_thermal_conductivity_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' irpeat::irp_dry_thermal_conductivity_1(
#'   irpeat_sample_data[1, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train",
#'   bulk_density = 0.04
#' )
#'
#' @export
irp_dry_thermal_conductivity_1 <- function(x, do_summary = FALSE, summary_function_mean = mean, summary_function_sd = stats::sd, check_prediction_domain = "train", return_as_list = FALSE, bulk_density = NULL) {

  # check additional packages
  if(! requireNamespace("brms", quietly = TRUE)) {
    rlang::abort("You have to install the 'brms' package to use this function.")
  }
  if(!(is.null(bulk_density) || length(bulk_density) == nrow(x))) {
    rlang::abort("`bulk_density` must either be `NULL` or a numeric vector with an element for each row in `x`.")
  }

  x_or <- x
  x_has_bulk_density <- any(colnames(x) == "bulk_density_1")

  # import data
  m_draws_kt_1 <- irpeatmodels::model_dry_thermal_conductivity_1_draws
  config_kt_1 <-  irpeatmodels::model_dry_thermal_conductivity_1_config

  # predict bulk density
  if(!is.null(bulk_density)) {
    if (requireNamespace("posterior", quietly = TRUE) && posterior::is_rvar(bulk_density)) {
      cond <- nrow(posterior::draws_of(bulk_density)) != nrow(m_draws_kt_1)
      if(cond) {
        rlang::abort(paste0("If `bulk_density` is provided as `rvar`, it must contain ", nrow(m_draws_kt_1), " draws."))
      }
      x[["bulk_density_1"]] <-
        posterior::draws_of(bulk_density) |>
        as.data.frame() |>
        purrr::map(function(.x) .x)
    } else if(is.list(bulk_density) && ! inherits(bulk_density, "rvar")) {
      cond <-
        purrr::map_lgl(bulk_density, function(.x) {
          length(.x) != nrow(m_draws_kt_1)
        })
      if(any(cond)) {
        rlang::abort(paste0("If `bulk_density` is provided as list, it must contain ", nrow(m_draws_kt_1), " values."))
      }
      x[["bulk_density_1"]] <- bulk_density
    } else {
      x[["bulk_density_1"]] <-
        purrr::map(bulk_density, function(.x) {
          rep(.x, nrow(m_draws_kt_1))
        })
    }
    x[["dry_thermal_conductivity_1_in_pd"]] <- NA
  } else {
  x <-
    x %>%
    dplyr::select(! dplyr::any_of(c("bulk_density_1", "bulk_density_1_in_pd"))) %>%
    irp_bulk_density_1(do_summary = FALSE, check_prediction_domain = check_prediction_domain, return_as_list = TRUE) %>%
    dplyr::rename(dry_thermal_conductivity_1_in_pd = "bulk_density_1_in_pd")
  }

  ## predict dry thermal conductivity

  # get predictor matrix
  X_bulk_density_1 <-
    as.data.frame(x$bulk_density_1) %>%
    as.matrix() %>%
    magrittr::subtract(config_kt_1$data_scale$x_center) %>%
    magrittr::divide_by(config_kt_1$data_scale$x_scale)
  logX_bulk_density_1 <-log(X_bulk_density_1)

  Intercept_mu <-
    m_draws_kt_1$b_Intercept +
    m_draws_kt_1$b_bulk_density * config_kt_1$model_scale$x_center[["b_bulkdensity"]] +
    m_draws_kt_1$b_logbulk_density * config_kt_1$model_scale$x_center[["b_bulkdensity"]]

  Intercept_mu_shape <-
    m_draws_kt_1$b_shape_Intercept +
    m_draws_kt_1$b_shape_bulk_density * config_kt_1$model_scale$x_center[["b_bulkdensity"]] +
    m_draws_kt_1$b_shape_logbulk_density * config_kt_1$model_scale$x_center[["b_bulkdensity"]]

  # linear predictor
  mu <-
    Intercept_mu +
    sweep(X_bulk_density_1, 1, (m_draws_kt_1 %>% dplyr::pull(.data$b_bulk_density)), FUN = "*") +
    sweep(logX_bulk_density_1, 1, (m_draws_kt_1 %>% dplyr::pull(.data$b_logbulk_density)), FUN = "*")

  mu_shape <-
    Intercept_mu_shape +
    sweep(X_bulk_density_1, 1, (m_draws_kt_1 %>% dplyr::pull(.data$b_shape_bulk_density)), FUN = "*") +
    sweep(logX_bulk_density_1, 1, (m_draws_kt_1 %>% dplyr::pull(.data$b_shape_logbulk_density)), FUN = "*")

  mu <- as.data.frame(config_kt_1$likelihood$linkinv(mu))
  shape <- as.data.frame(exp(mu_shape))

  # predictions
  res <-
    purrr::map_dfc(seq_len(ncol(mu)), function(i) {

      res <- stats::rgamma(n = nrow(mu), shape = shape[, i], rate = shape[, i]/mu[, i])

      # scale
      tibble::tibble(x = res * config_kt_1$data_scale$y_scale + config_kt_1$data_scale$y_center) %>%
        stats::setNames(nm = paste0("V", i))

    })

  # summarize and add unit
  res <-
    irp_summarize_predictions(
      x = res,
      x_unit = "W/(m * K)",
      do_summary = do_summary,
      return_as_list = return_as_list,
      summary_function_mean = summary_function_mean,
      summary_function_sd = summary_function_sd
    )

  res <-
    tibble::tibble(y = res) %>%
    stats::setNames(nm = "dry_thermal_conductivity_1")

  cbind(x_or, res, x %>% dplyr::select(.data$dry_thermal_conductivity_1_in_pd))

}


#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # microbial nitrogen content (note that the data are not ideal uses cases for
#' # the model (see the documentation). The following only demonstrates how to
#' # use the prediction function)
#' irpeat::irp_microbial_nitrogen_content_1(
#'   x = irpeat_sample_data[1, ],
#'   y = irpeat_sample_data[2, ],
#'   do_summary = TRUE,
#'   check_prediction_domain = "train"
#' )
#'
#' @export
irp_microbial_nitrogen_content_1 <- function(x, y, do_summary = FALSE, summary_function_mean = mean, summary_function_sd = stats::sd, check_prediction_domain = "train", return_as_list = FALSE) {

  # check additional packages
  if(! requireNamespace("brms", quietly = TRUE)) {
    rlang::abort("You have to install the 'brms' package to use this function.")
  }

  x_or <- x

  # import data
  m_mbn_1 <- irpeatmodels::model_microbial_nitrogen_content_1
  config <-  irpeatmodels::model_microbial_nitrogen_content_1_config
  prediction_domain <- irpeatmodels::model_microbial_nitrogen_content_1_prediction_domain
  target_variable_name <- "microbial_nitrogen_content_1"

  # check spectra
  x_flat <- ir::ir_flatten(x)
  if(x_flat$x[[1]] > config$irp_preprocess$clip_range$start) {
    rlang::warn(paste0("The minimum wavenumber value in `x` is ", x_flat$x[[1]], " , but should be ", config$irp_preprocess$clip_range$start, " or smaller."))
  }
  if(x_flat$x[[nrow(x_flat)]] < config$irp_preprocess$clip_range$end) {
    rlang::warn(paste0("The maximum wavenumber value in `x` is ", x_flat$x[[nrow(x)]], " , but should be ", config$irp_preprocess$clip_range$end, " or larger."))
  }
  y_flat <- ir::ir_flatten(y)
  if(y_flat$x[[1]] > config$irp_preprocess$clip_range$start) {
    rlang::warn(paste0("The minimum wavenumber value in `y` is ", y_flat$x[[1]], " , but should be ", config$irp_preprocess$clip_range$start, " or smaller."))
  }
  if(y_flat$x[[nrow(y_flat)]] < config$irp_preprocess$clip_range$end) {
    rlang::warn(paste0("The maximum wavenumber value in `y` is ", y_flat$x[[nrow(x)]], " , but should be ", config$irp_preprocess$clip_range$end, " or larger."))
  }

  # preprocessing
  x <- irp_preprocess_for(x = x, y = y, variable = target_variable_name)

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
        microbial_nitrogen_content_1_in_pd =
          x %>%
          irp_is_in_prediction_domain(prediction_domain = prediction_domain) %>%
          dplyr::pull(.data$is_in_prediction_domain)
      )
    } else {
      tibble::tibble(
        microbial_nitrogen_content_1_in_pd = rep(NA, nrow(x_or))
      )
    }

  # predict microbial N content
  res <-
    brms::posterior_predict(
      m_mbn_1,
      newdata =
        tibble::tibble(
          intensity_1220 =
            purrr::map_dbl(x$spectra, function(.x) .x$y)
        )
    ) %>%
    as.data.frame()

  # summarize and add unit
  res <-
    irp_summarize_predictions(
      x = res,
      x_unit = "g/g",
      do_summary = do_summary,
      return_as_list = return_as_list,
      summary_function_mean = summary_function_mean,
      summary_function_sd = summary_function_sd
    )

  res <-
    tibble::tibble(y = res) %>%
    stats::setNames(nm = "microbial_nitrogen_content_1")

  cbind(x_or, res, res_pd %>% dplyr::select(.data$microbial_nitrogen_content_1_in_pd))

}


#### eb1149 ####

#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # degree_of_decomposition_1
#' if(! requireNamespace("posterior", quietly = TRUE)) {
#' x <-
#'   irpeat::irp_degree_of_decomposition_1(
#'     irpeat_sample_data[1, ],
#'     do_summary = TRUE,
#'     check_prediction_domain = "train",
#'     summary_function_sd = posterior::sd
#'   )
#' }
#'
#' @export
irp_degree_of_decomposition_1 <-
  irp_function_factory_eb1149(
    model = irpeatmodels::model_degree_of_decomposition_1_brms,
    config = irpeatmodels::model_degree_of_decomposition_1_config,
    prediction_domain = irpeatmodels::model_degree_of_decomposition_1_prediction_domain,
    target_variable_name =
      "degree_of_decomposition_1",
    irpeatmodels_required_version =
      "0.0.0"
  )

#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # degree_of_decomposition_2
#' if(! requireNamespace("posterior", quietly = TRUE)) {
#' x <-
#'   irpeat::irp_degree_of_decomposition_2(
#'     irpeat_sample_data[1, ],
#'     do_summary = TRUE,
#'     check_prediction_domain = "train",
#'     summary_function_sd = posterior::sd
#'   )
#' }
#'
#' @export
irp_degree_of_decomposition_2 <-
  irp_function_factory_eb1149(
    model = irpeatmodels::model_degree_of_decomposition_2_brms,
    config = irpeatmodels::model_degree_of_decomposition_2_config,
    prediction_domain = irpeatmodels::model_degree_of_decomposition_2_prediction_domain,
    target_variable_name =
      "degree_of_decomposition_2",
    irpeatmodels_required_version =
      "0.0.0"
  )

#' @rdname irp-predict-transmission-mir
#'
#' @examples
#' # degree_of_decomposition_3
#' if(! requireNamespace("posterior", quietly = TRUE)) {
#' x <-
#'   irpeat::irp_degree_of_decomposition_3(
#'     irpeat_sample_data[1, ],
#'     do_summary = TRUE,
#'     check_prediction_domain = "train",
#'     summary_function_sd = posterior::sd
#'   )
#' }
#'
#' @export
irp_degree_of_decomposition_3 <-
  irp_function_factory_eb1149(
    model = irpeatmodels::model_degree_of_decomposition_3_brms,
    config = irpeatmodels::model_degree_of_decomposition_3_config,
    prediction_domain = irpeatmodels::model_degree_of_decomposition_3_prediction_domain,
    target_variable_name =
      "degree_of_decomposition_3",
    irpeatmodels_required_version =
      "0.0.0"
  )



