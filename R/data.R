#### Models ####

#' Linear model to predict Klason lignin content from MIR spectra from \insertCite{Hodgkins.2018;textual}{ir}
#'
#'
#' `irp_content_kl_hodgkins_model` is a linear model to predict
#' Klason lignin mass fractions \[g/g\] in samples based on
#' normalized and area corrected peak heights of peaks identified in
#' spectra. This model is meant for internal use in
#' [irp_content_klh_hodgkins()]. The model was originally
#' developed by \insertCite{Hodgkins.2018;textual}{ir} and is
#' here reproduced from the same data (available via
#' [ir::ir_sample_data]), but handling
#' mass fraction in g/g instead of as percentages.
#'
#' @format An object of class [stats::lm()].
#' @source The data set was derived from <https://www.nature.com/articles/s41467-018-06050-2>
#' and published by \insertCite{Hodgkins.2018;textual}{ir} under the CC BY 4.0 license <https://creativecommons.org/licenses/by/4.0/>.
#' \insertCite{Hodgkins.2018;textual}{ir} originally derived the data on Klason lignin content from
#' \insertCite{LaCruz.2016;textual}{ir} <https://www.liebertpub.com/doi/full/10.1089/ees.2014.0402>.
#' The linear model is reproduced in slightly modified form (see description) from \insertCite{Hodgkins.2018;textual}{ir}.
#' @references
#'   \insertAllCited{}
"irp_content_kl_hodgkins_model"


#' Linear model to predict holocellulose content from MIR spectra from \insertCite{Hodgkins.2018;textual}{ir}
#'
#'
#' `irp_content_h_hodgkins_model` is a linear model to predict
#' holocellulose mass fractions \[g/g\] in samples based on
#' normalized and area corrected peak heights of peaks identified in
#' spectra. This model is meant for internal use in
#' [irp_content_klh_hodgkins()]. The model was originally
#' developed by \insertCite{Hodgkins.2018;textual}{ir} and is
#' here reproduced from the same data (available via
#' [ir::ir_sample_data]), but handling
#' mass fraction in g/g instead of as percentages.
#'
#' @format An object of class [stats::lm()].
#' @source The data set was derived from <https://www.nature.com/articles/s41467-018-06050-2>
#' and published by \insertCite{Hodgkins.2018;textual}{ir} under the CC BY 4.0 license <https://creativecommons.org/licenses/by/4.0/>.
#' \insertCite{Hodgkins.2018;textual}{ir} originally derived the data on holocellulose content from
#' \insertCite{LaCruz.2016;textual}{ir} <https://www.liebertpub.com/doi/full/10.1089/ees.2014.0402>.
#' The linear model is reproduced in slightly modified form (see description) from \insertCite{Hodgkins.2018;textual}{ir}.
#' @references
#'   \insertAllCited{}
"irp_content_h_hodgkins_model"

#' Linear model to predict holocellulose content from MIR spectra from \[---todo:ref]
#'
#'
#' `model_holocellulose_2` is a linear model to predict
#' holocellulose mass fractions \[g/g\] in samples based on mid infrared spectra.
#' Predictions with this model can be generated with
#' [irp_holocellulose_2()].
#'
#' @details The model was trained on the data from \insertCite{Hodgkins.2018;textual}{ir}
#' (available via [ir::ir_sample_data()]) and is
#' described in \[---todo:ref]. See [model_holocellulose_2_config()]
#' for details on how the training spectra were preprocessed prior model
#' fitting.
#'
#' The model is an improved version of `irp_content_h_hodgkins_model`
#' \[---todo:ref]. It is a Bayesian beta regression model using all binned
#' spectral variables for prediction.
#'
#' @note Note that this is a preliminary model only which has not been fully
#' validated for peat samples yet and which has known limitations in predicting
#' contents for peat samples [--- todo: add reference].
#'
#' @format An object of class [brms::brmsfit-class()].
#'
#' @source The model is described in \[---todo:ref]. The data set was derived
#' from <https://www.nature.com/articles/s41467-018-06050-2> and published
#' by \insertCite{Hodgkins.2018;textual}{ir} under the CC BY 4.0 license
#' <https://creativecommons.org/licenses/by/4.0/>.
#' \insertCite{Hodgkins.2018;textual}{ir} originally derived the data on
#' holocellulose content from \insertCite{LaCruz.2016;textual}{ir}
#' <https://www.liebertpub.com/doi/full/10.1089/ees.2014.0402>.
#'
#' @references
#'   \insertAllCited{}
#'
"model_holocellulose_2"

#' Linear model to predict Klason lignin content from MIR spectra from \[---todo:ref]
#'
#'
#' `model_klason_lignin_2` is a linear model to predict
#' Klason lignin mass fractions \[g/g\] in samples based on mid infrared spectra.
#' Predictions with this model can be generated with
#' [irp_klason_lignin_2()].
#'
#' @details The model was trained on the data from \insertCite{Hodgkins.2018;textual}{ir}
#' (available via [ir::ir_sample_data()]) and is
#' described in \[---todo:ref]. See [model_klason_lignin_2_config()]
#' for details on how the training spectra were preprocessed prior model
#' fitting.
#'
#' The model is an improved version of `irp_content_kl_hodgkins_model`
#' \[---todo:ref]. It is a Bayesian beta regression model using all binned
#' spectral variables for prediction.
#'
#' @note Note that this is a preliminary model only which has not been fully
#' validated for peat samples yet and which has known limitations in predicting
#' contents for peat samples [--- todo: add reference].
#'
#' @format An object of class [brms::brmsfit-class()].
#'
#' @source The model is described in \[---todo:ref]. The data set was derived
#' from <https://www.nature.com/articles/s41467-018-06050-2> and published
#' by \insertCite{Hodgkins.2018;textual}{ir} under the CC BY 4.0 license
#' <https://creativecommons.org/licenses/by/4.0/>.
#' \insertCite{Hodgkins.2018;textual}{ir} originally derived the data on
#' Klason lignin content from \insertCite{LaCruz.2016;textual}{ir}
#' <https://www.liebertpub.com/doi/full/10.1089/ees.2014.0402>.
#'
#' @references
#'   \insertAllCited{}
#'
"model_klason_lignin_2"

#' Linear model to predict peat electron accepting capacities from mid infrared spectra from \insertCite{Teickner.2022;textual}{irpeat}
#'
#'
#' `model_eac_1` is a linear model to predict
#' electron accepting capacities \[\eqn{\mu}mol g\eqn{_\text{C}^{-1}}\] in peat
#' samples based on mid infrared spectra. Predictions with this model can be
#' generated with [irp_eac_1()].
#'
#' @details See \insertCite{Teickner.2022;textual}{irpeat} for a detailed description of the model.
#'
#' @note Note that this model still has several limitations described in
#' \insertCite{Teickner.2022;textual}{irpeat}.
#'
#' @format An object of class [`stanreg`][rstanarm::stanreg-objects].
#' @source The model is described in \insertCite{Teickner.2022;textual}{irpeat}.
#' @seealso [irp_eac_1()], [model_configuration()].
#' @references
#'   \insertAllCited{}
"model_eac_1"

#' Linear model to predict peat electron donating capacities from mid infrared spectra from \insertCite{Teickner.2022;textual}{irpeat}
#'
#'
#' `model_edc_1` is a linear model to predict
#' electron donating capacities \[\eqn{\mu}mol g\eqn{_\text{C}^{-1}}\] in peat
#' samples based on mid infrared spectra. Predictions with this model can be
#' generated with [irp_edc_1()].
#'
#' @details See \insertCite{Teickner.2022;textual}{irpeat} for a detailed description of the model.
#'
#' @note Note that this model is known to make biased predictions and has
#' several other limitations described in \insertCite{Teickner.2022;textual}{irpeat}.
#'
#' @format An object of class [`stanreg`][rstanarm::stanreg-objects].
#' @source The model is described in \insertCite{Teickner.2022;textual}{irpeat}.
#' @seealso [irp_edc_1()], [model_configuration()].
#' @references
#'   \insertAllCited{}
"model_edc_1"


#### Configuration ####

#' Configuration files defining preprocessing parameters.
#'
#' The different configuration files are used in the models with the respective
#' name. See [irp_preprocess()].
#'
#' @name model_configuration
#' @format A list with the following elements:
#' \describe{
#'    \item{`irp_preprocess`}{
#'      A list with the a subset of the following elements:
#'      \describe{
#'         \item{`do_interpolate`}{See [irp_preprocess()].}
#'         \item{`interpolate_start`}{See [irp_preprocess()].}
#'         \item{`interpolate_dw`}{See [irp_preprocess()].}
#'         \item{`do_clip`}{See [irp_preprocess()].}
#'         \item{`clip_range`}{See [irp_preprocess()].}
#'         \item{`do_interpolate_region`}{See [irp_preprocess()].}
#'         \item{`interpolate_region_range`}{See [irp_preprocess()].}
#'         \item{`do_bc`}{See [irp_preprocess()].}
#'         \item{`bc_method`}{See [irp_preprocess()].}
#'         \item{`bc_degree`}{See [irp_preprocess()].}
#'         \item{`bc_cutoff`}{See [irp_preprocess()].}
#'         \item{`do_smooth`}{See [irp_preprocess()].}
#'         \item{`smooth_method`}{See [irp_preprocess()].}
#'         \item{`smooth_p`}{See [irp_preprocess()].}
#'         \item{`smooth_n`}{See [irp_preprocess()].}
#'         \item{`smooth_m`}{See [irp_preprocess()].}
#'         \item{`smooth_ts`}{See [irp_preprocess()].}
#'         \item{`smooth_k`}{See [irp_preprocess()].}
#'         \item{`do_normalise`}{See [irp_preprocess()].}
#'         \item{`normalise_method`}{See [irp_preprocess()].}
#'         \item{`do_bin`}{See [irp_preprocess()].}
#'         \item{`bin_width`}{See [irp_preprocess()].}
#'         \item{`do_scale`}{See [irp_preprocess()].}
#'         \item{`scale_center`}{See [irp_preprocess()]. This is
#'         used to preprocess the training data before model computation. The
#'         values in `data_scale` are used to scale data used for
#'         prediction.}
#'         \item{`scale_scale`}{See [irp_preprocess()]. This is
#'         used to preprocess the training data before model computation. The
#'         values in `data_scale` are used to scale data used for
#'         prediction.}
#'      }
#'    }
#'    \item{`data_scale`}{
#'      Values used for scaling data used for prediction. A list with the
#'      following elements:
#'      \describe{
#'         \item{`y_center`}{A numeric value representing the value to
#'         subtract from the dependent variable of a model during
#'         scaling. Corresponds to argument `center` in
#'         [base::scale()].}
#'         \item{`y_scale`}{A numeric value representing the value by
#'         which the dependent variable of a model is divided during
#'         scaling. Corresponds to argument `scale` in
#'         [base::scale()].}
#'         \item{`x_center`}{The same as `y_center`, but for the
#'         independent variable(s). Can be a numeric vector with an element for
#'         each independent variable.}
#'         \item{`y_scale`}{The same as `y_center`, but for the
#'         independent variable(s). Can be a numeric vector with an element for
#'         each independent variable.}
#'      }
#'    }
#' }

#' @rdname model_configuration
"model_klason_lignin_2_config"

#' @rdname model_configuration
"model_holocellulose_2_config"

#' @rdname model_configuration
"model_eac_1_config"

#' @rdname model_configuration
"model_edc_1_config"

