#### Models ####

#' Linear model to predict Klason lignin content from MIR spectra from \insertCite{Hodgkins.2018;textual}{ir}.
#'
#'
#' \code{irp_content_kl_hodgkins_model} is a linear model to predict
#' Klason lignin mass fractions [g/g] in samples based on
#' normalised and area corrected peak heights of peaks identified in
#' spectra. This model is meant for internal use in
#' \code{\link{irp_content_klh_hodgkins}}. The model was originally
#' developed by \insertCite{Hodgkins.2018;textual}{ir} and is
#' here reproduced from the same data (available via
#' \code{\link[ir:ir_sample_data]{ir::ir_sample_data}}), but handling
#' mass fraction in g/g instead of as percentages.
#'
#' @format An object of class \code{\link[stats:lm]{lm}}.
#' @source The data set was derived from \url{https://www.nature.com/articles/s41467-018-06050-2}
#' and published by \insertCite{Hodgkins.2018;textual}{ir} under the CC BY 4.0 license \url{https://creativecommons.org/licenses/by/4.0/}.
#' \insertCite{Hodgkins.2018;textual}{ir} originally derived the data on Klason lignin content from
#' \insertCite{LaCruz.2016;textual}{ir} \url{https://www.liebertpub.com/doi/full/10.1089/ees.2014.0402}.
#' The linear model is reproduced in slightly modified form (see description) from \insertCite{Hodgkins.2018;textual}{ir}.
#' @references
#'   \insertAllCited{}
"irp_content_kl_hodgkins_model"


#' Linear model to predict holocellulose content from MIR spectra from \insertCite{Hodgkins.2018;textual}{ir}.
#'
#'
#' \code{irp_content_h_hodgkins_model} is a linear model to predict
#' holocellulose mass fractions [g/g] in samples based on
#' normalised and area corrected peak heights of peaks identified in
#' spectra. This model is meant for internal use in
#' \code{\link{irp_content_klh_hodgkins}}. The model was originally
#' developed by \insertCite{Hodgkins.2018;textual}{ir} and is
#' here reproduced from the same data (available via
#' \code{\link[ir:ir_sample_data]{ir::ir_sample_data}}), but handling
#' mass fraction in g/g instead of as percentages.
#'
#' @format An object of class \code{\link[stats:lm]{lm}}.
#' @source The data set was derived from \url{https://www.nature.com/articles/s41467-018-06050-2}
#' and published by \insertCite{Hodgkins.2018;textual}{ir} under the CC BY 4.0 license \url{https://creativecommons.org/licenses/by/4.0/}.
#' \insertCite{Hodgkins.2018;textual}{ir} originally derived the data on holocellulose content from
#' \insertCite{LaCruz.2016;textual}{ir} \url{https://www.liebertpub.com/doi/full/10.1089/ees.2014.0402}.
#' The linear model is reproduced in slightly modified form (see description) from \insertCite{Hodgkins.2018;textual}{ir}.
#' @references
#'   \insertAllCited{}
"irp_content_h_hodgkins_model"

#' Linear model to predict peat electron accepting capacities from mid infrared spectra from \insertCite{Teickner.submitted;textual}{irpeat}.
#'
#'
#' \code{model_eac_1} is a linear model to predict
#' electron accepting capacities [\eqn{\mu}mol g\eqn{_\text{C}^{-1}}] in peat
#' samples based on mid infrared spectra. Predictions with this model can be
#' generated with \code{\link{irp_eac_1}}.
#'
#' @details See \insertCite{Teickner.submitted;textual}{irpeat} for a detailed description of the model.
#'
#' @note Note that this model still has several limitations described in
#' \insertCite{Teickner.submitted;textual}{irpeat}.
#'
#' @format An object of class \code{\link[rstanarm:stanreg-objects]{stanreg}}.
#' @source The model is described in \insertCite{Teickner.submitted;textual}{irpeat}.
#' @seealso \code{\link{irp_eac_1}}, \code{\link{model_configuration}}.
#' @references
#'   \insertAllCited{}
"model_eac_1"

#' Linear model to predict peat electron donating capacities from mid infrared spectra from \insertCite{Teickner.submitted;textual}{irpeat}.
#'
#'
#' \code{model_edc_1} is a linear model to predict
#' electron donating capacities [\eqn{\mu}mol g\eqn{_\text{C}^{-1}}] in peat
#' samples based on mid infrared spectra. Predictions with this model can be
#' generated with \code{\link{irp_edc_1}}.
#'
#' @details See \insertCite{Teickner.submitted;textual}{irpeat} for a detailed description of the model.
#'
#' @note Note that this model is known to make biased predictions and has
#' several other limitations described in \insertCite{Teickner.submitted;textual}{irpeat}.
#'
#' @format An object of class \code{\link[rstanarm:stanreg-objects]{stanreg}}.
#' @source The model is described in \insertCite{Teickner.submitted;textual}{irpeat}.
#' @seealso \code{\link{irp_edc_1}}, \code{\link{model_configuration}}.
#' @references
#'   \insertAllCited{}
"model_edc_1"


#### Configuration ####

#' Configuration files defining preprocessing parameters.
#'
#' The different configuration files are used in the models with the respective
#' name. See \code{\link{irp_preprocess}}.
#'
#' @name model_configuration
#' @format A list with the following elements:
#' \describe{
#'    \item{\code{irp_preprocess}}{
#'      A list with the a subset of the following elements:
#'      \describe{
#'         \item{\code{do_interpolate}}{See \code{\link{irp_preprocess}}.}
#'         \item{\code{interpolate_start}}{See \code{\link{irp_preprocess}}.}
#'         \item{\code{interpolate_dw}}{See \code{\link{irp_preprocess}}.}
#'         \item{\code{do_clip}}{See \code{\link{irp_preprocess}}.}
#'         \item{\code{clip_range}}{See \code{\link{irp_preprocess}}.}
#'         \item{\code{do_interpolate_region}}{See \code{\link{irp_preprocess}}.}
#'         \item{\code{interpolate_region_range}}{See \code{\link{irp_preprocess}}.}
#'         \item{\code{do_bc}}{See \code{\link{irp_preprocess}}.}
#'         \item{\code{bc_method}}{See \code{\link{irp_preprocess}}.}
#'         \item{\code{bc_degree}}{See \code{\link{irp_preprocess}}.}
#'         \item{\code{bc_cutoff}}{See \code{\link{irp_preprocess}}.}
#'         \item{\code{do_smooth}}{See \code{\link{irp_preprocess}}.}
#'         \item{\code{smooth_method}}{See \code{\link{irp_preprocess}}.}
#'         \item{\code{smooth_p}}{See \code{\link{irp_preprocess}}.}
#'         \item{\code{smooth_n}}{See \code{\link{irp_preprocess}}.}
#'         \item{\code{smooth_m}}{See \code{\link{irp_preprocess}}.}
#'         \item{\code{smooth_ts}}{See \code{\link{irp_preprocess}}.}
#'         \item{\code{smooth_k}}{See \code{\link{irp_preprocess}}.}
#'         \item{\code{do_normalise}}{See \code{\link{irp_preprocess}}.}
#'         \item{\code{normalise_method}}{See \code{\link{irp_preprocess}}.}
#'         \item{\code{do_bin}}{See \code{\link{irp_preprocess}}.}
#'         \item{\code{bin_width}}{See \code{\link{irp_preprocess}}.}
#'         \item{\code{do_scale}}{See \code{\link{irp_preprocess}}.}
#'         \item{\code{scale_center}}{See \code{\link{irp_preprocess}}. This is
#'         used to preprocess the training data before model computation. The
#'         values in \code{data_scale} are used to scale data used for
#'         prediction.}
#'         \item{\code{scale_scale}}{See \code{\link{irp_preprocess}}. This is
#'         used to preprocess the training data before model computation. The
#'         values in \code{data_scale} are used to scale data used for
#'         prediction.}
#'      }
#'    }
#'    \item{\code{data_scale}}{
#'      Values used for scaling data used for prediction. A list with the
#'      following elements:
#'      \describe{
#'         \item{\code{y_center}}{A numeric value representing the value to
#'         subtract from the dependent variable of a model during
#'         scaling. Corresponds to argument \code{center} in
#'         \code{\link[base]{scale}}.}
#'         \item{\code{y_scale}}{A numeric value representing the value by
#'         which the dependent variable of a model is divided during
#'         scaling. Corresponds to argument \code{scale} in
#'         \code{\link[base]{scale}}.}
#'         \item{\code{x_center}}{The same as \code{y_center}, but for the
#'         independent variable(s). Can be a numeric vector with an element for
#'         each independent variable.}
#'         \item{\code{y_scale}}{The same as \code{y_center}, but for the
#'         independent variable(s). Can be a numeric vector with an element for
#'         each independent variable.}
#'      }
#'    }
#' }

#' @rdname model_configuration
"model_eac_1_config"

#' @rdname model_configuration
"model_edc_1_config"
