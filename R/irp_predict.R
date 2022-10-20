#' Batch-prediction of sample properties
#'
#' Wrapper function to batch-predict sample properties.
#'
#' @param x An object of class [`ir`][ir::ir_new_ir]. See the individual
#' prediction models for further data requirements.
#'
#' @param variable A character vector with one or more values that define for which
#' components contents are computed for the spectra in `x`. Currently allowed values
#' are:
#' \describe{
#'   \item{"all"}{
#'     `irp_content` computes all of the values below.
#'     }
#'   \item{"klason_lignin_content_1"}{
#'     Klason lignin mass fraction \[g/g\] as computed by
#'     [irp_content_klh_hodgkins()].
#'     }
#'   \item{"holocellulose_content_1"}{
#'     Holocellulose mass fraction \[g/g\] as computed by
#'     [irp_content_klh_hodgkins()].
#'     }
#'   \item{"klason_lignin_content_2"}{
#'     Klason lignin mass fraction \[g/g\] as computed by
#'     [irp_klason_lignin_content_2()].
#'     }
#'   \item{"holocellulose_content_2"}{
#'     Holocellulose mass fraction \[g/g\] as computed by
#'     [irp_holocellulose_content_2()].
#'     }
#'   \item{"eac_1"}{
#'     Electron accepting capacity as computed by [irp_eac_1()].
#'     }
#'   \item{"edc_1"}{
#'     Electron donating capacity as computed by [irp_edc_1()].
#'     }
#'   \item{"carbon_content_1"}{
#'     Carbon content as computed by [irp_carbon_content_1()].
#'   }
#'   \item{"nitrogen_content_1"}{
#'     Nitrogen content as computed by [irp_nitrogen_content_1()].
#'   }
#'   \item{"hydrogen_content_1"}{
#'     Hydrogen content as computed by [irp_hydrogen_content_1()].
#'   }
#'   \item{"oxygen_content_1"}{
#'     Oxygen content as computed by [irp_oxygen_content_1()].
#'   }
#'   \item{"phosphorous_content_1"}{
#'     Phosphorus content as computed by [irp_phosphorous_content_1()].
#'   }
#'   \item{"potassium_content_1"}{
#'     Potassium content as computed by [irp_potassium_content_1()].
#'   }
#'   \item{"sulfur_content_1"}{
#'     Sulfur content as computed by [irp_sulfur_content_1()].
#'   }
#'   \item{"titanium_content_1"}{
#'     Titanium content as computed by [irp_titanium_content_1()].
#'   }
#'   \item{"d13C_1"}{
#'     \eqn{\delta^{13}}C values as computed by [irp_d13C_1()].
#'   }
#'   \item{"d15N_1"}{
#'     \eqn{\delta^{15}}N values as computed by [irp_d15N_1()].
#'   }
#'   \item{"nosc_1"}{
#'     The nominal oxidation state of carbon as computed by [irp_nosc_1()].
#'   }
#'   \item{"dgf0_1"}{
#'     The standard Gibbs free energy of formation content as computed by
#'     [irp_dgf0_1()].
#'   }
#'   \item{"bulk_density_1"}{
#'     Bulk density as computed by [irp_bulk_density_1()].
#'   }
#'   \item{"O_to_C_1"}{
#'     O/C ratio as computed by [irp_O_to_C_1()].
#'   }
#'   \item{"H_to_C_1"}{
#'     H/C ratio as computed by [irp_H_to_C_1()].
#'   }
#'   \item{"C_to_N_1"}{
#'     C/N ratio as computed by [irp_C_to_N_1()].
#'   }
#'   \item{"volume_fraction_solids_1"}{
#'     Volume fraction of solids as computed by
#'     [irp_volume_fraction_solids_1()].
#'   }
#'   \item{"non_macroporosity_1"}{
#'     Non-macroporosity as computed by [irp_non_macroporosity_1()].
#'   }
#'   \item{"macroporosity_1"}{
#'     Macroporosity as computed by [irp_macroporosity_1()].
#'   }
#'   \item{"saturated_hydraulic_conductivity_1"}{
#'     Saturated hydraulic conductivity as computed by
#'     [saturated_hydraulic_conductivity_1()].
#'   }
#' }
#'
#' @param ... Further arguments passed to individual prediction functions.
#'
#' @note
#' \describe{
#'   \item{`value = "klason_lignin_content_1"` and
#'   `value = "holocellulose_content_1"`}{
#'     No warnings are shown and no values are exported to disk.
#'   }
#' }
#'
#'
#' @return An object of class [`ir`][ir::ir_new_ir] with additional columns
#' containing the predictions for the spectra in `x`.
#'
#' @examples
#' library(ir)
#'
#' irp_predict(
#'   ir::ir_sample_data[1, ],
#'   variable = "carbon_content_1",
#'   do_summary = TRUE
#' )
#'
#' \dontrun{
#' irp_predict(
#'   ir::ir_sample_data[1, ],
#'   variable = c("eac_1", "carbon_content_1", "nitrogen_content_1", "dgf0_1"),
#'   do_summary = TRUE
#' )
#' }
#'
#' @export
irp_predict <- function(x, variable, ...) {

  stopifnot(inherits(x, "ir"))
  if(!is.character(variable)) {
    rlang::abort(paste0("`variable` must be a character vector, not ", class(variable)[[1]], "."))
  }
  if(length(variable) == 0) {
    rlang::abort(paste0("`variable` must contain at least one element, not ", length(variable), " elements."))
  }
  variable_values <- c("klason_lignin_content_1", "holocellulose_content_1", "holocellulose_content_2", "klason_lignin_content_2", "eac_1", "edc_1", "carbon_content_1", "nitrogen_content_1", "hydrogen_content_1", "oxygen_content_1", "phosphorous_content_1", "potassium_content_1", "sulfur_content_1", "titanium_content_1", "d13C_1", "d15N_1", "nosc_1", "dgf0_1", "bulk_density_1", "O_to_C_1", "H_to_C_1", "C_to_N_1", "volume_fraction_solids_1", "non_macroporosity_1", "macroporosity_1", "saturated_hydraulic_conductivity_1")
  if("all" %in% variable) {
    variable <- variable_values
  }
  variable_values_present <- variable_values[variable_values %in% colnames(x)]
  variable_values_nonrequested <- variable_values[! variable_values %in% variable & ! variable_values %in% variable_values_present]

  # settings to avoid duplicate computation
  has_holocellulose_content_1 <- FALSE
  has_porosity_1 <- FALSE

  for(i in seq_along(variable)) {
    x <-
      switch(variable[[i]],
             "klason_lignin_content_1" = ,
             "holocellulose_content_1" =
               if(! has_holocellulose_content_1) {
                 has_holocellulose_content_1 <<- TRUE
                 irp_content_klh_hodgkins(
                   x = x,
                   export = NULL,
                   verbose = FALSE,
                   make_plots = FALSE
                 )
               } else {
                 x
               },
             "klason_lignin_content_2" =
               irp_klason_lignin_content_2(x = x, ...),
             "holocellulose_content_2" =
               irp_holocellulose_content_2(x = x, ...),
             "eac_1" =
               irp_eac_1(x = x, ...),
             "edc_1" =
               irp_edc_1(x = x, ...),
             "carbon_content_1" =
               irp_carbon_content_1(x = x, ...),
             "nitrogen_content_1" =
               irp_nitrogen_content_1(x = x, ...),
             "hydrogen_content_1" =
               irp_hydrogen_content_1(x = x, ...),
             "oxygen_content_1" =
               irp_oxygen_content_1(x = x, ...),
             "phosphorous_content_1" =
               irp_phosphorous_content_1(x = x, ...),
             "potassium_content_1" =
               irp_potassium_content_1(x = x, ...),
             "sulfur_content_1" =
               irp_sulfur_content_1(x = x, ...),
             "titanium_content_1" =
               irp_titanium_content_1(x = x, ...),
             "d13C_1" =
               irp_d13C_1(x = x, ...),
             "d15N_1" =
               irp_d15N_1(x = x, ...),
             "nosc_1" =
               irp_nosc_1(x = x, ...),
             "dgf0_1" =
               irp_dgf0_1(x = x, ...),
             "bulk_density_1" =
               irp_bulk_density_1(x = x, ...),
             "H_to_C_1" =
               irp_H_to_C_1(x = x, ...),
             "O_to_C_1" =
               irp_O_to_C_1(x = x, ...),
             "C_to_N_1" =
               irp_C_to_N_1(x = x, ...),
             "volume_fraction_solids_1" =,
             "non_macroporosity_1" =,
             "macroporosity_1" =
               irp_porosity_1(x = x, ...),
             "saturated_hydraulic_conductivity_1" = {
               if(! has_porosity_1) {
                 has_porosity_1 <<- TRUE
                 irp_saturated_hydraulic_conductivity_1(x = x, ...)
               } else {
                 x
               }
             },
             stop(paste0("Unknown value for `variable`: ", variable[[i]]))
      )
  }

  # delete not requested variables
  x %>%
    dplyr::select(!dplyr::any_of(c(variable_values_nonrequested, paste0(variable_values_nonrequested, "_in_pd"))))

}



#' @rdname irp_predict
#' @export
irp_content <- function(x, variable, ...) {

  lifecycle::deprecate_warn("0.2.0", "irp_content()", "irpeat::irp_predict()")

  irp_predict(x, variable, ...)

}
