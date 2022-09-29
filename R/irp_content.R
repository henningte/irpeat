#' Computes contents of components of peat from mid infrared spectra
#'
#' `irp_content` computes the contents of various components of peat
#' from mid infrared spectra of peat samples. Some functions may also work
#' for organic matter in general.
#' Note that this function is a wrapper function to several individual functions
#' in irpeat. In `irp_content` any warnings and messages regarding issues
#' with the data are
#' suppressed and nothing will be exported to disk. In case of doubt, use the
#' corresponding functions with warnings activated.
#'
#' @param x An object of class [`ir`][ir::ir_new_ir].
#'
#' @param variable A character vector with one or more values that define for which
#' components contents are computed for the spectra in `x`. Currently allowed values
#' are:
#' \describe{
#'   \item{"all"}{`irp_content` computes all of the values below.}
#'   \item{"klason_lignin_content_1"}{Klason lignin mass fraction \[g/g\]
#'     as computed by [irp_content_klh_hodgkins()].}
#'   \item{"holocellulose_content_1"}{Holocellulose mass fraction \[g/g\]
#'     as computed by [irp_content_klh_hodgkins()].}
#'   \item{"klason_lignin_content_2"}{Klason lignin mass fraction \[g/g\]
#'     as computed by [irp_klason_lignin_content_2()].}
#'   \item{"holocellulose_content_2"}{Holocellulose mass fraction \[g/g\]
#'     as computed by [irp_holocellulose_content_2()].}
#'   \item{"eac_1"}{Electron accepting capacity as computed by [irp_eac_1()].}
#'   \item{"edc_1"}{Electron donating capacity as computed by [irp_edc_1()].}
#' }
#'
#' @param ... Further arguments passed to individual prediction functions.
#'
#' @return An object of class [`ir`][ir::ir_new_ir] with additional columns
#' containing the predictions for the spectra in `x`.
#'
#' @examples
#' library(ir)
#'
#' irp_content(ir::ir_sample_data[1, ], variable = "klason_lignin_content_1")
#'
#' \dontrun{
#' irp_content(ir::ir_sample_data[1, ], variable = "klason_lignin_content_2")
#' irp_content(ir::ir_sample_data[1, ], variable = "holocellulose_content_2")
#' irp_content(ir::ir_sample_data[1, ], variable = "eac_1")
#' irp_content(ir::ir_sample_data[1, ], variable = "edc_1")
#' }
#'
#' @export
irp_content <- function(x, variable, ...) {

  stopifnot(inherits(x, "ir"))
  if(!is.character(variable)) {
    rlang::abort(paste0("`variable` must be a character vector, not ", class(variable)[[1]], "."))
  }
  if(length(variable) == 0) {
    rlang::abort(paste0("`variable` must contain at least one element, not ", length(variable), " elements."))
  }
  variable_values <- c("klason_lignin_content_1", "holocellulose_content_1", "holocellulose_content_2", "klason_lignin_content_2", "eac_1", "edc_1")
  variable_match <- variable %in% c(variable_values, "all")
  if(!all(variable_match)) {
    if(sum(variable_match) == 1) {
      rlang::abort(paste0("`variable` must contain elements that are contained in ", paste(variable_values, collapse = ", "), ". Element ", which(!variable_match), " of `variable` does not match."))
    } else {
      rlang::abort(paste0("`variable` must contain elements that are contained in ", paste(variable_values, collapse = ", "), ". Element ", paste(which(!variable_match), collapse = ", "), " of `variable` does not match."))
    }
  }
  if("all" %in% variable) {
    variable <- variable_values
  }

  klh_hodgkins_done <- FALSE
  klh_hodgkins_variable_values <- c("klason_lignin_content_1", "holocellulose_content_1")

  for(i in seq_along(variable)) {
    switch(variable[[i]],
           "klason_lignin_content_1" = ,
           "holocellulose_content_1" = {
             if(!klh_hodgkins_done) {
               klh_hodgkins_done <- TRUE
               x <- irp_content_klh_hodgkins(x = x,
                                             export = NULL,
                                             verbose = FALSE,
                                             make_plots = FALSE)
             }
           },
           "klason_lignin_content_2" = {
             x <- irp_klason_lignin_content_2(x = x, ...)
           },
           "holocellulose_content_2" = {
             x <- irp_holocellulose_content_2(x = x, ...)
           },
           "eac_1" = {
             x <- irp_eac_1(x = x, ...)
           },
           "edc_1" = {
             x <- irp_edc_1(x = x, ...)
           },
           stop("Unknown value for `variable`.")
          )
  }

  # delete not demanded variables
  if(sum(klh_hodgkins_variable_values %in% variable) == 1) {
    if(which(klh_hodgkins_variable_values %in% variable) == 1) {
      x <- x[, -match(klh_hodgkins_variable_values[[2]], colnames(x))]
    } else {
      x <- x[, -match(klh_hodgkins_variable_values[[1]], colnames(x))]
    }
  }

  x

}
