#' Computes contents of components of peat from mid infrared spectra.
#'
#' \code{irp_content} computes the contents of various components of peat
#' from mid infrared spectra of peat samples. Some functions may also work
#' for organic matter in general.
#' Note that this function is a wrapper function to several individual functions
#' in irpeat. In \code{irp_content} any warnings and messages regarding issues
#' with the data are
#' suppressed and nothing will be esported to disk. In case of doubt, use the
#' corresponding functions with warnings activated.
#'
#' @param x An object of class \code{\link[ir:ir_new_ir]{ir}}.
#' @param variable A character vector with one or more values that define for which
#' components contents are computed for the spectra in \code{x}. Currently allowed values
#' are:
#' \describe{
#'   \item{"all"}{\code{irp_content} computes all of the values below.}
#'   \item{"klason_lignin_hodgkins"}{Klason lignin mass fraction [unitless]
#'   as computed by \code{\link{irp_content_klh_hodgkins}}.}
#'   \item{"holocellulose_hodgkins"}{Holocellulose mass fraction [unitless]
#'   as computed by \code{\link{irp_content_klh_hodgkins}}.}
#' }
#' @return An object of class \code{\link[ir:ir_new_ir]{ir}} with additional
#' columns containing the computed compount contents for the spectra in
#' \code{x}.
#' @seealso
#' \code{\link{irp_content_klh_hodgkins}}.
#' @export
irp_content <- function(x,
                        variable) {

  ir::ir_check_ir(x)
  if(!is.character(variable)) {
    rlang::abort(paste0("`variable` must be a character vector, not ", class(variable)[[1]], "."))
  }
  if(length(variable) == 0) {
    rlang::abort(paste0("`variable` must contain at least one element, not ", length(variable), " elements."))
  }
  variable_values <- c("klason_lignin_hodgkins", "holocellulose_hodgkins")
  variable_match <- variable %in% c(variable_values, "all")
  if(!all(variable_match)) {
    if(sum(variable_match) == 1) {
      rlang::abort(paste0("`variable` must contain elements that are contained in c('klason_lignin_hodgkins', 'holocellulose_hodgkins'). Element ", which(!variable_match), " of `variable` does not match."))
    } else {
      rlang::abort(paste0("`variable` must contain elements that are contained in c('klason_lignin_hodgkins', 'holocellulose_hodgkins'). Element ", paste(which(!variable_match), collapse = ", "), " of `variable` does not match."))
    }
  }
  if("all" %in% variable) {
    variable <- variable_values
  }

  klh_hodgkins_done <- FALSE
  klh_hodgkins_variable_values <- c("klason_lignin_hodgkins", "holocellulose_hodgkins")

  for(i in seq_along(variable)) {
    switch(variable[[i]],
           "klason_lignin_hodgkins",
           "holocellulose_hodgkins" = {
             if(!klh_hodgkins_done) {
               klh_hodgkins_done <- TRUE
               x <- irp_content_klh_hodgkins(x = x,
                                             export = NULL,
                                             verbose = FALSE,
                                             make_plots = FALSE)
             }
           })
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
