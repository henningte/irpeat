% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/irp_content.R
\name{irp_content}
\alias{irp_content}
\title{Computes contents of components of peat from mid infrared spectra.}
\usage{
irp_content(x, variable)
}
\arguments{
\item{x}{An object of class \code{\link[ir:ir_new_ir]{ir}}.}

\item{variable}{A character vector with one or more values that define for which
components contents are computed for the spectra in \code{x}. Currently allowed values
are:
\describe{
  \item{"all"}{\code{irp_content} computes all of the values below.}
  \item{"klason_lignin_hodgkins"}{Klason lignin mass fraction [unitless]
  as computed by \code{\link{irp_content_klh_hodgkins}}.}
  \item{"holocellulose_hodgkins"}{Holocellulose mass fraction [unitless]
  as computed by \code{\link{irp_content_klh_hodgkins}}.}
}}
}
\value{
An object of class \code{\link[ir:ir_new_ir]{ir}} with additional
columns containing the computed compount contents for the spectra in
\code{x}.
}
\description{
\code{irp_content} computes the contents of various components of peat
from mid infrared spectra of peat samples. Some functions may also work
for organic matter in general.
Note that this function is a wrapper function to several individual functions
in irpeat. In \code{irp_content} any warnings and messages regarding issues
with the data are
suppressed and nothing will be esported to disk. In case of doubt, use the
corresponding functions with warnings activated.
}
\seealso{
\code{\link{irp_content_klh_hodgkins}}.
}
