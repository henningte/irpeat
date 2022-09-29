#' Computes contents of Klason lignin and holocellulose of peat from mid infrared spectra
#'
#' `irp_content_klh_hodgkins` computes the mass fraction of Klason lignin and
#' holocellulose in peat from mid infrared spectra of the peat samples. This
#' function may also work for organic matter in general. Note that the models
#' have been shown to be biased for peat samples
#' \insertCite{Teickner.2022a}{irpeat}.
#'
#' @param x An object of class [`ir`][ir::ir_new_ir].
#'
#' @param export Either a valid path to an existing directory where to
#' store the results of the original script (see section "Source")
#' `irp_content_klh_hodgkins` is based on (if this output should be exported) or
#' `NULL` (if nothing should be exported). Note that `export` has to be
#' specified in order to analyze any issues with the data (e.g. high clay content)
#' that may confound the computed values. Refer to the original publication
#' (\insertCite{Hodgkins.2018}{ir}) to get information on the meaning of
#' the exported objects.
#'
#' @param verbose A logical value indicating if messages should be printed
#' (`verbose = TRUE`) or not (`verbose = FALSE`).
#'
#' @param make_plots logical value indicating if plots should be printed
#' (`make_plots = TRUE`) or not (`make_plots = FALSE`).
#'
#' @note Note that the models have been shown to be biased for peat samples
#' \insertCite{Teickner.2022a}{irpeat}.
#'
#' @source `irp_content_klh_hodgkins` is a wrapper function to a script written by
#' Suzanne Hodgkins (<https://github.com/shodgkins/FTIRbaselines/blob/master/FTIRbaselines.R>).
#' and distributed under the GPL-3 (<https://www.gnu.org/licenses/gpl-3.0.html>) license.
#' The original script was modified for this purpose.
#'
#' @return An object of class [`ir`][ir::ir_new_ir] with additional
#' columns containing the computed mass fractions of Klason lignin and
#' holocellulose for the spectra in `x`.
#'
#' @examples
#' library(ir)
#'
#' irp_content_klh_hodgkins(
#'   ir::ir_sample_data[1:5, ],
#'   export = NULL,
#'   verbose = TRUE,
#'   make_plots = FALSE
#' )
#'
#' @references
#'   \insertAllCited{}
#'
#' @export
irp_content_klh_hodgkins <- function(x,
                                     export = NULL,
                                     verbose = FALSE,
                                     make_plots = FALSE) {

  check_irpeatmodels(version = "0.0.0")

  stopifnot(inherits(x, "ir"))
  if(!is.null(export)) {
    if(!is.character(export)) {
      rlang::abort(paste0("`export` must be a character value, not ", class(export)[[1]], "."))
    }
    if(length(export) != 1) {
      rlang::abort(paste0("`export` must be of length 1, not ", length(export), "."))
    }
    if(!dir.exists(export)) {
      rlang::abort("The directory specified by `export` does not exist.")
    }
  }
  if(length(verbose) != 1) {
    rlang::abort(paste0("`verbose` must be of length 1, not ", length(verbose), "."))
  }
  if(!is.logical(verbose)) {
    rlang::abort(paste0("`verbose` must be a character value, not ", class(verbose)[[1]], "."))
  }
  if(length(make_plots) != 1) {
    rlang::abort(paste0("`make_plots` must be of length 1, not ", length(make_plots), "."))
  }
  if(!is.logical(verbose)) {
    rlang::abort(paste0("`make_plots` must be a character value, not ", class(make_plots)[[1]], "."))
  }

  x_flat <- ir::ir_flatten(x)
  x_flat_empty <- ir::ir_flat_clean(x = x_flat, return_empty = TRUE)
  x_flat <- ir::ir_flat_clean(x = x_flat, return_empty = FALSE)

  x_flat_m <- irp_content_klh_hodgkins_prepare(x)

  res <- irp_content_klh_hodgkins_main(data = x_flat_m,
                                       export = export,
                                       verbose = verbose,
                                       make_plots = make_plots)

  if(verbose) {
    message("Computing Klason lignin and holocellulose contents.")
  }

  irp_content_klh_hodgkins_predict(x = x,
                                   x_flat = x_flat,
                                   x_flat_empty = x_flat_empty,
                                   res = res)

}

# Function irp_content_klh_hodgkins_predict
irp_content_klh_hodgkins_predict <- function(x,
                                             x_flat,
                                             x_flat_empty,
                                             res) {

  # utils::data("irp_content_h_hodgkins_model", envir=environment())
  # utils::data("irp_content_kl_hodgkins_model", envir=environment())

  # make sure that all values get their right row in x
  x_flat_new_names <- c(colnames(x_flat)[-1], colnames(x_flat_empty)[-1])
  x_flat_target_names <- c(as.character(seq_len(nrow(x))))
  index <- match(x_flat_target_names, x_flat_new_names)
  x_flat_empty_n <- ncol(x_flat_empty) - 1
  res_fake <- as.data.frame(matrix(nrow = x_flat_empty_n, ncol = ncol(res$norm.Acorr)))
  colnames(res_fake) <- colnames(res$norm.Acorr)
  res$norm.Acorr <- dplyr::bind_rows(res$norm.Acorr, res_fake)[index,]

  # holocellulose
  newdata_h_hodgkins <-
    data.frame(
      lm_x = res$norm.Acorr$carb,
      stringsAsFactors = FALSE
    )
  prediction_h_hodgkins <-
    as.data.frame(
      stats::predict(
        irpeatmodels::model_holocellulose_content_1,
        newdata = newdata_h_hodgkins,
        se.fit = TRUE)
    )
  prediction_h_hodgkins$se_pi <- sqrt(prediction_h_hodgkins$se.fit^2 + prediction_h_hodgkins$residual.scale^2)
  x$holocellulose_content_1 <-
    quantities::set_quantities(
      prediction_h_hodgkins$fit,
      unit = "g/g",
      errors = prediction_h_hodgkins$se_pi
    )

  # Klason lignin
  newdata_kl_hodgkins <-
    data.frame(
      lm_x = res$norm.Acorr$arom15 + res$norm.Acorr$arom16,
                                    stringsAsFactors = FALSE
      )
  prediction_kl_hodgkins <-
    as.data.frame(
      stats::predict(irpeatmodels::model_klason_lignin_content_1,
                     newdata = newdata_kl_hodgkins,
                     se.fit = TRUE
      )
    )
  prediction_kl_hodgkins$se_pi <- sqrt(prediction_kl_hodgkins$se.fit^2 + prediction_kl_hodgkins$residual.scale^2)
  x$klason_lignin_content <-
    quantities::set_quantities(
      prediction_kl_hodgkins$fit, unit = "g/g",
      errors = prediction_kl_hodgkins$se_pi
    )
  x

}

# Function to prepare objects of class ir for the main function
irp_content_klh_hodgkins_prepare <- function(x) {

  x_flat <- ir::ir_interpolate(x = x, start = NULL, dw = 1)
  x_flat <- ir::ir_flatten(x = x_flat, measurement_id = as.character(seq_len(nrow(x_flat))))
  x_flat <- ir::ir_flat_clean(x_flat, return_empty = FALSE)
  x_flat <- as.data.frame(x_flat[order(x_flat$x, decreasing = TRUE), ])
  rownames(x_flat) <- x_flat$x
  x_flat$x <- NULL
  x_flat

}

# Main function
irp_content_klh_hodgkins_main <- function(data,
                                          export = NULL,
                                          verbose = TRUE,
                                          make_plots = TRUE) {

  # The code in this function is a modified version of the script FTIRBaselines.R written by Suzanne
  # Hodgkins (Copyright 2017 Suzanne Hodgkins and Florida State University) and distributed under the GPL-3 license at https://github.com/shodgkins/FTIRbaselines.
  # Modifications include renaming functions, uncommenting interactive parts and other smaller modifications.
  # The modification were done on 2020-01-26.

  # Function for finding peak boundaries
  irp_content_klh_hodgkins_minima <- function(region,
                                              default = NULL) {

    # region = segment of xydata containing the peak boundary (see definition of xydata in next section)
    # default = wavenumber to use if a suitable minimum is not found

    success <- TRUE # changes to FALSE if the minimum resorts to default

    # First try to find a true local minimum, based on were dy crosses the x-axis and d2y > 0:
    offset.dy <- c(0, region$dy)[1:(length(region$dy))]
    just.crossed <- which(sign(region$dy*offset.dy) == -1 & region$d2y > 0)
    crossing <- c(just.crossed, just.crossed-1)
    if(length(crossing) > 0) {
      # Find the absolute minimum among the local minima candidates.
      index.min <- which(region$y==min(region$y[crossing]) & region$x %in% region$x[crossing]) # --- todo: to avoid cases, where there are multiple values with the same minimum intensity, let it select just the first value by default, by adding a "[[1]]" at the end.
      # Location of absolute minimum in region of index.min, for cases where multi-point derivative causes slight error.
      # Deactivate if this causes problems.
      index.min <- which(region$y==min(region$y[(index.min-min(m, index.min-1)):(index.min+min(m, nrow(region)-index.min))]) & region$x %in% region$x[(index.min-min(m, index.min-1)):(index.min+min(m, nrow(region)-index.min))])
      # wavenumber of local minimum
      minimum <- region$x[index.min]
    } else {
      # If the above doesn't work, use the maximum of the 2nd derivative if it is positive (i.e. "shoulder minimum"):
      if(max(region$d2y) > 0) {
        minimum <- region$x[which(region$d2y==max(region$d2y))]
      } else {
        # if there are no minima or "shoulder minima", use default
        minimum <- default
        success <- FALSE
      }
    }
    if(length(minimum)>1) {
      if(verbose) {
        message(paste(" Warning: multiple minima for", default, "peak boundary:", paste(minimum, collapse=" "), "\n"))
      }
    }
    return(list(minimum=minimum, success=success))

  }

  # Function for finding peaks, for use with acids peak
  irp_content_klh_hodgkins_maxima <- function(region,
                                              default = NULL) {

    # region = segment of xydata containing the peak boundary (see definition of xydata in next section)
    # default = wavenumber to use if a suitable minimum is not found

    success <- TRUE # changes to FALSE if the minimum resorts to default

    # First try to find a true local maximum, based on were dy crosses the x-axis and d2y < 0:
    offset.dy <- c(0, region$dy)[1:(length(region$dy))]
    just.crossed <- which(sign(region$dy*offset.dy) == -1 & region$d2y < 0)
    crossing <- c(just.crossed, just.crossed-1)
    if(length(crossing) > 0) {
      # Find the absolute maximum among the local maxima candidates.
      index.max <- which(region$y==max(region$y[crossing]) & region$x %in% region$x[crossing])
      # Location of absolute maximum in region of index.max, for cases where multi-point derivative causes slight error.
      # Deactivate if this causes problems.
      index.max <- which(region$y==max(region$y[(index.max-min(m, index.max-1)):(index.max+min(m, nrow(region)-index.max))]) & region$x %in% region$x[(index.max-min(m, index.max-1)):(index.max+min(m, nrow(region)-index.max))])
      # wavenumber of local maximum
      maximum <- region$x[index.max]
    } else {
      # If the above doesn't work, use the minimum of the 2nd derivative if it is negative (i.e. "shoulder maximum"):
      if(min(region$d2y) < 0) {
        maximum <- region$x[which(region$d2y==min(region$d2y))]
      } else {
        # if there are no maxima or "shoulder maxima", use default
        maximum <- default
        success <- FALSE
      }
    }
    if(length(maximum)>1) {
      if(verbose) {
        message(paste(" Warning: multiple maxima for", default, "peak:", paste(maximum, collapse=" "), "\n"))
      }
    }
    return(list(maximum=maximum, success=success))
  }

  # Startup notice ----
  # cat("FTIRbaselines: Copyright \u00a9 2017 Suzanne Hodgkins and Florida State University.
  #     This program comes with ABSOLUTELY NO WARRANTY. This is free software, and you
  #     are welcome to redistribute it under the terms of the GNU General Public
  #     License (version 3 or later). This license is included with this program and
  #     is also available at <http://www.gnu.org/licenses/>.")

  # Import data ----

  # filename <- readline("Enter name of csv with FTIR absorbance data: ")
  # data <- read.csv(filename, header=TRUE, row.names=1)

  # Check input data ----

  # if((as.numeric(row.names(data)[1]) - as.numeric(row.names(data)[2])) != 1) {
  #   stop("Input CSV must have wavenumbers in the first column in descending order.")
  # }

  # check for absorbance mode
  # if(max(data) > 5 | is.na(max(data))) {
  #   stop("Input data must be in Absorbance mode, with no extra rows except for the header row.")
  # }

  # Prepare to analyze data ----

  n <- 29  # number of points to differentiate over; must be odd
  m <- (n-1)/2

  peaks <- c("carb", "arom15", "arom16", "trough16", "acids", "aliph28", "trough28", "aliph29")

  # preallocate dataframes
  Wp <- data.frame(matrix(nrow=length(data),ncol=length(peaks),dimnames=list(names(data),peaks)))
  W1 <- Wp
  W2 <- Wp
  Ap <- Wp
  A1 <- Wp
  A2 <- Wp
  Ab <- Wp
  Acorr <- Wp
  success.W1 <- Wp
  success.W2 <- Wp
  area.wholepeak <- Wp
  area.corrpeak <- Wp
  notes <- data.frame(aliph.type=rep(NA,length(data)), acids.type=rep(NA,length(data)), notes=rep(NA,length(data)), row.names=names(data))

  # Find peaks and peak boundaries for each spectrum ----

  x <- as.numeric(row.names(data))

  for(i in seq_along(data)) {
    if(verbose) {
      message(paste0("Processing ", names(data)[i], " (", i, " of ", length(data), ")\n"))
    }

    y <- data[,i]
    xydata <- data.frame(x=x, y=y)

    # calculate first derivative
    dy <- rep(NA, length(y))
    for(j in (m+1):(length(dy)-m)) {
      dy[j] <- stats::coefficients(stats::lsfit(x[(j-m):(j+m)], y[(j-m):(j+m)]))[2]
    }
    xydata <- cbind(xydata, dy)

    # calculate second derivative
    d2y <- rep(NA, length(y))
    for(j in (2*m+1):(length(d2y)-2*m)) {
      d2y[j] <- stats::coefficients(stats::lsfit(x[(j-m):(j+m)], dy[(j-m):(j+m)]))[2]
    }
    xydata <- cbind(xydata, d2y)

    # Find w1, w2, and location of each peak (can alter ranges as needed).

    # 1030 peak ("carb") ----

    w1.result <- irp_content_klh_hodgkins_minima(xydata[x %in% 890:920,], default=905)
    w1 <- max(w1.result[[1]]) # the max accounts for the unlikely event where 2 minima were found with exactly the same y
    a1 <- y[x == w1]
    W1[i,"carb"] <- w1
    A1[i,"carb"] <- a1
    success.W1[i,"carb"] <- w1.result[[2]]

    #w2.result <- irp_content_klh_hodgkins_minima(xydata[x %in% 1100:1150,], default=1135) # lower w2 option
    w2.result <- irp_content_klh_hodgkins_minima(xydata[x %in% 1150:1210,], default=1185) # higher w2 option
    w2 <- min(w2.result[[1]]) # the min accounts for the unlikely event where 2 minima were found with exactly the same y
    a2 <- y[x == w2]
    W2[i,"carb"] <- w2
    A2[i,"carb"] <- a2
    success.W2[i,"carb"] <- w2.result[[2]]

    # baseline correction calculated for all wavenumbers in peak
    w.wholepeak <- w2:w1
    a.wholepeak <- y[x %in% w.wholepeak]
    a.corrpeak <- a.wholepeak-((a2-a1)*(w.wholepeak-w1)/(w2-w1)+a1)

    # peak areas
    area.wholepeak[i,"carb"] <- sum(a.wholepeak)
    area.corrpeak[i,"carb"] <- sum(a.corrpeak)

    # pick out peak location
    Acorr[i,"carb"] <- max(a.corrpeak)
    Wp[i,"carb"] <- min(w.wholepeak[a.corrpeak==Acorr[i,"carb"]]) # the min accounts for unlikely event of 2 peaks with exactly the same y; the same is true for all other peaks defined as min(w.wholepeak...) etc.
    Ap[i,"carb"] <- a.wholepeak[w.wholepeak==Wp[i,"carb"]]
    Ab[i,"carb"] <- Ap[i,"carb"] - Acorr[i,"carb"]

    rm(w1.result, w2.result, w1, w2, a1, a2, w.wholepeak, a.wholepeak, a.corrpeak)


    # 1510 peak ("arom15") *** NOTE THIS PEAK IS HIGHLY VARIABLE *** ----

    w1.result <- irp_content_klh_hodgkins_minima(xydata[x %in% 1470:1500,], default=1485)
    w1 <- max(w1.result[[1]])
    a1 <- y[x == w1]
    W1[i,"arom15"] <- w1
    A1[i,"arom15"] <- a1
    success.W1[i,"arom15"] <- w1.result[[2]]

    w2.result <- irp_content_klh_hodgkins_minima(xydata[x %in% 1515:1540,], default=1535)
    w2 <- min(w2.result[[1]])
    a2 <- y[x == w2]
    W2[i,"arom15"] <- w2
    A2[i,"arom15"] <- a2
    success.W2[i,"arom15"] <- w2.result[[2]]

    #  # redefine either w1 or w2 so that the baseline doesn't cut through part of the graph
    #  if((a2-a1)/(w2-w1) > 0) {
    #    # if baseline slope is positive, slightly increase w1
    #    w1 <- max(x[x %in% w1:(w1+20) & y <= (a1 + 0.005*max(y[x > 1800]))])
    #    w1.result <- list(minimum=w1, success=success.W1[i,"arom15"])
    #    a1 <- y[x == w1]
    #    W1[i,"arom15"] <- w1
    #    A1[i,"arom15"] <- a1
    #    success.W1[i,"arom15"] <- w1.result[[2]]
    #  } else {
    #    # if baseline slope is negative, slightly decrease w2
    #    w2 <- min(x[x %in% (w2-20):w2 & y <= (a2 + 0.005*max(y[x > 1800]))])
    #    w2.result <- list(minimum=w2, success=success.W2[i,"arom15"])
    #    a2 <- y[x == w2]
    #    W2[i,"arom15"] <- w2
    #    A2[i,"arom15"] <- a2
    #    success.W2[i,"arom15"] <- w2.result[[2]]
    #  }

    # baseline correction calculated for all wavenumbers in peak
    w.wholepeak <- w2:w1
    a.wholepeak <- y[x %in% w.wholepeak]
    a.corrpeak <- a.wholepeak-((a2-a1)*(w.wholepeak-w1)/(w2-w1)+a1)

    # peak areas
    area.wholepeak[i,"arom15"] <- sum(a.wholepeak)
    area.corrpeak[i,"arom15"] <- sum(a.corrpeak)

    # pick out peak location
    Acorr[i,"arom15"] <- max(a.corrpeak)
    if(Acorr[i,"arom15"]==0) {
      #what to do if there is no arom15 peak, i.e. a.corrpeak is negative except at endpoints
      Wp[i,"arom15"] <- 1510 # use this as a default value for Wp
      Acorr[i,"arom15"] <- a.corrpeak[w.wholepeak==Wp[i,"arom15"]] # redefine Acorr at Wp
      Ap[i,"arom15"] <- a.wholepeak[w.wholepeak==Wp[i,"arom15"]]
      Ab[i,"arom15"] <- Ap[i,"arom15"] - Acorr[i,"arom15"]
    } else {   # normal arom15 peak
      Wp[i,"arom15"] <- min(w.wholepeak[a.corrpeak==Acorr[i,"arom15"]])
      Ap[i,"arom15"] <- a.wholepeak[w.wholepeak==Wp[i,"arom15"]]
      Ab[i,"arom15"] <- Ap[i,"arom15"] - Acorr[i,"arom15"]
    }

    rm(w1.result, w2.result, w1, w2, a1, a2, w.wholepeak, a.wholepeak, a.corrpeak)


    # 1630 and 1720 peaks ("arom16" and "acids") ----

    # Processes of defining w1 and w2 are a bit different, but still use the same variable names:
    if (A2[i,"arom15"] != min(y[x %in% W2[i,"arom15"]:1590])) {
      # case where there is another local minimum in y after w2 of the arom15 peak
      w1.result <- list(minimum=x[y==min(y[x %in% W2[i,"arom15"]:1590])], success=TRUE)
    } else {
      # w1 for these peaks is the w2 of the arom15 peak
      w1.result <- list(minimum=W2[i,"arom15"], success=success.W2[i,"arom15"])
    }
    w1 <- max(w1.result[[1]])
    a1 <- y[x == w1]
    W1[i,"arom16"] <- w1
    A1[i,"arom16"] <- a1
    success.W1[i,"arom16"] <- w1.result[[2]]
    W1[i,"trough16"] <- w1
    A1[i,"trough16"] <- a1
    success.W1[i,"trough16"] <- w1.result[[2]]
    W1[i,"acids"] <- w1
    A1[i,"acids"] <- a1
    success.W1[i,"acids"] <- w1.result[[2]]

    # w2 for these peaks based on first wavenumber that is "close" to the absolute minimum
    w2 <- min(x[x %in% 1760:1850 & y <= (min(y[x %in% 1760:1850]) + 0.01*max(y[x > 1800]))])
    w2.result <- list(minimum=w2, success=TRUE)
    a2 <- y[x == w2]
    W2[i,"arom16"] <- w2
    A2[i,"arom16"] <- a2
    success.W2[i,"arom16"] <- w2.result[[2]]
    W2[i,"trough16"] <- w2
    A2[i,"trough16"] <- a2
    success.W2[i,"trough16"] <- w2.result[[2]]
    W2[i,"acids"] <- w2
    A2[i,"acids"] <- a2
    success.W2[i,"acids"] <- w2.result[[2]]

    # baseline correction calculated for all wavenumbers in peak
    w.wholepeak <- w2:w1
    a.wholepeak <- y[x %in% w.wholepeak]
    a.corrpeak <- a.wholepeak-((a2-a1)*(w.wholepeak-w1)/(w2-w1)+a1)

    # Now comes the hard part:
    # Separately find both peaks and the trough.
    # The arom16 peak is probably easiest, so find this first.

    # pick out 1630 peak location (only check wavenumbers up to a set point)
    Acorr[i,"arom16"] <- max(a.corrpeak[w.wholepeak < 1660])
    if(Acorr[i,"arom16"] == 0) {
      # what to do if there is no arom16 peak, i.e. a.corrpeak is negative except at endpoint
      Wp[i,"arom16"] <- 1630 # use this as a default value for Wp
      Acorr[i,"arom16"] <- a.corrpeak[w.wholepeak==Wp[i,"arom16"]] # redefine Acorr at Wp
      Ap[i,"arom16"] <- a.wholepeak[w.wholepeak==Wp[i,"arom16"]]
      Ab[i,"arom16"] <- Ap[i,"arom16"] - Acorr[i,"arom16"]
    } else {   # normal arom16 peak
      Wp[i,"arom16"] <- min(w.wholepeak[a.corrpeak==Acorr[i,"arom16"]])
      Ap[i,"arom16"] <- a.wholepeak[w.wholepeak==Wp[i,"arom16"]]
      Ab[i,"arom16"] <- Ap[i,"arom16"] - Acorr[i,"arom16"]
    }

    # find an initial guess for the trough using the minima function
    # use a range of arom16+20 --> min(1740, (w2-20)) (min with 1740 is new, MAY CHANGE THIS LATER)

    default.trough16 <- 1685
    default.acids <- 1725

    tr.ini <- irp_content_klh_hodgkins_minima(xydata[x %in% max(1650, Wp[i,"arom16"]+20):min(1740, (w2-20)),], default=default.trough16)[[1]]

    # pick out acids peak location
    if(all(dy[x %in% max(1650, Wp[i,"arom16"]+20):min(1740, (w2-20))] < 0)) {
      # weak or absent acids peak, i.e. dy is always negative
      acids.result <- irp_content_klh_hodgkins_maxima(xydata[x %in% max(tr.ini, 1685):w2,], default=default.acids)
      Wp[i,"acids"] <- min(acids.result[[1]])
      Acorr[i,"acids"] <- a.corrpeak[w.wholepeak==Wp[i,"acids"]]
      Ap[i,"acids"] <- a.wholepeak[w.wholepeak==Wp[i,"acids"]]
      Ab[i,"acids"] <- Ap[i,"acids"] - Acorr[i,"acids"]
      notes[i, "acids.type"] <- ifelse(acids.result[[2]], "shoulder", "no peak")
      rm(acids.result)
    } else {
      # case where there is a local maximum in region of acids peak, i.e. dy > 0 in part of this region
      if(max(a.corrpeak[w.wholepeak > max(tr.ini, 1685)]) != 0) {
        # case where there is a true/normal acids peak above the baseline
        Acorr[i,"acids"] <- max(a.corrpeak[w.wholepeak > max(tr.ini, 1685)])
        Wp[i,"acids"] <- min(w.wholepeak[a.corrpeak==Acorr[i,"acids"]])
        Ap[i,"acids"] <- a.wholepeak[w.wholepeak==Wp[i,"acids"]]
      } else {
        # case where there is a true acids peak, but it is below baseline
        # in this case, non-baseline-corrected absorbances are used to find the peak
        Ap[i,"acids"] <- max(a.wholepeak[w.wholepeak > max(tr.ini, 1685)])
        Wp[i,"acids"] <- min(w.wholepeak[a.wholepeak==Ap[i,"acids"]])
        Acorr[i,"acids"] <- a.corrpeak[w.wholepeak==Wp[i,"acids"]]
      }
      Ab[i,"acids"] <- Ap[i,"acids"] - Acorr[i,"acids"]
      notes[i, "acids.type"] <- "peak"
    }

    # pick out "true" trough location, defined based on whether acids is a real peak
    if(notes[i, "acids.type"] == "peak") {
      Acorr[i,"trough16"] <- min(a.corrpeak[w.wholepeak %in% Wp[i,"arom16"]:Wp[i,"acids"]])
      Wp[i,"trough16"] <- min(w.wholepeak[a.corrpeak==Acorr[i,"trough16"]])
    } else if (notes[i, "acids.type"] == "shoulder") {
      Wp[i,"trough16"] <- min(tr.ini)
      Acorr[i,"trough16"] <- a.corrpeak[w.wholepeak==Wp[i,"trough16"]]
    } else { # acids.type is "no peak"
      Wp[i,"trough16"] <- default.trough16
      Acorr[i,"trough16"] <- a.corrpeak[w.wholepeak==Wp[i,"trough16"]]
    }
    Ap[i,"trough16"] <- a.wholepeak[w.wholepeak==Wp[i,"trough16"]]
    Ab[i,"trough16"] <- Ap[i,"trough16"] - Acorr[i,"trough16"]

    # append "negative" to acids.type if Acorr is negative for either acids or trough16
    if(Acorr[i,"trough16"] < 0 | Acorr[i,"acids"] < 0) {
      notes[i, "acids.type"] <- paste0(notes[i, "acids.type"], ", negative")
    }

    # peak area of arom16 (<= trough16)
    area.wholepeak[i,"arom16"] <- sum(a.wholepeak[w.wholepeak %in% w1:Wp[i,"trough16"]])
    area.corrpeak[i,"arom16"] <- sum(a.corrpeak[w.wholepeak %in% w1:Wp[i,"trough16"]])
    # peak area of acids (> trough16)
    area.wholepeak[i,"acids"] <- sum(a.wholepeak[w.wholepeak %in% (Wp[i,"trough16"]+1):w2])
    area.corrpeak[i,"acids"] <- sum(a.corrpeak[w.wholepeak %in% (Wp[i,"trough16"]+1):w2])

    rm(w1.result, w2.result, w1, w2, a1, a2, w.wholepeak, a.wholepeak, a.corrpeak, tr.ini, default.trough16, default.acids)


    # 2850 and 2920 peaks ("aliph28" and "aliph29") ----

    # w1 for these peaks is a constant, which can be changed globally if needbe
    w1.result <- list(minimum=2750, success=TRUE)
    w1 <- max(w1.result[[1]])
    a1 <- y[x == w1]
    W1[i,"aliph28"] <- w1
    A1[i,"aliph28"] <- a1
    success.W1[i,"aliph28"] <- w1.result[[2]]
    W1[i,"trough28"] <- w1
    A1[i,"trough28"] <- a1
    success.W1[i,"trough28"] <- w1.result[[2]]
    W1[i,"aliph29"] <- w1
    A1[i,"aliph29"] <- a1
    success.W1[i,"aliph29"] <- w1.result[[2]]

    #w2 for these peaks is based on local minimum
    w2.result <- irp_content_klh_hodgkins_minima(xydata[x %in% 2950:3050,], default=3000)
    w2 <- min(w2.result[[1]])
    a2 <- y[x == w2]
    W2[i,"aliph28"] <- w2
    A2[i,"aliph28"] <- a2
    success.W2[i,"aliph28"] <- w2.result[[2]]
    W2[i,"trough28"] <- w2
    A2[i,"trough28"] <- a2
    success.W2[i,"trough28"] <- w2.result[[2]]
    W2[i,"aliph29"] <- w2
    A2[i,"aliph29"] <- a2
    success.W2[i,"aliph29"] <- w2.result[[2]]

    # baseline correction calculated for all wavenumbers in peak
    w.wholepeak <- w2:w1
    a.wholepeak <- y[x %in% w.wholepeak]
    a.corrpeak <- a.wholepeak-((a2-a1)*(w.wholepeak-w1)/(w2-w1)+a1)

    # pick out peak locations

    # The aliphatic region can be one of 3 basic shapes.
    # To test them, perform 2 checks:
    # (1) check if maximum is close to 2920
    # (2) if max is close to 2920, check if there is any local minimum in a.wholepeak or (if not) a.corrpeak between 2850 and 2920

    if(max(w.wholepeak[a.wholepeak==max(a.wholepeak)]) %in% 2915:w2) {
      # case where maximum is close to 2920 (or higher)
      if(!any(w.wholepeak[a.wholepeak==min(a.wholepeak[w.wholepeak %in% 2854:2913])] %in% c(2854,2913))) {
        # case where there is a local minimum in a.wholepeak between 2854 and 2913
        # (i.e., the absolute minimum of this range is not one of the endpoints)
        type <- "two separated peaks"
      } else if(sum(y[x %in% 2920:2850]-((y[x==2920]-y[x==2850])*((2920:2850)-2850)/(2920-2850)+y[x==2850])) < 0) {
        # case where there is no local minimum in a.wholepeak,
        # but absorbances between 2850 and 2920 are below a diagonal between these points
        type <- "two unseparated peaks"
      } else {
        # case where none of mid-region is concave up, indicating singlet with max at 2920
        type <- "one peak"
      }
    } else {
      # case where maximum is not close to 2920
      type <- "one peak"
    }
    notes[i,"aliph.type"] <- type

    if(type == "one peak") {
      # set aliph29 as peak maximum, regardless of where this maximum is
      # then set trough28 and aliph28 the same as aliph29
      Acorr[i,"aliph29"] <- max(a.corrpeak)
      Wp[i,"aliph29"] <- max(w.wholepeak[a.corrpeak==Acorr[i,"aliph29"]])
      Ap[i,"aliph29"] <- a.wholepeak[w.wholepeak==Wp[i,"aliph29"]]
      Ab[i,"aliph29"] <- Ap[i,"aliph29"] - Acorr[i,"aliph29"]

      Acorr[i,"trough28"] <- Acorr[i,"aliph29"]
      Wp[i,"trough28"] <- Wp[i,"aliph29"]
      Ap[i,"trough28"] <- Ap[i,"aliph29"]
      Ab[i,"trough28"] <- Ab[i,"aliph29"]

      Acorr[i,"aliph28"] <- Acorr[i,"aliph29"]
      Wp[i,"aliph28"] <- Wp[i,"aliph29"]
      Ap[i,"aliph28"] <- Ap[i,"aliph29"]
      Ab[i,"aliph28"] <- Ab[i,"aliph29"]

    } else if(type == "two separated peaks") {
      Acorr[i,"aliph29"] <- max(a.corrpeak[w.wholepeak %in% 2900:w2])
      Wp[i,"aliph29"] <- max(w.wholepeak[a.corrpeak==Acorr[i,"aliph29"]])
      Ap[i,"aliph29"] <- a.wholepeak[w.wholepeak==Wp[i,"aliph29"]]
      Ab[i,"aliph29"] <- Ap[i,"aliph29"] - Acorr[i,"aliph29"]

      Acorr[i,"trough28"] <- min(a.corrpeak[w.wholepeak %in% 2854:Wp[i,"aliph29"]])
      Wp[i,"trough28"] <- max(w.wholepeak[a.corrpeak==Acorr[i,"trough28"]])
      Ap[i,"trough28"] <- a.wholepeak[w.wholepeak== Wp[i,"trough28"]]
      Ab[i,"trough28"] <- Ap[i,"trough28"] - Acorr[i,"trough28"]

      Acorr[i, "aliph28"] <- max(a.corrpeak[w.wholepeak %in% w1:Wp[i,"trough28"]])
      Wp[i,"aliph28"] <- min(w.wholepeak[a.corrpeak==Acorr[i,"aliph28"]])
      Ap[i,"aliph28"] <- a.wholepeak[w.wholepeak==Wp[i,"aliph28"]]
      Ab[i,"aliph28"] <- Ap[i,"aliph28"] - Acorr[i,"aliph28"]

    } else { # type == "two unseparated peaks"
      Acorr[i,"aliph29"] <- max(a.corrpeak[w.wholepeak %in% 2900:w2])
      Wp[i,"aliph29"] <- max(w.wholepeak[a.corrpeak==Acorr[i,"aliph29"]])
      Ap[i,"aliph29"] <- a.wholepeak[w.wholepeak==Wp[i,"aliph29"]]
      Ab[i,"aliph29"] <- Ap[i,"aliph29"] - Acorr[i,"aliph29"]

      Acorr[i,"trough28"] <- min(a.corrpeak[w.wholepeak %in% 2854:Wp[i,"aliph29"]])
      Wp[i,"trough28"] <- max(w.wholepeak[a.corrpeak==Acorr[i,"trough28"]])
      Ap[i,"trough28"] <- a.wholepeak[w.wholepeak== Wp[i,"trough28"]]
      Ab[i,"trough28"] <- Ap[i,"trough28"] - Acorr[i,"trough28"]

      # manually specify aliph28
      # since two unseparated peaks otherwise produces unpredictable results
      Wp[i,"aliph28"] <- 2850
      Acorr[i,"aliph28"] <- a.corrpeak[w.wholepeak==Wp[i,"aliph28"]]
      Ap[i,"aliph28"] <- a.wholepeak[w.wholepeak==Wp[i,"aliph28"]]
      Ab[i,"aliph28"] <- Ap[i,"aliph28"] - Acorr[i,"aliph28"]
    }

    # peak area of aliph28 (<= trough28)
    area.wholepeak[i,"aliph28"] <- sum(a.wholepeak[w.wholepeak %in% w1:Wp[i,"trough28"]])
    area.corrpeak[i,"aliph28"] <- sum(a.corrpeak[w.wholepeak %in% w1:Wp[i,"trough28"]])
    # peak area of aliph29 (> trough28)
    area.wholepeak[i,"aliph29"] <- sum(a.wholepeak[w.wholepeak %in% (Wp[i,"trough28"]+1):w2])
    area.corrpeak[i,"aliph29"] <- sum(a.corrpeak[w.wholepeak %in% (Wp[i,"trough28"]+1):w2])

    rm(w1.result, w2.result, w1, w2, a1, a2, w.wholepeak, a.wholepeak, a.corrpeak, type)


    rm(xydata, y, dy, d2y)


  }

  # Show plots for each peak in each sample ----
  if(verbose) {
    message(paste("Spectra with baselines and peak locations will be shown for all", length(data), "samples."))
  }
  # plots.prompt <- readline("Prompt for notes on each spectrum, which will be stored in output? (y/n) ")
  # while(!any(plots.prompt == c('y','Y','n','N'))) {
  #   plots.prompt <- readline(paste("Invalid response. Prompt for notes on each spectrum? Type Y or N. "))
  # }

  col.peaks <- rep(NA, length(peaks))
  names(col.peaks) <- peaks
  col.peaks['carb'] <- 'blue'
  col.peaks['arom15'] <- 'red'
  col.peaks['arom16'] <- 'purple3'
  col.peaks['trough16'] <- 'deepskyblue3'
  col.peaks['acids'] <- 'forestgreen'
  col.peaks['aliph28'] <- 'darkorange3'
  col.peaks['trough28'] <- 'gold4'
  col.peaks['aliph29'] <- 'orangered3'

  lty.peaks <- ifelse(peaks=="trough16" | peaks=="trough28", "dashed", "solid")
  names(lty.peaks) <- peaks

  if(make_plots) {
    for(i in seq_along(data)) {
      graphics::plot(x, data[,i], type="l", xaxp=c(700,4000,33), xaxs='i', yaxs='i', main=names(data)[i], xlab="wavenumber", ylab="Absorbance")

      for(j in seq_along(peaks)) {
        graphics::segments(W1[i,j], A1[i,j], W2[i,j], A2[i,j], col="green")                        # baseline
        graphics::segments(Wp[i,j], Ap[i,j], Wp[i,j], Ab[i,j], col=col.peaks[j], lty=lty.peaks[j]) # Acorr
        graphics::segments(W1[i,j], A1[i,j], W1[i,j], 0, lty="dotted", col="gray50") # connect endpoints to x-axis
        graphics::segments(W2[i,j], A2[i,j], W2[i,j], 0, lty="dotted", col="gray50") # connect endpoints to x-axis
        #     points(c(W1[i,j], W2[i,j]), c(A1[i,j], A2[i,j]), col="gray50")
        #     points(c(Wp[i,j], Wp[i,j]), c(Ap[i,j], Ab[i,j]), col="gray50")
      }
      # if(plots.prompt == 'y' | plots.prompt == 'Y') {
      #   notes[i,"notes"] <- readline("Do these peak baseline corrections look OK? Type any notes, or just press ENTER. ")
      # }
    }
  }

  # Calculate areas and mineral indices (780 peak, for silicates; constant and has no baseline) ----
  area <- apply(data, 2, sum)
  silicate780 <- as.numeric(data[x==780,])
  norm.silicate780 <- silicate780/area
  area_and_silicate <- data.frame(area, silicate780, norm.silicate780)

  # Normalize to area ----
  norm.Ap <- Ap/area
  norm.Acorr <- Acorr/area
  norm.area.wholepeak <- area.wholepeak/area
  norm.area.corrpeak <- area.corrpeak/area

  # Export data ----

  dataset.name <- export

  if(!is.null(export)) {

    # write.csv(data, file.path(dataset.name, filename))
    utils::write.csv(area_and_silicate, file.path(dataset.name, "TotalArea.and.780.csv"))
    utils::write.csv(Wp, file.path(dataset.name, "Wp.csv"))
    utils::write.csv(W1, file.path(dataset.name, "W1.csv"))
    utils::write.csv(W2, file.path(dataset.name, "W2.csv"))
    utils::write.csv(Ap, file.path(dataset.name, "Heights_Raw.csv"))
    utils::write.csv(Acorr, file.path(dataset.name, "Heights_Corr.csv"))
    utils::write.csv(norm.Ap, file.path(dataset.name, "Heights_Norm.Raw.csv"))
    utils::write.csv(norm.Acorr, file.path(dataset.name, "Heights_Norm.Corr.csv"))
    utils::write.csv(area.wholepeak, file.path(dataset.name, "Areas_Raw.csv"))
    utils::write.csv(area.corrpeak, file.path(dataset.name, "Areas_Corr.csv"))
    utils::write.csv(norm.area.wholepeak, file.path(dataset.name, "Areas_Norm.Raw.csv"))
    utils::write.csv(norm.area.corrpeak, file.path(dataset.name, "Areas_Norm.Corr.csv"))
    utils::write.csv(success.W1, file.path(dataset.name, "success.W1.csv"))
    utils::write.csv(success.W2, file.path(dataset.name, "success.W2.csv"))
    utils::write.csv(notes, file.path(dataset.name, "Notes.csv"))

  }

  list(
    area_and_silicate = area_and_silicate,
    Wp = Wp,
    W1 = W1,
    W2 = W2,
    Ap = Ap,
    Acorr = Acorr,
    norm.Ap = norm.Ap,
    norm.Acorr = norm.Acorr,
    area.wholepeak = area.wholepeak,
    area.corrpeak = area.corrpeak,
    norm.area.wholepeak = norm.area.wholepeak,
    norm.area.corrpeak = norm.area.corrpeak,
    success.W1 = success.W1,
    success.W2 = success.W2
  )

}
