## code to prepare `model_holocellulose_2` dataset goes here

# packages
library(dplyr)
library(tibble)
library(ir)
library(brms)
library(magrittr)

# get data
d <- ir::ir_sample_data

#### preprocessing of holocellulose content data ####

#### merge data ####

#### prepare data for modeling ####

## parameters for MCMC sampling and Bayesian modeling

# define MCMC parameters
seed <- 1
chains <- 4
iter <- 3000
warmup <- iter %/% 2
control <- list(adapt_delta = 0.99, max_treedepth = 11)

## data

# index to filter data
index <- d$sample_type != "old magazines"

## spectral preprocessing

# store configurations in a list
model_holocellulose_2_config <-
  list(
    irp_preprocess =
      list(
        do_interpolate = TRUE,
        interpolate_start = NULL,
        interpolate_dw = 1,
        do_clip = FALSE,
        # clip_range = mir_target_range,
        do_interpolate_region = FALSE,
        # interpolate_region_range = mir_interpolation_range,
        do_bc = TRUE,
        bc_method = "rubberband",
        bc_cutoff = 0,
        do_smooth = FALSE,
        do_normalise = TRUE,
        normalise_method = "area",
        do_bin = TRUE,
        bin_width = 50, # ---todo: check
        do_scale = TRUE,
        scale_center = TRUE,
        scale_scale = TRUE
      )
  )

x <-
  d %>%
  irpeat::irp_preprocess(
    do_interpolate = model_holocellulose_2_config$irp_preprocess$do_interpolate,
    interpolate_start = model_holocellulose_2_config$irp_preprocess$interpolate_start,
    interpolate_dw = model_holocellulose_2_config$irp_preprocess$interpolate_dw,
    do_clip = model_holocellulose_2_config$irp_preprocess$do_clip,
    clip_range = model_holocellulose_2_config$irp_preprocess$clip_range,
    do_interpolate_region = model_holocellulose_2_config$irp_preprocess$do_interpolate_region,
    interpolate_region_range = model_holocellulose_2_config$irp_preprocess$interpolate_region_range,
    do_bc = model_holocellulose_2_config$irp_preprocess$do_bc,
    bc_method = model_holocellulose_2_config$irp_preprocess$bc_method,
    bc_cutoff = model_holocellulose_2_config$irp_preprocess$bc_cutoff,
    do_smooth = model_holocellulose_2_config$irp_preprocess$do_smooth,
    do_normalise = model_holocellulose_2_config$irp_preprocess$do_normalise,
    normalise_method = model_holocellulose_2_config$irp_preprocess$normalise_method,
    do_bin = model_holocellulose_2_config$irp_preprocess$do_bin,
    bin_width = model_holocellulose_2_config$irp_preprocess$bin_width,
    do_scale = model_holocellulose_2_config$irp_preprocess$do_scale,
    scale_center = model_holocellulose_2_config$irp_preprocess$scale_center,
    scale_scale = model_holocellulose_2_config$irp_preprocess$scale_scale
  )

d_cal <-
  d %>%
  dplyr::mutate(
    x = as.matrix(x),
    y =
      holocellulose %>%
      units::drop_units() %>%
      scale(center = FALSE, scale = FALSE)
  ) %>%
  dplyr::filter(index) %>%
  dplyr::select(y, x)

model_holocellulose_2_config$data_scale <-
  list(
    y_center = attr(d_cal$y, "scaled:center"),
    y_scale = attr(d_cal$y, "scaled:scale"),
    x_center = attr(x, "scaled:center"),
    x_scale = attr(x, "scaled:scale")
  )

#### model fit ####

# define priors
### START NEW
n <- nrow(d_cal$x)
p <- ncol(d_cal$x)
p0 <- 8L # prior guess for the number of relevant variables
tau0 <- p0/(p-p0) * 1/sqrt(n) # tau0

priors <-
  c(
    brms::prior_string(paste0("horseshoe(df = 1, scale_global = ", tau0, ", df_global = 1, autoscale = TRUE)"),
                       class = "b"),
    brms::prior_string("normal(0, 2.5)", class = "Intercept"),
    brms::prior_string("gamma(0.01, 0.01)", class = "phi")
  )
### END NEW

## fit the model
m <-
  brms::brm(y ~ .,
            data = d_cal,
            family = Beta(link = "logit", link_phi = "log"),
            prior = priors,
            cores = chains,
            chains = chains,
            warmup = warmup,
            iter = iter,
            seed = seed,
            control = control)

#### export preparation ####

# get the fitted model
model_holocellulose_2 <- m

# export
usethis::use_data(model_holocellulose_2, overwrite = TRUE)
usethis::use_data(model_holocellulose_2_config, overwrite = TRUE)
