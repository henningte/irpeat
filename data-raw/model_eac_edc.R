## code to prepare `model_eac_1` dataset goes here (including required datasets `y_eac_scale_attr_1` and `x_eac_scale_attr_1`)

# packages
library(redoxpeat)
library(dplyr)
library(errors)
library(quantities)
library(tibble)
library(ir)
library(elco)
library(magrittr)

# get data
d <- redoxpeat::d
d_mir <- redoxpeat::d_mir
el_t0 <- redoxpeat::el_t0
fe_t0 <- redoxpeat::fe_t0

#### preprocessing of electrochemical data ####

# bind to el_t0 to compute contribution to the EAC and EDC
el_t0 <-
  dplyr::left_join(el_t0,
                   fe_t0 %>%
                     dplyr::rename(mass_fe = mass) %>%
                     dplyr::select(id_90, fe2, fe3, mass_fe),
                   by = "id_90")

# compute the contribution of iron to the EAC and EDC
el_t0 <-
  el_t0 %>%
  dplyr::mutate(
    fe2_n = fe2/errors::drop_errors(elco::elco_drop_elco(C)),
    fe3_n = fe3/errors::drop_errors(elco::elco_drop_elco(C)),
    edc_c_uncorrected = edc_c,
    eac_c_uncorrected = eac_c,
    edc_c_fe_tot = edc_c - (fe2_n + fe3_n),
    eac_c_fe_tot = eac_c - (fe2_n + fe3_n),
    edc_c = edc_c - fe2_n,
    eac_c = eac_c - fe3_n,
    edc_m = edc_c * errors::drop_errors(elco::elco_drop_elco(C)),
    eac_m = eac_c * errors::drop_errors(elco::elco_drop_elco(C))
  )

# get an index for samples with too large contribution of iron to the EAC/EDC
p_fe_corrected_uncorrected_df <-
  dplyr::bind_rows(
    el_t0 %>%
      dplyr::group_by(id_90, site_label) %>%
      dplyr::select(id_90, site_label, eac_c_uncorrected, eac_c, eac_c_fe_tot) %>%
      dplyr::rename(uncorrected = eac_c_uncorrected,
                    corrected = eac_c,
                    corrected_fe_tot = eac_c_fe_tot) %>%
      dplyr::summarize_if(is.numeric,
                          .funs = list(median = mean, sd = sd),
                          na.rm = T) %>%
      dplyr::mutate(variable = "EAC"),
    el_t0 %>%
      dplyr::group_by(id_90, site_label) %>%
      dplyr::select(id_90, site_label, edc_c_uncorrected, edc_c, edc_c_fe_tot) %>%
      dplyr::rename(uncorrected = edc_c_uncorrected,
                    corrected = edc_c,
                    corrected_fe_tot = edc_c_fe_tot) %>%
      dplyr::summarize_if(is.numeric,
                          .funs = list(median = mean, sd = sd),
                          na.rm = T) %>%
      dplyr::mutate(variable = "EDC")
  ) %>%
  dplyr::mutate(
    uncorrected_lwr = as.numeric(uncorrected_median) - uncorrected_sd,
    uncorrected_upr = as.numeric(uncorrected_median) + uncorrected_sd,
    corrected_lwr = as.numeric(corrected_median) - corrected_sd,
    corrected_upr = as.numeric(corrected_median) + corrected_sd,
    corrected_fe_tot_lwr = as.numeric(corrected_fe_tot_median) - corrected_fe_tot_sd,
    corrected_fe_tot_upr = as.numeric(corrected_fe_tot_median) + corrected_fe_tot_sd,
    diff = uncorrected_median - corrected_median,
    diff_fe_tot = uncorrected_median - corrected_fe_tot_median
  )

index_fe_threshold <-
  tibble::tibble(
    id_90 = p_fe_corrected_uncorrected_df$id_90[p_fe_corrected_uncorrected_df$variable == "EAC"],
    index_eac = as.numeric(p_fe_corrected_uncorrected_df$diff_fe_tot)[p_fe_corrected_uncorrected_df$variable == "EAC"] <= 100,
    index_edc = as.numeric(p_fe_corrected_uncorrected_df$diff_fe_tot)[p_fe_corrected_uncorrected_df$variable == "EDC"] <= 100
  )

# discard outlier
index <- el_t0$eac_c > units::set_units(1500, "µmol/g")
el_t0 <-
  el_t0 %>%
  dplyr::mutate(eac = replace(eac, index, NA),
                eac_c = replace(eac_c, index, NA),
                eac_m = replace(eac_m, index, NA))

# get number edc measurements <=0
el_t0_negative <-
  el_t0 %>%
  dplyr::filter(el_t0$edc_c <= units::set_units(0, "µmol/g") & !is.na(el_t0$edc_c))

# set negative edc values to 0
el_t0$edc_c[el_t0$edc_c <= units::set_units(0, "µmol/g") & !is.na(el_t0$edc_c)] <- units::set_units(0, as.character(units(el_t0$edc_c)), mode = "standard")

# summarize
el_t0_summary <-
  el_t0 %>%
  dplyr::group_by(id_90) %>%
  dplyr::summarise(
    eac1 = quantities::set_quantities(mean(eac, na.rm = TRUE),
                                      unit = "µmol",
                                      errors = sd(eac, na.rm = TRUE)),
    eac_c1 = quantities::set_quantities(mean(eac_c, na.rm = TRUE),
                                        unit = "µmol/g",
                                        errors = sd(eac_c, na.rm = TRUE)),
    eac_m1 = quantities::set_quantities(mean(eac_m, na.rm = TRUE),
                                        unit = "µmol/g",
                                        errors = sd(eac_m, na.rm = TRUE)),
    edc1 = quantities::set_quantities(mean(edc, na.rm = TRUE),
                                      unit = "µmol",
                                      errors = sd(edc, na.rm = TRUE)),
    edc_c1 = quantities::set_quantities(mean(edc_c, na.rm = TRUE),
                                        unit = "µmol/g",
                                        errors = sd(edc_c, na.rm = TRUE)),
    edc_m1 = quantities::set_quantities(mean(edc_m, na.rm = TRUE),
                                        unit = "µmol/g",
                                        errors = sd(edc_m, na.rm = TRUE))
  )

#### mir preprocessing ####

# interpolation
# d_mir <-
#   d_mir %>%
#   ir::ir_interpolate(start = NULL, dw = 1)

# clipping to a common wavenumber range
# mir_target_range <-
#   tibble::tibble(
#     start = 800, # this range was chosen to remove artifacts from clay minerals which have a strong influence on the area normalized intensities
#     end = 4000
#   )

# d_mir <-
#   d_mir %>%
#   ir::ir_clip(mir_target_range)

# remove CO2 artifacts
# mir_interpolation_range <-
#   tibble::tibble(
#     start = 2290,
#     end = 2400
#   )

# d_mir <-
#   d_mir %>%
#   ir::ir_interpolate_region(range = mir_interpolation_range)

# baseline correction
# bc_range <-
#   tibble::tibble(
#     start = min(d_mir$spectra[[1]]$x) + 10,
#     end = max(d_mir$spectra[[1]]$x) - 10
#   )

# d_mir <-
#   d_mir %>%
#   ir::ir_bc(method = "rubberband") %>%
#   ir::ir_clip(range = bc_range) %>%
#   ir::ir_bc(method = "rubberband")

# normalization
# d_mir <-
#   d_mir %>%
#   ir::ir_normalize(method = "area")

# humification indices
# d_mir <-
#   d_mir %>%
#   irpeat::irp_hi()

#### merge data ####

# remove unneeded columns
d <-
  d %>%
  dplyr::select(site_label, site_name, site_order, depth_upper, depth_lower, id_90, id_ftir)
d_mir_combn <-
  d_mir %>%
  dplyr::rename(id_ftir = sample_id) %>%
  dplyr::mutate(id_ftir = as.character(id_ftir)) %>%
  dplyr::left_join(d %>%
                     dplyr::select(id_ftir, id_90),
                   by = c("id_ftir", "id_90"))

# combine all data for the subsequent analyses
d <-
  list(d, el_t0_summary, d_mir_combn) %>%
  purrr::reduce(dplyr::left_join, by = "id_90") %>%
  dplyr::select(-c(id_ftir.x, id_ftir.y))

# set units
d <-
  d %>%
  dplyr::mutate(
    depth_lower = units::set_units(depth_lower, "cm"),
    depth_upper = units::set_units(depth_upper, "cm")
  )

# define factors
d <-
  d %>%
  dplyr::mutate(
    site_label = factor(site_label, levels = unique(site_label[order(site_order)])),
    site_name = factor(site_name, levels = unique(site_name[order(site_order)]))
  )

# restore class
class(d) <- c("ir", class(d))

# add an index indicating which samples to use during the statistical analyses
d <-
  d %>%
  dplyr::mutate(
    index_eac = id_90 %in% {index_fe_threshold %>%
        dplyr::filter(index_eac) %>%
        dplyr::select(id_90) %>%
        unlist()},
    index_edc = id_90 %in% {index_fe_threshold %>%
        dplyr::filter(index_edc) %>%
        dplyr::select(id_90) %>%
        unlist()}
  )

#### prepare data for modeling ####

## parameters for MCMC sampling and Bayesian modeling

# define the number of expected non-zero coefficients
p0 <- 8

# define MCMC parameters
seed <- 1
chains <- 4
iter <- 2000
warmup <- 500

## data

## spectral preprocessing

# clipping range
mir_target_range <-
  tibble::tibble(
    start = 800, # this range was chosen to remove artifacts from clay minerals which have a strong influence on the area normalized intensities
    end = 4000
  )

# range where to remove CO2 artifacts by linear interpolation
mir_interpolation_range <-
  tibble::tibble(
    start = 2290,
    end = 2400
  )

# store configurations in a list
model_eac_1_config <-
  list(
    irp_preprocess =
      list(
        do_interpolate = TRUE,
        interpolate_start = NULL,
        interpolate_dw = 1,
        do_clip = TRUE,
        clip_range = mir_target_range,
        do_interpolate_region = TRUE,
        interpolate_region_range = mir_interpolation_range,
        do_bc = TRUE,
        bc_method = "rubberband",
        bc_cutoff = 10,
        do_smooth = FALSE,
        do_normalise = TRUE,
        normalise_method = "area",
        do_bin = TRUE,
        bin_width = 10,
        bin_new_x_type = "start",
        do_scale = TRUE,
        scale_center = TRUE,
        scale_scale = TRUE
      )
  )

d_cal_x <-
  d %>%
  dplyr::filter(index_eac) %>%
  irpeat::irp_preprocess(
    do_interpolate = model_eac_1_config$irp_preprocess$do_interpolate,
    interpolate_start = model_eac_1_config$irp_preprocess$interpolate_start,
    interpolate_dw = model_eac_1_config$irp_preprocess$interpolate_dw,
    do_clip = model_eac_1_config$irp_preprocess$do_clip,
    clip_range = model_eac_1_config$irp_preprocess$clip_range,
    do_interpolate_region = model_eac_1_config$irp_preprocess$do_interpolate_region,
    interpolate_region_range = model_eac_1_config$irp_preprocess$interpolate_region_range,
    do_bc = model_eac_1_config$irp_preprocess$do_bc,
    bc_method = model_eac_1_config$irp_preprocess$bc_method,
    bc_cutoff = model_eac_1_config$irp_preprocess$bc_cutoff,
    do_smooth = model_eac_1_config$irp_preprocess$do_smooth,
    do_normalise = model_eac_1_config$irp_preprocess$do_normalise,
    normalise_method = model_eac_1_config$irp_preprocess$normalise_method,
    do_bin = model_eac_1_config$irp_preprocess$do_bin,
    bin_width = model_eac_1_config$irp_preprocess$bin_width,
    bin_new_x_type = model_eac_1_config$irp_preprocess$bin_new_x_type,
    do_scale = model_eac_1_config$irp_preprocess$do_scale,
    scale_center = model_eac_1_config$irp_preprocess$scale_center,
    scale_scale = model_eac_1_config$irp_preprocess$scale_scale
  )

# bin the spectra
# d_cal <-
#   d %>%
#   ir::ir_bin(width = 10)

# EAC, non-derived spectra
d_cal_eac <-
  d %>%
  dplyr::filter(index_eac) %>%
  dplyr::mutate(x =
                  as.matrix(d_cal_x),
                y = quantities::drop_quantities(eac_c1),
                y_scaled = scale(y, center = TRUE, scale = TRUE),
                y_se = errors::errors(eac_c1)
  ) %>%
  dplyr::select(y, y_scaled, x)

model_eac_1_config$data_scale <-
  list(
    y_center = attr(d_cal_eac$y_scaled, "scaled:center"),
    y_scale = attr(d_cal_eac$y_scaled, "scaled:scale"),
    x_center = attr(d_cal_x, "scaled:center"),
    x_scale = attr(d_cal_x, "scaled:scale")
  )

#### model fit ####

x <- d_cal_eac

# define the prior distribution
cal_hs <-
  rstanarm::hs(df = 1,
               global_df = 1,
               global_scale = p0/(ncol(x$x) - p0) * 1/sqrt(nrow(x)),
               slab_df = 4,
               slab_scale = 1)

## compute the reference model
cal_ref_fit <-
  rstanarm::stan_glm(y_scaled ~ x,
                     family = gaussian(),
                     data = x,
                     prior = cal_hs,
                     seed = seed,
                     chains = chains,
                     iter = iter,
                     cores = chains,
                     warmup = warmup,
                     weights = NULL
  )


#### export preparation ####

# get the fitted model
model_eac_1 <- cal_ref_fit

# export
usethis::use_data(model_eac_1, overwrite = TRUE)
usethis::use_data(model_eac_1_config, overwrite = TRUE)

#####################################################

# store configurations in a list
model_edc_1_config <- model_eac_1_config

d_cal_x <-
  d %>%
  dplyr::filter(index_edc) %>%
  irpeat::irp_preprocess(
    do_interpolate = model_edc_1_config$irp_preprocess$do_interpolate,
    interpolate_start = model_edc_1_config$irp_preprocess$interpolate_start,
    interpolate_dw = model_edc_1_config$irp_preprocess$interpolate_dw,
    do_clip = model_edc_1_config$irp_preprocess$do_clip,
    clip_range = model_edc_1_config$irp_preprocess$clip_range,
    do_interpolate_region = model_edc_1_config$irp_preprocess$do_interpolate_region,
    interpolate_region_range = model_edc_1_config$irp_preprocess$interpolate_region_range,
    do_bc = model_edc_1_config$irp_preprocess$do_bc,
    bc_method = model_edc_1_config$irp_preprocess$bc_method,
    bc_cutoff = model_edc_1_config$irp_preprocess$bc_cutoff,
    do_smooth = model_edc_1_config$irp_preprocess$do_smooth,
    do_normalise = model_edc_1_config$irp_preprocess$do_normalise,
    normalise_method = model_edc_1_config$irp_preprocess$normalise_method,
    do_bin = model_edc_1_config$irp_preprocess$do_bin,
    bin_width = model_edc_1_config$irp_preprocess$bin_width,
    bin_new_x_type = model_edc_1_config$irp_preprocess$bin_new_x_type,
    do_scale = model_edc_1_config$irp_preprocess$do_scale,
    scale_center = model_edc_1_config$irp_preprocess$scale_center,
    scale_scale = model_edc_1_config$irp_preprocess$scale_scale
  )

# EDC, non-derived spectra
d_cal_edc <-
  d %>%
  dplyr::filter(index_edc) %>%
  dplyr::mutate(x =
                  as.matrix(d_cal_x),
                y = quantities::drop_quantities(edc_c1),
                y_scaled = scale(y, center = TRUE, scale = TRUE),
                y_se = errors::errors(edc_c1)
  ) %>%
  dplyr::select(y, y_scaled, x)

model_edc_1_config$data_scale <-
  list(
    y_center = attr(d_cal_edc$y_scaled, "scaled:center"),
    y_scale = attr(d_cal_edc$y_scaled, "scaled:scale"),
    x_center = attr(d_cal_x, "scaled:center"),
    x_scale = attr(d_cal_x, "scaled:scale")
  )

#### model fit ####

x <- d_cal_edc

# define the prior distribution
cal_hs <-
  rstanarm::hs(df = 1,
               global_df = 1,
               global_scale = p0/(ncol(x$x) - p0) * 1/sqrt(nrow(x)),
               slab_df = 4,
               slab_scale = 1)

## compute the reference model
cal_ref_fit <-
  rstanarm::stan_glm(y_scaled ~ x,
                     family = gaussian(),
                     data = x,
                     prior = cal_hs,
                     seed = seed,
                     chains = chains,
                     iter = iter,
                     cores = chains,
                     warmup = warmup,
                     weights = NULL
  )


#### export preparation ####

# get the fitted model
model_edc_1 <- cal_ref_fit

# export
usethis::use_data(model_edc_1, overwrite = TRUE)
usethis::use_data(model_edc_1_config, overwrite = TRUE)
