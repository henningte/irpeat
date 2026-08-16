# Predicts peat properties from transmission mid-infrared spectra

Functions to predict peat properties from transmission mid-infrared
spectra. All functions below have been computed using peat samples. For
detailed information on the underlying prediction models, see the
details section.

## Usage

``` r
irp_holocellulose_content_2(
  x,
  ...,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  return_as_list = FALSE,
  check_prediction_domain = "train"
)

irp_klason_lignin_content_2(
  x,
  ...,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  return_as_list = FALSE,
  check_prediction_domain = "train"
)

irp_eac_1(
  x,
  ...,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  return_as_list = FALSE,
  check_prediction_domain = "train"
)

irp_edc_1(
  x,
  ...,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  return_as_list = FALSE,
  check_prediction_domain = "train"
)

irp_carbon_content_1(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE
)

irp_nitrogen_content_1(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE
)

irp_hydrogen_content_1(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE
)

irp_oxygen_content_1(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE
)

irp_phosphorus_content_1(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE
)

irp_potassium_content_1(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE
)

irp_sulfur_content_1(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE
)

irp_titanium_content_1(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE
)

irp_silicon_content_1(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE
)

irp_calcium_content_1(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE
)

irp_d13C_1(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE
)

irp_d15N_1(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE
)

irp_nosc_1(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE
)

irp_dgf0_1(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE
)

irp_bulk_density_1(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE
)

irp_loss_on_ignition_1(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE
)

irp_O_to_C_1(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE
)

irp_C_to_N_1(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE
)

irp_H_to_C_1(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE
)

irp_volume_fraction_solids_1(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE,
  bulk_density = NULL
)

irp_non_macroporosity_1(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE,
  bulk_density = NULL
)

irp_macroporosity_1(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE,
  bulk_density = NULL
)

irp_saturated_hydraulic_conductivity_1(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE,
  bulk_density = NULL
)

irp_specific_heat_capacity_1(
  x,
  temperature = 273.15,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE,
  nitrogen_content = NULL
)

irp_dry_thermal_conductivity_1(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE,
  bulk_density = NULL
)

irp_microbial_nitrogen_content_1(
  x,
  y,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE
)

irp_degree_of_decomposition_1(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE
)

irp_degree_of_decomposition_2(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE
)

irp_degree_of_decomposition_3(
  x,
  do_summary = FALSE,
  summary_function_mean = mean,
  summary_function_sd = stats::sd,
  check_prediction_domain = "train",
  return_as_list = FALSE
)
```

## Source

- `irp_holocellulose_2()`, `irp_klason_lignin_2()`:

  Teickner and Knorr (2022) .

- `irp_eac_1()`, `irp_edc_1()`:

  Teickner et al. (2022) .

- `irp_microbial_nitrogen_content_1()`:

  Reuter et al. (2020)

- `irp_degree_of_decomposition_1()`, `irp_degree_of_decomposition_2()`,
  `irp_degree_of_decomposition_3()`:

  Teickner et al. (2025)

- All other models:

  Teickner and Knorr (2025)

## Arguments

- x:

  An object of class
  [`ir`](https://henningte.github.io/ir/reference/ir_new_ir.html) with
  transmission mid infrared spectra. Some tests are applied to check if
  the supplied spectra match the spectra used to fit the models (the
  spectral range is checked). The spectral resolution of the original
  spectral data should not be smaller than 4 cm\\^{-1}\\ and it is not
  checked if this assumption is met. For the following models, `x` has
  special meaning:

  `irp_microbial_nitrogen_content_1()`

  :   Here, `x` is a set of litter spectra after decomposition (and `y`
      is a set of litter spectra before decomposition). See the Details
      section.

- ...:

  Additional arguments passed to
  [`rstanarm::posterior_predict.stanreg()`](https://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.html)
  (`irp_eac_1()`,`irp_eac_2()`).

- do_summary:

  A logical value indicating if the predicted values should be returned
  in a summarized version (`TRUE`) or not (`FALSE`).

  - If `do_summary = FALSE`, a list column is returned and each element
    of the list column is a numeric vector, or an `rvar` object is
    returned as column in `x`, depending on the value of
    `return_as_list`. In both cases, the column contains draws from the
    posterior predictive distribution.

  - If `do_summary = TRUE`, each element is a
    [`quantities::quantities()`](https://r-quantities.github.io/quantities/reference/quantities.html)
    object with value and error summarized from posterior draws via
    `summary_function_mean` and `summary_function_sd`

- summary_function_mean:

  A function used to summarize the predicted values (average).

- summary_function_sd:

  A function used to summarize the predicted values (spread).

- return_as_list:

  Logical value. If set to `TRUE`, the result will be returned as list
  of draws, otherwise the result will be returned as `rvar` object. This
  is a new argument currently only implemented for models predicting the
  degree of decomposition.

- check_prediction_domain:

  A character value indicating if and how it should be checked whether
  the spectra in `x` are within the prediction domain of the model. One
  of:

  `"train"`

  :   It is checked whether the spectra in `x` are within the prediction
      domain formed by the training data for the model.

  `"test"`

  :   It is checked whether the spectra in `x` are within the prediction
      domain formed by the testing data for the model.

  `"none"`

  :   It is not checked whether the spectra in `x` are within the
      prediction domain for the model.

- bulk_density:

  For `irp_porosity_1()`, `irp_non_macroporosity_1()`,
  `irp_macroporosity_1()`, `irp_volume_fraction_of_solids_1()`,
  `irp_saturated_hydraulic_conductivity_1()`,
  `irp_dry_thermal_conductivity_1()`: One of:

  1\.

  :   A numeric vector with the same number of elements as spectra in
      `x` with values for the dry bulk density in g cm\\^{-3}\\. These
      values will be used to predict the peat property.

  2\.

  :   A list with the same number of elements as spectra in `x`. Each
      element must be a numeric vector with the same number of elements
      as there are MCMC draws in the corresponding model to use. Each
      numeric value is the the dry bulk density in g cm\\^{-3}\\. These
      values will be used to predict the peat property.

  3\.

  :   `NULL`: Dry bulk density will be estimated from the spectra in `x`
      and these estimates will be used to predict the peat property.

- temperature:

  For `irp_specific_heat_capacity_1()`: The temperature in K for which
  to predict the specific heat capacity.

- nitrogen_content:

  For `irp_specific_heat_capacity_1()`: One of:

  1\.

  :   A numeric vector with the same number of elements as spectra in
      `x` with values for the nitrogen content in g g\\^{-1}\\. These
      values will be used to predict the peat property.

  2\.

  :   A list with the same number of elements as spectra in `x`. Each
      element must be a numeric vector with the same number of elements
      as there are MCMC draws in the corresponding model to use. Each
      numeric value is the the nitrogen content in g g\\^{-1}\\. These
      values will be used to predict the peat property.

  3\.

  :   `NULL`: Nitrogen content will be estimated from the spectra in `x`
      and these estimates will be used to predict the peat property.

- y:

  An object of class
  [`ir`](https://henningte.github.io/ir/reference/ir_new_ir.html) with
  transmission mid infrared spectra. This argument is required for the
  following models which need more than one set of spectra to make
  predictions:

  `irp_microbial_nitrogen_content_1()`

  :   Here, `y` is a set of litter spectra before decomposition.

## Value

`x` with a new column with the predicted peat property and a new column
(with the same name as the predicted peat property and ending in `_pd`)
with value `TRUE` if the respective spectrum is within the prediction
domain for the model and `FALSE` if not (see argument
`check_prediction_domain` and section 2.6 in (Teickner and Knorr 2025)
for more details). If `check_prediction_domain = "none"`, all values in
this column are `NA`.

## Details

The models use the models of the same name in the 'irpeatmodels'
package. The 'irpeatmodels' package provides information on the models.

## Note

- `irp_eac_1()`, `irp_edc_1()`:

  The model still has a relatively large uncertainty because it is
  fitted with few samples (Teickner et al. 2022) . For further
  limitations, see Teickner et al. (2022) .

- `irp_microbial_nitrogen_content_1()`:

  Reuter et al. (2020) describes limitations and uncertainties: “Small
  method modifications should be considered for the applicability of the
  method in aerobic decomposition studies. These modifications include
  an optimization of the calibration curve, either through the addition
  of very low N litters to a decomposition study as calibration samples
  or through the artificial mixing of undecomposed litter with microbial
  biomass. Furthermore, the contribution of fungi must be considered,
  which we assumed to be negligible in anoxic soils. Differences in the
  amount of DNA per biomass units and in the C/N ratio should be
  considered for the decomposer biomass in aerobic systems. Finally, the
  applicability of the same calibration curve for decomposed litters of
  different plant species still has to be investigated.”

## References

Reuter H, Gensel J, Elvert M, Zak D (2020). “Evidence for Preferential
Protein Depolymerization in Wetland Soils in Response to External
Nitrogen Availability Provided by a Novel FTIR Routine.”
*Biogeosciences*, **17**(2), 499–514. ISSN 1726-4189.
[doi:10.5194/bg-17-499-2020](https://doi.org/10.5194/bg-17-499-2020) .
[2022-02-15](https://henningte.github.io/irpeat/reference/2022-02-15).  
  
Teickner H, Arsenault J, Gałka M, Knorr K (2025). “Estimation of the
Degree of Decomposition of Peat and Past Net Primary Production from
Mid-Infrared Spectra.”  
  
Teickner H, Gao C, Knorr K (2022). “Electrochemical Properties of Peat
Particulate Organic Matter on a Global Scale: Relation to Peat Chemistry
and Degree of Decomposition.” *Global Biogeochemical Cycles*, **36**(2),
e2021GB007160. ISSN 0886-6236, 1944-9224.
[doi:10.1029/2021GB007160](https://doi.org/10.1029/2021GB007160) .
[2022-02-03](https://henningte.github.io/irpeat/reference/2022-02-03).  
  
Teickner H, Knorr K (2022). “Improving Models to Predict Holocellulose
and Klason Lignin Contents for Peat Soil Organic Matter with
Mid-Infrared Spectra.” *SOIL*, **8**(2), 699–715.
[doi:10.5194/soil-8-699-2022](https://doi.org/10.5194/soil-8-699-2022)
.  
  
Teickner H, Knorr K (2025). “Prediction of Peat Properties from
Transmission Mid-Infrared Spectra.”

## See also

[`irp_predict()`](https://henningte.github.io/irpeat/reference/irp_predict.md)

## Examples

``` r
library(ir)

x <- ir::ir_sample_data[1, ]

## make predictions

# holocellulose content
x <- irpeat::irp_holocellulose_content_2(
  x,
  do_summary = TRUE,
  check_prediction_domain = "train"
)

# Klason lignin content
x <- irpeat::irp_klason_lignin_content_2(
  x,
  do_summary = TRUE,
  check_prediction_domain = "train"
)

# electron accepting capacity
x <- irpeat::irp_eac_1(
  x,
  do_summary = TRUE,
  check_prediction_domain = "train"
)

# electron donating capacity
x <- irpeat::irp_edc_1(
  x,
  do_summary = TRUE,
  check_prediction_domain = "train"
)

# carbon content
x <-
  irp_carbon_content_1(
    x,
    do_summary = TRUE,
    check_prediction_domain = "train"
  )
#> Warning: 650 selected instead of 645.

# nitrogen content
irpeat::irp_nitrogen_content_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train"
)
#> Warning: 650 selected instead of 645.
#> # A tibble: 1 × 13
#>   id_90 sample_id measurement_id spectra            C        H        N        O
#> * <int>     <int>          <int> <named lis> (err) [… (err) [… (err) [… (err) […
#> 1     1         1             23 <df>        0.47902… 0.05625… 0.00968… 0.39768…
#> # ℹ 5 more variables: S (err) [g/g], d15N <dbl>, d13C <dbl>,
#> #   nitrogen_content_1 (err) [g/g], nitrogen_content_1_in_pd <lgl>

# hydrogen content
irpeat::irp_hydrogen_content_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train"
)
#> Warning: 650 selected instead of 645.
#> # A tibble: 1 × 13
#>   id_90 sample_id measurement_id spectra            C        H        N        O
#> * <int>     <int>          <int> <named lis> (err) [… (err) [… (err) [… (err) […
#> 1     1         1             23 <df>        0.47902… 0.05625… 0.00968… 0.39768…
#> # ℹ 5 more variables: S (err) [g/g], d15N <dbl>, d13C <dbl>,
#> #   hydrogen_content_1 (err) [g/g], hydrogen_content_1_in_pd <lgl>

# oxygen content
irpeat::irp_oxygen_content_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train"
)
#> Warning: 650 selected instead of 645.
#> # A tibble: 1 × 13
#>   id_90 sample_id measurement_id spectra            C        H        N        O
#> * <int>     <int>          <int> <named lis> (err) [… (err) [… (err) [… (err) […
#> 1     1         1             23 <df>        0.47902… 0.05625… 0.00968… 0.39768…
#> # ℹ 5 more variables: S (err) [g/g], d15N <dbl>, d13C <dbl>,
#> #   oxygen_content_1 (err) [g/g], oxygen_content_1_in_pd <lgl>

# phosphorus content
irpeat::irp_phosphorus_content_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train"
)
#> Warning: 650 selected instead of 645.
#> # A tibble: 1 × 13
#>   id_90 sample_id measurement_id spectra            C        H        N        O
#> * <int>     <int>          <int> <named lis> (err) [… (err) [… (err) [… (err) […
#> 1     1         1             23 <df>        0.47902… 0.05625… 0.00968… 0.39768…
#> # ℹ 5 more variables: S (err) [g/g], d15N <dbl>, d13C <dbl>,
#> #   phosphorus_content_1 (err) [g/g], phosphorus_content_1_in_pd <lgl>

# potassium content
irpeat::irp_potassium_content_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train"
)
#> Warning: 650 selected instead of 645.
#> # A tibble: 1 × 13
#>   id_90 sample_id measurement_id spectra            C        H        N        O
#> * <int>     <int>          <int> <named lis> (err) [… (err) [… (err) [… (err) […
#> 1     1         1             23 <df>        0.47902… 0.05625… 0.00968… 0.39768…
#> # ℹ 5 more variables: S (err) [g/g], d15N <dbl>, d13C <dbl>,
#> #   potassium_content_1 (err) [g/g], potassium_content_1_in_pd <lgl>

# sulfur content
irpeat::irp_sulfur_content_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train"
)
#> Warning: 650 selected instead of 645.
#> # A tibble: 1 × 13
#>   id_90 sample_id measurement_id spectra            C        H        N        O
#> * <int>     <int>          <int> <named lis> (err) [… (err) [… (err) [… (err) […
#> 1     1         1             23 <df>        0.47902… 0.05625… 0.00968… 0.39768…
#> # ℹ 5 more variables: S (err) [g/g], d15N <dbl>, d13C <dbl>,
#> #   sulfur_content_1 (err) [g/g], sulfur_content_1_in_pd <lgl>

# titanium content
irpeat::irp_titanium_content_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train"
)
#> Warning: 650 selected instead of 645.
#> # A tibble: 1 × 13
#>   id_90 sample_id measurement_id spectra            C        H        N        O
#> * <int>     <int>          <int> <named lis> (err) [… (err) [… (err) [… (err) […
#> 1     1         1             23 <df>        0.47902… 0.05625… 0.00968… 0.39768…
#> # ℹ 5 more variables: S (err) [g/g], d15N <dbl>, d13C <dbl>,
#> #   titanium_content_1 (err) [g/g], titanium_content_1_in_pd <lgl>

# silicon content
irpeat::irp_silicon_content_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train"
)
#> Warning: 650 selected instead of 645.
#> # A tibble: 1 × 13
#>   id_90 sample_id measurement_id spectra            C        H        N        O
#> * <int>     <int>          <int> <named lis> (err) [… (err) [… (err) [… (err) […
#> 1     1         1             23 <df>        0.47902… 0.05625… 0.00968… 0.39768…
#> # ℹ 5 more variables: S (err) [g/g], d15N <dbl>, d13C <dbl>,
#> #   silicon_content_1 (err) [g/g], silicon_content_1_in_pd <lgl>

# calcium content
irpeat::irp_calcium_content_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train"
)
#> Warning: 650 selected instead of 645.
#> # A tibble: 1 × 13
#>   id_90 sample_id measurement_id spectra            C        H        N        O
#> * <int>     <int>          <int> <named lis> (err) [… (err) [… (err) [… (err) […
#> 1     1         1             23 <df>        0.47902… 0.05625… 0.00968… 0.39768…
#> # ℹ 5 more variables: S (err) [g/g], d15N <dbl>, d13C <dbl>,
#> #   calcium_content_1 (err) [g/g], calcium_content_1_in_pd <lgl>

# d13C values
irpeat::irp_d13C_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train"
)
#> Warning: 650 selected instead of 645.
#> # A tibble: 1 × 13
#>   id_90 sample_id measurement_id spectra            C        H        N        O
#> * <int>     <int>          <int> <named lis> (err) [… (err) [… (err) [… (err) […
#> 1     1         1             23 <df>        0.47902… 0.05625… 0.00968… 0.39768…
#> # ℹ 5 more variables: S (err) [g/g], d15N <dbl>, d13C <dbl>, d13C_1 (err) [1],
#> #   d13C_1_in_pd <lgl>

# d15N values
irpeat::irp_d15N_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train"
)
#> Warning: 650 selected instead of 645.
#> # A tibble: 1 × 13
#>   id_90 sample_id measurement_id spectra            C        H        N        O
#> * <int>     <int>          <int> <named lis> (err) [… (err) [… (err) [… (err) […
#> 1     1         1             23 <df>        0.47902… 0.05625… 0.00968… 0.39768…
#> # ℹ 5 more variables: S (err) [g/g], d15N <dbl>, d13C <dbl>, d15N_1 (err) [1],
#> #   d15N_1_in_pd <lgl>

# nominal oxidation state of carbon (NOSC)
irpeat::irp_nosc_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train"
)
#> Warning: 650 selected instead of 645.
#> # A tibble: 1 × 13
#>   id_90 sample_id measurement_id spectra            C        H        N        O
#> * <int>     <int>          <int> <named lis> (err) [… (err) [… (err) [… (err) […
#> 1     1         1             23 <df>        0.47902… 0.05625… 0.00968… 0.39768…
#> # ℹ 5 more variables: S (err) [g/g], d15N <dbl>, d13C <dbl>, nosc_1 (err) [1],
#> #   nosc_1_in_pd <lgl>

# Gibbs free energy of formation
irpeat::irp_dgf0_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train"
)
#> Warning: 650 selected instead of 645.
#> # A tibble: 1 × 13
#>   id_90 sample_id measurement_id spectra            C        H        N        O
#> * <int>     <int>          <int> <named lis> (err) [… (err) [… (err) [… (err) […
#> 1     1         1             23 <df>        0.47902… 0.05625… 0.00968… 0.39768…
#> # ℹ 5 more variables: S (err) [g/g], d15N <dbl>, d13C <dbl>,
#> #   dgf0_1 (err) [J/mol], dgf0_1_in_pd <lgl>

# dry bulk density
irpeat::irp_bulk_density_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train"
)
#> Warning: 650 selected instead of 645.
#> # A tibble: 1 × 13
#>   id_90 sample_id measurement_id spectra            C        H        N        O
#> * <int>     <int>          <int> <named lis> (err) [… (err) [… (err) [… (err) […
#> 1     1         1             23 <df>        0.47902… 0.05625… 0.00968… 0.39768…
#> # ℹ 5 more variables: S (err) [g/g], d15N <dbl>, d13C <dbl>,
#> #   bulk_density_1 (err) [g/cm^3], bulk_density_1_in_pd <lgl>

# loss on ignition
irp_loss_on_ignition_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train"
)
#> Warning: 650 selected instead of 645.
#> # A tibble: 1 × 13
#>   id_90 sample_id measurement_id spectra            C        H        N        O
#> * <int>     <int>          <int> <named lis> (err) [… (err) [… (err) [… (err) […
#> 1     1         1             23 <df>        0.47902… 0.05625… 0.00968… 0.39768…
#> # ℹ 5 more variables: S (err) [g/g], d15N <dbl>, d13C <dbl>,
#> #   loss_on_ignition_1 (err) [g/g], loss_on_ignition_1_in_pd <lgl>

# O/C
irpeat::irp_O_to_C_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train"
)
#> Warning: 650 selected instead of 645.
#> # A tibble: 1 × 13
#>   id_90 sample_id measurement_id spectra            C        H        N        O
#> * <int>     <int>          <int> <named lis> (err) [… (err) [… (err) [… (err) […
#> 1     1         1             23 <df>        0.47902… 0.05625… 0.00968… 0.39768…
#> # ℹ 5 more variables: S (err) [g/g], d15N <dbl>, d13C <dbl>,
#> #   O_to_C_1 (err) [g/g], O_to_C_1_in_pd <lgl>

# C/N
irpeat::irp_C_to_N_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train"
)
#> Warning: 650 selected instead of 645.
#> # A tibble: 1 × 13
#>   id_90 sample_id measurement_id spectra            C        H        N        O
#> * <int>     <int>          <int> <named lis> (err) [… (err) [… (err) [… (err) […
#> 1     1         1             23 <df>        0.47902… 0.05625… 0.00968… 0.39768…
#> # ℹ 5 more variables: S (err) [g/g], d15N <dbl>, d13C <dbl>,
#> #   C_to_N_1 (err) [g/g], C_to_N_1_in_pd <lgl>

# H/C
irpeat::irp_H_to_C_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train"
)
#> Warning: 650 selected instead of 645.
#> # A tibble: 1 × 13
#>   id_90 sample_id measurement_id spectra            C        H        N        O
#> * <int>     <int>          <int> <named lis> (err) [… (err) [… (err) [… (err) […
#> 1     1         1             23 <df>        0.47902… 0.05625… 0.00968… 0.39768…
#> # ℹ 5 more variables: S (err) [g/g], d15N <dbl>, d13C <dbl>,
#> #   H_to_C_1 (err) [g/g], H_to_C_1_in_pd <lgl>

# volume fraction of solids
irpeat::irp_volume_fraction_solids_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train"
)
#> Warning: 650 selected instead of 645.
#> # A tibble: 1 × 13
#>   id_90 sample_id measurement_id spectra            C        H        N        O
#>   <int>     <int>          <int> <named lis> (err) [… (err) [… (err) [… (err) […
#> 1     1         1             23 <df>        0.47902… 0.05625… 0.00968… 0.39768…
#> # ℹ 5 more variables: S (err) [g/g], d15N <dbl>, d13C <dbl>,
#> #   volume_fraction_solids_1 (err) [L/L], volume_fraction_solids_1_in_pd <lgl>

irpeat::irp_volume_fraction_solids_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train",
  bulk_density = 0.04
)
#> # A tibble: 1 × 13
#>   id_90 sample_id measurement_id spectra            C        H        N        O
#>   <int>     <int>          <int> <named lis> (err) [… (err) [… (err) [… (err) […
#> 1     1         1             23 <df>        0.47902… 0.05625… 0.00968… 0.39768…
#> # ℹ 5 more variables: S (err) [g/g], d15N <dbl>, d13C <dbl>,
#> #   volume_fraction_solids_1 (err) [L/L], volume_fraction_solids_1_in_pd <lgl>

# non-macroporosity
irpeat::irp_non_macroporosity_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train"
)
#> Warning: 650 selected instead of 645.
#> # A tibble: 1 × 13
#>   id_90 sample_id measurement_id spectra            C        H        N        O
#>   <int>     <int>          <int> <named lis> (err) [… (err) [… (err) [… (err) […
#> 1     1         1             23 <df>        0.47902… 0.05625… 0.00968… 0.39768…
#> # ℹ 5 more variables: S (err) [g/g], d15N <dbl>, d13C <dbl>,
#> #   non_macroporosity_1 (err) [L/L], non_macroporosity_1_in_pd <lgl>

irpeat::irp_non_macroporosity_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train",
  bulk_density = 0.04
)
#> # A tibble: 1 × 13
#>   id_90 sample_id measurement_id spectra            C        H        N        O
#>   <int>     <int>          <int> <named lis> (err) [… (err) [… (err) [… (err) […
#> 1     1         1             23 <df>        0.47902… 0.05625… 0.00968… 0.39768…
#> # ℹ 5 more variables: S (err) [g/g], d15N <dbl>, d13C <dbl>,
#> #   non_macroporosity_1 (err) [L/L], non_macroporosity_1_in_pd <lgl>

# macroporosity
irpeat::irp_macroporosity_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train"
)
#> Warning: 650 selected instead of 645.
#> # A tibble: 1 × 13
#>   id_90 sample_id measurement_id spectra            C        H        N        O
#>   <int>     <int>          <int> <named lis> (err) [… (err) [… (err) [… (err) […
#> 1     1         1             23 <df>        0.47902… 0.05625… 0.00968… 0.39768…
#> # ℹ 5 more variables: S (err) [g/g], d15N <dbl>, d13C <dbl>,
#> #   macroporosity_1 (err) [L/L], macroporosity_1_in_pd <lgl>

irpeat::irp_macroporosity_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train",
  bulk_density = 0.04
)
#> # A tibble: 1 × 13
#>   id_90 sample_id measurement_id spectra            C        H        N        O
#>   <int>     <int>          <int> <named lis> (err) [… (err) [… (err) [… (err) […
#> 1     1         1             23 <df>        0.47902… 0.05625… 0.00968… 0.39768…
#> # ℹ 5 more variables: S (err) [g/g], d15N <dbl>, d13C <dbl>,
#> #   macroporosity_1 (err) [L/L], macroporosity_1_in_pd <lgl>

# saturated hydraulic conductivity
irpeat::irp_saturated_hydraulic_conductivity_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train"
)
#> Warning: 650 selected instead of 645.
#> Error in inv_link(eta, link = slink): could not find function "inv_link"

irpeat::irp_saturated_hydraulic_conductivity_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train",
  bulk_density = 0.04
)
#> Error in inv_link(eta, link = slink): could not find function "inv_link"

# specific heat capacity
irpeat::irp_specific_heat_capacity_1(
  irpeat_sample_data[1, ],
  temperature = 290,
  do_summary = TRUE,
  check_prediction_domain = "train"
)
#> Warning: 650 selected instead of 645.
#> Error in inv_link(eta, link = slink): could not find function "inv_link"

irpeat::irp_specific_heat_capacity_1(
  irpeat_sample_data[1, ],
  temperature = 290,
  do_summary = TRUE,
  check_prediction_domain = "train",
  nitrogen_content = irpeat_sample_data$N[1]
)
#> Error in inv_link(eta, link = slink): could not find function "inv_link"

# dry thermal conductivity
irpeat::irp_dry_thermal_conductivity_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train"
)
#> Warning: 650 selected instead of 645.
#> Error in inv_link(eta, link = slink): could not find function "inv_link"

irpeat::irp_dry_thermal_conductivity_1(
  irpeat_sample_data[1, ],
  do_summary = TRUE,
  check_prediction_domain = "train",
  bulk_density = 0.04
)
#> Error in inv_link(eta, link = slink): could not find function "inv_link"

# microbial nitrogen content (note that the data are not ideal uses cases for
# the model (see the documentation). The following only demonstrates how to
# use the prediction function)
irpeat::irp_microbial_nitrogen_content_1(
  x = irpeat_sample_data[1, ],
  y = irpeat_sample_data[2, ],
  do_summary = TRUE,
  check_prediction_domain = "train"
)
#> # A tibble: 1 × 13
#>   id_90 sample_id measurement_id           C        H        N        O        S
#> * <int>     <int>          <int> (err) [g/g] (err) [… (err) [… (err) [… (err) […
#> 1     1         1             23 0.479025(0) 0.05625… 0.00968… 0.39768… 0.00395…
#> # ℹ 5 more variables: d15N <dbl>, d13C <dbl>,
#> #   microbial_nitrogen_content_1 (err) [g/g],
#> #   microbial_nitrogen_content_1_in_pd <lgl>, spectra <named list>

# degree_of_decomposition_1
if(! requireNamespace("posterior", quietly = TRUE)) {
x <-
  irpeat::irp_degree_of_decomposition_1(
    irpeat_sample_data[1, ],
    do_summary = TRUE,
    check_prediction_domain = "train",
    summary_function_sd = posterior::sd
  )
}

# degree_of_decomposition_2
if(! requireNamespace("posterior", quietly = TRUE)) {
x <-
  irpeat::irp_degree_of_decomposition_2(
    irpeat_sample_data[1, ],
    do_summary = TRUE,
    check_prediction_domain = "train",
    summary_function_sd = posterior::sd
  )
}

# degree_of_decomposition_3
if(! requireNamespace("posterior", quietly = TRUE)) {
x <-
  irpeat::irp_degree_of_decomposition_3(
    irpeat_sample_data[1, ],
    do_summary = TRUE,
    check_prediction_domain = "train",
    summary_function_sd = posterior::sd
  )
}
```
