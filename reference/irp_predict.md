# Batch-prediction of sample properties

Wrapper function to batch-predict sample properties.

## Usage

``` r
irp_predict(x, y = NULL, variable, ...)

irp_content(x, variable, ...)
```

## Arguments

- x:

  An object of class
  [`ir`](https://henningte.github.io/ir/reference/ir_new_ir.html). See
  the individual prediction models for further data requirements.

- y:

  An object of class
  [`ir`](https://henningte.github.io/ir/reference/ir_new_ir.html) or
  `NULL`. Needed for prediction functions using more than one set of
  spectra (`irp_microbial_nitrogen_content_1`). See the individual
  prediction models for further data requirements. For other functions,
  `y` is not used and should be set to `NULL` (the default).

- variable:

  A character vector with one or more values that define for which
  components contents are computed for the spectra in `x`. Currently
  allowed values are:

  "all"

  :   `irp_content` computes all of the values below.

  "klason_lignin_content_1"

  :   Klason lignin mass fraction \[g/g\] as computed by
      [`irp_content_klh_hodgkins()`](https://henningte.github.io/irpeat/reference/irp_content_klh_hodgkins.md).

  "holocellulose_content_1"

  :   Holocellulose mass fraction \[g/g\] as computed by
      [`irp_content_klh_hodgkins()`](https://henningte.github.io/irpeat/reference/irp_content_klh_hodgkins.md).

  "klason_lignin_content_2"

  :   Klason lignin mass fraction \[g/g\] as computed by
      [`irp_klason_lignin_content_2()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "holocellulose_content_2"

  :   Holocellulose mass fraction \[g/g\] as computed by
      [`irp_holocellulose_content_2()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "eac_1"

  :   Electron accepting capacity as computed by
      [`irp_eac_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "edc_1"

  :   Electron donating capacity as computed by
      [`irp_edc_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "carbon_content_1"

  :   Carbon content as computed by
      [`irp_carbon_content_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "nitrogen_content_1"

  :   Nitrogen content as computed by
      [`irp_nitrogen_content_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "hydrogen_content_1"

  :   Hydrogen content as computed by
      [`irp_hydrogen_content_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "oxygen_content_1"

  :   Oxygen content as computed by
      [`irp_oxygen_content_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "phosphorus_content_1"

  :   Phosphorus content as computed by
      [`irp_phosphorus_content_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "potassium_content_1"

  :   Potassium content as computed by
      [`irp_potassium_content_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "sulfur_content_1"

  :   Sulfur content as computed by
      [`irp_sulfur_content_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "titanium_content_1"

  :   Titanium content as computed by
      [`irp_titanium_content_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "silicon_content_1"

  :   Silicon content as computed by
      [`irp_silicon_content_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "calcium_content_1"

  :   Calcium content as computed by
      [`irp_calcium_content_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "d13C_1"

  :   \\\delta^{13}\\C values as computed by
      [`irp_d13C_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "d15N_1"

  :   \\\delta^{15}\\N values as computed by
      [`irp_d15N_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "nosc_1"

  :   The nominal oxidation state of carbon as computed by
      [`irp_nosc_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "dgf0_1"

  :   The standard Gibbs free energy of formation content as computed by
      [`irp_dgf0_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "bulk_density_1"

  :   Bulk density as computed by
      [`irp_bulk_density_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "loss_on_ignition_1"

  :   Loss on ignition as computed by
      [`irp_loss_on_ignition_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "O_to_C_1"

  :   O/C ratio as computed by
      [`irp_O_to_C_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "H_to_C_1"

  :   H/C ratio as computed by
      [`irp_H_to_C_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "C_to_N_1"

  :   C/N ratio as computed by
      [`irp_C_to_N_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "volume_fraction_solids_1"

  :   Volume fraction of solids as computed by
      [`irp_volume_fraction_solids_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "non_macroporosity_1"

  :   Non-macroporosity as computed by
      [`irp_non_macroporosity_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "macroporosity_1"

  :   Macroporosity as computed by
      [`irp_macroporosity_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "saturated_hydraulic_conductivity_1"

  :   Saturated hydraulic conductivity as computed by
      [`irp_saturated_hydraulic_conductivity_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "specific_heat_capacity_1"

  :   Specific heat capacity as computed by
      [`irp_specific_heat_capacity_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "dry_thermal_conductivity_1"

  :   Dry thermal conductivity as computed by
      [`irp_dry_thermal_conductivity_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "microbial_nitrogen_content_1"

  :   Microbial nitrogen content as computed by
      [`irp_microbial_nitrogen_content_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "degree_of_decomposition_1"

  :   Degree of decomposition as computed by
      [`irp_degree_of_decomposition_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "degree_of_decomposition_2"

  :   Degree of decomposition as computed by
      [`irp_degree_of_decomposition_2()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

  "degree_of_decomposition_3"

  :   Degree of decomposition as computed by
      [`irp_degree_of_decomposition_3()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md).

- ...:

  Further arguments passed to individual prediction functions.

## Value

An object of class
[`ir`](https://henningte.github.io/ir/reference/ir_new_ir.html) with
additional columns containing the predictions for the spectra in `x`.

## Note

- `value = "klason_lignin_content_1"` and
  `value = "holocellulose_content_1"`:

  No warnings are shown and no values are exported to disk.

## Examples

``` r
library(ir)

irp_predict(
  ir::ir_sample_data[1, ],
  variable = "carbon_content_1",
  do_summary = TRUE
)
#> Warning: 650 selected instead of 645.
#> # A tibble: 1 × 9
#>   id_measurement id_sample sample_type sample_comment       klason_lignin
#> *          <int> <chr>     <chr>       <chr>                          [1]
#> 1              1 GN 11-389 needles     Abies Firma Momi fir         0.360
#> # ℹ 4 more variables: holocellulose [1], spectra <named list>,
#> #   carbon_content_1 (err) [g/g], carbon_content_1_in_pd <lgl>

if (FALSE) { # \dontrun{
irp_predict(
  ir::ir_sample_data[1, ],
  variable = c("eac_1", "carbon_content_1", "nitrogen_content_1", "dgf0_1"),
  do_summary = TRUE
)
} # }
```
