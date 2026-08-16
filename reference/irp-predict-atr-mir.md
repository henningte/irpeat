# Predicts sample properties from attenuated total reflection mid-infrared spectra

Functions to predict peat properties from attenuated total reflection
(ATR) mid-infrared spectra. All functions below have been computed using
various organic matter samples. For detailed information on the
underlying prediction models, see the details section.

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
```

## Source

- `irp_holocellulose_2()`, `irp_klason_lignin_2()`:

  Teickner and Knorr (2022) .

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

  [`irp_microbial_nitrogen_content_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md)

  :   Here, `x` is a set of litter spectra after decomposition (and `y`
      is a set of litter spectra before decomposition). See the Details
      section.

- ...:

  Additional arguments passed to
  [`rstanarm::posterior_predict.stanreg()`](https://mc-stan.org/rstanarm/reference/posterior_predict.stanreg.html)
  ([`irp_eac_1()`](https://henningte.github.io/irpeat/reference/irp-predict-transmission-mir.md),`irp_eac_2()`).

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
  of draws, otherwise the result will be returned as `rvar` object.

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

## Details

The models use the models of the same name in the 'irpeatmodels'
package. The 'irpeatmodels' package provides information on the models.

## References

Teickner H, Knorr K (2022). “Improving Models to Predict Holocellulose
and Klason Lignin Contents for Peat Soil Organic Matter with
Mid-Infrared Spectra.” *SOIL*, **8**(2), 699–715.
[doi:10.5194/soil-8-699-2022](https://doi.org/10.5194/soil-8-699-2022) .

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
```
