# Creates an object of class `irp_prediction_domain`

Creates an object of class `irp_prediction_domain`

## Usage

``` r
new_irp_prediction_domain(x)
```

## Arguments

- x:

  A data frame with a row for each wavenumber value in the training data
  and the following columns:

  `x`

  :   A numeric value representing the wavenumber value \[cm\$^-1\$\].

  `ymin`

  :   A numeric value representing the minimum predictor variable value
      in the training data at this wavenumber value.

  `ymax`

  :   A numeric value representing the maximum predictor variable value
      in the training data at this wavenumber value.

  `x` may contain additional columns.

## Value

An object of class `irp_prediction_domain`. This is the same as `x`, but
with an additional class label.
