# Checks whether spectra are within a prediction domain

Checks whether spectra are within a prediction domain

## Usage

``` r
irp_is_in_prediction_domain(x, prediction_domain)
```

## Arguments

- x:

  An object of class `ir`. The wavenumber values in the spectra of `x`
  need to be identical.

- prediction_domain:

  An object of class `irp_prediction_domain`. The wavenumber values in
  `prediction_domain` need to be identical to those in `x`.

## Value

`x` with an additional column `is_in_prediction_domain` with value
`TRUE` if a spectrum in `x` has intensity values which are within the
prediction domain and `FALSE` if not.
