# Stores errors and units of error or units columns of objects of class `ir` in separate columns

**\[deprecated\]**

`irp_export_prepare` stores errors
([errors::errors](https://r-quantities.github.io/errors/reference/errors.html))
and units
([`units`](https://r-quantities.github.io/units/reference/units.html))
of error or units columns in an object of class `ir` in separate
columns.

## Usage

``` r
irp_export_prepare(x)
```

## Arguments

- x:

  An object of class
  [`ir`](https://henningte.github.io/ir/reference/ir_new_ir.html).

## Value

An object of class `ir` with a new column for each column that is of
class `units` and each column that is of class `errors`. New columns are
named `old_column_name_errors` and `old_column_name_units`,
respectively.

## Examples

``` r
library(ir)

irp_export_prepare(ir::ir_sample_data[1:5, ]) # now deprecated
#> Warning: `irp_export_prepare()` was deprecated in irpeat 0.1.0.
#> ℹ Please use `ir::ir_export_prepare()` instead.
#> # A tibble: 5 × 9
#>   id_measurement id_sample sample_type sample_comment              klason_lignin
#>            <int> <chr>     <chr>       <chr>                                 [1]
#> 1              1 GN 11-389 needles     Abies Firma Momi fir                0.360
#> 2              2 GN 11-400 needles     Cupressocyparis leylandii …         0.339
#> 3              3 GN 11-407 needles     Juniperus chinensis Chines…         0.268
#> 4              4 GN 11-411 needles     Metasequoia glyptostroboid…         0.350
#> 5              5 GN 11-416 needles     Pinus strobus Torulosa              0.331
#> # ℹ 4 more variables: holocellulose [1], spectra <named list>,
#> #   klason_lignin_units <chr>, holocellulose_units <chr>

# replacement for deprecated function:
if (FALSE) { # \dontrun{
ir::ir_export_prepare(ir::ir_sample_data[1:5, ], what = "metadata")
} # }
```
