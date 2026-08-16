# Computes contents of Klason lignin and holocellulose of peat from mid-infrared spectra

`irp_content_klh_hodgkins` computes the mass fraction of Klason lignin
and holocellulose in peat from mid-infrared spectra of the peat samples.
This function may also work for organic matter in general. Note that the
models have been shown to be biased for peat samples (Teickner and Knorr
2022) .

## Usage

``` r
irp_content_klh_hodgkins(x, export = NULL, verbose = FALSE, make_plots = FALSE)
```

## Source

`irp_content_klh_hodgkins` is a wrapper function to a script written by
Suzanne Hodgkins
(<https://github.com/shodgkins/FTIRbaselines/blob/master/FTIRbaselines.R>).
and distributed under the GPL-3
(<https://www.gnu.org/licenses/gpl-3.0.html>) license. The original
script was modified for this purpose.

## Arguments

- x:

  An object of class
  [`ir`](https://henningte.github.io/ir/reference/ir_new_ir.html) with
  all spectra being in the range 650 to 4000 cm\$\$^{-1}\$\$.

- export:

  Either a valid path to an existing directory where to store the
  results of the original script (see section "Source")
  `irp_content_klh_hodgkins` is based on (if this output should be
  exported) or `NULL` (if nothing should be exported). Note that
  `export` has to be specified in order to analyze any issues with the
  data (e.g. high clay content) that may confound the computed values.
  Refer to the original publication ((Hodgkins et al. 2018) ) to get
  information on the meaning of the exported objects.

- verbose:

  A logical value indicating if messages should be printed
  (`verbose = TRUE`) or not (`verbose = FALSE`).

- make_plots:

  logical value indicating if plots should be printed
  (`make_plots = TRUE`) or not (`make_plots = FALSE`).

## Value

An object of class
[`ir`](https://henningte.github.io/ir/reference/ir_new_ir.html) with
additional columns containing the computed mass fractions of Klason
lignin and holocellulose for the spectra in `x`.

## Note

Note that the models have been shown to be biased for peat samples
(Teickner and Knorr 2022) .

## References

Hodgkins SB, Richardson CJ, Dommain R, Wang H, Glaser PH, Verbeke B,
Winkler BR, Cobb AR, Rich VI, Missilmani M, Flanagan N, Ho M, Hoyt AM,
Harvey CF, Vining SR, Hough MA, Moore TR, Richard PJH, De La Cruz FB,
Toufaily J, Hamdan R, Cooper WT, Chanton JP (2018). “Tropical Peatland
Carbon Storage Linked to Global Latitudinal Trends in Peat
Recalcitrance.” *Nature Communications*, **9**(1), 3640. ISSN 2041-1723.
[doi:10.1038/s41467-018-06050-2](https://doi.org/10.1038/s41467-018-06050-2)
.  
  
Teickner H, Knorr K (2022). “Improving Models to Predict Holocellulose
and Klason Lignin Contents for Peat Soil Organic Matter with
Mid-Infrared Spectra.” *SOIL*, **8**(2), 699–715.
[doi:10.5194/soil-8-699-2022](https://doi.org/10.5194/soil-8-699-2022) .

## Examples

``` r
library(ir)

irp_content_klh_hodgkins(
  ir::ir_sample_data[1:5, ],
  export = NULL,
  verbose = TRUE,
  make_plots = FALSE
)
#> Processing 1 (1 of 5)
#> Processing 2 (2 of 5)
#> Processing 3 (3 of 5)
#> Processing 4 (4 of 5)
#> Processing 5 (5 of 5)
#> Spectra with baselines and peak locations will be shown for all 5 samples.
#> Computing Klason lignin and holocellulose contents.
#> # A tibble: 5 × 9
#>   id_measurement id_sample sample_type sample_comment              klason_lignin
#> *          <int> <chr>     <chr>       <chr>                                 [1]
#> 1              1 GN 11-389 needles     Abies Firma Momi fir                0.360
#> 2              2 GN 11-400 needles     Cupressocyparis leylandii …         0.339
#> 3              3 GN 11-407 needles     Juniperus chinensis Chines…         0.268
#> 4              4 GN 11-411 needles     Metasequoia glyptostroboid…         0.350
#> 5              5 GN 11-416 needles     Pinus strobus Torulosa              0.331
#> # ℹ 4 more variables: holocellulose [1], spectra <named list>,
#> #   holocellulose_content_1 (err) [g/g], klason_lignin_content_1 (err) [g/g]
```
