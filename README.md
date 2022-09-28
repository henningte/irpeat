
<!-- README.md is generated from README.Rmd. Please edit that file -->

# irpeat <img src='man/figures/logo-hex.png' align="right" height="139" alt="logo" style="float:right; height:200px;" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

## Overview

‘irpeat’ is an R package that contains simple functions to analyze
infrared spectra of peat samples. Some functions may also work with
organic matter samples in general.

Provided functions for analyzing infrared spectra of peat are:

1.  Computation of several humification indices.
2.  Klason lignin mass fraction (following Hodgkins et al. (2018) and
    Teickner and Knorr (2022)) (note that these models are not reliable
    for peat, see Teickner and Knorr (2022)).
3.  Holocellulose mass fraction (following Hodgkins et al. (2018) and
    Teickner and Knorr (2022)) (note that these models are not reliable
    for peat, see Teickner and Knorr (2022)).
4.  Peat electron accepting capacity (following Teickner, Gao, and Knorr
    (2022)).
5.  Peat electron donating capacity (following Teickner, Gao, and Knorr
    (2022)).

### How to install

You can install ‘irpeat’ from GitHub using R via:

``` r
remotes::install_github(repo = "henningte/irpeat")
```

‘irpeat’ relies on the R package [‘ir’](https://github.com/henningte/ir)
for handling infrared spectra.

If you want to use the prediction models, you have to install the
[‘irpeatmodels’](---todo:%20add%20url) package in addition to the
‘irpeat’ package:

``` r
remotes::install_url("---todo:add url", type = "source")
```

### How to use

You can load ‘irpeat’ in R with:

``` r
library(irpeat)

# load additional packages needed for this tutorial
library(ir)
library(ggplot2)
```

You can test ‘irpeat’ with sample data from the R package ‘ir’:

``` r
ir::ir_sample_data
#> # A tibble: 58 × 7
#>    id_measurement id_sample sample_type sample_comment             klason_lignin
#>             <int> <chr>     <chr>       <chr>                      <units>      
#>  1              1 GN 11-389 needles     Abies Firma Momi fir       0.359944     
#>  2              2 GN 11-400 needles     Cupressocyparis leylandii… 0.339405     
#>  3              3 GN 11-407 needles     Juniperus chinensis Chine… 0.267552     
#>  4              4 GN 11-411 needles     Metasequoia glyptostroboi… 0.350016     
#>  5              5 GN 11-416 needles     Pinus strobus Torulosa     0.331100     
#>  6              6 GN 11-419 needles     Pseudolarix amabili Golde… 0.279360     
#>  7              7 GN 11-422 needles     Sequoia sempervirens Cali… 0.329672     
#>  8              8 GN 11-423 needles     Taxodium distichum Cascad… 0.356950     
#>  9              9 GN 11-428 needles     Thuja occidentalis Easter… 0.369360     
#> 10             10 GN 11-434 needles     Tsuga caroliniana Carolin… 0.289050     
#> # … with 48 more rows, and 2 more variables: holocellulose <units>,
#> #   spectra <named list>
```

`ir::ir_sample_data` contains various ATR-MIR spectra of organic
reference material (e.g. newspaper, wood, grass).

A simple workflow could be, for example, to baseline correct the spectra
(using functions of the package ‘ir’) compute various humification
indices and Klason lignin and holocellulose mass fractions in the
samples. We use only the first few spectra from `ir::ir_sample_data` to
speed the computations a bit up.

``` r
x <- 
  ir::ir_sample_data[1:10, ] %>%                                # data
  ir::ir_bc(method = "rubberband") %>%                          # baseline correction
  irpeat::irp_hi() %>%                                          # humification indices
  irpeat::irp_klason_lignin_2(do_summary = TRUE)                # Klason lignin content
```

`x` is identical to `ir::ir_sample_data[1:10, ]`, but contains
additional columns for the computed humification indices (`h1`, `h2`,
`h3`, `h4`) and the computed Klason lignin content (`klason_lignin_2`)

``` r
x
#> # A tibble: 10 × 12
#>    id_measurement id_sample sample_type sample_comment             klason_lignin
#>  *          <int> <chr>     <chr>       <chr>                                [1]
#>  1              1 GN 11-389 needles     Abies Firma Momi fir               0.360
#>  2              2 GN 11-400 needles     Cupressocyparis leylandii…         0.339
#>  3              3 GN 11-407 needles     Juniperus chinensis Chine…         0.268
#>  4              4 GN 11-411 needles     Metasequoia glyptostroboi…         0.350
#>  5              5 GN 11-416 needles     Pinus strobus Torulosa             0.331
#>  6              6 GN 11-419 needles     Pseudolarix amabili Golde…         0.279
#>  7              7 GN 11-422 needles     Sequoia sempervirens Cali…         0.330
#>  8              8 GN 11-423 needles     Taxodium distichum Cascad…         0.357
#>  9              9 GN 11-428 needles     Thuja occidentalis Easter…         0.369
#> 10             10 GN 11-434 needles     Tsuga caroliniana Carolin…         0.289
#> # … with 7 more variables: holocellulose [1], spectra <list>, hi1 <dbl>,
#> #   hi2 <dbl>, hi3 <dbl>, hi4 <dbl>, klason_lignin_2 (err) [g/g]
```

Plot of the humification index (ratio of the intensities at 1420 and
1090 cm<sup>-1</sup> (Broder et al. 2012)) versus the Klason lignin
content:

``` r
ggplot2::ggplot(x, aes(x = quantities::drop_quantities(klason_lignin_2) * 100, y = hi1)) + 
  ggplot2::geom_point() +
  ggplot2::labs(x = "Klason lignin content [mass-%]", 
                y = expression("Ratio of the intensities at"~1420~and~1090~cm^{-1}))
```

![](man/figures/README-x_plot-1.png)<!-- -->

All computed quantities come with units and standard errors (thanks to
the [quantities](https://github.com/r-quantities/quantities) package):

``` r
x$klason_lignin_2
#> Units: [g/g]
#> Errors: 0.05305209 0.03890806 0.03432737 0.03900097 0.03664285 ...
#>        V1        V2        V3        V4        V5        V6        V7        V8 
#> 0.3762468 0.3420424 0.2542360 0.3085094 0.2976808 0.2766546 0.3137988 0.3524349 
#>        V9       V10 
#> 0.3391504 0.2918670
```

### Future development

Henning Teickner plans, as part of his PhD project, to extensively
extent ‘irpeat’ by developing a set of calibration models that can
predict various peat physicochemical properties from mid infrared
spectra. These models should be finished by September 2024. Currently, a
data compendium ([pmird](https://henningte.github.io/pmird/index.html))
is developed to collect the data required for this task.

### How to cite

Please cite this R package as:

> Henning Teickner, Suzanne B. Hodgkins (2022). *irpeat: Functions to
> Analyze Mid Infrared Spectra of Peat Samples*. Accessed 2022-09-28.
> Online at <https://github.com/henningte/irpeat>.

### Licenses

**Text and figures :**
[CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/)

**Code :** See the [DESCRIPTION](DESCRIPTION) file

**Data :** [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
attribution requested in reuse. See the sources section for data sources
and how to give credit to the original author(s) and the source.

### Contributions

We welcome contributions from everyone. Before you get started, please
see our [contributor guidelines](CONTRIBUTING.md). Please note that this
project is released with a [Contributor Code of Conduct](CONDUCT.md). By
participating in this project you agree to abide by its terms.

### Sources

The data and prediction models for holocellulose and Klason lignin
(`irp_content_h_hodgkins_model`, `irp_content_kl_hodgkins_model`) are
derived from Hodgkins et al. (2018) and were restructured to match the
requirements of ir. The original article containing the data can be
downloaded from <https://www.nature.com/articles/s41467-018-06050-2> and
is distributed under the Creative Commons Attribution 4.0 International
License (<https://creativecommons.org/licenses/by/4.0/>). The data on
Klason lignin and holocellulose content was originally derived from De
La Cruz, Florentino B., Osborne, and Barlaz (2016).

Modified prediction models for holocellulose and Klason lignin
(`model_holocellulose_2`, `model_klason_lignin_2`) are derived from
Teickner and Knorr (2022).

Data and models for the electrochemical accepting and donating
capacities (EAC, EDC) of peat were derived from Teickner, Gao, and Knorr
(2022) and Teickner, Gao, and Knorr (2021)

This packages was developed in R (R version 4.2.0 (2022-04-22 ucrt)) (R
Core Team 2019) using functions from devtools (Wickham, Hester, and
Chang 2019), usethis (Wickham and Bryan 2019), rrtools (Marwick 2019)
and roxygen2 (Wickham et al. 2019).

### References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Broder.2012" class="csl-entry">

Broder, T., C. Blodau, H. Biester, and K. H. Knorr. 2012. “<span
class="nocase">Peat decomposition records in three pristine ombrotrophic
bogs in southern Patagonia</span>.” *Biogeosciences* 9 (4): 1479–91.
<https://doi.org/10.5194/bg-9-1479-2012>.

</div>

<div id="ref-LaCruz.2016" class="csl-entry">

De La Cruz, Florentino B., Jason Osborne, and Morton A. Barlaz. 2016.
“<span class="nocase">Determination of Sources of Organic Matter in
Solid Waste by Analysis of Phenolic Copper Oxide Oxidation Products of
Lignin</span>.” *Journal of Environmental Engineering* 142 (2):
04015076. <https://doi.org/10.1061/(ASCE)EE.1943-7870.0001038>.

</div>

<div id="ref-Hodgkins.2018" class="csl-entry">

Hodgkins, Suzanne B., Curtis J. Richardson, René Dommain, Hongjun Wang,
Paul H. Glaser, Brittany Verbeke, B. Rose Winkler, et al. 2018. “<span
class="nocase">Tropical peatland carbon storage linked to global
latitudinal trends in peat recalcitrance</span>.” *Nature
Communications* 9 (1): 3640.
<https://doi.org/10.1038/s41467-018-06050-2>.

</div>

<div id="ref-Marwick.2019" class="csl-entry">

Marwick, Ben. 2019. “<span class="nocase">rrtools: Creates a
Reproducible Research Compendium</span>.”
<https://github.com/benmarwick/rrtools>.

</div>

<div id="ref-RCoreTeam.2019" class="csl-entry">

R Core Team. 2019. “<span class="nocase">R: A Language and Environment
for Statistical Computing</span>.” Vienna, Austria: R Foundation for
Statistical Computing. <https://www.R-project.org/>.

</div>

<div id="ref-Teickner.2021c" class="csl-entry">

Teickner, Henning, Chuanyu Gao, and Klaus-Holger Knorr. 2021.
“Reproducible Research Compendium with R Code and Data for:
’Electrochemical Properties of Peat Particulate Organic Matter on a
Global Scale: Relation to Peat Chemistry and Degree of Decomposition’.”
Zenodo. <https://doi.org/10.5281/zenodo.5792970>.

</div>

<div id="ref-Teickner.2022" class="csl-entry">

———. 2022. “Electrochemical Properties of Peat Particulate Organic
Matter on a Global Scale: Relation to Peat Chemistry and Degree of
Decomposition.” *Global Biogeochemical Cycles*, February.
<https://doi.org/10.1029/2021GB007160>.

</div>

<div id="ref-Teickner.2022a" class="csl-entry">

Teickner, Henning, and Klaus-Holger Knorr. 2022. “Improving Models to
Predict Holocellulose and Klason Lignin Contents for Peat Soil Organic
Matter with Mid Infrared Spectra.” Preprint. Soil and methods.
<https://doi.org/10.5194/soil-2022-27>.

</div>

<div id="ref-Wickham.2019b" class="csl-entry">

Wickham, Hadley, and Jennifer Bryan. 2019. “<span
class="nocase">usethis: Automate Package and Project Setup</span>.”
<https://CRAN.R-project.org/package=usethis>.

</div>

<div id="ref-Wickham.2019c" class="csl-entry">

Wickham, Hadley, Peter Danenberg, Gábor Csárdi, and Manuel Eugster.
2019. “<span class="nocase">roxygen2: In-Line Documentation for
R</span>.” <https://CRAN.R-project.org/package=roxygen2>.

</div>

<div id="ref-Wickham.2019" class="csl-entry">

Wickham, Hadley, Jim Hester, and Winston Chang. 2019. “<span
class="nocase">devtools: Tools to Make Developing R Packages
Easier</span>.” <https://CRAN.R-project.org/package=devtools>.

</div>

</div>
