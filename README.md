
<!-- README.md is generated from README.Rmd. Please edit that file -->

# irpeat <img src='man/figures/logo-hex.png' align="right" height="139" alt="logo" style="float:right; height:200px;" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

## Overview

‘irpeat’ is an R package which contains functions to analyze infrared
spectra of peat samples. These functions are functions to compute
humification indices and functions to predict peat properties. Some
functions may also work with organic matter samples in general.

The following peat properties can currently be predicted:

-   Elemental contents (C, H, N, O, S, P, K, Ti)
-   isotope values
    (![\delta^{13}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cdelta%5E%7B13%7D "\delta^{13}")C
    and
    ![\delta^{15}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cdelta%5E%7B15%7D "\delta^{15}")N)
-   physical properties (bulk density, volume fraction of solids,
    non-macroporosity, macroporosity saturaed hydraulic conductivity)
-   standard Gibbs free energy of formation
    (![\Delta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CDelta "\Delta")G![\_\text{f}^0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;_%5Ctext%7Bf%7D%5E0 "_\text{f}^0"))
-   electrochemical properties (electron accepting capacity, electron
    donating capacity)

The package also contains functions to predict holocellulose and Klason
lignin contents \[Hodgkins et al. (2018); Teickner.2022a\], but these
models have been shown to be biased (Teickner and Knorr 2022).

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

`irpeat_sample_data` contains transmission mid infrared spectra of peat
different samples (See Teickner, Gao, and Knorr (2022) and Teickner,
Gao, and Knorr (2021) for details).

A simple workflow could be, for example, to baseline correct the spectra
(using functions of the package ‘ir’) compute various humification
indices and Klason lignin and holocellulose mass fractions in the
samples. We use only the first few spectra from `ir::ir_sample_data` to
speed the computations a bit up.

``` r
x <- 
  irpeat_sample_data %>%                                  # data
  dplyr::mutate(
    hi_1630_1090 =
      irpeat_sample_data %>%
      ir::ir_bc(method = "rubberband") %>%                # baseline correction
      ir::ir_interpolate(start = NULL, dw = 1) %>%        # interpolation
      irp_hi(x1 = 1630, x2 = 1090) %>%                    # humification index
      dplyr::pull(hi_1630_1090)
  ) %>%
  irpeat::irp_carbon_content_1(do_summary = TRUE) %>%     # C content
  irpeat::irp_bulk_density_1(do_summary = TRUE)           # bulk density
```

`x` is identical to `ir::ir_sample_data`, but contains additional
columns for the computed humification indices (`h1`, `h2`, `h3`, `h4`)
and the computed carbon content (`carbon_content_1`)

``` r
x
#> # A tibble: 59 × 9
#>    id_90 sample_id measurement_id hi_1630_1090 carbon_content_1 carbon_content_…
#>  * <int>     <int>          <int> <numeric>         (err) [g/g] <logical>       
#>  1     1         1             23 0.6190412             0.43(2) FALSE           
#>  2     2         2             32 0.4393524             0.39(2) FALSE           
#>  3     3         3             38 0.5568726             0.43(2) FALSE           
#>  4     5         5             52 0.6164633             0.42(2) FALSE           
#>  5     6         6             54 0.9495150             0.46(2) FALSE           
#>  6     7         7             55 0.8387912             0.45(2) FALSE           
#>  7     8         8             56 0.7038940             0.42(2) FALSE           
#>  8     9         9             57 0.8111736             0.46(2) FALSE           
#>  9    10        10             24 0.5434465             0.38(2) FALSE           
#> 10    11        11             25 0.6056091             0.38(2) FALSE           
#> # … with 49 more rows, and 3 more variables: bulk_density_1 (err) [g/cm^3],
#> #   bulk_density_1_in_pd <logical>, spectra <named list>
```

Plot of the humification index (ratio of the intensities at 1420 and
1090 cm<sup>-1</sup> (Broder et al. 2012)) versus the Klason lignin
content:

``` r
ggplot2::ggplot(
  x, 
  aes(x = quantities::drop_quantities(carbon_content_1) * 100, y = hi_1630_1090)
) + 
  ggplot2::geom_point() +
  ggplot2::labs(
    x = "Carbon content [mass-%]", 
    y = expression("Ratio of the intensities at"~1630~and~1090~cm^{-1})
  )
```

![](man/figures/README-x_plot-1.png)<!-- -->

All computed quantities come with units and standard errors (thanks to
the [quantities](https://github.com/r-quantities/quantities) package):

``` r
x$carbon_content_1[1:5]
#> Units: [g/g]
#> Errors: 0.01611127 0.01714743 0.01647694 0.01669999 0.01584673
#>        V1        V2        V3        V4        V5 
#> 0.4317706 0.3906542 0.4286142 0.4227900 0.4623638
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
> Analyze Mid Infrared Spectra of Peat Samples*. Accessed 2022-09-29.
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
