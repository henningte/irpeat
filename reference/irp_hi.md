# Compute humification indices from mid-infrared spectra

`irp_hi` computes either custom humification indices or a predefined set
of humification indices as reported by Broder et al. (2012) for
mid-infrared spectra with the x values representing wavenumber values
(no checks are performed) (see the details section for the humification
indices computed by default). A humification index is the ratio of the
intensity values at two different x axis values (e.g. wavenumbers) and
defined as: \$\$ HI = \frac{y_x1}{y_x2} \$\$, where \\y_x1\\ is the
intensity of the spectrum at x axis value \\x1\\ and \\y_x2\\ is the
intensity of the spectrum at x axis value \\x2\\.

## Usage

``` r
irp_hi(x, x1 = NULL, x2 = NULL)
```

## Arguments

- x:

  An object of class
  [`ir`](https://henningte.github.io/ir/reference/ir_new_ir.html).

- x1:

  A numeric vector with values representing x axis values in the spectra
  of `x`. This is \\x1\\ in the equation displayed above. If multiple
  humification indices for the same `x2` should be computed, `x1` can
  contain multiple values. If `x1 = NULL`, the default humification
  indices will be computed (see details section).

- x2:

  A numeric value representing an x axis value in the spectra of `x`.
  This is \\x2\\ in the equation displayed above. If `x2 = NULL`, the
  default humification indices will be computed (see details section).

## Value

An object of class
[`ir`](https://henningte.github.io/ir/reference/ir_new_ir.html) with
additional columns for additional humification indices. If `x1 = NULL`
and `x2 = NULL`, these are four new columns (hi1, hi2, hi3 and hi4) that
correspond to the humification indices defined in the details section.
If `x1` and `x2` are not `NULL`, the columns have names `"hi_x1_x2"`
where `x1` and `x2` are replaced by the respective values.

## Details

The following humification indices are computed by default (if
`x1 = NULL` and `x2 = NULL`) (values represent wavenumbers
\[cm\\^{-1}\\\]) (Broder et al. 2012) :

- hi1: 1420/1090:

  OH and CO of phenols or CH of CH\\\_1\\ and CH\\\_3\\ groups (phenolic
  and aliphatic structures)/polysaccharides

- hi2: 1510/1090:

  Aromatic C=C or C=O of amides/polysaccharides

- hi3: 1630/1090:

  Aromatic C=C and COO\\^{-1}\\ (aromatics and aromatic or aliphatic
  carboxylates)/polysaccharides

- hi4: 1720/1090:

  carbonylic and carboxylic C=O (carboxylic acids and aromatic
  esters)/polysaccharides

## References

Broder T, Blodau C, Biester H, Knorr KH (2012). “Peat Decomposition
Records in Three Pristine Ombrotrophic Bogs in Southern Patagonia.”
*Biogeosciences*, **9**(4), 1479–1491. ISSN 1726-4189.
[doi:10.5194/bg-9-1479-2012](https://doi.org/10.5194/bg-9-1479-2012) .
[2021-06-02](https://henningte.github.io/irpeat/reference/2021-06-02).

## Examples

``` r
# get sample data
library(ir)

# compute default humification indices
d <- ir_sample_data[1:5, ]
d <- irp_hi(d)

# compute custom humification index
d <- ir_sample_data
d <- irp_hi(d, x1 = 2900, x2 = 1090)

# compute custom humification indices
d <- ir_sample_data
d <- irp_hi(d, x1 = c(2900, 1630), x2 = 1090)
```
