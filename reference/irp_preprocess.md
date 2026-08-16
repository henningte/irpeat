# General purpose preprocessing of spectra for predictions

`irp_preprocess` is a function that provides a general-purpose
preprocessing workflow for spectra. The workflow comprises
interpolation, clipping, baseline correction, smoothing (including
optionally derivatization), normalization, binning, and scaling. All
these steps are optionally, but occur on a fixed order that cannot be
changed.

## Usage

``` r
irp_preprocess(
  x,
  do_interpolate = TRUE,
  interpolate_start = NULL,
  interpolate_dw = 1,
  do_clip = TRUE,
  clip_range,
  do_interpolate_region = FALSE,
  interpolate_region_range,
  do_bc = TRUE,
  bc_method = "rubberband",
  bc_degree = 2,
  bc_cutoff = 0,
  bc_do_impute = FALSE,
  do_smooth = TRUE,
  smooth_method = "sg",
  smooth_p = 3,
  smooth_n = smooth_p + 3 - smooth_p%%2,
  smooth_m = 1,
  smooth_ts = 0,
  smooth_k = 111,
  do_normalise = TRUE,
  normalise_method = "area",
  do_bin = TRUE,
  bin_width = 10,
  bin_new_x_type = "start",
  do_scale = TRUE,
  scale_center = TRUE,
  scale_scale = TRUE,
  do_return_as_ir = FALSE
)
```

## Arguments

- x:

  An object of class
  [`ir`](https://henningte.github.io/ir/reference/ir_new_ir.html) (it is
  assumed that `x` is not yet preprocessed).

- do_interpolate:

  A logical value indicating if spectra should be interpolated using
  [`ir::ir_interpolate()`](https://henningte.github.io/ir/reference/ir_interpolate.html).

- interpolate_start:

  See `ir_interpolate` (parameter `start`).

- interpolate_dw:

  See `ir_interpolate` (parameter `dw`).

- do_clip:

  A logical value indicating if spectra should be clipped using
  [`ir::ir_clip()`](https://henningte.github.io/ir/reference/ir_clip.html).

- clip_range:

  See `ir_clip` (parameter `range`).

- do_interpolate_region:

  A logical value indicating if spectra should be linearly interpolated
  in selected regions using
  [`ir::ir_interpolate_region()`](https://henningte.github.io/ir/reference/ir_interpolate_region.html).

- interpolate_region_range:

  See `ir_interpolate_region` (parameter `range`).

- do_bc:

  A logical value indicating if spectra should be baseline corrected
  using
  [`ir::ir_bc()`](https://henningte.github.io/ir/reference/ir_bc.html).

- bc_method:

  See `ir_bc` (parameter `method`).

- bc_degree:

  See `ir_bc` (parameter `degree`).

- bc_cutoff:

  A numeric value representing the wavenumber units to remove at the
  start and end of each spectrum in `x` during baseline correction. This
  may be done to remove artifacts due to baseline correction.

- bc_do_impute:

  See `ir_bc` (parameter `do_impute`).

- do_smooth:

  A logical value indicating if spectra should be smoothed using
  [`ir::ir_smooth()`](https://henningte.github.io/ir/reference/ir_smooth.html).

- smooth_method:

  See `ir_smooth` (parameter `method`).

- smooth_p:

  See `ir_smooth` (parameter `p`).

- smooth_n:

  See `ir_smooth` (parameter `n`).

- smooth_m:

  See `ir_smooth` (parameter `m`).

- smooth_ts:

  See `ir_smooth` (parameter `ts`).

- smooth_k:

  See `ir_smooth` (parameter `k`).

- do_normalise:

  A logical value indicating if spectra should be normalized using
  [`ir::ir_normalize()`](https://henningte.github.io/ir/reference/ir_normalize.html).

- normalise_method:

  See `ir_normalize` (parameter `method`).

- do_bin:

  A logical value indicating if spectra should be binned using
  [`ir::ir_bin()`](https://henningte.github.io/ir/reference/ir_bin.html).

- bin_width:

  See `ir_bin` (parameter `width`).

- bin_new_x_type:

  See `ir_bin` (parameter `new_x_type`).

- do_scale:

  A logical value indicating if spectral variables should be scaled
  using [`base::scale()`](https://rdrr.io/r/base/scale.html).

- scale_center:

  See `scale` (parameter `center`). To scale each spectral variable
  independently, provide a vector with length equal to the number of
  spectral variables returned after preprocessing.

- scale_scale:

  See `scale` (parameter `scale`). To scale each spectral variable
  independently, provide a vector with length equal to the number of
  spectral variables returned after preprocessing.

- do_return_as_ir:

  Logical value indicating whether the preprocessed spectra should be
  returned as
  [`ir`](https://henningte.github.io/ir/reference/ir_new_ir.html) object
  (`TRUE`) or as data frame with the same structure as
  [`ir_flat`](https://henningte.github.io/ir/reference/ir_new_ir_flat.html)
  object (`FALSE`).

## Value

A data frame with spectra in rows and a column for each spectral
variable after preprocessing.

## Examples

``` r
# get sample data
x <- ir::ir_sample_data[1:5, ]

# example preprocessing
res <-
  irpeat::irp_preprocess(
    x,
    do_interpolate = TRUE,
    interpolate_start = NULL,
    interpolate_dw = 1,
    do_bc = TRUE,
    do_clip = FALSE,
    do_interpolate_region = FALSE,
    bc_method = "rubberband",
    bc_cutoff = 10,
    bc_do_impute = FALSE,
    do_smooth = FALSE,
    do_normalise = TRUE,
    normalise_method = "area",
    do_bin = TRUE,
    bin_width = 10,
    bin_new_x_type = "start",
    do_scale = TRUE,
    scale_center = TRUE,
    scale_scale = TRUE,
    do_return_as_ir = FALSE
  )
```
