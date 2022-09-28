# irpeat (development version)

* Performance improvement for `irp_preprocess()`.
* Updated `irp_preprocess()`, model preprocessing configurations, and prediction functions to consider the new arguments (1) `new_x_type` from `ir::ir_bin()`, and (2) `bd_do_impute` from `ir::ir_bc()`.
* Update the documentation: Adding explicit warnings to `irp_content_klh_hodgkins()`, `irp_content_kl_hodgkins_model`, and `irp_content_h_hodgkins_model` that the related models are not reliable for peat.

# irpeat 0.1.0

* Added a `NEWS.md` file to track changes to the package.
* Update the documentation (adding missing examples).
* Add pkgdown site.

#### Bug fixes

* Fix bugs caused by changes in 'ir' 0.2.0 (`ir::ir_check_ir()` no longer is exported, `measurement_id` is no longer a required column in `ir` objects).
* Fix bug in `irp_content()` if `variable = "klason_lignin_hodgkins"` (no values were computed).

#### New models

* New models to predict holocellulose (`model_holocellulose_2`) and Klason lignin (`model_klason_lingin_2`) contents.

#### Deprecated functions
