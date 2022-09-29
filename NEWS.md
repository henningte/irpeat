# irpeat (development version)

* Performance improvement for `irp_preprocess()`.
* Updated `irp_preprocess()`, model preprocessing configurations, and prediction functions to consider the new arguments (1) `new_x_type` from `ir::ir_bin()`, (2) `bc_do_impute` from `ir::ir_bc()`, (3) `do_return_as_ir` (allows to return the preprocessed spectra as `ir` object).
* Update the documentation: Adding explicit warnings to `irp_content_klh_hodgkins()`, `irp_content_kl_hodgkins_model`, and `irp_content_h_hodgkins_model` that the related models are not reliable for peat.
* Correct a typo which made it impossible to compute `holocellulose_2` with `irp_content()` (now: `irp_predict()`, see below).
* Major restructuration: The models originally included with the 'irpeat' package are now stored in a separate data package ['irpeatmodels'](---todo: add url) to reduce the size of the 'irpeat' package. The 'irpeatmodels' package can be installed from Zenodo.
* New models: `irp_carbon_content_1()`, `irp_nitrogen_content_1()`, `irp_oxygen_content_1()`, `irp_hydrogen_content_1()`, `irp_phosphorous_content_1()`, `irp_potassium_content_1()`, `irp_sulfur_content_1()`, `irp_titanium_content_1()`, `irp_d13C_1()`, `irp_d15N_1()`, `irp_nosc_1()`, `irp_dgf0_1()`, `irp_bulk_density_1()`, `irp_O_to_C_1()`, `irp_H_to_C_1()`, `irp_C_to_N_1()`.
* Add a new class `irp_prediction_domain`: This class stores the prediction domain for a model. Methods available are: Conversion from `ir` and `ir_flat` objects, plotting, check whether spectra in an `ir` object are within a prediction domain.
* Rename `irp_content()` to `irp_predict()`.

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
