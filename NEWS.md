# irpeat (development version)

#### Major changes

* Major restructuration: The models originally included with the 'irpeat' package are now stored in a separate data package ['irpeatmodels'](---todo: add url) to reduce the size of the 'irpeat' package. The 'irpeatmodels' package can be installed from Zenodo.
* Add a new class `irp_prediction_domain`: This class stores the prediction domain for a model. Methods available are: Conversion from `ir` and `ir_flat` objects, plotting, check whether spectra in an `ir` object are within a prediction domain.
* Rename:
    1. `irp_content()` to `irp_predict()`.
    2. `irp_holocellulose_2()` to `irp_holocellulose_content_2()`.
    3. `irp_klason_lignin_2()` to `irp_klason_lignin_content_2()`.
    4. Output from `irp_hkl_hodgkins()`: `holocellulose_hodgkins` --> `holocellulose_content_1`, `klason_lignin_hodgkins` --> `klason_lignin_content_1`.
    5. Output from `irp_eac_1()`: `eac` --> `eac_1`.
    6. Output from `irp_edc_1()`: `edc` --> `edc_1`.
* For existing and new models: Training and testing prediction domains were added and can now be used to check whether the models cannot be used to make reliable predictions.

#### New functions

* Add `irp_preprocess_for()`: Allows to extract the spectra after automated preprocessing as they would be used to make predictions.
* Add `irp_get_prediction_domain_for()`: Allows to extract the prediction domain for a specific model.

#### New data

* Add sample data (transmission mid infrared spectra for peat samples from the 'redoxpeat' R package): `irpeat_sample_data`.

#### New models

* New models: `irp_carbon_content_1()`, `irp_nitrogen_content_1()`, `irp_oxygen_content_1()`, `irp_hydrogen_content_1()`, `irp_phosphorus_content_1()`, `irp_potassium_content_1()`, `irp_sulfur_content_1()`, `irp_titanium_content_1()`, `irp_d13C_1()`, `irp_d15N_1()`, `irp_nosc_1()`, `irp_dgf0_1()`, `irp_bulk_density_1()`, `irp_O_to_C_1()`, `irp_H_to_C_1()`, `irp_C_to_N_1()`, `irp_volume_fraction_solids_1()`, `irp_non_macroporosity_1()`, `irp_macroporosity_1()`, `irp_saturated_hydraulic_conductivity_1()`, `irp_specific_heat_capacity_1()`, `irp_dry_thermal_conductivity_1()`, `irp_microbial_nitrogen_content_1()`.


#### Improvements

* Performance improvement for `irp_preprocess()`.
* Updated `irp_preprocess()`, model preprocessing configurations, and prediction functions to consider the new arguments (1) `new_x_type` from `ir::ir_bin()`, (2) `bc_do_impute` from `ir::ir_bc()`, (3) `do_return_as_ir` (allows to return the preprocessed spectra as `ir` object).
* Update the documentation: Adding explicit warnings to `irp_content_klh_hodgkins()`, `irp_content_kl_hodgkins_model`, and `irp_content_h_hodgkins_model` that the related models are not reliable for peat.
* Correct a typo which made it impossible to compute `holocellulose_2` with `irp_content()` (now: `irp_predict()`, see below).
* `irp_content_klh_hodgkins_predict()`: The function now checks that spectra have the correct wavenumber range.
* `irp_eac_1()`, `irp_edc_1()`: Both functions now also accept `ir` objects with empty spectra as arguments.


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
