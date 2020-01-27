## code to prepare `irp_content_h_hodgkins_model` dataset goes here

devtools::load_all()

x <- ir::ir_sample_data

x_flat <- irp_content_klh_hodgkins_prepare(x)

res <- irp_content_klh_hodgkins_main(data = x_flat,
                                     export = export,
                                     verbose = verbose,
                                     make_plots = make_plots)

lm_x <- res$norm.Acorr$carb[x$sample_type != "old magazines"]
lm_y <- x$holocellulose[x$sample_type != "old magazines"]

irp_content_h_hodgkins_model <- lm(lm_y~lm_x)

usethis::use_data(irp_content_h_hodgkins_model, overwrite = TRUE)
