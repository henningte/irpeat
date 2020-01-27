## code to prepare `irp_content_kl_hodgkins_model` dataset goes here

devtools::load_all()

x <- ir::ir_sample_data

x_flat <- irp_content_klh_hodgkins_prepare(x)

res <- irp_content_klh_hodgkins_main(data = x_flat,
                                     export = export,
                                     verbose = verbose,
                                     make_plots = make_plots)

lm_x <- res$norm.Acorr$arom15[x$sample_type != "office paper"] + res$norm.Acorr$arom16[x$sample_type != "office paper"]
lm_y <- x$klason_lignin[x$sample_type != "office paper"]

irp_content_kl_hodgkins_model <- lm(lm_y~lm_x)

usethis::use_data(irp_content_kl_hodgkins_model, overwrite = TRUE)
