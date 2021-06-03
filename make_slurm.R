library(targets)

tar_make_future(names = c(sim_fit_dtmc_gapprox, sim_fits_hanks), workers = 400)