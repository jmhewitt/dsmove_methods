dir_targets = list(
  
  # location to store simulated data
  tar_target(
    name = sim_dir,
    command = {
      f = file.path('output', 'simulation')
      dir.create(f, recursive = TRUE, showWarnings = FALSE)
      f
    }
  ),
  
  # location for storing mcmc samples
  tar_target(
    name = mcmc_sample_dir, 
    command = {
      f = file.path('output', 'mcmc')
      dir.create(f, recursive = TRUE, showWarnings = FALSE)
      f
    }
  )
  
)
