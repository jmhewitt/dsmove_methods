library(targets)
library(future)
library(future.batchtools)

# basic check for existence of SLURM job submission command to determine if 
# running on a SLURM-enabled server
if(system('command -v sbatch') == 0) {
  nnodes = as.numeric(Sys.getenv('SLURM_JOB_NUM_NODES'))
  if(is.na(nnodes)) {
    nnodes = 1
  }
  if(nnodes > 1) {
    plan(cluster, workers = snow::getMPIcluster())
  } else {
    plan(batchtools_slurm, template = file.path("hpc", "slurm_batchtools.tmpl"))
  }
} else {
  plan(multisession)
}

# set packages to load
tar_option_set(
  # packages = c('dplyr', 'lubridate', 'ggplot2', 'ggthemes', 'stringr', 
  #              'nimble', 'expm', 'pryr', 'suncalc', 'tarchetypes', 'coda',
  #              'tidyr', 'future', 'future.batchtools', 'viridis'),
  packages = c('dsmovetools', 'ctmcmove', 'coda', 'ggplot2', 'ggthemes', 
               'dplyr', 'bisque', 'sp', 'lubridate', 'ggnewscale', 'viridis',
               'metR', 'fields', 'dsmovetools2d'),
  imports = 'dsmovetools',
  deployment = 'main'
)

## load R files and workflows
lapply(list.files("R", full.names = TRUE, recursive = TRUE), source)

# assemble workflow
c(
  dir_targets,
  simulation_targets,
  exact_targets,
  whale_targets,
  fastloc_targets,
  whale_targets_multires
)