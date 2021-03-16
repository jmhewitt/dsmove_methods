## Load your packages, e.g. library(drake).
source("./packages.R")

## Load your R files and subplans
invisible(lapply(list.files("./R", full.names = TRUE, recursive = TRUE), source))

source('R/plan.R')

# options(clustermq.scheduler = "multicore")

# make(the_plan, lock_envir = FALSE, targets = 'gibbs_fits_init_fits_impute_segments_sim_obs_0.25_sim_trajectory_sim_params_1_prior_params_0.125')

make(the_plan, lock_envir = FALSE, 
     targets = grep('sim_obs_empirical_durations', the_plan$target, value = TRUE))

make(the_plan, lock_envir = FALSE, parallelism = 'clustermq', jobs = 5)

# build graph components
graph = vis_drake_graph(the_plan, targets_only = TRUE)

# view graph
visNetwork::visHierarchicalLayout(graph, direction = "LR",
                                  edgeMinimization = FALSE)

# r_make()

build_times(starts_with('proposed_path'), digits = 8)
