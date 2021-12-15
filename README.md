# Reproducibility workflow for state space estimation of CTDS models

Code and workflows used to conduct analyses and create figures for the introduction of state space estimation of CTDS models.  The workflow is designed to be run using the [targets workflow manager](https://books.ropensci.org/targets/) for [R](https://cran.r-project.org), partly on the [Duke University Compute Cluster](https://rc.duke.edu/dcc/) via [SLURM](https://slurm.schedmd.com) cluster management software, and partly on personal computers.  The cluster configurations are difficult to make completely reproducible.


# Workflow

1. Install the `dsmovetools2d` package in `packages/dsmovetools2d`.  The package is primarily used for evaluating the likelihood and for providing additional utility functions; the package is not intended for distribution via CRAN.
2. Build targets that contain posterior surfaces and comparisons.  Targets can be built by using the `targets::tar_make` function within R.  Setting the `tar_make` option `callr_function = NULL` within an interactive R environment can let users explore code as it is running if `browser()` calls are added to the workflow scripts:
    - Posterior distribution targets for application
        - whale_ll_approx: evaluates likelihood at gridded parameter values
            - Save results to file `whale_ll_subsets.rds` as some targets require the output in this format for subsequent processing
        - whale_marginal_approx: evalutes marginal distributions for locations at gridded parameter values
        - whale_marginal_location_post, whale_marginal_additional_location_post, whale_marginal_additional_location_post2: build posterior distributions for a whale's location
    - Additional plots for application
        - zc095_post_loc_comparison, zc095_post_uncertainty_comparison: build plots that compare uncertainty in estimates of a whale's location under different scenarios
    - Simulation targets
        - simulation_results_combined_univariate, simulation_results_combined_bivariate: builds tables and figures that compare parameter estimates across computational methods
            - Running these targets should automatically run all required simulations and estimations.
3. Run [tmp_crawl.R](tmp_crawl.R) to generate AID-based estimates of a whale's location for comparison.
4. Run [distance_from_seafloor.R](distance_from_seafloor.R) to make and compare estimates of seafloor depth at specific times, and summarize predictions in a figure.

## Modifications

The function `targets::tar_make_future` can be used to build targets instead of using `targets::tar_make` if one has access to parallel computing resources.  See the targets manual chapter on [High performance computing (HPC)](https://books.ropensci.org/targets/hpc.html) for more details about supported features.  HPC workflows are difficult to make reproducible because system configurations
