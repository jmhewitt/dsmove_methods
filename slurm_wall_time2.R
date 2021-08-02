library(dplyr)
library(ggplot2)
library(ggthemes)

sim = function(njobs, ncores, time_to_submit, time_to_complete) {
  # Given njobs to run on at most ncores, simulate the times at which jobs
  # begin and end as well as the nodes that the jobs run on.  The simulation
  # accounts for possibly variable job run times (via time_to_complete function)
  # and the time required to send a job to a core to run (via time_to_submit).  
  # The time_to_submit parameter is important for compute clusters since 
  # this is a communication bottleneck and is unlikely to be small enough so as 
  # to have a negligible effect on the overall computation time and effort.
  #
  # All "time" values are unitless in the simulation, but intended to represent
  # seconds.
  #
  # Parameters:
  #   njobs - number of jobs to complete
  #   ncores - number of cores available for running jobs
  #   time_to_submit - (seconds) required to send a job to a worker
  #   time_to_complete - function returning (seconds) required to finish a job
  #
  # Return:
  #   data.frame with one line for each job, indicating the job number, the 
  #   time the job started and finished, and the node the job ran on.
  
  # total array of jobs to complete
  job_list = data.frame(
    jobid = 1:njobs,
    start_time = NA,
    end_time = NA,
    core = NA
  )
  
  # seed the job list
  job_list$start_time[1] = time_to_submit
  job_list$end_time[1] = complete_time() + job_list$start_time[1]
  job_list$core[1] = 1
  
  # process all remaining jobs
  for(i in 2:nrow(job_list)) {
    # list of active jobs after the previous job is submitted
    active_jobs = job_list$start_time[i-1] < job_list$end_time
    # number of jobs active after the previous job is submitted
    jobs_active = sum(active_jobs, na.rm = TRUE)
    # find next time available to submit a job
    if(jobs_active >= ncores) {
      # wait for a job to finish before submitting a new job
      t_avail_to_submit = min(job_list$end_time[active_jobs], na.rm = TRUE)
    } else {
      # submit a new job immediately
      t_avail_to_submit = job_list$start_time[i-1]
    }
    # number of jobs active when the current job is submitted
    active_jobs = t_avail_to_submit < job_list$end_time
    # job characteristics
    job_list$start_time[i] = t_avail_to_submit + time_to_submit
    job_list$end_time[i] = complete_time() + job_list$start_time[i]
    job_list$core[i] = min(setdiff(1:ncores, job_list$core[active_jobs]))
  }
  
  # time required to finish computations
  wall_time = max(job_list$end_time)
  
  # package results
  list(
    # data frame with details of completed jobs
    job_history = job_list,
    # time required to finish computations
    wall_time = wall_time,
    # most number of cores actually used for computation
    max_cores = max(job_list$core),
    # batch efficiency (i.e., active compute time / total CPU time)
    efficiency = sum(job_list$end_time - job_list$start_time) / 
      (wall_time * ncores)
  )
}

# timple for a single movement diffusion likelihood evaluaton
complete_time = function(){
  # runif(n = 1, min = 19125.94, max = 22444.54)  # fastloc whale data
  # runif(n = 1, min = 1100, max = 2000) # zc093 data
  runif(n = 1, min = 2300, max = 3600) # zc093 data with no gps truncation
}

# Monte Carlo simulations of runs when varying number of cores available
sims = lapply(seq(from = 100, to = 300, by = 10), function(ncores) {
  list(
    ncores = ncores,
    simulations = replicate(
      n = 1e2, 
      expr = list(sim(
        njobs = 1600, ncores = ncores, 
        time_to_submit = 60/3.7,
        time_to_complete = complete_time
      ))
    )
  )
})

# aggregate Monte Carlo output
df = do.call(rbind, lapply(sims, function(batch) {
  cbind(
    ncores = batch$ncores,
    do.call(rbind, lapply(batch$simulations, function(res) {
      data.frame(efficiency = res$efficiency, wall_time = res$wall_time,
                 max_cores = res$max_cores)
    }))
  )
}))

# view computation time against node use efficiency
ggplot(df, aes(x = wall_time/3600, y = efficiency, col = factor(ncores))) + 
  geom_point() + 
  xlab('Wall time (h)') +
  stat_smooth(se = FALSE, inherit.aes = FALSE, 
              mapping = aes(x = wall_time/3600, y = efficiency)) + 
  theme_few()

# view computation time against number of cores available
ggplot(df, aes(x = factor(ncores), y = wall_time/3600)) + 
  geom_boxplot() + 
  ylab('Wall time (h)') +
  theme_few()

# view node use efficiency against number of cores available
ggplot(df, aes(x = factor(ncores), y = efficiency)) + 
  geom_boxplot() + 
  theme_few()

# view number of cores that actually processed jobs vs 
# number of cores available
ggplot(df, aes(x = factor(ncores), y = max_cores)) + 
  geom_boxplot() + 
  theme_few()

# view one of the Monte Carlo simulations, where horizontal bars
# represent active computation time, and the y-axis value is the 
# core on which computation occurs
sim_batch = 10
ggplot(sims[[sim_batch]]$simulations[[1]]$job_history, 
       aes(xmin = start_time/3600, xmax = end_time/3600, y = core)) + 
  geom_linerange(alpha = .5) + 
  geom_hline(yintercept = sims[[sim_batch]]$ncores, lty = 3) + 
  ylim(0,sims[[sim_batch]]$ncores) + 
  theme_few() + 
  xlab('Time (h)') + 
  ggtitle(paste('MC sample when', sims[[sim_batch]]$ncores, 
                'nodes available for computation.'))
  