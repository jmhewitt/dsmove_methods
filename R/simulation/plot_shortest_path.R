plot_shortest_path = function(ctds_struct, trajectory, obs, plot_dir, 
                              beta_loc) {
  
  # load simulated trajectory
  ctds_sim = readRDS(trajectory)
  
  # load observations
  ctds_obs = readRDS(obs)
  
  # extract timestep
  tstep = round(mean(diff(ctds_obs$times)), 2)
  
  # impute path from observation
  imputed = ctds.shortest_impute(states = ctds_obs$states, 
                                 times = ctds_obs$times, 
                                 ctds_struct = ctds_struct)
  
  # coordinates associated with simulated trajectory
  sim_coords = data.frame(ctds_struct$coords[ctds_sim$states, ], 
                          time = ctds_sim$times,
                          Type = 'Truth', rep = 1)
  
  # coordinates associated with imputed trajectory
  imputation_coords = data.frame(
    ctds_struct$coords[imputed$states[1:length(imputed$times)], ],
    time = imputed$times, 
    Type = 'Imputation'
  )
  
  # coordinates associated with observed trajectory
  observation_coords = data.frame(
    ctds_struct$coords[ctds_obs$states, ],
    time = ctds_obs$times,
    Type = 'Observation'
  )
  
  # compute empirical speed between observations or transitions
  empirical_speed = function(dat) {
    data.frame(do.call(rbind, lapply(2:nrow(dat), function(i) {
      ds = sqrt(sum(apply(dat[c(i, i-1), c('s1', 's2')], 2, diff)^2))
      dt = diff(dat$time[c(i-1, i)])
      c(speed = ds/dt, time = dat$time[i])
    })))
  }
  
  # empirical speed associated with trajectories
  sim_speed = empirical_speed(sim_coords)
  obs_speed = empirical_speed(observation_coords)
  imputation_speed = empirical_speed(imputation_coords)
  
  pl.x = ggplot(imputation_coords, 
                aes(x = time, y = s1, col = Type)) + 
    geom_step() + 
    geom_step(data = sim_coords) + 
    # geom_point(data = observation_coords, alpha = .5, col = 'black') + 
    scale_color_brewer('', type = 'qual', palette = 'Dark2') + 
    ggtitle(paste('Time between observations:', round(tstep, 2))) + 
    ylab('x') + 
    theme_few() + 
    theme(panel.border = element_blank(), 
          axis.title.y = element_text(angle = 0, vjust = 0.5), 
          axis.title.x = element_blank())
  
  pl.y = ggplot(imputation_coords, 
                aes(x = time, y = s2, col = Type)) + 
    geom_step() + 
    geom_step(data = sim_coords) + 
    # geom_point(data = observation_coords, alpha = .5, col = 'black') + 
    scale_color_brewer('', type = 'qual', palette = 'Dark2') + 
    xlab('Time') + 
    ylab('y') + 
    theme_few() + 
    theme(panel.border = element_blank(), 
          axis.title.y = element_text(angle = 0, vjust = 0.5))
  
  df = rbind(
    cbind(imputation_speed, Type = 'Imputation'),
    cbind(obs_speed, Type = 'Observed'),
    cbind(sim_speed, Type = 'Simulation')
  )
  
  pl.dsdt = ggplot(df, aes(x = time, y = speed, col = Type)) + 
    geom_line(alpha = .2) + 
    geom_point() + 
    geom_point(data = df %>% dplyr::filter(is.infinite(speed))) + 
    geom_hline(yintercept = exp(beta_loc), lty = 2) + 
    scale_color_brewer(type = 'qual', palette = 'Dark2') + 
    scale_y_log10() + 
    xlab('Time') + 
    ylab('Speed') + 
    ggtitle('Empirical speeds (Sim. parameter at dotted line)') + 
    theme_few() + 
    theme(panel.border = element_blank())
  
  df = rbind(df,
             data.frame(Type = 'Reference',
                        time = 1,
                        speed = 1/rexp(n = 1e4, rate = exp(beta_loc))))
  
  pl.dsdt_summary = ggplot(df, aes(x = speed, col = Type)) + 
    stat_density(geom = 'line', trim = TRUE) + 
    geom_vline(xintercept = exp(beta_loc), lty = 2) + 
    geom_vline(data = df %>% 
                 dplyr::group_by(Type) %>%
                 dplyr::summarise(speed = mean(speed)), 
               mapping = aes(xintercept = speed, col = Type)) + 
    scale_x_log10() +
    xlab('Speed') + 
    ylab('Density') + 
    ggtitle('Distribution of empirical speeds (Sim. parameter at dotted line)', 
            subtitle = '(Means at solid lines)') + 
    theme_few() + 
    theme(panel.border = element_blank())
  
  pl = ggarrange(pl.x, pl.y, nrow = 2, ncol = 1, common.legend = TRUE,
                 legend = 'right')
  
  ggsave(pl, filename = file.path(plot_dir, 
                                  paste('shortest_path_imputations_vs_truth_', 
                                        round(tstep, 2), '.pdf', sep = ''))
         )
  
  ggsave(pl.dsdt, filename = file.path(plot_dir, 
                                  paste('shortest_path_speeds_vs_truth_', 
                                        round(tstep, 2), '.pdf', sep = ''))
  )
  
  ggsave(pl.dsdt_summary, filename = file.path(plot_dir, 
                                       paste('shortest_path_speed_dists', 
                                             round(tstep, 2), '.pdf', sep = ''))
  )
  
  pl
}
