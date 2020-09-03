plot_shortest_path = function(ctds_struct, trajectory, obs, plot_dir) {
  
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
  
  pl.x = ggplot(imputation_coords, 
                aes(x = time, y = s1, col = Type)) + 
    geom_line() + 
    geom_line(data = sim_coords) + 
    scale_color_brewer('', type = 'qual', palette = 'Dark2') + 
    ggtitle(paste('Time between observations:', round(tstep, 2))) + 
    ylab('x') + 
    theme_few() + 
    theme(panel.border = element_blank(), 
          axis.title.y = element_text(angle = 0, vjust = 0.5), 
          axis.title.x = element_blank())
  
  pl.y = ggplot(imputation_coords, 
                aes(x = time, y = s2, col = Type)) + 
    geom_line() + 
    geom_line(data = sim_coords) + 
    scale_color_brewer('', type = 'qual', palette = 'Dark2') + 
    xlab('Time') + 
    ylab('y') + 
    theme_few() + 
    theme(panel.border = element_blank(), 
          axis.title.y = element_text(angle = 0, vjust = 0.5))
  
  pl = ggarrange(pl.x, pl.y, nrow = 2, ncol = 1, common.legend = TRUE,
                 legend = 'right')
  
  ggsave(pl, filename = file.path(plot_dir, 
                                  paste('shortest_path_imputations_vs_truth_', 
                                        round(tstep, 2), '.pdf', sep = ''))
         )
  
  pl
}
