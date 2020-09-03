plot_imputations = function(ctds_struct, trajectory, imputations, 
                            raster.coords, tstep, plot_dir) {

  n_imputations = length(imputations)
  
  # load simulated trajectory
  ctds_sim = readRDS(trajectory)
  
  # coordinates associated with simulated trajectory
  sim_coords = data.frame(ctds_struct$coords[ctds_sim$states, ], 
                          time = ctds_sim$times,
                          Type = 'Truth', rep = 1)
  
  # coordinates associated with imputed trajectories
  imputation_coords = do.call(rbind, lapply(1:n_imputations, function(i) {
    imputed = imputations[[i]]
    data.frame(raster.coords[imputed$ec, ],
               time = imputed$trans.times,
               Type = 'Imputation', rep = i)
  }))
  colnames(imputation_coords)[1:2] = c('s1', 's2')
  
  pl.x = ggplot(imputation_coords, 
                aes(x = time, y = s1, col = Type, group = rep)) + 
    geom_line(alpha = .3) + 
    geom_line(data = sim_coords) + 
    scale_color_brewer('', type = 'qual', palette = 'Dark2') + 
    ggtitle(paste('Time between observations:', round(tstep, 2))) + 
    ylab('x') + 
    theme_few() + 
    theme(panel.border = element_blank(), 
          axis.title.y = element_text(angle = 0, vjust = 0.5), 
          axis.title.x = element_blank())
  
  pl.y = ggplot(imputation_coords, 
                aes(x = time, y = s2, col = Type, group = rep)) + 
    geom_line(alpha = .3) + 
    geom_line(data = sim_coords) + 
    scale_color_brewer('', type = 'qual', palette = 'Dark2') + 
    xlab('Time') + 
    ylab('y') + 
    theme_few() + 
    theme(panel.border = element_blank(), 
          axis.title.y = element_text(angle = 0, vjust = 0.5))
  
  pl = ggarrange(pl.x, pl.y, nrow = 2, ncol = 1, common.legend = TRUE,
                 legend = 'right')
  
  ggsave(pl, filename = file.path(plot_dir, paste('imputations_vs_truth_', 
                                                  round(tstep, 2),
                                                  '.pdf', sep = ''))
         )
  
  pl
}
