summarize_fit_hanks = function(fit_output) {
  
  fit = fit_output$fit
  
  if(inherits(fit$fit, 'list')) {
    
    o = fit$fit
    
    se = sqrt(diag(solve(-o$hessian)))
    
    df = data.frame(
      param = c('(Intercept)', 'crw', 'shape'),
      tstep = fit$tstep,
      Estimate = o$par, 
      lwr = o$par - 1.96 * se,
      upr = o$par + 1.96 * se,
      truth = c(fit_output$params$beta_loc, fit_output$params$beta_ar, 
                fit_output$params$weibull_shape),
      weibull_est = TRUE,
      weibull_truth = fit_output$params$weibull_shape
    ) %>% 
      dplyr::mutate(covered = (lwr <= truth) & (truth <= upr))
    
  } else {
    
    s = summary(fit$fit)
    
    df = data.frame(
      param = c(rownames(s$coefficients)),
      tstep = fit$tstep,
      Estimate = s$coefficients[, 'Estimate'], 
      lwr = s$coefficients[, 'Estimate'] - 
        1.96 * s$coefficients[, 'Std. Error'],
      upr = s$coefficients[, 'Estimate'] + 
        1.96 * s$coefficients[, 'Std. Error'],
      truth = c(fit_output$params$beta_loc, fit_output$params$beta_ar),
      weibull_est = FALSE,
      weibull_truth = fit_output$params$weibull_shape
    ) %>% 
      dplyr::mutate(covered = (lwr <= truth) & (truth <= upr))
  }
  
  df
}
