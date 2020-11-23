summarize_fit_hanks = function(fit_output) {
  
  fit = fit_output$fit

  o = fit$fit
  
  se = sqrt(diag(solve(-o$hessian)))
  
  df = data.frame(
    param = c('(Intercept)', 'crw'),
    tstep = fit$tstep,
    Estimate = o$par, 
    lwr = o$par - 1.96 * se,
    upr = o$par + 1.96 * se,
    truth = c(fit_output$params$beta_loc, fit_output$params$beta_ar)
  ) %>% 
    dplyr::mutate(covered = (lwr <= truth) & (truth <= upr))
  
  df
}
