summarize_fit_hanks = function(fit, beta_loc, beta_ar) {
  
  s = summary(fit$fit)
  
  data.frame(
    param = rownames(s$coefficients),
    tstep = fit$tstep,
    s$coefficients, 
    lwr = s$coefficients[, 'Estimate'] - 1.96 * s$coefficients[, 'Std. Error'],
    upr = s$coefficients[, 'Estimate'] + 1.96 * s$coefficients[, 'Std. Error'],
    truth = c(beta_loc, beta_ar)
  ) %>% 
    dplyr::mutate(covered = (lwr <= truth) & (truth <= upr))
  
}
