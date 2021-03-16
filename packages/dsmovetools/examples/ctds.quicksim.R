beta_loc = 1

x = ctds.quicksim(dims = c(100, 100), beta_loc = beta_loc, beta_ar = 1, 
                  v0 = c(50,50), t0 = 0, tf = 1e2, max.steps = 1e4, 
                  v0.last = c(49,50))

# view locations in 2D 
plot(x$states[,1], x$states[,2])
# view progression in 1D, removing v0.last from plot
par(mfrow = c(2,1))
plot(x$times, x$states[-1,1], type = 'l', xlab = 'Time', ylab = 'x')
plot(x$times, x$states[-1,2], type = 'l', xlab = 'Time', ylab = 'y')

# look at density of durations
par(mfrow = c(1,1))
plot(density(x$durations))
curve(dexp(x, rate = exp(beta_loc)), add = TRUE, col = 2)
