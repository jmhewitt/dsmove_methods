# simulation transition rate parameter
beta_loc = 1

# simulate trajectory
x = ctds.quicksim(dims = c(100, 100), beta_loc = beta_loc, beta_ar = 1, 
                  v0 = c(50,50), t0 = 0, tf = 1e2, max.steps = 1e4, 
                  v0.last = c(49,50))

# observe trajectory, making sure to remove v0.last from ctds.observe input
tobs = seq(from = x$times[1], to = tail(x$times, 1), length.out = 100)
x.obs = ctds.observe(states = x$states[-1,], times = x$times, t.obs = tobs)

# view simulated locations in 2D 
plot(x$states[,1], x$states[,2])
# highlight observations
points(x.obs$states[,1], x.obs$states[,2], col = 2, pch = 4)

# view progression in 1D
par(mfrow = c(2,1))
plot(stepfun(x = x$times, y = x$states[,1]), pch = '.',
     xlab = 'Time', ylab = 'x', main = '')
points(x.obs$times, x.obs$states[,1], col = 2)
plot(stepfun(x = x$times, y = x$states[,2]), pch = '.',
     xlab = 'Time', ylab = 'y', main = '')
points(x.obs$times, x.obs$states[,2], col = 2)

# look at density of durations
par(mfrow = c(1,1))
plot(density(x$durations))
curve(dexp(x, rate = exp(beta_loc)), add = TRUE, col = 2)
