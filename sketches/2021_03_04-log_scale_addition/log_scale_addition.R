logAdd = function(log_a, log_b) {
  # Evaluates log(a + b) when given log(a) and log(b)
  # 
  # Parameters: 
  #  log_a - (vector of) constants to add
  #  log_b - (vector of) constants to add
  
  mapply(function(log_a, log_b) {
    
    x = log_a - log_b
    exp_x = exp(x)
    
    if(exp_x == Inf) {
      log_a
    } else if (exp_x == 0) {
      log_b
    } else {
      log_b + log(1 + exp_x)
    }
    
  }, log_a, log_b)
}

#
# standard test
#

# test values - standard scale
xseq = seq(from = 1e-100, to = 1, length.out = 100)

# data frame of test value pairs
df = expand.grid(a = xseq, b = xseq)
df$log_a = log(df$a)
df$log_b = log(df$b)

# direct addition
df$c = df$a + df$b
df$log_c = log(df$c)

# indirect addition
df$log_c_alt = logAdd(log_a = df$log_a, log_b = df$log_b)
df$c_alt = exp(df$log_c_alt)

# accuracy
max(abs(df$log_c - df$log_c_alt))
max(abs(df$c - df$c_alt))

# distribution of numbers
plot(density(xseq))
plot(density(log(xseq)))


#
# log-test
#

# test values - log scale
xseq_log = seq(from = -900, to = 0, length.out = 100)

# data frame of test value pairs
df = expand.grid(log_a = xseq_log, log_b = xseq_log)
df$a = exp(df$log_a)
df$b = exp(df$log_b)

# direct addition
df$c = df$a + df$b
df$log_c = log(df$c)

# indirect addition
df$log_c_alt = logAdd(log_a = df$log_a, log_b = df$log_b)
df$c_alt = exp(df$log_c_alt)

# accuracy (removing -Inf's due to underflow using standard methods)
# we can see that the direct method is losing accuracy as underflow occurs.
# this can be tested by changing  the minimum value for xseq_log
finite_direct = (df$a != 0) & (df$b != 0)
max(abs(df$log_c[finite_direct] - df$log_c_alt[finite_direct]))
max(abs(df$c - df$c_alt))

# distribution of numbers
plot(density(exp(xseq_log)))
plot(density(xseq_log))

