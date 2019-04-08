# pragmatic
Pragmatic hypothesis calculation

```{r}
# Generate original data
data <- rnorm(100)

# Define the log-likelihood function
log_f = function(x, mu, sigma) {
  output = 0
    for(i in 1:length(x))
      output = output + log((sigma * sqrt(2*pi))^(-1) * exp(-((x[i]-mu)^2/2*sigma^2)))
  return(output)
}

# Estimate the variance
sigma = var(data)

# Define a secondary function that calculates the log-likelihood considering estimated variance
log_f_sigma = function(x, mu) {
  output = log_f(x, mu, sigma)
}

#Define a function that generate samples considering estimated variance
samples = function(B, mu) {
  return(rnorm(B, mu, sigma))
}

# Calculate the pragmatic hypothesis and summarize it
results = pragmatic(null = 0, epsilon = 0.6, log_f = log_f_sigma, generate_samples = samples, B = 10000, par_grid = seq(-1, 1, length = 201))
summary(results)
```

[1] "Pragmatic hypothesis:"
[1] "(-0.02, 0.02)"
