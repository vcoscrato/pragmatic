# pragmatic
Pragmatic hypothesis calculation

## Installation

You can install pragmatic directly into R from github with devtools:

``` r
if(!("devtools" %in% rownames(installed.packages())))
  install.packages('devtools')
devtools::install_github("vcoscrato/pragmatic")
```

## Example
This example performs pragmatic hypothesis calculation on for a t-test.

```{r}
# Generate original data
data <- rnorm(100)

# Define the log-likelihood function
log_f = function(x, mu, sigma) {
  sum(dnorm(x, mu, sigma, log = TRUE))
}

# Estimate the data's standard deviation
sigma_data = sd(data)

# Define a secondary function that calculates the log-likelihood considering estimated deviation
log_f_sigma = function(x, mu) {
  log_f(x, mu, sigma_data)
}

#Define a function that generate samples considering estimated variance
samples = function(B, mu) {
  rnorm(B, mu, sigma_data)
}

# Calculate the pragmatic hypothesis
pragmatic(null = 0, epsilon = 0.6, log_f = log_f_sigma, generate_samples = samples, symmetrical = TRUE)
```

[1] "(-0.1844, 0.1844)"
