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
# Generate data
data <- rnorm(100)

# Define the log-likelihood function
log_f = function(x, mu, sigma) {
  sum(dnorm(x, mu, sigma, log = TRUE))
}

# Estimate the data's standard deviation
sigma_data = sd(data)

# Define a secondary function that calculates the log-likelihood considering estimated variance
log_f_sigma = function(x, mu) {
  output = log_f(x, mu, sigma_data)
}

#Define a function that generates samples considering estimated variance
samples = function(B, mu) {
  return(rnorm(B, mu, sigma_data))
}

# Calculate the pragmatic hypothesis and summarize it
results = pragmatic(null = 0, epsilon = 0.6, log_f = log_f_sigma, 
                    generate_samples = samples, B = 10000, 
                    par_grid = seq(-1, 1, length = 201))
summary(results)
```

[1] "Pragmatic hypothesis:"

[1] "(-0.01, 0)"
