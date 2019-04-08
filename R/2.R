#' Pragmatic hypothesis calculation
#'
#' @description The pragmatic function calculates a pragmatic hypothesis basic on a simple one.
#'
#' @param null The original null hypothesis point.
#' @param epsilon Dissimilarity threshold.
#' @param log_f A function that computes the log-likelihood of the new experiment.
#' @param generate_samples A funtion that generates new samples from the original distribution.
#' @param B Number of Monte Carlo simulations.
#' @param par_grid Grid of parameters to test.
#' @param grid_lim If par_grid is not defined, this controls the limits of an internal constructed grid.
#' @param grid_dots If par_grid is not defined, this controls the length of an internal constructed grid.
#'
#' @return A list containing: The used parameter grid, the used threshold, and a logical vector indicating which of the parameters from the grid belongs to the pragmatic hypothesis.
#' @export pragmatic
#' @S3method summary,pragmatic
#'
#' @examples
#'
#' # Generate original data
#' data <- rnorm(100)
#'
#' # Define the log-likelihood function
#' log_f = function(x, mu, sigma) {
#'   output = 0
#'     for(i in 1:length(x))
#'       output = output + log((sigma * sqrt(2*pi))^(-1) * exp(-((x[i]-mu)^2/2*sigma^2)))
#'   return(output)
#' }
#'
#' # Estimate the variance
#' sigma = var(data)
#'
#' # Define a secondary function that calculates the log-likelihood considering estimated variance
#' log_f_sigma = function(x, mu) {
#'   output = log_f(x, mu, sigma)
#' }
#'
#' #Define a function that generate samples considering estimated variance
#' samples = function(B, mu) {
#'   return(rnorm(B, mu, sigma))
#' }
#'
#' # Calculate the pragmatic hypothesis and summarize it
#' results = pragmatic(null = 0, epsilon = 0.6, log_f = log_f_sigma, generate_samples = samples, B = 10000, par_grid = seq(-1, 1, length = 201))
#' summary(results)
#'

pragmatic <- function(null, epsilon, log_f, generate_samples, B, par_grid = NULL, grid_lim = NULL, grid_dots = NULL) {

  if(is.null(par_grid) & (is.null(grid_lim) | is.null(grid_dots))) {

    stop('Either par_grid or both grid_lim and grid_dots must be specified')

  } else if(is.null(par_grid)) {

    par_grid <- seq(grid_lim[1], grid_lim[2], length = par_dots)

  } else if(!is.null(grid_lim) | !is.null(grid_dots)) {

    warning('par_grid already defined, ignoring grid_lim and grid_dots')

  }

  z0 <- generate_samples(B, null)

  pragmatic <- logical(length = length(par_grid))

  for (i in 1:length(par_grid)) {

    z1 <- generate_samples(B, par_grid[i])

    theta0 <- mean(log_f(z0, null) > log_f(z0, par_grid[i]))
    theta1 <- mean(log_f(z1, par_grid[i]) > log_f(z1, null))
    dist <- (theta0 + theta1)/2

    pragmatic[i] <- (dist < epsilon)

  }

  output <- list(pragmatic = pragmatic, par_grid = par_grid, epsilon = epsilon)
  class(output) <- 'pragmatic'

  return(output)
}



summary.pragmatic <- function(output) {

  if(sum(output$pragmatic) == 0) {
    stop('Pragmatic is empty')
  }

  first <- min(which(output$pragmatic == TRUE))
  last <- min(which(output$pragmatic[-(1:first)] == FALSE)) + first - 1

  limits <- list(c(output$par_grid[first], output$par_grid[last]))

  repeat {

    first <- suppressWarnings(min(which(output$pragmatic[-(1:(last))] == TRUE)) + last)

    if(first == Inf)
      break

    last <- min(which(output$pragmatic[-(1:first)] == FALSE)) + first - 1
    limits[[length(limits)+1]] <- c(output$par_grid[first], output$par_grid[last])

  }

  print('Pragmatic hypothesis:')
  if(length(limits) > 1) {
    for (i in 1:(length(limits)-1)) {
      print(paste0('(', round(limits[[i]][1], 4), ', ', round(limits[[i]][2], 4), ') ', 'U'))
    }
    print(paste0('(', round(limits[[i+1]][1], 4), ', ', round(limits[[i+1]][2], 4), ')'))
  } else {
    print(paste0('(', round(limits[[1]][1], 4), ', ', round(limits[[1]][2], 4), ')'))
  }
}


