#' Pragmatic hypothesis calculation
#'
#' @description The pragmatic function calculates a pragmatic hypothesis basic on a simple one.
#'
#' @param null Numeric. The original null hypothesis point.
#' @param epsilon Numeric. The desired dissimilarity threshold.
#' @param log_f Function. A function that computes the log-likelihood of a new experiment. See also 'Examples'.
#' @param generate_samples Function. A funtion that generates new samples from the original distribution. See also 'Examples'.
#' @param B Numeric. Number of Monte Carlo simulations. Higher values might take longer to run but lead to more accurate results. Default set to 1000 simulations.
#' @param connected Boolean. For default, the pragmatic function assumes the real pragmatic is connect (a convex set), if this does not suit a particular example, set this parameter to FALSE.
#' @param par_grid Numeric vector. If connected = FALSE, par_grid have to be a grid of points to test belonginess to the pragmatic. Ignored if connected = TRUE.
#' @param symmetrical Boolean. If one knows in advance that the real pragmatic is symmetrical around the null hypothesis, setting this to TRUE might speed up calculations and lead to more accurate results.
#' @param ... Additional arguments to be passed to lower level functions.
#'
#' @details The pragmatic function performs a grid search for the pragmatic hypothesis related to a simple one. For default, the function assumes thats the real pragmatic is a convex set, and them, performs a conservative binary search to find the set limits. When there is no information weather the pragmatic is a convex set, a user defined grid is needed to perform the search.
#'
#' @return A string determining the pragmatic hypothesis set.
#' @export pragmatic
#'
#' @examples
#'
#' # Generate original data
#' data <- rnorm(100)
#'
#' # Define the log-likelihood function
#' log_f = function(x, mu, sigma) {
#'   sum(dnorm(x, mu, sigma, log = TRUE))
#' }
#'
#' # Estimate the data's standard deviation
#' sigma_data = sd(data)
#'
#' # Define a secondary function that calculates the log-likelihood considering estimated deviation
#' log_f_sigma = function(x, mu) {
#'   log_f(x, mu, sigma_data)
#' }
#'
#' #Define a function that generate samples considering estimated variance
#' samples = function(mu) {
#'   rnorm(10, mu, sigma_data) # 10 is the number of new samples Z
#' }
#'
#' # Calculate the pragmatic hypothesis
#' pragmatic(null = 0, epsilon = 0.6, log_f = log_f_sigma, generate_samples = samples, symmetrical = TRUE)

pragmatic <- function(null, epsilon, log_f, generate_samples, B = 1000, connected = TRUE, par_grid = NULL, symmetrical = FALSE, ...) {

  if(connected) {

    if(!is.null(par_grid))
      warning('connected = TRUE, ignoring par_grid.')

    results = .connected_pragmatic(null, epsilon, log_f, generate_samples, B, symmetrical, ...)
    string = paste(c('(', paste(round(results[1], 4), round(results[2], 4), sep = ', '), ')'), collapse = '')

  } else {

    if(is.null(par_grid))
      stop('when connected = FALSE, par_grid must be given.')

    output = .unconnected_pragmatic(null, epsilon, log_f, generate_samples, B, par_grid)

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

    string = paste(c('(', paste(round(limits[[1]][1], 4), round(limits[[1]][2], 4), sep = ', '), ')'), collapse = '')

    if(length(limits) > 1) {

      for (i in 2:length(limits)) {
        string = paste(string, paste(c('(', paste(round(limits[[i]][1], 4), round(limits[[i]][2], 4), sep = ', '), ')'), collapse = ''), sep = ' U ')
      }
    }
  }

  return(string)
}


.unconnected_pragmatic <- function(null, epsilon, log_f, generate_samples, B, par_grid) {

  browser()

  z0 <- replicate(B, generate_samples(null))

  log_f_z0_null <- apply(z0, 2, log_f, mu = null)

  pragmatic <- logical(length = length(par_grid))

  for (i in 1:length(par_grid)) {

    log_f_z0_par <- apply(z0, 2, log_f, mu = par_grid[i])

    z1 <- replicate(B, generate_samples(par_grid[i]))

    log_f_z1_null <- apply(z1, 2, log_f, mu = null)
    log_f_z1_par <- apply(z1, 2, log_f, mu = par_grid[i])

    theta0 <- mean(log_f_z0_null > log_f_z0_par)
    theta1 <- mean(log_f_z1_par > log_f_z1_null)
    dist <- (theta0 + theta1)/2

    pragmatic[i] <- (dist < epsilon)

  }

  output <- list(pragmatic = pragmatic, par_grid = par_grid, epsilon = epsilon)

  return(output)
}


.connected_pragmatic <- function(null, epsilon, log_f, generate_samples, B, symmetrical = FALSE, step = 0.999, start = 0.2) {

  z0 <- replicate(B, generate_samples(null))
  log_f_z0_null <- apply(z0, 2, log_f, mu = null)

  pragmatic <- c(null, null)

  dot = null + start

  dec = FALSE

  repeat {

    log_f_z0_par <- apply(z0, 2, log_f, mu = dot)

    z1 <- replicate(B, generate_samples(dot))

    log_f_z1_null <- apply(z1, 2, log_f, mu = null)
    log_f_z1_par <- apply(z1, 2, log_f, mu = dot)

    theta0 <- mean(log_f_z0_null > log_f_z0_par)
    theta1 <- mean(log_f_z1_par > log_f_z1_null)
    dist <- (theta0 + theta1)/2

    if(dist < epsilon) {

      pragmatic[2] = dot

      if(dec) {

        break

      } else {

        dot = dot + 1

      }

    } else {

      if(dot > 0) {

        dot = step*dot

      } else {

        dot = dot*(2-step)

      }

      dec = TRUE

    }

  }

    if(symmetrical) {

      pragmatic[1] <- null - (dot - null)

    } else {

      dot = null - start

      inc = FALSE

      repeat {

        log_f_z0_par <- apply(z0, 2, log_f, mu = dot)

        z1 <- replicate(B, generate_samples(dot))

        log_f_z1_null <- apply(z1, 2, log_f, mu = null)
        log_f_z1_par <- apply(z1, 2, log_f, mu = dot)

        theta0 <- mean(log_f_z0_null > log_f_z0_par)
        theta1 <- mean(log_f_z1_par > log_f_z1_null)
        dist <- (theta0 + theta1)/2

        if(dist < epsilon) {

          pragmatic[1] = dot

          if(inc) {

            break

          } else {

            dot = dot - 1

          }

        } else {

          if(dot < 0) {

            dot = step*dot

          } else {

            dot = dot*(2-step)

          }

          inc = TRUE

        }

    }

  }

  return(pragmatic)
}
