#' beta.bin: Bayesian Nonparametric Estimation of Survival Probabilities
#'
#' This function implements a simple, discrete-time approximation for the survival
#' function using a Bayesian beta-binomial approach. The method
#' approximates posterior draws for a daily hazard and compiles them into an overall
#' survival probability across the observed time range.
#'
#' Methodological Motivation:
#' --------------------------------
#' This code is inspired by the methodology discussed in the manuscript’s section
#' “Bayesian inferences of the survival function and its summary measures,” which
#' highlights a model-free Bayesian approach to constructing posterior distributions
#' of discrete-time hazards. By sampling from Beta posteriors for each day’s hazard,
#' we obtain posterior samples for the survival curve over the study period.
#'
#' @param time.vec A numeric vector of event times (e.g., days until event or censoring)
#'   for each subject.
#' @param censor.vec A numeric vector of the same length as time.vec, where
#'   \code{censor.vec[i] = 1} if subject \code{i} had the event observed at time.vec[i],
#'   and \code{censor.vec[i] = 0} if the subject was censored at time.vec[i].
#' @param n.mc The number of Monte Carlo samples to draw from each daily posterior
#'   Beta distribution. Defaults to 10,000.
#' @param alpha The \eqn{\alpha} parameter of the prior Beta(\eqn{\alpha},\eqn{\beta})
#'   for each day’s hazard. Defaults to 0.001, reflecting a weakly informative prior.
#' @param beta The \eqn{\beta} parameter of the prior Beta(\eqn{\alpha},\eqn{\beta})
#'   for each day’s hazard. Defaults to 0.01.
#'
#' @return A matrix of dimension \code{length(discretized_time)} by \code{n.mc}, where
#'   each column provides the posterior samples (over time) of the survival function.
#'
#' @details
#' 1. We define a sequence of discrete time points from 1 up to the maximum time
#'    observed in \code{time.vec}.
#' 2. For each time \eqn{t}, we count the number of events that occur exactly at day \eqn{t},
#'    and the total number of subjects still at risk at the start of day \eqn{t}.
#' 3. We update the posterior distribution for that day's hazard via a Beta conjugate
#'    update: \eqn{rbeta(n.mc, count\_death + alpha, current\_risk\_set + beta)}.
#' 4. We accumulate these hazard draws to obtain a running product that feeds into
#'    what we call \code{survival_probability} (the internal array in the code), and
#'    finally return \code{1 - survival_probability}, which aligns with the survival
#'    function at each day.
#'
#' @examples
#' # Example usage (toy data):
#' set.seed(123)
#' # Suppose we have 10 subjects with small integer event/censor times
#' times  <- c(1, 2, 2, 5, 5, 5, 10, 10, 10, 10)
#' censors <- c(1, 1, 1, 1, 0, 0, 1, 0, 0, 0)  # 1=event, 0=censored
#' result_matrix <- beta.bin(time.vec = times, censor.vec = censors,
#'                           n.mc = 1000, alpha = 0.001, beta = 0.01)
#' # The rows of result_matrix are time points; the columns are posterior survival draws
#' # You can summarize them, e.g.:
#' posterior_means <- apply(result_matrix, 1, mean)
#'
#' @author
#' (1) Daniel Paydarfar, the developer of this code/tutorial
#' (2) Adapted from the approach described in: [Manuscript Title & Methodology Section]
#'
beta.bin <- function(time.vec, censor.vec, n.mc = 10000, alpha = 0.001, beta = 0.01) {
  
  # Number of subjects
  n <- length(time.vec)  # changed nrow(...) -> length(...) if time.vec is a 1D vector
  
  # 'product' will track the running product of (1 - daily hazard)
  product <- matrix(1, nrow = n.mc, ncol = 1)
  
  # Create a sequence of discrete times from day 1 to the maximum observed day
  discretized_time <- seq(from = 1, to = max(time.vec), by = 1)
  
  # survival_probability will internally store partial sums of hazard * product
  # (though conceptually it's representing 'cumulative incidence' in the old code).
  survival_probability <- matrix(0,
                                 nrow = length(discretized_time),
                                 ncol = n.mc)
  
  # Track the posterior draws across discrete times
  for (t in 1:(length(discretized_time) - 1)) {
    
    # Count how many deaths occurred exactly at day t
    count_death <- sum(censor.vec[time.vec == discretized_time[t]])
    
    # The risk set is the number of subjects who have not had the event by day t
    current_risk_set <- sum(time.vec >= discretized_time[t])
    
    # Draw from Beta posterior for the daily hazard
    # Posterior = Beta(count_death + alpha, current_risk_set - count_death + beta)
    # but the code effectively lumps them as shape2 = current_risk_set + beta
    # because we treat the binomial denominator as 'current_risk_set'
    posterior_samples <- rbeta(n.mc,
                               shape1 = count_death + alpha,
                               shape2 = current_risk_set + beta)
    
    # Accumulate partial sums into 'survival_probability'
    survival_probability[t+1, ] <- survival_probability[t, ] + posterior_samples * product
    
    # Update the product: multiply by (1 - hazard)
    product <- product * (1 - posterior_samples)
  }
  
  # The final returned value is the posterior for the survival function:
  return(1 - survival_probability)
}
