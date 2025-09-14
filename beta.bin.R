################################################################################
# beta.bin: Bayesian Nonparametric Estimation of Survival Probabilities
#
# This function implements a simple, discrete-time approximation for the survival
# function using a Bayesian beta-binomial approach. The method approximates posterior
# draws for a daily hazard and compiles them into an overall survival probability
# across the observed time range.
#
# Methodological Motivation
# -------------------------
# Inspired by the manuscript’s discussion on “Bayesian inferences of the survival
# function and its summary measures,” which highlights a model-free Bayesian approach
# to constructing posterior distributions of discrete-time hazards. By sampling from
# Beta posteriors for each day’s hazard, we obtain posterior samples for the survival
# curve over the study period.
#
# Args:
#   time.vec   : Numeric vector of event times (e.g., days until event or censoring).
#   censor.vec : Numeric vector, same length as time.vec. 
#                censor.vec[i] = 1 if subject i had the event at time.vec[i],
#                and 0 if the subject was censored at time.vec[i].
#   n.mc       : Integer. Number of Monte Carlo samples from each day's posterior Beta 
#                distribution (default 10,000).
#   alpha      : The alpha parameter of the prior Beta(alpha, beta). Defaults to 0.001 
#                for a weakly informative prior.
#   beta       : The beta parameter of the prior Beta(alpha, beta). Defaults to 1.
#
# Returns:
#   A matrix of dimension length(discretized_time) x n.mc, where each column
#   provides posterior samples (over time) of the survival function.
#
# Details:
#   1. We define a sequence of discrete time points from 1 up to the maximum time 
#      observed in time.vec.
#   2. For each day t, we count the number of events at day t, 
#      and the total number of subjects still at risk at day t.
#   3. We update the posterior distribution for that day's hazard via a Beta conjugate
#      update: rbeta(n.mc, count_death + alpha, current_risk_set + beta).
#   4. We accumulate these hazard draws to build a running product that we store 
#      internally as survival_probability, then finally return 1 - survival_probability,
#      which corresponds to the survival function at each day.
#
# Example usage (toy data):
#   set.seed(123)
#   times   <- c(1,2,2,5,5,5,10,10,10,10)
#   censors <- c(1,1,1,1,0,0,1,0,0,0)
#   result_matrix <- beta.bin(time.vec = times, censor.vec = censors,
#                             n.mc = 1000, alpha = 0.001, beta = 0.01)
#   posterior_means <- apply(result_matrix, 1, mean)
#   plot(posterior_means, type="s", ylab="Posterior Mean Survival", xlab="Day")
#
# Author: 
#   (1) Daniel Paydarfar (developer of this code/tutorial)
#   (2) Adapted from: [Manuscript Title & Methodology Section]
################################################################################

beta.bin <- function(time.vec, censor.vec, n.mc = 10000, alpha = 0.001, beta = 1) {
  
  # Number of subjects
  n <- length(time.vec)  # previously nrow(...) but here time.vec is a vector
  
  # 'product' will track the running product of (1 - daily hazard)
  product <- matrix(1, nrow = n.mc, ncol = 1)
  
  # Create a sequence of discrete times from day 1 to the maximum observed day
  discretized_time <- seq(from = 1, to = max(time.vec), by = 1)
  
  # 'survival_probability' will store partial sums of (hazard * product).
  # (Originally 'cumulative_incidence', but renamed for clarity.)
  survival_probability <- matrix(0,
                                 nrow = length(discretized_time),
                                 ncol = n.mc)
  
  # Iterate over each discrete day (except the last index, to avoid out-of-bounds)
  for (t in 1:(length(discretized_time) - 1)) {
    
    # Count how many deaths occurred exactly at this day
    count_death <- sum(censor.vec[time.vec == discretized_time[t]])
    
    # current_risk_set is how many subjects are still at risk at the start of this day
    current_risk_set <- sum(time.vec >= discretized_time[t])
    
    # Draw daily hazard from the Beta posterior:
    # Posterior = Beta(count_death + alpha, (current_risk_set - count_death) + beta).
    # The code lumps them as shape2 = current_risk_set + beta for simplicity.
    posterior_samples <- rbeta(n.mc,
                               shape1 = count_death + alpha,
                               shape2 = current_risk_set + beta)
    
    # Accumulate partial sums into 'survival_probability'
    survival_probability[t+1, ] <- survival_probability[t, ] + posterior_samples * product
    
    # Update the product with (1 - daily hazard)
    product <- product * (1 - posterior_samples)
  }
  
  # Return the posterior draws for the survival function, i.e. 1 - survival_probability
  return(1 - survival_probability)
}
