# Bayesian Stratified Survival Analysis Tutorial

This tutorial demonstrates how to conduct a **stratified** Bayesian survival analysis using a discrete-time **beta-binomial** approach, as described in the paper:
> Bayesian and Frequentist Stratified Analysis of Treatment Effects with Survival Data in Comparative Trials.

We use the **beta.bin** function, which implements the basic building block for obtaining posterior samples of the survival function from a single group (or stratum). Once we have posterior draws for each stratum, we can combine them (via stratified weights) to derive an overall survival curve, restricted mean survival time (RMST), or other summary measures, as discussed in the manuscript.

--------------------------------------------------------------------------------
## 1. Overview

### 1.1 What Is beta.bin?

The beta.bin function computes the **posterior draws** for the survival function in a single population (or stratum) by treating each discrete day as a binomial trial and assigning a Beta prior to that day’s hazard. Because the Beta prior is conjugate to the binomial likelihood, the posterior for each day’s hazard is also Beta, making sampling straightforward and computationally efficient.

- Inputs:
    - time.vec: numeric vector of observed times (e.g., days)
    - censor.vec: numeric vector (same length) with 1 = event at time.vec[i], 0 = censored
    - n.mc: number of Monte Carlo draws (e.g., 10,000)
    - alpha, beta: Beta prior hyperparameters (weakly informative defaults: 0.001, 1)
- Output:
    - A matrix whose rows correspond to time points (from day 1 up to the max of time.vec), and columns correspond to Monte Carlo draws of the survival function.

### 1.2 Stratified Analysis Steps

In a **stratified** design, you have multiple strata (groups) defined by combinations of baseline characteristics—like sex, disease severity, or other covariates. Each stratum’s data is analyzed separately, and then combined in a **weighted** fashion.

1. Identify Strata: Partition your data by stratum (e.g., stratum 1 = “male, EF ≤ 20%”, stratum 2 = “male, EF > 20%”, etc.).
2. Apply beta.bin to Each Stratum: Obtain posterior draws for each stratum’s survival curve.
3. Weight and Combine: Multiply each stratum’s posterior survival by that stratum’s proportion of the overall population, then sum (or otherwise combine) across strata to get an overall survival curve.

--------------------------------------------------------------------------------
## 2. Installing / Loading the Function

### 2.1 Get the Code
1. Download or clone this repository.
2. Look for the file beta.bin.R.

### 2.2 Load into R

    # Set your working directory to the repository folder
    setwd("path/to/your/repo")

    # Source the R script
    source("beta.bin.R")

Alternatively, you can load the function directly into R via the following command in R:

```r
# Source directly from GitHub
source("https://raw.githubusercontent.com/danali-bst/bayesian-survival-tutorial/main/beta.bin.R")
```

Now `beta_bin.R` is available in your R session.

Any of these work. That’s it—beta.bin is now available in your R session.

--------------------------------------------------------------------------------
## 3. Single-Stratum Example

Below we illustrate how to use beta.bin on a **single stratum** of fictitious survival data. You can apply the same procedure to each stratum in a real dataset.

    # Toy data: times + censors for 10 subjects in one stratum
    time.vec   <- c(1, 2, 2, 5, 5, 5, 10, 10, 10, 10)
    censor.vec <- c(1, 1, 1, 1, 0, 0, 1, 0, 0, 0)  # 1=event, 0=censored

    # Run the beta-binomial approach, sampling daily hazards
    set.seed(123)
    result_matrix <- beta.bin(
      time.vec   = time.vec,
      censor.vec = censor.vec, 
      n.mc       = 2000,
      alpha      = 0.001, 
      beta       = 1
    )

    dim(result_matrix)
    # e.g. 10 x 2000 if max(time.vec) = 10

    # Posterior mean survival curve
    posterior_means <- apply(result_matrix, 1, mean)

    # Quick plot
    plot(
      1:length(posterior_means), 
      posterior_means, 
      type = "s",
      xlab = "Days",
      ylab = "Posterior Mean Survival Probability",
      main = "Single-Stratum Bayesian Survival (Toy Data)"
    )

Interpretation: Each row in result_matrix is the time index; each column is a posterior draw. Taking row-wise means or quantiles yields estimates of the survival function at each day.

--------------------------------------------------------------------------------
## 4. Stratified Analysis: Multi-Stratum Demonstration

### 4.1 Simulate (or Load) Data

    # Suppose we have two strata: S1 and S2
    # Each stratum has (time, censor) data for some subset of subjects.

    # Stratum 1
    time_s1   <- c(1,2,3,5,5,7,10,15)
    censor_s1 <- c(1,1,1,0,1,1, 0,1)
    n_s1      <- length(time_s1)

    # Stratum 2
    time_s2   <- c(2,2,2,4,6,6,8,12,15)
    censor_s2 <- c(1,0,1,1,1,0,1, 0, 0)
    n_s2      <- length(time_s2)


### 4.2 Apply beta.bin to Each Stratum

    set.seed(123)
    res_s1 <- beta.bin(time.vec = time_s1,
                       censor.vec = censor_s1,
                       n.mc = 2000)
    res_s2 <- beta.bin(time.vec = time_s2,
                       censor.vec = censor_s2,
                       n.mc = 2000)

    dim(res_s1)  # rows = max(time_s1), columns = 2000
    dim(res_s2)  # rows = max(time_s2), columns = 2000

### 4.3 Combine Posterior Estimates Across Strata

To get an overall survival curve for the entire population (both strata):

1. Weight each curve by the stratum’s proportion in the overall sample.
2. Sum those weighted survival curves at each time point.

```r
#Follow-up time used for plotting / RMST (if follow-up times between strata not the same,
#see note in discussion)
max_day = max(time_s1, time_s2)

#Weighted average of the posterior survival
prop_s1 <- n_s1 / (n_s1 + n_s2)
prop_s2 <- n_s2 / (n_s1 + n_s2)

overall_posterior <- prop_s1 * res_s1 + prop_s2 * res_s2
overall_mean <- apply(overall_posterior, 1, mean)
```

    # Quick plot
    plot(
      1:max_day,
      overall_mean,
      type = "s",
      xlab = "Days",
      ylab = "Posterior Mean Survival (Combined Strata)",
      main = "Stratified Bayesian Survival (2-Stratum Example)"
    )

Result: overall_posterior is a matrix of size [max_day x n.mc], combining the draws from both strata, weighted by their relative sizes.

### 4.4 Additional Summaries (e.g., RMST)

    # Example: compute restricted mean survival time (RMST) up to day tau
    tau <- 12
    rmst_draws <- colSums(overall_posterior[1:tau, ])
    rmst_mean  <- mean(rmst_draws)
    rmst_ci    <- quantile(rmst_draws, c(0.025, 0.975))

    cat("Posterior Mean RMST up to day", tau, "=", rmst_mean, "\n")
    cat("95% Credible Interval =", rmst_ci, "\n")

--------------------------------------------------------------------------------
## 5. Interpretation & Discussion

- **Probabilistic Inference**: With Bayesian posterior draws, a 95% credible interval can be read as “there’s a 95% chance the true survival or RMST lies in this range,” in contrast to frequentist confidence intervals.
- **Stratified Weights**: By weighting each stratum’s survival draws by stratum size, we get an overall measure that respects the original stratification.
- **Sensitivity Analysis**: If you suspect different priors, you can change alpha/beta and see how results shift.

--------------------------------------------------------------------------------
## 6. Methodology References

1. Manuscript section “Bayesian inferences of the survival function and its summary measures.”
2. Beta-Binomial Background: Hjort (1990)
3. RMST Discussion: Various references in manuscript about the advantages over hazard ratios.

--------------------------------------------------------------------------------
## 7. Closing Notes

- This method does not rely on proportional hazards or a hazard ratio.
- It is straightforward, requiring only binomial likelihoods and Beta priors.
- For real-world data with more strata, repeat the same pattern: call beta.bin for each stratum, align row lengths, then weight and sum the posterior draws.
- If strata have different follow-up, use the minimum of the maximimum follow-up times across strata as your truncation time for RMST and as your last observation in time.vec for generating survival curve, to ensure there is no extrapolation. 

**We hope this tutorial helps you explore a model-free Bayesian approach to stratified survival analysis!** For questions or contributions, please open an issue or pull request in this repository, or contact Daniel Paydarfar at danielpaydarfar@fas.harvard.edu.
