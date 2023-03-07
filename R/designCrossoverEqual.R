#' Power Calculations for Two-Group Crossover Designs with Equal Variances
#'
#' Approximates the power of equivalence and noninferiority tests for two-sequence, two-period crossover designs with equal variances. One can either find sample sizes that achieve desired statistical power or estimate the power for given sample sizes. More information about the notation and terminology for `diff` and `sigma` is provided in Chapter 10 of Chow et al. (2008).
#'
#' @param diff the anticipated effect size for average (bio)equivalence. Let \eqn{\mu_{jk}} be the population mean of the normal observations in the \eqn{k}th sequence at the \eqn{j}th period. Then, the effect size `diff` is \eqn{\frac{1}{2}(\mu_{11} - \mu_{21}) - \frac{1}{2}(\mu_{12} - \mu_{22})}.
#' @param sigma the anticipated standard deviation for intra-subject comparison, which is assumed to be the same for both groups. That is, \eqn{\sigma} is the standard deviation of all \eqn{y_{i1k} - y_{i2k}} observations, where \eqn{i} denotes the \eqn{i}th subject in the \eqn{k}th sequence.
#' @param deltaL the lower bound for the interval of equivalence (can be `-Inf` for noninferiority test).
#' @param deltaU the upper bound for the interval of equivalence (can be `Inf` for noninferiority test).
#' @param alpha the significance level \eqn{\alpha \in (0,1)}.
#' @param targetPower the desired statistical power of the equivalence or noninferiority test, must be a single number between 0 and 1 (exclusive). Exactly one of the following pairs must be specified: `targetPower` and `q` or `n1` and `n2`. Specify both `targetPower` and `q` to find sample sizes `n1` and `n2` that yield desired power such that \eqn{n_{2} \approx q n_{1}}.
#' @param q the multiplicative constant for imbalanced two-group sample size determination (\eqn{n_{2} = q n_{1}}).
#' @param n1 the sample size for group 1. This is the number of subjects assigned to the first sequence and must be a single integer such that \eqn{n_{1} \ge 2}. Exactly one of the following pairs must be specified: `targetPower` and `q` or `n1` and `n2`. Specify both `n1` and `n2` to estimate statistical power for these sample sizes.
#' @param n2 the sample size for group 2. This is the number of subjects assigned to the first sequence and must be a single integer such that \eqn{n_{2} \ge 2}.
#' @param seed if provided, a single positive integer is used to ensure reproducibility when randomizing the Sobol' sequence via `sobol()` in the `qrng` package (Hofert and Lemieux, 2020).
#' @param sobol one of the following integers: \eqn{s \in \{0, 1, 2, 3, 4 \}}. When approximating the power curve using `targetPower` and `q`, \eqn{2^{s + 10}} points are generated from the Sobol' sequence. When estimating power for given sample sizes `n1` and `n2`, \eqn{2^{s + 16}} points are generated from the Sobol' sequence. The default setting is \eqn{s = 0}, which ensures that each function call should take less than two seconds. As \eqn{s} increases, the sample size calculation takes longer to complete. However, all function calls should still take less than 30 seconds when \eqn{s = 4}.
#'
#' @examples
#' # specify targetPower and q to obtain sample sizes n1 and n2
#' DesignCrossoverEqual(diff = 0.05, sigma = 0.4, deltaL = -0.223, deltaU = 0.223,
#' targetPower = 0.8, q = 1, alpha = 0.05, seed = 1, sobol = 0)
#'
#' # specify n1 and n2 to estimate power for this design
#' DesignCrossoverEqual(diff = 0.05, sigma = 0.4, deltaL = -0.223, deltaU = 0.223,
#' n1 = 18, n2 = 18, alpha = 0.05, seed = 1, sobol = 0)
#'
#' @references
#' Chow, S. C., J. Shao, and H. Wang. (2008). *Sample size calculations in clinical research*. Chapman and Hall/CRC.
#'
#' Hofert, M. and C. Lemieux. (2020). *qrng: (Randomized) Quasi-Random Number Generators*. R package version 0.0-8.
#'
#' @return The sample sizes or power estimate are returned as a list with supplementary information. If `targetPower` and `q` are specified to find sample sizes `n1` and `n2`, a plot of the approximated power curve will also appear in the plot pane.
#' @export
DesignCrossoverEqual <- function(diff = NULL, sigma = NULL, deltaL = -Inf,
                                  deltaU = Inf, alpha = NULL, targetPower = NULL, q = 1, n1 = NULL,
                                  n2 = NULL, seed = NULL, sobol = 0){

  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  ## error checking
  if(!is.numeric(diff) | length(diff) != 1) {
    return("Error: Please specify a valid number for diff.")}
  if(!is.numeric(sigma) | length(sigma) != 1){
    return("Error: Please specify a valid number for sigma.")}
  else if (sigma <= 0){
    return("Error: Please specify a valid number for sigma.")}
  if(!is.numeric(deltaL) | length(deltaL) != 1){
    return("Error: Please specify a valid number for deltaL.")}
  if(!is.numeric(deltaU) | length(deltaU) != 1){
    return("Error: Please specify a valid number for deltaU.")}
  if(deltaL == -Inf & deltaU == Inf){
    return("Error: Please specify valid interval endpoints deltaL and deltaU.")}
  if (deltaL >= deltaU){
    return("Error: Please specify valid interval endpoints deltaL and deltaU.")}
  if(!is.numeric(alpha) | length(alpha) != 1) {
    return("Error: Please specify a valid number for alpha.")}
  if (is.numeric(alpha)){
    if (alpha <= 0 | alpha >= 1){
      return("Error: Please specify a valid number for alpha.")}
  }
  if(sum(c((length(targetPower) != 1 | length(q) != 1), (length(n1) != 1 | length(n2) != 1))) != 1){
    return("Error: Please specify valid inputs for one of the following pairs: targetPower and q or n1 and n2.")}
  if(length(targetPower) == 1 & length(q) == 1){
    if(!is.numeric(targetPower)) {
      return("Error: Please specify a valid number for targetPower.")}
    # if(is.null(targetPower) | !is.numeric(targetPower)) {
    #   return("Error: Please specify a valid number for targetPower.")}
    if (is.numeric(targetPower)){
      if (targetPower <= 0 | targetPower >= 1){
        return("Error: Please specify a valid number for targetPower.")}
      if(!is.numeric(q)) {
        return("Error: Please specify a valid number for q.")}
      else if (is.numeric(q)){
        if (q <= 0) {
          return("Error: Please specify a valid number for q.")}
      }
      if (diff >= deltaU | diff <= deltaL){
        return("Error: Please ensure diff is between deltaL and deltaU.")
      }
    }
  }
  if(length(n1) == 1 & length(n2) == 1){
    if(!is.numeric(n1)) {
      return("Error: Please specify a valid integer for n1.")}
    else if (n1 < 2){
      return("Error: Please specify a valid integer for n1.")}
    else if (n1%%1 != 0){
      return("Error: Please specify valid a integer for n1.")
    }
    if(!is.numeric(n2)) {
      return("Error: Please specify valid a integer for n2.")}
    else if (n2 < 2){
      return("Error: Please specify valid a integer for n2.")}
    else if (n2%%1 != 0){
      return("Error: Please specify valid a integer for n2.")
    }
  }
  if(!is.null(seed) & (!is.numeric(seed) | length(seed) != 1)) {
    return("Error: Please specify a valid seed for random number generation.")}
  else if (!is.null(seed)){
    if (seed%%1 != 0){
      return("Error: Please specify a valid seed for random number generation.")}
  }
  if(!is.numeric(sobol) | length(sobol) != 1) {
    return("Error: Please specify a valid integer between 0 and 4 to initialize the Sobol' sequence.")}
  else if (sobol < 0){
    return("Error: Please specify a valid integer between 0 and 4 to initialize the Sobol' sequence.")}
  else if (!(sobol %in% c(0,1,2,3,4))){
    return("Error: Please specify a valid integer between 0 and 4 to initialize the Sobol' sequence.")}

  sigma <- sigma/2

  if (length(targetPower) == 1 & length(q) == 1){
    targetPowerfn <- function(u, deltaL, deltaU, diff, sigma, n_val, alpha, q){
      n1 <- n_val[1]; n2 <- max(2, q*n_val)
      x <- stats::qchisq(u[1], n1 + n2 - 2)
      z <- stats::qnorm(u[2], diff, sigma*sqrt(1/n1 + 1/n2))

      sdv <- sqrt(x*sigma^2/(n1 + n2 - 2))*sqrt(1/n1 + 1/n2)

      if (deltaU == Inf){
        thres <- (z - deltaL)/stats::qt(1-alpha, n1 + n2 - 2)
      }
      else if (deltaL == -Inf){
        thres <- (deltaU - z)/stats::qt(1-alpha, n1 + n2 - 2)
      }
      else{
        thresUp <- (deltaU - z)/stats::qt(1-alpha, n1 + n2 - 2)
        thresLow <- (z - deltaL)/stats::qt(1-alpha, n1 + n2 - 2)
        thres <- pmin(thresUp, thresLow)
      }

      return(thres - sdv)
    }

    uu <- function(f, lower, upper, tol = 1e-4, maxiter =1000L, ...) {
      f.lower <- f(lower, ...)
      f.upper <- f(upper, ...)
      val <- .External2(stats:::C_zeroin2, function(arg) f(arg, ...),
                        lower, upper, f.lower, f.upper, tol, as.integer(maxiter))
      return(val[1])
    }

    if (is.null(seed)){
      sob <- qrng::sobol(2^(sobol + 10), d = 2, randomize = "digital.shift")
    }
    else{
      sob <- qrng::sobol(2^(sobol + 10), d = 2, randomize = "digital.shift", seed = seed)
    }

    upper_val <- (stats::qnorm(alpha/2) - stats::qnorm(0.9999))^2*(2*sigma^2)/(min(abs(diff-deltaL), abs(diff-deltaU)))^2
    endpoints0_vec <- NULL
    endpoints1_vec <- NULL
    endpoints2_vec <- NULL
    endpoints3_vec <- NULL
    for (i in 1:nrow(sob)){
      endpoints0_vec[i] <- targetPowerfn(n_val = 2,
                                         deltaL = deltaL, deltaU = deltaU, diff = diff, sigma = sigma,
                                         u = sob[i,], alpha = alpha, q = q)
      endpoints1_vec[i] <- targetPowerfn(n_val = ceiling(upper_val/3),
                                         deltaL = deltaL, deltaU = deltaU, diff = diff, sigma = sigma,
                                         u = sob[i,], alpha = alpha, q = q)
      endpoints2_vec[i] <- targetPowerfn(n_val = ceiling(2*upper_val/3),
                                         deltaL = deltaL, deltaU = deltaU, diff = diff, sigma = sigma,
                                         u = sob[i,], alpha = alpha, q = q)
      endpoints3_vec[i] <- targetPowerfn(n_val = ceiling(upper_val),
                                         deltaL = deltaL, deltaU = deltaU, diff = diff, sigma = sigma,
                                         u = sob[i,], alpha = alpha, q = q)
    }

    endpoints_cat <- ifelse(endpoints0_vec >= 0, 0,
                            ifelse(endpoints1_vec >= 0, 1,
                                   ifelse(endpoints2_vec >= 0, 2,
                                          ifelse(endpoints3_vec >= 0, 3, 4))))

    last_group <- which(endpoints_cat == 4)
    if (length(last_group) == 0){
      upper_c <- 2
    }
    else{
      upper_c <- 1
      while(length(last_group) > 0){
        upper_c <- 2*upper_c
        endpoints4_vec <- NULL
        for (i in 1:length(last_group)){
          endpoints4_vec[i] <- targetPowerfn(n_val = ceiling(upper_c*upper_val),
                                             deltaL = deltaL, deltaU = deltaU, diff = diff, sigma = sigma,
                                             u = sob[last_group[i],], alpha = alpha, q = q)
        }
        keep_vec <- ifelse(endpoints4_vec >= 0, FALSE, TRUE)
        last_group <- last_group[keep_vec]
      }
    }

    samps <- NULL
    for (i in 1:nrow(sob)){
      if (endpoints_cat[i] == 0){
        samps[i] <- 2
      }
      else if (endpoints_cat[i] == 1){
        samps[i] <- uu(targetPowerfn, lower =2,
                       upper = ceiling(upper_val/3),
                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma = sigma,
                       u = sob[i,], alpha = alpha, q = q)
      }
      else if (endpoints_cat[i] == 2){
        samps[i] <- uu(targetPowerfn, lower = ceiling(upper_val/3),
                       upper = ceiling(2*upper_val/3),
                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma = sigma,
                       u = sob[i,], alpha = alpha, q = q)
      }
      else if (endpoints_cat[i] == 3){
        samps[i] <- uu(targetPowerfn, lower = ceiling(2*upper_val/3),
                       upper = ceiling(upper_val),
                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma = sigma,
                       u = sob[i,], alpha = alpha, q = q)
      }
      else{
        samps[i] <- uu(targetPowerfn, lower = ceiling(upper_val),
                       upper = ceiling(upper_c*upper_val),
                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma = sigma,
                       u = sob[i,], alpha = alpha, q = q)
      }
    }
    plot(stats::ecdf(samps), main = "Estimated Power Curve", col = "blue")
    funecdf <- stats::ecdf(samps)

    ecdf_root <- function(quant, pwr){return(funecdf(quant) - pwr)}
    n_rough <- uu(ecdf_root, lower = stats::quantile(samps, targetPower*0.5), upper = stats::quantile(samps, targetPower + 0.5*(1 - targetPower)),
                  pwr = targetPower)

    pwrs <- NULL
    n1 <- c(ceiling(n_rough)-1, ceiling(n_rough))
    n2 <- pmax(2, round(q*n1))

    for (j in 1:length(n1)){
      n1_temp <- n1[j]; n2_temp <- n2[j]

      x <- stats::qchisq(sob[,1], n1_temp + n2_temp - 2)
      z <- stats::qnorm(sob[,2], diff, sigma*sqrt(1/n1_temp + 1/n2_temp))

      sdv <- sqrt(x*sigma^2/(n1_temp + n2_temp - 2))*sqrt(1/n1_temp + 1/n2_temp)

      if (deltaU == Inf){
        thresUp <- rep(Inf, length(z))
      }
      else{
        thresUp <- (deltaU - z)/stats::qt(1-alpha, n1_temp + n2_temp - 2)
      }

      if (deltaL == -Inf){
        thresLow <- rep(Inf, length(z))
      }
      else{
        thresLow <- (z - deltaL)/stats::qt(1-alpha, n1_temp + n2_temp - 2)
      }

      thres <- pmin(thresUp, thresLow)

      pwrs <- c(pwrs, mean(ifelse(sdv <= thres,1,0)))
    }

    df_samps <- data.frame(n_plot = samps)

    n_plot <- NULL
    plot_pwr <- ggplot2::ggplot(df_samps, ggplot2::aes(x = n_plot)) +
      ggplot2::stat_ecdf(geom = "step", pad = FALSE, colour = cbbPalette[6], size = 2) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 13)) +
      ggplot2::theme(axis.text.x =  ggplot2::element_text(size = 13)) +
      ggplot2::labs(title = "Approximated Power Curve") +
      ggplot2::labs(y = "Power", x=bquote(italic(n)[1]*"  ("*italic(n)[2]*" = "*.(round(q,3))*italic(n)[1]*")")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,size=20,face="bold",
                                                          margin= ggplot2::margin(0,0,10,0))) +
      ggplot2::theme(axis.title.x = ggplot2::element_text(size = 16, margin= ggplot2::margin(10,0,0,0))) +
      ggplot2::theme(axis.title.y = ggplot2::element_text(size = 16, margin= ggplot2::margin(0,10,0,0)))

    print(plot_pwr)

    type <- ifelse(is.finite(deltaL) & is.finite(deltaU), "a",
                   ifelse(!is.finite(deltaL), "b", "c"))

    METHOD <- switch(type, a = "Equivalence test power calculation (TOST)",
                     b = "Noninferiority test power calculation (for Group 2)",
                     c = "Noninferiority test power calculation (for Group 1)")

    NOTE <- paste0("As a check, power was estimated to be ", pwrs[1], " for n1 = ", n1[1], " and n2 = ", n2[1],".",
                   "\n", "Particularly if q is not 1, you may need to explore sample sizes near n1 and n2.",
                   "\n", "n1 and n2 are the number of subjects assigned to sequences 1 and 2.")

    results <- structure(list(n1 = n1[2], n2 = n2[2], q = q, diff = diff, sigma = 2*sigma,
                              sig.level = alpha, power = pwrs[2], bounds = c(deltaL, deltaU),
                              note = NOTE,
                              method = METHOD), class = "power.htest")
    return(results)
  }

  if (length(n1) == 1 & length(n2) == 1){
    if (is.null(seed)){
      sob <- qrng::sobol(2^(sobol+16), d = 2, randomize = "digital.shift")
    }
    else{
      sob <- qrng::sobol(2^(sobol+16), d = 2, randomize = "digital.shift", seed = seed)
    }

    ## generate two sds and one mean
    x <- stats::qchisq(sob[,1], n1 + n2 - 2)
    z <- stats::qnorm(sob[,2], diff, sigma*sqrt(1/n1 + 1/n2))

    sdv <- sqrt(x*sigma^2/(n1 + n2 - 2))*sqrt(1/n1 + 1/n2)

    if (deltaU == Inf){
      thresUp <- rep(Inf, length(z))
    }
    else{
      thresUp <- (deltaU - z)/stats::qt(1-alpha, n1 + n2 - 2)
    }

    if (deltaL == -Inf){
      thresLow <- rep(Inf, length(z))
    }
    else{
      thresLow <- (z - deltaL)/stats::qt(1-alpha, n1 + n2 - 2)
    }

    thres <- pmin(thresUp, thresLow)

    pwrs <- mean(ifelse(sdv <= thres,1,0))

    type <- ifelse(is.finite(deltaL) & is.finite(deltaU), "a",
                   ifelse(!is.finite(deltaL), "b", "c"))

    METHOD <- switch(type, a = "Equivalence test power calculation (TOST)",
                     b = "Noninferiority test power calculation (for Group 2)",
                     c = "Noninferiority test power calculation (for Group 1)")

    NOTE <- paste0("diff is not in (deltaL, deltaU), so we compute type I error rate.",
                   "\n", "n1 and n2 are the number of subjects assigned to sequences 1 and 2.")


    if(diff >= deltaU | diff <= deltaL){

      NOTE <- paste0("diff is not in (deltaL, deltaU), so we compute type I error rate.",
                     "\n", "n1 and n2 are the number of subjects assigned to sequences 1 and 2.")

      results <- structure(list(n1 = n1, n2 = n2, diff = diff, sigma = 2*sigma,
                                sig.level = alpha, type.I.error = pwrs, bounds = c(deltaL, deltaU),
                                note = NOTE,
                                method = METHOD), class = "power.htest")
    }
    else{
      NOTE <- "n1 and n2 are the number of subjects assigned to sequences 1 and 2."
      results <- structure(list(n1 = n1, n2 = n2, diff = diff, sigma = 2*sigma,
                                sig.level = alpha, power = pwrs, bounds = c(deltaL, deltaU),
                                note = NOTE,
                                method = METHOD), class = "power.htest")
    }

    return(results)
  }
}
