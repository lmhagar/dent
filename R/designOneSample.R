#' Power Calculations for One-Group Designs
#'
#' Approximates the power of equivalence and noninferiority tests for one sample. One can either find a sample size that achieves desired statistical power or estimate the power for a given sample size.
#'
#' @param mu the anticipated mean of the normal population (\eqn{\mu}).
#' @param sigma the anticipated population standard deviation (\eqn{\sigma}).
#' @param deltaL the lower bound for the interval of equivalence (can be `-Inf` for nonsuperiority test).
#' @param deltaU the upper bound for the interval of equivalence (can be `Inf` for noninferiority test).
#' @param alpha the significance level \eqn{\alpha \in (0,1)}.
#' @param targetPower the desired statistical power of the equivalence or noninferiority test, must be a single number between 0 and 1 (exclusive). Exactly one of the following inputs must be specified: `targetPower` or `n`. Specify `targetPower` to find a sample size `n` that yields desired power.
#' @param n the sample size, must be a single integer such that \eqn{n \ge 2}. Exactly one of the following inputs must be specified: `targetPower` or `n`. Specify `n` to estimate statistical power for this sample size.
#' @param plot a logical variable indicating whether to return a plot of the power curve. If `n` is specified instead of `targetPower`, this variable is automatically set to `FALSE`. If you wish to approximate many power curves, suppressing the plots will expedite this process.
#' @param seed if provided, a single positive integer is used to ensure reproducibility when randomizing the Sobol' sequence via `sobol()` in the `qrng` package.
#' @param sobol one of the following integers: \eqn{s \in \{0, 1, 2, 3, 4 \}}. When approximating the power curve using `targetPower`, \eqn{2^{s + 10}} points are generated from the Sobol' sequence. When estimating power for a given sample size `n`, \eqn{2^{s + 16}} points are generated from the Sobol' sequence. The default setting is \eqn{s = 0}, which ensures that each function call should take less than two seconds.  As \eqn{s} increases, the sample size calculation is less sensitive to simulation variability but takes longer to complete. However, all function calls should still take less than 30 seconds when \eqn{s = 4}.
#'
#' @examples
#' # specify targetPower to obtain sample size n
#' DesignOneSample(mu = -4, sigma = 15, deltaL = -19.2, deltaU = 19.2,
#' targetPower = 0.8, alpha = 0.05, plot = TRUE, seed = 1, sobol = 0)
#'
#' # specify n to estimate power for this design
#' DesignOneSample(mu = -4, sigma = 15, deltaL = -19.2, deltaU = 19.2,
#' n = 17, alpha = 0.05, plot = FALSE, seed = 1, sobol = 0)
#'
#' @return The sample size or power estimate are returned as a list with supplementary information. If `targetPower` is specified to find sample size `n`, a plot of the approximated power curve will also appear in the plot pane if `plot = TRUE`. To confirm the sample size recommendation, power will be approximated at sample size `n - 1`. This power estimate should be *less* than `targetPower`. To find a sample size that corresponds to a different `targetPower`, save this function's output to an object and use the `UpdateTargetPower()` function.
#' @export
DesignOneSample <- function(mu = NULL, sigma = NULL, deltaL = -Inf,
                                  deltaU = Inf, alpha = NULL, targetPower = NULL,
                                  n = NULL, plot = TRUE, seed = NULL, sobol = 0){

  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  ## error checking
  if(!is.numeric(mu) | length(mu) != 1) {
    stop("Please specify a valid number for mu.")}
  if(!is.numeric(sigma) | length(sigma) != 1){
    stop("Please specify a valid number for sigma.")}
  else if (sigma <= 0 | !is.finite(sigma)){
    stop("Please specify a valid number for sigma.")}
  if(!is.numeric(deltaL) | length(deltaL) != 1){
    stop("Please specify a valid number for deltaL.")}
  if(!is.numeric(deltaU) | length(deltaU) != 1){
    stop("Please specify a valid number for deltaU.")}
  if(deltaL == -Inf & deltaU == Inf){
    stop("Please specify valid interval endpoints deltaL and deltaU.")}
  if (deltaL >= deltaU){
    stop("Please specify valid interval endpoints deltaL and deltaU.")}
  if(!is.numeric(alpha) | length(alpha) != 1) {
    stop("Please specify a valid number for alpha.")}
  if (is.numeric(alpha)){
    if (alpha <= 0 | alpha >= 1){
      stop("Please specify a valid number for alpha.")}
  }
  if(sum(c(length(targetPower) != 1, length(n) != 1)) != 1){
    stop("Please specify valid inputs for one of the following: targetPower or n.")}
  if(length(targetPower) == 1){
    if(!is.numeric(targetPower)) {
      stop("Please specify a valid number for targetPower.")}
    # if(is.null(targetPower) | !is.numeric(targetPower)) {
    #   stop("Please specify a valid number for targetPower.")}
    if (is.numeric(targetPower)){
      if (targetPower <= 0 | targetPower >= 1){
        stop("Please specify a valid number for targetPower.")}
      if (mu >= deltaU | mu <= deltaL){
        stop("Please ensure mu is between deltaL and deltaU.")
      }
    }
  }
  if(length(n) == 1){
    if(!is.numeric(n)) {
      stop("Please specify a valid integer for n.")}
    else if (n < 2){
      stop("Please specify a valid integer for n.")}
    else if (n%%1 != 0){
      stop("Please specify a valid integer for n.")
    }
  }
  if(!is.null(seed) & (!is.numeric(seed) | length(seed) != 1)) {
    stop("Please specify a valid seed for random number generation.")}
  else if (!is.null(seed)){
    if (seed%%1 != 0){
      stop("Please specify a valid seed for random number generation.")}
  }
  if(!is.numeric(sobol) | length(sobol) != 1) {
    stop("Please specify a valid integer between 0 and 4 to initialize the Sobol' sequence.")}
  else if (sobol < 0){
    stop("Please specify a valid integer between 0 and 4 to initialize the Sobol' sequence.")}
  else if (!(sobol %in% c(0,1,2,3,4))){
    stop("Please specify a valid integer between 0 and 4 to initialize the Sobol' sequence.")}

  if (length(plot) != 1 | !(plot %in% c(TRUE, FALSE))){
    stop("Please provide valid logical input for plot.")
  }
  if (length(n) == 1){
    plot <- FALSE
  }

  if (length(targetPower) == 1){
    targetPowerfn <- function(u, deltaL, deltaU, mu, sigma, n_val, alpha){
      n <- n_val
      x <- stats::qchisq(u[1], n - 1)
      z <- stats::qnorm(u[2], mu, sigma*sqrt(1/n))

      sdv <- sqrt(x*sigma^2/(n - 1))*sqrt(1/n)

      if (deltaU == Inf){
        thres <- (z - deltaL)/stats::qt(1-alpha, n - 1)
      }
      else if (deltaL == -Inf){
        thres <- (deltaU - z)/stats::qt(1-alpha, n - 1)
      }
      else{
        thresUp <- (deltaU - z)/stats::qt(1-alpha, n - 1)
        thresLow <- (z - deltaL)/stats::qt(1-alpha, n - 1)
        thres <- pmin(thresUp, thresLow)
      }

      return(thres - sdv)
    }

    uu <- function (fun, lower, upper, maxiter = 1000, tol = 1e-4, ...)
    {
      f <- function(x) fun(x, ...)
      x1 <- lower
      f1 <- f(x1)
      x2 <- upper
      f2 <- f(x2)
      x3 <- 0.5 * (lower + upper)
      niter <- 1
      while (niter <= maxiter) {
        f3 <- f(x3)
        if (abs(f3) < tol) {
          x0 <- x3
          return(x0)
        }
        if (f1 * f3 < 0) {
          upper <- x3}
        else {lower <- x3}
        if ((upper - lower) < tol * max(abs(upper), 1)) {
          x0 <- 0.5 * (lower + upper)
          return(x0)
        }
        denom <- (f2 - f1) * (f3 - f1) * (f2 - f3)
        numer <- x3 * (f1 - f2) * (f2 - f3 + f1) + f2 * x1 *
          (f2 - f3) + f1 * x2 * (f3 - f1)
        if (denom == 0) {
          dx <- upper - lower
        }
        else {
          dx <- f3 * numer/denom
        }
        x <- x3 + dx
        if ((upper - x) * (x - lower) < 0) {
          dx <- 0.5 * (upper - lower)
          x <- lower + dx
        }
        if (x1 < x3) {
          x2 <- x3
          f2 <- f3
        }
        else {
          x1 <- x3
          f1 <- f3
        }
        niter <- niter + 1
        if (abs(x - x3) < tol) {
          x0 <- x
          return(x0)
        }
        x3 <- x
      }
      return(x0)
    }

    if (is.null(seed)){
      seed <- ceiling(1000*stats::runif(1))
    }
    sob <- qrng::sobol(2^(sobol + 10), d = 2, randomize = "digital.shift", seed = seed)

    upper_val <- (stats::qnorm(alpha/2) - stats::qnorm(0.9999))^2*(2*sigma^2)/(min(abs(mu-deltaL), abs(mu-deltaU)))^2
    endpoints0_vec <- NULL
    endpoints1_vec <- NULL
    endpoints2_vec <- NULL
    endpoints3_vec <- NULL
    for (i in 1:nrow(sob)){
      endpoints0_vec[i] <- targetPowerfn(n_val = 2,
                                         deltaL = deltaL, deltaU = deltaU, mu = mu, sigma = sigma,
                                         u = sob[i,], alpha = alpha)
      endpoints1_vec[i] <- targetPowerfn(n_val = ceiling(upper_val/3),
                                         deltaL = deltaL, deltaU = deltaU, mu = mu, sigma = sigma,
                                         u = sob[i,], alpha = alpha)
      endpoints2_vec[i] <- targetPowerfn(n_val = ceiling(2*upper_val/3),
                                         deltaL = deltaL, deltaU = deltaU, mu = mu, sigma = sigma,
                                         u = sob[i,], alpha = alpha)
      endpoints3_vec[i] <- targetPowerfn(n_val = ceiling(upper_val),
                                         deltaL = deltaL, deltaU = deltaU, mu = mu, sigma = sigma,
                                         u = sob[i,], alpha = alpha)
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
                                             deltaL = deltaL, deltaU = deltaU, mu = mu, sigma = sigma,
                                             u = sob[last_group[i],], alpha = alpha)
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
                       deltaL = deltaL, deltaU = deltaU, mu = mu, sigma = sigma,
                       u = sob[i,], alpha = alpha)
      }
      else if (endpoints_cat[i] == 2){
        samps[i] <- uu(targetPowerfn, lower = ceiling(upper_val/3),
                       upper = ceiling(2*upper_val/3),
                       deltaL = deltaL, deltaU = deltaU, mu = mu, sigma = sigma,
                       u = sob[i,], alpha = alpha)
      }
      else if (endpoints_cat[i] == 3){
        samps[i] <- uu(targetPowerfn, lower = ceiling(2*upper_val/3),
                       upper = ceiling(upper_val),
                       deltaL = deltaL, deltaU = deltaU, mu = mu, sigma = sigma,
                       u = sob[i,], alpha = alpha)
      }
      else{
        samps[i] <- uu(targetPowerfn, lower = ceiling(upper_val),
                       upper = ceiling(upper_c*upper_val),
                       deltaL = deltaL, deltaU = deltaU, mu = mu, sigma = sigma,
                       u = sob[i,], alpha = alpha)
      }
    }
    # plot(stats::ecdf(samps), main = "Estimated Power Curve", col = "blue")
    funecdf <- stats::ecdf(samps)

    ecdf_root <- function(quant, pwr){return(funecdf(quant) - pwr)}
    n_rough <- uu(ecdf_root, lower = stats::quantile(samps, targetPower*0.5), upper = stats::quantile(samps, targetPower + 0.5*(1 - targetPower)),
                  pwr = targetPower)

    pwrs <- NULL
    n <- c(ceiling(n_rough)-1, ceiling(n_rough))

    for (j in 1:length(n)){
      n_temp <- n[j]

      x <- stats::qchisq(sob[,1], n_temp - 1)
      z <- stats::qnorm(sob[,2], mu, sigma*sqrt(1/n_temp))

      sdv <- sqrt(x*sigma^2/(n_temp - 1))*sqrt(1/n_temp)

      if (deltaU == Inf){
        thresUp <- rep(Inf, length(z))
      }
      else{
        thresUp <- (deltaU - z)/stats::qt(1-alpha, n_temp - 1)
      }

      if (deltaL == -Inf){
        thresLow <- rep(Inf, length(z))
      }
      else{
        thresLow <- (z - deltaL)/stats::qt(1-alpha, n_temp - 1)
      }

      thres <- pmin(thresUp, thresLow)

      pwrs <- c(pwrs, mean(ifelse(sdv <= thres,1,0)))
    }

    df_samps <- data.frame(n_plot = samps)

    n_plot <- NULL
    if (plot == TRUE){
    plot_pwr <- ggplot2::ggplot(df_samps, ggplot2::aes(x = n_plot)) +
      ggplot2::stat_ecdf(geom = "step", pad = FALSE, colour = cbbPalette[6], linewidth = 2) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 13)) +
      ggplot2::theme(axis.text.x =  ggplot2::element_text(size = 13)) +
      ggplot2::labs(title = "Approximated Power Curve") +
      ggplot2::labs(y = "Power", x=bquote(italic(n))) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,size=20,face="bold",
                                                          margin= ggplot2::margin(0,0,10,0))) +
      ggplot2::theme(axis.title.x = ggplot2::element_text(size = 16, margin= ggplot2::margin(10,0,0,0))) +
      ggplot2::theme(axis.title.y = ggplot2::element_text(size = 16, margin= ggplot2::margin(0,10,0,0))) +
      ggplot2::geom_segment(ggplot2::aes(x = n_rough, y = 0, xend = n_rough, yend = targetPower), linetype="dashed", color = "black")

    plot_min <- ggplot2::ggplot_build(plot_pwr)$layout$panel_params[[1]]$x$breaks[1]
    if (is.na(plot_min)){
      plot_min <- floor(ggplot2::ggplot_build(plot_pwr)$layout$panel_params[[1]]$x$continuous_range[1])
    }

    plot_pwr <- plot_pwr +
      ggplot2::geom_segment(ggplot2::aes(x = plot_min, y = targetPower, xend = n_rough, yend = targetPower), linetype="dashed", color = "black")

    print(plot_pwr)
    }

    type <- ifelse(is.finite(deltaL) & is.finite(deltaU), "a",
                   ifelse(!is.finite(deltaL), "b", "c"))

    METHOD <- switch(type, a = "Equivalence test power calculation (TOST)",
                     b = "Nonsuperiority test power calculation",
                     c = "Noninferiority test power calculation")

    NOTE <- paste0("As a check, power was estimated to be ", pwrs[1], " for n = ", n[1],".")

    results <- structure(list(n = n[2], mu = mu, sigma = sigma,
                              sig.level = alpha, power = pwrs[2], bounds = c(deltaL, deltaU),
                              note = NOTE,
                              method = METHOD,
                              samps = samps, seed = seed, sobol = sobol, design = "OneSample"), class = "power.en.test")
    return(results)
  }

  if (length(n) == 1){
    if (is.null(seed)){
      seed <- ceiling(1000*stats::runif(1))
    }
    sob <- qrng::sobol(2^(sobol+16), d = 2, randomize = "digital.shift", seed = seed)

    ## generate two sds and one mean
    x <- stats::qchisq(sob[,1], n  - 1)
    z <- stats::qnorm(sob[,2], mu, sigma*sqrt(1/n))

    sdv <- sqrt(x*sigma^2/(n - 1))*sqrt(1/n)

    if (deltaU == Inf){
      thresUp <- rep(Inf, length(z))
    }
    else{
      thresUp <- (deltaU - z)/stats::qt(1-alpha, n  - 1)
    }

    if (deltaL == -Inf){
      thresLow <- rep(Inf, length(z))
    }
    else{
      thresLow <- (z - deltaL)/stats::qt(1-alpha, n  - 1)
    }

    thres <- pmin(thresUp, thresLow)

    pwrs <- mean(ifelse(sdv <= thres,1,0))

    type <- ifelse(is.finite(deltaL) & is.finite(deltaU), "a",
                   ifelse(!is.finite(deltaL), "b", "c"))

    METHOD <- switch(type, a = "Equivalence test power calculation (TOST)",
                     b = "Nonsuperiority test power calculation",
                     c = "Noninferiority test power calculation")

    NOTE <- "mu is not in (deltaL, deltaU), so we compute type I error rate."

    if(mu >= deltaU | mu <= deltaL){
      results <- structure(list(n = n, mu = mu, sigma = sigma,
                                sig.level = alpha, type.I.error = pwrs, bounds = c(deltaL, deltaU),
                                note = NOTE,
                                method = METHOD), class = "power.en.test")
    }
    else{
      results <- structure(list(n = n, mu = mu, sigma = sigma,
                                sig.level = alpha, power = pwrs, bounds = c(deltaL, deltaU),
                                method = METHOD), class = "power.en.test")
    }

    return(results)
  }
}
