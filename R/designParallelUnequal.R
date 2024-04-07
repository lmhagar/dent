#' Power Calculations for Two-Group Parallel Designs with Unequal Variances
#'
#' Approximates the power of equivalence and noninferiority tests for two independent samples with unequal variances. One can either find sample sizes that achieve desired statistical power or estimate the power for given sample sizes.
#'
#' @param diff the anticipated difference between the group means (\eqn{\mu_{1} - \mu_{2}}).
#' @param sigma1 the anticipated within-group standard deviation for group 1 (\eqn{\sigma_{1}}).
#' @param sigma2 the anticipated within-group standard deviation for group 2 (\eqn{\sigma_{2}}).
#' @param deltaL the lower bound for the interval of equivalence (can be `-Inf` for noninferiority test).
#' @param deltaU the upper bound for the interval of equivalence (can be `Inf` for noninferiority test).
#' @param alpha the significance level \eqn{\alpha \in (0,1)}.
#' @param targetPower the desired statistical power of the equivalence or noninferiority test, must be a single number between 0 and 1 (exclusive). Exactly one of the following pairs must be specified: `targetPower` and `q` or `n1` and `n2`. Specify both `targetPower` and `q` to find sample sizes `n1` and `n2` that yield desired power such that \eqn{n_{2} \approx q n_{1}}.
#' @param q the multiplicative constant for imbalanced two-group sample size determination (\eqn{n_{2} = q n_{1}}). The default value is 1.
#' @param n1 the sample size for group 1, must be a single integer such that \eqn{n_{1} \ge 2}. Exactly one of the following pairs must be specified: `targetPower` and `q` or `n1` and `n2`. Specify both `n1` and `n2` to estimate statistical power for these sample sizes.
#' @param n2 the sample size for group 2, must be a single integer such that \eqn{n_{2} \ge 2}.
#' @param plot a logical variable indicating whether to return a plot of the power curve. If `n1` and `n2` are specified instead of `q` and `targetPower`, this variable is automatically set to `FALSE`. If you wish to approximate many power curves, suppressing the plots will expedite this process.
#' @param seed if provided, a single positive integer is used to ensure reproducibility when randomizing the Sobol' sequence via `sobol()` in the `qrng` package.
#' @param sobol one of the following integers: \eqn{s \in \{0, 1, 2, 3, 4 \}}. When approximating the power curve using `targetPower` and `q`, \eqn{2^{s + 10}} points are generated from the Sobol' sequence. When estimating power for given samples sizes `n1` and `n2`, \eqn{2^{s + 16}} points are generated from the Sobol' sequence. The default setting is \eqn{s = 0}, which ensures that each function call should take less than two seconds. As \eqn{s} increases, the sample size calculation is less sensitive to simulation variability but takes longer to complete. However, all function calls should still take less than 30 seconds when \eqn{s = 4}.
#'
#' @examples
#' # specify targetPower and q to obtain sample sizes n1 and n2
#' DesignParallelUnequal(diff = -4, sigma1 = 18, sigma2 = 15, deltaL = -19.2, deltaU = 19.2,
#' alpha = 0.05, targetPower = 0.8, q = 1, plot = TRUE, seed = 1, sobol = 0)
#'
#' # specify n1 and n2 to estimate power for this design
#' DesignParallelUnequal(diff = -4, sigma1 = 18, sigma2 = 15, deltaL = -19.2, deltaU = 19.2,
#' alpha = 0.05, n1 = 17, n2 = 17, plot = FALSE, seed = 1, sobol = 0)
#'
#' @return The sample sizes or power estimate are returned as a list with supplementary information. If `targetPower` and `q` are specified to find sample sizes `n1` and `n2`, a plot of the approximated power curve will appear on the plot pane if `plot = TRUE`. To find a sample size that corresponds to a different `targetPower`, save this function's output to an object and use the `UpdateTargetPower()` function.
#' @export
DesignParallelUnequal <- function(diff = NULL, sigma1 = NULL, sigma2 = NULL, deltaL = -Inf,
                                  deltaU = Inf, alpha = NULL, targetPower = NULL, q = 1, n1 = NULL,
                                  n2 = NULL, plot = TRUE, seed = NULL, sobol = 0){

  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  ## error checking
  if(!is.numeric(diff) | length(diff) != 1) {
    stop("Please specify a valid number for diff.")}
  if(!is.numeric(sigma1) | length(sigma1) != 1){
    stop("Please specify a valid number for sigma1.")}
  else if (sigma1 <= 0 | !is.finite(sigma1)){
    stop("Please specify a valid number for sigma1.")}
  if(!is.numeric(sigma2) | length(sigma2) != 1){
    stop("Please specify a valid number for sigma2.")}
  else if (sigma2 <= 0 | !is.finite(sigma2)){
    stop("Please specify a valid number for sigma2.")}
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
  if(sum(c((length(targetPower) != 1 | length(q) != 1), (length(n1) != 1 | length(n2) != 1))) != 1){
    stop("Please specify valid inputs for one of the following pairs: targetPower and q or n1 and n2.")}
  if(length(targetPower) == 1 & length(q) == 1){
    if(!is.numeric(targetPower)) {
      stop("Please specify a valid number for targetPower.")}
    # if(is.null(targetPower) | !is.numeric(targetPower)) {
    #   stop("Please specify a valid number for targetPower.")}
    if (is.numeric(targetPower)){
      if (targetPower <= 0 | targetPower >= 1){
        stop("Please specify a valid number for targetPower.")}
      if(!is.numeric(q)) {
        stop("Please specify a valid number for q.")}
      else if (is.numeric(q)){
        if (q <= 0) {
          stop("Please specify a valid number for q.")}
      }
      if (diff >= deltaU | diff <= deltaL){
        stop("Please ensure diff is between deltaL and deltaU.")
      }
    }
  }
  if(length(n1) == 1 & length(n2) == 1){
    if(!is.numeric(n1)) {
      stop("Please specify a valid integer for n1.")}
    else if (n1 < 2){
      stop("Please specify a valid integer for n1.")}
    else if (n1%%1 != 0){
      stop("Please specify a valid a integer for n1.")
    }
    if(!is.numeric(n2)) {
      stop("Please specify a valid integer for n2.")}
    else if (n2 < 2){
      stop("Please specify a valid integer for n2.")}
    else if (n2%%1 != 0){
      stop("Please specify a valid integer for n2.")
    }
  }
  if(!is.null(seed) & (!is.numeric(seed) | length(seed) != 1)) {
    stop("Please specify a valid seed for random number generation.")}
  else if (!is.null(seed)){
    if (seed%%1 != 0 | seed < 0){
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
  if (length(n1) == 1 & length(n2) == 1){
    plot <- FALSE
  }

  if (length(targetPower) == 1 & length(q) == 1){
    targetPowerfn <- function(u, deltaL, deltaU, diff, sigma1, sigma2, n_val, alpha, q){
      n1 <- n_val[1]; n2 <- max(2, q*n_val)
      x <- stats::qchisq(u[1], n1 - 1)
      y <- stats::qchisq(u[2], n2 - 1)
      z <- stats::qnorm(u[3], diff, sqrt(sigma1^2/n1 + sigma2^2/n2))

      sdv <- sqrt(sigma1^2*x/((n1 - 1)*n1) + sigma2^2*y/(n2*(n2 - 1)))

      dfw <- sdv^4/(sigma1^4*x^2/((n1 - 1)^3*n1^2) + sigma2^4*y^2/((n2 - 1)^3*n2^2))

      if (deltaU == Inf){
        thres <- (z - deltaL)/stats::qt(1-alpha, dfw)
      }
      else if (deltaL == -Inf){
        thres <- (deltaU - z)/stats::qt(1-alpha, dfw)
      }
      else{
        thresUp <- (deltaU - z)/stats::qt(1-alpha, dfw)
        thresLow <- (z - deltaL)/stats::qt(1-alpha, dfw)
        thres <- pmin(thresUp, thresLow)
      }

      return(thres - sdv)
    }

    ## this is a leaner inplementation of uniroot in base R
    uu <- function (fun, lower, upper, f_lower = NULL, f_upper = NULL, maxiter = 1000, tol = 1e-4, tol2 = 0.0005, ...)
    {
      f <- function(x) fun(x, ...)
      x1 <- lower
      if (!is.null(f_lower)){
        f1 <- f_lower
      }
      else{
        f1 <- f(x1)
      }
      if (f1 > 0){return(x1)}
      x2 <- upper
      if (!is.null(f_upper)){
        f2 <- f_upper
      }
      else{
        f2 <- f(x2)
      }
      f2 <- f(x2)
      if (f2 < 0){return(x2)}
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
        if ((upper - lower) < tol2 * max(abs(upper), 1)) {
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
        if (abs(x - x3) < tol2) {
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
    sob <- qrng::sobol(2^(sobol + 10), d = 3, randomize = "digital.shift", seed = seed)

    ## find starting point for root-finding algorithm using normal approximations
    if (!is.finite(deltaU)){
      mid_val <- ((stats::qnorm(targetPower) + stats::qnorm(1 - alpha))*sqrt(sigma1^2 + sigma2^2/q)/(diff - deltaL))^2
      upper_val <- ((stats::qnorm((3 +targetPower)/4) + stats::qnorm(1 - alpha))*sqrt(sigma1^2 + sigma2^2/q)/(diff - deltaL))^2
    }
    else if (!is.finite(deltaL)){
      mid_val <- ((stats::qnorm(targetPower) + stats::qnorm(1 - alpha))*sqrt(sigma1^2 + sigma2^2/q)/(deltaU - diff))^2
      upper_val <- ((stats::qnorm((3 +targetPower)/4) + stats::qnorm(1 - alpha))*sqrt(sigma1^2 + sigma2^2/q)/(deltaU - diff))^2
    }
    else{
      ## find starting point for root-finding algorithm using normal approximations to
      ## the t-distribution
      a_cons <- (deltaU - diff)/sqrt(sigma1^2 + sigma2^2/q)
      b_cons <- (deltaL - diff)/sqrt(sigma1^2 + sigma2^2/q)
      c_cons <- stats::qnorm(1-alpha)
      ## lower bound for root-finding algorithm
      lower_cons <- 2*c_cons*sqrt(sigma1^2 + sigma2^2/q)/(deltaU - deltaL)
      upper_cons <- lower_cons

      fn_ci = function(n_sq, a, b, c, pwr){
        return(stats::pnorm(a*n_sq - c) - stats::pnorm(b*n_sq + c) - pwr)}

      upper_large <- FALSE
      while(upper_large == FALSE){
        upper_cons <- 10*upper_cons
        upper_check <- fn_ci(n_sq = sqrt(upper_cons), a = a_cons, b = b_cons, c = c_cons, pwr = targetPower)
        if (upper_check > 0){
          upper_large <- TRUE
        }
      }

      ## mid_val should be close to the final sample size
      mid_val <- (uu(fn_ci, a = a_cons, b = b_cons, c = c_cons, pwr = targetPower,
                     lower = lower_cons, upper = upper_cons))^2

      upper_large <- FALSE
      upper_cons <- mid_val
      while(upper_large == FALSE){
        upper_cons <- 10*upper_cons
        upper_check <- fn_ci(n_sq = sqrt(upper_cons), a = a_cons, b = b_cons, c = c_cons, pwr = (3 + targetPower)/4)
        if (upper_check > 0){
          upper_large <- TRUE
        }
      }

      upper_val <- (uu(fn_ci, a = a_cons, b = b_cons, c = c_cons, pwr = (3 + targetPower)/4,
                       lower = sqrt(mid_val), upper = sqrt(upper_cons)))^2
    }

    ## upper_val and lower_val will be the next sample sizes explored by the root-finding algorithm
    ## depending on whether the point corresponds to the rejection region at mid_val
    mid_val <- max(mid_val, 10)
    lower_val <- 0.5*mid_val

    endpoints_vec <- rep(0, nrow(sob))
    samps <- NULL
    for (i in 1:nrow(sob)){
      power1 <- targetPowerfn(n_val = ceiling(mid_val),
                              deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                              u = sob[i,], alpha = alpha, q = q)
      if (power1 >= 0){
        power1b <- targetPowerfn(n_val = ceiling(lower_val),
                                 deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                                 u = sob[i,], alpha = alpha, q = q)
        if (power1b >= 0){
          samps[i] <- uu(targetPowerfn, lower =2,
                         upper = ceiling(lower_val), f_upper = power1b,
                         deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                         u = sob[i,], alpha = alpha, q = q)
        }
        else{
          samps[i] <- uu(targetPowerfn, lower = ceiling(lower_val), f_lower = power1b,
                         upper = ceiling(mid_val), f_upper = power1,
                         deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                         u = sob[i,], alpha = alpha, q = q)
        }
      }
      else{
        power2 <- targetPowerfn(n_val = ceiling(upper_val),
                                deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                                u = sob[i,], alpha = alpha, q = q)
        if (power2 >= 0){
          samps[i] <- uu(targetPowerfn, lower = ceiling(mid_val), f_lower = power1,
                         upper = ceiling(upper_val), f_upper = power2,
                         deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                         u = sob[i,], alpha = alpha, q = q)
        }
        else{
          endpoints_vec[i] <- 1
        }
      }
    }

    last_group <- which(endpoints_vec == 1)
    if (length(last_group) == 0){
      upper_c <- 2
    } else{
      upper_c <- 1
      while(length(last_group) > 0){
        if (upper_c > 32){
          last_group <- NULL
        }
        upper_c <- 2*upper_c
        endpoints1_vec <- NULL
        for (i in 1:length(last_group)){
          endpoints1_vec[i] <- targetPowerfn(n_val = ceiling(upper_c*upper_val),
                                             deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                                             u = sob[last_group[i],], alpha = alpha, q = q)
        }
        keep_vec <- ifelse(endpoints1_vec >= 0, FALSE, TRUE)
        ## only keep points that still do not satisfy power criterion after increasing
        ## the upper bound for the sample size
        last_group <- last_group[keep_vec]
      }
    }

    ## implement the root-finding algorithm for each point in the Sobol' sequence
    ## that required a large upper bound (i.e., those in last_group)
    for (i in 1:nrow(sob)){
      if (endpoints_vec[i] == 1){
        samps[i] <- uu(targetPowerfn, lower = ceiling(upper_val),
                       upper = ceiling(upper_c*upper_val),
                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                       u = sob[i,], alpha = alpha, q = q)
      }
    }

    n_rough <- stats::quantile(samps, targetPower)

    pwrs <- NULL
    n1 <- n_rough
    n2 <- pmax(2, q*n1)

    ## explore sampling distribution at the preliminary sample size recommendation
    n1_temp <- n1; n2_temp <- n2
    ## generate two sds and one mean
    x <- stats::qchisq(sob[,1], n1_temp - 1)
    y <- stats::qchisq(sob[,2], n2_temp - 1)
    z <- stats::qnorm(sob[,3], diff, sqrt(sigma1^2/n1_temp + sigma2^2/n2_temp))

    sdv <- sqrt(sigma1^2*x/((n1_temp - 1)*n1_temp) + sigma2^2*y/(n2_temp*(n2_temp - 1)))

    dfw <- sdv^4/(sigma1^4*x^2/((n1_temp - 1)^3*n1_temp^2) + sigma2^4*y^2/((n2_temp - 1)^3*n2_temp^2))

    if (deltaU == Inf){
      thresUp <- rep(Inf, length(z))
    }
    else{
      thresUp <- (deltaU - z)/stats::qt(1-alpha, dfw)
    }

    if (deltaL == -Inf){
      thresLow <- rep(Inf, length(z))
    }
    else{
      thresLow <- (z - deltaL)/stats::qt(1-alpha, dfw)
    }

    thres <- pmin(thresUp, thresLow)

    pwrs <- thres - sdv

    consistency <- ifelse(samps < n_rough, round(pwrs,2) >= 0, round(pwrs,2) <= 0)
    consistency <- ifelse(consistency == 1, 1, as.numeric(abs(pwrs) < 0.02))

    ## for any points where the root-finding algorithm has caused issues,
    ## re-run the root-finding algorithm starting at n_*
    inconsistent <- which(consistency == 0)
    if (length(inconsistent) > 0){
      for (i in 1:length(inconsistent)){
        if (pwrs[inconsistent[i]] < 0){
          samps[inconsistent[i]] <- uu(targetPowerfn, lower = n_rough, upper = ceiling(upper_c*upper_val),
                                       f_lower = pwrs[inconsistent[i]],
                                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                                       u = sob[inconsistent[i],], alpha = alpha, q = q)
        }
        else{
          samps[inconsistent[i]] <- uu(targetPowerfn, lower = 2, upper = n_rough,
                                       f_upper = pwrs[inconsistent[i]],
                                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                                       u = sob[inconsistent[i],], alpha = alpha, q = q)
        }
      }
    }
    funecdf <- stats::ecdf(samps)

    ecdf_root <- function(quant, pwr){return(funecdf(quant) - pwr)}
    n_rough <- uu(ecdf_root, lower = stats::quantile(samps, targetPower*0.5), upper = stats::quantile(samps, targetPower + 0.5*(1 - targetPower)),
                  pwr = targetPower)

    ## get confirmatory power estimate for output
    n1_temp <- ceiling(n_rough); n2_temp <- pmax(2, ceiling(q*n_rough))
    ## generate two sds and one mean
    x <- stats::qchisq(sob[,1], n1_temp - 1)
    y <- stats::qchisq(sob[,2], n2_temp - 1)
    z <- stats::qnorm(sob[,3], diff, sqrt(sigma1^2/n1_temp + sigma2^2/n2_temp))

    sdv <- sqrt(sigma1^2*x/((n1_temp - 1)*n1_temp) + sigma2^2*y/(n2_temp*(n2_temp - 1)))

    dfw <- sdv^4/(sigma1^4*x^2/((n1_temp - 1)^3*n1_temp^2) + sigma2^4*y^2/((n2_temp - 1)^3*n2_temp^2))

    if (deltaU == Inf){
      thresUp <- rep(Inf, length(z))
    }
    else{
      thresUp <- (deltaU - z)/stats::qt(1-alpha, dfw)
    }

    if (deltaL == -Inf){
      thresLow <- rep(Inf, length(z))
    }
    else{
      thresLow <- (z - deltaL)/stats::qt(1-alpha, dfw)
    }

    thres <- pmin(thresUp, thresLow)

    pwrs <- thres - sdv

    df_samps <- data.frame(n_plot = samps)

    n_plot <- NULL
    if (plot == TRUE){
      plot_pwr <- ggplot2::ggplot(df_samps, ggplot2::aes(x = n_plot)) + ggplot2::theme_bw() +
        ggplot2::stat_ecdf(geom = "step", pad = FALSE, colour = cbbPalette[6], linewidth = 2) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 13)) +
        ggplot2::theme(axis.text.x =  ggplot2::element_text(size = 13)) +
        ggplot2::labs(title = "Approximated Power Curve") +
        ggplot2::labs(y = "Power", x=bquote(italic(n)[1]*"  ("*italic(n)[2]*" = "*.(round(q,3))*italic(n)[1]*")")) +
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
                     b = "Noninferiority test power calculation (for Group 2)",
                     c = "Noninferiority test power calculation (for Group 1)")

    results <- structure(list(n1 = n1_temp, n2 = n2_temp, q = q, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                              sig.level = alpha, power = mean(pwrs >= 0), bounds = c(deltaL, deltaU),
                              method = METHOD, upper.bound = ceiling(upper_c*upper_val),
                              samps = samps, seed = seed, sobol = sobol, design = "ParallelUnequal"), class = "power.en.test")

    if (sigma1 == sigma2){
      message("This function works with equal population variances, but did you mean to use DesignParallelEqual()?")
    }
    return(results)
  }

  if (length(n1) == 1 & length(n2) == 1){
    if (is.null(seed)){
      seed <- ceiling(1000*stats::runif(1))
    }
    sob <- qrng::sobol(2^(sobol+16), d = 3, randomize = "digital.shift", seed = seed)

    ## generate two sds and one mean
    x <- stats::qchisq(sob[,1], n1 - 1)
    y <- stats::qchisq(sob[,2], n2 - 1)
    z <- stats::qnorm(sob[,3], diff, sqrt(sigma1^2/n1 + sigma2^2/n2))

    sdv <- sqrt(sigma1^2*x/((n1 - 1)*n1) + sigma2^2*y/(n2*(n2 - 1)))

    dfw <- sdv^4/(sigma1^4*x^2/((n1 - 1)^3*n1^2) + sigma2^4*y^2/((n2 - 1)^3*n2^2))

    if (deltaU == Inf){
      thresUp <- rep(Inf, length(z))
    }
    else{
      thresUp <- (deltaU - z)/stats::qt(1-alpha, dfw)
    }

    if (deltaL == -Inf){
      thresLow <- rep(Inf, length(z))
    }
    else{
      thresLow <- (z - deltaL)/stats::qt(1-alpha, dfw)
    }

    thres <- pmin(thresUp, thresLow)

    pwrs <- mean(ifelse(sdv <= thres,1,0))

    type <- ifelse(is.finite(deltaL) & is.finite(deltaU), "a",
                   ifelse(!is.finite(deltaL), "b", "c"))

    METHOD <- switch(type, a = "Equivalence test power calculation (TOST)",
                     b = "Noninferiority test power calculation (for Group 2)",
                     c = "Noninferiority test power calculation (for Group 1)")

    NOTE <- "diff is not in (deltaL, deltaU), so we compute the type I error rate."

    if(diff >= deltaU | diff <= deltaL){
      results <- structure(list(n1 = n1, n2 = n2, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                                sig.level = alpha, type.I.error = pwrs, bounds = c(deltaL, deltaU),
                                note = NOTE,
                                method = METHOD), class = "power.en.test")
    }
    else{
      results <- structure(list(n1 = n1, n2 = n2, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                                sig.level = alpha, power = pwrs, bounds = c(deltaL, deltaU),
                                method = METHOD), class = "power.en.test")
    }

    if (sigma1 == sigma2){
      message("This function works with equal population variances, but did you mean to use DesignParallelEqual()?")
    }
    return(results)
  }
}

#' @title Printing power analysis for equivalence and noninferiority tests
#' @name power.en.test-methods
#'
#' @description Helper function to print power.en.test objects
#'
#' @param x an power.en.test object
#' @param ... further arguments to be passed through (see `print()` function)
#'
#' @rdname power.en.test-methods
#' @method print power.en.test
#' @keywords internal
#' @export
print.power.en.test <- function(x, ...){

  names_remove <- c("samps", "seed", "sobol", "design", "upper.bound")

  inds_remove <- which(names_remove %in% names(x))

  if (length(inds_remove) > 0){
    purge <- which(unlist(lapply(lapply(names(x), function(z) which(z == names_remove)), function(y) as.numeric((length(y) > 0)))) > 0)
    output_temp <- structure(x[-c(purge)], class = "power.htest")

  } else {
    output_temp <- structure(x, class = "power.htest")
  }


  # if (!(is.null(x$q))){
  #   output_temp <- structure(list(n1 = x$n1, n2 = x$n2, q = x$q, diff = x$diff, sigma1 = x$sigma1, sigma2 = x$sigma2,
  #                           sig.level = x$sig.level, power = x$power, bounds = x$bounds,
  #                           note = x$note,
  #                           method = x$method), class = "power.htest")
  # } else if (!(is.null(x$note))) {
  #   output_temp <- structure(list(n1 = x$n1, n2 = x$n2, diff = x$diff, sigma1 = x$sigma1, sigma2 = x$sigma2,
  #                                 sig.level = x$sig.level, type.I.error = x$type.I.error, bounds = x$bounds,
  #                                 note = x$note,
  #                                 method = x$method), class = "power.htest")
  # } else{
  #   output_temp <- structure(list(n1 = x$n1, n2 = x$n2, diff = x$diff, sigma1 = x$sigma1, sigma2 = x$sigma2,
  #                                 sig.level = x$sig.level, power = x$power, bounds = x$bounds,
  #                                 method = x$method), class = "power.htest")
  # }

  print(output_temp)
}

#' Updated Power Calculations for Equivalence and Noninferiority Tests
#'
#' Performs sample size calculations for equivalence and noninferiority tests using output saved from a previous function call with an updated target for statistical power.
#'
#' @param power.en.test an object to which the output from a `Design...()` function call was saved. The input `targetPower` must have been specified in this function call so that the power curve was approximated. This function will return an error if the original `Design...()` function call specified `n`, `n1`, or `n2`. All other inputs from the initial `Design()` function call will be passed to this function.
#' @param targetPower the desired statistical power of the equivalence or noninferiority test, must be a single number between 0 and 1 (exclusive).
#' @param plot a logical variable indicating whether to return a plot of the power curve.
#'
#' @examples
#' # save output from `Design...()` function call with previous `targetPower` specification
#' # to object `design.result`.
#' design.result <- DesignParallelUnequal(diff = -4, sigma1 = 15, sigma2 = 18, deltaL = -19.2,
#'                  deltaU = 19.2, alpha = 0.05, targetPower = 0.8, q = 1, seed = 1, sobol = 0)
#'
#' UpdateTargetPower(power.en.test = design.result, targetPower = 0.7, plot = TRUE)
#'
#' @return The sample sizes or power estimate are returned as a list with supplementary information. A plot of the approximated power curve will appear on the plot pane if `plot = TRUE`.
#' @export
UpdateTargetPower <- function(power.en.test = NULL, targetPower = NULL, plot = TRUE){

  if(!(methods::is(power.en.test, "power.en.test"))){
    stop("Please input valid power.en.test object.")
  }

  deltaL <- power.en.test$bounds[1]
  deltaU <- power.en.test$bounds[2]
  alpha <- power.en.test$sig.level
  seed <- power.en.test$seed
  sobol <- power.en.test$sobol
  design <- power.en.test$design
  upper.bound <- power.en.test$upper.bound

  if (length(design) != 1){
    stop("Please input a power.en.test object generated by specifying targetPower with a valid design type.")
  } else if (!(design %in% c("ParallelEqual", "ParallelUnequal", "OneSample", "PairedSample",
                      "Crossover2x2Unequal", "Crossover2x2Equal", "CrossoverDualUnequal",
                      "CrossoverDualEqualComp", "CrossoverDualEqualNoComp"))){
    stop("Please input a power.en.test object generated by specifying targetPower with a valid design type.")
  }

  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  ## error checking
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
  if(!is.null(seed) & (!is.numeric(seed) | length(seed) != 1)) {
    stop("Please specify a valid seed for random number generation.")}
  else if (!is.null(seed)){
    if (seed%%1 != 0 | seed < 0){
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

  if(!is.numeric(upper.bound)){
    stop("Please specify a valid number for the upper bound.")
  }

  uu <- function (fun, lower, upper, f_lower = NULL, f_upper = NULL, maxiter = 1000, tol = 1e-4, tol2 = 0.0005, ...)
  {
    f <- function(x) fun(x, ...)
    x1 <- lower
    if (!is.null(f_lower)){
      f1 <- f_lower
    }
    else{
      f1 <- f(x1)
    }
    if (f1 > 0){return(x1)}
    x2 <- upper
    if (!is.null(f_upper)){
      f2 <- f_upper
    }
    else{
      f2 <- f(x2)
    }
    f2 <- f(x2)
    if (f2 < 0){return(x2)}
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
      if ((upper - lower) < tol2 * max(abs(upper), 1)) {
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
      if (abs(x - x3) < tol2) {
        x0 <- x
        return(x0)
      }
      x3 <- x
    }
    return(x0)
  }

  if (design == "ParallelUnequal"){

    diff <- power.en.test$diff
    sigma1 <- power.en.test$sigma1
    sigma2 <- power.en.test$sigma2
    q <- power.en.test$q

    if(!is.numeric(diff) | length(diff) != 1) {
      stop("Please specify a valid number for diff.")}
    if(!is.numeric(sigma1) | length(sigma1) != 1){
      stop("Please specify a valid number for sigma1.")}
    else if (sigma1 <= 0 | !is.finite(sigma1)){
      stop("Please specify a valid number for sigma1.")}
    if(!is.numeric(sigma2) | length(sigma2) != 1){
      stop("Please specify a valid number for sigma2.")}
    else if (sigma2 <= 0 | !is.finite(sigma2)){
      stop("Please specify a valid number for sigma2.")}

    if(sum(c((length(targetPower) != 1 | length(q) != 1)))){
      stop("Please specify a valid input for targetPower.")}
    if(length(targetPower) == 1 & length(q) == 1){
      if(!is.numeric(targetPower)) {
        stop("Please specify a valid number for targetPower.")}
      # if(is.null(targetPower) | !is.numeric(targetPower)) {
      #   stop("Please specify a valid number for targetPower.")}
      if (is.numeric(targetPower)){
        if (targetPower <= 0 | targetPower >= 1){
          stop("Please specify a valid number for targetPower.")}
        if(!is.numeric(q)) {
          stop("Please specify a valid number for q.")}
        else if (is.numeric(q)){
          if (q <= 0) {
            stop("Please specify a valid number for q.")}
        }
        if (diff >= deltaU | diff <= deltaL){
          stop("Please ensure diff is between deltaL and deltaU.")
        }
      }
    }

    if (is.null(seed)){
      seed <- ceiling(1000*stats::runif(1))
    }
    sob <- qrng::sobol(2^(sobol + 10), d = 3, randomize = "digital.shift", seed = seed)

    samps <- power.en.test$samps
    funecdf <- stats::ecdf(samps)

    ecdf_root <- function(quant, pwr){return(funecdf(quant) - pwr)}
    n_rough <- uu(ecdf_root, lower = stats::quantile(samps, targetPower*0.5), upper = stats::quantile(samps, targetPower + 0.5*(1 - targetPower)),
                  pwr = targetPower)

    pwrs <- NULL
    n1 <- n_rough
    n2 <- pmax(2, q*n1)

    ## explore sampling distribution at the preliminary sample size recommendation
    n1_temp <- n1; n2_temp <- n2
    ## generate two sds and one mean
    x <- stats::qchisq(sob[,1], n1_temp - 1)
    y <- stats::qchisq(sob[,2], n2_temp - 1)
    z <- stats::qnorm(sob[,3], diff, sqrt(sigma1^2/n1_temp + sigma2^2/n2_temp))

    sdv <- sqrt(sigma1^2*x/((n1_temp - 1)*n1_temp) + sigma2^2*y/(n2_temp*(n2_temp - 1)))

    dfw <- sdv^4/(sigma1^4*x^2/((n1_temp - 1)^3*n1_temp^2) + sigma2^4*y^2/((n2_temp - 1)^3*n2_temp^2))

    if (deltaU == Inf){
      thresUp <- rep(Inf, length(z))
    }
    else{
      thresUp <- (deltaU - z)/stats::qt(1-alpha, dfw)
    }

    if (deltaL == -Inf){
      thresLow <- rep(Inf, length(z))
    }
    else{
      thresLow <- (z - deltaL)/stats::qt(1-alpha, dfw)
    }

    thres <- pmin(thresUp, thresLow)

    pwrs <- thres - sdv

    consistency <- ifelse(samps < n_rough, round(pwrs,2) >= 0, round(pwrs,2) <= 0)
    consistency <- ifelse(consistency == 1, 1, as.numeric(abs(pwrs) < 0.02))

    ## for any points where the root-finding algorithm has caused issues,
    ## re-run the root-finding algorithm starting at n_*
    inconsistent <- which(consistency == 0)
    if (length(inconsistent) > 0){

      targetPowerfn <- function(u, deltaL, deltaU, diff, sigma1, sigma2, n_val, alpha, q){
        n1 <- n_val[1]; n2 <- max(2, q*n_val)
        x <- stats::qchisq(u[1], n1 - 1)
        y <- stats::qchisq(u[2], n2 - 1)
        z <- stats::qnorm(u[3], diff, sqrt(sigma1^2/n1 + sigma2^2/n2))

        sdv <- sqrt(sigma1^2*x/((n1 - 1)*n1) + sigma2^2*y/(n2*(n2 - 1)))

        dfw <- sdv^4/(sigma1^4*x^2/((n1 - 1)^3*n1^2) + sigma2^4*y^2/((n2 - 1)^3*n2^2))

        if (deltaU == Inf){
          thres <- (z - deltaL)/stats::qt(1-alpha, dfw)
        }
        else if (deltaL == -Inf){
          thres <- (deltaU - z)/stats::qt(1-alpha, dfw)
        }
        else{
          thresUp <- (deltaU - z)/stats::qt(1-alpha, dfw)
          thresLow <- (z - deltaL)/stats::qt(1-alpha, dfw)
          thres <- pmin(thresUp, thresLow)
        }

        return(thres - sdv)
      }

      for (i in 1:length(inconsistent)){
        if (pwrs[inconsistent[i]] < 0){
          if (n_rough > 0.5*upper.bound){
            upper.bound <- 2*n_rough
          }
          samps[inconsistent[i]] <- uu(targetPowerfn, lower = n_rough, upper = upper.bound,
                                       f_lower = pwrs[inconsistent[i]],
                                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                                       u = sob[inconsistent[i],], alpha = alpha, q = q)
        }
        else{
          samps[inconsistent[i]] <- uu(targetPowerfn, lower = 2, upper = n_rough,
                                       f_upper = pwrs[inconsistent[i]],
                                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                                       u = sob[inconsistent[i],], alpha = alpha, q = q)
        }
      }
    }
    funecdf <- stats::ecdf(samps)

    ecdf_root <- function(quant, pwr){return(funecdf(quant) - pwr)}
    n_rough <- uu(ecdf_root, lower = stats::quantile(samps, targetPower*0.5), upper = stats::quantile(samps, targetPower + 0.5*(1 - targetPower)),
                  pwr = targetPower)

    ## get confirmatory power estimate for output
    n1_temp <- ceiling(n_rough); n2_temp <- pmax(2, ceiling(q*n_rough))
    ## generate two sds and one mean
    x <- stats::qchisq(sob[,1], n1_temp - 1)
    y <- stats::qchisq(sob[,2], n2_temp - 1)
    z <- stats::qnorm(sob[,3], diff, sqrt(sigma1^2/n1_temp + sigma2^2/n2_temp))

    sdv <- sqrt(sigma1^2*x/((n1_temp - 1)*n1_temp) + sigma2^2*y/(n2_temp*(n2_temp - 1)))

    dfw <- sdv^4/(sigma1^4*x^2/((n1_temp - 1)^3*n1_temp^2) + sigma2^4*y^2/((n2_temp - 1)^3*n2_temp^2))

    if (deltaU == Inf){
      thresUp <- rep(Inf, length(z))
    }
    else{
      thresUp <- (deltaU - z)/stats::qt(1-alpha, dfw)
    }

    if (deltaL == -Inf){
      thresLow <- rep(Inf, length(z))
    }
    else{
      thresLow <- (z - deltaL)/stats::qt(1-alpha, dfw)
    }

    thres <- pmin(thresUp, thresLow)

    pwrs <- thres - sdv

    df_samps <- data.frame(n_plot = samps)

    n_plot <- NULL
    if (plot == TRUE){
      plot_pwr <- ggplot2::ggplot(df_samps, ggplot2::aes(x = n_plot)) + ggplot2::theme_bw() +
        ggplot2::stat_ecdf(geom = "step", pad = FALSE, colour = cbbPalette[6], linewidth = 2) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 13)) +
        ggplot2::theme(axis.text.x =  ggplot2::element_text(size = 13)) +
        ggplot2::labs(title = "Approximated Power Curve") +
        ggplot2::labs(y = "Power", x=bquote(italic(n)[1]*"  ("*italic(n)[2]*" = "*.(round(q,3))*italic(n)[1]*")")) +
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
                     b = "Noninferiority test power calculation (for Group 2)",
                     c = "Noninferiority test power calculation (for Group 1)")

    results <- structure(list(n1 = n1_temp, n2 = n2_temp, q = q, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                              sig.level = alpha, power = mean(pwrs >= 0), bounds = c(deltaL, deltaU),
                              method = METHOD, upper.bound = upper.bound,
                              samps = samps, seed = seed, sobol = sobol, design = "ParallelUnequal"), class = "power.en.test")

    if (sigma1 == sigma2){
      message("This function works with equal population variances, but did you mean to use DesignParallelEqual()?")
    }
    return(results)
  }

  if (design == "ParallelEqual"){

    diff <- power.en.test$diff
    sigma <- power.en.test$sigma
    q <- power.en.test$q

    if(!is.numeric(diff) | length(diff) != 1) {
      stop("Please specify a valid number for diff.")}

    if(!is.numeric(sigma) | length(sigma) != 1 ){
      stop("Please specify a valid number for sigma.")}
    else if (sigma <= 0 | !is.finite(sigma)){
      stop("Please specify a valid number for sigma.")}

    if(sum(c((length(targetPower) != 1 | length(q) != 1)))){
      stop("Please specify a valid input for targetPower.")}
    if(length(targetPower) == 1 & length(q) == 1){
      if(!is.numeric(targetPower)) {
        stop("Please specify a valid number for targetPower.")}
      if (is.numeric(targetPower)){
        if (targetPower <= 0 | targetPower >= 1){
          stop("Please specify a valid number for targetPower.")}
        if(!is.numeric(q)) {
          stop("Please specify a valid number for q.")}
        else if (is.numeric(q)){
          if (q <= 0) {
            stop("Please specify a valid number for q.")}
        }
        if (diff >= deltaU | diff <= deltaL){
          stop("Please ensure diff is between deltaL and deltaU.")
        }
      }
    }

    if (is.null(seed)){
      seed <- ceiling(1000*stats::runif(1))
    }
    sob <- qrng::sobol(2^(sobol + 10), d = 2, randomize = "digital.shift", seed = seed)

    samps <- power.en.test$samps
    funecdf <- stats::ecdf(samps)

    ecdf_root <- function(quant, pwr){return(funecdf(quant) - pwr)}
    n_rough <- uu(ecdf_root, lower = stats::quantile(samps, targetPower*0.5), upper = stats::quantile(samps, targetPower + 0.5*(1 - targetPower)),
                  pwr = targetPower)

    pwrs <- NULL
    n1 <- n_rough
    n2 <- pmax(2, q*n1)

    ## explore sampling distribution at the preliminary sample size recommendation
    n1_temp <- n1; n2_temp <- n2

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

    pwrs <- thres - sdv

    consistency <- ifelse(samps < n_rough, round(pwrs,2) >= 0, round(pwrs,2) <= 0)
    consistency <- ifelse(consistency == 1, 1, as.numeric(abs(pwrs) < 0.02))

    ## for any points where the root-finding algorithm has caused issues,
    ## re-run the root-finding algorithm starting at n_*
    inconsistent <- which(consistency == 0)
    if (length(inconsistent) > 0){

      targetPowerfnEq <- function(u, deltaL, deltaU, diff, sigma, n_val, alpha, q){
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

      for (i in 1:length(inconsistent)){
        if (pwrs[inconsistent[i]] < 0){
          if (n_rough > 0.5*upper.bound){
            upper.bound <- 2*n_rough
          }
          samps[inconsistent[i]] <- uu(targetPowerfnEq, lower = n_rough, upper = upper.bound,
                                       f_lower = pwrs[inconsistent[i]],
                                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma = sigma,
                                       u = sob[inconsistent[i],], alpha = alpha, q = q)
        }
        else{
          samps[inconsistent[i]] <- uu(targetPowerfnEq, lower = 2, upper = n_rough,
                                       f_upper = pwrs[inconsistent[i]],
                                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma = sigma,
                                       u = sob[inconsistent[i],], alpha = alpha, q = q)
        }
      }
    }
    funecdf <- stats::ecdf(samps)

    ecdf_root <- function(quant, pwr){return(funecdf(quant) - pwr)}
    n_rough <- uu(ecdf_root, lower = stats::quantile(samps, targetPower*0.5), upper = stats::quantile(samps, targetPower + 0.5*(1 - targetPower)),
                  pwr = targetPower)

    ## get confirmatory power estimate for output
    n1_temp <- ceiling(n_rough); n2_temp <- pmax(2, ceiling(q*n_rough))

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

    pwrs <- thres - sdv

    df_samps <- data.frame(n_plot = samps)

    n_plot <- NULL
    if (plot == TRUE){
      plot_pwr <- ggplot2::ggplot(df_samps, ggplot2::aes(x = n_plot)) + ggplot2::theme_bw() +
        ggplot2::stat_ecdf(geom = "step", pad = FALSE, colour = cbbPalette[6], linewidth = 2) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 13)) +
        ggplot2::theme(axis.text.x =  ggplot2::element_text(size = 13)) +
        ggplot2::labs(title = "Approximated Power Curve") +
        ggplot2::labs(y = "Power", x=bquote(italic(n)[1]*"  ("*italic(n)[2]*" = "*.(round(q,3))*italic(n)[1]*")")) +
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
                     b = "Noninferiority test power calculation (for Group 2)",
                     c = "Noninferiority test power calculation (for Group 1)")

    results <- structure(list(n1 = n1_temp, n2 = n2_temp, q = q, diff = diff, sigma = sigma,
                              sig.level = alpha, power = mean(pwrs >= 0), bounds = c(deltaL, deltaU),
                              method = METHOD, upper.bound = upper.bound,
                              samps = samps, seed = seed, sobol = sobol, design = "ParallelEqual"), class = "power.en.test")
    return(results)
  }

  if (design == "OneSample"){

    mu <- power.en.test$mu
    sigma <- power.en.test$sigma

    if(!is.numeric(mu) | length(mu) != 1) {
      stop("Please specify a valid number for mu.")}

    if(!is.numeric(sigma) | length(sigma) != 1){
      stop("Please specify a valid number for sigma.")}
    else if (sigma <= 0 | !is.finite(sigma)){
      stop("Please specify a valid number for sigma.")}

    if(length(targetPower) != 1){
      stop("Please specify a valid input for targetPower.")}
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

    if (is.null(seed)){
      seed <- ceiling(1000*stats::runif(1))
    }
    sob <- qrng::sobol(2^(sobol + 10), d = 2, randomize = "digital.shift", seed = seed)

    samps <- power.en.test$samps
    funecdf <- stats::ecdf(samps)

    ecdf_root <- function(quant, pwr){return(funecdf(quant) - pwr)}
    n_rough <- uu(ecdf_root, lower = stats::quantile(samps, targetPower*0.5), upper = stats::quantile(samps, targetPower + 0.5*(1 - targetPower)),
                  pwr = targetPower)

    pwrs <- NULL
    n <- n_rough

    n_temp <- n_rough

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

    pwrs <- thres - sdv

    consistency <- ifelse(samps < n_rough, round(pwrs,2) >= 0, round(pwrs,2) <= 0)
    consistency <- ifelse(consistency == 1, 1, as.numeric(abs(pwrs) < 0.02))

    ## for any points where the root-finding algorithm has caused issues,
    ## re-run the root-finding algorithm starting at n_*
    inconsistent <- which(consistency == 0)
    if (length(inconsistent) > 0){

      targetPowerfnOne <- function(u, deltaL, deltaU, mu, sigma, n_val, alpha){
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

      for (i in 1:length(inconsistent)){
        if (pwrs[inconsistent[i]] < 0){
          if (n_rough > 0.5*upper.bound){
            upper.bound <- 2*n_rough
          }
          samps[inconsistent[i]] <- uu(targetPowerfnOne, lower = n_rough, upper = upper.bound,
                                       f_lower = pwrs[inconsistent[i]],
                                       deltaL = deltaL, deltaU = deltaU, mu = mu, sigma = sigma,
                                       u = sob[inconsistent[i],], alpha = alpha)
        }
        else{
          samps[inconsistent[i]] <- uu(targetPowerfnOne, lower = 2, upper = n_rough,
                                       f_upper = pwrs[inconsistent[i]],
                                       deltaL = deltaL, deltaU = deltaU, mu = mu, sigma = sigma,
                                       u = sob[inconsistent[i],], alpha = alpha)
        }
      }
    }
    funecdf <- stats::ecdf(samps)

    ecdf_root <- function(quant, pwr){return(funecdf(quant) - pwr)}
    n_rough <- uu(ecdf_root, lower = stats::quantile(samps, targetPower*0.5), upper = stats::quantile(samps, targetPower + 0.5*(1 - targetPower)),
                  pwr = targetPower)

    ## get confirmatory power estimate for output
    n_temp <- ceiling(n_rough)

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

    pwrs <- thres - sdv

    df_samps <- data.frame(n_plot = samps)

    n_plot <- NULL
    if (plot == TRUE){
      plot_pwr <- ggplot2::ggplot(df_samps, ggplot2::aes(x = n_plot)) + ggplot2::theme_bw() +
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

    results <- structure(list(n = n_temp, mu = mu, sigma = sigma,
                              sig.level = alpha, power = mean(pwrs >= 0), bounds = c(deltaL, deltaU),
                              upper.bound = upper.bound,
                              method = METHOD,
                              samps = samps, seed = seed, sobol = sobol, design = "OneSample"), class = "power.en.test")
    return(results)
  }

  if (design == "PairedSample"){

    diff <- power.en.test$diff
    sigma <- power.en.test$sigma

    if(!is.numeric(diff) | length(diff) != 1) {
      stop("Please specify a valid number for diff.")}

    if(!is.numeric(sigma) | length(sigma) != 1){
      stop("Please specify a valid number for sigma.")}
    else if (sigma <= 0 | !is.finite(sigma)){
      stop("Please specify a valid number for sigma.")}

    if(length(targetPower) != 1){
      stop("Please specify a valid input for targetPower.")}
    if(length(targetPower) == 1){
      if(!is.numeric(targetPower)) {
        stop("Please specify a valid number for targetPower.")}
      # if(is.null(targetPower) | !is.numeric(targetPower)) {
      #   stop("Please specify a valid number for targetPower.")}
      if (is.numeric(targetPower)){
        if (targetPower <= 0 | targetPower >= 1){
          stop("Please specify a valid number for targetPower.")}
        if (diff >= deltaU | diff <= deltaL){
          stop("Please ensure diff is between deltaL and deltaU.")
        }
      }
    }

    if (is.null(seed)){
      seed <- ceiling(1000*stats::runif(1))
    }
    sob <- qrng::sobol(2^(sobol + 10), d = 2, randomize = "digital.shift", seed = seed)

    samps <- power.en.test$samps
    funecdf <- stats::ecdf(samps)

    ecdf_root <- function(quant, pwr){return(funecdf(quant) - pwr)}
    n_rough <- uu(ecdf_root, lower = stats::quantile(samps, targetPower*0.5), upper = stats::quantile(samps, targetPower + 0.5*(1 - targetPower)),
                  pwr = targetPower)

    pwrs <- NULL
    n <- n_rough

    n_temp <- n_rough

    x <- stats::qchisq(sob[,1], n_temp - 1)
    z <- stats::qnorm(sob[,2], diff, sigma*sqrt(1/n_temp))

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

    pwrs <- thres - sdv

    consistency <- ifelse(samps < n_rough, round(pwrs,2) >= 0, round(pwrs,2) <= 0)
    consistency <- ifelse(consistency == 1, 1, as.numeric(abs(pwrs) < 0.02))

    ## for any points where the root-finding algorithm has caused issues,
    ## re-run the root-finding algorithm starting at n_*
    inconsistent <- which(consistency == 0)
    if (length(inconsistent) > 0){

      targetPowerfnPair <- function(u, deltaL, deltaU, diff, sigma, n_val, alpha){
        n <- n_val
        x <- stats::qchisq(u[1], n - 1)
        z <- stats::qnorm(u[2], diff, sigma*sqrt(1/n))

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

      for (i in 1:length(inconsistent)){
        if (pwrs[inconsistent[i]] < 0){
          if (n_rough > 0.5*upper.bound){
            upper.bound <- 2*n_rough
          }
          samps[inconsistent[i]] <- uu(targetPowerfnPair, lower = n_rough, upper = upper.bound,
                                       f_lower = pwrs[inconsistent[i]],
                                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma = sigma,
                                       u = sob[inconsistent[i],], alpha = alpha)
        }
        else{
          samps[inconsistent[i]] <- uu(targetPowerfnPair, lower = 2, upper = n_rough,
                                       f_upper = pwrs[inconsistent[i]],
                                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma = sigma,
                                       u = sob[inconsistent[i],], alpha = alpha)
        }
      }
    }
    funecdf <- stats::ecdf(samps)

    ecdf_root <- function(quant, pwr){return(funecdf(quant) - pwr)}
    n_rough <- uu(ecdf_root, lower = stats::quantile(samps, targetPower*0.5), upper = stats::quantile(samps, targetPower + 0.5*(1 - targetPower)),
                  pwr = targetPower)

    ## get confirmatory power estimate for output
    n_temp <- ceiling(n_rough)

    x <- stats::qchisq(sob[,1], n_temp - 1)
    z <- stats::qnorm(sob[,2], diff, sigma*sqrt(1/n_temp))

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

    pwrs <- thres - sdv

    df_samps <- data.frame(n_plot = samps)

    n_plot <- NULL
    if (plot == TRUE){
      plot_pwr <- ggplot2::ggplot(df_samps, ggplot2::aes(x = n_plot)) + ggplot2::theme_bw() +
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
                     b = "Noninferiority test power calculation (for Group 2)",
                     c = "Noninferiority test power calculation (for Group 1)")

    results <- structure(list(n = n_temp, diff = diff, sigma = sigma,
                              sig.level = alpha, power = mean(pwrs >= 0), bounds = c(deltaL, deltaU),
                              upper.bound = upper.bound,
                              method = METHOD,
                              samps = samps, seed = seed, sobol = sobol, design = "PairedSample"), class = "power.en.test")
    return(results)
  }

  if (design == "Crossover2x2Unequal"){

    diff <- power.en.test$diff
    sigma1 <- power.en.test$sigma1/2
    sigma2 <- power.en.test$sigma2/2
    q <- power.en.test$q

    if(!is.numeric(diff) | length(diff) != 1) {
      stop("Please specify a valid number for diff.")}
    if(!is.numeric(sigma1) | length(sigma1) != 1){
      stop("Please specify a valid number for sigma1.")}
    else if (sigma1 <= 0 | !is.finite(sigma1)){
      stop("Please specify a valid number for sigma1.")}
    if(!is.numeric(sigma2) | length(sigma2) != 1){
      stop("Please specify a valid number for sigma2.")}
    else if (sigma2 <= 0 | !is.finite(sigma2)){
      stop("Please specify a valid number for sigma2.")}

    if(sum(c((length(targetPower) != 1 | length(q) != 1)))){
      stop("Please specify a valid input for targetPower.")}
    if(length(targetPower) == 1 & length(q) == 1){
      if(!is.numeric(targetPower)) {
        stop("Please specify a valid number for targetPower.")}
      # if(is.null(targetPower) | !is.numeric(targetPower)) {
      #   stop("Please specify a valid number for targetPower.")}
      if (is.numeric(targetPower)){
        if (targetPower <= 0 | targetPower >= 1){
          stop("Please specify a valid number for targetPower.")}
        if(!is.numeric(q)) {
          stop("Please specify a valid number for q.")}
        else if (is.numeric(q)){
          if (q <= 0) {
            stop("Please specify a valid number for q.")}
        }
        if (diff >= deltaU | diff <= deltaL){
          stop("Please ensure diff is between deltaL and deltaU.")
        }
      }
    }

    if (is.null(seed)){
      seed <- ceiling(1000*stats::runif(1))
    }
    sob <- qrng::sobol(2^(sobol + 10), d = 3, randomize = "digital.shift", seed = seed)

    samps <- power.en.test$samps
    funecdf <- stats::ecdf(samps)

    ecdf_root <- function(quant, pwr){return(funecdf(quant) - pwr)}
    n_rough <- uu(ecdf_root, lower = stats::quantile(samps, targetPower*0.5), upper = stats::quantile(samps, targetPower + 0.5*(1 - targetPower)),
                  pwr = targetPower)

    pwrs <- NULL
    n1 <- n_rough
    n2 <- pmax(2, q*n1)

    ## explore sampling distribution at the preliminary sample size recommendation
    n1_temp <- n1; n2_temp <- n2
    ## generate two sds and one mean
    x <- stats::qchisq(sob[,1], n1_temp - 1)
    y <- stats::qchisq(sob[,2], n2_temp - 1)
    z <- stats::qnorm(sob[,3], diff, sqrt(sigma1^2/n1_temp + sigma2^2/n2_temp))

    sdv <- sqrt(sigma1^2*x/((n1_temp - 1)*n1_temp) + sigma2^2*y/(n2_temp*(n2_temp - 1)))

    dfw <- sdv^4/(sigma1^4*x^2/((n1_temp - 1)^3*n1_temp^2) + sigma2^4*y^2/((n2_temp - 1)^3*n2_temp^2))

    if (deltaU == Inf){
      thresUp <- rep(Inf, length(z))
    }
    else{
      thresUp <- (deltaU - z)/stats::qt(1-alpha, dfw)
    }

    if (deltaL == -Inf){
      thresLow <- rep(Inf, length(z))
    }
    else{
      thresLow <- (z - deltaL)/stats::qt(1-alpha, dfw)
    }

    thres <- pmin(thresUp, thresLow)

    pwrs <- thres - sdv

    consistency <- ifelse(samps < n_rough, round(pwrs,2) >= 0, round(pwrs,2) <= 0)
    consistency <- ifelse(consistency == 1, 1, as.numeric(abs(pwrs) < 0.02))

    ## for any points where the root-finding algorithm has caused issues,
    ## re-run the root-finding algorithm starting at n_*
    inconsistent <- which(consistency == 0)
    if (length(inconsistent) > 0){

      targetPowerfnCross <- function(u, deltaL, deltaU, diff, sigma1, sigma2, n_val, alpha, q){
        n1 <- n_val[1]; n2 <- max(2, q*n_val)
        x <- stats::qchisq(u[1], n1 - 1)
        y <- stats::qchisq(u[2], n2 - 1)
        z <- stats::qnorm(u[3], diff, sqrt(sigma1^2/n1 + sigma2^2/n2))

        sdv <- sqrt(sigma1^2*x/((n1 - 1)*n1) + sigma2^2*y/(n2*(n2 - 1)))

        dfw <- sdv^4/(sigma1^4*x^2/((n1 - 1)^3*n1^2) + sigma2^4*y^2/((n2 - 1)^3*n2^2))

        if (deltaU == Inf){
          thres <- (z - deltaL)/stats::qt(1-alpha, dfw)
        }
        else if (deltaL == -Inf){
          thres <- (deltaU - z)/stats::qt(1-alpha, dfw)
        }
        else{
          thresUp <- (deltaU - z)/stats::qt(1-alpha, dfw)
          thresLow <- (z - deltaL)/stats::qt(1-alpha, dfw)
          thres <- pmin(thresUp, thresLow)
        }

        return(thres - sdv)
      }

      for (i in 1:length(inconsistent)){
        if (pwrs[inconsistent[i]] < 0){
          if (n_rough > 0.5*upper.bound){
            upper.bound <- 2*n_rough
          }
          samps[inconsistent[i]] <- uu(targetPowerfnCross, lower = n_rough, upper = upper.bound,
                                       f_lower = pwrs[inconsistent[i]],
                                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                                       u = sob[inconsistent[i],], alpha = alpha, q = q)
        }
        else{
          samps[inconsistent[i]] <- uu(targetPowerfnCross, lower = 2, upper = n_rough,
                                       f_upper = pwrs[inconsistent[i]],
                                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                                       u = sob[inconsistent[i],], alpha = alpha, q = q)
        }
      }
    }
    funecdf <- stats::ecdf(samps)

    ecdf_root <- function(quant, pwr){return(funecdf(quant) - pwr)}
    n_rough <- uu(ecdf_root, lower = stats::quantile(samps, targetPower*0.5), upper = stats::quantile(samps, targetPower + 0.5*(1 - targetPower)),
                  pwr = targetPower)

    ## get confirmatory power estimate for output
    n1_temp <- ceiling(n_rough); n2_temp <- pmax(2, ceiling(q*n_rough))
    ## generate two sds and one mean
    x <- stats::qchisq(sob[,1], n1_temp - 1)
    y <- stats::qchisq(sob[,2], n2_temp - 1)
    z <- stats::qnorm(sob[,3], diff, sqrt(sigma1^2/n1_temp + sigma2^2/n2_temp))

    sdv <- sqrt(sigma1^2*x/((n1_temp - 1)*n1_temp) + sigma2^2*y/(n2_temp*(n2_temp - 1)))

    dfw <- sdv^4/(sigma1^4*x^2/((n1_temp - 1)^3*n1_temp^2) + sigma2^4*y^2/((n2_temp - 1)^3*n2_temp^2))

    if (deltaU == Inf){
      thresUp <- rep(Inf, length(z))
    }
    else{
      thresUp <- (deltaU - z)/stats::qt(1-alpha, dfw)
    }

    if (deltaL == -Inf){
      thresLow <- rep(Inf, length(z))
    }
    else{
      thresLow <- (z - deltaL)/stats::qt(1-alpha, dfw)
    }

    thres <- pmin(thresUp, thresLow)

    pwrs <- thres - sdv

    df_samps <- data.frame(n_plot = samps)

    n_plot <- NULL
    if(plot == TRUE){
      plot_pwr <- ggplot2::ggplot(df_samps, ggplot2::aes(x = n_plot)) + ggplot2::theme_bw() +
        ggplot2::stat_ecdf(geom = "step", pad = FALSE, colour = cbbPalette[6], linewidth = 2) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 13)) +
        ggplot2::theme(axis.text.x =  ggplot2::element_text(size = 13)) +
        ggplot2::labs(title = "Approximated Power Curve") +
        ggplot2::labs(y = "Power", x=bquote(italic(n)[1]*"  ("*italic(n)[2]*" = "*.(round(q,3))*italic(n)[1]*")")) +
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
                     b = "Noninferiority test power calculation (for Group 2)",
                     c = "Noninferiority test power calculation (for Group 1)")

    NOTE <- paste0("n1 and n2 are the # of subjects in sequences 1 and 2, and sigma1 and sigma2 denote intra-subject sd.")

    results <- structure(list(n1 = n1_temp, n2 = n2_temp, q = q, diff = diff, sigma1 = 2*sigma1, sigma2 = 2*sigma2,
                              sig.level = alpha, power = mean(pwrs >= 0), bounds = c(deltaL, deltaU),
                              note = NOTE,
                              method = METHOD, upper.bound = upper.bound,
                              samps = samps, seed = seed, sobol = sobol, design = "Crossover2x2Unequal"), class = "power.en.test")
    if (sigma1 == sigma2){
      message("This function works with equal population variances, but did you mean to use DesignCrossover2x2Equal()?")
    }
    return(results)
  }

  if (design == "Crossover2x2Equal"){

    diff <- power.en.test$diff
    sigma <- power.en.test$sigma/2
    q <- power.en.test$q

    if(!is.numeric(diff) | length(diff) != 1) {
      stop("Please specify a valid number for diff.")}

    if(!is.numeric(sigma) | length(sigma) != 1){
      stop("Please specify a valid number for sigma.")}
    else if (sigma <= 0 | !is.finite(sigma)){
      stop("Please specify a valid number for sigma.")}

    if(sum(c((length(targetPower) != 1 | length(q) != 1)))){
      stop("Please specify a valid input for targetPower.")}
    if(length(targetPower) == 1 & length(q) == 1){
      if(!is.numeric(targetPower)) {
        stop("Please specify a valid number for targetPower.")}
      if (is.numeric(targetPower)){
        if (targetPower <= 0 | targetPower >= 1){
          stop("Please specify a valid number for targetPower.")}
        if(!is.numeric(q)) {
          stop("Please specify a valid number for q.")}
        else if (is.numeric(q)){
          if (q <= 0) {
            stop("Please specify a valid number for q.")}
        }
        if (diff >= deltaU | diff <= deltaL){
          stop("Please ensure diff is between deltaL and deltaU.")
        }
      }
    }

    if (is.null(seed)){
      seed <- ceiling(1000*stats::runif(1))
    }
    sob <- qrng::sobol(2^(sobol + 10), d = 2, randomize = "digital.shift", seed = seed)

    samps <- power.en.test$samps
    funecdf <- stats::ecdf(samps)

    ecdf_root <- function(quant, pwr){return(funecdf(quant) - pwr)}
    n_rough <- uu(ecdf_root, lower = stats::quantile(samps, targetPower*0.5), upper = stats::quantile(samps, targetPower + 0.5*(1 - targetPower)),
                  pwr = targetPower)

    pwrs <- NULL
    n1 <- n_rough
    n2 <- pmax(2, q*n1)

    ## explore sampling distribution at the preliminary sample size recommendation
    n1_temp <- n1; n2_temp <- n2

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

    pwrs <- thres - sdv

    consistency <- ifelse(samps < n_rough, round(pwrs,2) >= 0, round(pwrs,2) <= 0)
    consistency <- ifelse(consistency == 1, 1, as.numeric(abs(pwrs) < 0.02))

    ## for any points where the root-finding algorithm has caused issues,
    ## re-run the root-finding algorithm starting at n_*
    inconsistent <- which(consistency == 0)
    if (length(inconsistent) > 0){

      targetPowerfnCrossEq <- function(u, deltaL, deltaU, diff, sigma, n_val, alpha, q){
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

      for (i in 1:length(inconsistent)){
        if (pwrs[inconsistent[i]] < 0){
          if (n_rough > 0.5*upper.bound){
            upper.bound <- 2*n_rough
          }
          samps[inconsistent[i]] <- uu(targetPowerfnCrossEq, lower = n_rough, upper = upper.bound,
                                       f_lower = pwrs[inconsistent[i]],
                                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma = sigma,
                                       u = sob[inconsistent[i],], alpha = alpha, q = q)
        }
        else{
          samps[inconsistent[i]] <- uu(targetPowerfnCrossEq, lower = 2, upper = n_rough,
                                       f_upper = pwrs[inconsistent[i]],
                                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma = sigma,
                                       u = sob[inconsistent[i],], alpha = alpha, q = q)
        }
      }
    }
    funecdf <- stats::ecdf(samps)

    ecdf_root <- function(quant, pwr){return(funecdf(quant) - pwr)}
    n_rough <- uu(ecdf_root, lower = stats::quantile(samps, targetPower*0.5), upper = stats::quantile(samps, targetPower + 0.5*(1 - targetPower)),
                  pwr = targetPower)

    ## get confirmatory power estimate for output
    n1_temp <- ceiling(n_rough); n2_temp <- pmax(2, ceiling(q*n_rough))

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

    pwrs <- thres - sdv

    df_samps <- data.frame(n_plot = samps)

    n_plot <- NULL
    if (plot == TRUE){
      plot_pwr <- ggplot2::ggplot(df_samps, ggplot2::aes(x = n_plot)) + ggplot2::theme_bw() +
        ggplot2::stat_ecdf(geom = "step", pad = FALSE, colour = cbbPalette[6], linewidth = 2) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 13)) +
        ggplot2::theme(axis.text.x =  ggplot2::element_text(size = 13)) +
        ggplot2::labs(title = "Approximated Power Curve") +
        ggplot2::labs(y = "Power", x=bquote(italic(n)[1]*"  ("*italic(n)[2]*" = "*.(round(q,3))*italic(n)[1]*")")) +
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
                     b = "Noninferiority test power calculation (for Group 2)",
                     c = "Noninferiority test power calculation (for Group 1)")

    NOTE <- paste0("n1 and n2 are the # of subjects in sequences 1 and 2, and sigma denotes intra-subject sd.")

    results <- structure(list(n1 = n1_temp, n2 = n2_temp, q = q, diff = diff, sigma = 2*sigma,
                              sig.level = alpha, power = mean(pwrs >= 0), bounds = c(deltaL, deltaU),
                              note = NOTE,
                              method = METHOD, upper.bound = upper.bound,
                              samps = samps, seed = seed, sobol = sobol, design = "Crossover2x2Equal"), class = "power.en.test")
    return(results)
  }

  if (design == "CrossoverDualUnequal"){

    diff <- power.en.test$diff
    sigma1 <- power.en.test$sigma1/sqrt(8)
    sigma2 <- power.en.test$sigma2/sqrt(8)
    q <- power.en.test$q

    if(!is.numeric(diff) | length(diff) != 1) {
      stop("Please specify a valid number for diff.")}
    if(!is.numeric(sigma1) | length(sigma1) != 1){
      stop("Please specify a valid number for sigma1.")}
    else if (sigma1 <= 0 | !is.finite(sigma1)){
      stop("Please specify a valid number for sigma1.")}
    if(!is.numeric(sigma2) | length(sigma2) != 1){
      stop("Please specify a valid number for sigma2.")}
    else if (sigma2 <= 0 | !is.finite(sigma2)){
      stop("Please specify a valid number for sigma2.")}

    if(sum(c((length(targetPower) != 1 | length(q) != 1)))){
      stop("Please specify a valid input for targetPower.")}
    if(length(targetPower) == 1 & length(q) == 1){
      if(!is.numeric(targetPower)) {
        stop("Please specify a valid number for targetPower.")}
      # if(is.null(targetPower) | !is.numeric(targetPower)) {
      #   stop("Please specify a valid number for targetPower.")}
      if (is.numeric(targetPower)){
        if (targetPower <= 0 | targetPower >= 1){
          stop("Please specify a valid number for targetPower.")}
        if(!is.numeric(q)) {
          stop("Please specify a valid number for q.")}
        else if (is.numeric(q)){
          if (q <= 0) {
            stop("Please specify a valid number for q.")}
        }
        if (diff >= deltaU | diff <= deltaL){
          stop("Please ensure diff is between deltaL and deltaU.")
        }
      }
    }

    if (is.null(seed)){
      seed <- ceiling(1000*stats::runif(1))
    }
    sob <- qrng::sobol(2^(sobol + 10), d = 3, randomize = "digital.shift", seed = seed)

    samps <- power.en.test$samps
    funecdf <- stats::ecdf(samps)

    ecdf_root <- function(quant, pwr){return(funecdf(quant) - pwr)}
    n_rough <- uu(ecdf_root, lower = stats::quantile(samps, targetPower*0.5), upper = stats::quantile(samps, targetPower + 0.5*(1 - targetPower)),
                  pwr = targetPower)

    pwrs <- NULL
    n1 <- n_rough
    n2 <- pmax(2, q*n1)

    ## explore sampling distribution at the preliminary sample size recommendation
    n1_temp <- n1; n2_temp <- n2
    ## generate two sds and one mean
    x <- stats::qchisq(sob[,1], n1_temp - 1)
    y <- stats::qchisq(sob[,2], n2_temp - 1)
    z <- stats::qnorm(sob[,3], diff, sqrt(sigma1^2/n1_temp + sigma2^2/n2_temp))

    sdv <- sqrt(sigma1^2*x/((n1_temp - 1)*n1_temp) + sigma2^2*y/(n2_temp*(n2_temp - 1)))

    dfw <- sdv^4/(sigma1^4*x^2/((n1_temp - 1)^3*n1_temp^2) + sigma2^4*y^2/((n2_temp - 1)^3*n2_temp^2))

    if (deltaU == Inf){
      thresUp <- rep(Inf, length(z))
    }
    else{
      thresUp <- (deltaU - z)/stats::qt(1-alpha, dfw)
    }

    if (deltaL == -Inf){
      thresLow <- rep(Inf, length(z))
    }
    else{
      thresLow <- (z - deltaL)/stats::qt(1-alpha, dfw)
    }

    thres <- pmin(thresUp, thresLow)

    pwrs <- thres - sdv

    consistency <- ifelse(samps < n_rough, round(pwrs,2) >= 0, round(pwrs,2) <= 0)
    consistency <- ifelse(consistency == 1, 1, as.numeric(abs(pwrs) < 0.02))

    ## for any points where the root-finding algorithm has caused issues,
    ## re-run the root-finding algorithm starting at n_*
    inconsistent <- which(consistency == 0)
    if (length(inconsistent) > 0){

      targetPowerfnDual <- function(u, deltaL, deltaU, diff, sigma1, sigma2, n_val, alpha, q){
        n1 <- n_val[1]; n2 <- max(2, q*n_val)
        x <- stats::qchisq(u[1], n1 - 1)
        y <- stats::qchisq(u[2], n2 - 1)
        z <- stats::qnorm(u[3], diff, sqrt(sigma1^2/n1 + sigma2^2/n2))

        sdv <- sqrt(sigma1^2*x/((n1 - 1)*n1) + sigma2^2*y/(n2*(n2 - 1)))

        dfw <- sdv^4/(sigma1^4*x^2/((n1 - 1)^3*n1^2) + sigma2^4*y^2/((n2 - 1)^3*n2^2))

        if (deltaU == Inf){
          thres <- (z - deltaL)/stats::qt(1-alpha, dfw)
        }
        else if (deltaL == -Inf){
          thres <- (deltaU - z)/stats::qt(1-alpha, dfw)
        }
        else{
          thresUp <- (deltaU - z)/stats::qt(1-alpha, dfw)
          thresLow <- (z - deltaL)/stats::qt(1-alpha, dfw)
          thres <- pmin(thresUp, thresLow)
        }

        return(thres - sdv)
      }

      for (i in 1:length(inconsistent)){
        if (pwrs[inconsistent[i]] < 0){
          if (n_rough > 0.5*upper.bound){
            upper.bound <- 2*n_rough
          }
          samps[inconsistent[i]] <- uu(targetPowerfnDual, lower = n_rough, upper = upper.bound,
                                       f_lower = pwrs[inconsistent[i]],
                                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                                       u = sob[inconsistent[i],], alpha = alpha, q = q)
        }
        else{
          samps[inconsistent[i]] <- uu(targetPowerfnDual, lower = 2, upper = n_rough,
                                       f_upper = pwrs[inconsistent[i]],
                                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                                       u = sob[inconsistent[i],], alpha = alpha, q = q)
        }
      }
    }
    funecdf <- stats::ecdf(samps)

    ecdf_root <- function(quant, pwr){return(funecdf(quant) - pwr)}
    n_rough <- uu(ecdf_root, lower = stats::quantile(samps, targetPower*0.5), upper = stats::quantile(samps, targetPower + 0.5*(1 - targetPower)),
                  pwr = targetPower)

    ## get confirmatory power estimate for output
    n1_temp <- ceiling(n_rough); n2_temp <- pmax(2, ceiling(q*n_rough))
    ## generate two sds and one mean
    x <- stats::qchisq(sob[,1], n1_temp - 1)
    y <- stats::qchisq(sob[,2], n2_temp - 1)
    z <- stats::qnorm(sob[,3], diff, sqrt(sigma1^2/n1_temp + sigma2^2/n2_temp))

    sdv <- sqrt(sigma1^2*x/((n1_temp - 1)*n1_temp) + sigma2^2*y/(n2_temp*(n2_temp - 1)))

    dfw <- sdv^4/(sigma1^4*x^2/((n1_temp - 1)^3*n1_temp^2) + sigma2^4*y^2/((n2_temp - 1)^3*n2_temp^2))

    if (deltaU == Inf){
      thresUp <- rep(Inf, length(z))
    }
    else{
      thresUp <- (deltaU - z)/stats::qt(1-alpha, dfw)
    }

    if (deltaL == -Inf){
      thresLow <- rep(Inf, length(z))
    }
    else{
      thresLow <- (z - deltaL)/stats::qt(1-alpha, dfw)
    }

    thres <- pmin(thresUp, thresLow)

    pwrs <- thres - sdv

    df_samps <- data.frame(n_plot = samps)

    n_plot <- NULL
    if(plot == TRUE){
      plot_pwr <- ggplot2::ggplot(df_samps, ggplot2::aes(x = n_plot)) + ggplot2::theme_bw() +
        ggplot2::stat_ecdf(geom = "step", pad = FALSE, colour = cbbPalette[6], linewidth = 2) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 13)) +
        ggplot2::theme(axis.text.x =  ggplot2::element_text(size = 13)) +
        ggplot2::labs(title = "Approximated Power Curve") +
        ggplot2::labs(y = "Power", x=bquote(italic(n)[1]*"  ("*italic(n)[2]*" = "*.(round(q,3))*italic(n)[1]*")")) +
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
                     b = "Noninferiority test power calculation (for Group 2)",
                     c = "Noninferiority test power calculation (for Group 1)")

    NOTE <- paste0("n1 and n2 are the # of subjects in sequences 1 and 2, and sigma1 and sigma2 denote intra-subject sd.")

    results <- structure(list(n1 = n1_temp, n2 = n2_temp, q = q, diff = diff, sigma1 = sqrt(8)*sigma1, sigma2 = sqrt(8)*sigma2,
                              sig.level = alpha, power = mean(pwrs >= 0), bounds = c(deltaL, deltaU),
                              note = NOTE, upper.bound = upper.bound,
                              method = METHOD,
                              samps = samps, seed = seed, sobol = sobol, design = "CrossoverDualUnequal"), class = "power.en.test")
    if (sigma1 == sigma2){
      message("This function works with equal population variances, but did you mean to use DesignCrossoverDualEqual()?")
    }
    return(results)
  }

  if (design == "CrossoverDualEqualNoComp"){

    diff <- power.en.test$diff
    sigma <- power.en.test$sigma/sqrt(8)
    q <- power.en.test$q

    if(!is.numeric(diff) | length(diff) != 1) {
      stop("Please specify a valid number for diff.")}

    if(!is.numeric(sigma) | length(sigma) != 1){
      stop("Please specify a valid number for sigma.")}
    else if (sigma <= 0 | !is.finite(sigma)){
      stop("Please specify a valid number for sigma.")}

    if(sum(c((length(targetPower) != 1 | length(q) != 1)))){
      stop("Please specify a valid input for targetPower.")}
    if(length(targetPower) == 1 & length(q) == 1){
      if(!is.numeric(targetPower)) {
        stop("Please specify a valid number for targetPower.")}
      if (is.numeric(targetPower)){
        if (targetPower <= 0 | targetPower >= 1){
          stop("Please specify a valid number for targetPower.")}
        if(!is.numeric(q)) {
          stop("Please specify a valid number for q.")}
        else if (is.numeric(q)){
          if (q <= 0) {
            stop("Please specify a valid number for q.")}
        }
        if (diff >= deltaU | diff <= deltaL){
          stop("Please ensure diff is between deltaL and deltaU.")
        }
      }
    }

    if (is.null(seed)){
      seed <- ceiling(1000*stats::runif(1))
    }
    sob <- qrng::sobol(2^(sobol + 10), d = 2, randomize = "digital.shift", seed = seed)

    samps <- power.en.test$samps
    funecdf <- stats::ecdf(samps)

    ecdf_root <- function(quant, pwr){return(funecdf(quant) - pwr)}
    n_rough <- uu(ecdf_root, lower = stats::quantile(samps, targetPower*0.5), upper = stats::quantile(samps, targetPower + 0.5*(1 - targetPower)),
                  pwr = targetPower)

    pwrs <- NULL
    n1 <- n_rough
    n2 <- pmax(2, q*n1)

    ## explore sampling distribution at the preliminary sample size recommendation
    n1_temp <- n1; n2_temp <- n2
    df_temp <- n1_temp + n2_temp - 2

    x <- stats::qchisq(sob[,1], df_temp)
    z <- stats::qnorm(sob[,2], diff, sigma*sqrt(1/n1_temp + 1/n2_temp))

    sdv <- sqrt(x*sigma^2/(df_temp))*sqrt(1/n1_temp + 1/n2_temp)

    if (deltaU == Inf){
      thresUp <- rep(Inf, length(z))
    }
    else{
      thresUp <- (deltaU - z)/stats::qt(1-alpha, df_temp)
    }

    if (deltaL == -Inf){
      thresLow <- rep(Inf, length(z))
    }
    else{
      thresLow <- (z - deltaL)/stats::qt(1-alpha, df_temp)
    }

    thres <- pmin(thresUp, thresLow)

    pwrs <- thres - sdv

    consistency <- ifelse(samps < n_rough, round(pwrs,2) >= 0, round(pwrs,2) <= 0)
    consistency <- ifelse(consistency == 1, 1, as.numeric(abs(pwrs) < 0.02))

    ## for any points where the root-finding algorithm has caused issues,
    ## re-run the root-finding algorithm starting at n_*
    inconsistent <- which(consistency == 0)
    if (length(inconsistent) > 0){

      targetPowerfnNoComp <- function(u, deltaL, deltaU, diff, sigma, n_val, alpha, q){
        n1 <- n_val[1]; n2 <- max(2, q*n_val)
        df_use <- n1 + n2 - 2

        x <- stats::qchisq(u[1], df_use)
        z <- stats::qnorm(u[2], diff, sigma*sqrt(1/n1 + 1/n2))

        sdv <- sqrt(x*sigma^2/(df_use))*sqrt(1/n1 + 1/n2)

        if (deltaU == Inf){
          thres <- (z - deltaL)/stats::qt(1-alpha, df_use)
        }
        else if (deltaL == -Inf){
          thres <- (deltaU - z)/stats::qt(1-alpha, df_use)
        }
        else{
          thresUp <- (deltaU - z)/stats::qt(1-alpha, df_use)
          thresLow <- (z - deltaL)/stats::qt(1-alpha, df_use)
          thres <- pmin(thresUp, thresLow)
        }

        return(thres - sdv)
      }

      for (i in 1:length(inconsistent)){
        if (pwrs[inconsistent[i]] < 0){
          if (n_rough > 0.5*upper.bound){
            upper.bound <- 2*n_rough
          }
          samps[inconsistent[i]] <- uu(targetPowerfnNoComp, lower = n_rough, upper = upper.bound,
                                       f_lower = pwrs[inconsistent[i]],
                                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma = sigma,
                                       u = sob[inconsistent[i],], alpha = alpha, q = q)
        }
        else{
          samps[inconsistent[i]] <- uu(targetPowerfnNoComp, lower = 2, upper = n_rough,
                                       f_upper = pwrs[inconsistent[i]],
                                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma = sigma,
                                       u = sob[inconsistent[i],], alpha = alpha, q = q)
        }
      }
    }
    funecdf <- stats::ecdf(samps)

    ecdf_root <- function(quant, pwr){return(funecdf(quant) - pwr)}
    n_rough <- uu(ecdf_root, lower = stats::quantile(samps, targetPower*0.5), upper = stats::quantile(samps, targetPower + 0.5*(1 - targetPower)),
                  pwr = targetPower)

    ## get confirmatory power estimate for output
    n1_temp <- ceiling(n_rough); n2_temp <- pmax(2, ceiling(q*n_rough))
    df_temp <- n1_temp + n2_temp - 2

    x <- stats::qchisq(sob[,1], df_temp)
    z <- stats::qnorm(sob[,2], diff, sigma*sqrt(1/n1_temp + 1/n2_temp))

    sdv <- sqrt(x*sigma^2/(df_temp))*sqrt(1/n1_temp + 1/n2_temp)

    if (deltaU == Inf){
      thresUp <- rep(Inf, length(z))
    }
    else{
      thresUp <- (deltaU - z)/stats::qt(1-alpha, df_temp)
    }

    if (deltaL == -Inf){
      thresLow <- rep(Inf, length(z))
    }
    else{
      thresLow <- (z - deltaL)/stats::qt(1-alpha, df_temp)
    }

    thres <- pmin(thresUp, thresLow)

    pwrs <- thres - sdv

    df_samps <- data.frame(n_plot = samps)

    n_plot <- NULL
    if (plot == TRUE){
      plot_pwr <- ggplot2::ggplot(df_samps, ggplot2::aes(x = n_plot)) + ggplot2::theme_bw() +
        ggplot2::stat_ecdf(geom = "step", pad = FALSE, colour = cbbPalette[6], linewidth = 2) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 13)) +
        ggplot2::theme(axis.text.x =  ggplot2::element_text(size = 13)) +
        ggplot2::labs(title = "Approximated Power Curve") +
        ggplot2::labs(y = "Power", x=bquote(italic(n)[1]*"  ("*italic(n)[2]*" = "*.(round(q,3))*italic(n)[1]*")")) +
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
                     b = "Noninferiority test power calculation (for Group 2)",
                     c = "Noninferiority test power calculation (for Group 1)")

    NOTE <- paste0("n1 and n2 are the # of subjects in sequences 1 and 2, and sigma denotes intra-subject sd.")

    results <- structure(list(n1 = n1_temp, n2 = n2_temp, q = q, diff = diff, sigma = sqrt(8)*sigma,
                              sig.level = alpha, power = mean(pwrs>= 0), bounds = c(deltaL, deltaU),
                              note = NOTE,
                              method = METHOD, upper.bound = upper.bound,
                              samps = samps, seed = seed, sobol = sobol, design = "CrossoverDualEqualNoComp"), class = "power.en.test")
    return(results)
  }

  if (design == "CrossoverDualEqualComp"){

    diff <- power.en.test$diff
    sigma <- power.en.test$sigma/sqrt(8)
    q <- power.en.test$q

    if(!is.numeric(diff) | length(diff) != 1) {
      stop("Please specify a valid number for diff.")}

    if(!is.numeric(sigma) | length(sigma) != 1){
      stop("Please specify a valid number for sigma.")}
    else if (sigma <= 0 | !is.finite(sigma)){
      stop("Please specify a valid number for sigma.")}

    if(sum(c((length(targetPower) != 1 | length(q) != 1)))){
      stop("Please specify a valid input for targetPower.")}
    if(length(targetPower) == 1 & length(q) == 1){
      if(!is.numeric(targetPower)) {
        stop("Please specify a valid number for targetPower.")}
      if (is.numeric(targetPower)){
        if (targetPower <= 0 | targetPower >= 1){
          stop("Please specify a valid number for targetPower.")}
        if(!is.numeric(q)) {
          stop("Please specify a valid number for q.")}
        else if (is.numeric(q)){
          if (q <= 0) {
            stop("Please specify a valid number for q.")}
        }
        if (diff >= deltaU | diff <= deltaL){
          stop("Please ensure diff is between deltaL and deltaU.")
        }
      }
    }

    if (is.null(seed)){
      seed <- ceiling(1000*stats::runif(1))
    }
    sob <- qrng::sobol(2^(sobol + 10), d = 2, randomize = "digital.shift", seed = seed)

    samps <- power.en.test$samps
    funecdf <- stats::ecdf(samps)

    ecdf_root <- function(quant, pwr){return(funecdf(quant) - pwr)}
    n_rough <- uu(ecdf_root, lower = stats::quantile(samps, targetPower*0.5), upper = stats::quantile(samps, targetPower + 0.5*(1 - targetPower)),
                  pwr = targetPower)

    pwrs <- NULL
    n1 <- n_rough
    n2 <- pmax(2, q*n1)

    ## explore sampling distribution at the preliminary sample size recommendation
    n1_temp <- n1; n2_temp <- n2
    df_temp <- 2*(n1_temp + n2_temp - 2)

    x <- stats::qchisq(sob[,1], df_temp)
    z <- stats::qnorm(sob[,2], diff, sigma*sqrt(1/n1_temp + 1/n2_temp))

    sdv <- sqrt(x*sigma^2/(df_temp))*sqrt(1/n1_temp + 1/n2_temp)

    if (deltaU == Inf){
      thresUp <- rep(Inf, length(z))
    }
    else{
      thresUp <- (deltaU - z)/stats::qt(1-alpha, df_temp)
    }

    if (deltaL == -Inf){
      thresLow <- rep(Inf, length(z))
    }
    else{
      thresLow <- (z - deltaL)/stats::qt(1-alpha, df_temp)
    }

    thres <- pmin(thresUp, thresLow)

    pwrs <- thres - sdv

    consistency <- ifelse(samps < n_rough, round(pwrs,2) >= 0, round(pwrs,2) <= 0)
    consistency <- ifelse(consistency == 1, 1, as.numeric(abs(pwrs) < 0.02))

    ## for any points where the root-finding algorithm has caused issues,
    ## re-run the root-finding algorithm starting at n_*
    inconsistent <- which(consistency == 0)
    if (length(inconsistent) > 0){

      targetPowerfnComp <- function(u, deltaL, deltaU, diff, sigma, n_val, alpha, q){
        n1 <- n_val[1]; n2 <- max(2, q*n_val)
        df_use <- 2*(n1 + n2 - 2)

        x <- stats::qchisq(u[1], df_use)
        z <- stats::qnorm(u[2], diff, sigma*sqrt(1/n1 + 1/n2))

        sdv <- sqrt(x*sigma^2/(df_use))*sqrt(1/n1 + 1/n2)

        if (deltaU == Inf){
          thres <- (z - deltaL)/stats::qt(1-alpha, df_use)
        }
        else if (deltaL == -Inf){
          thres <- (deltaU - z)/stats::qt(1-alpha, df_use)
        }
        else{
          thresUp <- (deltaU - z)/stats::qt(1-alpha, df_use)
          thresLow <- (z - deltaL)/stats::qt(1-alpha, df_use)
          thres <- pmin(thresUp, thresLow)
        }

        return(thres - sdv)
      }

      for (i in 1:length(inconsistent)){
        if (pwrs[inconsistent[i]] < 0){
          if (n_rough > 0.5*upper.bound){
            upper.bound <- 2*n_rough
          }
          samps[inconsistent[i]] <- uu(targetPowerfnComp, lower = n_rough, upper = upper.bound,
                                       f_lower = pwrs[inconsistent[i]],
                                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma = sigma,
                                       u = sob[inconsistent[i],], alpha = alpha, q = q)
        }
        else{
          samps[inconsistent[i]] <- uu(targetPowerfnComp, lower = 2, upper = n_rough,
                                       f_upper = pwrs[inconsistent[i]],
                                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma = sigma,
                                       u = sob[inconsistent[i],], alpha = alpha, q = q)
        }
      }
    }
    funecdf <- stats::ecdf(samps)

    ecdf_root <- function(quant, pwr){return(funecdf(quant) - pwr)}
    n_rough <- uu(ecdf_root, lower = stats::quantile(samps, targetPower*0.5), upper = stats::quantile(samps, targetPower + 0.5*(1 - targetPower)),
                  pwr = targetPower)

    ## get confirmatory power estimate for output
    n1_temp <- ceiling(n_rough); n2_temp <- pmax(2, ceiling(q*n_rough))
    df_temp <- 2*(n1_temp + n2_temp - 2)

    x <- stats::qchisq(sob[,1], df_temp)
    z <- stats::qnorm(sob[,2], diff, sigma*sqrt(1/n1_temp + 1/n2_temp))

    sdv <- sqrt(x*sigma^2/(df_temp))*sqrt(1/n1_temp + 1/n2_temp)

    if (deltaU == Inf){
      thresUp <- rep(Inf, length(z))
    }
    else{
      thresUp <- (deltaU - z)/stats::qt(1-alpha, df_temp)
    }

    if (deltaL == -Inf){
      thresLow <- rep(Inf, length(z))
    }
    else{
      thresLow <- (z - deltaL)/stats::qt(1-alpha, df_temp)
    }

    thres <- pmin(thresUp, thresLow)

    pwrs <- thres - sdv

    df_samps <- data.frame(n_plot = samps)

    n_plot <- NULL
    if (plot == TRUE){
      plot_pwr <- ggplot2::ggplot(df_samps, ggplot2::aes(x = n_plot)) + ggplot2::theme_bw() +
        ggplot2::stat_ecdf(geom = "step", pad = FALSE, colour = cbbPalette[6], linewidth = 2) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 13)) +
        ggplot2::theme(axis.text.x =  ggplot2::element_text(size = 13)) +
        ggplot2::labs(title = "Approximated Power Curve") +
        ggplot2::labs(y = "Power", x=bquote(italic(n)[1]*"  ("*italic(n)[2]*" = "*.(round(q,3))*italic(n)[1]*")")) +
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
                     b = "Noninferiority test power calculation (for Group 2)",
                     c = "Noninferiority test power calculation (for Group 1)")

    NOTE <- paste0("n1 and n2 are the # of subjects in sequences 1 and 2, and sigma denotes intra-subject sd.")

    results <- structure(list(n1 = n1_temp, n2 = n2_temp, q = q, diff = diff, sigma = sqrt(8)*sigma,
                              sig.level = alpha, power = mean(pwrs >= 0), bounds = c(deltaL, deltaU),
                              note = NOTE,
                              method = METHOD, upper.bound = upper.bound,
                              samps = samps, seed = seed, sobol = sobol, design = "CrossoverDualEqualComp"), class = "power.en.test")
    return(results)
  }
}
