% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/designPairedSample.R
\name{DesignPairedSample}
\alias{DesignPairedSample}
\title{Power Calculations for Designs with Paired Data}
\usage{
DesignPairedSample(
  diff = NULL,
  sigma1 = NULL,
  sigma2 = NULL,
  rho = NULL,
  deltaL = -Inf,
  deltaU = Inf,
  alpha = NULL,
  targetPower = NULL,
  n = NULL,
  plot = TRUE,
  seed = NULL,
  copies = 8,
  sobol = 0
)
}
\arguments{
\item{diff}{the anticipated difference between the group means (\eqn{\mu_{1} - \mu_{2}}).}

\item{sigma1}{the anticipated within-group standard deviation for group 1 (\eqn{\sigma_{1}}).}

\item{sigma2}{the anticipated within-group standard deviation for group 2 (\eqn{\sigma_{2}}).}

\item{rho}{the anticipated correlation between the pair of observations \eqn{\rho \in [0,1]}.}

\item{deltaL}{the lower bound for the interval of equivalence (can be \code{-Inf} for noninferiority test).}

\item{deltaU}{the upper bound for the interval of equivalence (can be \code{Inf} for noninferiority test).}

\item{alpha}{the significance level \eqn{\alpha \in (0,1)}.}

\item{targetPower}{the desired statistical power of the equivalence or noninferiority test, must be a single number between 0 and 1 (exclusive). Exactly one of the following inputs must be specified: \code{targetPower} or \code{n}. Specify \code{targetPower} to find a sample size \code{n} that yields desired power.}

\item{n}{the sample size specified in number of pairs, must be a single integer such that \eqn{n \ge 2}. Exactly one of the following inputs must be specified: \code{targetPower} or \code{n}. Specify \code{n} to estimate statistical power for this sample size.}

\item{plot}{a logical variable indicating whether to return a plot of the power curve. If \code{n} is specified instead of \code{targetPower}, this variable is automatically set to \code{FALSE}. If you wish to approximate many power curves, suppressing the plots will expedite this process.}

\item{seed}{if provided, a single positive integer is used to ensure reproducibility when randomizing Sobol' sequences via \code{sobol()} in the \code{qrng} package.}

\item{copies}{an integer between 4 and 10 (inclusive) denoting how many independent copies of the Sobol' sequence are used to construct confidence intervals for power. The default value is 8.}

\item{sobol}{an integer \eqn{s} between 0 and 4 (inclusive). When approximating the power curve using \code{targetPower}, \eqn{\texttt{copies} \times 2^{s + 7}} simulation repetitions are used. When estimating power for given a samples size \code{n}, \eqn{\texttt{copies} \times 2^{s + 13}} simulation repetitions are used. The default setting is \eqn{s = 0}, which ensures that each function call should take less than two seconds. As \eqn{s} increases, the sample size calculation is less sensitive to simulation variability but takes longer to complete. However, function calls should still take less than 30 seconds when \eqn{s = 4}.}
}
\value{
The sample size or power estimate are returned as a list with supplementary information.
If \code{targetPower} is specified to find sample size \code{n}, a plot of the approximated power curve will also appear in the plot pane if \code{plot = TRUE}. To find a sample size that corresponds to a different \code{targetPower}, save this function's output to an object and use the \code{UpdateTargetPower()} function.
}
\description{
Approximates the power of equivalence and noninferiority tests with paired data. One can either find a sample size (i.e., number of pairs) that achieves desired statistical power or estimate the power for a given sample size.
}
\examples{
# specify targetPower to obtain sample size n
DesignPairedSample(diff = -4, sigma1 = 15, sigma2 = 18, rho = 0.25, deltaL = -19.2,
deltaU = 19.2, targetPower = 0.8, alpha = 0.05, plot = TRUE, seed = 1)

# specify n to estimate power for this design
DesignPairedSample(diff = -4, sigma1 = 15, sigma2 = 18, rho = 0.25, deltaL = -19.2,
deltaU = 19.2, n = 17, alpha = 0.05, plot = FALSE, seed = 1)

}
