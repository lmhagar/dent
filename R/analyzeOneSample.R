#' Equivalence-Based Tests for One-Group Studies
#'
#' Conducts equivalence and noninferiority tests with one sample.
#'
#' @param x a vector of data values with at least two distinct observations. Please exclude missing data and `NA` values.
#' @param deltaL the lower bound for the interval of equivalence (can be `-Inf` for nonsuperiority test).
#' @param deltaU the upper bound for the interval of equivalence (can be `Inf` for noninferiority test).
#' @param alpha the significance level \eqn{\alpha \in (0,1)}.
#'
#' @examples
#' # analyze data from one group
#' AnalyzeOneSample(x = rnorm(15), deltaL = -1, deltaU = 1, alpha = 0.05)
#'
#' @return The test results are returned as a list with supplementary information. Equivalence, noninferiority, or nonsuperiority is concluded when all *p*-values in the Results table are less than `alpha`. Otherwise, we do not reject the null hypotheses of inequivalence, inferiority, or superiority.
#' @export
AnalyzeOneSample <- function(x = NULL, deltaL = -Inf,
                                  deltaU = Inf, alpha = NULL){

  ## error checking
  if(!is.numeric(x) | length(x) < 2) {
    stop("Please specify a valid input for data x.")}
  else if (sum(is.na(x)) > 0){
    stop("Please specify a valid input for data x.")}
  else if (stats::var(x) == 0){
    stop("Please specify a valid input for data x.")}
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

  if (is.finite(deltaL) & is.finite(deltaU)){
    t1 <- stats::t.test(x, mu = deltaL, alternative = "greater")
    t2 <- stats::t.test(x, mu = deltaU, alternative = "less")
    tab <- data.frame(t = c(t1$statistic, t2$statistic),
                      SE = c(t1$stderr, t2$stderr),
                      df = c(t1$parameter, t2$parameter),
                      p.value= c(t1$p.value, t2$p.value))
    row.names(tab) <- c("Lower Bound", "Upper Bound")
  }
  else if (is.finite(deltaL)){
    t1 <- stats::t.test(x, mu = deltaL, alternative = "greater")
    tab <- data.frame(t = t1$statistic, SE = t1$stderr, df = t1$parameter, p.value = t1$p.value)
    row.names(tab) <- "Lower Bound"
  }
  else {
    t2 <- stats::t.test(x, mu = deltaU, alternative = "less")
    tab <- data.frame(t = t2$statistic, SE = t2$stderr, df = t2$parameter, p.value = t2$p.value)
    row.names(tab) <- "Upper Bound"
  }

  METHOD <- "One Sample t-test"

  type <- ifelse(is.finite(deltaL) & is.finite(deltaU), "a",
                 ifelse(!is.finite(deltaL), "b", "c"))


  test <- paste0(switch(type, a = "Equivalence test with (deltaL, deltaU) = (",
                 b = "Nonsuperiority test with (deltaL, deltaU) = (",
                 c = "Noninferiority test with (deltaL, deltaU) = ("), round(deltaL, 4), ",", round(deltaU,4), ")")

  summary <- paste0(ifelse(max(tab$p.value) <= alpha, "Conclude ", "Do not conclude "),
                           switch(type, a = "equivalence.",
                        b = "nonsuperiority.",
                        c = "noninferiority."))

  results <- structure(list(test = test, table = tab, summary = summary,
                              sig.level = round(alpha,4),
                              method = METHOD), class = "en.test")
  return(results)
}


#' @title Printing equivalence and noninferiority test results
#' @name en.test-methods
#'
#' @description Helper function to print en.test objects
#'
#' @param x an en.test object
#' @param ... further arguments to be passed through (see `print()` function)
#'
#' @rdname en.test-methods
#' @method print en.test
#' @keywords internal
#' @export
print.en.test <- function(x, ...){
  cat(x$method, "\n",
      x$test, "\n",
      "Significance Level: ", x$sig.level, "\n", "\n", "Results", "\n", sep = "")
  print(as.matrix(x$table),quote=F)
  cat("\n", "Summary", "\n", x$summary, sep = "")
}
