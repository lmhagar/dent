#' Equivalence-Based Tests for Two-Group Studies
#'
#' Conducts equivalence and noninferiority tests with two samples (including those with paired data).
#'
#' @param x a vector of data values with at least two distinct observations from Group 1. Please exclude missing data and `NA` values.
#' @param y a vector of data values with at least two distinct observations from Group 2. Please exclude missing data and `NA` values. `x` and `y` must contain the same number of observations if `type = "paired"`.
#' @param type a character string specifying the type of *t*-test(s) to be used, must be one of `"welch"` (default, for independent samples without equal variance assumption), `"student"` (independent sample with equal variance assumption), or `"paired"` (for paired data).
#' @param deltaL the lower bound for the interval of equivalence (can be `-Inf` for noninferiority test).
#' @param deltaU the upper bound for the interval of equivalence (can be `Inf` for noninferiority test).
#' @param alpha the significance level \eqn{\alpha \in (0,1)}.
#'
#' @examples
#' # analyze data from two groups
#' AnalyzeTwoSample(x = rnorm(15), y = rnorm(15), type = "welch", deltaL = -1,
#' deltaU = 1, alpha = 0.05)
#'
#' @return The test results are returned as a list with supplementary information. Equivalence or noninferiority is concluded when all *p*-values in the Results table are less than `alpha`. Otherwise, we do not reject the null hypotheses of inequivalence or inferiority.
#' @export
AnalyzeTwoSample <- function(x = NULL, y = NULL, type = "welch", deltaL = -Inf,
                                  deltaU = Inf, alpha = NULL){

  ## error checking
  if(!is.numeric(x) | length(x) < 2) {
    stop("Please specify a valid input for data x.")}
  else if (sum(is.na(x)) > 0){
    stop("Please specify a valid input for data x.")}
  else if (stats::var(x) == 0){
    stop("Please specify a valid input for data x.")}
  if(!is.numeric(y) | length(y) < 2) {
    stop("Please specify a valid input for data y.")}
  else if (sum(is.na(y)) > 0){
    stop("Please specify a valid input for data y.")}
  else if (stats::var(y) == 0){
    stop("Please specify a valid input for data y.")}
  if(length(type) != 1 | !(type %in% c("welch", "student", "paired"))){
    stop("Please specify a valid type for the t-test(s).")}
  if (length(x) != length(y) & type == "paired"){
    stop("Please ensure x and y have the same number of observations.")}
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
    t1 <- stats::t.test(x, y, mu = deltaL, alternative = "greater", var.equal = ifelse(type == "student", TRUE, FALSE),
                        paired = ifelse(type == "paired", TRUE, FALSE))
    t2 <- stats::t.test(x, y, mu = deltaU, alternative = "less", var.equal = ifelse(type == "student", TRUE, FALSE),
                        paired = ifelse(type == "paired", TRUE, FALSE))
    tab <- data.frame(t = c(t1$statistic, t2$statistic),
                      SE = c(t1$stderr, t2$stderr),
                      df = c(t1$parameter, t2$parameter),
                      p.value= c(t1$p.value, t2$p.value))
    row.names(tab) <- c("Lower Bound", "Upper Bound")
  }
  else if (is.finite(deltaL)){
    t1 <- stats::t.test(x, y, mu = deltaL, alternative = "greater", var.equal = ifelse(type == "student", TRUE, FALSE),
                        paired = ifelse(type == "paired", TRUE, FALSE))
    tab <- data.frame(t = t1$statistic, SE = t1$stderr, df = t1$parameter, p.value = t1$p.value)
    row.names(tab) <- "Lower Bound"
  }
  else {
    t2 <- stats::t.test(x, y, mu = deltaU, alternative = "less", var.equal = ifelse(type == "student", TRUE, FALSE),
                        paired = ifelse(type == "paired", TRUE, FALSE))
    tab <- data.frame(t = t2$statistic, SE = t2$stderr, df = t2$parameter, p.value = t2$p.value)
    row.names(tab) <- "Upper Bound"
  }

  METHOD <-  switch(type, welch = "Welch's Two-Sample t-test",
                    student = "Student's Two-Sample t-test",
                    paired = "Paired Two-Sample t-test")

  typehyp <- ifelse(is.finite(deltaL) & is.finite(deltaU), "a",
                 ifelse(!is.finite(deltaL), "b", "c"))


  test <- paste0(switch(typehyp, a = "Equivalence test with (deltaL, deltaU) = (",
                 b = "Noninferiority test for group 2 with (deltaL, deltaU) = (",
                 c = "Noninferiority test for group 1 with (deltaL, deltaU) = ("), round(deltaL, 4), ",", round(deltaU,4), ")")

  summary <- paste0(ifelse(max(tab$p.value) <= alpha, "Conclude ", "Do not conclude "),
                           switch(typehyp, a = "equivalence.",
                        b = "noninferiority for group 2.",
                        c = "noninferiority for group 1."))

  results <- structure(list(test = test, table = tab, summary = summary,
                              sig.level = round(alpha,4),
                              method = METHOD), class = "en.test")
  return(results)
}
