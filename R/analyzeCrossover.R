#' Equivalence-Based Tests for Crossover Studies
#'
#' Conducts equivalence and noninferiority tests with two samples (including those with paired data).
#'
#' @param x a vector of data values with at least two distinct observations from Group 1. The observations should be the intra-subject *differences* in the first sequence/group (i.e., of the form \eqn{x_{i} = x_{i1} - x_{i2}}, where \eqn{i} denotes the \eqn{i}th subject in the 1st sequence and the second subscript denotes the period). Please exclude missing data and `NA` values.
#' @param y a vector of data values with at least two distinct observations from Group 2. The observations should be the intra-subject *differences* in the second sequence/group (i.e., of the form \eqn{y_{i} = y_{i1} - y_{i2}}, where \eqn{i} denotes the \eqn{i}th subject in the 2nd sequence and the second subscript denotes the period). Please exclude missing data and `NA` values.
#' @param type a character string specifying the type of *t*-test(s) to be used, must be one of `"welch"` (default, does not assume intra-subject variance is the same in both groups) or `"student"` (assumes intra-subject variance is the same in both groups).
#' @param deltaL the lower bound for the interval of equivalence (can be `-Inf` for noninferiority test).
#' @param deltaU the upper bound for the interval of equivalence (can be `Inf` for noninferiority test).
#' @param alpha the significance level \eqn{\alpha \in (0,1)}.
#'
#' @examples
#' # specify targetPower to obtain sample sizes n
#' AnalyzeCrossover(x = rnorm(15), y = rnorm(15), type = "welch", deltaL = -1,
#' deltaU = 1, alpha = 0.05)
#'
#' @return The test results are returned as a list with supplementary information. Equivalence or noninferiority is concluded when all *p*-values in the Results table are less than `alpha`. Otherwise, we do not reject the null hypotheses of inequivalence or inferiority.
#' @export
AnalyzeCrossover <- function(x = NULL, y = NULL, type = "welch", deltaL = -Inf,
                                  deltaU = Inf, alpha = NULL){

  ## error checking
  if(!is.numeric(x) | length(x) < 2) {
    return("Error: Please specify a valid input for data x.")}
  else if (sum(is.na(x)) > 0){
    return("Error: Please specify a valid input for data x.")}
  else if (stats::var(x) == 0){
    return("Error: Please specify a valid input for data x.")}
  if(!is.numeric(y) | length(y) < 2) {
    return("Error: Please specify a valid input for data y.")}
  else if (sum(is.na(y)) > 0){
    return("Error: Please specify a valid input for data y.")}
  else if (stats::var(y) == 0){
    return("Error: Please specify a valid input for data y.")}
  if (length(x) != length(y)){
    return("Error: Please ensure x and y have the same number of observations.")}
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
  if(!(type %in% c("welch", "student"))){
    return("Error: Please specify a valid type for the t-test(s).")}

  x <- x/2; y <- y/2

  if (is.finite(deltaL) & is.finite(deltaU)){
    t1 <- stats::t.test(x, y, mu = deltaL, alternative = "greater", var.equal = ifelse(type == "student", TRUE, FALSE))
    t2 <- stats::t.test(x, y, mu = deltaU, alternative = "less", var.equal = ifelse(type == "student", TRUE, FALSE))
    tab <- data.frame(t = c(t1$statistic, t2$statistic),
                      SE = c(t1$stderr, t2$stderr),
                      df = c(t1$parameter, t2$parameter),
                      p.value= c(t1$p.value, t2$p.value))
    row.names(tab) <- c("Lower Bound", "Upper Bound")
  }
  else if (is.finite(deltaL)){
    t1 <- stats::t.test(x, y, mu = deltaL, alternative = "greater", var.equal = ifelse(type == "student", TRUE, FALSE))
    tab <- data.frame(t = t1$statistic, SE = t1$stderr, df = t1$parameter, p.value = t1$p.value)
    row.names(tab) <- "Lower Bound"
  }
  else {
    t2 <- stats::t.test(x, y, mu = deltaU, alternative = "less", var.equal = ifelse(type == "student", TRUE, FALSE))
    tab <- data.frame(t = t2$statistic, SE = t2$stderr, df = t2$parameter, p.value = t2$p.value)
    row.names(tab) <- "Upper Bound"
  }

  METHOD <-  switch(type, welch = "Welch's Two-Sample t-test",
                    student = "Student's Two-Sample t-test")

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
