#' Equivalence-Based Tests for 2\eqn{\times}2 Crossover Studies
#'
#' Conducts equivalence and noninferiority tests for 2\eqn{\times}2 crossover studies.
#'
#' @param X a matrix of data values with at least two rows of observations from Group 1. The observations in each row should be of the form \eqn{\boldsymbol{x}_{i} = (x_{i1}, x_{i2})}, where \eqn{i} denotes the \eqn{i}th subject in the 1st sequence and the second subscript denotes the period. That is, column \eqn{j} corresponds to the \eqn{j}th period of sequence 1. The relevant intra-subject contrasts will be computed, and there should be at least two distinct values for these contrasts. Please exclude missing data and `NA` values.
#' @param Y a matrix of data values with at least two rows of observations from Group 2. The observations in each row should be of the form \eqn{\boldsymbol{y}_{i} = (y_{i1}, y_{i2})}, where \eqn{i} denotes the \eqn{i}th subject in the 2nd sequence and the second subscript denotes the period. That is, column \eqn{j} corresponds to the \eqn{j}th period of sequence 2. The relevant intra-subject contrasts will be computed, and there should be at least two distinct values for these contrasts. Please exclude missing data and `NA` values.
#' @param type a character string specifying the type of *t*-test(s) to be used, must be one of `"welch"` (default, does not assume intra-subject variance is the same in both groups) or `"student"` (assumes intra-subject variance is the same in both groups).
#' @param deltaL the lower bound for the interval of equivalence (can be `-Inf` for noninferiority test).
#' @param deltaU the upper bound for the interval of equivalence (can be `Inf` for noninferiority test).
#' @param alpha the significance level \eqn{\alpha \in (0,1)}.
#'
#' @examples
#' # analyze data from a 2x2 crossover trial
#' AnalyzeCrossover2x2(X = cbind(rnorm(15), rnorm(15, 0.5)), Y = cbind(rnorm(10, 0.5), rnorm(10)),
#' type = "welch", deltaL = -1, deltaU = 1, alpha = 0.05)
#'
#' @return The test results are returned as a list with supplementary information. Equivalence or noninferiority is concluded when all *p*-values in the Results table are less than `alpha`. Otherwise, we do not reject the null hypotheses of inequivalence or inferiority.
#' @export
AnalyzeCrossover2x2 <- function(X = NULL, Y = NULL, type = "welch", deltaL = -Inf,
                             deltaU = Inf, alpha = NULL){

  ## error checking
  if(!is.numeric(X[,1]) | !is.numeric(X[,2]) | nrow(X) < 2) {
    stop("Please specify a valid input for data X.")}
  else if (sum(is.na(X)) > 0){
    stop("Please specify a valid input for data x.")}
  else if (stats::var(X[,1] - X[,2]) == 0){
    stop("All intra-subject contrasts are the same for data X.")}
  if(!is.numeric(Y[,1]) | !is.numeric(Y[,2]) | length(Y) < 2) {
    stop("Please specify a valid input for data Y.")}
  else if (sum(is.na(Y)) > 0){
    stop("Please specify a valid input for data Y.")}
  else if (stats::var(Y[,1] - Y[,2]) == 0){
    stop("All intra-subject contrasts are the same for data Y.")}
  # if (length(x) != length(y)){
  #   stop("Please ensure x and y have the same number of observations.")}
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
  if(!(type %in% c("welch", "student"))){
    stop("Please specify a valid type for the t-test(s).")}

  x <- 0.5*(X[,1] - X[,2]); y <- 0.5*(Y[,1] - Y[,2])

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

#' Equivalence-Based Tests for Dual Crossover Studies
#'
#' Conducts equivalence and noninferiority tests for two-sequence dual crossover designs.
#' @param X a matrix of data values with at least two rows of observations from Group 1. The observations in each row should be of the form \eqn{\boldsymbol{x}_{i} = (x_{i1}, x_{i2}, x_{i3})}, where \eqn{i} denotes the \eqn{i}th subject in the 1st sequence and the second subscript denotes the period. That is, column \eqn{j} corresponds to the \eqn{j}th period of sequence 1. The relevant intra-subject contrasts will be computed, and there should be at least two distinct values for these contrasts. Please exclude missing data and `NA` values.
#' @param Y a matrix of data values with at least two rows of observations from Group 2. The observations in each row should be of the form \eqn{\boldsymbol{y}_{i} = (y_{i1}, y_{i2}, y_{i3})}, where \eqn{i} denotes the \eqn{i}th subject in the 2nd sequence and the second subscript denotes the period. That is, column \eqn{j} corresponds to the \eqn{j}th period of sequence 2. The relevant intra-subject contrasts will be computed, and there should be at least two distinct values for these contrasts. Please exclude missing data and `NA` values.
#' @param type a character string specifying the type of *t*-test(s) to be used, must be one of `"welch"` (default, does not assume intra-subject variance is the same in both groups) or `"student"` (assumes intra-subject variance is the same in both groups).
#' @param deltaL the lower bound for the interval of equivalence (can be `-Inf` for noninferiority test).
#' @param deltaU the upper bound for the interval of equivalence (can be `Inf` for noninferiority test).
#' @param alpha the significance level \eqn{\alpha \in (0,1)}.
#' @param compSymm a logical variable indicating whether the covariance matrix for the three responses from the same subject can be assumed to take a compound symmetric structure. If so, the relevant Student's *t*-tests are conducted with \eqn{2(n_1 + n_2 - 2)} degrees of freedom. If not, the *t*-tests are conducted with \eqn{n_1 + n_2 - 2} degrees of freedom. The default setting is `FALSE`. Please see Chapter 9 of Chow and Liu (2008) for more information.
#'
#' @examples
#' # analyze data from a two-sequence dual crossover trial
#' AnalyzeCrossoverDual(X = cbind(rnorm(15), rnorm(15, 0.5), rnorm(15, 0.25)),
#' Y = cbind(rnorm(10, 0.5), rnorm(10), rnorm(10, 0.75)), type = "welch",
#' deltaL = -1, deltaU = 1, alpha = 0.05, compSymm = FALSE)
#'
#' @references
#' Chow, S. C. and J.P. Liu. (2008). *Design and analysis of bioavailability and bioequivalence studies*. Chapman and Hall/CRC.
#'
#' @return The test results are returned as a list with supplementary information. Equivalence or noninferiority is concluded when all *p*-values in the Results table are less than `alpha`. Otherwise, we do not reject the null hypotheses of inequivalence or inferiority.
#' @export
AnalyzeCrossoverDual <- function(X = NULL, Y = NULL, type = "welch", deltaL = -Inf,
                                  deltaU = Inf, alpha = NULL, compSymm = FALSE){

  ## error checking
  if(!is.numeric(X[,1]) | !is.numeric(X[,2]) | !is.numeric(X[,3]) | nrow(X) < 2) {
    stop("Please specify a valid input for data X.")}
  else if (sum(is.na(X)) > 0){
    stop("Please specify a valid input for data x.")}
  else if (stats::var(2*X[,1] - X[,2] - X[,3]) == 0){
    stop("All intra-subject contrasts are the same for data X.")}
  if(!is.numeric(Y[,1]) | !is.numeric(Y[,2]) | !is.numeric(Y[,3]) | length(Y) < 2) {
    stop("Please specify a valid input for data Y.")}
  else if (sum(is.na(Y)) > 0){
    stop("Please specify a valid input for data Y.")}
  else if (stats::var(2*Y[,1] - Y[,2] - Y[,3]) == 0){
    stop("All intra-subject contrasts are the same for data Y.")}
  # if (length(x) != length(y)){
  #   stop("Please ensure x and y have the same number of observations.")}
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
  if(!(type %in% c("welch", "student"))){
    stop("Please specify a valid type for the t-test(s).")}

  if (length(compSymm) != 1 | !(compSymm %in% c(TRUE, FALSE))){
    stop("Please provide valid logical input for compSymm.")
  }
  else if (compSymm == TRUE & type == "welch"){
    stop("type must be 'student' when compSymm is TRUE.")
  }

  if (compSymm == FALSE){
    x <- 0.25*(2*X[,1] - X[,2] - X[,3]); y <- 0.25*(2*Y[,1] - Y[,2] - Y[,3])

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

  if (compSymm == TRUE){

    n1 <- nrow(X); n2 <- nrow(Y)
    mydata11 <- data.frame(Y = X[,1],Subjects=seq(1, n1),Time=1,Groups="TRR",Treatment=2, Carryover = 0)
    mydata21 <- data.frame(Y = X[,2],Subjects=seq(1, n1),Time=2,Groups="TRR",Treatment=1, Carryover = 2)
    mydata31 <- data.frame(Y = X[,3],Subjects=seq(1, n1),Time=3,Groups="TRR",Treatment=1, Carryover = 1)
    mydata12 <- data.frame(Y = Y[,1],Subjects=seq(n1 + 1, n1 + n2),Time=1,Groups="RTT",Treatment=1, Carryover = 0)
    mydata22 <- data.frame(Y = Y[,2],Subjects=seq(n1 + 1, n1 + n2),Time=2,Groups="RTT",Treatment=2, Carryover = 1)
    mydata32 <- data.frame(Y = Y[,3],Subjects=seq(n1 + 1, n1 + n2),Time=3,Groups="RTT",Treatment=2, Carryover = 2)
    mydata <- rbind(mydata11,mydata21,mydata31,mydata12,mydata22,mydata32)
    mydata$Subjects<-factor(mydata$Subjects)
    mydata$Time<-factor(mydata$Time)

    aov.1 <- stats::aov(Y~Groups+Time+Treatment + Carryover + Error(Subjects), data=mydata)
    seCompSymm <- sqrt((3/8)*summary(aov.1)[2][[1]][[1]][4,3]*(1/n1 + 1/n2))

    x <- 0.25*(2*X[,1] - X[,2] - X[,3]); y <- 0.25*(2*Y[,1] - Y[,2] - Y[,3])
    eff <- mean(x) - mean(y)

    if (is.finite(deltaL) & is.finite(deltaU)){
      # t1 <- stats::t.test(x, y, mu = deltaL, alternative = "greater", var.equal = ifelse(type == "student", TRUE, FALSE))
      # t2 <- stats::t.test(x, y, mu = deltaU, alternative = "less", var.equal = ifelse(type == "student", TRUE, FALSE))
      tab <- data.frame(t = c((eff - deltaL)/seCompSymm, (eff - deltaU)/seCompSymm),
                        SE = c(seCompSymm, seCompSymm),
                        df = c(2*(n1 + n2 - 2), 2*(n1 + n2 - 2)),
                        p.value= c(stats::pt((eff - deltaL)/seCompSymm, 2*(n1 + n2 - 2), lower.tail = FALSE),
                                   stats::pt((eff - deltaU)/seCompSymm, 2*(n1 + n2 - 2))))
      row.names(tab) <- c("Lower Bound", "Upper Bound")
    }
    else if (is.finite(deltaL)){
      # t1 <- stats::t.test(x, y, mu = deltaL, alternative = "greater", var.equal = ifelse(type == "student", TRUE, FALSE))
      tab <- data.frame(t = (eff - deltaL)/seCompSymm, SE = seCompSymm, df = 2*(n1 + n2 - 2),
                        p.value = stats::pt((eff - deltaL)/seCompSymm, 2*(n1 + n2 - 2), lower.tail = FALSE))
      row.names(tab) <- "Lower Bound"
    }
    else {
      # t2 <- stats::t.test(x, y, mu = deltaU, alternative = "less", var.equal = ifelse(type == "student", TRUE, FALSE))
      tab <- data.frame(t = (eff - deltaU)/seCompSymm, SE = seCompSymm, df = 2*(n1 + n2 - 2),
                        p.value = stats::pt((eff - deltaU)/seCompSymm, 2*(n1 + n2 - 2)))
      row.names(tab) <- "Upper Bound"
    }

    METHOD <-  "Student's Two-Sample t-test"

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
}
