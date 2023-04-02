test_that("invalid data",{
  expect_error(AnalyzeOneSample(x = c(2,2,NA,3)), "Please specify a valid input for data x.")
})

test_that("no endpoints one sample",{
  expect_error(AnalyzeOneSample(x = rnorm(10)), "Please specify valid interval endpoints deltaL and deltaU.")
})

test_that("one endpoint invalid one sample",{
  expect_error(AnalyzeOneSample(x = rnorm(10), deltaL = "a", deltaU = 6),
               "Please specify a valid number for deltaL.")
})

test_that("invalid alpha one sample",{
  expect_error(AnalyzeOneSample(x = rnorm(10), deltaL = -19.2, deltaU = 19.2, alpha = -2),
               "Please specify a valid number for alpha.")
})

test_that("working example one sample equivalence",{
  expect_equal(round(AnalyzeOneSample(x = seq(1,7, 0.5), deltaL = -3, deltaU = 5, alpha = 0.05)$table$p.value[2],4), 0.0444)
})

test_that("working example one sample noninferiority",{
  expect_equal(AnalyzeOneSample(x = seq(1,7, 0.5), deltaL = -3, alpha = 0.05)$summary, "Conclude noninferiority.")
})

test_that("invalid data two sample",{
  expect_error(AnalyzeTwoSample(x = c(2,2,4,3), y = NULL), "Please specify a valid input for data y.")
})

test_that("not paired",{
  expect_error(AnalyzeTwoSample(x = c(2,2,4,3), y = c(5,6,7), deltaL = -2, alpha = 0.05, type = "paired"),
               "Please ensure x and y have the same number of observations.")
})

test_that("not paired",{
  expect_error(AnalyzeTwoSample(x = c(2,2,4,3), y = c(5,6,7), deltaL = -2, alpha = 0.05, type = "paired"),
               "Please ensure x and y have the same number of observations.")
})

test_that("invalid type two sample",{
  expect_equal(tryCatch(AnalyzeTwoSample(x = c(2,2,4,3), y = c(5,6,7,8), type = "a", deltaL = -2, deltaU = 2, alpha = 0.05),
                        error = function(x){return("Please specify a valid type for the t-test(s).")}),
               "Please specify a valid type for the t-test(s).")
})

test_that("working example two sample equivalence",{
  expect_equal(round(AnalyzeTwoSample(x = c(2,2,4,3), y = c(5,6,7,8), deltaL = -2, alpha = 0.05, type = "paired")$table$p.value, 4), 0.9823)
})

set.seed(1)
test_that("working example 2x2 crossover unequal noninferiority",{
  expect_equal(
    round(AnalyzeCrossover2x2(X = cbind(rnorm(15), rnorm(15, 0.5)), Y = cbind(rnorm(10, 0.5), rnorm(10)),
                              type = "welch", deltaL = -1, alpha = 0.05)$table$p.value, 4)
    , 0.0095)
})

test_that("same intrasubject contrasts",{
  expect_error(
    AnalyzeCrossover2x2(X = cbind(rnorm(15), rnorm(15, 0.5)), Y = cbind(seq(1,10), seq(1,10)),
                              type = "welch", deltaL = -1, alpha = 0.05)
    , "All intra-subject contrasts are the same for data Y.")
})

set.seed(2)
test_that("working example 2x2 crossover equal",{
  expect_equal(
    round(AnalyzeCrossover2x2(X = cbind(rnorm(15), rnorm(15, 0.5)), Y = cbind(rnorm(10, 0.5), rnorm(10)),
                              type = "student", deltaL = -0.223, deltaU = 0.223, alpha = 0.05)$table$p.value[1], 4)
    , 0.6420)
})

set.seed(3)
test_that("working example dual crossover unequal",{
  expect_equal(
    round(AnalyzeCrossoverDual(X = cbind(rnorm(15), rnorm(15, 0.5), rnorm(15, 0.25)),
                               Y = cbind(rnorm(10, 0.5), rnorm(10), rnorm(10, 0.75)), type = "welch",
                               deltaL = -1, deltaU = 1, alpha = 0.05, compSymm = FALSE)$table$t[1], 4)
    , 4.7274)
})

test_that("welch compSymm invalid combination",{
  expect_error(
    AnalyzeCrossoverDual(X = cbind(rnorm(15), rnorm(15, 0.5), rnorm(15, 0.25)),
                         Y = cbind(rnorm(10, 0.5), rnorm(10), rnorm(10, 0.75)), type = "welch",
                         deltaL = -1, deltaU = 1, alpha = 0.05, compSymm = TRUE)
    , "type must be 'student' when compSymm is TRUE.")
})

set.seed(4)
test_that("working example dual crossover equal no compSymm",{
  expect_equal(
    round(AnalyzeCrossoverDual(X = cbind(rnorm(15), rnorm(15, 0.5), rnorm(15, 0.25)),
                               Y = cbind(rnorm(10, 0.5), rnorm(10), rnorm(10, 0.75)), type = "student",
                               deltaL = -1, deltaU = 1, alpha = 0.05, compSymm = FALSE)$table$df[1], 4)
    , 23)
})

set.seed(5)
test_that("working example dual crossover equal compSymm",{
  expect_equal(
    round(AnalyzeCrossoverDual(X = cbind(rnorm(15), rnorm(15, 0.5), rnorm(15, 0.25)),
                               Y = cbind(rnorm(10, 0.5), rnorm(10), rnorm(10, 0.75)), type = "student",
                               deltaL = -1, deltaU = 1, alpha = 0.05, compSymm = TRUE)$table$t[1], 4)
    , 1.4494)
})

