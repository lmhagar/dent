test_that("invalid diff",{
  expect_error(DesignParallelUnequal(diff = "a"), "Please specify a valid number for diff.")
})

test_that("invalid sigma1",{
  expect_error(DesignParallelUnequal(diff = -4, sigma1 = -2), "Please specify a valid number for sigma1.")
})

test_that("invalid sigma2",{
  expect_error(DesignParallelUnequal(diff = -4, sigma1 = 18, sigma2 = c(15,12)), "Please specify a valid number for sigma2.")
})

test_that("no endpoints",{
  expect_error(DesignParallelUnequal(diff = -4, sigma1 = 18, sigma2 = 15), "Please specify valid interval endpoints deltaL and deltaU.")
})

test_that("invalid endpoints",{
  expect_error(DesignParallelUnequal(diff = -4, sigma1 = 18, sigma2 = 15, deltaL = 8, deltaU = 6),
               "Please specify valid interval endpoints deltaL and deltaU.")
})

test_that("one endpoint invalid",{
  expect_error(DesignParallelUnequal(diff = -4, sigma1 = 18, sigma2 = 15, deltaL = "a", deltaU = 6),
               "Please specify a valid number for deltaL.")
})

test_that("invalid alpha",{
  expect_error(DesignParallelUnequal(diff = -4, sigma1 = 18, sigma2 = 15, deltaL = -19.2, deltaU = 19.2, alpha = 2),
               "Please specify a valid number for alpha.")
})

test_that("invalid targetPower",{
  expect_error(DesignParallelUnequal(diff = -4, sigma1 = 18, sigma2 = 15, deltaL = -19.2, deltaU = 19.2, alpha = 0.05, targetPower = 0),
               "Please specify a valid number for targetPower.")
})

test_that("both types of calculations",{
  expect_error(DesignParallelUnequal(diff = -4, sigma1 = 18, sigma2 = 15, deltaL = -19.2, deltaU = 19.2, alpha = 0.05, targetPower = 0.8,
                                     n1 = 10, n2 = 10),
               "Please specify valid inputs for one of the following pairs: targetPower and q or n1 and n2.")
})

test_that("invalid q",{
  expect_error(DesignParallelUnequal(diff = -4, sigma1 = 18, sigma2 = 15, deltaL = -19.2, deltaU = 19.2, alpha = 0.05, targetPower = 0.8, q = -2),
               "Please specify a valid number for q.")
})

test_that("invalid true difference",{
  expect_error(DesignParallelUnequal(diff = -20, sigma1 = 18, sigma2 = 15, deltaL = -19.2, deltaU = 19.2, alpha = 0.05, targetPower = 0.8),
               "Please ensure diff is between deltaL and deltaU.")
})

test_that("invalid n1",{
  expect_error(DesignParallelUnequal(diff = -4, sigma1 = 18, sigma2 = 15, deltaL = -19.2, deltaU = 19.2, alpha = 0.05, n1 = 1, n2 = 3),
               "Please specify a valid integer for n1.")
})

test_that("invalid n2",{
  expect_error(DesignParallelUnequal(diff = -4, sigma1 = 18, sigma2 = 15, deltaL = -19.2, deltaU = 19.2, alpha = 0.05, n1 = 4, n2 = -3),
               "Please specify a valid integer for n2.")
})

test_that("invalid seed",{
  expect_error(DesignParallelUnequal(diff = -4, sigma1 = 18, sigma2 = 15, deltaL = -19.2, deltaU = 19.2, alpha = 0.05, n1 = 4, n2 = 3, seed = c(2,3)),
               "Please specify a valid seed for random number generation.")
})

test_that("invalid sobol",{
  expect_error(DesignParallelUnequal(diff = -4, sigma1 = 18, sigma2 = 15, deltaL = -19.2, deltaU = 19.2, alpha = 0.05, n1 = 4, n2 = 3, sobol = 6),
               "Please specify a valid integer between 0 and 4 to initialize the Sobol' sequence.")
})

test_that("invalid sobol",{
  expect_error(DesignParallelUnequal(diff = -4, sigma1 = 18, sigma2 = 15, deltaL = -19.2, deltaU = 19.2, alpha = 0.05, n1 = 4, n2 = 3, plot = 8),
               "Please provide valid logical input for plot.")
})

test_that("working example 1",{
  expect_equal(round(DesignParallelUnequal(diff = -4, sigma1 = 18, sigma2 = 15, q = 1, deltaL = -19.2, deltaU = 19.2,
                                           targetPower = 0.8, alpha = 0.05, seed = 1, sobol = 0)$power,4), 0.8232)
})

test_that("working example 2",{
  expect_equal(round(DesignParallelUnequal(diff = -4, sigma1 = 18, sigma2 = 15, q = 1, deltaL = -19.2, deltaU = 19.2,
                                           n1 = 10, n2 = 10, alpha = 0.05, seed = 2, sobol = 0)$power,4), 0.5369)
})

test_that("working example 3",{
  expect_equal(round(DesignParallelUnequal(diff = -19.2, sigma1 = 18, sigma2 = 15, q = 1, deltaL = -19.2, deltaU = Inf,
                                           n1 = 10, n2 = 10, alpha = 0.05, seed = 2, sobol = 0)$type.I.error,4), 0.0491)
})

test_that("working example equal 1",{
  expect_equal(round(DesignParallelEqual(diff = -4, sigma = 15, q = 1, deltaL = -19.2, deltaU = 19.2,
                                           targetPower = 0.8, alpha = 0.05, seed = 1, sobol = 0)$power,4), 0.8242)
})

test_that("working example equal 2",{
  expect_equal(round(DesignParallelEqual(diff = -4, sigma = 15, q = 1, deltaL = -19.2, deltaU = 19.2,
                                         n1 =10, n2 = 12, alpha = 0.05, seed = 1, sobol = 0)$power,4), 0.7063)
})

test_that("working example equal 3",{
  expect_equal(round(DesignParallelEqual(diff = 19.2, sigma = 15, q = 1, deltaL = -19.2, deltaU = 19.2,
                                         n1 =10, n2 = 12, alpha = 0.05, seed = 2, sobol = 0)$type.I.error,4), 0.0500)
})

test_that("invalid mu",{
  expect_error(DesignOneSample(mu = "a"), "Please specify a valid number for mu.")
})

test_that("both calculations one sample",{
  expect_error(DesignOneSample(mu = -4, sigma = 15, deltaL = -19.2, deltaU = 19.2, targetPower = 0.8, alpha = 0.05, n = 10),
               "Please specify valid inputs for one of the following: targetPower or n.")
})

test_that("working example one 1",{
  expect_equal(round(DesignOneSample(mu = -4, sigma = 15, deltaL = -19.2, deltaU = 19.2,
                                     targetPower = 0.8, alpha = 0.05, plot = FALSE, seed = 1, sobol = 0)$power,4), 0.8115)
})

test_that("working example one 2",{
  expect_equal(round(DesignOneSample(mu = -4, sigma = 15, deltaL = -19.2, deltaU = 19.2,
                                     n = 17, alpha = 0.05, plot = FALSE, seed = 1, sobol = 0)$power,4), 0.9905)
})

test_that("working example one 3",{
  expect_equal(round(DesignOneSample(mu = -19.2, sigma = 15, deltaL = -19.2, deltaU = 19.2,
                                     n = 17, alpha = 0.05, plot = FALSE, seed = 1, sobol = 0)$type.I.error,4), 0.0501)
})

test_that("invalid rho",{
  expect_error(DesignPairedSample(diff = -4, sigma1 = 15, sigma2 = 18, rho = 1.25, deltaL = -19.2,
                                  deltaU = 19.2, targetPower = 0.8, alpha = 0.05, plot = FALSE, seed = 1, sobol = 0),
               "Please specify a valid number for rho.")
})

test_that("working example paired 1",{
  expect_equal(round(DesignPairedSample(diff = -4, sigma1 = 15, sigma2 = 18, rho = 0.25, deltaL = -19.2,
                                        deltaU = 19.2, targetPower = 0.8, alpha = 0.05, plot = FALSE, seed = 1, sobol = 0)$power,4), 0.8047)
})

test_that("working example paired 2",{
  expect_equal(round(DesignPairedSample(diff = -4, sigma1 = 15, sigma2 = 18, rho = 0.25, deltaL = -19.2,
                                        deltaU = 19.2, n = 17, alpha = 0.05, plot = FALSE, seed = 1, sobol = 0)$power,4), 0.9011)
})

test_that("working example paired 3",{
  expect_equal(round(DesignPairedSample(diff = -19.2, sigma1 = 15, sigma2 = 18, rho = 0.25, deltaL = -19.2,
                                        deltaU = 19.2, n = 17, alpha = 0.05, plot = FALSE, seed = 2, sobol = 0)$type.I.error,4), 0.0500)
})

test_that("working example crossover 2x2 unequal 1",{
  expect_equal(round(DesignCrossover2x2Unequal(diff = 0.05, sigma1 = 0.4, sigma2 = 0.3, deltaL = -0.223,
                                               deltaU = 0.223, targetPower = 0.8, q = 1, alpha = 0.05, plot = FALSE, seed = 1, sobol = 0)$power,4), 0.8213)
})

test_that("working example crossover 2x2 unequal 2",{
  expect_equal(round(DesignCrossover2x2Unequal(diff = 0.05, sigma1 = 0.4, sigma2 = 0.3, deltaL = -0.223,
                                               deltaU = 0.223, n1 = 14, n2 = 12, alpha = 0.05, plot = FALSE, seed = 2, sobol = 0)$power,4), 0.7721)
})

test_that("working example crossover 2x2 unequal 3",{
  expect_equal(round(DesignCrossover2x2Unequal(diff = -0.223, sigma1 = 0.4, sigma2 = 0.3, deltaL = -0.223,
                                               deltaU = 0.223, n1 = 14, n2 = 12, alpha = 0.05, plot = FALSE, seed = 2, sobol = 0)$type.I.error,4), 0.0496)
})

test_that("working example crossover 2x2 equal 1",{
  expect_equal(round(DesignCrossover2x2Equal(diff = 0.05, sigma = 0.4, deltaL = -0.223,
                                             deltaU = 0.223, targetPower = 0.8, q = 1.5, alpha = 0.05, plot = FALSE, seed = 1, sobol = 0)$power,4), 0.8027)
})

test_that("working example crossover 2x2 equal 2",{
  expect_equal(round(DesignCrossover2x2Equal(diff = 0.05, sigma = 0.4, deltaL = -0.223,
                                             deltaU = 0.223, n1 = 8, n2 = 12, alpha = 0.05, plot = FALSE, seed = 2, sobol = 0)$power,4), 0.4621)
})

test_that("working example crossover 2x2 equal 3",{
  expect_equal(round(DesignCrossover2x2Equal(diff = -0.223, sigma = 0.4, deltaL = -0.223,
                                             deltaU = 0.223, n1 = 8, n2 = 12, alpha = 0.05, plot = FALSE, seed = 2, sobol = 0)$type.I.error,4), 0.0489)
})

test_that("working example crossover dual unequal 1",{
  expect_equal(round(DesignCrossoverDualUnequal(diff = 0.05, sigma1 = sqrt(2*0.2^2 + 0.25^2), sigma2 = sqrt(0.2^2 + 2*0.25^2), deltaL = -0.223,
                                                deltaU = 0.223, targetPower = 0.8, q = 1.5, alpha = 0.05, plot = FALSE, seed = 1, sobol = 0)$power,4), 0.8330)
})

test_that("working example crossover dual unequal 2",{
  expect_equal(round(DesignCrossoverDualUnequal(diff = 0.05, sigma1 = sqrt(2*0.2^2 + 0.25^2), sigma2 = sqrt(0.2^2 + 2*0.25^2), deltaL = -0.223,
                                                deltaU = 0.223, n1 = 8, n2 = 12, alpha = 0.05, plot = FALSE, seed = 2, sobol = 0)$power,4), 0.8294)
})

test_that("working example crossover dual unequal 3",{
  expect_equal(round(DesignCrossoverDualUnequal(diff = 0.223, sigma1 = sqrt(2*0.2^2 + 0.25^2), sigma2 = sqrt(0.2^2 + 2*0.25^2), deltaL = -0.223,
                                                deltaU = 0.223, n1 = 8, n2 = 12, alpha = 0.05, plot = FALSE, seed = 2, sobol = 0)$type.I.error,4), 0.0494)
})

test_that("working example crossover dual equal 1",{
  expect_equal(round(DesignCrossoverDualEqual(diff = 0.05, sigma = sqrt(3)*0.25, compSymm = FALSE, deltaL = -0.223,
                                              deltaU = 0.223, targetPower = 0.8, q = (2/3), alpha = 0.05, plot = FALSE, seed = 1, sobol = 0)$power,4), 0.8086)
})

test_that("working example crossover dual equal 2",{
  expect_equal(round(DesignCrossoverDualEqual(diff = 0.05, sigma = sqrt(3)*0.25, compSymm = FALSE, deltaL = -0.223,
                                              deltaU = 0.223, n1 = 15, n2 = 10, alpha = 0.05, plot = FALSE, seed = 1, sobol = 0)$power,4), 0.8460)
})

test_that("working example crossover dual equal 3",{
  expect_equal(round(DesignCrossoverDualEqual(diff = 0.223, sigma = sqrt(3)*0.25, compSymm = FALSE, deltaL = -0.223,
                                              deltaU = 0.223, n1 = 15, n2 = 10, alpha = 0.05, plot = FALSE, seed = 2, sobol = 0)$type.I.error,4), 0.0500)
})

test_that("working example crossover dual equal 4",{
  expect_equal(round(DesignCrossoverDualEqual(diff = 0.05, sigma = sqrt(3)*0.25, compSymm = TRUE, deltaL = -0.223,
                                              deltaU = 0.223, targetPower = 0.8, q = (2/3), alpha = 0.05, plot = FALSE, seed = 1, sobol = 0)$power,4), 0.8125)
})

test_that("working example crossover dual equal 5",{
  expect_equal(round(DesignCrossoverDualEqual(diff = 0.05, sigma = sqrt(3)*0.25, compSymm = TRUE, deltaL = -0.223,
                                              deltaU = 0.223, n1 = 15, n2 = 10, alpha = 0.05, plot = FALSE, seed = 2, sobol = 0)$power,4), 0.8564)
})

test_that("working example crossover dual equal 6",{
  expect_equal(round(DesignCrossoverDualEqual(diff = 0.223, sigma = sqrt(3)*0.25, compSymm = TRUE, deltaL = -0.223,
                                              deltaU = 0.223, n1 = 15, n2 = 10, alpha = 0.05, plot = FALSE, seed = 2, sobol = 0)$type.I.error,4), 0.0501)
})
