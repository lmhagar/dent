test_that("function works",{
  expect_equal(round(DesignParallelUnequal(diff = -4, sigma1 = 15, sigma2 = 18, q = 1, deltaL = -19.2, deltaU = 19.2, targetPower = 0.8, alpha = 0.05, seed = 1, sobol = 0)$power,4), 0.8184)
})
