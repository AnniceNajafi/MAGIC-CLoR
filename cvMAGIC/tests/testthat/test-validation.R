test_that("validate_split_p accepts interior values and rejects boundaries/garbage", {
  expect_equal(validate_split_p(0.5), 0.5)
  expect_equal(validate_split_p(0.1), 0.1)
  expect_error(validate_split_p(0), "between 0 and 1")
  expect_error(validate_split_p(1), "between 0 and 1")
  expect_error(validate_split_p(-0.1), "between 0 and 1")
  expect_error(validate_split_p(1.5), "between 0 and 1")
  expect_error(validate_split_p(NA_real_), "between 0 and 1")
  expect_error(validate_split_p(c(0.3, 0.5)), "between 0 and 1")
  expect_error(validate_split_p("0.5"), "between 0 and 1")
})

test_that("validate_t_values canonicalises valid grids and rejects invalid ones", {
  expect_equal(validate_t_values(0:10), 0:10)
  expect_equal(validate_t_values(c(3, 1, 1, 2)), c(1L, 2L, 3L))
  expect_equal(validate_t_values(c(2, 0)), c(0L, 2L))
  expect_error(validate_t_values(c(-1, 2)), "non-negative integers")
  expect_error(validate_t_values(2.5), "non-negative integers")
  expect_error(validate_t_values(c(1, NA)), "finite")
  expect_error(validate_t_values(integer(0)), "non-empty")
  expect_error(validate_t_values("1"), "non-empty|integers")
})

test_that("public MCV entry points fail fast on invalid split_p and t_values", {
  set.seed(1)
  X <- matrix(rpois(60 * 40, 2), nrow = 60, ncol = 40)
  gp <- list(npca = 5L, k = 8L, ka = 4L)

  # valid defaults still run end-to-end
  r <- mcv_select_t(X, t_values = 0:4, npca = gp$npca, k = gp$k, ka = gp$ka)
  expect_true(r$t %in% 0:4)

  expect_error(mcv_select_t(X, t_values = 0:4, npca = gp$npca, k = gp$k,
                            ka = gp$ka, split_p = 0), "between 0 and 1")
  expect_error(mcv_select_t(X, t_values = -1:3, npca = gp$npca, k = gp$k,
                            ka = gp$ka), "non-negative integers")
  expect_error(mcv_select_t_tolerant(X, t_values = 0:4, npca = gp$npca,
                                     k = gp$k, ka = gp$ka, split_p = 1),
               "between 0 and 1")
  expect_error(magic_per_gene(X, t_values = c(1, 2.5), npca = gp$npca,
                              k = gp$k, ka = gp$ka), "non-negative integers")
  expect_error(magic_modules(X, n_modules = 3L, t_values = 1:4, npca = gp$npca,
                             k = gp$k, ka = gp$ka, split_p = 0),
               "between 0 and 1")
})
