test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("unicode micro works", {
	u<-unit.from.string("µM")
	expect_equal(u$exponent,c(1,-1))
})

test_that("fractions work", {
	u<-unit.from.string("km/h")
	expect_equal(u,data.frame(scale=c(3,0),multiplier=c(1,60),exponent=c(1,-1)))
})
