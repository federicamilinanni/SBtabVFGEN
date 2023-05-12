test_that("unicode micro works", {
	u<-unit.from.string("µM")
	expect_equal(u$exponent,c(1,-1))
})

test_that("fractions work", {
	u<-unit.from.string("km/h")
	U <- data.frame(scale=c(3,0),multiplier=c(1,60),exponent=c(1,-1),kind=c('metre','second'))
	comment(U) <- "km/h"
	attr(U,'id') <- "km_per_h"
	expect_equal(u,U)
})
