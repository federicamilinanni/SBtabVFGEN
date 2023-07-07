require(rgsl)
require(SBtabVFGEN)

test_that("event import works",{
	modelDir <- system.file("extdata/CaMKII", package = "SBtabVFGEN")
	model.tsv <- dir(modelDir,pattern="[.]tsv$",full.names=TRUE)
	model <- sbtab_from_tsv(model.tsv)
	expect_type(model,"list")
	expect_equal(class(model[[1]]),"data.frame")
	p <- as.matrix(sbtab_quantity(model$Parameter))
	simExperiments <- sbtab.data(model)
	expect_type(simExperiments,"list")
	vf.c <- dir(
		system.file("extdata/CaMKII", package = "SBtabVFGEN"),
		full.names=TRUE,pattern="_gvf[.]c"
	)
	expect_true(file.exists(vf.c))
	so <- rgsl::model.so(vf.c)
	print(so)
	expect_type(so,"character")
	expect_true(file.exists(so))
	modelName <- comment(model)
	comment(modelName) <- so
	## simulate
	yf <- rgsl::r_gsl_odeiv2_outer(modelName,simExperiments,p)
	expect_type(yf,"list")
	expect_length(yf,length(simExperiments))
	## plot
	oTime <- simExperiment$SpikeSeries$outputTimes
	Ca <- yf[["SpikeSeries"]]$func[5,]
	expect_equal(Ca,length(oTime))
	plot(oTime,Ca)
})
