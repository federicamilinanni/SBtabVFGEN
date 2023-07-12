require(rgsl)
require(SBtabVFGEN)

test_that("event import works",{
	sink("/dev/null")
	modelDir <- system.file("extdata/CaMKII", package = "SBtabVFGEN")
	model.tsv <- dir(modelDir,pattern="[.]tsv$",full.names=TRUE)
	CaMKII <- sbtab_from_tsv(model.tsv)
	expect_type(CaMKII,"list")
	expect_equal(class(CaMKII[[1]]),"data.frame")
	source(system.file("extdata/CaMKII/CaMKIIs.R", package = "SBtabVFGEN"))
	n <- nrow(CaMKII$Parameter)
	p <- as.matrix(model$par()[1:n])
	simExperiments <- sbtab.data(CaMKII)
	expect_type(simExperiments,"list")
	vf.c <- dir(
		system.file("extdata/CaMKII", package = "SBtabVFGEN"),
		full.names=TRUE,pattern="_gvf[.]c"
	)
	expect_true(file.exists(vf.c))
	print(getwd())
	so <- rgsl::model.so(vf.c)
	print(so)
	expect_type(so,"character")
	expect_true(file.exists(so))
	modelName <- comment(CaMKII)
	comment(modelName) <- so
	## simulate
	yf <- rgsl::r_gsl_odeiv2_outer(modelName,simExperiments,p)
	file.remove(so)
	expect_type(yf,"list")
	expect_length(yf,length(simExperiments))
	## plot
	oTime <- simExperiments$SpikeSeries$outputTimes
	Ca <- yf[["SpikeSeries"]]$func[5,,1]
	expect_equal(length(Ca),length(oTime))
	sink()
	plot(oTime,Ca)
})
