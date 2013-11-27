library(testthat)
library(VariantAnnotation)
library(Rsamtools)

sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
        if(trace) cat(nm,":")
        source(file.path(path, nm), ...)
        if(trace) cat("\n")
    }
}

sourceDir("../../R/", trace=FALSE)

context("BamCounter class")

test_that("TestBamCounterInstanceCreation",{
    file=system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE)
    param=ScanBamParam(tag=c("NM", "H1"), what="flag")
    bc = BamCounter(file=file, param=param)

    expect_match(class(bc),"BamCounter")
    expect_match(class(bc@file),"BamFile")
    expect_match(class(bc@param),"ScanBamParam")
    expect_match(class(bc@res),"list")
    expect_true(length(bc@res[["tag"]]$NM) > 0)
    expect_match(class(bc@counts),"data.frame")
})

test_that("TestBamCounterSetCounts",{
	file=system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE)
    param=ScanBamParam(tag=c("NM", "H1"), what="flag")
    bc = BamCounter(file=file, param=param)
	tags=NULL
	setCounts(bc) <- countMismatches(bc,tags)
	
	expect_match(class(bc@counts),"data.frame")
	expect_true(length(row.names(bc@counts)) > 0) 
	cat("\n")
	h2=bc@counts[1:2,]
	print(h2, zero.print = ".")
})

test_that("TestCountMismatchesOnly",{
	file=system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE)
    param=ScanBamParam(tag=c("NM", "H1"), what="flag")
    bc = BamCounter(file=file, param=param)
	tags=NULL
	df <- countMismatches(bc,tags)
	
	expect_match(class(df),"data.frame")	
	cat("\n")
	h3=df[1:3,]
	print(h3, zero.print = ".")
})

test_that("TestCountMismatchesCrossedTags",{
	file=system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE)
    param=ScanBamParam(tag=c("NM", "H1"), what="flag")
    bc = BamCounter(file=file, param=param)
	tags=c("H1")
	df <- countMismatches(bc,tags)
	
	expect_match(class(df),"data.frame")	
	cat("\n")
	h4=df[1:4,]
	print(h4, zero.print = ".")
})


