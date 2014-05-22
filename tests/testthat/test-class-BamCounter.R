library(testthat)
library(VariantAnnotation)
library(Rsamtools)
library(plyr)

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

context("Parallel clusterMap")

test_that("TestParallelClusterMapCountMismatches",{
	file=system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE)
	p1=ScanBamParam(tag=c("NM", "H1"), what="flag")
	p2=ScanBamParam(tag=c("NM", "H1"), what="flag", flag=scanBamFlag(isProperPair=TRUE))
	bc1 = BamCounter(file=file, param=p1)
	bc2 = BamCounter(file=file, param=p2)
	tags=c("H1")

	bcl <- clusterMapCountMismatches(list(bc1,bc2),list(tags,tags))

	#print(paste("bcl length: ",length(bcl), sep=""), zero.print = ".")
	expect_true(length(bcl)==2)
	llply(seq_len(length(bcl)), function(i){
											expect_match(class(bcl[[i]]),"BamCounter")
											expect_true(length(row.names(bcl[[i]]@counts)) > 0)
											cat("\n")
											h5=bcl[[i]]@counts[1:5,]
											print(h5, zero.print = ".")
										}
	)
})

context("Merging counts dataframes")

test_that("TestJoiningCounts",{
	file=system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE)
	p1=ScanBamParam(tag=c("NM", "H1"), what="flag")
	p2=ScanBamParam(tag=c("NM", "H1"), what="flag", flag=scanBamFlag(isProperPair=TRUE))
	bc1 = BamCounter(file=file, param=p1)
	bc2 = BamCounter(file=file, param=p2)
	tags=c("H1")

	bcl <- clusterMapCountMismatches(list(bc1,bc2),list(tags,tags))

	dfj <- joinCounts(bcl,c("AllAln_Freq","PropPairAln_Freq"),by=c("NM","H1"),type="left",match="first")

	expect_match(class(dfj),"data.frame")
	expect_true(dim(dfj)[1]>0L)
	expect_true(length(colnames(dfj))==4L)

	cat("\n")
	h6=dfj[1:6,]
	print(h6, zero.print = ".")
})

context("Filtering tags")

test_that("TestFilteringTagsHI:i:1",{
	file="/home/ldap/users/jtran/dev/R/projects/variantutils/inst/extdata/sample_NHI5.sorted.bam"
	#file="../../data/sample_NHI5.sorted.bam"
	if (!file.exists(file))
		warning(paste(file, "file does not exist.", sep=" "))

	withOnlyFirstHit=TRUE
	HI=ifelse(withOnlyFirstHit,"HI",NULL)
    taglist=c("NM", "NH", HI)
    p1=ScanBamParam(tag=taglist, what="flag")
	bc1 = BamCounter(file=file, param=p1)
	bc1filt <- filterTag(bc1, "HI", 1)
	
	str(bc1filt)
	HI1len=length(bc1filt$tag[["HI"]])
	print(HI1len)
	expect_true(HI1len==19L)
})
