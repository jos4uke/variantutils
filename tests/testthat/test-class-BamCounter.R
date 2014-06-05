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
	mismatchTag="NM"
	setCounts(bc) <- countMismatches(bc,mismatchTag,tags)
	
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
	mismatchTag="NM"
	df <- countMismatches(bc,mismatchTag,tags)
	
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
	mismatchTag="NM"
	df <- countMismatches(bc,mismatchTag,tags)
	
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

	mismatchTag="NM"
	bcl <- clusterMapCountMismatches(list(bc1,bc2),mismatchTag,list(tags,tags))

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

#test_that("TestCountMismatchesTagsEquality: NM==XW+XV",{
#    file="/home/ldap/users/jtran/dev/R/projects/variantutils/inst/extdata/test_XW_sorted.bam"
#	#file=system.file("extdata", "test_XW_sorted.bam", package="VariantUtils", mustWork=TRUE)
#	p1=ScanBamParam(tag=c("NM", "XW", "XV","HI"), what="flag")
#    bc1 = BamCounter(file=file, param=p1)
#    tags=c("XW","XV")
#
#    mismatchTag="NM"
#    bcl <- clusterMapCountMismatches(list(bc1),mismatchTag,list(tags))
#
#    #print(paste("bcl length: ",length(bcl), sep=""), zero.print = ".")
#    expect_true(length(bcl)==1)
#    llply(seq_len(length(bcl)), function(i){
#                                            expect_match(class(bcl[[i]]),"BamCounter")
#                                            expect_true(length(row.names(bcl[[i]]@counts)) > 0)
#                                            cat("\n")
#                                            h6=bcl[[i]]@counts[1:6,]
#                                            print(h6, zero.print = ".")
#                                        }
#    )
#	bc_one<-bcl[[1]]
#	llply(seq_len(length(bc_one@counts$NM)), function(i){
#													expect_equal(as.numeric(bc_one@counts$NM[i]),sum(as.numeric(bc_one@counts$XW[i]),as.numeric(bc_one@counts$XV[i])))
#	}
#	  )
#})


context("Merging counts dataframes")

test_that("TestJoiningCounts",{
	file=system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE)
	p1=ScanBamParam(tag=c("NM", "H1"), what="flag")
	p2=ScanBamParam(tag=c("NM", "H1"), what="flag", flag=scanBamFlag(isProperPair=TRUE))
	bc1 = BamCounter(file=file, param=p1)
	bc2 = BamCounter(file=file, param=p2)
	tags=c("H1")
	mismatchTag="NM"
	bcl <- clusterMapCountMismatches(list(bc1,bc2),mismatchTag,list(tags,tags))

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
	#file=system.file("extdata", "sample_NHI5.sorted.bam", package="VariantUtils", mustWork=TRUE)
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

context("Independent Events")

test_that("TestCountIndels",{
	cigar="3S6M1D2M3I4M4S"
	matches=gregexpr(pattern="([0-9]+[ID])+", text=cigar)
	indels=unlist(regmatches(cigar, m=matches))
	indelsSizeSum=sum(as.numeric(substr(indels,1,nchar(indels)-1))) 
	indelsCountRef=c(length(indels),indelsSizeSum)

	indelsCount=countIndels(cigar)
	
	print(paste("Expected: ", paste(indelsCountRef, collapse=" ")))
	print(paste("Calculated: ", paste(indelsCount, collapse= " ")))
	
	expect_true(isTRUE(all.equal(indelsCount,indelsCountRef)))
})

#test_that("TestCalculateIndependentEventsWithoutSNP", {
#	file="/home/ldap/users/jtran/dev/R/projects/variantutils/inst/extdata/test_XW_sorted.bam"
#	#file=system.file("extdata", "test_XW_sorted.bam", package="VariantUtils", mustWork=TRUE)
#	ie_ref_file="/home/ldap/users/jtran/dev/R/projects/variantutils/inst/extdata/test_XW_sorted_countIE.txt"
#	#ie_ref_file=system.file("extdata", "", package="VariantUtils", mustWork=TRUE)
#	
#	if (!file.exists(file))
#		warning(paste(file, "file does not exist.", sep=" "))
#	if (!file.exists(ie_ref_file))
#		warning(paste(ie_ref_file, "file does not exist.", sep=" "))
#
#	withOnlyFirstHit=TRUE
#	HI=ifelse(withOnlyFirstHit,"HI",NULL)
#	taglist=c("NM", "NH", HI)
#	flags=c("cigar")
#	p1=ScanBamParam(tag=taglist, what=flags)
#	bc1 = BamCounter(file=file, param=p1)
#	bc1filt <- filterTag(bc1, "HI", 1)
#	bc1@res<-bc1filt
#	
#	ieTag="IE"
#	bc1@res$tag[[ieTag]] <- calculateIndependentEvents(bc1, ieTag="IE")
#	HI1len=length(bc1filt$tag[["HI"]])
#	IElen=length(bc1@res$tag[[ieTag]])
#	
#	print(IElen)
#	expect_true(HI1len==IElen)
#
#	# TODO: check value using the output of bash function calcIE
#	ie_ref=read.table(file=ie_ref_file, header=TRUE, sep="\t")
#	expect_true(IElen==length(na.omit(ie_ref$IE)))
#	expect_true(isTRUE(all.equal(bc1@res$tag[[ieTag]],na.omit(ie_ref$IE))))
#})
#
#
