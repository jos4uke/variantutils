### ======================================
### Utils functions
### ======================================

clusterMapCountMismatches <- function(bamcounter_list,crosstags_list) {
	# check validity of arguments
	bclen <- length(bamcounter_list)
	if (bclen==0L)
		stop("bamcounter list length must be greater than 0")
	if (bclen != length(crosstags_list))
		stop("bamcounter object and bam tags lists must be of same length")
	if (!(any(llply(bamcounter_list,class)=="BamCounter"))) 
		stop("all elements in bamcounter list must be of class BamCounter")
	if (!(any(llply(crosstags_list, function(t){
											  if (is.null(class(t))) return(TRUE) 
											  ifelse(class(t)!="character", FALSE, TRUE)}
					)==TRUE)))
		stop("all elements in crosstags list must be NULL or of class character")

	# parallel clusterMap
	limit_cores=detectCores()/6
	cs <- ifelse(bclen<=limit_cores, bclen, limit_cores)
	cl <- makeCluster(cs, type="FORK")
	bamcounter_list <- BiocGenerics::clusterApplyLB(cl, 1:bclen, function(i){
																			bci <- bamcounter_list[[i]]
																			xtagi <- crosstags_list[[i]]
																			setCounts(bci) <- countMismatches(bci,xtagi)
																			bci
																			})
	stopCluster(cl)
	bamcounter_list
}

joinCounts <- function(bamcounter_list,freq_labels,by=NULL,type="left",match="first") {
	# check validity of arguments
	bclen <- length(bamcounter_list)
	if (bclen==0L)
		stop("bamcounter list length must be greater than 0")
	if (bclen != length(freq_labels))
		stop("bamcounter object and frequency labels lists must be of same length")
	if (!(any(llply(freq_labels, function(t){
					ifelse(class(t)!="character", FALSE, TRUE)}
	)==TRUE)))
		stop("all elements in frequency labels list must be of class character")
	if (!(any(llply(seq_len(bclen), function(i){
					ifelse(dim(bamcounter_list[[i]]@counts)[1]==0L,FALSE,TRUE)}
	)==TRUE)))
		stop("all dataframes in bamcounter objects must not be void")

	# rename "Freq" column name in each dataframe with the right label
	dfs <- llply(seq_len(bclen), function(i){
											  bci <- bamcounter_list[[i]]
											  colnames(bci@counts) <- replace(colnames(bci@counts), which(colnames(bci@counts) == "Freq"), freq_labels[[i]])
											  bci@counts
											}
	)
																			
	# join dataframes
	dfj <- join_all(dfs, by=by, type=type, match=match)
	dfj
}

filterTag <- function(bamcounter, tag, value, skip=FALSE) {
    # check validity of arguments
    if (!(class(bamcounter)=="BamCounter"))
        stop("provided bamcounter value is not of class BamCounter")
    bamlist<-bamcounter@res
    bamlistlen <- length(bamlist)
    if (bamlistlen==0L)
        stop("provided bam list length must be greater than 0")
    if (!(any(llply(list(tag, value), function(t){
                          if (is.null(class(t))) return(TRUE)
                          ifelse(!(class(t) %in% c("character","numeric")), FALSE, TRUE)}

        )==TRUE)))
        stop("provided tag, value pair must be NULL or of class character or numeric")
    if (!(is.logical(skip)))
        stop("provided skip value is not logical/boolean")

    # filter for tag/value
    if (!(is.null(class(tag))) && !(is.null(class(value)))) {
        if (!(is.null(bamlist$tag[[tag]]))) {
            # 
            fields<-names(bamlist)[which(!names(bamlist) == "tag")]
            # output only alignments with the given tag/value pair
            if (!(skip==TRUE)) {
                # idx       
                idx<-which(bamlist$tag[[tag]]==value)

            # skip alignments with the given tag/value pair
            } else {
                # idx       
                idx<-which(bamlist$tag[[tag]]!=value)
            }
            # filter fields
            bamlistfilt <- lapply(bamlist[fields], function(v) {
                                        if (!(is.vector(v)))
                                            stop("should be a vector")
                                        v[idx]
                })
            # filter tag list
            bamlistfilt$tag <- lapply(bamlist$tag, function(v) {
                                          if (!(is.vector(v)))
                                              stop("should be a vector")
                                          v[idx]
                })
            bamlistfilt
        }
    }
}

#countMismatchesAndJoin <- function(bam){
#	if (!(file.exists(bam)))
#		stop("Bam file does not exist. Check file name/path.")
#	
#	# count and join
#	p1=ScanBamParam(tag=c("NM", "NH"), what="flag")
#	p2=ScanBamParam(tag=c("NM", "NH"), what="flag", flag=scanBamFlag(isProperPair=TRUE))
#	bc1 = BamCounter(file=bam, param=p1)
#	bc2 = BamCounter(file=bam, param=p2)
#	xtags=c("NH")
#	bcl <- clusterMapCountMismatches(list(bc1,bc2),list(xtags,xtags))
#	dfj <- joinCounts(bcl,c("AllAln_Freq","PropPairAln_Freq"),by=c("NM","NH"),type="left",match="first")
#
#	write.table(dfj,file=paste(dirname(path.expand(bam)),"/",basename(bam),"_countMM.tab",sep=""),sep="\t",row.names = FALSE)
#}

countMismatchesAndJoin <- function(bam, withOnlyFirstHit=TRUE, by="NH"){
	# check validity of arguments
    if (!(file.exists(bam)))
        stop("provided Bam file does not exist. Check file name/path.")
	if (!(is.logical(withOnlyFirstHit)))
		stop("provided withOnlyFirstHit value is not logical/boolean")
	if (!(class(by)=="character") && !(nchar(by)==2))
		stop("provided 'by' value is not a 2-character string")
 
    # count and join
	HI=ifelse(withOnlyFirstHit,"HI",NULL)
	taglist=c("NM", by, HI)
    p1=ScanBamParam(tag=taglist, what="flag")
    p2=ScanBamParam(tag=taglist, what="flag", flag=scanBamFlag(isProperPair=TRUE))
    bc1 = BamCounter(file=bam, param=p1)
    bc2 = BamCounter(file=bam, param=p2)
	if (withOnlyFirstHit) {
		bc1filt <- filterTag(bc1, "HI", 1)
    	bc2filt <- filterTag(bc2, "HI", 1)
		bc1@res<-bc1filt
		bc2@res<-bc2filt
		countColnames=c("AllAln_WithOnlyFirstHit_Freq","ProperPairAln_WithOnlyFirstHit_Freq")
		outFile=paste(dirname(path.expand(bam)),"/",basename(bam),"_countMM_withOnlyFirstHit.tab",sep="")
    } else {
		countColnames=c("AllAln_WithAllHits_Freq","ProperPairAln_WithAllhHits_Freq")
		outFile=paste(dirname(path.expand(bam)),"/",basename(bam),"_countMM_withAllHits.tab",sep="")
	}
	bcl=list(bc1, bc2)
	xtags=c(by)
    bcl <- clusterMapCountMismatches(bcl,list(xtags,xtags))
    dfj <- joinCounts(bcl,countColnames,by=c("NM",by),type="left",match="first")

    write.table(dfj,file=outFile,sep="\t",row.names = FALSE)
}


