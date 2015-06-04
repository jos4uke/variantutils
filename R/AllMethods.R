###############
### BamCounter
###############

setMethod("initialize",
          signature(.Object="BamCounter"),
          function(.Object, file, param){ 
            .Object@file = BamFile(file=file) 
            .Object@param = param
            .Object@res = scanBam(file=.Object@file, param=.Object@param)[[1]] 
            validObject(.Object) 
			return(.Object)
          }
)

setMethod(
		  f = "countMismatches",
		  signature = "BamCounter",
		  definition = function(object, mismatchTag="NM",tags){
			 if (!(is.null(tags))) {
				 if (!(class(mismatchTag)=="character") && !(nchar(mismatchTag)==2))
					stop("provided 'mismatchTag' value is not a 2-character string")
				 if (!(is(tags,"character") && any(nchar(tags)>0)))
					stop("tags must be a character vector and each tag element must be non-empty string")

				# cross NM tag with given other tags
				df <- as.data.frame(xtabs(as.formula(paste(paste("~",mismatchTag,sep=""),tags, sep="+",collapse="+")), data=object@res[["tag"]]))
				return(df)
			 } else {
				# only NM tag
				df <- as.data.frame(xtabs(paste("~",mismatchTag,sep=""), data=object@res[["tag"]]))
				return(df)
			 }
		  })

setMethod(
		  f = "countPrimaryTag",
		  signature = "BamCounter",
		  definition = function(object, primaryTag="IE", tags){
			if (!(is.null(tags))) {
                 if (!(class(primaryTag)=="character") && !(nchar(primaryTag)==2))
                    stop("provided 'primaryTag' value is not a 2-character string")
                 if (!(is(tags,"character") && any(nchar(tags)>0)))
                    stop("tags must be a character vector and each tag element must be non-empty string")
				 if (is.null(object@res[["tag"]][[primaryTag]]))
					stop(paste(primaryTag, "tag does not exist in BamCounter object", sep=" "))

				 # cross primaryTag with given other tags
				 df <- as.data.frame(xtabs(as.formula(paste(paste("~",primaryTag,sep=""),tags, sep="+",collapse="+")), data=object@res[["tag"]]))
				 df
			} else {
				# only primaryTag
				df <- as.data.frame(xtabs(paste("~",primaryTag,sep=""), data=object@res[["tag"]]))
				df
			}
		  })

setReplaceMethod(
				 f = "setCounts",
				 signature = "BamCounter",
				 definition = function(object,value){
					 object@counts = value
					 validObject(object)
					 return(object)
				 }
				)
