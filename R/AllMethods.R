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

setReplaceMethod(
				 f = "setCounts",
				 signature = "BamCounter",
				 definition = function(object,value){
					 object@counts = value
					 validObject(object)
					 return(object)
				 }
				)
