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
		  definition = function(object, tags){
			 if (!(is.null(tags))) {
				if (!(is(tags,"character") && any(nchar(tags)>0)))
					stop("tags must be a character vector and each tag element must be non-empty string")
				# cross NM tag with given other tags
				df <- as.data.frame(xtabs(as.formula(paste("~NM",tags, sep="+")), data=object@res[["tag"]]))
				return(df)
			 } else {
				# only NM tag
				df <- as.data.frame(xtabs(~NM, data=object@res[["tag"]]))
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
