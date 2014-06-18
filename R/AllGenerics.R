#setGeneric(
#    name="countMismatches",
#    def = function( object, tags){standardGeneric("countMismatches")}
#    )

setGeneric(
	name = "countMismatches",
	def = function(object, tags, ...) {
		df <- standardGeneric("countMismatches")
		if(!(is(df, "data.frame"))) 
			stop("countMismatches methods must return data frame objects")
		df
	})

setGeneric(
    name = "countPrimaryTag",
    def = function(object, tags, ...) {
        df <- standardGeneric("countPrimaryTag")
        if(!(is(df, "data.frame")))
            stop("countPrimaryTag method must return data frame objects")
        df
    })

setGeneric(
		   name = "setCounts<-",
		   def = function(object,value){standardGeneric("setCounts<-")})
