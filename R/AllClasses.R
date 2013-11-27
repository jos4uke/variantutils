setClass("BamCounter",
		representation(
			file="BamFile",
			param="ScanBamParam",
			res="list",
			counts="data.frame"
			),
		 sealed = FALSE
		)

BamCounter<-function(file=NA_character_, param=NA_ScanBamParam_){ 
    new("BamCounter", file=file, param=param)                                                                                                                        }
