#' @include  comradesDataSet.R comradesFoldedDataSet.R comradesClusteredDataSet.R
NULL


##################################################
###############      METHODS       ###############
##################################################

##################################################
###############      Show          ###############
##################################################


setMethod("show", "comradesDataSet", function(object) {
    cat("comradesDataSet Object \n")
    cat("RNAs Analysed - ",rnas(object), "\n")
    cat("Samples Analysed - ",sampleNames(object), "\n")
    types = c()
    for(i in names(hybFiles(object)[[rnas(object)[1]]])){
        types = c(types  , i)
    }
    cat("Raw data  - ", types, "\n")

})


setMethod("show", "comradesClusteredDataSet", function(object) {
    cat("comradesClusteredDataSet Object \n")
    cat("RNAs Analysed - ",rnas(object), "\n")
    cat("Samples Analysed - ",sampleNames(object), "\n")
    types = c()
    for(i in names(hybFiles(object)[[rnas(object)[1]]])){
        types = c(types  , i)
    }
    cat("Raw data  - ", types, "\n")
    types = c()
    for(i in names(matrixList(object)[[rnas(object)[1]]])){
        types = c(types  , i)
    }
    cat("Matrix Types - ", types, "\n")
    types = c()
    for(i in names(clusterTableList(object)[[rnas(object)[1]]])){
        types = c(types  , i)
    }
    cat("Cluster Types - ", types, "\n")
})


setMethod("show", "comradesFoldedDataSet", function(object) {
    cat("comradesFoldedDataSet Object \n")
    cat("RNAs Analysed - ",rnas(object), "\n")
    cat("Samples Analysed - ",sampleNames(object), "\n")
    types = c()
    for(i in names(hybFiles(object)[[rnas(object)[1]]])){
        types = c(types  , i)
    }
    cat("Raw data Types - ", types, "\n")
    types = c()
    for(i in names(matrixList(object)[[rnas(object)[1]]])){
        types = c(types  , i)
    }
    cat("Matrix Types - ", types, "\n")
    types = c()
    for(i in names(clusterTableList(object)[[rnas(object)[1]]])){
        types = c(types  , i)
    }
    cat("Cluster Types - ", types, "\n")
})



##################################################
###############  Accessors         ###############
##################################################
# Functions to access the attributes

# comradesDataSet
setGeneric("rnas", function(x) standardGeneric("rnas"))
setMethod("rnas", "comradesDataSet", function(x)  x@rnas)

setGeneric("hybDir", function(x) standardGeneric("hybDir"))
setMethod("hybDir", "comradesDataSet", function(x)   x@hybDir)

setGeneric("sampleTable", function(x) standardGeneric("sampleTable"))
setMethod("sampleTable", "comradesDataSet", function(x)   x@sampleTable)

setGeneric("hybFiles", function(x) standardGeneric("hybFiles"))
setMethod("hybFiles", "comradesDataSet", function(x)   x@hybFiles)

setGeneric("matrixList", function(x) standardGeneric("matrixList"))
setMethod("matrixList", "comradesDataSet", function(x)   x@matrixList)

setGeneric("group", function(x) standardGeneric("group"))
setMethod("group", "comradesDataSet", function(x)   x@group)

setGeneric("sampleNames", function(x) standardGeneric("sampleNames"))
setMethod("sampleNames", "comradesDataSet", function(x)   x@"sampleNames")


# comradesClusteredDataSet
setGeneric("clusterGrangesList", function(x) standardGeneric("clusterGrangesList"))
setMethod("clusterGrangesList", "comradesClusteredDataSet", function(x)  x@clusterGrangesList)

setGeneric("clusterTableList", function(x) standardGeneric("clusterTableList"))
setMethod("clusterTableList", "comradesClusteredDataSet", function(x)  x@clusterTableList)


# comradesFoldedDataSet
setGeneric("clusterTableFolded", function(x) standardGeneric("clusterTableFolded"))
setMethod("clusterTableFolded", "comradesFoldedDataSet", function(x)  x@clusterTableFolded)



##################################################
###############  Setters           ###############
##################################################
# for slots that need to be altered during processing


#comradesDataSet
setGeneric("matrixList<-", function(x, value) standardGeneric("matrixList<-"))
setMethod("matrixList<-", "comradesDataSet", function(x, value) {
    x@matrixList  = value
})

#comradesClusteredDataSet
setGeneric("clusterGrangesList<-", function(x, value) standardGeneric("clusterGrangesList<-"))
setMethod("clusterGrangesList<-", "comradesClusteredDataSet", function(x, value) {
    x@clusterGrangesList  = value
})

setGeneric("clusterTableList<-", function(x, value) standardGeneric("clusterTableList<-"))
setMethod("clusterTableList<-", "comradesClusteredDataSet", function(x, value) {
    x@clusterTableList  = value
})









