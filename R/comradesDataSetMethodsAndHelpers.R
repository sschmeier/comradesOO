#' @include  comradesDataSet.R comradesFoldedDataSet.R comradesClusteredDataSet.R
NULL



################################################################################
# Methods functions and methods relating to creation and use of 
# comradesDataSets
################################################################################


#' topTranscripts
#'
#' This method prints the top transcripts that have the most duplexes 
#' assigned
#'
#' @param cds a \code{comradesClusteredDataSet} object
#' 
#' 
#' @return Returns a \code{comradesClusteredDataSet} object
#' 
#' The 3 attributes matrixList, clusterTableList and clusterGrangesList 
#' will gain the \code{types} "superClusters" and "trimmedClusters"
#' 
#' @export
#' 
setGeneric("topTranscripts",
           function(cds,ntop=10, ...) standardGeneric("topTranscripts" ) )

setMethod("topTranscripts",
          "comradesDataSet",
          function(cds,ntop= 10)  {
              c = group(cds)[["s"]]
              vect = c()
              for(i in c){
                  vect = c(vect, 
                           hybFiles(cds)[["all"]][["all"]][[i]]$V4,
                           hybFiles(cds)[["all"]][["all"]][[i]]$V10)
              }
              x = table(vect)[order(table(vect), decreasing = T)]
              x[1:ntop]
          })

################################################################################
# Swapping and subsetting
################################################################################




#' Swap the table to ensure that the rna of interest is on the left of the table
#'
#' Swap the table to ensure that 3 prime most duplex side is ont he left of the table
#' used to make one sides heatmaps and other reasons where having the left of the table
#' coming after the right side is a problem. Different from swapHybs as it
#' ensure that BOTH duplex sides originate from the RNA of interest.
#' @param hybList the original hybList created with readHybFiles or subsetHybList
#' @param rna The rna of interest
#' @return A list of "swapped" hyb datas
#' @examples
#'  hybListSwapped =  swapHybs2(hybList,"myRNA");
#' @export
swapHybs = function(hybList, rna){
    
    hybListSwap = hybList
    for(hyb in 1:length(hybList)){
        
        hybListS  =  hybList[[hyb]][hybList[[hyb]]$V4 == rna | hybList[[hyb]]$V10 == rna,]
        
        
        
        comb = hybListS
        
        
        
        
        hybListSwap[[hyb]] = comb
    }
    
    return(hybListSwap)
}



#' Swap the table to ensure that 3 prime most duplex side is on the left of the table
#'
#' Swap the table to ensure that 3 prime most duplex side is on the left of the table
#' used to make one sides heatmaps and other reasons where having the left of the table
#' coming after the right side is a problem. Different from swapHybs as it
#' ensure that BOTH duplex sides originate from the RNA of interest.
#' @param hybList the original hybList created with readHybFiles or subsetHybList
#' @param rna The rna of interest
#' @return A list of "swapped" hyb data
#' @examples
#'  hybListSwapped =  swapHybs2(hybList,"myRNA");
#' @export
swapHybs2 = function(hybList, rna){
    
    for(hyb in 1:length(hybList)){
        
        hybList18S = hybList[[hyb]][as.character(hybList[[hyb]]$V4) == rna & as.character(hybList[[hyb]]$V10) == rna,]
        
        
        tmp1 = hybList18S[hybList18S$V7 < hybList18S$V13,]
        tmp2 = hybList18S[hybList18S$V7 > hybList18S$V13,]
        tmp2 = tmp2[,c("V1","V2","V3",
                       "V4","V11","V12",
                       "V13","V14","V15",
                       "V10","V5","V6",
                       "V7","V8","V9")]
        colnames(tmp2) = colnames(tmp1)
        comb = rbind.data.frame(tmp1, tmp2)
        
        
        
        
        hybList[[hyb]] = comb
    }
    
    return(hybList)
}




#' Swap the table to ensure that the rna of interest is on the left of the table
#'
#' Swap the table to ensure that 3 prime most duplex side is ont he left of the table
#' used to make one sides heatmaps and other reasons where having the left of the table
#' coming after the right side is a problem. Different from swapHybs as it
#' ensure that BOTH duplex sides originate from the RNA of interest.
#' @param hybList the original hybList created with readHybFiles or subsetHybList
#' @param rna The rna of interest
#' @return A list of "swapped" hyb datas
#' @examples
#'  hybListSwapped =  swapHybs2(hybList,"myRNA");
#' @export
swapHybs3 = function(hybList, rna){
    
    hybListSwap = hybList
    for(hyb in 1:length(hybList)){
        
        hybListS  =  hybList[[hyb]][hybList[[hyb]]$V4 == rna | hybList[[hyb]]$V10 == rna,]
        
        
        tmp1 = hybListS[hybListS$V4 == rna,]
        tmp2 = hybListS[hybListS$V4 != rna,]
        tmp2 = tmp2[,c("V1","V2","V3",
                       "V10","V11","V12",
                       "V13","V14","V15",
                       "V4","V5","V6",
                       "V7","V8","V9")]
        colnames(tmp2) = colnames(tmp1)
        comb = rbind.data.frame(tmp1, tmp2, stringsAsFactors = F)
        
        
        
        
        hybListSwap[[hyb]] = comb
    }
    
    return(hybListSwap)
}








