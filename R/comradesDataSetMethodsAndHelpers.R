#' @include  comradesDataSet.R comradesFoldedDataSet.R comradesClusteredDataSet.R
NULL


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


swapHybs = function(hybList, rna){
    
    hybListSwap = hybList
    for(hyb in 1:length(hybList)){
        
        hybListS  =  hybList[[hyb]][hybList[[hyb]]$V4 == rna | hybList[[hyb]]$V10 == rna,]
        
        
        
        comb = hybListS
        
        
        
        
        hybListSwap[[hyb]] = comb
    }
    
    return(hybListSwap)
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







#' Subset a list of hyb files
#'
#' Function used to subset a list of hyb data created by readHybFiles
#' This function produces the same size list as before but with an
#' it returns ONLY the rna of interest and also
#' the number parameter specifies the number of chimeras to return
#' this is helpful when you want to directly compare libraries with
#' different sizes.
#' @param hybList the original hybList created with readHybFiles
#' @param min the rna of interest that you want to subset
#' @param max The number of randomly subsetted chimeric reads you need
#' @param length The number of randomly subsetted chimeric reads you need
#' @return A list of subsetted hyb files
#' @examples
#' hybListSubset <- subsetHybList(hybList, "myRNA" ,  10000);
#' @export



subsetHybList2 = function(hybList, min, max, length){
    longDistHyb = list()
    for (i in 1:length(hybList)){
        hybList[[i]]$dist = hybList[[i]]$V13 -  hybList[[i]]$V8
        longDistHyb[[i]] = hybList[[i]][hybList[[i]]$dist < max & hybList[[i]]$dist >= min,]
        leftLength = longDistHyb[[i]]$V6 - longDistHyb[[i]]$V5
        rightLength = longDistHyb[[i]]$V12 - longDistHyb[[i]]$V11
        longDistHyb[[i]] = longDistHyb[[i]][leftLength < length & rightLength < length,]
    }
    return(longDistHyb)
}
