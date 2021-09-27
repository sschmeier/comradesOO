#' @include  comradesDataSet.R comradesFoldedDataSet.R comradesClusteredDataSet.R
NULL



################################################################################
# Methods functions and methods relating to creation and use of 
# comradesDataSets
################################################################################


#' featureInfo
#'
#' this method prints some feature level statistics
#'
#' @param cds a \code{comradesClusteredDataSet} object
#' @param file a folder location to write the table of statistics for RNAS
#' 
#' @export
#' 
setGeneric("featureInfo",
           function(cds, file, ...) standardGeneric("featureInfo" ) )

setMethod("featureInfo",
          "comradesDataSet",
          function(cds,file)  {
              
              
              alteredHybList = list()
              
              for(TE in rnas(cds)){
                  
                  hybList = getData(cds,TE, data = "hybFiles", "original" )
                  
                  for(hyb in 1:length(hybList)){
                      
                      controlHyb = hybList[[hyb]]
                      #get the right columns
                      controlHybO  = as.data.frame(cbind(as.character(controlHyb$V4),
                                                         as.character(controlHyb$V10)))
                      
                      #get only those lines with the RNA
                      controlHyb = controlHybO[controlHybO$V1 == TE | controlHybO$V2 == TE,]
                      
                      #remove the factor levels 
                      controlHyb$V1 = as.character(controlHyb$V1)
                      controlHyb$V2 = as.character(controlHyb$V2)
                      
                      
                      #Swap the columns around that do not have zika in column 2 
                      controlHybtmp1 = controlHyb[controlHyb$V1 == TE,]
                      controlHybtmp2 = controlHyb[!(controlHyb$V1 == TE),]
                      controlHybtmp2 = rev(controlHybtmp2)
                      colnames(controlHybtmp2) = c("V1","V2")
                      controlHyb = as.data.frame(rbind(controlHybtmp1, controlHybtmp2))
                      
                      # check to see if the unique stuff has reduced 
                      print(unique(sort(controlHyb$V1)))
                      print(length(unique(sort(controlHyb$V2))))
                      
                      
                      #add to list
                      alteredHybList[[TE]][[hyb]] = controlHyb
                  }
              }
              
              
              
              
              ##############################
              # Now Combine the two  
              ##############################
              
              
              
              df = list()

              for(TE in rnas(cds)){
                  df[[TE]] = data.frame()
                  aggList = list()
                  totalNames = c()
                  for(hyb in 1:length(alteredHybList[[TE]])){
                      #aggList[[TE]] = list()
                      sampleHyb = alteredHybList[[TE]][[hyb]]
                      sData = sampleHyb[sampleHyb$V1 == TE,]
                      freqSample2 = aggregate(sData$V1, by = list(sData$V2), FUN = length)
                      freqSample = freqSample2$x
                      names(freqSample) = freqSample2$Group.1
                      aggList[[TE]][[hyb]] = freqSample
                      # Get the total features that exist in the dataset 
                      totalNames = unique(sort(c(totalNames, names(freqSample))))
                  }
                  
                  tmpMat = matrix(0, nrow = length(totalNames), ncol = length(hybList))
                  row.names(tmpMat) = totalNames
                  colnames(tmpMat) = sampleNames(cds)
                  
                  
                  
                  for(i in 1:nrow(tmpMat)){
                      for(j in 1:ncol(tmpMat)){
                          tmpMat[i,j] = aggList[[TE]][[j]][row.names(tmpMat)[i]]
                      }
                      #print(i)
                  }
                  
                  tmpMat[is.na(tmpMat)] <- 0.00001
                  write.table(tmpMat, file = paste(file,"/",TE,"mRNA_featureStats.txt",sep = ""), quote = F, sep = "\t")
                  df[[TE]] = as.data.frame(tmpMat)
              }
              
              
              
              
              
              
              
              
              
              
              
              statsList = list()
              
              for(TE in rnas(cds)){
                  statsList[[TE]] = df[[TE]]
                  statsList[[TE]]$ID = sapply(row.names(statsList[[TE]]), function(x) strsplit(x, "_")[[1]][4], USE.NAMES=FALSE)
                  statsList[[TE]]
              }
              
              
              # Now get the stats from aggregate 
              
              aggList = list()
              aggList2 = list()
              
              for(TE in rnas(cds)){
                  
                  aggList[[TE]] =   aggregate(.~ID , statsList[[TE]], sum)
                  
                  
                  #get sample Table
                  st = sampleTable(cds)
                  
                  #find out which samples have controls
                  print("finding sample with controls")
                  samples = st[which(duplicated(st$sample)), "sample"]
                  
                  
                  
                  for(i in samples){
                      
                      # now get the control and sample index
                      control =  as.numeric(which(st$sample == i & st$group == "c")) +1
                      sample = as.numeric(which(st$sample == i & st$group == "s"))+1
                      
                      
                      factor1 = sum(aggList[[TE]][,sample]) /sum(aggList[[TE]][,control])
                      
                      aggList[[TE]][,paste("norm",i, sep = "_")] = aggList[[TE]][,sample] / (aggList[[TE]][,control]* factor1)
                      
                      
                  }
                  
                  data = melt(aggList[[TE]],id.vars = c("ID") )
                  
                  samples = paste("norm",samples, sep = "_")
                  

                  
                  plot(ggplot() + 
                           geom_boxplot(data = data[data$variable %in% samples,], aes(x = reorder(ID, value, FUN = median), y = log2(value))) +
                           geom_point(data = data[data$variable %in% samples,], aes(x = reorder(ID, value, FUN = median), y = log2(value), fill = ID)) +
                           geom_hline(yintercept = 0, colour = "darkred")+
                           theme_classic() +
                           theme(axis.text.x = element_text(angle = 90, hjust = 1)))

                  
              }
              
              
              
              
          })





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








