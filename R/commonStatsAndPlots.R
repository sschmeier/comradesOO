#' @include  comradesDataSet.R comradesFoldedDataSet.R comradesClusteredDataSet.R
NULL


#' getClusterStats
#'
#' 
#' @export
#' 
setGeneric("getClusterStats", function(knowClusteredCds, ...) standardGeneric("getClusterStats"))

setMethod("getClusterStats", "comradesClusteredDataSet", function(knowClusteredCds, rna)  {
    
    sampleTable = sampleTable(knowClusteredCds)
    clusters = clusterTableList(knowClusteredCds)
    
    
    
    
    for(rna in rnas(knowClusteredCds)){
        sampleTable$rna = rna
        for(type in 1:length(clusters[[rna]])){
            print(type)
            typeID = names(clusters[[rna]])[type]
            sampleTable[,typeID] = 0
            
            v = c()
            for(sample in 1:length(clusters[[rna]][[type]])){
                sampleID = sampleNames(knowClusteredCds)[sample]
                
                v = c(v, as.numeric(nrow(clusters[[rna]][[type]][[sample]])) )
                
                
            }
            sampleTable[ , typeID] = v
            
        }
    }
    return(sampleTable)
    
})





#' getDuplexStats
#'
#' 
#' @export
#' 
setGeneric("getDuplexStats", function(knowClusteredCds,rna, ...) standardGeneric("getDuplexStats"))

setMethod("getDuplexStats", "comradesDataSet", function(knowClusteredCds,rna)  {
    
    sampleTable = sampleTable(knowClusteredCds)
    clusters = hybFiles(knowClusteredCds)
    #get numbers of original hyb files
    for(rna in c(rnas(knowClusteredCds),"all")){
        sampleTable$rna = rna
        for(type in 1:length(clusters[[rna]])){
            print(type)
            typeID = names(clusters[[rna]])[type]
            sampleTable[,typeID] = 0
            
            v = c()
            for(sample in 1:length(clusters[[rna]][[type]])){
                sampleID = sampleNames(knowClusteredCds)[sample]
                
                v = c(v, as.numeric(nrow(clusters[[rna]][[type]][[sample]])) )
                
                
            }
            sampleTable[ , typeID] = v
            
        }
    }
    
    if(class(knowClusteredCds)[1] == "comradesClusteredDataSet"){
        clusters = clusterGrangesList(knowClusteredCds)
        
        for(rna in rnas(knowClusteredCds)){
            sampleTable$rna = rna
            for(type in 1:length(clusters[[rna]])){
                print(type)
                typeID = names(clusters[[rna]])[type]
                sampleTable[,typeID] = 0
                
                v = c()
                for(sample in 1:length(clusters[[rna]][[type]])){
                    sampleID = sampleNames(knowClusteredCds)[sample]
                    if(length(clusters[[rna]][[type]][[sample]]) == 2){
                        v = c(v, as.numeric(length(unlist(clusters[[rna]][[type]][[sample]]))) )
                    } else{
                        v = c(v, as.numeric(length((clusters[[rna]][[type]][[sample]]))) )
                    }
                    
                }
                sampleTable[ , typeID] = v
                
            }
        }
    }
    return(sampleTable)
    
})















# Plot Matrices
#' @export
setGeneric("plotMatrices", function(cds,type, directory,a,b,c,d,h, ...) standardGeneric("plotMatrices"))

setMethod("plotMatrices", "comradesDataSet", function(cds,type, directory,a,b,c,d,h)  {
    
    hybMatList = matrixList(cds)
    rnaS = rnas(cds)
    sampleNames = names(hybMatList[[1]][[type]])
    if(is.null(sampleNames)){
        sampleNames = 1:length(hybMatList[[1]][[type]])
    }
    
    print(sampleNames)
    
    for(sample in  1:length(sampleNames)  ){
        
        for(rna in c(rnaS)){
            
            print(paste("*** doing ", sampleNames[sample], rna, "  ***"))
            
            myCol = colorRampPalette(c("black","black",brewer.pal(9,"YlOrRd")))(13)
            
            cols = log2(max(hybMatList[[rna]][[type]][[sample]][a:b,c:d]+1))
            
            myCol = myCol[1:cols]
            
            pdf(paste(directory,"/",rna ,"_", sampleNames[sample], "-",type ,".pdf", sep = ""),
                height = h,
                width = h)
            heatmap3((log2(t(hybMatList[[rna]][[type]][[sample]][a:b,c:d]+1))),
                     col=myCol,
                     scale="none" ,
                     Rowv = NA,
                     Colv = NA,
                     useRaster = T)
            dev.off()
            
        }
    }
    
    
})



# Plot Matrices
#' @export
setGeneric("plotMatricesAverage", function(cds,type, directory,a,b,c,d,h, ...) standardGeneric("plotMatricesAverage"))

setMethod("plotMatricesAverage", "comradesDataSet", function(cds,type, directory,a,b,c,d,h)  {
    
    for(rna in rnas(cds)) {
        
        hybMatList = matrixList(cds)
        hybMatList2 = hybMatList
        for(i in c("c", "s")){
            c = 1
            for(j in group(cds)[[i]] ){
                if(length( group(cds)[[i]] ) < 2 | c == 1 ){
                    sum(hybMatList[[rna]][[type]][[ j ]])
                    hybMatList2[[rna]][[type]][[ i ]] =   hybMatList[[rna]][[type]][[ j ]]
                    print(sum(hybMatList2[[rna]][[type]][[ i ]]))
                    print(c)
                    # hybMatList2[[rna]][[type]][["s"]] =   hybMatList[[rna]][[type]][[ j ]]
                }else{
                    hybMatList2[[rna]][[type]][[ i ]] =  hybMatList2[[rna]][[type]][[ i ]] +
                        hybMatList[[rna]][[type]][[ j ]]
                    print(sum(hybMatList2[[rna]][[type]][[ i  ]] ))
                    print(c)
                    # hybMatList2[[rna]][[type]][["s"]] =  Reduce('+', hybMatList[[rna]][[type]][[ j ]] )
                }
                c = c+1
            }
            
        }
        
        sampleNames = c("s", "c")
        #c = 1
        for(sample in  sampleNames ){
            
            
            
            print(paste("*** doing ",sample, rna, "  ***"))
            
            myCol = colorRampPalette(c("black","black",brewer.pal(9,"YlOrRd")))(13)
            
            cols = log2(max(hybMatList2[[rna]][[type]][[sample]][a:b,c:d]+1))-2
            
            myCol = myCol[1:cols]
            
            pdf(paste(directory,"/",rna ,"_", sample , "-",type ,".pdf", sep = ""),
                height = h,
                width = h)
            heatmap3((log2(t(hybMatList2[[rna]][[type]][[sample]][a:b,c:d]+1))),
                     col=myCol,
                     scale="none" ,
                     Rowv = NA,
                     Colv = NA,
                     useRaster = T)
            dev.off()
            
            #  c = c + 1
        }
    }
    
    
})

#plot




