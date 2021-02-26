

#' hybToGRanges
#'
#' This function is useful to turn a list of hyb data into lists of GRanges
#' It creates a list for each sample one for the left side one for the right
#' side and one for the gap in the middle.
#' @param hybList the original hybList created with readHybFiles or subsetHybList
#' @param rna The rna of interest
#' @return A list of GRanges data in hyb format
#' @examples
#'  hybGranges =  hybToGRanges(hybList,"myRNA");
#' @export
hybToGRanges = function(hybList, rna){
    seqName = rna
    hybOutput2 = hybList
    gList = list()
    for(i in 1:length(hybOutput2)){
        hybOutput = hybOutput2[[i]]
        gList[[i]] = GRangesList()
        #make a Granges from the left
        left <- GRanges(seqnames=seqName,
                        IRanges(
                            start=hybOutput$V7,
                            end=hybOutput$V8
                        ))
        names(left) <- hybOutput$V1
        
        
        
        
        #make a GRanges from the right
        right <-GRanges(seqnames=seqName,
                        IRanges(
                            start=hybOutput$V13,
                            end=hybOutput$V14
                        ))
        names(right) <- hybOutput$V1
        
        
        
        distances = GRanges(seqnames=seqName,
                            IRanges(
                                start=end(left),
                                end=start(right)
                            ))
        names(distances) <- hybOutput$V1
        
        gList[[i]][["left"]] = left
        gList[[i]][["right"]] = right
        gList[[i]][["gap"]] = distances
    }
    return(gList)
}




#' sampleChimeras
#'
#' This function samples chimeras
#' 
#' @param chimeraList list of chimeras
#' 
#' @examples
#' @export
sampleChimeras = function(chimeraList){
    
    chimeraListSampled =list()
    
    for(i in 1:length(chimeraList)){
        max =  length(chimeraList[[i]][["left"]])
        seq = c(1,max)
        if(max > 40000){seq = seq(1,max,by = 20000)
        seq = c(seq,max)}
        
        print(seq)
        chimeraListSampled[[i]] = list()
        for(j in c("left","right","gap")){
            chimeraListSampled[[i]][[j]] = list()
            for(k in 1:(length(seq)-1)){
                sample = seq[k]:seq[k+1]
                chimeraListSampled[[i]][[j]][[k]] = chimeraList[[i]][[j]][sample]
            }
        }
        
    }
    return(chimeraListSampled)
}










#' Makes a table with the coordinates of the clusters
#'
#' Does the same as printClusters but is a lot faster and does not create plots
#' of each cluster
#' @param  dir the directory that contains the *hybrids.hyb files
#' @param  clustering The output from the iGraph function cluster_walktrap for the (made with adjacency matrix input)
#' @param  highest_clusters The cluster you are interested in keeping
#' @param  left list created with hybToGRanges (but just the left section of the list)
#' @param  right list created with hybToGRanges (but just the right section of the list)
#' @return A table of clusters and coordinates
#' @examples
#' plottingList[[i]][[k]]  = printClustersFast(dir,
#' clustering, 
#' highest_clusters, 
#' chimeraList[["sample"]][["left"]], 
#' chimeraList[["sample"]][["right"]])
#' @export
printClustersFast = function(dir, clustering, highest_clusters, left, right ){
    
    plotting = GRanges()
    for(i in highest_clusters ){#:(max(as.numeric((names(table(membership(cluster3))))))-300)){
        c1c1 = names(membership(clustering)[membership(clustering) == i])
        plotting2 = GRanges()
        plotting2 = addCluster(left, c1c1, plotting2,i  ,"left")
        plotting2 = addCluster(right, c1c1, plotting2,i ,"right" )
        plotting = addCluster(left, c1c1, plotting,i  ,"left")
        plotting = addCluster(right, c1c1, plotting,i ,"right" )
    }
    return(plotting)
}



#helper function for printClusters
addCluster = function(granges, indexes, prev, cluster, type){
    x = granges[as.numeric(indexes)]
    x$cluster = cluster
    x$type = type
    prev= c(prev, x)
}

