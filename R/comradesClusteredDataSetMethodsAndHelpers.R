#' @include  comradesDataSet.R comradesFoldedDataSet.R comradesClusteredDataSet.R
NULL

#' trimClusters
#'
#' This method trims the clusters
#'
#' @param clusteredCds a \code{comradesClusteredDataSet} object
#' 
#' 
#' @return Returns a \code{comradesClusteredDataSet} object
#' 
#' The 3 attributes matrixList, clusterTableList and clusterGrangesList 
#' will gain the \code{types} "superClusters" and "trimmedClusters"
#' 
#' @export
#' 
setGeneric("trimClusters",
           function(clusteredCds, ...) standardGeneric("trimClusters" ) )

setMethod("trimClusters",
          "comradesClusteredDataSet",
          function(clusteredCds)  {
              
              ##############################
              # set up variables
              allChimerasForSuperClustersPlotting = list()
              for(rna in rnas(clusteredCds)){  ## for each RNA
                  
                  # size of rna
                  rnaSize = ncol(matrixList(clusteredCds)[[rna]][["noHost"]][[1]])
                  # original clusters
                  originalClusters =  clusterTableList(clusteredCds)[[rna]][["original"]]
                  
                  
                  ##############################
                  # Now cluster the clusters
                  #get original tables
                  clusterPositionsList = clusterTableList(clusteredCds)[[rna]][["original"]]
                  #get original gRanges
                  combinedPlotting     = clusterGrangesList(clusteredCds)[[rna]][["original"]]
                  
                  ##############################
                  # Set up new tables, matric and granges lists
                  superclustersPoisitonList = list()
                  superclustersPlotting = list()
                  matList = list()
                  
                  for(z in 1:length(sampleNames(clusteredCds))){
                      
                      clusterPositions = clusterPositionsList[[z]]
                      ##############################
                      # changes coordinates of clusters where the 2 sides overlap
                      clusterPositions2 = clusterPositions
                      for(i in 1:nrow(clusterPositions )){
                          if(clusterPositions$le[i] > clusterPositions$rs[i]){
                              clusterPositions2[i,"rs"] =     clusterPositions[i,"le"] +1
                          }
                      }
                      ##############################
                      #make Granges left right and gap
                      left = GRanges(seqnames=rna,
                                     IRanges(start=clusterPositions2$ls,
                                             end=clusterPositions2$le))
                      names(left) <- clusterPositions2$id
                      right= GRanges(seqnames=rna,
                                     IRanges(start=clusterPositions2$rs,
                                             end=clusterPositions2$re))
                      names(right) <- clusterPositions2$id
                      distances = GRanges(seqnames=rna,
                                          IRanges(start=clusterPositions2$le,
                                                  end=clusterPositions2$rs))
                      names(distances) <- clusterPositions2$id
                      
                      
                      ##############################
                      # Now make super clusters
                      ##############################
                      
                      # from the gaps, make a adjacancy matrix
                      adjacancyMat = getAdjacancyMat(distances,"nucleotide", 35)
                      # create Graph
                      net = graph_from_adjacency_matrix(adjacancyMat,
                                                        mode = "undirected",
                                                        weighted = T)
                      # clusterGraph
                      clustering = cluster_walktrap(net,steps = 1)
                      
                      ##############################
                      # Store Super-clusters
                      highest_clusters = names(table(membership(clustering)))
                      # printClustersFast function creates a the standard clustering
                      # table from the iGraph output
                      superclustersPlotting[[z]]  = printClustersFast("../clustering/combined/",clustering, highest_clusters, left, right)
                      plottingListFull = superclustersPlotting[[z]]
                      
                      ##############################
                      # Identify orphan clusters that missed with superclustering
                      missing = as.character(clusterPositions$id[which( !(as.character(clusterPositions$id) %in% unique(names(plottingListFull)) ) )])
                      clusterPositionsmissing = clusterPositions[clusterPositions$id %in% missing,]
                      
                      ##############################
                      #  Get super cluster and cluster identity
                      cluster = mcols(plottingListFull)$cluster
                      names(cluster)= names(plottingListFull)
                      #  add this super cluster membership to the clustering table
                      clusterPositions$superCluster = cluster[as.character((clusterPositions$id))]
                      clusterPositions = clusterPositions[!is.na(clusterPositions$superCluster),]
                      # no find the number of chimeras in each supercluster
                      lengths = aggregate(clusterPositions$size.x, by = list(clusterPositions$superCluster), FUN = sum)
                      row.names(lengths) = lengths$Group.1
                      
                      ##############################
                      # make Table
                      #for each cluster get the min start and max end
                      plottingSplit = split(plottingListFull, paste(mcols(plottingListFull)$cluster, mcols(plottingListFull)$type))
                      
                      #returns the min start and max ends of each cluster
                      minStarts = unlist(lapply(plottingSplit, function(x) {return(min(start(x)))  }))
                      maxEnd = unlist(lapply(plottingSplit, function(x) {return(max(end(x)))  }))
                      # Make the clustering table
                      clusterPositionsCombined = data.frame("id" = names(maxEnd)[seq(1,length(minStarts),2)],
                                                            "ls" = minStarts[seq(1,length(minStarts),2)],
                                                            "le" = maxEnd[seq(1,length(maxEnd),2)],
                                                            "rs" = minStarts[seq(2,length(minStarts),2)],
                                                            "re" = maxEnd[seq(2,length(maxEnd),2)],
                                                            "size" = lengths[as.numeric(sub("\\s.*","",names(maxEnd)[seq(1,length(minStarts),2)])),])
                      
                      # Make the clustering table
                      superclustersPoisitonList[[z]] = rbind.data.frame(clusterPositionsmissing,
                                                                        clusterPositionsCombined,
                                                                        stringsAsFactors = F)
                      
                      ##############################
                      # Make matrices of superclusters
                      clusterPositions = superclustersPoisitonList[[z]]
                      mat = matrix(0,nrow = rnaSize, ncol = rnaSize)
                      for(i in 1:nrow(clusterPositions)){
                          mat[clusterPositions[i,"ls"]:clusterPositions[i,"le"],
                              clusterPositions[i,"rs"]:clusterPositions[i,"re"]] =    mat[clusterPositions[i,"ls"]:clusterPositions[i,"le"],
                                                                                          clusterPositions[i,"rs"]:clusterPositions[i,"re"]] + clusterPositions[i, "size.x"]
                          
                      }
                      matList[[z]] = mat
                  }
                  ##############################
                  # Save new Table and matrix list -  super clusters
                  print("saving")
                  print("saving mat list")
                  ml = matrixList(clusteredCds)
                  ml[[rna]][["superClusters"]] = matList
                  
                  print("saving table list")
                  ctl = clusterTableList(clusteredCds)
                  ctl[[rna]][["superClusters"]] =   superclustersPoisitonList
              }
              
              
              ##############################
              # Get Granges List for the super clusters (containing original duplexes)
              
              # **
              allChimerasForSuperClustersPlotting = list()
              combinedPlottingSplit = list()
              combinedPlottingUnlist = list()
              
              for(b in 1:length(sampleNames(clusteredCds))){
                  
                  combinedPlottingUnlist = unlist(combinedPlotting[[b]])
                  combinedPlottingSplit = split(combinedPlottingUnlist,
                                                mcols(combinedPlottingUnlist)$cluster)
                  superClusterArray =  sub("\\s.*","",names(superclustersPlotting[[b]][duplicated(names(superclustersPlotting[[b]]))]))
                  names(superClusterArray) = superclustersPlotting[[b]][duplicated(names(superclustersPlotting[[b]]))]$cluster
                  x = superclustersPoisitonList[[b]][ grep("bin", row.names(superclustersPoisitonList[[b]])),]
                  names = c(names(superClusterArray), row.names(x))
                  superClusterArray = c(superClusterArray,row.names(x))
                  names(superClusterArray) = names
                  combinedPlottingUnlist$superCluster = "X"
                  for( z in 1:length(superClusterArray)){
                      supercluster = names(superClusterArray)[z]
                      cluster = unique(sub("\\s.*","",superClusterArray[z]))
                      combinedPlottingUnlist[combinedPlottingUnlist$k == cluster,]$superCluster = supercluster
                  }
                  allChimerasForSuperClustersPlotting[[b]] = combinedPlottingUnlist
              }
              ##############################
              # Save new Granges super clusters
              print("saving granges list")
              cgr = clusterGrangesList(clusteredCds)
              cgr[[rna]][["superClusters"]]   =  allChimerasForSuperClustersPlotting
              
              
              
              
              
              
              ##############################
              # Now Trim the clusters
              # the new granges list
              allChimerasForSuperClustersPlottingTrimmed = list()
              # for each sample
              for(i in 1:length(sampleNames(clusteredCds))){
                  allChimerasForSuperClustersPlottingTrimmed[[i]] = GRanges()
                  #for each cluster
                  for(cluster in unique(allChimerasForSuperClustersPlotting[[i]]$superCluster)){
                      cluster2 = sub("\\s.*","" ,cluster)
                      lefty = list()
                      
                      #for the left and right sides for each cluster
                      # cut the ends based mean and sd of evidence
                      for(l in c("left","right")){
                          clusterrange = allChimerasForSuperClustersPlotting[[i]][allChimerasForSuperClustersPlotting[[i]]$superCluster == cluster & allChimerasForSuperClustersPlotting[[i]]$type == l  ,]
                          min = min(start(clusterrange[clusterrange$superCluster == cluster,]))
                          max = max(end(clusterrange[clusterrange$superCluster == cluster,]))
                          s = (start(clusterrange[clusterrange$superCluster == cluster &clusterrange$type == l  ,]))
                          e = (end(clusterrange[clusterrange$superCluster == cluster &clusterrange$type == l  ,]))
                          
                          # function that vectorises Seq
                          seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
                          # get a vector of each cluster and side from start to end
                          # the more eveidence the more times a number appears
                          x = unlist(c(seq2(from = (s), to = e)))
                          xt = table(x)
                          x1 = x
                          # find rhe mean + and - one sd and
                          # make a GRanges with this value
                          removal =  mean(x) + sd(x)*3
                          included  = GRanges(seqnames=rna,
                                              IRanges(
                                                  start=rep(min, length(clusterrange)),
                                                  end=rep(removal, length(clusterrange))
                                              ))
                          if( l == "left"){
                              removal =  mean(x) -  sd(x)*3
                              included = GRanges(seqnames=rna,
                                                 IRanges(
                                                     start=rep(removal, length(clusterrange)),
                                                     end=rep(max, length(clusterrange))
                                                 ))
                          }
                          
                          # now remove any bases that overlap with that
                          t = pintersect( clusterrange,included)
                          
                          allChimerasForSuperClustersPlottingTrimmed[[i]] = c(  allChimerasForSuperClustersPlottingTrimmed[[i]],t)
                          # un comment to print the views of the trimming
                          # s = (start(t))
                          #    e = (end(t))
                          #    x = unlist(c(seq2(from = s, to = e)))
                          #    if( l == "left"){
                          #        lefty[[1]] = x
                          #            lefty[[2]] = x1
                          #        }
                          # }
                          
                          #print(cluster )
                          #tbl1 = data.frame(table(c(x1,lefty[[2]])))
                          #tbl2 = data.frame(table(c(x,lefty[[1]])))
                          #plot(ggplot(mapping =  aes(x = Var1, y = as.numeric(as.character(Freq))))+
                          #       geom_bar(data = tbl1, stat = "identity")+
                          #       geom_bar(data = tbl2, stat = "identity", colour = "firebrick") +
                          #       theme_classic())
                      }
                  }
              }
              
              
              
              
              
              ##############################
              # From the trimmed Granges make a cluster table
              # foir the trimmed super clusters
              matListTrimmed = list()
              clusterPositionsListTrimmed = list()
              for(j in 1:length(sampleNames(clusteredCds))){
                  
                  plotting  = allChimerasForSuperClustersPlottingTrimmed[[j]]
                  lengths = aggregate(mcols(plotting)$superCluster, by = list(mcols(plotting)$superCluster), FUN = length)
                  row.names(lengths) = lengths$Group.1
                  plottingSplit = split(plotting, paste(mcols(plotting)$superCluster, mcols(plotting)$type))
                  
                  minStarts = unlist(lapply(plottingSplit, function(x) {return(min(start(x)))  }))
                  maxEnd = unlist(lapply(plottingSplit, function(x) {return(max(end(x)))  }))
                  x = sub("\\sleft","",names(maxEnd)[seq(1,length(minStarts),2)])
                  x = sub("\\sright","",x,2)
                  clusterPositionsListTrimmed[[j]] = data.frame("id" = names(maxEnd)[seq(1,length(minStarts),2)],
                                                                "ls" = minStarts[seq(1,length(minStarts),2)],
                                                                "le" = maxEnd[seq(1,length(maxEnd),2)],
                                                                "rs" = minStarts[seq(2,length(minStarts),2)],
                                                                "re" = maxEnd[seq(2,length(maxEnd),2)],
                                                                "size" = lengths[x,])
                  
                  ###################################
                  # make the matrices
                  matListTrimmed[[j]] = matrix(0,nrow = rnaSize, ncol = rnaSize)
                  for(i in 1:nrow(clusterPositionsListTrimmed[[j]])){
                      matListTrimmed[[j]][clusterPositionsListTrimmed[[j]][i,"ls"]:clusterPositionsListTrimmed[[j]][i,"le"],
                                          clusterPositionsListTrimmed[[j]][i,"rs"]:clusterPositionsListTrimmed[[j]][i,"re"]] =     matListTrimmed[[j]][clusterPositionsListTrimmed[[j]][i,"ls"]:clusterPositionsListTrimmed[[j]][i,"le"],
                                                                                                                                                       clusterPositionsListTrimmed[[j]][i,"rs"]:clusterPositionsListTrimmed[[j]][i,"re"]] + clusterPositionsListTrimmed[[j]][i, "size.x"]
                  }
              }
              
              
              ###################################
              # And save
              print("saving")
              print("saving mat list")
              ml[[rna]][["trimmedClusters"]] =     matListTrimmed
              
              print("saving granges list")
              cgr[[rna]][["trimmedClusters"]]   =  allChimerasForSuperClustersPlottingTrimmed
              
              print("saving table list")
              ctl[[rna]][["trimmedClusters"]] =   clusterPositionsListTrimmed
              
              
              ###################################
              # Re-make the object
              object  = new("comradesClusteredDataSet",
                            rnas = rnas(clusteredCds),
                            rnaSize = rnaSize(clusteredCds),
                            sampleTable = sampleTable(clusteredCds),
                            hybFiles = hybFiles(clusteredCds),
                            matrixList = ml,
                            group = group(clusteredCds),
                            sampleNames = sampleNames(clusteredCds),
                            clusterTableList = ctl,
                            clusterGrangesList = cgr
              )
              return(object)
          })








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







#' compareKnown
#'
#' This method compares the current object to a know structure.run 
#' \code{trimClusters()} on the  \code{comradesClusteredDataSet} first
#'
#' @param trimmedClusters a \code{comradesClusteredDataSet} object, 
#' run \code{trimClusters()} on the  \code{comradesClusteredDataSet} first
#' 
#' @param knownMat Matrix - A marix(ncol = lengthRNA,nrow = lengthRNA) where a
#' value in matrix[x,y] would indicate a known interation between nucleotide 
#' x and nucleotide y 
#' 
#' @slot rna string - a single RNA to analyse - must be present in \code{rnas(cdsObject)} 
#' 
#' @param type string - the type of clusters you would like to compare you can find 
#' available types by just running the objects name
#' 
#' 
#' @return Returns a \code{comradesClusteredDataSet} object
#' 
#' The 3 attributes matrixList, clusterTableList and clusterGrangesList 
#' will gain the \code{types} "known" and "novel" and "knownAndNovel"
#' 
#' @export
setGeneric("compareKnown", 
           function(trimmedClusters, knownMat,rna,  type, ...) standardGeneric("compareKnown"))

setMethod("compareKnown", "comradesClusteredDataSet", function(trimmedClusters, knownMat,rna, type)  {
    
    ###################################
    # Inputs
    k18Smat = knownMat
    type = type
    trimmedClusters = trimmedClusters
    sampleNames = sampleNames(trimmedClusters)
    ml = matrixList(trimmedClusters)
    rnaSize = ncol(matrixList(trimmedClusters)[[rna]][["noHost"]][[1]])
    
    ###################################
    # set up variables
    ml[[rna]][["KnownAndNovel"]] = list()
    novelClusters = list()
    novelClustersMat = list()
    novelClustersMat2 = list()
    cannonicalClusters = list()
    cannonicalClustersMat = list()
    cannonicalClustersMat2 = list()
    
    
    clusterPositionsListTrimmed = clusterTableList(trimmedClusters)[[rna]][[type]]
    # for each sample
    for(i in 1:length(clusterPositionsListTrimmed)){
        
        # set up the matrix for this sample and the cluster positions
        clusters = clusterPositionsListTrimmed[[i]]
        
        ###################################
        # test each cluster against the known interactions
        #this a matrix of known positions:
        k18Smat2 = k18Smat
        tf = c()
        for(j in 1:nrow(clusters)){
            #For each cluster, make a individual matrix
            clusterMat = matrix(0, nrow = rnaSize, ncol = rnaSize)
            # add 10 to the positions of this cluster
            clusterMat[ clusters$ls[j]:clusters$le[j] ,  clusters$rs[j]:clusters$re[j] ] =
                clusterMat[  clusters$ls[j]:clusters$le[j] ,  clusters$rs[j]:clusters$re[j] ] + 10
            
            # Add that to the matric of cannonical interactions
            k18Smat2 = k18Smat + clusterMat
            # find those when the two annotations overlap.
            tf = c(tf, all(k18Smat2 < 11))
        }
        
        ###################################
        # Get novel and cannonical tables and matrices
        print(which((tf == F)))
        novelClusters[[i]] = clusterPositionsListTrimmed[[i]][which((tf == T)),]
        novelClustersMat[[i]] = matrix(0, nrow = rnaSize, ncol = rnaSize)
        if(nrow(novelClusters[[i]]) > 0){
            for(j in 1:nrow(novelClusters[[i]])){
                
                novelClustersMat[[i]][novelClusters[[i]]$ls[j]:novelClusters[[i]]$le[j] ,
                                      novelClusters[[i]]$rs[j]:novelClusters[[i]]$re[j] ] =
                    novelClustersMat[[i]][ novelClusters[[i]]$ls[j]:novelClusters[[i]]$le[j] ,
                                           novelClusters[[i]]$rs[j]:novelClusters[[i]]$re[j] ] +  novelClusters[[i]]$size.x[j]
                
            }
        }
        
        print(which((tf == T)))
        cannonicalClusters[[i]] = clusterPositionsListTrimmed[[i]][which((tf == F)),]
        cannonicalClustersMat[[i]] = matrix(0, nrow = rnaSize, ncol = rnaSize)
        for(j in 1:nrow(cannonicalClusters[[i]])){
            
            cannonicalClustersMat[[i]][cannonicalClusters[[i]]$ls[j]:cannonicalClusters[[i]]$le[j] ,
                                       cannonicalClusters[[i]]$rs[j]:cannonicalClusters[[i]]$re[j] ] =
                cannonicalClustersMat[[i]][ cannonicalClusters[[i]]$ls[j]:cannonicalClusters[[i]]$le[j] ,
                                            cannonicalClusters[[i]]$rs[j]:cannonicalClusters[[i]]$re[j] ] + cannonicalClusters[[i]]$size.x[j]
        }
        
        ###################################
        # Add the known interactions to the known and novel matrices
        cannonicalClustersMat2[[i]] = cannonicalClustersMat[[i]] + knownMat*30000
        novelClustersMat2[[i]] = novelClustersMat[[i]] + knownMat*30000
        ml[[rna]][["KnownAndNovel"]][[i]] = ml[[rna]][[type]][[i]] + knownMat*30000
    }
    
    
    
    ###################################
    # add to the lists for the object
    ml[[rna]][["novel"]] = novelClustersMat2
    ml[[rna]][["known"]] = cannonicalClustersMat2
    ctl = clusterTableList(trimmedClusters)
    
    ctl[[rna]][["novel"]] = novelClusters
    ctl[[rna]][["known"]] = cannonicalClusters
    cgl = clusterGrangesList(trimmedClusters)
    
    
    
    ###################################
    # create object
    object  = new("comradesClusteredDataSet",
                  rnas = rnas(clusteredCds),
                  rnaSize = rnaSize(clusteredCds),
                  sampleTable = sampleTable(clusteredCds),
                  hybFiles = hybFiles(clusteredCds),
                  matrixList = ml,
                  group = group(clusteredCds),
                  sampleNames = sampleNames(clusteredCds),
                  clusterTableList = ctl,
                  clusterGrangesList = cgl
    )
    
    return(object)
    
})







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








#' Subset a list of hyb files
#'
#' Function used to subset a list of hyb data created by readHybFiles
#' This function produces the same size list as before but 
#' it returns ONLY the rna of interest and also
#' Choose duplexes where the nt difference in position between the
#' one side and other side of an interaction is between min and max
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

