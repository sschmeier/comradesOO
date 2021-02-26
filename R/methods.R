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











##################################################
###########  Other methods         ###############
##################################################






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
                          removal =  mean(x) + sd(x)*1
                          included  = GRanges(seqnames=rna,
                                              IRanges(
                                                  start=rep(min, length(clusterrange)),
                                                  end=rep(removal, length(clusterrange))
                                              ))
                          if( l == "left"){
                              removal =  mean(x) -  sd(x)*1
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
                            hybDir = hybDir(clusteredCds),
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
    type = "trimmedClusters"
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
        for(j in 1:nrow(novelClusters[[i]])){

            novelClustersMat[[i]][novelClusters[[i]]$ls[j]:novelClusters[[i]]$le[j] ,
                                  novelClusters[[i]]$rs[j]:novelClusters[[i]]$re[j] ] =
                novelClustersMat[[i]][ novelClusters[[i]]$ls[j]:novelClusters[[i]]$le[j] ,
                                       novelClusters[[i]]$rs[j]:novelClusters[[i]]$re[j] ] +  novelClusters[[i]]$size.x[j]

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
                  hybDir = hybDir(clusteredCds),
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







# Make stats from cds object
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
# get stats

# Make stats from cds object
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
setGeneric("plotMatrices", function(cds,type, directory,a,b,c,d,h, ...) standardGeneric("plotMatrices"))

setMethod("plotMatrices", "comradesDataSet", function(cds,type, directory,a,b,c,d,h)  {

    hybMatList = matrixList(cds)
    rnaS = rnas(cds)
    sampleNames = sampleNames(cds)


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
setGeneric("plotMatricesAverage", function(cds,type, directory,a,b,c,d,h, ...) standardGeneric("plotMatricesAverage"))

setMethod("plotMatricesAverage", "comradesDataSet", function(cds,type, directory,a,b,c,d,h)  {

    for(rna in rnas(cds)) {

        hybMatList = matrixList(cds)
        hybMatList2 = hybMatList
        c = 1
        for(i in "s"){
            for(j in group(cds)[[i]] ){
                if(length( group(cds)[[i]] ) < 2 | c == 1 ){
                    sum(hybMatList[[rna]][[type]][[ j ]])
                    hybMatList2[[rna]][[type]][[ i ]] =   hybMatList[[rna]][[type]][[ j ]]
                    sum(hybMatList2[[rna]][[type]][[ i ]])
                    print("first")
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

        sampleNames = c("s")
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







##################################################
########### Helper Functions       ###############
##################################################



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






#' Make a matrix of contact interactions
#'
#' Function used to create a list of matrices for plotting with
#' plotMatrixList or plotMatrixListFull, the output list will be same as the
#' input except for an extra list layer for the specific RNA
#' @param hybList the original hybList created with readHybFiles or subsetHybList
#' @param rna the rna of interest that you want to subset
#' @return A list of matrices
#' @examples
#' hybMatList <- subsetHybList(hybList, "myRNA");
#' @export



getMatrices = function(hybList, rna, size){
    hybMatList = list()
    for(hyb in 1:length(hybList)){

        hybOutputO = hybList[[hyb]]



        hybOutput =  hybOutputO[as.character(hybOutputO$V4) == rna & as.character(hybOutputO$V10) == rna,]
        hybOutput = unique(hybOutput)

        startsends = hybOutput[,c(7,8,13,14)]



        mat = matrix(0, ncol = size, nrow =size)
        for(i in 1:nrow(startsends)){
            data = startsends[i,]
            xData = seq(data$V7, data$V8)
            yData = seq(data$V14, data$V13)
            mat[xData,yData] = mat[xData,yData] +1


        }

        hybMatList[[hyb]] = mat

    }
    return(hybMatList)
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











#' Makes and adjacency matrix list (for clustering)
#'
#' Makes and adjacency matrix list (for clustering)
#' @param hybGranges list created with hybToGRanges (but just the gap section of the list)
#' @param  nucletideOrPerc measure difference by percentage or nucleotides
#' @param  cutoff The maximum difference before giving these two gaps 0
#' @return A list of Adjacancy matrices
#' @examples
#'  adjcancyMat =  getAdjacancyMat(hybGranges[["sample"]][["gap"]],"nucleotide", 5);
#' @export



getAdjacancyMat = function(hybGranges, nucletideOrPerc, cutoff){
    distances = hybGranges
    max = max(width(distances))

    #get overlapping
    hits <- findOverlaps(distances, drop.self=T, drop.redundant=F)
    # get the relative overlap for the weights
    x <- distances[queryHits(hits)]
    y <- distances[subjectHits(hits)]

    relative_overlap <-  width(pintersect(x, y)) / pmax(width(x), width(y))

    hitsWithOverlap = hits
    # parameter for relative overlap

    #print(length(hitsWithOverlap))

    if(nucletideOrPerc == "none"){
        hitsWithOverlap = hits
    } else if(nucletideOrPerc == "nucleotide"){
        relative_overlap =     pmax(width(x), width(y)) - width(pintersect(x, y))
        relative_overlap = cutoff - relative_overlap


        hitsWithOverlap = hits[relative_overlap <= cutoff & relative_overlap >= 0 ]
        relative_overlap = relative_overlap[relative_overlap <= cutoff  & relative_overlap >= 0 ]
    } else if(nucletideOrPerc == "perc"){

        relative_overlap = (1- (width(pintersect(x, y)) / max))

        hitsWithOverlap = hits[relative_overlap >= cutoff]
        if(length(hitsWithOverlap ) ==0 ){return(0)}else{}
        # print(length(hitsWithOverlap))
        relative_overlap = relative_overlap[relative_overlap >= cutoff]

    }

    hitsMat = as.data.frame(hitsWithOverlap)
    hitsMat$weight = relative_overlap


    testLong = dcast(hitsMat, queryHits ~ subjectHits, value.var = "weight")
    row.names(testLong) = testLong$queryHits
    testLong = testLong[,-1]
    testLong = as.matrix(testLong)
    testLong[is.na(testLong)] =0

    return(testLong)

}




