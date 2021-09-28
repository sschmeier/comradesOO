#' @include comradesOO.R comradesDataSet.R






#' comradesClusteredDataSet
#'
#' This object is a parent of comradesDataSet with 2 extra attributes
#' (clusterTableList and clusterGrangesList) it takes as input a 
#' comradesDataSet object and performs clustering of the duplexes. 
#'
#'
#'
#' @rdname comradesClusteredDataSet
#' @export
#' 
setClass("comradesClusteredDataSet",
         contains = "comradesDataSet",
         slots = c(
             clusterGrangesList = "list",
             clusterTableList = "list"
         ),
         prototype = list(

         ))

setValidity("comradesClusteredDataSet", function(object) {

    #test to make sure the object is O.K

})




#' comradesClusteredDataSet
#'
#' This object is a parent of comradesDataSet with 2 extra attributes
#' (clusterTableList and clusterGrangesList) it takes as input a 
#' comradesDataSet object and performs clustering of the duplexes. 
#'
#'
#' @param cds comradesDataSet object created with comradesDataSet
#' @param cores numeric - The number of cores to use 
#'
#'
#' @slot clusterTableList List - Follows the pattern for list slots of comradesDataSet
#' objects, \code{matrixList(cds)[[rna]][[type]][[sample]]}. contains a table
#' with coordinates and information about the clusters identified
#' @slot clusterGrangesList List - Follows the pattern for list slots of comradesDataSet
#' objects, \code{matrixList(cds)[[rna]][[type]][[sample]]}. contains GRanges 
#' objects of the original suplexes with their cluster membership
#' #' @slot sampleTable table - Column names; fileName, group (s or c),
#'  sample (1,2,3, etc), sampleName (must be unique)
#' @slot rnas vector - A vector of RNA names to analyse 
#' @slot group list - This is made from the a sample table during object 
#' creation, it is a list with two vector elements ("c","s") containing the 
#' indexes of the sampleTable that have "c" or "s" in the group column.
#' @slot matrixList List - Follows the pattern for list slots of comradesDataSet
#' objects, \code{matrixList(cds)[[rna]][[type]][[sample]]}. Contains a set
#' of contact matrices, each cell contains the number of duplexes identified 
#' for position x,y.
#' @slot hybFiles List - Follows the pattern for list slots of comradesDataSet
#' objects, \code{hybFiles(cds)[[rna]][[type]][[sample]]}. Contains a set of 
#' tables, these are the original Hyb files that were read in. 
#'
#' @rdname comradesClusteredDataSet
#' @export

comradesClusteredDataSet <- function(cds,
                                     cores) {


    ###########################################################
    # Create clusters from the cds object
    ###########################################################


    ##############################
    # Set up variables
    clusters = list()
    clusterTables = list()
    for(rna in rnas(cds)){


        ##############################
        # Set up variables
        clusters[[rna]] = list() #will contain the clusterGranges
        clusterTables[[rna]] = list() #will contain the cluster tables
        rnaSize = ncol(matrixList(cds)[[rna]][["noHost"]][[1]]) #calculate size of rna
        print(paste("*** clustering ", rnaSize," nt ", rna, " ***"))



        ##############################
        # get long and short range interactions

        ##############################
        # short Long interactions
        print(paste("*** doing Long Range***"))

        #subset hybLists2
        longDistHyb = subsetHybList2(hybFiles(cds)[[rna]][[ "noHost" ]],
                                     10,rnaSize,length = 800)

        chimeraList = hybToGRanges(longDistHyb,rna)

        print(paste("*** sampling Long Range***"))
        # Reduce the GRanges to make them smaller
        chimeraListSampled = sampleChimeras(chimeraList)

        matrixList = matrixList(cds)
        print(paste("*** clustering Long Range***"))
        registerDoParallel(cores)
        plottingList = list()

        for(i in 1:length(sampleNames(cds))){

            plottingList[[i]] = list()
            foreach (k=1:length( chimeraListSampled[[i]][["gap"]])) %do% {
                adjacancyMat = getAdjacancyMat(chimeraListSampled[[i]][["gap"]][[k]],"nucleotide", 15)
                net = graph_from_adjacency_matrix(adjacancyMat,
                                                  mode = "undirected",
                                                  weighted = T)
                clustering = cluster_walktrap(net,steps = 2)

                highest_clusters = names(table(membership(clustering))[table(membership(clustering)) > 10])
                plottingList[[i]][[k]]  = printClustersFast(dir,clustering, highest_clusters,
                                                            chimeraListSampled[[i]][["left"]][[k]],
                                                            chimeraListSampled[[i]][["right"]][[k]])

                print(k)
            }

        }


        longRange = plottingList

        ##############################
        # short Long interactions
        print(paste("*** doing Short Range***"))
        longDistHyb = subsetHybList2(hybFiles(cds)[[rna]][[ "noHost" ]],
                                     1,9,length = 800)

        chimeraList = hybToGRanges(longDistHyb,rna)
        # Reduce the GRanges to make them smaller
        print(paste("*** sampling short Range***"))
        chimeraListSampled = sampleChimeras(chimeraList)

        print(paste("*** clustering short Range***"))
        plottingList = list()
        for(i in 1:length(sampleNames(cds))){
            plottingList[[i]] = list()
            foreach (k=1:length( chimeraListSampled[[i]][["gap"]])) %do% {
                adjacancyMat = getAdjacancyMat(chimeraListSampled[[i]][["gap"]][[k]],"nucleotide", 15)
                net = graph_from_adjacency_matrix(adjacancyMat, mode = "undirected", weighted = T)
                clustering = cluster_walktrap(net,steps = 2)
                highest_clusters = names(table(membership(clustering))[table(membership(clustering)) > 10])
                plottingList[[i]][[k]]  = printClustersFast(dir,clustering, highest_clusters,
                                                            chimeraListSampled[[i]][["left"]][[k]],
                                                            chimeraListSampled[[i]][["right"]][[k]])


            }

        }


        shortRange = plottingList



        ##############################
        #make the IDs for the table and combine
        combinedPlotting = shortRange
        for(i in 1:length(shortRange)){
            plotting = GRangesList(shortRange[[i]])
            for(j in 1:length(plotting)){
                plotting[[j]]$k = paste(plotting[[j]]$cluster,"binShort",j, sep = ":")
            }
            combinedPlotting[[i]] =  plotting
        }

        longRange
        for(i in 1:length(longRange)){
            plotting = GRangesList(longRange[[i]])
            for(j in 1:length(plotting)){
                plotting[[j]]$k = paste(plotting[[j]]$cluster,"binLong",j, sep = ":")
            }
            combinedPlotting[[i]] =  c(combinedPlotting[[i]],plotting)
        }



        clusters[[rna]][["original"]] = combinedPlotting




        ##############################
        # Make a matrix and table for the clusters
        clusterPositionsList = list()
        matList = list()
        for(j in 1:length(sampleNames(cds))){


            plotting =  unlist(combinedPlotting[[j]])

            lengths = aggregate(mcols(plotting)$cluster, by = list(mcols(plotting)$k), FUN = length)
            row.names(lengths) = lengths$Group.1
            #for each cluster get the min start and max end
            plottingSplit = split(plotting, paste(mcols(plotting)$k, mcols(plotting)$type))
            minStarts = unlist(lapply(plottingSplit, function(x) {return(min(start(x)))  }))
            maxEnd = unlist(lapply(plottingSplit, function(x) {return(max(end(x)))  }))
            clusterPositionsList[[j]] = data.frame("id" = names(maxEnd)[seq(1,length(minStarts),2)],
                                                   "ls" = minStarts[seq(1,length(minStarts),2)],
                                                   "le" = maxEnd[seq(1,length(maxEnd),2)],
                                                   "rs" = minStarts[seq(2,length(minStarts),2)],
                                                   "re" = maxEnd[seq(2,length(maxEnd),2)],
                                                   "size" = lengths[sub("\\s.*","",names(maxEnd)[seq(1,length(minStarts),2)]),])



            matList[[j]] = matrix(0,nrow = rnaSize, ncol = rnaSize)


            for(i in 1:nrow(clusterPositionsList[[j]])){
                matList[[j]][clusterPositionsList[[j]][i,"ls"]:clusterPositionsList[[j]][i,"le"],
                             clusterPositionsList[[j]][i,"rs"]:clusterPositionsList[[j]][i,"re"]] =    matList[[j]][clusterPositionsList[[j]][i,"ls"]:clusterPositionsList[[j]][i,"le"],
                                                                                                                    clusterPositionsList[[j]][i,"rs"]:clusterPositionsList[[j]][i,"re"]] + clusterPositionsList[[j]][i, "size.x"]

            }
        }
        names(matList) = sampleNames(cds)
        names(clusterPositionsList) = sampleNames(cds)

        clusterTables[[rna]][["original"]] =clusterPositionsList
        matrixList[[rna]][["originalClusters"]] = matList
        clusterGranges = clusters
    }



    ###########################################################
    # Make object
    ###########################################################
    print(" *** Creating object ***")
    #create comrades dataset object
    object  = new("comradesClusteredDataSet",
                  rnas = rnas(cds),
                  rnaSize = rnaSize(cds),
                  sampleTable = sampleTable(cds),
                  hybFiles = hybFiles(cds),
                  matrixList = matrixList,
                  group = group(cds),
                  sampleNames = sampleNames(cds),
                  clusterTableList = clusterTables,
                  clusterGrangesList = clusterGranges
    )

    return(object)

}
