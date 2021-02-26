#' @include comradesOO.R comradesDataSet.R comradesClusteredDataSet.R



#' comradesFoldedDataSet
#'
#' Desc
#'
#' @rdname comradesFoldedDataSet
#' @export
setClass("comradesFoldedDataSet",
         contains = "comradesClusteredDataSet",
         slots = c(
             clusterTableFolded = "data.frame"
         ),
         prototype = list(

         ))

setValidity("comradesFoldedDataSet", function(object) {

    #test to make sure the object is O.K

})





#' comradesFoldedDataSet
#'
#'
#' This object is a parent of comradesClusteredDataSet with one extra attribute
#' clusterTableFolded.
#'
#' @param cdsObject comradesDataSet object created with comradesDataSet
#' @param rna string - a single RNA to analyse  
#' @param rnarefs named List - a list with named elements that correspond to the 
#'     .RNA of interest. The element of the list must be a fasta file that has 
#'     been read with \code{seqinr::read.fasta()}
#'
#' @slot clusterTableFolded table - a table similar to the \code{clusterTableList}
#' it contains coordinates of the clusters along with vienna format fold and 
#' RNA sequences for each cluster
#' @slot clusterTableList List - Follows the pattern for list slots of comradesDataSet
#' objects, \code{matrixList(cds)[[rna]][[type]][[sample]]}. contains a table
#' with coordinates and information about the clusters identified
#' @slot clusterGrangesList List - Follows the pattern for list slots of comradesDataSet
#' objects, \code{matrixList(cds)[[rna]][[type]][[sample]]}. contains GRanges 
#' objects of the original duplexes with their cluster membership
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
#' @rdname comradesFoldedDataSet
#' @export

comradesFoldedDataSet <- function(cdsObject,
                                  rna,
                                  rnaRefs) {


    ###########################################################
    # Fold each cluster and add vienna to the table
    ###########################################################


    ##############################
    # get trimmed cluster tables

    clusterPositionsListTrimmed = clusterTableList(cdsObject)[[rna]][["trimmedClusters"]]


    ##############################
    #make combined tables for the samples
    clusterPositionsListTrimmedSarsCombined = list()
    for(j in group(cdsObject)[["s"]]){

        #  clusterPositionsListTrimmedSarsCombined[[j]]$sample =
        clusterPositionsListTrimmed[[j]]$sample = sampleTable(cdsObject)[j,"sampleName"]
        clusterPositionsListTrimmedSarsCombined = rbind.data.frame(clusterPositionsListTrimmedSarsCombined,
                                                                   clusterPositionsListTrimmed[[j]],
                                                                   stringsAsFactors = F)
    }
    ##############################
    # add the sequences to the table
    # seq1 = left seq2 = right, type = short or long (range interaction)
    colnames(clusterPositionsListTrimmedSarsCombined )
    clusterPositionsListTrimmedSarsCombined$type = ""
    clusterPositionsListTrimmedSarsCombined$seq1 = ""
    clusterPositionsListTrimmedSarsCombined$seq2 = ""
    # add the sequences tot he table
    for(i in 1:nrow(clusterPositionsListTrimmedSarsCombined)){
        x = getClusterClusterShortRangeWhole(clusterPositionsListTrimmedSarsCombined[i,],rnaRefs[[rna]] )
        clusterPositionsListTrimmedSarsCombined$type[i]= x[[1]]
        clusterPositionsListTrimmedSarsCombined$seq2[i] = x[[3]]
        clusterPositionsListTrimmedSarsCombined$seq1[i] = x[[2]]
        print(i)
    }


    ##############################
    # now Fold

    tableAll = data.frame()


    for(i in 1:nrow(clusterPositionsListTrimmedSarsCombined)){
        table = c()

        if(clusterPositionsListTrimmedSarsCombined[i,"type"] == "short"){

            table = findBasePairsRNAfold(startPos = clusterPositionsListTrimmedSarsCombined[i,"ls"],
                                         endPos = clusterPositionsListTrimmedSarsCombined[i,"re"],
                                         seqs = clusterPositionsListTrimmedSarsCombined[i,"seq1"],
                                         fasta = rnaRefs)
        }else{
            table = findBasePairsRNAcoFold(startPos1 = clusterPositionsListTrimmedSarsCombined[i,"ls"],
                                           endPos1 = clusterPositionsListTrimmedSarsCombined[i,"le"],
                                           seq1 = clusterPositionsListTrimmedSarsCombined[i,"seq1"],
                                           startPos2 = clusterPositionsListTrimmedSarsCombined[i,"rs"],
                                           endPos2 = clusterPositionsListTrimmedSarsCombined[i,"re"],
                                           seq2 = clusterPositionsListTrimmedSarsCombined[i,"seq2"],
                                           fasta = rnaRefs)
        }
        print(i)
        table$cluster = clusterPositionsListTrimmedSarsCombined[i,"id"]
        table$evidence = clusterPositionsListTrimmedSarsCombined[i,"size.x"]
        table$sample = clusterPositionsListTrimmedSarsCombined[i,"sample"]
        table$type = clusterPositionsListTrimmedSarsCombined[i,"type"]
        clusterPositionsListTrimmedSarsCombined$vienna[i] = table$Group.5[1]
        clusterPositionsListTrimmedSarsCombined$seq1new[i] = table$Group.6[1]
        clusterPositionsListTrimmedSarsCombined$seq2new[i] = table$Group.7[1]
        tableAll = rbind.data.frame(table,tableAll)
        print("done")
    }


    colnames(tableAll) = c("p1","p2","nt1","nt2","vienna", "seq","x1","x2","cluster","evidence","sample","type")
    clusterPositionsListTrimmedSarsCombinedWithStructures = clusterPositionsListTrimmedSarsCombined



    ###########################################################
    # Make object
    ###########################################################
    print(" *** Creating object ***")
    #create comrades dataset object
    object  = new("comradesFoldedDataSet",
                  rnas = rnas(cdsObject),
                  hybDir = hybDir(cdsObject),
                  sampleTable = sampleTable(cdsObject),
                  hybFiles = hybFiles(cdsObject),
                  matrixList = matrixList(cdsObject),
                  group = group(cdsObject),
                  sampleNames = sampleNames(cdsObject),
                  clusterTableList = clusterTableList(cdsObject),
                  clusterGrangesList = clusterGrangesList(cdsObject),
                  clusterTableFolded = clusterPositionsListTrimmedSarsCombinedWithStructures
    )

    return(object)

}

