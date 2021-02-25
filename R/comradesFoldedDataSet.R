#' @include comradesOO.R comradesDataSet.R comradesClusteredDataSet.R



#' comradesFoldedDataSet
#'
#' Desc
#'
#' @inheritParams comradesClusteredDataSet
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
#' Desc
#'
#' @inheritParams comradesClusteredDataSet
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

