#' @include  comradesDataSet.R comradesFoldedDataSet.R comradesClusteredDataSet.R
NULL


#' comradesFoldedDataSet
#'
#' Desc
#'
#' @rdname comradesFoldedDataSet
#' @export
setClass("comradesFoldedDataSet",
         contains = "comradesClusteredDataSet",
         slots = c(
             clusterTableFolded = "data.frame",
             interactionTable = "data.frame",
             viennaStructures = "character"
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
#' @param cdsObject comradesClusteredDataSet object created with comradesClusteredDataSet
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
#' @slot rna string - a single RNA to analyse - must be present in \code{rnas(cdsObject)} 
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


    
    
    ############################################################################
    # Fold the whole molecule
    ############################################################################
    
    # get interactions
    interactionTable = tableAll
    
    
    
    # Just get the columns needed
    interactionTable = interactionTable[,c(1,2,3,4,9,10,11,12)]
    
    
    # Aggregate the table to combine interactions by sample 
    interactionTable2 = aggregate(interactionTable$evidence, by = list(interactionTable$p1,  interactionTable$p2, interactionTable$nt1 ,
                                                                       interactionTable$nt2, interactionTable$sample), FUN = sum)
    colnames(interactionTable2) = c("p1","p2","nt1","nt2","sample","evidence")
    interactionTable3 = aggregate(interactionTable2$sample, by = list(interactionTable2$p1,  interactionTable2$p2, interactionTable2$nt1 ,
                                                                      interactionTable2$nt2), FUN = paste, collapse=",")
    interactionTable4 = aggregate(interactionTable2$evidence, by = list(interactionTable2$p1,  interactionTable2$p2, interactionTable2$nt1 ,
                                                                        interactionTable2$nt2), FUN = paste, collapse=",")
    
    interactionTable5 = aggregate(interactionTable2$evidence, by = list(interactionTable2$p1,  interactionTable2$p2, interactionTable2$nt1 ,
                                                                        interactionTable2$nt2), FUN = sum)
    
    # Check if the interactions agree with canonical sites
    
    #interactionTable3$cannonical = 0
    
#    for(i in 1:nrow(interactionTable3)){
#        if(paste(interactionTable3$Group.1[i], interactionTable3$Group.2[i]) %in%
#           paste(known18S$V1, known18S$V2)){
#            interactionTable3$cannonical[i] = 1
#    }
#        }
    
#    unique(interactionTable3$cannonical)
#    nrow(known18S)
    
    interactionTable3$evidence = interactionTable5$x
    interactionTable3$evidence2 = interactionTable4$x
    
    colnames(interactionTable3) = c("p1", "p2", "nt1", "nt2", "samples" ,"evidence", "evidence2")
    
    
    # SUBSET the interaction based on thei evidence and appearing in all the samples 
    
    interactionTable3_sub = interactionTable3[interactionTable3$evidence > 1000 &
                                                  interactionTable3$samples == "s1,s2,s3"   , ]
    
    unique(interactionTable3_sub$samples)

    
    
    
    head(interactionTable3_sub)
    nrow(interactionTable3_sub)
    
    # find the probability of each constraint
    normalized_evidence = interactionTable3_sub$evidence / sum(interactionTable3_sub$evidence)
    

    
    # this is how you sample the "bag"
    #sample(1:nrow(interactionTable3_sub),1,prob = normalized_evidence)
    
    # run the structures 
    prevConstraints = 0
    
    viennas = c()
    
    for(j in 1:100){
        goodvienna = ""
        for(i in 1:10) {
            # pull constraints and re pick if constraint breaks the structure
            constraints = c(prevConstraints, sample(1:nrow(interactionTable3_sub),1,prob = normalized_evidence) )
            prevConstraints = constraints
            #write contraint to file
            
            constraints = unique(constraints)
            
            constraintFile = interactionTable3_sub[constraints, c("p1","p2")]
            #F i j k
            constraintFile$F = "F"
            constraintFile$k = 1
            constraintFile = constraintFile[,c(3,1,2,4)]
            print(constraintFile)
            write.table(constraintFile, file = "constraints.txt", quote = F, row.names = F, col.names = F)
            
            
            
            
            
            table = data.frame(stringsAsFactors = FALSE)
            command = paste("echo \">",rna,"\n",paste(rnaRefs[[1]]$NR_003286.4, collapse = ""),"\" | RNAfold  --noPS --constraint=constraints.txt ", sep = "")
            x = system(command,intern = T)
            
            if(!(grepl("\\(\\(", x[3])) ){
                prevConstraints = prevConstraints[-1]
                next;
            }else{
                goodvienna = x[3]
            }
            print(prevConstraints)
            print(x)
            
            
            
        }
        viennas = c(viennas,sub("\\s.*","",goodvienna))
        prevConstraints = 0
    }
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    ###########################################################
    # Make object
    ###########################################################
    print(" *** Creating object ***")
    #create comrades dataset object
    object  = new("comradesFoldedDataSet",
                  rnas = rnas(cdsObject),
                  rnaSize = rnaSize(cdsObject),
                  sampleTable = sampleTable(cdsObject),
                  hybFiles = hybFiles(cdsObject),
                  matrixList = matrixList(cdsObject),
                  group = group(cdsObject),
                  sampleNames = sampleNames(cdsObject),
                  clusterTableList = clusterTableList(cdsObject),
                  clusterGrangesList = clusterGrangesList(cdsObject),
                  clusterTableFolded = clusterPositionsListTrimmedSarsCombinedWithStructures,
                  interactionTable = tableAll,
                  viennaStructures = viennas
    )

    return(object)

}

