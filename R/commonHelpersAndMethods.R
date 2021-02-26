#' @include  comradesDataSet.R comradesFoldedDataSet.R comradesClusteredDataSet.R
NULL


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








##################################################
########### Helper Functions       ###############
##################################################



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




