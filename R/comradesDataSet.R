#' @include comradesOO.R


#' comradesDataSet
#'
#' An S4 class to represent a COMRADES dataset
#'
#'
#'
#'
#' @slot sampleTable Table File Name -  Column names - fileName, group (s or c),
#'  sample (1,2,3, etc), sampleName (must be unique)
#' @slot rnas Vector - A vecotor of RNA names to analyse. (must be in the Hyb
#' input files)
#' @slot group This is made from the a sample table
#' @slot matrixList List - Follows the pattern for list slots of comradesDataSet
#' objects, \code{matrixList(cds)[[rna]][[type]][[sample]]}. Contains a set
#' of contact matrices, each cell contains the number of duplexes identified 
#' for position x,y.
#' @slot hybFiles List - Follows the pattern for list slots of comradesDataSet
#' objects, \code{hybFiles(cds)[[rna]][[type]][[sample]]}. Contains a set of 
#' tables, these are the original Hyb files that were read in. 
#'
#'
#' @rdname comradesDataSet
#' @export
#'
setClass("comradesDataSet",
         slots = c(
             rnas = "character",
             sampleTable = "data.frame",    # meta data
             hybFiles = "list",       # data tables
             matrixList = "list",
             group = "list",
             sampleNames = "character" # Column containing control or sample
         ),
         prototype = list(

         ))

setValidity("comradesDataSet", function(object) {

    #test to make sure the object is O.K

})




#' comradesDataSet
#'
#' \code{comradesDataSet} objects are used to store the input meta-data, data and
#' create a framework for the storage of results. Whilst creating the object, 
#' the original hyb files are also filtered for the RNA of interest. 
#'
#'
#' @param rnas vector - The names of the RNA interest, these must be displayed
#' the same way as in the input Hyb Files. 
#' @param rnaSize named list - The sizes (nt) of the RNAs of interest, the list
#'  elements must have same names as the \code{rnas} vector and each each contain 
#'  one numeric value. 
#' @param sampletable string - The address of the sample table, the sample table
#'  must have 4 columns, fileName (the full path and file name of the input 
#'  hyb file for each sample ), group ("s" - sample or "c" - control), 
#'  sample (1,2,3, etc), sampleName (must be unique).
#'  
#' @return A comradesDataSet object.
#' 
#' 
#' 
#' @slot sampleTable table - Column names; fileName, group (s or c),
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
#' @references See \url{https://github.com/gkudla/hyb} for hyb
#'
#' @docType class
#' 
#' @rdname comradesDataSet
#'
#' @examples
#'
#'
#' @export

comradesDataSet <- function(rnas,
                            rnaSize,
                            sampleTable) {
    ###########################################################
    # Read in the sample table
    ###########################################################
    # check the inputs here, stop if wrong
    print(" ***** ******************* ****** ")
    print(" ***** Reading SampleTable ****** ")

    # Read in sample table
    sampleTable = read.table(sampleTable,
                             header = T, stringsAsFactors = F)
    #check for more than two samples
    if( nrow(sampleTable) < 2 ){
        stop( "The sample Table must contain at least 1
              sample and 1 control" )
    }
    print(paste("***  detected ",nrow(sampleTable), " samples  ***"))


    #check column names of sampleTable
    colnamesST = c("file", "group", "sample", "sampleName")
    if(all(colnames( sampleTable ) != colnamesST)){
        stop( "Column names of metaData table should be :
              file, group, sample, sampleNames" )
    }

    ###########################################################
    # Get the comparison groups
    # check group has the c and s
    if(! ( unique(as.character( sampleTable$group ) )[1] %in% c("c", "s") &
           unique(as.character( sampleTable$group ) )[2] %in% c("c", "s") ) ) {
        stop( "Groups should be c and s" )
    }

    # Make group into a list with control and sample
    group = sampleTable[,"group"]
    group2 = list()
    group2[["c"]] = which(group == "c")
    group2[["s"]] = which(group == "s")

    group = group2
    print(paste("*** detected group c::",  paste(group[["c"]], collapse = " ") , "***"))
    print(paste("*** detected group s::",  paste(group[["s"]], collapse = " ") , "***"))



    ###########################################################
    # Get the sampleNames
    sampleNames = c()
    if( is.null(sampleTable$sampleName) ){
        stop( "The sample Table must have a column named sampleName" )
    }else if( length(unique( sampleTable$sampleName)) != length( sampleTable$sampleName) ){
        stop( "SampleName column must be unique" )
    }else{
        sampleNames = as.character( sampleTable$sampleName )

        print(paste("*** detected ", paste(sampleNames, collapse = " "), " sample Names ***"))
    }



    ###########################################################
    # Read in the  hyb files
    ###########################################################
    #load the files into a list
    print(" ***** Reading Hyb Files ******")

    hybFiles = list()
    hybFiles[[ "all" ]] = list()
    hybFiles[[ "all" ]][[ "all" ]] = list()

    for(i in 1:nrow( sampleTable )){

        #get file and path
        file =  as.character( sampleTable$file[i] )
        print( file )

        # Read in
        sampleHyb = read.table( file ,
                                header = F,
                                stringsAsFactors = F )

        #check the hyb file column names
        colnamesHyb = c("V1", "V2", "V3", "V4", "V5",
                        "V6", "V7", "V8", "V9", "V10",
                        "V11", "V12", "V13", "V14", "V15")

        if( !( identical(colnames(sampleHyb), colnamesHyb) ) ){
            stop(" The input hyb files do not look they are produced with the
                 hyb program. ")
        }

        # Store
        hybFiles[[ "all" ]][[ "all" ]][[ sampleNames[i] ]] = unique( sampleHyb )
        print(nrow(hybFiles[[ "all" ]][[ "all" ]][[ sampleNames[i] ]] ))
    }


    ###########################################################
    # Change the hyb files to have specific rna of interest
    # and with host and without
    ###########################################################
    print(" ***** Getting RNAs of Interest ******")
    # Get the rna of interest it comes in two ways, with host and without
    for(i in 1:length(rnas)){
        
        print(" *** RNA of interest + Host RNA ***")
        hybFiles[[ rnas[ i ] ]][[ "original"]] = swapHybs(hybList = hybFiles[[ "all" ]][[ "all" ]],
                                                          rna = rnas[ i ] )

        hybFiles[[ rnas[ i ] ]][[ "host"]] = swapHybs3(hybList = hybFiles[[ "all" ]][[ "all" ]],
                                                       rna = rnas[ i ] )
        print(" *** RNA of interest Alone ***")
        hybFiles[[ rnas[ i ] ]][[ "noHost"]] = swapHybs2(hybList = hybFiles[[ "all" ]][[ "all" ]],
                                                         rna = rnas[ i ] )
    }


    ###########################################################
    # Make matrices of the specific RNA without host
    ###########################################################
    print(" ***** Making Matrices ******")

    matrixList = list()
    c = 1
    for(i in rnas){
        rnaSize2 =   rnaSize[[i]]
        print(i)
        matrixList[[i]][[ "noHost" ]] = list()
        matrixList[[i]][[ "noHost" ]] = getMatrices(hybFiles[[ i  ]][[ "noHost"]],
                                                    i, rnaSize2)
        names(matrixList[[i]][[ "noHost" ]]) = sampleNames

        matrixList[[i]][[ "original" ]] = list()
        matrixList[[i]][[ "original" ]] = getMatrices(hybFiles[[ i  ]][[ "original"]],
                                                      i, rnaSize2)
        names(matrixList[[i]][[ "original" ]]) = sampleNames
        c = c +1
    }



    ###########################################################
    # Make object
    ###########################################################
    print(" *** Creating object ***")
    #create comrades dataset object
    object  = new("comradesDataSet",
                  rnas = rnas,
                  sampleTable = sampleTable,
                  hybFiles = hybFiles,
                  matrixList = matrixList,
                  group = group,
                  sampleNames = sampleNames)

    return(object)


}



