#' @include  comradesDataSet.R comradesFoldedDataSet.R comradesClusteredDataSet.R
NULL

#' findBasePairsRNAcoFold
#'
#' Folds the clusters using Vienna RNA duplex
#' 
#' @param  startPos1  Start of the cluster side x
#' @param  endPos1  End of the cluster side x
#' @param  seq1 Sequence of x
#' @param  startPos2 Start of the cluster side y
#' @param  startPos2 End of the cluster side y
#' @param  seq2 Sequence of y
#' @param  fasta \code{rnaRefs}
#' 
#' @return A table of clusters and coordinates with folds
#' 
#' @examples
#' 
#' @export
findBasePairsRNAcoFold = function(startPos1, endPos1,seq1,startPos2, endPos2,seq2, fasta){
    
    
    table = data.frame(stringsAsFactors = FALSE)
    #create vienna command and run capturing output
    command = paste("echo \">",startPos1,"-", endPos1,"\n",seq1,"\n>",startPos2,"-",endPos2,"\n",seq2,"\" | RNAduplex", sep = "")
    x = system(command,intern = T)
    
    #extract the vienna and the start locations
    vienna = sub(" .*","",x[3])
    startl = as.numeric(gsub(".* (\\d{1,3}),\\d{1,3}\\s+:\\s+(\\d{1,3}),\\d{1,3} .*", "\\1", x[3], perl = T))
    endl = as.numeric(gsub(".* \\d{1,3},(\\d{1,3})\\s+:\\s+\\d{1,3},(\\d{1,3}) .*", "\\1", x[3], perl = T))
    startr = as.numeric(gsub(".* (\\d{1,3}),\\d{1,3}\\s+:\\s+\\d{1,3},(\\d{1,3}) .*", "\\2", x[3], perl = T))
    endr = as.numeric(gsub(".* \\d{1,3},(\\d{1,3})\\s+:\\s+(\\d{1,3}),\\d{1,3} .*", "\\2", x[3], perl = T))
    #print(nchar(vienna))
    seq1New = substr(seq1,startl,endl)
    seq2New = substr(seq2,endr,startr)
    
    viennar =   startPos2 + (startr-1)
    viennal =   startPos1 + (startl-1)
    
    additionl = 0
    additionr = 0
    for(i in 1:nchar(vienna)){
        j = i+ additionl
        k = i + additionr
        l = substr(vienna,j,j)
        r = substr(vienna,(nchar(vienna)+1)-k,(nchar(vienna)+1)-k)
        #print(paste(j,k,l,r))
        
        
        if(l == "(" & r == ")"){
            posl = (j-1) + viennal
            posr =  viennar - (k-1)
            p1 = fasta[[1]]$NR_003286.4[posl]
            p2 = fasta[[1]]$NR_003286.4[posr]
            #print(p1)
            #table = rbind.data.frame(table,c(j,k,posl,posr,fasta$RN18S1[posl],fasta$RN18S1[posr]))
            table = rbind.data.frame(table,c(as.numeric(j),as.numeric(k),as.numeric(posl),as.numeric(posr),as.character(p1),as.character(p2)),stringsAsFactors = FALSE)
        }
        
        if(l == "." & r == ")"){
            additionr = additionr -1
        }
        if(l == "(" & r == "."){
            additionl = additionl -1
        }
        if(l == "&" || r == "&"){
            break
        }
    }
    #}
    #}
    table$veinna = vienna
    table$seq1New = seq1New
    table$seq2New = seq2New
    colnames(table) = c("l","r","rl","rr","pl","pr","vienna","seq1new","seq2new")
    #print(table)
    table$check = 1
    aggTable = aggregate(table$check, by=list(table$rl,table$rr,table$pl,table$pr, table$vienna, table$seq1new, table$seq2new), FUN = sum )
    return(aggTable)
}




#' findBasePairsRNAcoFold2
#'
#' Folds the clusters using Vienna RNA duplex
#' 
#' @param  startPos1  Start of the cluster side x
#' @param  endPos1  End of the cluster side x
#' @param  seq1 Sequence of x
#' @param  startPos2 Start of the cluster side y
#' @param  startPos2 End of the cluster side y
#' @param  seq2 Sequence of y
#' @param  fasta \code{rnaRefs}
#' 
#' @return A table of clusters and coordinates with folds
#' 
#' @examples
#' 
#' @export
findBasePairsRNAcoFold2 = function(startPos1, endPos1,seq1,startPos2, endPos2,seq2, fasta){
    
    
    table = data.frame(stringsAsFactors = FALSE)
    
    As = 15
    #make a new sequence with AAA as the middle 
    seq = paste(seq1, paste(rep("A", As), collapse = ""), seq2, sep = "")
    
    #make the constraints, you want the 40 A's not to be paired
    start = nchar(seq1) + 1

    #P i 0 k
    row = c("P", start, 0, As)
    write.table(t(as.data.frame(row)), file = "constraints.txt", quote = F, row.names = F, col.names = F)
    
    
    command = paste("echo \">",startPos1,"-", endPos1,"\n",seq,"\" | RNAfold  --noPS --constraint=constraints.txt ", sep = "")
    
    
    #create vienna command and run capturing output
    #command = paste("echo \">",startPos1,"-", endPos1,"\n",seq1,"\n>",startPos2,"-",endPos2,"\n",seq2,"\" | RNAduplex", sep = "")
    x = system(command,intern = T)
    
    
    #extract the vienna and the start locations
    vienna = sub(" .*","",x[3])
    
    helix = viennaToHelix(vienna )
    
    # to anything below a certain number we want to add startPos1 -1 
    # and anything over a certain number we want to add startpos2 -2 
    cutoff = nchar(seq1) + As
    helix2 = helix
    for(i in c(1,2)){
        for(j in 1:nrow(helix)){
            if(helix[j,i] > cutoff){
                helix2[j,i] = helix2[j,i] + startPos2 - 1
            }else{
                helix2[j,i] = helix2[j,i] + startPos1 - 1
            }
        
        }
    }
    helix = helix2
    helix = helix[,-c(3,4)]
    helix$rl = helix$i 
    helix$rr = helix$j 
    helix$pl = fasta[[1]]$NR_003286.4[helix$rl]
    helix$pr = fasta[[1]]$NR_003286.4[helix$rr]
    helix$veinna = vienna
    helix$seq1new = seq1
    helix$seq2new = seq2
    colnames(helix) = c("l","r","rl","rr","pl","pr","vienna","seq1new","seq2new")
    helix$check = 1
    aggTable = aggregate(helix$check, by=list(helix$rl,helix$rr,helix$pl,helix$pr,helix$vienna,helix$seq1new,helix$seq2new), FUN = sum )
    return(aggTable)
}











#' findBasePairsRNAfold
#'
#' Folds the clusters using Vienna RNA duplex
#' 
#' @param  startPos1  Start of the cluster side x
#' @param  endPos1  End of the cluster side x
#' @param  seq1 Sequence of x
#' @param  fasta \code{rnaRefs}
#' 
#' @return A table of clusters and coordinates with folds
#' 
#' @examples
#' 
#' @export
findBasePairsRNAfold = function(startPos, endPos, seqs,fasta){
    
    
    table = data.frame(stringsAsFactors = FALSE)
    #  for(i in 1:10){
    #command = paste("echo \">",startPos,"-",endPos,"\n",seqs,"\" | RNAfold --ImFeelingLucky", sep = "")
    command = paste("echo \">",startPos,"-",endPos,"\n",seqs,"\" | RNAfold  --noPS", sep = "")
    x = system(command,intern = T)
    
    #extract the vienna and the start locations
    vienna = sub(" .*","",x[3])
    
    helix = viennaToHelix(vienna )
    
    helix = helix[,-c(3,4)]
    helix$rl = helix$i + startPos -1
    helix$rr = helix$j + startPos  -1
    helix$pl = fasta[[1]]$NR_003286.4[helix$rl]
    helix$pr =fasta[[1]]$NR_003286.4[helix$rr]
    helix$veinna = vienna
    helix$seq1new = seqs
    helix$seq2new = ""
    colnames(helix) = c("l","r","rl","rr","pl","pr","vienna","seq1new","seq2new")
    helix$check = 1
    aggTable = aggregate(helix$check, by=list(helix$rl,helix$rr,helix$pl,helix$pr,helix$vienna,helix$seq1new,helix$seq2new), FUN = sum )
    return(aggTable)
}




#' getClusterClusterShortRangeWhole
#'
#' Decides if a cluster is long or short range
#' then either grabs the whole sequence or the sequence of the two sides of the 
#' interaction separately.
#' 
#' @return The same table with an extra column
#' 
#' @examples
#' 
#' @export
getClusterClusterShortRangeWhole = function(cluster, seq){
    fasta = seq
    if(cluster$rs - cluster$le <= 40){
        coords = list()
        coords[[1]] = cluster$ls
        coords[[2]] = cluster$re
        seqs1 = paste(fasta[[1]][coords[[1]]:coords[[2]]], collapse = "")
        type = "short"
        seqs2 = ""
        return(list(type, seqs1, seqs2, cluster))
    }else{
        coords = list()
        coords[[1]] = cluster$ls
        coords[[2]] = cluster$le
        coords[[3]] = cluster$rs
        coords[[4]] = cluster$re
        seqs1 = paste(fasta[[1]][coords[[1]]:coords[[2]]], collapse = "")
        seqs2 = paste(fasta[[1]][coords[[3]]:coords[[4]]], collapse = "")
        type = "long"
        return(list(type, seqs1, seqs2, cluster))
    }
}






#' findBasePairsRNAfold
#'
#' Folds the clusters using Vienna RNA duplex
#' 
#' @param  startPos1  Start of the cluster side x
#' @param  endPos1  End of the cluster side x
#' @param  seq1 Sequence of x
#' @param  fasta \code{rnaRefs}
#' 
#' @return A table of clusters and coordinates with folds
#' 
#' @examples
#' 
#' @export
findBasePairsRNAfold2 = function(startPos, endPos, seqs,fasta){
    
    
    table = data.frame(stringsAsFactors = FALSE)
    #  for(i in 1:10){
    #command = paste("echo \">",startPos,"-",endPos,"\n",seqs,"\" | RNAfold --ImFeelingLucky", sep = "")
    command = paste("echo \">",startPos,"-",endPos,"\n",seqs,"\" | RNAfold  --noPS", sep = "")
    x = system(command,intern = T)
    
    #extract the vienna and the start locations
    vienna = sub(" .*","",x[3])
    
    helix = viennaToHelix(vienna )
    
    helix = helix[,-c(3,4)]
    helix$rl = helix$i + startPos -1
    helix$rr = helix$j + startPos  -1
    helix$pl = fasta[[1]]$NR_003286.4[helix$rl]
    helix$pr =fasta[[1]]$NR_003286.4[helix$rr]
    helix$veinna = vienna
    helix$seq1new = seqs
    helix$seq2new = ""
    colnames(helix) = c("l","r","rl","rr","pl","pr","vienna","seq1new","seq2new")
    helix$check = 1
    aggTable = aggregate(helix$check, by=list(helix$rl,helix$rr,helix$pl,helix$pr,helix$vienna,helix$seq1new,helix$seq2new), FUN = sum )
    return(aggTable)
}





