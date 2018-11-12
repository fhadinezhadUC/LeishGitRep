# This script will take the integrated gene file (output of Integrate_Tse_Ara script), 
# removes the introns found by TSE, removes the variable arm
# output the result in file Integrated_Genes_NoVarIntron.fasta 
removeIntronVararm <- function(integrated_tse_ara) {
  
  # we use genes that are found by both tse and ara and their identity matches. 
  both <- integrated_tse_ara$foundby == "both"
  genefile <- integrated_tse_ara[both, ]
  isvagueIdentity <- (genefile$tseidentity != genefile$araidentity)
  genefile <- genefile[!isvagueIdentity,]
  nonpseudo <- genefile$tsenote=="notfound"
  genefile <- genefile[nonpseudo,]
  
  # add one column to the genefile dataframe and call it ClassID 
  genefile$classID <- genefile$tseidentity 
  for (i in 1:nrow(genefile)) {
    if(genefile$classID[i] == "Ala")
      genefile$classID[i] = "A"
    if(genefile$classID[i] == "Arg")
      genefile$classID[i] = "R"
    if(genefile$classID[i] == "Asn")
      genefile$classID[i] = "N"
    if(genefile$classID[i] == "Asp")
      genefile$classID[i] = "D"
    if(genefile$classID[i] == "Cys")
      genefile$classID[i] = "C"
    if(genefile$classID[i] == "Gln")
      genefile$classID[i] = "Q"
    if(genefile$classID[i] == "Glu")
      genefile$classID[i] = "E"
    if(genefile$classID[i] == "Gly")
      genefile$classID[i] = "G"
    if(genefile$classID[i] == "His")
      genefile$classID[i] = "H"
    if(genefile$classID[i] == "Ile")
      genefile$classID[i] = "I"
    if(genefile$classID[i] == "Leu")
      genefile$classID[i] = "L"
    if(genefile$classID[i] == "Lys")
      genefile$classID[i] = "K"
    if(genefile$classID[i] == "Met")
      genefile$classID[i] = "M"
    if(genefile$classID[i] == "Phe")
      genefile$classID[i] = "F"
    if(genefile$classID[i] == "Pro")
      genefile$classID[i] = "P"
    if(genefile$classID[i] == "SeC")
      genefile$classID[i] = "Z"
    if(genefile$classID[i] == "Ser")
      genefile$classID[i] = "S"
    if(genefile$classID[i] == "Thr")
      genefile$classID[i] = "T"
    if(genefile$classID[i] == "Trp")
      genefile$classID[i] = "W"
    if(genefile$classID[i] == "Tyr")
      genefile$classID[i] = "Y"
    if(genefile$classID[i] == "Val")
      genefile$classID[i] = "V"
  }
  
  # mark initiators with X 
  # pass the genefile to the initiatorDetection function to mark the initiators 
  genefile <- initiatorDetecting(genefile)
  
  genefile_Nointron <-
    data.frame(paste(genefile$geneid,"_",genefile$classID,sep = ""),
               genefile$tsegeneseq,
               genefile$tsegeness,
               genefile$arageness)
  names(genefile_Nointron) <- c("GeneID", "GeneSeq", "SSTSE", "SSARA")
  
  genefile$arageneseq <- as.character(genefile$arageneseq)
  genefile$tsegeneseq <- as.character(genefile$tsegeneseq)
  genefile_Nointron$GeneID <- as.character(genefile_Nointron$GeneID)
  genefile_Nointron$GeneSeq <-
    as.character(genefile_Nointron$GeneSeq)
  genefile_Nointron$SSARA <- as.character(genefile_Nointron$SSARA)
  genefile_Nointron$SSTSE <- as.character(genefile_Nointron$SSTSE)
  
  for (i in 1:nrow(genefile)) {
    if (genefile$tseintronbegin[i] != "notfound" &
        genefile$tseintronbegin[i] != "0")
    {
      # find the location on sequence
      start <-
        as.integer(genefile$tseintronbegin[i]) - as.integer(genefile$tsebegin[i]) + 1
      firstchunk <- substring(genefile$tsegeneseq[i], 1, start - 1)
      SSTSE_firstchunk <- substring(genefile$tsegeness[i], 1, start - 1)
      
      intronlen <-
        as.integer(genefile$tseintronend[i]) - as.integer(genefile$tseintronbegin[i]) + 1
      secondchunk <-
        substring(genefile$tsegeneseq[i],
                  start + intronlen,
                  nchar(genefile$tsegeneseq[i]))
      SSTSE_secondchunk <-
        substring(genefile$tsegeness[i],
                  start + intronlen,
                  nchar(genefile$tsegeness[i]))
      
      genefile_Nointron$GeneSeq[i] <-
        paste(firstchunk, secondchunk, sep = "")
      genefile_Nointron$SSTSE[i] <-
        paste(SSTSE_firstchunk, SSTSE_secondchunk, sep = "")
      
    }
    else if (genefile$araintronbegin[i] != "nointron" &
             genefile$araintronbegin[i]  != "notfound" &
             genefile$tseintronbegin[i] != "0")
    {
      print("hmm!")
      print(i)
      start <- as.integer(genefile$araintronbegin[i])
      firstchunk <- substring(genefile$arageneseq[i], 1, start - 1)
      SSARA_firstchunk <- substring(genefile$arageness[i], 1, start - 1)
      
      intronlen <-
        as.integer(genefile$araintronend[i]) - as.integer(genefile$araintronbegin[i])
      secondchunk <-
        substring(genefile$arageneseq[i],
                  start + intronlen,
                  nchar(genefile$arageneseq[i]))
      SSARA_secondchunk <-
        substring(genefile$arageness[i],
                  start + intronlen,
                  nchar(genefile$arageness[i]))
      
      genefile_Nointron$GeneSeq[i] <-
        paste(firstchunk, secondchunk, sep = "")
      genefile_Nointron$SSARA[i] <-
        paste(SSARA_firstchunk, SSARA_secondchunk, sep = "")
    }
    
    else if ((genefile$tseintronbegin[i] == "notfound" |
              genefile$tseintronbegin[i] == "0") &
             (genefile$araintronbegin[i] == "nointron" |
              genefile$araintronbegin[i]  == "notfound"))
    {
      genefile_Nointron$GeneSeq[i] <- genefile$tsegeneseq[i]
      genefile_Nointron$SSTSE[i] <-  genefile$tsegeness[i]
      if (genefile$tsegeneseq[i] == "notfound")
      {
        genefile_Nointron$GeneSeq[i] <- genefile$arageneseq[i]
        genefile_Nointron$SSARA[i] <-  genefile$arageness[i]
        print("hey!")
      }
    }
  }
  
  genefile_Nointron <- removeVariableArm(genefile_Nointron)
  write.fwf(
    data.frame(
      paste(">", genefile_Nointron$GeneID, sep = ""),
      genefile_Nointron$GeneSeq
    ),
    file = "/home/fatemeh/Leishmania_Aug2018/phyloclassification/Integrated_Genes_NoVarIntron.fasta",
    sep = "\n",
    colnames = FALSE
  )
  
  
}

initiatorDetecting <- function(geneDF){
  
  isnotcat <- geneDF$araac!="(cat)" & geneDF$tseac!="CAT"
  notcatDF <- geneDF[isnotcat,]
  
  iscat <- geneDF$araac=="(cat)" & geneDF$tseac=="CAT"
  catDF <- geneDF[iscat,]
  geneseqs <- character(length = nrow(catDF))
  geneseqs <- catDF$tsegeneseq
  
  m <- matrix(nrow = nrow(catDF), ncol = nrow(catDF))
  distanceDF <- as.data.frame(m)
  for (i in 1:length(geneseqs)) {
    for (j in 1:length(geneseqs)) {
      distanceDF[i, j] <- adist(geneseqs[i], geneseqs[j])
    }
  }
  
  
  #library(tidyverse)  # data manipulation
  library(cluster)    # clustering algorithms
  library(factoextra) # clustering visualization
  library(dendextend) # for comparing two dendrograms

  hc <- hclust(as.dist(distanceDF), method='ward.D2')
  plot(hc, cex = 0.6, hang = -1)
  

  rect.hclust(hc, k = 3, border = 2:5)
  sub_groups <- cutree(hc, k = 3)
  fviz_cluster(list(data = distanceDF, cluster = sub_groups))
  
  subg1 <- sub_groups == 1 
  subg2 <- sub_groups == 2
  subg3 <- sub_groups == 3

  clus1DF <- catDF[subg1,]
  clus2DF <- catDF[subg2,]
  clus3DF <- catDF[subg3,]
  
  clus1DF$classID = "X"
  
  Editedgenefile <- rbind(clus1DF,clus2DF,clus3DF,notcatDF)
  Editedgenefile
  
}
removeVariableArm <- function(genefile_Nointron) {
  for (i in 1:nrow(genefile_Nointron)) {
    # count number of "><"
    numarms <- countarms(genefile_Nointron[i, ]$SSTSE)
    if (numarms == 4)
    {
      varcor <-
        removearm(genefile_Nointron[i, ]$GeneSeq, genefile_Nointron[i, ]$SSTSE)
      firstchunk <-
        substring(genefile_Nointron$GeneSeq[i], 1, varcor[1] - 1)
      SSTSE_firstchunk <-
        substring(genefile_Nointron$SSTSE[i], 1, varcor[1] - 1)
      
      secondchunk <-
        substring(genefile_Nointron$GeneSeq[i],
                  varcor[2] + 1,
                  nchar(genefile_Nointron$GeneSeq[i]))
      SSTSE_secondchunk <-
        substring(genefile_Nointron$SSTSE[i],
                  varcor[2] + 1,
                  nchar(genefile_Nointron$SSTSE[i]))
      
      genefile_Nointron$GeneSeq[i] <-
        paste(firstchunk, secondchunk, sep = "")
      genefile_Nointron$SSTSE[i] <-
        paste(SSTSE_firstchunk, SSTSE_secondchunk, sep = "")
    }
    if (numarms < 3)
    {
      print("We have genes with unusual structure!")
    }
  }
  genefile_Nointron
}

countarms <- function(geneSS) {
  geneSSarr <-
    substring(geneSS, seq(1, nchar(geneSS), 1), seq(1, nchar(geneSS), 1))
  counter = 0
  flag = 0
  vararmend = 0
  for (i in 1:length(geneSSarr)) {
    if (geneSSarr[i] == "<" & flag == 0)
    {
      counter <- counter + 1
      flag = 1
      
    }
    if (geneSSarr[i] == ">")
      flag = 0
  }
  
  counter
}


removearm <- function(geneSeq, geneSS)
{
  flag = 0
  counter = 0
  forward = 1
  varbeg = 0
  varend = 0
  geneSSarr <-
    substring(geneSS, seq(1, nchar(geneSS), 1), seq(1, nchar(geneSS), 1))
  for (i in 1:length(geneSSarr)) {
    if (geneSSarr[i] == ">" & flag == 0)
    {
      counter <- counter + 1
      flag = 1
      if (counter == 3)
      {
        # this is the begining of the variable arm
        # count howmany > you have
        # go forward untill you see that many <
        forward = 0
        varbeg = i
        while (geneSSarr[i] != "<")
        {
          if (geneSSarr[i] == ">")
            forward = forward + 1
          i = i + 1
        }
        while (forward != 0) {
          if (geneSSarr[i] == "<")
            forward = forward - 1
          
          i = i + 1
        }
        varend = i - 1
      }
    }
    if (geneSSarr[i] == "<")
      flag = 0
  }
  
  print(varend)
  print(varbeg)
  varcor <- c(varbeg, varend)
  varcor
}