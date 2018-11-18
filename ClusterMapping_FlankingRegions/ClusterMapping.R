ClusterMapping <- function() {
  library(ggrepel)
  library(gridExtra)
  library(ggplot2)
  # Assign clusters to gene file
  # ClusterDistance: Minimum distance between clusters
  ClusterDistance <- 3500
  
  # FlankingDistance: This is the distance between the end of left flanking region and the begining of the cluster
  # And, the distance between the end of the cluster and the begining of the right flanking region
  FlankingDistance <- 1000
  # teh width of flanking regions around the custers
  Flankingwidth <- 2000
  
  # The integrated genefile found by tse and ara
  geneDF <- integrated_tse_ara
  
  # assignGeneID() :
  # 1. assign two columns "tseGeneID" and "araGeneID" to the genefile
  # 2.
  geneDF <- assignGeneID(geneDF)
  
  # split the genefile into a list of Data frames. One DF for each genome tRNAs
  # Assign cluster number to tRNAs in each Dataframe
  genome_list = split(geneDF, f = geneDF$sourceOrg)
  for (i in 1:length(genome_list)) {
    genome_list[[i]] <- assignCluster(genome_list[[i]], ClusterDistance)
  }
  
  # Make a list of cluster DataFrames. one DF for each genome. Each DF will have the fallowing columns:
  # GenomesCluster_list: "cluster" "tsetRNAIDs" "aratRNAIDs" "seqname" "LeftFlank_start" "LeftFlank_end" "RightFlank_start" "RightFlank_end"
  # If the left flanking region is out of bound, set it to the begining of the sequence till the begining of the cluster.
  GenomesCluster_list <- list()
  for (i in 1:length(genome_list)) {
    GenomesCluster_list[[i]] <-
      makeClusterDf(genome_list[[i]], FlankingDistance, Flankingwidth)
  }
  names(GenomesCluster_list) <- names(genome_list)
  
  # Write each genome's Flanking regions in a file with pattern: "FlankingR_<GenomeName>.txt"
  # (For each of these files we will to run a bash script to align all the flanking regions to other genomes)
  for (i in 1:length(GenomesCluster_list)) {
    # write the DF in a file to be processed by the bash
    filepath <-
      paste(
        "/home/fatemeh/Leishmania_Aug2018/ClusteringAndGenomeSynteny/ClusterMapping_FlankingRegions/GenomesFlankingRegions/",
        "FlankingR_",
        names(GenomesCluster_list[i]),
        ".txt",
        sep = ""
      )
    write.csv(
      GenomesCluster_list[[i]],
      quote = FALSE,
      row.names = FALSE,
      file = filepath
    )
  }
  
  # ________________________________________________________________________________________________________________________________________________________
  #
  # Bash script gives us us a folder for each genomes. each folder has all the files for "pairwise flanking region alignment" with other genomes
  # in function genomeDistance we will read all the file and match the clusters between each pair of genomes
  #
  # ________________________________________________________________________________________________________________________________________________________
  genomeDistance(GenomesCluster_list)
  
  
  
}


clustergenome <- function(GenomeDistance_Average) {
  library(cluster)    # clustering algorithms
  library(factoextra) # clustering visualization
  library(dendextend) # for comparing two dendrograms
  distanceDF <- as.data.frame(GenomeDistance_Average)
  names(distanceDF) <- genomenames
  varzero <- as.data.frame(which(apply(distanceDF, 2, var) == 0))
  distanceDF2 <- distanceDF[-varzero[, 1], -varzero[, 1]]
  genomenames2 <- genomenames[!(genomenames %in% rownames(varzero))]
  #names(distanceDF2) <- genomenames
  rownames(distanceDF2) <- colnames(distanceDF2)
  for (i in 1:ncol(distanceDF2)) {
    distanceDF2[, i] <- (100 - distanceDF2[, i])
  }
  
  hc <- hclust(as.dist(distanceDF2), method = 'ward.D2')
  plot(
    hc,
    cex = 0.6,
    hang = -1,
    sub = NULL,
    xlab = ""
  )
  
  # The height of the cut to the dendrogram controls the number of clusters obtained. It plays the same role as the k in k-means clustering.
  #
  rect.hclust(hc, k = 9, border = 2:5)
  sub_groups <- cutree(hc, k = 9)
  fviz_cluster(list(data = distanceDF2, cluster = sub_groups))
  
}
assignGeneID <- function(geneDF) {
  geneDF$tseID <- character(length = nrow(geneDF))
  geneDF$araID <- character(length = nrow(geneDF))
  for (i in 1:nrow(geneDF)) {
    id <- tolower(geneDF$tseidentity[i])
    if (id == "ala")
      geneDF$tseID[i] <- paste("A", geneDF$tseac[i], sep = "")
    if (id == "arg")
      geneDF$tseID[i] <- paste("R", geneDF$tseac[i], sep = "")
    if (id == "asn")
      geneDF$tseID[i] <- paste("N", geneDF$tseac[i], sep =  "")
    if (id == "asp")
      geneDF$tseID[i] <- paste("D", geneDF$tseac[i], sep =  "")
    if (id == "cys")
      geneDF$tseID[i] <- paste("C", geneDF$tseac[i], sep =  "")
    if (id == "gln")
      geneDF$tseID[i] <- paste("Q", geneDF$tseac[i], sep =  "")
    if (id == "glu")
      geneDF$tseID[i] <- paste("E", geneDF$tseac[i], sep = "")
    if (id == "gly")
      geneDF$tseID[i] <- paste("G", geneDF$tseac[i], sep = "")
    if (id == "his")
      geneDF$tseID[i] <- paste("H", geneDF$tseac[i], sep =  "")
    if (id == "ile")
      geneDF$tseID[i] <- paste("l", geneDF$tseac[i], sep =  "")
    if (id == "start")
      geneDF$tseID[i] <- paste("start", geneDF$tseac[i], sep = "")
    if (id == "leu")
      geneDF$tseID[i] <- paste("L", geneDF$tseac[i], sep = "")
    if (id == "lys")
      geneDF$tseID[i] <-  paste("K", geneDF$tseac[i], sep = "")
    if (id == "met")
      geneDF$tseID[i] <-  paste("M", geneDF$tseac[i], sep =  "")
    if (id == "phe")
      geneDF$tseID[i] <-  paste("F", geneDF$tseac[i], sep =  "")
    if (id == "pro")
      geneDF$tseID[i] <-  paste("P", geneDF$tseac[i], sep =  "")
    if (id == "ser")
      geneDF$tseID[i] <-  paste("S", geneDF$tseac[i], sep =  "")
    if (id == "thr")
      geneDF$tseID[i] <-  paste("T", geneDF$tseac[i], sep =  "")
    if (id == "trp")
      geneDF$tseID[i] <- paste("W", geneDF$tseac[i], sep =  "")
    if (id == "tyr")
      geneDF$tseID[i] <-
        paste("Y", geneDF$tseac[i], sep = "")
    if (id == "val")
      geneDF$tseID[i] <-
        paste("V", geneDF$tseac[i], sep =  "")
    if (id == "sec")
      geneDF$tseID[i] <-
        paste("U", geneDF$tseac[i], sep =  "")
    if (id == "stop")
      geneDF$tseID[i] <-
        paste("stop", geneDF$tseac[i], sep = "")
  }
  
  #______________________________________________________
  for (i in 1:nrow(geneDF)) {
    id <- tolower(geneDF$araidentity[i])
    if (id == "ala")
      geneDF$araID[i] <- paste("A", geneDF$araac[i], sep = "")
    if (id == "arg")
      geneDF$araID[i] <- paste("R", geneDF$araac[i], sep = "")
    if (id == "asn")
      geneDF$araID[i] <- paste("N", geneDF$araac[i], sep =  "")
    if (id == "asp")
      geneDF$araID[i] <- paste("D", geneDF$araac[i], sep =  "")
    if (id == "cys")
      geneDF$araID[i] <- paste("C", geneDF$araac[i], sep =  "")
    if (id == "gln")
      geneDF$araID[i] <- paste("Q", geneDF$araac[i], sep =  "")
    if (id == "glu")
      geneDF$araID[i] <- paste("E", geneDF$araac[i], sep = "")
    if (id == "gly")
      geneDF$araID[i] <- paste("G", geneDF$araac[i], sep = "")
    if (id == "his")
      geneDF$araID[i] <- paste("H", geneDF$araac[i], sep =  "")
    if (id == "ile")
      geneDF$araID[i] <- paste("l", geneDF$araac[i], sep =  "")
    if (id == "start")
      geneDF$araID[i] <- paste("start", geneDF$araac[i], sep = "")
    if (id == "leu")
      geneDF$araID[i] <- paste("L", geneDF$araac[i], sep = "")
    if (id == "lys")
      geneDF$araID[i] <-  paste("K", geneDF$araac[i], sep = "")
    if (id == "met")
      geneDF$araID[i] <-  paste("M", geneDF$araac[i], sep =  "")
    if (id == "phe")
      geneDF$araID[i] <-  paste("F", geneDF$araac[i], sep =  "")
    if (id == "pro")
      geneDF$araID[i] <-  paste("P", geneDF$araac[i], sep =  "")
    if (id == "ser")
      geneDF$araID[i] <-  paste("S", geneDF$araac[i], sep =  "")
    if (id == "thr")
      geneDF$araID[i] <-  paste("T", geneDF$araac[i], sep =  "")
    if (id == "trp")
      geneDF$araID[i] <- paste("W", geneDF$araac[i], sep =  "")
    if (id == "tyr")
      geneDF$araID[i] <-
        paste("Y", geneDF$araac[i], sep = "")
    if (id == "val")
      geneDF$araID[i] <-
        paste("V", geneDF$araac[i], sep =  "")
    if (id == "sec")
      geneDF$araID[i] <-
        paste("U", geneDF$araac[i], sep =  "")
    if (id == "stop")
      geneDF$araID[i] <-
        paste("stop", geneDF$araac[i], sep = "")
    
  }
  geneDF
}
assignCluster <- function(geneDF, clusterdistance) {
  isboth <- geneDF$foundby == "both"
  tRNAGenesDF <- geneDF[isboth,]
  tRNAGenesDF$tsebegin <- as.integer(tRNAGenesDF$tsebegin)
  tRNAGenesDF$tseend <- as.integer(tRNAGenesDF$tseend)
  
  # they should be sorted based on seqname and start
  # remove the spaces in seqname
  
  tRNAGenesDF <-
    tRNAGenesDF[order(tRNAGenesDF$sourceseq, tRNAGenesDF$tsebegin),]
  
  tRNAGenesDF$cluster <- 0
  tRNAGenesDF$cluster[1] <- 1
  setnumber <- 1
  for (i in 1:(nrow(tRNAGenesDF) - 1)) {
    geneDist <-
      abs(tRNAGenesDF$tsebegin[i + 1] - tRNAGenesDF$tsebegin[i])
    if ((geneDist < clusterdistance) &&
        (tRNAGenesDF$sourceOrg[i] == tRNAGenesDF$sourceOrg[i + 1]) &&
        (tRNAGenesDF$sourceseq[i + 1] == tRNAGenesDF$sourceseq[i]))
      tRNAGenesDF$cluster[i + 1] <- setnumber
    else
    {
      setnumber <- setnumber + 1
      tRNAGenesDF$cluster[i + 1] <- setnumber
    }
  }
  # generating breaks
  tRNAGenesDF
}
makeClusterDf <- function(tRNAdf, FlankingDistance, Flankingwidth) {
  m <- matrix(nrow = tRNAdf$cluster[nrow(tRNAdf)], ncol = 10)
  outputdf <- as.data.frame(m)
  names(outputdf) <-
    c(
      "cluster",
      "tsetRNAIDs",
      "aratRNAIDs",
      "seqname",
      "clusterStart",
      "clusterEnd",
      "LeftFlank_start",
      "LeftFlank_end",
      "RightFlank_start",
      "RightFlank_end"
    )
  outputdf$cluster <- seq(1, tRNAdf$cluster[nrow(tRNAdf)])
  
  for (i in 1:tRNAdf$cluster[nrow(tRNAdf)]) {
    mydf <- tRNAdf[tRNAdf$cluster == i, ]
    outputdf$tsetRNAIDs[i] <- paste(mydf$tseID, collapse = " ")
    outputdf$aratRNAIDs[i] <- paste(mydf$araID, collapse = " ")
    outputdf$seqname[i] <- mydf$sourceseq[1]
    outputdf$clusterStart[i] <- mydf$tsebegin[1]
    outputdf$clusterEnd[i] <- mydf$tseend[nrow(mydf)]
    mydf <- mydf[order(mydf$tsebegin), ]
    LeftFlankS <-
      mydf$tsebegin[1] - (FlankingDistance + Flankingwidth)
    LeftFlankE <- LeftFlankS + Flankingwidth
    if (LeftFlankS < 0)
    {
      LeftFlankS <- 1
      LeftFlankE <- mydf$tsebegin[1] - 10
    }
    RightFlank <-
      mydf$tseend[nrow(mydf)] + FlankingDistance
    # if(RightFlank <- seqlength)
    
    outputdf$LeftFlank_start[i] <- LeftFlankS
    outputdf$LeftFlank_end[i] <- LeftFlankE
    
    outputdf$RightFlank_start[i] <- RightFlank
    outputdf$RightFlank_end[i] <- RightFlank + Flankingwidth
  }
  
  outputdf
}
matchqryclusters <- function(AlignmentDF, qrygenomeclusters){
  naends <-
    (is.na(AlignmentDF$EndStart) & !is.na(AlignmentDF$LeftStart))
  
  AlignmentDF[naends, ]$EndStart <-
    AlignmentDF[naends, ]$LeftStart + 5000
  AlignmentDF$QryCluster <- 0
  AlignmentDF$QrytsetRNAs <- " "
  AlignmentDF$QryaratRNAs <- " "
  for (j in 1:nrow(AlignmentDF)) {
    for (i in 1:nrow(qrygenomeclusters)) {
      if (qrygenomeclusters$seqname[i] == AlignmentDF$QrySeqname[j])
      {
        if (qrygenomeclusters$clusterStart[i] > AlignmentDF$LeftStart[j])
        {
          if (qrygenomeclusters$clusterEnd[i] < AlignmentDF$EndStart[j])
          {
            AlignmentDF$QryCluster[j] <-  qrygenomeclusters$cluster[i]
            AlignmentDF$QrytsetRNAs[j] <-
              qrygenomeclusters$tsetRNAIDs[i]
            AlignmentDF$QryaratRNAs[j] <-
              qrygenomeclusters$aratRNAIDs[i]
          }
        }
      }
    }
  }
  AlignmentDF
}

genomeDistance <- function() {
  #  in every round AlignmentDF has the matched clusters
  GenomeDistance <- matrix(nrow = 46, ncol = 46)
  dirpath <-
    "/home/fatemeh/Leishmania_Aug2018/ClusteringAndGenomeSynteny/ClusterMapping_FlankingRegionsCOPY/Alignments"
  Alignmentdirs <- list.dirs(path = dirpath)
  Alignmentdirs <- Alignmentdirs[2:length(Alignmentdirs)]
  for (d in 1:length(Alignmentdirs)) {
    print("**d**")
    print(d)
    filepath <- Alignmentdirs[d]
    tempf <- unlist(strsplit(Alignmentdirs[d], split = "_"))
    refgenome <- tempf[length(tempf)]
    
    pairwisefiles_Left <-
      list.files(path = filepath, pattern = "LeftF_AlTo")
    pairwisefiles_Right <-
      gsub("LeftF", "RightF", pairwisefiles_Left)
    
    pairwisefiles_Left <-
      paste(filepath, "/", pairwisefiles_Left, sep = "")
    pairwisefiles_Right <-
      paste(filepath, "/", pairwisefiles_Right, sep = "")
    
    for (i in 1:length(pairwisefiles_Left)) {
      # print("**i**")
      # print(i)
      info1 = file.info(pairwisefiles_Left[i])
      info2 = file.info(pairwisefiles_Right[i])
      if (info1$size <= 60 | info2$size <= 60)
      {
        print("emptyfile")
        GenomeDistance[d, i] <- -1
      }
      else{
        LeftLine <- readLines(pairwisefiles_Left[i])
        RightLine <- readLines(pairwisefiles_Right[i])
        AlignmentDF <-
          data.frame(character(length(LeftLine)), integer(length(LeftLine)))
        names(AlignmentDF) <- c("QrySeqname", "LeftStart")
        AlignmentDF_Right <-
          data.frame(character(length(RightLine)), integer(length(RightLine)))
        names(AlignmentDF_Right) <- c("QrySeqname", "RightStart")
        AlignmentDF_Right$QrySeqname <-
          as.character(AlignmentDF_Right$QrySeqname)
        AlignmentDF$QrySeqname <-
          as.character(AlignmentDF$QrySeqname)
        # their default is NA
        AlignmentDF_Right$RightStart <- NA
        AlignmentDF$LeftStart <- NA
        
        for (j in 1:length(LeftLine)) {
          if (LeftLine[j] == "")
            AlignmentDF$QrySeqname[j] <- "N"
          if (RightLine[j] == "")
            AlignmentDF_Right$QrySeqname[j] <- "N"
          if (LeftLine[j] != "")
            AlignmentDF[j,] <-
              unlist(strsplit(LeftLine[j], split = " "))
          if (RightLine[j] != "")
            AlignmentDF_Right[j,] <-
              unlist(strsplit(RightLine[j], split = " "))
        }
        
        AlignmentDF_Right$RightStart <-
          as.integer(AlignmentDF_Right$RightStart)
        AlignmentDF$LeftStart <- as.integer(AlignmentDF$LeftStart)
        AlignmentDF$EndStart <- AlignmentDF_Right$RightStart
        # if start_Left > EndStart then switch them
        
        # find all the tRNAs from Qry that are in between the flanking positions
        qrygenomename <- gsub(filepath, "", pairwisefiles_Left[i])
        qrygenomename <- gsub(".txt", "", qrygenomename)
        temp <- unlist(strsplit(qrygenomename, split = "_"))
        qrygenomename <- temp[length(temp)]
        qrygenomeclusters <-
          GenomesCluster_list[names(GenomesCluster_list) == qrygenomename][[1]]
        AlignmentDF <-
          matchqryclusters(AlignmentDF, qrygenomeclusters)
        refgenomeclusters <-
          GenomesCluster_list[names(GenomesCluster_list) == refgenome][[1]]
        #add three columns to refgenomeclusters one for cluster number and one for tsetRNAIDS, and one for aratRNAIDs
        column1 <- paste(qrygenomename, "cluster", sep = "")
        column2 <- paste(qrygenomename, "tsetRNAIDS", sep = "")
        column3 <- paste(qrygenomename, "aratRNAIDs", sep = "")
        refgenomeclusters$column1 <- AlignmentDF$QryCluster
        refgenomeclusters$column2 <- AlignmentDF$QrytsetRNAs
        refgenomeclusters$column3 <- AlignmentDF$QryaratRNAs
        names(refgenomeclusters)[ncol(refgenomeclusters)] <-
          paste(qrygenomename, "aratRNAIDs", sep = "")
        names(refgenomeclusters)[ncol(refgenomeclusters) - 1] <-
          paste(qrygenomename, "tsetRNAIDS", sep = "")
        names(refgenomeclusters)[ncol(refgenomeclusters) - 2] <-
          paste(qrygenomename, "cluster", sep = "")
        
        # number of matched clusters for reference
        
        GenomeDistance[d, i] <-
          (sum(AlignmentDF$QryCluster != 0) / nrow(AlignmentDF)) * 100
      }
    }
  }
  
  # extract the genomenames from the Alignmentdirs which has all the folders of pairwise alignment for each genome ____
  
  genomenames <- character(length = length(Alignmentdirs))
  for (d in 1:length(Alignmentdirs)) {
    filepath <- Alignmentdirs[d]
    tempf <- unlist(strsplit(Alignmentdirs[d], split = "_"))
    refgenome <- tempf[length(tempf)]
    genomenames[d] <- refgenome
  }
  tempm <- GenomeDistance
  
  # for each pair of genomes i and j take the average of teh distance between i,j and j,i
  GenomeDistance_Average <- matrix(nrow = 46,
                                   ncol = 46,
                                   data = 0)
  for (j in 1:nrow(GenomeDistance_Average)) {
    for (i in 1:nrow(GenomeDistance_Average)) {
      if (tempm[j, i] == -1)
        tempm[j, i] = 0
      if (tempm[i, j] == -1)
        tempm[i, j] = 0
      GenomeDistance_Average[j, i] <- (tempm[j, i] + tempm[j, i]) / 2
      GenomeDistance_Average[i, j] <- (tempm[j, i] + tempm[j, i]) / 2
    }
    # pass teh distance matrix to clustergenome to plot the cluster of genomes
    clustergenome(GenomeDistance_Average)
  }
}
