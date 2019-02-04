split_tRNAgene_with_genomename <- function() {
  # This function reads in the .fasta file containg aligned tRNA genes of HomoC and TryTryp genomes
  # Extract the functional classes from Headers
  # Makes a dataframe of "genomename", "sequence", "headers", "funclass"
  # Splits the dataframe into a list of data frames based on genomename
  # Writes each dataframe as a fasta file in a file with genome's names
  # At the end it calls function vidualization passing the list of data frames
  # Need to run script splitFuncClass.sh to split each file based on functional classes
  
  library(gdata)
  fastafile <-
    read.table(
      "/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfminput_2019/HomoTryTryp_EditedCovea.fasta",
      sep = "\n"
    )
  fastafile2 <- as.character(fastafile$V1)
  sequences <-
    as.character(fastafile2[seq(2, length(fastafile2), 2)])
  headers <- as.character(fastafile2[seq(1, length(fastafile2), 2)])
  
  genomenames2 <-
    lapply(X = headers , function(X)
      gsub("HOMO", "HOMO_", X, fixed = TRUE))
  genomenames2 <- as.character(genomenames2)
  genomenames3 <-
    lapply(X = genomenames2 , function(X)
      gsub(">", "", X, fixed = TRUE))
  genomenames3 <- as.character(genomenames3)
  genomenames4 <-
    lapply(X = genomenames3, function(X)
      unlist(strsplit(X, split = "_"))[1])
  genomenames4 <- as.character(genomenames4)
  
  headers <- gsub("\\s", "", headers)
  functionalclasses <-
    lapply(X = headers, function(X)
      unlist(strsplit(X, split = "_"))[length(unlist(strsplit(X, split = "_")))])
  
  functionalclasses <- as.character(functionalclasses)
  
  tRNAdf <-
    data.frame(genomenames4, sequences, genomenames2, functionalclasses)
  names(tRNAdf) <-
    c("genomename", "sequence", "headers", "funclass")
  
  genome_list = split(tRNAdf, f = tRNAdf$genomename)
  homo <-
    genome_list[names(genome_list) == "HOMO"][[1]]
  homo$funclass <-
    as.character(lapply(X = homo$headers , function(X)
      substr(X, 7, 7)))
  homo$funclass <- gsub("\\s", "", homo$funclass)
  homo$headers <- gsub("\\s", "", homo$headers)
  
  homoheader <- paste(homo$headers, homo$funclass, sep = "")
  genome_list[names(genome_list) == "HOMO"][[1]]$headers <-
    homoheader
  
  # make a file for each genome and put the genomes in there
  for (i in 1:length(genome_list)) {
    # write the DF in a file to be processed by the bash
    filepath <-
      paste(
        "/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfminput_2019/tRNADB/",
        names(genome_list[i]),
        ".fasta",
        sep = ""
      )
    
    write.fwf(
      data.frame(genome_list[[i]]$headers,
                 genome_list[[i]]$sequence),
      filepath,
      sep = "\n",
      colnames = FALSE
    )
  }
  vidualization(genome_list)
}
split_tRNAgene_with_ultiplegenomename <- function() {
  
  library(gdata)
  fastafile <-
    read.table(
      "/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfminput_2019/HomoTryTryp_EditedCovea.fasta",
      sep = "\n"
    )
  fastafile2 <- as.character(fastafile$V1)
  sequences <-
    as.character(fastafile2[seq(2, length(fastafile2), 2)])
  headers <- as.character(fastafile2[seq(1, length(fastafile2), 2)])

  genomenames2 <-
    lapply(X = headers , function(X)
      gsub("HOMO", "HOMO_", X, fixed = TRUE))
  genomenames2 <- as.character(genomenames2)
  genomenames3 <-
    lapply(X = genomenames2 , function(X)
      gsub(">", "", X, fixed = TRUE))
  genomenames3 <- as.character(genomenames3)
  genomenames4 <-
    lapply(X = genomenames3, function(X)
      unlist(strsplit(X, split = "_"))[1])
  genomenames4 <- as.character(genomenames4)
  
  # we will split genomes based on clusters of genomes
  # add another column "clustername" to tRNAdf to assign the new names to each cluster of genomes:
  # AfricanTrypanosome: TbruceigambienseDAL972, TbruceiLister427, TbruceiTREU927, TevansiSTIB805, TcongolenseIL3000, TvivaxY486
  # AmericanTrypanosome: TrangeliSC58, TcruziCLBrener, TcruziCLBrenerEsmeraldo-like, TcruziCLBrenerNon-Esmeraldo-like, TcruzicruziDm28c, TcruziDm28c, TcruziEsmeraldo, TcruziJRcl4, TcruzimarinkelleiB7, TcruziSylvioX10-1, TcruziSylvioX10-1-2012, TcruziTulacl2
  # Leishmania1: CfasciculataCfCl, LseymouriATCC30220LpyrrhocorisH10, LspMARLEM2494
  # Leishmania2: LmajorFriedlin, LmajorLV39c5, LmajorSD75, LdonovaniBHU1220, LdonovaniBPK282A1, LinfantumJPCM5, LturanicaLEM423, LarabicaLEM1108, LamazonensisMHOMBR71973M2269, LmexicanaMHOMGT2001U1103, LtropicaL590, LaethiopicaL147, LgerbilliLEM452
  # Lvianna: LbraziliensisMHOMBR75M2904, LbraziliensisMHOMBR75M2903, LpanamensisMHOMPA94PSC1, LpanamensisMHOMCOL81L13
  clusNames <- character(length = length(headers))
  clusNames <- genomenames4
  for (i in 1:length(clusNames)) {
    clusternames <- clusNames[i]
    if (clusternames == "TbruceigambienseDAL972" |
        clusternames == "TbruceiLister427" |
        clusternames == "TbruceiTREU927" |
        clusternames == "TevansiSTIB805" |
        clusternames == "TcongolenseIL3000" |
        clusternames == "TvivaxY486")
      clusNames[i] <- "AfricanTrypanosome"
    if (clusternames == "TgrayiANR4" |
        clusternames == "TrangeliSC58" |
        clusternames == "TcruziCLBrener" |
        clusternames == "TcruziCLBrenerEsmeraldo-like" |
        clusternames == "TcruziCLBrenerNon-Esmeraldo-like" |
        clusternames == "TcruzicruziDm28c" |
        clusternames == "TcruziDm28c" |
        clusternames == "TcruziEsmeraldo" |
        clusternames == "TcruziJRcl4" |
        clusternames == "TcruzimarinkelleiB7" |
        clusternames == "TcruziSylvioX10-1" |
        clusternames == "TcruziSylvioX10-1-2012" |
        clusternames == "TcruziTulacl2")
    clusNames[i] <- "AmericanTrypanosome"
    if (clusternames == "CfasciculataCfCl" |
        clusternames == "LseymouriATCC30220" |
        clusternames == "LpyrrhocorisH10" |
        clusternames == "LspMARLEM2494")
      clusNames[i] <- "Leishmania1"
    if (clusternames == "LmajorFriedlin" |
        clusternames == "LmajorLV39c5" |
        clusternames == "LmajorSD75" |
        clusternames == "LdonovaniBHU1220" |
        clusternames == "LdonovaniBPK282A1" |
        clusternames == "LinfantumJPCM5" |
        clusternames == "LturanicaLEM423" |
        clusternames == "LarabicaLEM1108" |
        clusternames == "LamazonensisMHOMBR71973M2269" |
        clusternames == "LmexicanaMHOMGT2001U1103" |
        clusternames == "LtropicaL590" |
        clusternames == "LaethiopicaL147" |
        clusternames == "LgerbilliLEM452")
    clusNames[i] <- "Leishmania2"
    if (clusternames == "LbraziliensisMHOMBR75M2904" |
        clusternames == "LbraziliensisMHOMBR75M2903" |
        clusternames == "LpanamensisMHOMPA94PSC1" |
        clusternames == "LpanamensisMHOMCOL81L13")
      clusNames[i] <- "Lvianna"
  }
  
  headers <- gsub("\\s", "", headers)
  functionalclasses <-
    lapply(X = headers, function(X)
      unlist(strsplit(X, split = "_"))[length(unlist(strsplit(X, split = "_")))])
  
  functionalclasses <- as.character(functionalclasses)
  
  tRNAdf <-
    data.frame(genomenames4, sequences, genomenames2, functionalclasses,clusNames)
  names(tRNAdf) <-
    c("genomename", "sequence", "headers", "funclass","clusNames")
  
  genome_list = split(tRNAdf, f = tRNAdf$clusNames)
  
  homo <-
    genome_list[names(genome_list) == "HOMO"][[1]]
  homo$funclass <-
    as.character(lapply(X = homo$headers , function(X)
      substr(X, 7, 7)))
  homo$funclass <- gsub("\\s", "", homo$funclass)
  homo$headers <- gsub("\\s", "", homo$headers)
  
  homoheader <- paste(homo$headers, homo$funclass, sep = "")
  genome_list[names(genome_list) == "HOMO"][[1]]$headers <-
    homoheader
  
  # make a file for each genome and put the genomes in there
  for (i in 1:length(genome_list)) {
    # write the DF in a file to be processed by the bash
    filepath <-
      paste(
        "/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfminput_2019/tRNADBClustered/",
        names(genome_list[i]),
        ".fasta",
        sep = ""
      )
    
    write.fwf(
      data.frame(genome_list[[i]]$headers,
                 genome_list[[i]]$sequence),
      filepath,
      sep = "\n",
      colnames = FALSE
    )
  }
  vidualization(genome_list)
   
}
# some statistics on TryTryp genes _________________________________________________________
vidualization <- function(genome_list) {
  # This function will vidualize:
  # 1. number of tRNA genes in HomoC and TryTryp genomes
  # 2. Percentage of 21 tRNA functional classes covered by each genome
  plotpath <-
    "/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfminput_2019/"
  tRNAcountdf <- data.frame(names(genome_list))
  # I used genomeorder <- tRNAcountdf$genome from vidualization of integrate_Te_Ara.R to be able to compare them
  #tRNAcountdf <- data.frame(genomeorder)
  tRNAcountdf$counts <- 0
  names(tRNAcountdf) <- c("genome", "tRNAcounts")
  for (i in 1:length(genome_list)) {
    #tRNAcountdf$genome[i] <- names(genome_list[i])
    tRNAcountdf[tolower(tRNAcountdf$genome) == tolower(names(genome_list[i])), ]$tRNAcounts <-
      nrow(genome_list[[i]])
  }
  
  tRNAcountdf <- tRNAcountdf[order(tRNAcountdf$tRNAcounts),]
  tRNAcountdf$genome <-
    factor(tRNAcountdf$genome, levels = tRNAcountdf$genome)
  library(ggplot2)
  p <- ggplot(data = tRNAcountdf, aes(x = genome, y = tRNAcounts)) +
    geom_bar(stat = "identity", fill = "steelblue") + geom_text(
      aes(label = tRNAcounts),
      vjust = 1.2,
      color = "white",
      position = position_dodge(0.9),
      size = 3.5
    ) + theme(axis.text.x = element_text(angle = 90, hjust = 1),
              plot.title = element_text(hjust = 0.5)) +
    labs(title =
           "Number of tRNA genes in HomoC and TryTryp clusters",
         x =
           "Genome", y = "Number of tRNAs")
  p
  ggsave(paste(plotpath, "tRNAcounts_clustered.png", sep = ""),
         width = 14,
         height = 7)
  
  
  # plot genomes vs percentage of 21 functional classes they contain
  
  tRNAcountdf$genome <-
    as.character(tRNAcountdf$genome, levels = tRNAcountdf$genome)
  
  func_classes <-
    c(
      "A",
      "R",
      "N" ,
      "D",
      "C" ,
      "Q",
      "E",
      "G",
      "H",
      "I",
      "L",
      "K",
      "M" ,
      "X",
      "F",
      "P",
      "S",
      "T",
      "W",
      "Y" ,
      "V"
    )
  tRNAcountdf$percent_21 <- 0
  tRNAcountdf$missing <- ''
  for (i in 1:length(genome_list)) {
    print(i)
    if (tolower(names(genome_list[i])) != "homo")
    {
      # remove function Z for now
      # isZ <- genome_list[[i]]$funclass == "Z"
      # genome_list[[i]] <- genome_list[[i]][!isZ,]
      # countsdf <-
      #   as.data.frame(table(as.character(genome_list[[i]]$funclass)))
      # fcount <- nrow(countsdf)
      # tRNAcountdf[tRNAcountdf$genome == names(genome_list[i]),]$percent_21 <-
      #   (fcount / 21) * 100
      # # make a vector of missing functional classes
      # tRNAcountdf[tRNAcountdf$genome == names(genome_list[i]),]$missing <-
      #   paste(setdiff(func_classes, as.character(countsdf$Var1)), collapse = " ")
      #
      #
      countsdf <-
        as.data.frame(table(as.character(genome_list[[i]]$funclass)))
      tRNAcountdf[tRNAcountdf$genome == names(genome_list[i]),]$missing <-
        paste(setdiff(func_classes, as.character(countsdf$Var1)), collapse = " ")
      diffset <-
        setdiff(func_classes, as.character(countsdf$Var1))
      fcount <- 21 - length(diffset)
      tRNAcountdf[tRNAcountdf$genome == names(genome_list[i]),]$percent_21 <-
        (fcount / 21) * 100
      # make a vector of missing functional classes
      
    }
    else
      tRNAcountdf[tolower(tRNAcountdf$genome) == tolower(names(genome_list[i])),]$percent_21 <-
        100
    
  }
  tRNAcountdf$genome <-
    factor(tRNAcountdf$genome, levels = tRNAcountdf$genome)
  p <- ggplot(data = tRNAcountdf, aes(x = genome, y = percent_21)) +
    geom_bar(stat = "identity", fill = "#56B4E9")  + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                                                           plot.title = element_text(hjust = 0.5)) +
    labs(title =
           "Percentage of 21 tRNA functional classes covered by each cluster",
         x =
           "Genome", y = "Percentage of tRNA functional classes(total=21)") + geom_text(
             aes(label =
                   missing),
             vjust = 0.5,
             angle = 90,
             hjust = -0.1
           )
  p
  ggsave(paste(plotpath, "funcPerc_clustered.png", sep = ""),
         width = 14,
         height = 7)
  
  # some statistics on TryTryp genes _________________________________________________________
  #
  # nonhomo <- (tRNAdf$genomename != "HOMO")
  # TryTrypfunccounts_df <-
  #   as.data.frame(table(as.character(tRNAdf$funclass[nonhomo])))
  # names(TryTrypfunccounts_df) <- c("func", "freq")
  # plot(TryTrypfunccounts_df$func, TryTrypfunccounts_df$freq)
  # barplot(table(as.character(tRNAdf$funclass[nonhomo])), xlab = "functional class", ylab = "frequency in TryTryp Genomes")
  #_____________________________________________________________________________________________
}