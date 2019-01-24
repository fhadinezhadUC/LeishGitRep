split_tRNAgene_with_genomename <- function() {
  # This function reads in the .fasta file containg aligned tRNA genes of HomoC and TryTryp genomes
  # Extract the functional classes from Headers
  # Makes a dataframe of "genomename", "sequence", "headers", "funclass"
  # Splits the dataframe into a list of data frames based on genomename
  # Writes each dataframe as a fasta file in a file with genome's names
  # At the end it call function vidualization passing the list of data frames
  
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
# some statistics on TryTryp genes _________________________________________________________
vidualization <- function(genome_list) {
  # This function will vidualize:
  # 1. number of tRNA genes in HomoC and TryTryp genomes
  # 2. Percentage of 21 tRNA functional classes covered by each genome
  plotpath <-
    "/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfminput_2019/"
  tRNAcountdf <- data.frame(names(genome_list))
  tRNAcountdf$counts <- 0
  names(tRNAcountdf) <- c("genome", "tRNAcounts")
  for (i in 1:length(genome_list)) {
    tRNAcountdf$genome[i] <- names(genome_list[i])
    tRNAcountdf$tRNAcounts[i] <- nrow(genome_list[[i]])
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
           "Number of tRNA genes in HomoC and TryTryp genomes",
         x =
           "Genome", y = "Number of tRNAs")
  p
  ggsave(paste(plotpath, "tRNAcounts.png", sep = ""),
         width = 14,
         height = 7)
  
  
  # plot genomes vs percentage of 21 functional classes they contain
  
  tRNAcountdf$genome <-
    as.character(tRNAcountdf$genome, levels = tRNAcountdf$genome)
  
  tRNAcountdf$percent_21 <- 0
  for (i in 1:length(genome_list)) {
    print(i)
    if (names(genome_list[i]) != "HOMO")
    {
      # remove function Z for now
      isZ <- genome_list[[i]]$funclass=="Z"
      genome_list[[i]] <- genome_list[[i]][!isZ,]
      countsdf <-
        as.data.frame(table(as.character(genome_list[[i]]$funclass)))
      fcount <- nrow(countsdf)
      tRNAcountdf[tRNAcountdf$genome == names(genome_list[i]),]$percent_21 <-
        (fcount / 21) * 100
    }
    else
      tRNAcountdf[tRNAcountdf$genome == names(genome_list[i]),]$percent_21 <-
        100
    
  }
  tRNAcountdf$genome <-
    factor(tRNAcountdf$genome, levels = tRNAcountdf$genome)
  p <- ggplot(data = tRNAcountdf, aes(x = genome, y = percent_21)) +
    geom_bar(stat = "identity", fill = "#56B4E9")  + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                                                           plot.title = element_text(hjust = 0.5)) +
    labs(title =
           "Percentage of 21 tRNA functional classes covered by each genome",
         x =
           "Genome", y = "Percentage of tRNA functional classes(total=21)")
  p
  ggsave(paste(plotpath, "funcPerc.png", sep = ""),
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