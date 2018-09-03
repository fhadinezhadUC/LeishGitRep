Tse1.3_VS_2_Comparison <- function() {
  # read the file files for two tse version output and make the tse dataframe for each of them
  library(IRanges)
  library(stringr)
  library(miceadds)
  library(GenomicRanges)
  library(gdata)
  library(VennDiagram)
  library(data.table)
  library(ggpubr)
  tse2_filenames <-
    list.files(
      "/home/fatemeh/Leishmania_Aug2018/genefindersScriptAndOutput/TseOutput/tse.out/",
      pattern = "*.tse.out",
      full.names = TRUE
    )
  tse2_ss_filenames <-
    list.files(
      "/home/fatemeh/Leishmania_Aug2018/genefindersScriptAndOutput/TseOutput/ss.tse.out/",
      pattern = "*.SS.tse.out",
      full.names = TRUE
    )
  
  tse1_filenames <- tse2_filenames
  tse1_ss_filenames <- tse2_filenames
  for (i in 1:length(tse2_filenames)) {
    namearray <- unlist(strsplit(tse2_filenames[i], split = "/"))
    nameorg <- namearray[length(namearray)]
    nameorg <- unlist(strsplit(nameorg, split = ".tse.out"))[1]
    
    tse1_filenames[i] = paste(
      "/home/fatemeh/Leishmania_Aug2018/genefindersScriptAndOutput/TseOutput/L/tse.out/",
      nameorg,
      ".tse.out",
      sep = ""
    )
    tse1_ss_filenames[i] = paste(
      "/home/fatemeh/Leishmania_Aug2018/genefindersScriptAndOutput/TseOutput/L/ss.tse.out/",
      nameorg,
      ".SS.tse.out",
      sep = ""
    )
  }
  
  tse2df <- making_tse_df(tse2_filenames[1], tse2_ss_filenames[1])
  tse1df <- making_tse_df(tse1_filenames[1], tse1_ss_filenames[1])
  # make the tsedf dataframe for each version
  for (i in 2:length(tse2_filenames)) {
    temp2 <- making_tse_df(tse2_filenames[i], tse2_ss_filenames[i])
    temp1 <- making_tse_df(tse1_filenames[i], tse1_ss_filenames[i])
    tse1df <- rbind(tse1df, temp1)
    tse2df <- rbind(tse2df, temp2)
  }
  
  
  visualize(tse1df, tse2df)
  
  
}

visualize <- function(tse1df, tse2df) {
  
  tse1df$tsesourceseq <- as.character(tse1df$tsesourceseq)
  tse2df$tsesourceseq <- as.character(tse2df$tsesourceseq)
  tse1df$tsebegin <- as.integer(as.character(tse1df$tsebegin))
  tse2df$tsebegin <- as.integer(as.character(tse2df$tsebegin))
  tse1df$tseend <- as.integer(as.character(tse1df$tseend))
  tse2df$tseend <- as.integer(as.character(tse2df$tseend))
  tse1df$tsedirection <- factor(tse1df$tsedirection)
  tse2df$tsedirection <- factor(tse2df$tsedirection)
  tse1df$tsedirection = as.character(tse1df$tsedirection)
  tse2df$tsedirection = as.character(tse2df$tsedirection)
  tse1df$tsesourceseq <- factor(tse1df$tsesourceseq)
  tse2df$tsesourceseq <- factor(tse2df$tsesourceseq)
  tse1df$tsesourceseq = as.character(tse1df$tsesourceseq)
  tse2df$tsesourceseq = as.character(tse2df$tsesourceseq)
  
  #################################################################################################################
  # venn diagram for 4 sets of genes from tse1 and tse2
  
  ispseudo1 <- tse1df$note == "pseudo"
  ispseudo2 <- tse2df$note == "pseudo"
  pseudodf1 <- tse1df[ispseudo1, ]
  pseudodf2 <- tse2df[ispseudo2, ]
  nonpseudodf1 <- tse1df[!ispseudo1, ]
  nonpseudodf2 <- tse2df[!ispseudo2, ]
  
  grid.newpage()
  draw.quad.venn(
    nrow(nonpseudodf1),
    nrow(nonpseudodf2),
    nrow(pseudodf1),
    nrow(pseudodf2),
    nrow(overlap_nonpseudodf),
    0,
    nrow(overlap_pseudo2_nonpseudo1df),
    nrow(overlap_pseudo1_nonpseudo2df),
    0,
    nrow(overlap_pseudodf),
    0,
    0,
    0,
    0,
    0
    ,
    category = c(
      "Nonpseudo Genes TSE1",
      "Nonpseudo Genes TSE2",
      "Pseudo Genes TSE1",
      "Pseudo Genes TSE2"
    ),
    lty = "blank",
    fill = c("skyblue", "pink1",
             "mediumorchid", "orange")
  )
  
  overlapdf <- findoverlapdf(tse1df, tse2df)
  
  #################################################################################################################
  # Corelation of gene scores found by two tse1 and tse2.
  # Group colored by being pseudo or non-pseudo
  
  tse1score <- as.integer(as.character(overlapdf$tsescore1))
  tse2score <- as.integer(as.character(overlapdf$tsescore2))
  my_data <- cbind(tse1score, tse2score)
  my_data <- as.data.frame(my_data)
  #http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r#what-is-correlation-test
  # Add one column to overlapdf for their color :
  # 4 colors for 4 cases:
  # Non pseudo in both, psudo in both,  pseudo1-nonpseudo2, nonpseudo1-pseudo2
  bothnonpseudo = overlapdf$note1 != "pseudo" &
    overlapdf$note2 != "pseudo"
  bothpseudo = overlapdf$note1 == "pseudo" &
    overlapdf$note2 == "pseudo"
  tse1pseudo = overlapdf$note1 == "pseudo" &
    overlapdf$note2 != "pseudo"
  tse2pseudo = overlapdf$note1 != "pseudo" &
    overlapdf$note2 == "pseudo"
  
  overlapdf$color <- "Both_nonpseudo"
  overlapdf[bothpseudo,]$color <- "Both_pseudo"
  overlapdf[tse1pseudo,]$color <- "Tse1_pseudo"
  overlapdf[tse2pseudo,]$color <- "Tse2_pseudo"
  my_data$color <- overlapdf$color
  ggplot(data = my_data, aes(x = tse1score, y = tse2score)) + geom_point(aes(colour = factor(color)),
                                                                         size = 2 ,
                                                                         show.legend = TRUE) +
    geom_smooth(
      data = my_data,
      method = "lm",
      se = FALSE,
      color = "grey"
    ) +  labs(colour = "Gene Type") +
    ggtitle("Corelation of gene covariance scores found by tRNAscan1.3 and tRNAscan2.0") +
    theme(plot.title = element_text(hjust = 0.5)) + xlab("TSE1.3 covariance score") + ylab("TSE2.0 covariance score")
  
}
findoverlapdf <- function(tse1df, tse2df) {
  tse1 <-
    GRanges(
      seqnames = tse1df$tsesourceseq,
      ranges = IRanges(tse1df$tsebegin, tse1df$tseend),
      strand = Rle(strand(tse1df$tsedirection))
    )
  tse2 <-
    GRanges(
      seqnames = tse2df$tsesourceseq,
      ranges = IRanges(tse2df$tsebegin, tse2df$tseend),
      strand = Rle(strand(tse2df$tsedirection))
    )
  
  overlaps <- as.data.frame(findOverlaps(tse1, tse2))
  Exactmatch <-
    as.data.frame(findOverlaps(tse1, tse2), type = "equal")
  
  sprintf("Total number of nonpseudo-genes found by TSE2.0 is: %d",
          nrow(tse2df))
  sprintf("total number of nonpseudo-genes found by TSE1.3 is: %d",
          nrow(tse1df))
  sprintf("total number of nonpseudo-genes found by both is: %d",
          nrow(overlaps))
  if (nrow(overlaps) == nrow(Exactmatch))
    print("genes found by both TSE1.3 and TSE2.0 have exact same coordinate.")
  
  # compare the score of the common found genes. using their coordinate
  
  names(overlaps) <- c("tse1", "tse2")
  m <- matrix(ncol = ncol(tse1df) + ncol(tse2df), nrow = 1)
  overlapdf <- as.data.frame(m)
  names(tse1df) <- paste(names(tse1df), "1", sep = "")
  names(tse2df) <- paste(names(tse2df), "2", sep = "")
  names(overlapdf) <- c(names(tse1df), names(tse2df))
  # remove the last one!
  overlapdf <- overlapdf[1:nrow(overlapdf) - 1 , ]
  
  for (i in 1:nrow(overlaps)) {
    tempbind <-
      as.data.frame(c(tse1df[overlaps$tse1[i], ], tse2df[overlaps$tse2[i], ]))
    overlapdf <- rbind(tempbind, overlapdf)
  }
  
  visualize_overlapdfs(overlapdf)
  #visualize_nonoverlapdfs(overlaps, tse1df, tse2df)
  overlapdf
  
}
making_tse_df <- function(tse_filename, tse_ss_filename) {
  # extracting the sourceOrg from filemame
  namearray <- unlist(strsplit(tse_filename, split = "/"))
  nameorg <- namearray[length(namearray)]
  nameorg <- unlist(strsplit(nameorg, split = "\\."))[1]
  
  tsegenename = ""
  tsesourceOrg =  nameorg
  tseidentity = ""
  tsedirection = ""
  tsebegin = ""
  tseend = ""
  tseac = ""
  tsesourceseq = ""
  tsescore = ""
  tsegeneseq = ""
  tsegeness = ""
  tseintronbegin = ""
  tseintronend = ""
  tseacloc = ""
  note = ""
  
  geneinfo <-
    data.frame(
      tsegenename,
      tsesourceOrg,
      tseidentity,
      tsedirection,
      tsebegin,
      tseend,
      tseac,
      tsesourceseq,
      tsescore,
      tsegeneseq,
      tsegeness,
      tseintronbegin,
      tseintronend,
      tseacloc,
      note
    )
  
  con = file(tse_filename, "r")
  con2 = file(tse_ss_filename, "r") # read 6 lines for each tRNA gene
  
  # # # skip the header
  readLines(con, n = 1)
  readLines(con, n = 1)
  readLines(con, n = 1)
  
  while (TRUE) {
    line = readLines(con, n = 1)
    if (length(line) == 0) {
      break
    }
    
    line_array <- unlist(strsplit(line, split = "\\s+"))
    tsesourceseq = line_array[1]
    tsegenename <-
      paste(tsesourceOrg, tsesourceseq, line_array[2], sep = "_")
    tseidentity = line_array[5]
    if (as.integer(line_array[3]) > as.integer(line_array[4]))
    {
      tsedirection = "-"
    }
    else
    {
      tsedirection = "+"
    }
    tsebegin = line_array[3]
    tseend = line_array[4]
    tseac = line_array[6]
    tsescore = line_array[9]
    tseintronbegin = line_array[7]
    tseintronend = line_array[8]
    if (!is.na(line_array[10]))
      note = line_array[10]
    else
      note = "notfound"
    
    # read the first two lines
    ss_line1 = readLines(con2, n = 1)
    ss_line2 = readLines(con2, n = 1)
    tseacloc = unlist(strsplit(ss_line2, split = "\\s+"))[6]
    l = "."
    while (TRUE) {
      l = readLines(con2, n = 1)
      if (l == "")
        break
      if (unlist(strsplit(l, split = "\\s+"))[1] == "Seq:")
        tsegeneseq <- unlist(strsplit(l, split = "\\s+"))[2]
      if (unlist(strsplit(l, split = "\\s+"))[1] == "Str:")
        tsegeness <- unlist(strsplit(l, split = "\\s+"))[2]
    }
    
    tempdf <-
      data.frame(
        tsegenename,
        tsesourceOrg,
        tseidentity,
        tsedirection,
        tsebegin,
        tseend,
        tseac,
        tsesourceseq,
        tsescore,
        tsegeneseq,
        tsegeness,
        tseintronbegin,
        tseintronend,
        tseacloc,
        note
      )
    
    geneinfo <- rbind(geneinfo, tempdf)
  }
  close(con)
  close(con2)
  geneinfo <- geneinfo[2:nrow(geneinfo),]
  
  geneinfo$tsebegin <- as.integer(as.character(geneinfo$tsebegin))
  geneinfo$tseend <- as.integer(as.character(geneinfo$tseend))
  bigger <- geneinfo$tsebegin > geneinfo$tseend
  for (i in 1:length(bigger)) {
    if (bigger[i])
    {
      temp <- geneinfo$tsebegin[i]
      geneinfo$tsebegin[i] <- geneinfo$tseend[i]
      geneinfo$tseend[i] <- temp
    }
  }
  geneinfo
}

# visualize_nonoverlapdfs <- function(overlaps, tse1df, tse2df) {
#   diff_tse1_tse2 <-
#     setdiff(seq(1, nrow(tse1df), 1), overlaps$tse1)
#   diff_tse2_tse1 <-
#     setdiff(seq(1, nrow(tse2df), 1), overlaps$tse2)
#
#   if (length(diff_tse1_tse2) > 0)
#     diff_tse1_tse2_df <- tse1df[diff_tse1_tse2, ]
#   if (length(diff_tse2_tse1) > 0)
#     diff_tse2_tse1_df <- tse2df[diff_tse2_tse1, ]
#
#   # vertical axis is score, X axis is genefinder tse1, tse2
#   # make the 2 column dataframe score genefnder
#   tse1data <-
#     data.frame(diff_tse1_tse2_df$tsescore1, rep("tse1", length(diff_tse1_tse2_df$tsescore1)))
#   tse2data <-
#     data.frame(diff_tse2_tse1_df$tsescore2, rep("tse2", length(diff_tse2_tse1_df$tsescore2)))
#   names(tse1data) <- c("score", "genefinder")
#   names(tse2data) <- c("score", "genefinder")
#   boxplotdata <- rbind(tse1data, tse2data)
#   boxplot(as.integer(as.character(score)) ~ as.character(genefinder), data = boxplotdata,
#           main = "Gene-scores Found by One of tRNAscan2.0 or tRNAscan1.3")
#
# }
