Integrate_Tse_Ara <- function() {
  # check for this ADWP02026113 	1	9   	76  	Undet	NNN	0	0	22.0	pseudo
  #http://www.bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.pdf
  #https://support.bioconductor.org/p/56123/
  library(IRanges)
  library(stringr)
  library(miceadds)
  library(GenomicRanges)
  library(gdata)
  # read ara
  ara_filenames <-
    list.files(
      "/home/fatemeh/Leishmania_Aug2018/genefindersScriptAndOutput/AraOutput/AraOutput_i30_ps95_rp100/",
      pattern = "*.ara.out",
      full.names = TRUE
    )
  #just take the ara_filenames and make the tse file for it to read.
  tse_filenames <- ara_filenames
  tse_ss_filenames <- ara_filenames
  for (i in 1:length(ara_filenames)) {
    namearray <- unlist(strsplit(ara_filenames[i], split = "/"))
    nameorg <- namearray[length(namearray)]
    nameorg <- unlist(strsplit(nameorg, split = ".ara.out"))[1]
    
    tse_filenames[i] = paste(
      "/home/fatemeh/Leishmania_Aug2018/genefindersScriptAndOutput/TseOutput/tse.out/",
      nameorg,
      ".tse.out",
      sep = ""
    )
    tse_ss_filenames[i] = paste(
      "/home/fatemeh/Leishmania_Aug2018/genefindersScriptAndOutput/TseOutput/ss.tse.out/",
      nameorg,
      ".SS.tse.out",
      sep = ""
    )
  }
  
  
  numberofgenes = 0
  
  aradf <- making_ara_df(ara_filenames[1])
  tsedf <-
    making_tse_df(tse_filenames[1], tse_ss_filenames[1])
  integrated_tse_ara <- integrate(aradf, tsedf)
  
  for (i in 2:length(ara_filenames)) {
    #length(ara_filenames)
    
    aradf <- making_ara_df(ara_filenames[i])
    tsedf <-
      making_tse_df(tse_filenames[i], tse_ss_filenames[i])
    temp_integrated <- integrate(aradf, tsedf)
    integrated_tse_ara <-
      rbind(integrated_tse_ara, temp_integrated)
    numberofgenes = nrow(integrated_tse_ara)
    # every round we need to rbind the old df with the new one
  }
  print(numberofgenes)
  formatoutput(integrated_tse_ara)
  
  
}

formatoutput <- function(integrated_tse_ara) {
  # writing the dataframe in two files :
  # file1 for coordinates: geneid, sourceOrg, sourceseq, sourceSO, direction, tse/aracoordinate /
  # file2for SS:          geneid:\n tse: tsegeneseq \n tseSS/ arageneseq, araSS /
  # file3 for identities:  geneid, tse/araidentity, tse/araac, tseacloc, arascore/tsescore
  # file4 for intron:      geneid, tseintroncoordinate
  
  n <- data.frame("geneid",
                  "tseintronbegin",
                  "tseintronend")
  names(n) <- c("geneid",
                "tseintronbegin",
                "tseintronend")
  introndf <-  data.frame(
    integrated_tse_ara$geneid,
    integrated_tse_ara$tseintronbegin,
    integrated_tse_ara$tseintronend
  )
  names(introndf) <- names(n)
  write.fwf(
    rbind(n, introndf),
    width = c(60, 20, 20),
    colnames = FALSE,
    file = "./introns.txt"
  )
  
  
  n <- data.frame(
    "geneid",
    "sourceOrg",
    "sourceseq",
    "sourceSO",
    "direction",
    "tsebegin",
    "tseend",
    "arabegin",
    "araend"
  )
  
  names(n) <- c(
    "geneid",
    "sourceOrg",
    "sourceseq",
    "sourceSO",
    "direction",
    "tsebegin",
    "tseend",
    "arabegin",
    "araend"
  )
  coordinatedf <- data.frame(
    integrated_tse_ara$geneid,
    integrated_tse_ara$sourceOrg,
    integrated_tse_ara$sourceseq,
    integrated_tse_ara$sourceSO,
    integrated_tse_ara$direction,
    integrated_tse_ara$tsebegin,
    integrated_tse_ara$tseend,
    integrated_tse_ara$arabegin,
    integrated_tse_ara$araend
  )
  names(coordinatedf) <- names(n)
  write.fwf(
    rbind(n, coordinatedf),
    colnames = FALSE,
    width = c(60, 35, 35, 13, 10, 10, 10, 10, 10),
    file = "./coordinates.txt"
  )
  
  
  n <-
    data.frame(
      "geneid",
      "tseidentity",
      "araidentity",
      "tseac",
      "araac",
      "tseacloc",
      "arascore",
      "tsescore",
      "foundby"
    )
  names(n) <-
    c(
      "geneid",
      "tseidentity",
      "araidentity",
      "tseac",
      "araac",
      "tseacloc",
      "arascore",
      "tsescore",
      "foundby"
    )
  identitydf <- data.frame(
    integrated_tse_ara$geneid,
    integrated_tse_ara$tseidentity,
    integrated_tse_ara$araidentity,
    integrated_tse_ara$tseac,
    integrated_tse_ara$araac,
    integrated_tse_ara$tseacloc,
    integrated_tse_ara$arascore,
    integrated_tse_ara$tsescore,
    integrated_tse_ara$foundby
  )
  names(identitydf) <- names(n)
  write.fwf(
    rbind(n, identitydf),
    width = c(60, 12, 12, 15, 15, 12, 12, 12, 10),
    file = "./identities.txt",
    colnames = FALSE
  )
  
  write.fwf(
    data.frame(
      paste("geneid: ", integrated_tse_ara$geneid, sep = ""),
      paste("tseseq: ", integrated_tse_ara$tsegeneseq, sep = ""),
      paste("tsess:  ", integrated_tse_ara$tsegeness, sep = ""),
      paste("araseq: ", integrated_tse_ara$arageneseq, sep = ""),
      paste("arass:  ", integrated_tse_ara$arageness, sep = ""),
      character(length = length(integrated_tse_ara$arageness))
    ),
    file = "./secondaryS.txt",
    sep = "\n",
    colnames = FALSE
  )
  
  
  
}

integrate <- function(aradf, tsedf) {
  aradf$arabegin <- as.integer(as.character(aradf$arabegin))
  aradf$araend <- as.integer(as.character(aradf$araend))
  tsedf$tsebegin <- as.integer(as.character(tsedf$tsebegin))
  tsedf$tseend <- as.integer(as.character(tsedf$tseend))
  tsedf$tsesourceseq <- as.character(tsedf$tsesourceseq)
  aradf$arasourceseq <- as.character(aradf$arasourceseq)
  
  # for tse output if it is on reverse strand replace the begin and end
  bigger <- tsedf$tsebegin > tsedf$tseend
  for (i in 1:length(bigger)) {
    if (bigger[i])
    {
      temp <- tsedf$tsebegin[i]
      tsedf$tsebegin[i] <- tsedf$tseend[i]
      tsedf$tseend[i] <- temp
    }
  }
  
  aradf$aradirection = as.character(aradf$aradirection)
  tsedf$tsedirection = as.character(tsedf$tsedirection)
  aradf$arasourceseq <- as.character(aradf$arasourceseq)
  tsedf$tsesourceseq <- as.character(tsedf$tsesourceseq)
  
  ara <-
    GRanges(
      seqnames = aradf$arasourceseq ,
      ranges = IRanges(aradf$arabegin, aradf$araend),
      strand = Rle(strand(aradf$aradirection))
    )
  
  tse <-
    GRanges(
      seqnames = tsedf$tsesourceseq,
      ranges = IRanges(tsedf$tsebegin, tsedf$tseend),
      strand = Rle(strand(tsedf$tsedirection))
    )
  overlapps <- as.data.frame(findOverlaps(ara, tse))
  names(overlapps) <- c("ararecord", "tserecord")
  # making a dataframe for binding two dataframes
  m <- matrix(ncol = ncol(aradf) + ncol(tsedf), nrow = 1)
  overlapdf <- as.data.frame(m)
  names(overlapdf) <- c(names(aradf), names(tsedf))
  for (i in 1:nrow(overlapps)) {
    tempbind <-
      as.data.frame(c(aradf[overlapps$ararecord[i],], tsedf[overlapps$tserecord[i],]))
    overlapdf <- rbind(tempbind, overlapdf)
  }
  integrated_gene_file <- overlapdf
  integrated_gene_file <-
    integrated_gene_file[1:(nrow(integrated_gene_file) - 1), ]
  integrated_gene_file$foundby <- "both"
  
  # dealing with those that do not have overlapps
  diff_ara_tse <-
    setdiff(seq(1, nrow(aradf), 1), overlapps$ararecord)
  diff_tse_ara <-
    setdiff(seq(1, nrow(tsedf), 1), overlapps$tserecord)
  
  if (length(diff_ara_tse) > 0)
    for (i in 1:length(diff_ara_tse)) {
      tempbind <-
        as.data.frame(c(aradf[diff_ara_tse[i],], rep("notfound", ncol(tsedf)), "ara"))
      names(tempbind) <- names(integrated_gene_file)
      tempbind$tsebegin = -1
      tempbind$tseend = -1
      integrated_gene_file <- rbind(tempbind, integrated_gene_file)
    }
  if (length(diff_tse_ara) > 0)
    for (i in 1:length(diff_tse_ara)) {
      tempbind <-
        as.data.frame(c(rep("notfound", ncol(aradf)), tsedf[diff_tse_ara[i],], "tse"))
      names(tempbind) <- names(integrated_gene_file)
      tempbind$arabegin = -1
      tempbind$araend = -1
      integrated_gene_file <- rbind(tempbind, integrated_gene_file)
    }
  
  # now we make integrated_tse_ara with unified geneid for tse and ara
  m <-
    matrix(nrow = nrow(integrated_gene_file),
           ncol = 27)
  integrated_tse_ara <- as.data.frame(m)
  names(integrated_tse_ara) <- c(
    "geneid",
    "sourceOrg",
    "tseidentity",
    "araidentity",
    "direction",
    "tsebegin",
    "tseend",
    "arabegin",
    "araend",
    "tseac",
    "araac",
    "sourceseq",
    "tsescore",
    "arascore",
    "tsegeneseq",
    "arageneseq",
    "tsegeness",
    "arageness",
    "tseintronbegin",
    "tseintronend",
    "araintronbegin",
    ####
    "araintronend",
    #####
    "tseacloc",
    "sourceSO",
    "tsenote",
    ####
    "aranote",
    #####
    "foundby"
  )
  
  
  # make the geneid from sourceOrg+sourceseq+index within that sourceOrg(which is the rownumber of our dataframe)
  
  for (i in 1:ncol(integrated_gene_file)) {
    integrated_gene_file[, i] <- as.character(integrated_gene_file[, i])
  }
  
  settse <-
    (integrated_gene_file$foundby == "tse" |
       integrated_gene_file$foundby == "both")
  integrated_tse_ara$sourceOrg[settse] = integrated_gene_file$tsesourceOrg[settse]
  integrated_tse_ara$sourceseq[settse] = integrated_gene_file$tsesourceseq[settse]
  integrated_tse_ara$direction[settse] = integrated_gene_file$tsedirection[settse]
  setara <- (integrated_gene_file$foundby == "ara")
  integrated_tse_ara$sourceOrg[setara] =  integrated_gene_file$arasourceOrg[setara]
  integrated_tse_ara$sourceseq[setara] = integrated_gene_file$arasourceseq[setara]
  integrated_tse_ara$direction[setara] = integrated_gene_file$aradirection[setara]
  
  integrated_tse_ara$tseidentity = integrated_gene_file$tseidentity
  integrated_tse_ara$araidentity = integrated_gene_file$araidentity
  integrated_tse_ara$tsebegin = integrated_gene_file$tsebegin
  integrated_tse_ara$tseend = integrated_gene_file$tseend
  integrated_tse_ara$arabegin = integrated_gene_file$arabegin
  integrated_tse_ara$araend = integrated_gene_file$araend
  integrated_tse_ara$tseac = integrated_gene_file$tseac
  integrated_tse_ara$araac = integrated_gene_file$araac
  integrated_tse_ara$tsescore = integrated_gene_file$tsescore
  integrated_tse_ara$arascore = integrated_gene_file$arascore
  integrated_tse_ara$tsegeneseq = integrated_gene_file$tsegeneseq
  integrated_tse_ara$arageneseq = integrated_gene_file$arageneseq
  integrated_tse_ara$tsegeness = integrated_gene_file$tsegeness
  integrated_tse_ara$arageness = integrated_gene_file$arageness
  integrated_tse_ara$tseintronbegin = integrated_gene_file$tseintronbegin
  integrated_tse_ara$tseintronend = integrated_gene_file$tseintronend
  integrated_tse_ara$tseacloc = integrated_gene_file$tseacloc
  integrated_tse_ara$sourceSO = integrated_gene_file$arasourceSO
  integrated_tse_ara$foundby = integrated_gene_file$foundby
  integrated_tse_ara$araintronbegin = integrated_gene_file$araintronbegin
  integrated_tse_ara$araintronend = integrated_gene_file$araintronend
  integrated_tse_ara$tsenote = integrated_gene_file$note
  integrated_tse_ara$aranote = integrated_gene_file$aranote
  
  integrated_tse_ara <- integrated_tse_ara[order(
    integrated_tse_ara$sourceOrg,
    integrated_tse_ara$sourceseq,
    integrated_tse_ara$arabegin
  ), ]
  
  for (i in 1:nrow(integrated_tse_ara)) {
    integrated_tse_ara$geneid[i] = paste(integrated_tse_ara$sourceOrg[i],
                                         integrated_tse_ara$sourceseq[i],
                                         i,
                                         sep = "_")
  }
  
  integrated_tse_ara
  # change the sourceorg for ara and make it from file name
  
}
making_ara_df <- function(arafilename) {
  namearray <- unlist(strsplit(arafilename, split = "/"))
  nameorg <- namearray[length(namearray)]
  nameorg <- unlist(strsplit(nameorg, split = "\\."))[1]
  
  aragenename = ""
  arasourceOrg = nameorg
  araidentity = ""
  aradirection = ""
  arabegin = ""
  araend = ""
  araintronbegin = ""
  araintronend = ""
  araac = ""
  arasourceSO = ""
  arasourceseq = ""
  arascore = ""
  arageneseq = ""
  arageness = ""
  aranote = ""
  
  geneinfo <-
    data.frame(
      aragenename ,
      arasourceOrg ,
      araidentity ,
      aradirection ,
      arabegin ,
      araend ,
      araintronbegin,
      araintronend,
      araac ,
      arasourceSO ,
      arasourceseq ,
      arascore ,
      arageneseq ,
      arageness,
      aranote
    )
  
  con = file(arafilename, "r")
  tmrinaline = ""
  checktmrna = ""
  while (TRUE) {
    #for (x in 1:54) {
    aranote = ""
    line1 = readLines(con, n = 1)
    
    #read two lines each loop
    if (length(line1) == 0) {
      break
    }
    line2 = readLines(con, n = 1)
    if (length(line2) == 0) {
      break
    }
    # if line1 starts with > and line2 does not start with 0
    temp1 <- grep("^>+", line1, value = TRUE)
    temp2 <- grep("^0+", line2, value = TRUE)
    if (length(temp1) > 0 & length(temp2) == 0)
    {
      # read the first element of line2 to see how many genes are found for this sequence (count)
      # read the next count * 3 lines
      # extract the sourceOrg, sourceseq, sourceSO
      line1array <-
        unlist(strsplit(line1, split = " | ", fixed = TRUE))
      arasourceseq <- substring(line1array[1], 2)
      #arasourceOrg <- substring(line1array[2], 10)
      arasourceSO <- substring(line1array[5], 4)
      
      genecount <-
        as.integer(unlist(strsplit(
          line2, split = " ", fixed = TRUE
        ))[1])
      for (j in 1:genecount) {
        #making the geneID with sourceorg+sourceseq+genenum
        aragenename <-
          paste(arasourceOrg, arasourceseq, j, sep = "_")
        
        trna_line1 <- readLines(con, n = 1)
        trna_line2 <- readLines(con, n = 1)
        trna_line3 <- readLines(con, n = 1)
        
        arageneseq <- trna_line2
        arageness <- trna_line3
        
        trna_line1_array <-
          unlist(strsplit(trna_line1, split = "\\s+"))
        
        #if we had ")i(" then we have the location for intron
        araactemp <- trna_line1_array[6]
        if (length(grep("+[)]i[(]+", araactemp)) != 0)
        {
          araac <- unlist(strsplit(araactemp, split = "i"))[1]
          intr <-
            unlist(strsplit(unlist(
              strsplit(araactemp, split = "i")
            )[2], split = ","))
          araintronbegin <- substring(intr[1], 2)
          araintronend <-
            as.integer(as.character(araintronbegin)) + as.integer(as.character(substring(intr[2], 1, nchar(intr[2]) -
                                                                                           1)))
        }
        else
        {
          araac <- araactemp
          araintronbegin = "nointron"
          araintronend = "nointron"
        }
        araidentity <- substring(trna_line1_array[2], 6)
        if (substring(araidentity, nchar(araidentity)) == "*")
        {
          aranote <- "pseudo"
          araidentity <-
            substring(trna_line1_array[2], 6, nchar(trna_line1_array[2]) - 1)
        }
        else
          aranote = ""
        coordinate <-
          unlist(str_extract_all(trna_line1_array[3], "\\d+"))
        arabegin <- coordinate[1]
        araend <- coordinate[2]
        firstchar = substring(trna_line1_array[3], 1, 1)
        if (firstchar == "c")
        {
          aradirection <- "-"
        }
        else
        {
          aradirection <- "+"
        }
        
        arascore <- trna_line1_array[4]
        
        # make one record in dataframe
        tempdf <- data.frame(
          aragenename ,
          arasourceOrg ,
          araidentity ,
          aradirection ,
          arabegin ,
          araend ,
          araintronbegin,
          araintronend,
          araac,
          arasourceSO ,
          arasourceseq ,
          arascore ,
          arageneseq ,
          arageness,
          aranote
        )
        geneinfo$araintronend <-
          as.character(geneinfo$araintronend)
        geneinfo <- rbind(geneinfo, tempdf)
        
      }
    }
  }
  close(con)
  geneinfo <- geneinfo[2:nrow(geneinfo),]
  geneinfo
}

# tse also reports introns. a possible intron
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
  geneinfo <- geneinfo[2:nrow(geneinfo), ]
  geneinfo
}

summery_table() <- function() {
  both <- integrated_tse_ara$foundby == "both"
  tseonly <- integrated_tse_ara$foundby == "tse"
  araonly <- integrated_tse_ara$foundby == "ara"
  
  # genes found by tse
  istse <- (tseonly | both)
  isara <- (araonly | both)
  isintersect_tse_ara <- both
  
  tse <- integrated_tse_ara[istse, ]
  ara <- integrated_tse_ara[isara, ]
  intersect_tse_ara <- integrated_tse_ara[isintersect_tse_ara, ]
  union_tse_ara <- integrated_tse_ara
  
  # making a 4 * 32 size dataframe
  m <- matrix(nrow = 4,
              ncol = 33,
              data = 0)
  summerytable <- as.data.frame(m)
  names(summerytable) <-
    c(
      "GeneSet",
      "#tRNA",
      "#nucleotides",
      "N/T",
      "lengthrange",
      "%G",
      "%C",
      "%T",
      "%A",
      "%intronContaining",
      "A",
      "C",
      "D",
      "E",
      "F",
      "G",
      "H",
      "I",
      "K",
      "L",
      "M",
      "N",
      "P",
      "Q",
      "R",
      "S",
      "T",
      "V",
      "W",
      "Y",
      "X",
      "Z",
      "$"
    )
  summerytable[, 1] <- c("TSE2", "ARA", "UNION", "INTERSECTION")
  summerytable[, 2] <-
    c(nrow(tse),
      nrow(ara),
      nrow(union_tse_ara),
      nrow(intersect_tse_ara))
  
  # sum of the gene length for union
  unionnucleotides <-
    sum(nchar(tse$tsegeneseq)) + sum(nchar(ara$arageneseq)) - sum(nchar(intersect_tse_ara$tsegeneseq))
  summerytable[, 3] <-
    c(sum(nchar(tse$tsegeneseq)) , sum(nchar(ara$arageneseq)) , unionnucleotides , sum(nchar(intersect_tse_ara$tsegeneseq)))
  summerytable[, 4] <-
    c(
      summerytable[1, 3] / summerytable[1, 2],
      summerytable[2, 3] / summerytable[2, 2],
      summerytable[3, 3] / summerytable[3, 2],
      summerytable[4, 3] / summerytable[4, 2]
    )
  tserange <-
    paste(min(nchar(tse$tsegeneseq)) , max(nchar(tse$tsegeneseq)) , sep = "-")
  ararange <-
    paste(min(nchar(ara$arageneseq)) , max(nchar(ara$arageneseq)) , sep = "-")
  unionrange <-
    paste(min(min(nchar(tse$tsegeneseq)), min(nchar(ara$arageneseq))) ,  max(max(nchar(tse$tsegeneseq)), max(nchar(ara$arageneseq))) , sep = "-")
  intersectionrange <-
    paste(min(nchar(intersect_tse_ara$tsegeneseq)) , max(nchar(intersect_tse_ara$tsegeneseq)) , sep =  "-")
  summerytable[, 5] <-
    c(tserange, ararange, unionrange, intersectionrange)
  
  # A T C G percentage
  tse$tsegeneseq <- tolower(tse$tsegeneseq)
  count <-
    sapply(c("a", "g", "c", "t"), function(nuc)
      str_count(tse$tsegeneseq, fixed(nuc)))
  countdf <- as.data.frame(count)
  A <- sum(countdf$a)
  G <- sum(countdf$g)
  C <- sum(countdf$c)
  TT <- sum(countdf$t)
  total <- sum(A, G, C, TT)
  summerytable$`%A`[1] <- A * 100 / total
  summerytable$`%G`[1] <- G * 100 / total
  summerytable$`%C`[1] <- C * 100 / total
  summerytable$`%T`[1] <- TT * 100 / total
  
  ara$arageneseq <- tolower(ara$arageneseq)
  count <-
    sapply(c("a", "g", "c", "t"), function(nuc)
      str_count(as.character(ara$arageneseq), fixed(nuc)))
  countdf <- as.data.frame(count)
  A <- sum(countdf$a)
  G <- sum(countdf$g)
  C <- sum(countdf$c)
  TT <- sum(countdf$t)
  total <- sum(A, G, C, TT)
  summerytable$`%A`[2] <- A * 100 / total
  summerytable$`%G`[2] <- G * 100 / total
  summerytable$`%C`[2] <- C * 100 / total
  summerytable$`%T`[2] <- TT * 100 / total
  
  intersect_tse_ara$tsegeneseq <-
    tolower(intersect_tse_ara$tsegeneseq)
  count <-
    sapply(c("a", "g", "c", "t"), function(nuc)
      str_count(as.character(intersect_tse_ara$tsegeneseq), fixed(nuc)))
  countdf <- as.data.frame(count)
  A <- sum(countdf$a)
  G <- sum(countdf$g)
  C <- sum(countdf$c)
  TT <- sum(countdf$t)
  total <- sum(A, G, C, TT)
  summerytable$`%A`[4] <- A * 100 / total
  summerytable$`%G`[4] <- G * 100 / total
  summerytable$`%C`[4] <- C * 100 / total
  summerytable$`%T`[4] <- TT * 100 / total
  
  seq <- integrated_tse_ara$tsegeneseq
  
  for (i in 1:length(seq)) {
    if (integrated_tse_ara$foundby[i] == "ara")
      seq[i] <- integrated_tse_ara$arageneseq[i]
  }
  
  seq <- tolower(seq)
  count <-
    sapply(c("a", "g", "c", "t"), function(nuc)
      str_count(as.character(seq), fixed(nuc)))
  countdf <- as.data.frame(count)
  A <- sum(countdf$a)
  G <- sum(countdf$g)
  C <- sum(countdf$c)
  TT <- sum(countdf$t)
  total <- sum(A, G, C, TT)
  summerytable$`%A`[3] <- A * 100 / total
  summerytable$`%G`[3] <- G * 100 / total
  summerytable$`%C`[3] <- C * 100 / total
  summerytable$`%T`[3] <- TT * 100 / total
  
  #_____________________ intron containing% ________________________________
  
  # araintronbegin is the location of the intron in the tRNA gene not genome!
  arahasintron <- ara$araintronbegin != "nointron"
  nrow(ara[arahasintron, ])
  
  tsehasintron <-
    (tse$tseintronbegin != "0" & tse$tseintronbegin != "notfound")
  nrow(tse[tsehasintron, ])
  
  intersecthasintron <-
    ((intersect_tse_ara$tseintronbegin == "0") &
       (intersect_tse_ara$araintronbegin == "nointron")
    )
  nrow(intersect_tse_ara[!intersecthasintron, ])
  
  
  unionhasintron <-
    (((union_tse_ara$tseintronbegin == "0") |
        (union_tse_ara$tseintronbegin == "notfound")
    ) &
      (union_tse_ara$araintronbegin == "nointron"))
  nrow(union_tse_ara[!unionhasintron,])
  
  summerytable[, 10] <-
    c(
      nrow(tse[tsehasintron,]) * 100 / nrow(tse),
      nrow(ara[arahasintron,]) * 100 / nrow(ara),
      nrow(union_tse_ara[!unionhasintron,]) * 100 / nrow(union_tse_ara),
      nrow(intersect_tse_ara[!intersecthasintron, ]) * 100 / nrow(intersect_tse_ara)
    )
  
  #______________________ pseudo frequency ______________________________________
  
  tsepseudo <- tse[tse$tsenote == "pseudo", ]
  arapseudo <- ara[ara$aranote == "pseudo", ]
  intersectpseudo <-
    intersect_tse_ara[intersect_tse_ara$tsenote == "pseudo" |
                        intersect_tse_ara$aranote == "pseudo",]
  unionpseudo <-
    union_tse_ara[union_tse_ara$tsenote == "pseudo" |
                    union_tse_ara$aranote == "pseudo",]
  summerytable[, 33] <-
    c(nrow(tsepseudo),
      nrow(arapseudo),
      nrow(unionpseudo),
      nrow(intersectpseudo))
  
  
  
  #___________________ Aminoacid frequency ________________________________________
  aminoacid <- character(length = nrow(tse))
  # frequency 23 classes
  for (i in 1:nrow(tse)) {
    id <- tolower(tse$tseidentity[i])
    
    if (id == "ala")
      aminoacid[i] <- "A"
    if (id == "arg")
      aminoacid[i] <- "R"
    if (id == "asn")
      aminoacid[i] <- "N"
    if (id == "asp")
      aminoacid[i] <- "D"
    if (id == "cys")
      aminoacid[i] <- "C"
    if (id == "gln")
      aminoacid[i] <- "Q"
    if (id == "glu")
      aminoacid[i] <- "E"
    if (id == "gly")
      aminoacid[i] <- "G"
    if (id == "his")
      aminoacid[i] <- "H"
    if (id == "ile")
      aminoacid[i] <- "I"
    if (id == "start")
      aminoacid[i] <- "start"
    if (id == "leu")
      aminoacid[i] <- "L"
    if (id == "lys")
      aminoacid[i] <- "K"
    if (id == "met")
      aminoacid[i] <- "M"
    if (id == "phe")
      aminoacid[i] <- "F"
    if (id == "pro")
      aminoacid[i] <- "P"
    if (id == "ser")
      aminoacid[i] <- "S"
    if (id == "thr")
      aminoacid[i] <- "T"
    if (id == "trp")
      aminoacid[i] <- "W"
    if (id == "tyr")
      aminoacid[i] <-
        "Y"
    if (id == "val")
      aminoacid[i] <-
        "V"
    if (id == "sec")
      aminoacid[i] <-
        "Z"
    if (id == "stop")
      aminoacid[i] <-
        "stop"
    
  }
  
  tseaminoacid <- data.frame(table(aminoacid))
  summerytable$A[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "A",]$Freq
  summerytable$C[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "C",]$Freq
  summerytable$D[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "D",]$Freq
  summerytable$E[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "E",]$Freq
  summerytable$F[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "F",]$Freq
  summerytable$G[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "G",]$Freq
  summerytable$H[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "H",]$Freq
  summerytable$I[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "I",]$Freq
  summerytable$K[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "K",]$Freq
  summerytable$L[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "L",]$Freq
  summerytable$M[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "M",]$Freq
  summerytable$N[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "N",]$Freq
  summerytable$P[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "P",]$Freq
  summerytable$Q[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "Q",]$Freq
  summerytable$R[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "R",]$Freq
  summerytable$S[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "S",]$Freq
  summerytable$T[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "T",]$Freq
  summerytable$V[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "V",]$Freq
  summerytable$W[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "W",]$Freq
  summerytable$Y[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "Y",]$Freq
  #summerytable$X[1]<-tseaminoacid[tseaminoacid$aminoacid=="X",]$Freq
  summerytable$Z[1] <-
    tseaminoacid[tseaminoacid$aminoacid == "Z",]$Freq
  
  
  
  
  #####################################################################
  aminoacid <- character(length = nrow(ara))
  # frequency 23 classes
  for (i in 1:nrow(ara)) {
    id <- tolower(ara$araidentity[i])
    
    if (id == "ala")
      aminoacid[i] <- "A"
    if (id == "arg")
      aminoacid[i] <- "R"
    if (id == "asn")
      aminoacid[i] <- "N"
    if (id == "asp")
      aminoacid[i] <- "D"
    if (id == "cys")
      aminoacid[i] <- "C"
    if (id == "gln")
      aminoacid[i] <- "Q"
    if (id == "glu")
      aminoacid[i] <- "E"
    if (id == "gly")
      aminoacid[i] <- "G"
    if (id == "his")
      aminoacid[i] <- "H"
    if (id == "ile")
      aminoacid[i] <- "I"
    if (id == "start")
      aminoacid[i] <- "start"
    if (id == "leu")
      aminoacid[i] <- "L"
    if (id == "lys")
      aminoacid[i] <- "K"
    if (id == "met")
      aminoacid[i] <- "M"
    if (id == "phe")
      aminoacid[i] <- "F"
    if (id == "pro")
      aminoacid[i] <- "P"
    if (id == "ser")
      aminoacid[i] <- "S"
    if (id == "thr")
      aminoacid[i] <- "T"
    if (id == "trp")
      aminoacid[i] <- "W"
    if (id == "tyr")
      aminoacid[i] <-
        "Y"
    if (id == "val")
      aminoacid[i] <-
        "V"
    if (id == "sec")
      aminoacid[i] <-
        "Z"
    if (id == "stop")
      aminoacid[i] <-
        "stop"
    
  }
  
  tseaminoacid <- data.frame(table(aminoacid))
  summerytable$A[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "A", ]$Freq
  summerytable$C[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "C", ]$Freq
  summerytable$D[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "D", ]$Freq
  summerytable$E[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "E", ]$Freq
  summerytable$F[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "F", ]$Freq
  summerytable$G[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "G", ]$Freq
  summerytable$H[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "H", ]$Freq
  summerytable$I[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "I", ]$Freq
  summerytable$K[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "K", ]$Freq
  summerytable$L[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "L", ]$Freq
  summerytable$M[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "M", ]$Freq
  summerytable$N[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "N", ]$Freq
  summerytable$P[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "P", ]$Freq
  summerytable$Q[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "Q", ]$Freq
  summerytable$R[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "R", ]$Freq
  summerytable$S[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "S", ]$Freq
  summerytable$T[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "T", ]$Freq
  summerytable$V[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "V", ]$Freq
  summerytable$W[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "W", ]$Freq
  summerytable$Y[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "Y", ]$Freq
  #summerytable$X[1]<-tseaminoacid[tseaminoacid$aminoacid=="X",]$Freq
  summerytable$Z[2] <-
    tseaminoacid[tseaminoacid$aminoacid == "Z", ]$Freq
  
  
  
  #####################################################################
  aminoacid <- character(length = nrow(integrated_tse_ara))
  # frequency 23 classes
  for (i in 1:nrow(intersect_tse_ara)) {
    id <- tolower(intersect_tse_ara$tseidentity[i])
    if (id == "undet")
      id <- tolower(intersect_tse_ara$araidentity[i])
    if (id == "ala")
      aminoacid[i] <- "A"
    if (id == "arg")
      aminoacid[i] <- "R"
    if (id == "asn")
      aminoacid[i] <- "N"
    if (id == "asp")
      aminoacid[i] <- "D"
    if (id == "cys")
      aminoacid[i] <- "C"
    if (id == "gln")
      aminoacid[i] <- "Q"
    if (id == "glu")
      aminoacid[i] <- "E"
    if (id == "gly")
      aminoacid[i] <- "G"
    if (id == "his")
      aminoacid[i] <- "H"
    if (id == "ile")
      aminoacid[i] <- "I"
    if (id == "start")
      aminoacid[i] <- "start"
    if (id == "leu")
      aminoacid[i] <- "L"
    if (id == "lys")
      aminoacid[i] <- "K"
    if (id == "met")
      aminoacid[i] <- "M"
    if (id == "phe")
      aminoacid[i] <- "F"
    if (id == "pro")
      aminoacid[i] <- "P"
    if (id == "ser")
      aminoacid[i] <- "S"
    if (id == "thr")
      aminoacid[i] <- "T"
    if (id == "trp")
      aminoacid[i] <- "W"
    if (id == "tyr")
      aminoacid[i] <-
        "Y"
    if (id == "val")
      aminoacid[i] <-
        "V"
    if (id == "sec")
      aminoacid[i] <-
        "Z"
    if (id == "stop")
      aminoacid[i] <-
        "stop"
    
  }
  
  tseaminoacid <- data.frame(table(aminoacid))
  summerytable$A[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "A", ]$Freq
  summerytable$C[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "C", ]$Freq
  summerytable$D[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "D", ]$Freq
  summerytable$E[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "E", ]$Freq
  summerytable$F[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "F", ]$Freq
  summerytable$G[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "G", ]$Freq
  summerytable$H[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "H", ]$Freq
  summerytable$I[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "I", ]$Freq
  summerytable$K[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "K", ]$Freq
  summerytable$L[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "L", ]$Freq
  summerytable$M[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "M", ]$Freq
  summerytable$N[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "N", ]$Freq
  summerytable$P[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "P", ]$Freq
  summerytable$Q[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "Q", ]$Freq
  summerytable$R[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "R", ]$Freq
  summerytable$S[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "S", ]$Freq
  summerytable$T[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "T", ]$Freq
  summerytable$V[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "V", ]$Freq
  summerytable$W[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "W", ]$Freq
  summerytable$Y[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "Y", ]$Freq
  #summerytable$X[1]<-tseaminoacid[tseaminoacid$aminoacid=="X",]$Freq
  summerytable$Z[4] <-
    tseaminoacid[tseaminoacid$aminoacid == "Z", ]$Freq
  
  
  
  #####################################################################
  
  aminoacid <- character(length = nrow(union_tse_ara))
  unionidentity <- union_tse_ara$tseidentity
  for (i in 1:nrow(union_tse_ara)) {
    if (union_tse_ara$foundby[i] == "ara")
      unionidentity[i] <- union_tse_ara$araidentity[i]
  }
  # frequency 23 classes
  for (i in 1:length(unionidentity)) {
    id <- tolower(unionidentity[i])
    if (id == "undet" & union_tse_ara$foundby[i] == "both")
      id <- tolower(union_tse_ara$araidentity[i])
    if (id == "ala")
      aminoacid[i] <- "A"
    if (id == "arg")
      aminoacid[i] <- "R"
    if (id == "asn")
      aminoacid[i] <- "N"
    if (id == "asp")
      aminoacid[i] <- "D"
    if (id == "cys")
      aminoacid[i] <- "C"
    if (id == "gln")
      aminoacid[i] <- "Q"
    if (id == "glu")
      aminoacid[i] <- "E"
    if (id == "gly")
      aminoacid[i] <- "G"
    if (id == "his")
      aminoacid[i] <- "H"
    if (id == "ile")
      aminoacid[i] <- "I"
    if (id == "start")
      aminoacid[i] <- "start"
    if (id == "leu")
      aminoacid[i] <- "L"
    if (id == "lys")
      aminoacid[i] <- "K"
    if (id == "met")
      aminoacid[i] <- "M"
    if (id == "phe")
      aminoacid[i] <- "F"
    if (id == "pro")
      aminoacid[i] <- "P"
    if (id == "ser")
      aminoacid[i] <- "S"
    if (id == "thr")
      aminoacid[i] <- "T"
    if (id == "trp")
      aminoacid[i] <- "W"
    if (id == "tyr")
      aminoacid[i] <-
        "Y"
    if (id == "val")
      aminoacid[i] <-
        "V"
    if (id == "sec")
      aminoacid[i] <-
        "Z"
    if (id == "stop")
      aminoacid[i] <-
        "stop"
    
  }
  
  tseaminoacid <- data.frame(table(aminoacid))
  summerytable$A[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "A", ]$Freq
  summerytable$C[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "C", ]$Freq
  summerytable$D[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "D", ]$Freq
  summerytable$E[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "E", ]$Freq
  summerytable$F[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "F", ]$Freq
  summerytable$G[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "G", ]$Freq
  summerytable$H[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "H", ]$Freq
  summerytable$I[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "I", ]$Freq
  summerytable$K[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "K", ]$Freq
  summerytable$L[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "L", ]$Freq
  summerytable$M[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "M", ]$Freq
  summerytable$N[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "N", ]$Freq
  summerytable$P[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "P", ]$Freq
  summerytable$Q[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "Q", ]$Freq
  summerytable$R[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "R", ]$Freq
  summerytable$S[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "S", ]$Freq
  summerytable$T[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "T", ]$Freq
  summerytable$V[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "V", ]$Freq
  summerytable$W[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "W", ]$Freq
  summerytable$Y[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "Y", ]$Freq
  #summerytable$X[1]<-tseaminoacid[tseaminoacid$aminoacid=="X",]$Freq
  summerytable$Z[3] <-
    tseaminoacid[tseaminoacid$aminoacid == "Z", ]$Freq
  
}