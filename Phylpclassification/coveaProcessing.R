coveaProcessing <- function() {
  readSeqsIntoDf()
  library(readr)
  seqDB <-
    read_csv("/home/fatemeh/Leishmania_Aug2018/phyloclassification/Leishmania_HomoC/HomoTryTryp_coveaDF.txt")
  SSDB <-
    read_csv("/home/fatemeh/Leishmania_Aug2018/phyloclassification/Leishmania_HomoC/HomoTryTryp_coveaDF_SS.txt")
  editAlignment(seqDB, SSDB)
  
  
}
#####################################################################################
# # of columns in the dataframe is: # of sequences + 1 CS line
# name of the columns are "CS","seqname1", "seqname2", ...
readSeqsIntoDf <- function()
{
  CS <-
    grep(
      "#=CS +",
      readLines(
        "/home/fatemeh/Leishmania_Aug2018/phyloclassification/Leishmania_HomoC/LeishmaniaHomoC.covea"#"/home/fatemeh/Leishmania_Aug2018/phyloclassification/Trytryp_genes_NoVarIntron.covea" #
      ),
      value = TRUE
    )
  CStemp <-
    unlist(strsplit(CS, split = "\\s+"))
  CS <- CStemp[seq(2, length(CStemp), 2)]
  CS <- paste(CS, collapse = '')
  CSarr <- substring(CS, seq(1, nchar(CS), 1), seq(1, nchar(CS), 1))
  
  # read the name of the sequences into variable "seqnames"
  SQs <- grep(
    "#=SQ +",
    readLines(
      "/home/fatemeh/Leishmania_Aug2018/phyloclassification/Leishmania_HomoC/LeishmaniaHomoC.covea"#"/home/fatemeh/Leishmania_Aug2018/phyloclassification/Trytryp_genes_NoVarIntron.covea" #
    ),
    value = TRUE
  )
  seqnames <- character(length = length(SQs))
  for (i in 1:length(SQs)) {
    seqnames[i] <- unlist(strsplit(SQs[i], split = " "))[2]
  }
  
  # define the main data frame as "seqDB" and assign names to it
  m <- matrix(ncol = (length(seqnames) + 1), nrow = length(CSarr))
  seqDB  <- as.data.frame(m)
  names(seqDB) <- c(seqnames, "CS")
  
  for (i in 1:length(seqnames)) {
    pat <- paste("^", seqnames[i], " +", sep = "")
    myseq <-
      grep(
        pattern = pat,
        readLines(
          "/home/fatemeh/Leishmania_Aug2018/phyloclassification/Leishmania_HomoC/LeishmaniaHomoC.covea" #"/home/fatemeh/Leishmania_Aug2018/phyloclassification/Trytryp_genes_NoVarIntron.covea"#
        ),
        value = TRUE
      )
    temp <-
      unlist(strsplit(myseq, split = "\\s+"))
    myseq <- temp[seq(2, length(temp), 2)]
    myseq <- paste(myseq, collapse = '')
    myseqarr <-
      substring(myseq, seq(1, nchar(myseq), 1), seq(1, nchar(myseq), 1))
    seqDB[, i] <- myseqarr
  }
  seqDB$CS <- CSarr
  
  #_____________________ reading the #=SS lines into a nother data frame as SSDB ____________
  # define the main data frame as "seqDB" and assign names to it
  
  SSDB = seqDB
  
  SSs <- grep(
    "#=SS +",
    readLines(
      "/home/fatemeh/Leishmania_Aug2018/phyloclassification/Leishmania_HomoC/LeishmaniaHomoC.covea"#"/home/fatemeh/Leishmania_Aug2018/phyloclassification/Trytryp_genes_NoVarIntron.covea"#
    ),
    value = TRUE
  )
  
  # number of sections is CStemp/2 = 3
  numsec <- length(CStemp) / 2
  end <- length(SSs) / numsec
  for (i in 1:end) {
    if (numsec == 3)
      myseq <-
        paste(unlist(strsplit(SSs[i], split = "\\s+"))[2],
              unlist(strsplit(SSs[end + i], split = "\\s+"))[2],
              unlist(strsplit(SSs[(2 * end) + i], split = "\\s+"))[2],
              sep = "")
    if (numsec == 2)
      myseq <-
        paste(unlist(strsplit(SSs[i], split = "\\s+"))[2], unlist(strsplit(SSs[end +
                                                                                 i], split = "\\s+"))[2], sep = "")
    
    myseqarr <-
      substring(myseq, seq(1, nchar(myseq), 1), seq(1, nchar(myseq), 1))
    SSDB[, i] <- myseqarr
  }
  
  library(readr)
  write_csv(seqDB, path  = "/home/fatemeh/Leishmania_Aug2018/phyloclassification/Leishmania_HomoC/HomoTryTryp_coveaDF.txt")
  write_csv(SSDB, path  = "/home/fatemeh/Leishmania_Aug2018/phyloclassification/Leishmania_HomoC/HomoTryTryp_coveaDF_SS.txt")
  
  
}
#####################################################################################
writeCovea <- function(SSDB, seqDB) {
  coveaseqs <- character(length = ncol(seqDB) - 1)
  coveass <- character(length = ncol(seqDB) - 1)
  
  seqDB <- seqDB[, -ncol(seqDB)]
  SSDB <- SSDB[, -ncol(SSDB)]
  for (i in 1:(length(coveaseqs))) {
    # coveaseqs[i] <- paste(seqDB[, i], collapse = '')
    # coveass[i] <- paste(SSDB[, i], collapse = '')
    print(i)
    coveaseqs[i] <-
      paste((as.data.frame(seqDB[, i])[, 1]), collapse = '')
    coveass[i] <- paste((as.data.frame(SSDB[, i])[, 1]), collapse = '')
    
  }
  
  mynames <- names(seqDB)[!names(seqDB) %in% "CS"]
  
  library(gdata)
  write.fwf(
    data.frame(mynames,
               coveaseqs,
               coveass),
    file = "/home/fatemeh/Leishmania_Aug2018/phyloclassification/Leishmania_HomoC/HomoTryTryp_EditedCovea.covea",
    sep = "\n",
    colnames = FALSE
  )
  
  # remove "."s from sequences and write them in a fasta file to run with covea
  library(Hmisc)
  for (i in 1:length(coveaseqs)) {
    coveaseqs[i] <- translate(coveaseqs[i], "[.]", "-")
    #coveaseqs[i] <- translate(coveaseqs[i], "[n]", "-")
    #coveaseqs[i] <- translate(coveaseqs[i], "[N]", "-")
  }
  write.fwf(
    data.frame(paste(">", mynames, sep = ""),
               coveaseqs),
    file = "/home/fatemeh/Leishmania_Aug2018/phyloclassification/Leishmania_HomoC/HomoTryTryp_EditedCovea.fasta",
    sep = "\n",
    colnames = FALSE
  )
  
  
}
#####################################################################################
editAlignment <- function(seqDB, SSDB) {
  # removing sites that have more than 98% gap (removing rows that have more than 99% ".")
  # seqDB = mydb
  delpos = " "
  for (i in 1:nrow(seqDB)) {
    temp <-
      data.frame(table(
        seqDB[i, ] == "." |
          seqDB[i, ] == "t" |
          seqDB[i, ] == "g" | seqDB[i, ] == "c" | seqDB[i, ] == "a"
      ))
    if (length(temp[temp$Var1 == "FALSE", 2]) != 0)
    {
      if (temp[temp$Var1 == "FALSE", 2] != ncol(seqDB))
      {
        gapperc <-
          (temp[temp$Var1 == "TRUE", 2] / (temp[temp$Var1 == "TRUE", 2] + temp[temp$Var1 ==
                                                                                 "FALSE", 2])) * 100
        if (gapperc > 97)
        {
          delpos = c(delpos, i)
        }
        
      }
    }
    else{
      if (temp[temp$Var1 == "TRUE", 2] == ncol(seqDB))
        delpos = c(delpos, i)
    }
  }
  delpos <- delpos[2:length(delpos)]
  delpos = as.integer(delpos)
  seqDB <- seqDB[-delpos, ]
  SSDB <- SSDB[-delpos, ]
  
  cs <- paste(seqDB$CS, collapse = '')
  
  # remove sequences with more than 8 gaps!
  delpos = " "
  flag=0
  for (i in 1:ncol(seqDB)) {
    temp <- data.frame(table(seqDB[, i] == "."))
    if (length(temp[temp$Var1 == "TRUE", 2]) != 0)
      if (temp[temp$Var1 == "TRUE", 2] > 8)
      {
        if (names(seqDB[, i]) != "CS")
          delpos = c(delpos, i)
      }
    for (x in 1:length(unlist(seqDB[,i]))) {
      if(unlist(seqDB[,i])[x]=="n" | unlist(seqDB[,i])[x]=="N" )
        flag=1
    }
    if(flag == 1)
    {
      print("shit")
      delpos = c(delpos, i)
      flag = 0
    }
  }
  if (length(delpos) > 1)
  {
    delpos <- delpos[2:length(delpos)]
    delpos = as.integer(delpos)
    seqDB <- seqDB[, -delpos]
    SSDB <- SSDB[, -delpos]
  }
  
  # remove sites that are lowercase or . in 99% of sequences
  
  # remove sequences with gap in their anticodon
  
  writeCovea(SSDB, seqDB)
  
}
#####################################################################################
# partition the fasta file based on their functional classes
# partition the EditedCovea.fasta files based on the sourceOrg 
