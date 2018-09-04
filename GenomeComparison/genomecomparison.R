# read all the genomes from file
# for each genome count number of sequences by counting number of headlines
# calculate the average length of the sequences
# plot a diagram with x= name of the genome, seqcountcomparison = (number of sequence/length of sequence)*average length of all genomes
# seqlencomparison = (average length of sequences / length of genome) * 100
# I can also read the gene file and sort it based on organism name
# I cannot use SourceSeq because many of the genes do not have this field
genomecomparison <- function() {
  library(gsubfn)
  library(ggplot2)
  filenames <-
    list.files(
      "/home/fatemeh/Leishmania_Aug2018/tritrypdb.org/",
      pattern = "*.fasta",
      full.names = TRUE
    )
  #filenames <- filenames[[1]]
  # array of seqlencomparison seqcountcomparison
  seqlencomp <- integer(length = length(filenames))
  seqcountcomp <- integer(length = length(filenames))
  genomelength <- integer(length = length(filenames))
  for (i in 1:length(filenames)) {
    headers <- read.pattern(filenames[i], "^>+.*")
    genomelen = 0
    myseqlen = 0
    # extract the length of the sequence from this line by taking the field after "to"
    # sum the length of all sequences as the length of the genome
    
    for (r in 1:nrow(headers)) {
      fields <- unlist(strsplit(as.character(headers[r,1]), " "))
      for (j in 1:length(fields)) {
        if(length(grep("length=*",fields[j])) != 0)
        {
          myseqlen <- substring(fields[j],8)
          myseqlen <- as.integer(myseqlen)
          break
        }
      }
      genomelen <- genomelen + myseqlen
    }
    genomelength[i] <- genomelen
    seqcountcomp[i] <- nrow(headers) / genomelen
    seqaveragelen <- genomelen / nrow(headers)
    seqlencomp[i] <- ( seqaveragelen / genomelen ) * 100000
  }
  seqcountcomp <- ave(genomelength) * seqcountcomp
  genomenames <- substring(filenames,55)
  info <- data.frame(genomenames,seqcountcomp,seqlencomp)
  names(info) <- c("genome","seqcountcomp", "seqlencomp")
  # ggplot(info, aes(y=info$seqcountcomp, x=info$genome)) + 
  #        geom_bar(stat="identity")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  mymatrix <- rbind(seqcountcomp,seqlencomp)
  gnames <- lapply(info$genome,function(X){gsub(x = X,pattern = "*.fasta",replacement = "")})
  names <- as.character(gnames)
  gnames <- unlist(gnames)
  colnames(mymatrix) <- seq(1,length(filenames),1)
  barplot(mymatrix,beside = TRUE,col = c("blue","red"),las=2,cex.names=0.6)
  gnames <- paste(seq(1,length(filenames),1),gnames)
  legend(2,25000,legend = gnames,col = seq(1,length(filenames),1),cex=0.6)
  legend(40,25000,legend = c("number of sequences","length of sequences"),col = c("blue","red"),cex=0.6,lty = 1)

}