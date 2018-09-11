# end-analysis of genes found by both aragorn and TSE
# venndiagram of ARA and TSE (pseudo and non-pseudo)

visualization <- function() {
  library(VennDiagram)
  library(ggplot2)
  library(dplyr)
  
  # source('/home/fatemeh//Leishmania_Aug2018/Tse_Ara_Integration/Integrate_Tse_Ara.R')
  # integrated_tse_ara <- Integrate_Tse_Ara()
  geneDF <- integrated_tse_ara
  isboth <- geneDF$foundby == "both"
  istseonly <- geneDF$foundby == "tse"
  isaraonly <- geneDF$foundby == "ara"
  # I know that aragorn reported 0 pseudo, so the pseudos are for tse
  ispseudo <- geneDF$tsenote == "pseudo"
  both <- geneDF[isboth,]
  tseonly <- geneDF[istseonly,]
  araonly <- geneDF[isaraonly,]
  tsepseudo <- geneDF[ispseudo,]
  ispseudo_ara <- tsepseudo$foundby == "both"
  pseudo_ara <- tsepseudo[ispseudo_ara,]
  
  grid.newpage()
  draw.triple.venn(
    area1 = nrow(tseonly) + nrow(both) ,
    area2 = nrow(araonly) + nrow(both),
    area3 = nrow(tsepseudo),
    n12 = nrow(both),
    n23 = nrow(pseudo_ara),
    n13 = nrow(tsepseudo),
    n123 = nrow(pseudo_ara),
    category = c("tRNAscan2.0", "Aragorn", "Pseudo"),
    lty = "blank",
    fill = c("skyblue", "pink1", "mediumorchid")
  )
  
  #___________________________________________________________________________________________________________
  Fiveprimend = as.integer(both$tsebegin) - as.integer(both$arabegin)
  bad <- is.na(Fiveprimend)
  Fiveprimend <- Fiveprimend[!bad]
  #Difference between tse1end and tse2end
  Threeprimend = as.integer(both$tseend) - as.integer(both$araend)
  bad <- is.na(Threeprimend)
  Threeprimend <- Threeprimend[!bad]
  
  par(mfrow = c(1, 1))
  X = Fiveprimend
  Y = Threeprimend
  data <- data.frame(X, Y)
  data = group_by(data, X, Y)
  data = summarize(data, frequency = n())
  data$frequency = data$frequency #/ sims
  xbreaks <- as.integer(levels(factor(X)))
  ybreaks <- as.integer(levels(factor(Y)))
  total <- sum(data$frequency)
  plottitle <-
    paste(
      "Empirical joint distribution of end displacements of ",
      nrow(both),
      " TryTryp pseudo genes found by Aragorn and tRNAscan2.0",
      sep = ""
    )
  ggplot(data = data, aes(X, Y)) +
    geom_tile(aes(fill = frequency), color = "white") + geom_text(aes(label = frequency)) +
    scale_fill_gradient(low = "lightblue", high = "darkred", name = "frequency") +
    ggtitle(plottitle) +
    theme(plot.margin = unit(c(1.8, .5, 1.75, 1.55), "cm")) + theme(plot.title = element_text(
      color = "#383838",
      face = "bold",
      size = 17,
      vjust = 4,
      hjust = 0.5
    )) +
    theme(
      axis.title.x = element_text(
        color = "#383838",
        face = "bold",
        size = 13,
        vjust = -1.5
      ),
      axis.title.y = element_text(
        color = "#383838",
        face = "bold",
        size = 13,
        vjust = 2
      )
    ) +
    theme(
      axis.text.x = element_text(
        face = "bold",
        color = "#383838",
        size = 10
      ),
      axis.text.y = element_text(
        face = "bold",
        color = "#383838",
        size = 10
      )
    ) +
    xlab("5' (Tse-Ara)") + ylab("3' (Tse-Ara)") + scale_y_continuous(breaks = ybreaks) +
    scale_x_continuous(breaks = xbreaks)
}

#________________________________________________________________________________________________
# We will keep the genes found by both genefinders 
# investigate genes that are not common between genefinders
# first check the 750 genes that are found by just aragorn:
# 218 of them have no intron 
# non of the pseudo genes incommon with ara has any intron
# one of the pseudo genes found by just tse have an intron 
# check how many of those that do not have intron 
# all the genes that have undertermind identity have no introns 


