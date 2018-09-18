initiatorDetecting <- function(){
 
  # the input is the filtered gene file from Integrate_Tse_Ara script. (geneDF6)
  iscat <- integrated_tse_ara$araac=="(cat)" & integrated_tse_ara$tseac=="CAT"
  catDF <- integrated_tse_ara[iscat,]
  # we have 173 cat genes found by both ara and tse
  # we have 4 cat ara genes not detected by tse
  # we have 2 genes found by both ara and tse but they are cat in tse and cta in ara
  # what to do with these 6 genes ?
  geneseqs <- character(length = nrow(catDF))
  for (i in 1:nrow(catDF)) {
    if (catDF$foundby[i] == "both")
      geneseqs[i] <- catDF$tsegeneseq[i]
    else if (catDF$foundby[i] == "tse")
      geneseqs[i] <- catDF$tsegeneseq[i]
    else if (catDF$foundby[i] == "ara")
      geneseqs[i] <- catDF$arageneseq[i]
  }
  m <- matrix(nrow = nrow(catDF), ncol = nrow(catDF))
  distanceDF <- as.data.frame(m)
  for (i in 1:length(geneseqs)) {
    for (j in 1:length(geneseqs)) {
      distanceDF[i, j] <- adist(geneseqs[i], geneseqs[j])
    }
  }
  

  ####################################################### Hierarchical Cluster Analysis ######################################
  # Hierarchical clustering is an alternative approach to k-means clustering for identifying groups in the dataset.
  # we used Agglomerative clustering which works in a bottom-up manner 
  # we measure the dissimilarity between two clusters of observations with Ward.D2 which It minimizes the total within-cluster variance. 
  # At each step the pair of clusters with minimum between-cluster distance are merged
  
  # libraries 
  #library(tidyverse)  # data manipulation
  library(cluster)    # clustering algorithms
  library(factoextra) # clustering visualization
  library(dendextend) # for comparing two dendrograms
  
  # The height of the fusion, provided on the vertical axis, indicates the (dis)similarity between two observations.
  # Note that, conclusions about the proximity of two observations can be drawn only based on the height where branches 
  # containing those two observations first are fused. We cannot use the proximity of two observations along the horizontal axis as a criteria of their similarity
  hc <- hclust(as.dist(distanceDF), method='ward.D2')
  plot(hc, cex = 0.6, hang = -1)

  # The height of the cut to the dendrogram controls the number of clusters obtained. It plays the same role as the k in k-means clustering.
  # 
  rect.hclust(hc, k = 3, border = 2:5)
  sub_groups <- cutree(hc, k = 3)
  fviz_cluster(list(data = distanceDF, cluster = sub_groups))
  
  
  # extracting minimum and maximum 
  subg1 <- sub_groups == 1 
  subg2 <- sub_groups == 2
  subg3 <- sub_groups == 3
  
  cluster1 <- paste(min(distanceDF[subg1,subg1]),max(distanceDF[subg1,subg1]),sep = "-")
  cluster2 <- paste(min(distanceDF[subg2,subg2]),max(distanceDF[subg2,subg2]),sep = "-")
  cluster3 <- paste(min(distanceDF[subg3,subg3]),max(distanceDF[subg3,subg3]),sep = "-")
  
  clusterDF <- data.frame(table(sub_groups))
  clusterDF$editdistanceRange <- c(cluster1,cluster2,cluster3)
  
  
  
  ####################################################### 
  # clustering 8:30 am monday november 
  mydata <- distanceDF
  # K-Means Clustering with 2 clusters
  fit <- kmeans(mydata, 3)
  
  # Cluster Plot against 1st 2 principal components
  # vary parameters for most readable graph
  library(cluster) 
  clusplot(mydata, fit$cluster, color=TRUE, shade=TRUE, 
           labels=2, lines=0)
  
  # Centroid Plot against 1st 2 discriminant functions
  library(fpc)
  plotcluster(mydata, fit$cluster)
  
  #cluster number 2 is the one that their distances are lowest and hopefuly they are initiators
  isinitiator <- fit$cluster==fit$cluster[2]
  initiatorDF <- catDF[isinitiator,]
  
}