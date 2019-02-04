## THIS R SCRIPT DRAWS TRNA INFORMATION BUBBLEPLOTS BASED ON TABULAR
## OUTPUT m , FROM THE PERL SCRIPT logofun_eps2table WHICH IS LOADED AND
## PASSED TO THE ENCLOSED FUNCTION all.bubble.wg 

## COPYLEFT 2011 DAVID H. ARDELL
## ALL WRONGS REVERSED
## REVISED ON FEB 3 2019
bubble.main <- function(){
library("heR.Misc"); ## THIS IS REQUIRED FOR THE BUBBLEPLOT FUNCTION AND MUST BE DOWNLOADED FROM
## http://exposurescience.org/her.html

# THESE ARE THE COORDINATES FOR THE TRNA STRUCTURE BACKBONE IN THE FIGURES

# THESE ARE THE COORDINATES FOR THE TRNA STRUCTURE BACKBONE IN THE FIGURES

# positions x= 8.75 and y= 10.5 are manually added

line.x <- c(6.875,6.500,6.125,5.750,5.375,5.000,4.625,4.625,5.000,5.000,2.375
            ,2.750,2.375,2.750,2.500,2.875,3.250,2.875,2.500,2.125,1.750,1.375,1.000
            ,0.625,0.250,0.625,1.000,1.500,1.125,1.500,1.125,1.500,1.125,1.500,1.125
            ,1.500,1.125,0.625,1.000,1.375,1.750,2.125,2.500,2.875,2.375,2.750,2.375
            ,2.750,2.375,4.250,4.250,3.875,4.250,3.875,4.250,4.250,3.875,3.500,3.125
            ,2.750,2.250,1.875,1.500,1.125,1.500,1.875,2.250,2.750,3.125,3.500,3.875
            ,4.250,4.625,5.000,5.375,5.750,6.125,6.500,6.875,7.250,7.625,8.000,8.375,8.75);

line.y <- c(8.875,8.500,8.875,8.500,8.875,8.500,8.875,7.375,7.000
            ,3.500,3.500,3.875,4.250,4.625,5.125,5.500,5.875,6.250,6.625
            ,7.000,7.375,7.000,6.625,6.250,5.875,5.500,5.125,4.625,4.250
            ,3.875,3.500,3.125,2.750,2.375,2.000,1.625,1.250,0.750,0.375
            ,0.000,-0.375,0.000,0.375,0.750,1.250,1.625,2.000,2.375,2.750
            ,2.750,3.875,4.250,4.625,5.000,5.375,8.500,8.875,8.500,8.875
            ,8.500,8.000,8.375,8.750,9.125,9.500,9.875,10.250,9.750,10.125
            ,9.750,10.125,9.750,10.125,9.750,10.125,9.750,10.125,9.750,10.125
            ,9.750,10.125,9.750,10.125,10.5);

## THESE ARE FOR THE SPRINZL COORD LABELS
coord.labels <- c("1","5","10","14A","18","20C","25","30","35","40","45","50","55","60","65","70");
coord.labels.x <- c(6.875,5.375,2.375,2.875,1.750,0.625,1.125,1.500,1.750,2.750,3.875,3.875,1.875,2.250,4.250,6.125);
coord.labels.y <- c(8.875,8.875,3.500,5.500,7.375,5.500,3.500,1.625,-0.375,1.625,4.250,8.875,8.375,10.250,9.750,10.125);

xbump <- 0.35;
ybump <- 0.4;

up <- c(5,13,14,15,16,17);
up.coord.labels <- coord.labels[up];
up.coord.labels.x <- coord.labels.x[up];
up.coord.labels.y <- coord.labels.y[up] + ybump;

dn <- c(1,2,9,11,12);
dn.coord.labels <- coord.labels[dn];
dn.coord.labels.x <- coord.labels.x[dn];
dn.coord.labels.y <- coord.labels.y[dn] - ybump;

lt <- c(3,4,7);
lt.coord.labels <- coord.labels[lt];
lt.coord.labels.x <- coord.labels.x[lt] - xbump;
lt.coord.labels.y <- coord.labels.y[lt];

rt <- c(6,8,10);
rt.coord.labels <- coord.labels[rt];
rt.coord.labels.x <- coord.labels.x[rt] + xbump;
rt.coord.labels.y <- coord.labels.y[rt];




tRNA_L_skel_df <-
  read.table(
    "/home/fatemeh/Leishmania_Aug2018/tsfmVisualization/tRNA_L_skel.txt",
    sep = ",",
    header = FALSE
  )
names(tRNA_L_skel_df) <- c("sprinzl", "x", "y", "sprinzl2")
dirpath <-
  "/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfminput_2019/tRNADBClustered/Logos/"
clusterdir <- list.dirs(path = dirpath, recursive = FALSE)
for (i in 1:length(clusterdir)) {
  splitedpath <- unlist(strsplit(clusterdir[i], split = "/"))
  clus_name <- splitedpath[length(splitedpath)]
  if (clus_name == "HOMO")
    next
  tablepath <-
    paste(
      "/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfminput_2019/tRNADBClustered/Logos/",
      clus_name,
      "/Bubble/",
      clus_name,
      "_Table.txt",
      sep = ""
    )
  df <- read.table(tablepath, header = TRUE)
  df <- match_bubble_coords(df, tRNA_L_skel_df)
  outputpath <-
    paste(
      "/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfminput_2019/tRNADBClustered/Logos/",
      clus_name,
      "/Bubble/",
      sep = ""
    )
  all.bubble(
    df,
    name = "bubble",
    clus_name,
    outputpath,
    alpha = 0.5,
    fact = 0.5,
    area = TRUE,
    legend = FALSE
  )
}

}
all.bubble <- function(df,name="bubble",clus_name,outputpath,alpha=0.5,fact=0.5,area=TRUE,legend=FALSE) { 
  gains <- (df$gainbits * df$gainfht);
  convs <- (df$convbits * df$convfht);
  map2rgb <- function (c) { rgb(t(col2rgb(c))/255,alpha=alpha);}
  colormap <- function (g,c) { 
     y <- rep(0,length(g));
     y[g <  0.48            & c < 0.44]             <- rgb(t(col2rgb("white"))/255,alpha=alpha);
     y[g >= 0.48 & g < 0.95 & c < 0.44]             <- rgb(t(col2rgb("darkred"))/255,alpha=alpha);
     y[g >= 0.95            & c < 0.44]             <- rgb(t(col2rgb("red"))/255,alpha=alpha);
     y[g <  0.48            & c >= 0.44 & c < 0.70] <- rgb(t(col2rgb("darkblue"))/255,alpha=alpha);
     y[g >= 0.48 & g < 0.95 & c >= 0.44 & c < 0.70] <- rgb(t(col2rgb("darkmagenta"))/255,alpha=alpha);
     y[g >= 0.95            & c >= 0.44 & c < 0.70] <- map2rgb("deeppink");
     y[g <  0.48            & c >= 0.70]            <- map2rgb("blue");
     y[g >= 0.48 & g < 0.95 & c >= 0.70]            <- map2rgb("blueviolet");
     y[g >= 0.95            & c >= 0.70]            <- map2rgb("magenta");
     y;
  }
  colors <- colormap(gains,convs);
  
  widthmap <- function (g,c) {
     y <- rep(1,length(g));
     y[g < 0.48 & c < 0.44] <- 1;
     y[g >= 0.48 | c >= 0.44] <- 2;
     y[g >= 0.95 | c >= 0.70] <- 3;	
     y;
  }
  widths <- widthmap(gains,convs);
  
  for (class in levels(df$aa)) {
    # add the first part of filename 
    filenm <- paste(outputpath,clus_name,"_",name,"_",class,".pdf",sep="");
    pdf(file=filenm,version="1.4");
    bubbleplot(
             df$x[df$aa == class],
             df$y[df$aa == class],
	     (df$fbits[df$aa == class] * df$fht[df$aa == class]),
	     fact=fact, #0.265165 = sqrt(2*0.375^2)/2
	     area=area, 
	     fg = rgb(t(col2rgb("black")/255)),
	     bg = colors[df$aa == class],
	     box=FALSE,
	     axes=FALSE,
	     lwd=2
	     );
    lines(line.x,line.y);
    text(labels=up.coord.labels,x=up.coord.labels.x,y=up.coord.labels.y);
    text(labels=dn.coord.labels,x=dn.coord.labels.x,y=dn.coord.labels.y);
    text(labels=lt.coord.labels,x=lt.coord.labels.x,y=lt.coord.labels.y);
    text(labels=rt.coord.labels,x=rt.coord.labels.x,y=rt.coord.labels.y);
    if (legend) {
      legend.x <- rep(df$x[df$aa == "X" &  df$state == "A" & 
                    (df$sprinzl == "68" |  df$sprinzl == "70" |
                     df$sprinzl == "72" )] + c(0.2,0.3,0.4),3);  
      legend.y <- c(rep(df$y[df$aa == "X" &  df$state == "A" & df$sprinzl == "31"],3),
      	       	    rep(df$y[df$aa == "X" &  df$state == "A" & df$sprinzl == "27"],3), 
                    rep(df$y[df$aa == "X" &  df$state == "A" & df$sprinzl == "23"],3));
      legend.z <- rep(2.2,9);
      legend.c <- colormap(c(0,0,0,0.5,0.5,0.5,1,1,1),c(0,0.5,1,0,0.5,1,0,0.5,1));
      bubbleplot(legend.x, legend.y, legend.z,fact=fact,area=area,bg=legend.c,
	     add=TRUE,
	     box=FALSE,axes=FALSE,lwd=2);

    } 
    prime.x <- df$x[df$aa == "X" &  df$state == "A" &
                (df$sprinzl == "1" |  df$sprinzl == "75")];
    prime.y <- df$y[df$aa == "X" &  df$state == "A" &
                (df$sprinzl == "1" |  df$sprinzl == "75")];
    text(labels=c("5'","3'"),x=prime.x,y=prime.y,adj=c(-0.9,0));	
    dev.off();
  }
}

match_bubble_coords <- function(df, tRNA_L_skel_df) {
  for (i in 1:nrow(tRNA_L_skel_df)) {
    samecoord = df$coord == tRNA_L_skel_df$sprinzl[i]
    df$x[samecoord] = tRNA_L_skel_df$x[i]
    df$y[samecoord] = tRNA_L_skel_df$y[i]
    df$sprinzl[samecoord] = tRNA_L_skel_df$sprinzl[i]
  }
  df$x[df$coord == 74] = 8.75
  df$y[df$coord == 74] = 10.5
  df
}
  
  
