# Load required libraries
# Advanced image processing package Magick
library(magick)
library("scatterplot3d")
# this library is used to train a perceptron
library(nnet)
library(neuralnet)
library(jpeg)
library(ggplot2)
library(extrafont)

# returns a list of individual bands in jpeg format from the labels in clusteredImg (by kMeans)
getBands <- function(file) {
  
  # read the image 
  img <- readJPEG(file)
  # structure of the image
  str(img)
  # get the dimensions of an image
  imgDm <- dim(img)
  # Assign RGB channels to data frame
  imgRGB <- data.frame(
    x = rep(1:imgDm[2], each = imgDm[1]),
    y = rep(imgDm[1]:1, imgDm[2]),
    R = as.vector(img[,,1]),
    G = as.vector(img[,,2]),
    B = as.vector(img[,,3])
  )
  # ggplot theme to plot image 
  plotTheme <- function() {
    theme(
      panel.background = element_rect(
        size = 3,
        colour = "pink",
        fill = "white"),
      axis.ticks = element_line(
        size = 2),
      panel.grid.major = element_line(
        colour = "blue",
        linetype = "dotted"),
      panel.grid.minor = element_line(
        colour = "red",
        linetype = "dashed"),
      axis.title.x = element_text(
        size = rel(1.2),
        face = "bold"),
      axis.title.y = element_text(
        size = rel(1.2),
        face = "bold"),
      plot.title = element_text(
        size = 20,
        face = "bold",
        vjust = 1.5)
    )
  }
  
  # plot image 
  ggplot(data = imgRGB, aes(x = x, y = y)) + 
    geom_point(colour = rgb(imgRGB[c("R", "G", "B")])) +
    labs(title = "C band  image") +
    xlab("x") +
    ylab("y") + 
    plotTheme()
  
  # number of clusters
  kClusters <- 9
  
  # KMeans algorithm and assign colors to each cluster and plot the clusters
  kMeans <- kmeans(imgRGB[, c("R", "G", "B")], centers = kClusters)
  kColours <- rgb(kMeans$centers[kMeans$cluster,])
  
  gg <- ggplot(data = imgRGB, aes(x = x, y = y)) + 
    geom_point(colour = kColours) +
    labs(title = paste("k-Means Clustering of", kClusters, "Colours")) +
    xlab("x") +
    ylab("y") + 
    plotTheme()
  
  gg
  
  # Save Last plot 
  # ggsave(file="s-arm-cluster.jpeg")
  
  # get bands from clusters
  imgRGB$Cluster <- as.factor(kMeans$cluster)
  imgRGB$ClusterColor <- kColours[kMeans$cluster]
  
  # Get unique cluster labels
  uniqueLabels <- unique(imgRGB$Cluster)
  # Sort labels in ascending order
  sortedLabels <- sort(uniqueLabels)
  
  # Create a list to store individual bands
  bandList <- list()
  h <- imgDm[1]
  w <- imgDm[2]
  
  # Loop through each cluster label
  for (label in sortedLabels) {
    # Create a blank image with the same dimensions as the clustered image
    bandImg <- array(img[0,0,], dim = c(h, w, 3))
    
    for (i in 1:w) {
      for (j in 1:h) {
        # Extract RGB values of the pixel
        pixel <- img[j, i, ]
        
        # Set pixels corresponding to the current label to the pixel values, others remain white
        if (imgRGB$Cluster[(i-1) * h + j] == label) {
          bandImg[j, i,] <- pixel
        }
      }
    }
    
    # Save the image as a JPEG file
    filename <- paste("band_", label, ".jpeg", sep = "")
    writeJPEG(bandImg, target = filename)
    
    # Add the JPEG file to the band list
    bandList[[label]] <- filename
  }
  
  # Return the list of individual bands
  return(bandList)
}


# returns the shape of img in a dataframe consisting of pairs of coordinates of  
# selected pixels on the border of the image at quartile distances from its  centroid
getFeats <- function(img){ 
  # convert image to rgb and copy image
  imgGray <- image_convert(img0, format = "rgb")
  # copy of image gray to recolor it after clustering
  imgC<-imgGray
  # # copy of image gray to recolor it after clustering
  # imgC<-imgGray
  # # # rotate image2bands-C3-C.pngp the image)
  imgI <-  image_data(imgGray)
  image_read(imgI)
  # # copy image array into another array imgI (intensity)
  # imgI <-  image_data(imgGray)
  # image_data(imgGray)
  # imgI<-imgFlip[[1]]
  # image information (height and x`` width) for image function
  imginfo<-image_info(imgGray)
  w <- image_info(imgGray)$width
  h <- image_info(imgGray)$height
  # calculate and draw the centroid of nonwhite pixels. Create xCoordinate and yCoordinate variables as a list
  xCoords <- list()
  yCoords <- list()
  # plot( as.raster(imgGray))
  for(i in 1:w){
    for(j in 1:h){
      if(imgI[1,i,j] != "ff")
      {
        xCoords <- c(xCoords, i)
        yCoords <- c(yCoords, h-j)
      }
    }
  }
  plot(cbind(xCoords,yCoords))
  # Find convex shape and coordinates
  getconvex <- function (xCoord,yCoord){
    pts <- chull(xCoord,yCoord); # df
    mpts <- cbind(xCoord[pts], yCoord[pts])
    # str(dfhull)
    plot(mpts)
    # c.hull
    chull0 <- chull(mpts)
    chull <- mpts[c(chull0, chull0[1]),]
    plot(mpts, pch=19)
    chull
    
  }
  
  # Draw convex haul 
  conv <- getconvex(xCoords,yCoords)
  # Separate coordinates from convex haul
  conXCoord <- as.numeric(conv[,1])
  conYCoord <- as.numeric(conv[,2])
  # Lenght of convex haul coordinates
  lenX <- length(conXCoord)
  # centroid and data extraction 
  centroidX <- ceiling(mean(as.numeric(conXCoord)))        
  centroidY <- ceiling(mean(as.numeric(conYCoord)))
  centroid <- c(centroidX, centroidY)
  plot(image_read(imgI))
  points(centroidX, centroidY, col="blue", pch=20, cex=2)
  # Next calculate Eucleadian distance from the centroid 
  distFromC <-array(0,dim = c(w,h))
  for( i in 1:lenX){
    x <- as.numeric(conXCoord[i])
    y <- as.numeric(conYCoord[i])
    distFromC[x,y] <- ceiling(sqrt((conXCoord[i]-centroidX)^2 + (conYCoord[i]-centroidY)^2))
  }
  # Store non zero distances in dist variable 
  dists <- c()
  for(i in 1:w)
  {
    for( j in 1:h){
      if( distFromC[i,j] != 0 )
      {
        dists <- c(dists, distFromC[i,j])
      }
    }
  }
  # calculate qurtiles 
  qs <- ceiling(quantile(unlist(dists),type = 1))
  plot(image_read(imgI))
  points(centroidX, centroidY, col="blue", pch=20, cex=2)
  # plot nonwhite pixel at max (last quartile) distance  
  print(qs)
  coordns <- c()
  coordns <- which(distFromC==qs[5], arr.ind = T )
  points(coordns,col="red",pch=20, cex=2)
  # c1 <- cbind(coordns[1,]) x1 <- coordns[1,1]
  y1 <- coordns[1,2]
  
  # plot nonwhite pixel at fourth quartile  distance
  coordns <- c()
  coordns <- which( distFromC == qs[4], arr.ind = T )
  # coordns <- which( distFromC<0, arr.ind = T )
  points(coordns, col="green",pch=20, cex=2)
  c2 <- cbind(coordns[1,])
  x2 <- coordns[1,1]
  y2 <- coordns[1,2]
  
  # plot nonwhite pixel at third quartile  distance
  coordns <- c()
  coordns <- which(distFromC == qs[3], arr.ind = T )
  points(coordns, col="green",pch=20, cex=2)
  c3 <- cbind(coordns[1,])
  x3 <- coordns[1,1]
  y3 <- coordns[1,2]
  # str(c3)
  
  # plot nonwhite pixel at second  quartile  distance
  coordns <- c()
  coordns <- which( distFromC==qs[2], arr.ind = T )
  points(coordns,col="green",pch=20, cex=2)
  c4 <- cbind(coordns[1,])
  x4 <- coordns[1,1]
  y4 <- coordns[1,2]
  
  # plot nonwhite pixel at first or minimum quartile  distance
  coordns <- c()
  coordns <- which(distFromC==qs[1], arr.ind = T )
  points(coordns,col="green",pch=20, cex=2)
  c5 <- cbind(coordns[1,])
  x5 <- coordns[1,1]
  y5 <- coordns[1,2]
  
  # create data frame for each band
  # df = rbind(c1,c2,c3,c4,c5)
  # Create the data frame
  df <- data.frame(bandName, colorOFBand, labelB, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5)
  
  
  # append data into CSV file 
  write.table(df , file="data-s-arm-df4.csv", append = T,sep = ",", row.names = FALSE, col.names = FALSE)
} # end getFeats


# main function to start the program 
# Set a working directory
setwd("/Users/naganaresh/Desktop/Personal/MastersProject/Completed/input/")
file <- "chrosome-nkip.jpeg"

# get a list of 15 candidate bands for the bimarkers (name of images)
candBands <- getBands(file)

# get the labels for the bands for the chromosome decomposition.
chrDec <- list()
