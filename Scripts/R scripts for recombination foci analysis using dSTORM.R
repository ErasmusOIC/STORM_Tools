#############################################################
### Script: "R Scripts for recombination foci analysis "  ###
### Authors: Lieke Koornneef & Johan Slotman              ###
### Affiliation: Erasmus MC, Rotterdam, The Netherlands   ###
### Contact: j.slotman@erasmusmc.nl  						          ###
### License: LGPLv3                                       ###
### Date: 01-02-2023                                      ###
#############################################################

# Installation SMoLR package ----------------------------------------------
install.packages("devtools")  
library(devtools)  
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install("EBImage") 
install_github("ErasmusOIC/SMoLR", build_vignettes = TRUE)
library(SMoLR) 


# Create dSTORM image per channel ------------------------------------------

#Load packages
library(SMoLR)
library(data.table)

#Select documents
dir <- "D:/"
nr_ch <- 3       #number of channels
name <- "Exp1"   #name of experiment

#Create Images 
roi <- paste(dir, name, "_Roi.roi", sep ="")
roi_chr <- paste(dir, name, "_Roi_characteristics.csv", sep ="")
roi_size <- read.csv(roi_chr)
X <- roi_size$BX*200
width <- roi_size$Width*200
Y <- roi_size$BY*200
height <- roi_size$Height*200
profile <- SMOLR_PROFILE("Elyra",c("Index","First Frame","Number Frames","Frames Missing","Position X [nm]","Position Y [nm]","Precision [nm]","Number Photons","Background variance","Chi square","PSF half width [nm]","Channel","Z Slice"),c(5,6,7,12),skip=0)

file1 <- paste(dir, name, "_Ch1_g.txt", sep ="")
image1 <- paste(dir, name, "_Ch1.tif", sep ="")
data1 <- SMOLR_FAST_IMPORT(file1,profile)
data1 <- data.frame(data1)
data_sub <- IJROI_subset(x = data1,file = roi,pxsize = 5)
rm(data1)
smolr <- SMOLR(data_sub,output = "tiff",file = image1, xlim=c(X*5,(X+width)*5), ylim=c(Y*5,(Y+height)*5))
rm(data_sub)
print("Image of channel 1 is generated")

if(nr_ch>1) {
  file2 <- paste(dir, name, "_Ch2_gt.txt", sep ="")
  image2 <- paste(dir, name, "_Ch2.tif", sep ="")
  data2 <- SMOLR_FAST_IMPORT(file2,profile)
  data2 <- data.frame(data2)
  data_sub <- IJROI_subset(x = data2,file = roi,pxsize = 5)
  rm(data2)
  smolr <- SMOLR(data_sub,output = "tiff",file = image2, xlim=c(X*5,(X+width)*5), ylim=c(Y*5,(Y+height)*5))
  rm(data_sub)
  print("Image of channel 2 is generated")
}
if(nr_ch>2) {
  file3 <- paste(dir, name, "_Ch3_gt.txt", sep ="")
  image3 <- paste(dir, name, "_Ch3.tif", sep ="")
  data3 <- SMOLR_FAST_IMPORT(file3,profile)
  data3 <- data.frame(data3)
  data_sub <- IJROI_subset(x = data3,file = roi,pxsize = 5)
  rm(data3)
  smolr <- SMOLR(data_sub,output = "tiff",file = image3, xlim=c(X*5,(X+width)*5), ylim=c(Y*5,(Y+height)*5))
  rm(data_sub)
  print("Image of channel 3 is generated")
}



# Combine channels --------------------------

#Load packages
library(SMoLR)
library(data.table)

#Select documents
dir <- "D:/"
nr_ch = 3
name = "Exp1"

#Create total localization file 
file1 <- paste(dir, name, "_Ch1_g.txt", sep ="")
data1 <- SMOLR_FAST_IMPORT(file1,profile)
data1 <- data.frame(data1)

if(nr_ch>1) {
  file2 <- paste(dir, name, "_Ch2_gt.txt", sep ="")
  data2 <- SMOLR_FAST_IMPORT(file2,profile)
  data2 <- data.frame(data2)
  data2$Channel <- 2
}
if(nr_ch>2) {
  file3 <- paste(dir, name, "_Ch3_gt.txt", sep ="")
  data3 <- SMOLR_FAST_IMPORT(file3,profile)
  data3 <- data.frame(data3)
  data3$Channel <- 3
}
if(nr_ch == 1){
  alldata <- data1
}
if(nr_ch == 2){
  alldata <- rbind(data1,data2)
}
if(nr_ch == 3){
  alldata <- rbind(data1,data2,data3)
}

names(alldata)[1:13] <- c("Index", "First Frame", "Number Frames", "Frames Missing", "Position X [nm]", "Position Y [nm]", "Precision [nm]", "Number Photons", "Background variance", "Chi square", "PSF width [nm]", "Channel", "Z Slice")            
lastrow <- c("",0,0,0,0,0,0,0,0,0,0,0,0)
lastrow <- as.numeric(lastrow)
alldata <- rbind(alldata,lastrow)
write.table(alldata, file = paste(dir,name,"_alldata.txt", sep=""), quote = FALSE, sep= "\t", row.names = FALSE, col.names = TRUE)



# Analysis of ROIs --------------------------------------------------------

#Load packages
library(SMoLR)
library(data.table)
library(EBImage)
library(plyr)

#Select documents
dir <- "D:/"
name = "Exp1"
thr_kde = 0.15   #Threshold of KDE in localizations / nm^2
thr_area = 50    #Threshold of nanofoci in pixels

#Analysis
alldata <- read.delim(paste(dir,name,"_alldata.txt", sep=""), header = TRUE)
names(alldata)[1:13] <- c("Index", "First_Frame", "Number_Frames", "Frames_Missing", "X", "Y", "Precision", "Number_Photons", "Background_variance", "X_i_square", "PSF_half_width_nm", "Channel", "Z_Slice")            
roi <- paste(dir,name,"_RoiSet_t.zip", sep = "")
data_sub <- IJROI_subset(x = alldata,file = roi, pxsize = 5)

#Create statistic_file (sf)
sf <- read.csv(paste(dir,name, "_RoiSet_t_characteristics.csv", sep = ""), header = T)   
sf <- cbind(name, c(1:nrow(sf)), sf)
colnames(sf) <- c("name", "ROI_ID", "X", "Y", "BX", "BY", "Width", "Height")

#Create KDE file 
kde <- lapply(1:length(data_sub), function(x){SMOLR_KDE(data_sub[[x]], threshold = thr_kde, xlim = c(sf$BX[x], sf$BX[x]+sf$Width[x]), ylim = c(sf$BY[x], sf$BY[x]+sf$Height[x]) )} )
features <- SMOLR_FEATURES(kde, filter = "x.0.s.area",filter_value = thr_area)
kde <- apply_filtered_features(features = features, kde=kde)
catagories <- ldply(lapply(features, FUN = function(x){
  if(!length(x$parameters[x$parameters[,1]==3,3]) == "0") {d_number <- x$parameters[x$parameters[,1]==3,3]}  else {(d_number <- 0) }
  if(!length(x$parameters[x$parameters[,1]==2,3]) == "0") {r_number <- x$parameters[x$parameters[,1]==2,3]}  else {(r_number <- 0) }
  c(d_number, r_number)
}))
sf <- cbind(sf,catagories)
sf <- cbind(sf, apply(catagories,1,FUN=function(x){paste(x[1],x[2],sep="")}))
names(sf)[9:11] <- c("DMC1","RAD51","Cat")
table <- sort(table(sf$Cat), decreasing = T)

#Plot
barplot((table)/sum(table)*100, ylim = c(0,70), xlab = "DxRy", ylab = "Frequency (%)")

#Save RData files
save(data_sub, file = paste(dir,name,"_data_sub.RData", sep = "")) 
save(kde, file = paste(dir,name,"_kde.RData", sep = ""))
save(sf, file = paste(dir,name,"_sf.RData", sep = ""))


# Generate binary DxRy images ---------------------------------------------

#Load packages
library(SMoLR)
#library(EBImage)

#Select documents
dir <- "D:/"
name = "Exp1"
thr_kde = 0.15   #Threshold of KDE in localizations / nm^2

#Load data
load(paste(dir,name,"_data_sub.RData", sep = ""))  
load(paste(dir,name,"_sf.RData", sep = "")) #sf


#Analysis
for(i in 1:length(data_sub)){
  SMOLR_KDE(data_sub[[i]], threshold = 0.15, output = "tiff",file = paste(dir,name, "_binair_", i, ".tif", sep = ""), 
            xlim = c(sf$BX[i], sf$BX[i]+sf$Width[i]),         
            ylim = c(sf$BY[i], sf$BY[i]+sf$Height[i]) )
}


