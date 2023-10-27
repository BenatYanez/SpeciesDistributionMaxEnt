library(terra)
library(geodata)
library(predicts)
library(rJava)
occdata <- geodata::sp_occurrence("Gypaetus", "barbatus*", geo=FALSE,removeZeros=TRUE,start=1,end=10000)
occdata[1:10,]

tiff_files <- list.files(path="C:/Users/benat/OneDrive/Documents/Biology/MastersProgramme/BiodiversityUnderPressure/Practical4-SDMPressence/data/wc2.1_10m",pattern = "\\.tif$", full.names = TRUE)
raster_list <- lapply(tiff_files, rast)
predictors <- c(raster_list[[1]],raster_list[[12]],raster_list[[13]],raster_list[[14]],raster_list[[15]],raster_list[[16]],raster_list[[17]],raster_list[[18]],raster_list[[19]],raster_list[[2]],raster_list[[3]],raster_list[[4]],raster_list[[5]],raster_list[[6]],raster_list[[7]],raster_list[[8]],raster_list[[9]],raster_list[[10]],raster_list[[11]])

plot(predictors$wc2.1_10m_bio_1, xlim=c(-180,180),ylim=c(-80,80),col="light yellow",border="light gray")
points(occdata$lon, occdata$lat,col="blue",pch=20)
#Remove duplicated data
dups <- duplicated(occdata[, c('lon', 'lat')])
sum(dups)
occ  <- occdata[!dups,]
#The distribution makes sense, there seems to be sampling bias with more samples in Europe than Africa
summary(occ$lon)
summary(occ$lat)

e <- ext(-10, 120, -35, 67)
predictors <- crop(predictors, e)
names(predictors) <- substring(names(predictors),11,16) #Shorten names of predictors
plot(predictors,1:9)

plot(predictors,1)
points(occ$lon,occ$lat, col='blue',pch=16)

#Generate background data and sample the background data
bg<-spatSample(predictors,5000,"random", na.rm=TRUE, as.points=TRUE,ext=e)
#Plot background points on a map of variable 1
plot(predictors,1)
points(bg,cex=0.1)
#Match climate and occurrence data
occlatlon<-cbind(occ$lon,occ$lat)
presvals <- extract(predictors, occlatlon) #Climate data for the regions where the species is present 
backvals <- values(bg) #The climate data for the background
pb <- c(rep(1, nrow(presvals)), rep(0, nrow(backvals))) #First collumn of the dataset is a vector of 1 for pressence and 0 for background
sdmdata  <- data.frame(cbind(pb,rbind(presvals,backvals)))
#examine how correlated the predictors are
pairs(sdmdata[,2:5], cex=0.1) #bio3 and bio 4 appear highly correlated

#Apply the Maxent approach
model<-MaxEnt(sdmdata[,-1],sdmdata[,1],removeDuplicates=TRUE)
#The model ahs a bad fit, probably because climatic variables are not what limits teh range of the vulture, 
#Their habitat ranges is limited by mountain ranges so tehya re only found above 2000m elevation, so the model is poor at predicting their location
#The more important variables are bio1,bio3, bio19 Howerever it doesnt say how good the fit is or what shape it has, 
#AUC was 0.733
plot(model)

predictedocc <- predict(model, predictors, args=c("outputformat=raw")) 
par(mfrow=c(2,1))
plot(predictedocc)
plot(predictedocc)
points(occlatlon,pch=".")
#Download future climate data
#bio_fut  <- cmip6_world(model="ACCESS-ESM1-5",ssp="245",time="2041-2060",var="bioc",download=F,res=10,path="C:/Users/benat/OneDrive/Documents/Biology/MastersProgramme/BiodiversityUnderPressure/Practical4-SDMPressence/data/")
#fut_predictors<-crop(bio_fut,e)
bio_fut<-rast("data/wc2.1_10m_bioc_ACCESS-ESM1-5_ssp245_2041-2060.tif")

dim(bio_fut)

fut_predictors<-crop(bio_fut,e)

plot(predictors,2)
names(fut_predictors)<-names(predictors)
#Generate future predictions
futpredictedocc <- predict(model, fut_predictors, args=c("outputformat=raw")) 
par(mfrow=c(2,1))
plot(predictedocc,main="current")
plot(futpredictedocc,main="Future") #The climate suitability of the range appears to decrease

#Run for different species
occdata <- geodata::sp_occurrence("Antilocapra", "americana*", geo=FALSE,removeZeros=TRUE,start=1,end=10000)
occdata[1:10,]
tiff_files <- list.files(path="C:/Users/benat/OneDrive/Documents/Biology/MastersProgramme/BiodiversityUnderPressure/Practical4-SDMPressence/data/wc2.1_10m",pattern = "\\.tif$", full.names = TRUE)
raster_list <- lapply(tiff_files, rast)
predictors <- c(raster_list[[1]],raster_list[[12]],raster_list[[13]],raster_list[[14]],raster_list[[15]],raster_list[[16]],raster_list[[17]],raster_list[[18]],raster_list[[19]],raster_list[[2]],raster_list[[3]],raster_list[[4]],raster_list[[5]],raster_list[[6]],raster_list[[7]],raster_list[[8]],raster_list[[9]],raster_list[[10]],raster_list[[11]])

plot(predictors$wc2.1_10m_bio_1, xlim=c(-180,180),ylim=c(-80,80),col="light yellow",border="light gray")
points(occdata$lon, occdata$lat,col="blue",pch=20)
occdata<-subset(occdata,lon<0)
occdata<-subset(occdata,lat<80)

dups <- duplicated(occdata[, c('lon', 'lat')])
sum(dups)
occ  <- occdata[!dups,]
#The distribution makes sense, there seems to be sampling bias with more samples in Europe than Africa
summary(occ$lon)
summary(occ$lat)
e <- ext(-130, -50, 20, 90)

predictors <- crop(predictors, e)

names(predictors)<-substring(names(predictors),11,16)
plot(predictors,1:9)

plot(predictors,1)
points(occ$lon,occ$lat, col='blue',pch=16)

#Generate background data and sample the abckground data
bg<-spatSample(predictors,5000,"random", na.rm=TRUE, as.points=TRUE,ext=e)
#Plot background points on a map of variable 1
plot(predictors,1)
points(bg,cex=0.1)
#Match climate and occurrence data
occlatlon<-cbind(occ$lon,occ$lat)
presvals <- extract(predictors, occlatlon) #Climate data for the regions where the species is present 
backvals <- values(bg) #The climate data for the background
pb <- c(rep(1, nrow(presvals)), rep(0, nrow(backvals))) #First collumn of the dataset is a vector of 1 for pressence and 0 for background
sdmdata  <- data.frame(cbind(pb,rbind(presvals,backvals)))
model<-MaxEnt(sdmdata[,-1],sdmdata[,1],removeDuplicates=TRUE)
#The model ahs a bad fit, probably because climatic variables are not what limits teh range of the vulture, 
#Their habitat ranges is limited by mountain ranges so tehya re only found above 2000m elevation, so the model is poor at predicting their location
#The more important variables are bio1,bio3, bio19 Howerever it doesnt say how good the fit is or what shape it has, 
#AUC was 0.733
plot(model)
predictedocc <- predict(model, predictors, args=c("outputformat=raw")) 
par(mfrow=c(2,1))
plot(predictedocc)
plot(predictedocc)
points(occlatlon,pch=".")
#Download future climate data
#bio_fut  <- cmip6_world(model="ACCESS-ESM1-5",ssp="245",time="2041-2060",var="bioc",download=F,res=10,path="C:/Users/benat/OneDrive/Documents/Biology/MastersProgramme/BiodiversityUnderPressure/Practical4-SDMPressence/data/")
#fut_predictors<-crop(bio_fut,e)
bio_fut<-rast("data/wc2.1_10m_bioc_ACCESS-ESM1-5_ssp245_2041-2060.tif")

dim(bio_fut)

fut_predictors<-crop(bio_fut,e)

plot(predictors,2)
names(fut_predictors)<-names(predictors)
#Generate future predictions
futpredictedocc <- predict(model, fut_predictors, args=c("outputformat=raw")) 
par(mfrow=c(2,1))
plot(predictedocc,main="current")
plot(futpredictedocc,main="Future")