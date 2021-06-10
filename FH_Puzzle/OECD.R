##########################################################################
# Figure 8: London network for inverse distance and 7 nearest neighbours #
# (weight matrices)                                                      #
##########################################################################

library("rworldmap")
library("R.matlab")

#############load data and weight matrices #########

FHD <- readMat("Data/FHD.mat")
Wn <-  readMat("Data/matrices.mat")

#### plot the ISO3 names of the countries in red
theCountries <- c('AUS',
                  'AUT',
                  'BEL',
                  'CAN',
                  'CHE',
                  'DNK',
                  'ESP',
                  'FIN',
                  'FRA',
                  'GBR',
                  'GRC',
                  'IRL',
                  'ISL',
                  'ITA',
                  'JPN',
                  'KOR',
                  'MEX',
                  'NLD',
                  'NOR',
                  'NZL',
                  'PRT',
                  'SWE',
                  'TUR',
                  'USA')

# malDF is a data.frame with the ISO3 country names plus a variable to
# merge to the map data
ocedDF <- data.frame(country = c('AUS',
                                'AUT',
                                'BEL',
                                'CAN',
                                'CHE',
                                'DNK',
                                'ESP',
                                'FIN',
                                'FRA',
                                'GBR',
                                'GRC',
                                'IRL',
                                'ISL',
                                'ITA',
                                'JPN',
                                'KOR',
                                'MEX',
                                'NLD',
                                'NOR',
                                'NZL',
                                'PRT',
                                'SWE',
                                'TUR',
                                'USA'),
                     OECD= rep("OECD",24))


#### join malDF data.frame to the country map data
malMap <- joinCountryData2Map(ocedDF, joinCode = "ISO3",
                              nameJoinColumn = "country")
#### London Inverse Distance network

new_world <- subset(malMap, continent != "Antarctica")
mal <- mapCountryData(new_world, nameColumnToPlot="OECD", catMethod = "categorical",
                      missingCountryCol = NA,
                      addLegend=F, mapTitle = " ",xlim=c(-153,190),ylim= range(FHD$lat))

points(FHD$long, FHD$lat,pch=16, cex = .6)

###### w.o. map #####

plot(FHD$long, FHD$lat,pch=16, cex = .6, ylab= " ", xlab = " ",frame.plot = FALSE, yaxt = "n", xaxt = "n")

for(j in 1:24){
    segments(FHD$long[10], FHD$lat[10],FHD$long[j], FHD$lat[j],col="blue", lwd=0.6)
}
text(FHD$long+1, FHD$lat+2,
     labels = c('Sydney','Vienna','Brussels','Ottawa','Bern','Copenhagen','Madrid','Helsinki','Paris','London','Athens','Dublin','Reykjavík','Rome','Tokyo','Seoul','Mexico City','Amsterdam','Oslo','Wellington','Lisbon','Stockholm','Ankara','Washington, D.C.'),
     cex=0.5)


#### London 7NN network

mal <- mapCountryData(malMap, nameColumnToPlot="OECD", catMethod = "categorical",
                      missingCountryCol = gray(.8),addLegend=F, mapTitle = " ",
                      xlim=c(-6, 20), ylim=c(42.5, 70))
points(FHD$long, FHD$lat,pch=16, cex = 1)


W.a <- as.matrix(Wn$W)
a <- which(W.a[10,]!=0)
plot(FHD$long[a], FHD$lat[a],pch=16, cex = .6, ylab= " ", xlab = " ",frame.plot = FALSE, yaxt = "n", xaxt = "n")
points(FHD$long[10], FHD$lat[10],pch=16, cex = .6)
for(j in 1:7){
  segments(FHD$long[10], FHD$lat[10],FHD$long[a[j]], FHD$lat[a[j]],col = "blue")
}
text(FHD$long[10], FHD$lat[10]+1,"London",cex=0.6)
text(FHD$long[3]+1.6, FHD$lat[3],"Brussels",cex=0.6)
text(FHD$long[5]+1, FHD$lat[5],"Bern",cex=0.6)
text(FHD$long[6]-2, FHD$lat[6]+1,"Copenhagen",cex=0.6)
text(FHD$long[9], FHD$lat[9]-1,"Paris",cex=0.6)
text(FHD$long[12]+1, FHD$lat[12]+1,"Dublin",cex=0.6)
text(FHD$long[18]+2.1, FHD$lat[18],"Amsterdam",cex=0.6)
text(FHD$long[19]+1, FHD$lat[19],"Oslo",cex=0.6)
#text(FHD$long[a], FHD$lat[a]-1,labels = c('Brussels','Bern','Copenhagen','Paris','Dublin','Amsterdam','Oslo'),cex = 0.6)


######### 

plot(FHD$long, FHD$lat,pch=16, cex = .6, ylab= " ", xlab = " ",frame.plot = FALSE, yaxt = "n", xaxt = "n")
for (i in 1:24) {
  a <- which(W.a[i,]!=0)
  for (j in 1:7) {
    segments(FHD$long[i], FHD$lat[i],FHD$long[a[j]], FHD$lat[a[j]],col = "blue",lwd = 0.6)
  }
}

####### 7NN Istanbul
a <- c(which(W.a[23,]!=0),23)
plot(FHD$long[a], FHD$lat[a],pch=16, cex = .6, ylab= " ", xlab = " ",frame.plot = FALSE, yaxt = "n", xaxt = "n")
#points(FHD$long[23], FHD$lat[23],pch=16, cex = .6)
for(j in 1:7){
  segments(FHD$long[23], FHD$lat[23],FHD$long[a[j]], FHD$lat[a[j]],col = "blue")
}
text(FHD$long[23], FHD$lat[23]-1,"Istanbul",cex=0.6)
text(FHD$long[a[1]], FHD$lat[a[1]]+1,"Vienna",cex=0.6)
text(FHD$long[a[2]], FHD$lat[a[2]]+1,"Bern",cex=0.6)
text(FHD$long[a[3]], FHD$lat[a[3]]+1,"Copenhagen",cex=0.6)
text(FHD$long[a[4]]+1, FHD$lat[a[4]]+0.6,"Helsinki",cex=0.6)
text(FHD$long[a[5]], FHD$lat[a[5]]+1,"Athens",cex=0.6)
text(FHD$long[a[6]], FHD$lat[a[6]]+1,"Rome",cex=0.6)
text(FHD$long[a[7]], FHD$lat[a[7]]+1,"Stockholm",cex=0.6)