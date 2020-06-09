##########################################################################
# Figure 8: London network for inverse distance and 7 nearest neighbours #
# (weight matrices)                                                      #
##########################################################################

library("rworldmap")
library("R.matlab")

#############load data and weight matrices #########

FHD <- readMat("FHD.mat")
Wn <-  readMat("matrices.mat")

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
                      missingCountryCol = gray(.8),addLegend=F, mapTitle = " ",xlim=c(-153,190),ylim= range(FHD$lat))

points(FHD$long, FHD$lat,pch=16, cex = .6)

for(j in 1:24){
    segments(FHD$long[10], FHD$lat[10],FHD$long[j], FHD$lat[j],col="blue")
}


#### London 7NN network

mal <- mapCountryData(malMap, nameColumnToPlot="OECD", catMethod = "categorical",
                      missingCountryCol = gray(.8),addLegend=F, mapTitle = " ",
                      xlim=c(-6, 20), ylim=c(42.5, 70))
points(FHD$long, FHD$lat,pch=16, cex = 1)
W.a <- as.matrix(Wn$W)
a <- which(W.a[10,]!=0)
for(j in 1:7){
  segments(FHD$long[10], FHD$lat[10],FHD$long[a[j]], FHD$lat[a[j]],col = "blue")
}
text(FHD$long[10], FHD$lat[10]+2,"London")

