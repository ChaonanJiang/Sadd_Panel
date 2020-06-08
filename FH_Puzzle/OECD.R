library(rworldmap)
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
# These are the ISO3 names of the countries you'd like to plot in red

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
# malDF is a data.frame with the ISO3 country names plus a variable to
# merge to the map data

malMap <- joinCountryData2Map(ocedDF, joinCode = "ISO3",
                              nameJoinColumn = "country")
# This will join your malDF data.frame to the country map data
new_world <- subset(malMap, continent != "Antarctica")
mal <- mapCountryData(new_world, nameColumnToPlot="OECD", catMethod = "categorical",
                      missingCountryCol = gray(.8),addLegend=F, mapTitle = " ",xlim=c(-153,190),ylim= range(FHD$lat))

points(FHD$long, FHD$lat,pch=16, cex = .6)



#### London
points(FHD$long[5], FHD$lat[5],pch=16, cex = .6)
#text(FHD$long[10], FHD$lat[10]+20,"London")
for(j in 1:24){
    segments(FHD$long[5], FHD$lat[5],FHD$long[j], FHD$lat[j],col="blue")
}
a <- which(W.a[10,]!=0)
for(j in 1:7){
  segments(FHD$long[10], FHD$lat[10],FHD$long[a[j]], FHD$lat[a[j]],col = "blue")
}
text(FHD$long[10], FHD$lat[10]+2,"London")


plot(FHD$long, FHD$lat,pch=16, cex =1,xlab = " ", ylab = " ", main = "Inverse distance",axes=FALSE,col= "red")
for(i in 1:24){
  for(j in (i+1):24){
    segments(FHD$long[i], FHD$lat[i],FHD$long[j], FHD$lat[j])
  }
  
}

plot(FHD$long, FHD$lat,pch=16, cex =1,xlab = " ", ylab = " ", main = "7 nearest neighbours",axes=FALSE,col= "red")

W.a <- as.matrix(Wn$W)
for(i in 1:24){
  a <- which(W.a[i,]!=0)
  for(j in 1:7){
    segments(FHD$long[i], FHD$lat[i],FHD$long[a[j]], FHD$lat[a[j]])
  }
  
  
}








#######
theCountries <- c(
                  'AUT',
                  'BEL',
                 
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
                  
                  'NLD',
                  'NOR',
                  'NZL',
                  'PRT',
                  'SWE',
                  'TUR'
                  )
# These are the ISO3 names of the countries you'd like to plot in red

ocedDF <- data.frame(country = c(
                                 'AUT',
                                 'BEL',
                                
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
                                
                                 'NLD',
                                 'NOR',
                                 'NZL',
                                 'PRT',
                                 'SWE',
                                 'TUR'
                                ),
                     OECD= rep("OECD",18))
# malDF is a data.frame with the ISO3 country names plus a variable to
# merge to the map data

malMap <- joinCountryData2Map(ocedDF, joinCode = "NAME",
                              nameJoinColumn = "country")
# This will join your malDF data.frame to the country map data

mal <- mapCountryData(malMap, nameColumnToPlot="OECD", catMethod = "categorical",
                      missingCountryCol = gray(.8),addLegend=F, mapTitle = " ",
                      xlim=c(-6, 20), ylim=c(42.5, 70))
points(FHD$long, FHD$lat,pch=16, cex = 1)

for(i in c(2,3,5:14,18,19,21:23)){
  a <- which(W.a[i,]!=0)
  for(j in 1:7){
    segments(FHD$long[i], FHD$lat[i],FHD$long[a[j]], FHD$lat[a[j]],col="blue")
  }
  
  
}

