library(rgbif)
library(sdm)
library(usdm)
library(CoordinateCleaner)
library(maptools)

# number of occurrences reported for Malawi:
occ_count(country = isocodes[grep("Malawi", isocodes$name), "code"])

# number of observations reported for Malawi:
occ_count(country = isocodes[grep("Malawi", isocodes$name), "code"],
          basisOfRecord = 'OBSERVATION')

# Get species data for Black rhinoceros. Check for synonyms
name_suggest(q = 'Diceros bicornis', rank = 'species')

# Get species data for Cheetah
name_suggest(q = 'Acinonyx jubatus', rank = 'species')

# Get species data for fresh water snail
name_suggest(q = 'Bulinus nyassanus', rank = 'species')

# Check number of records - here filtered to those with coordinate information
occ_search(scientificName = "Diceros bicornis", hasCoordinate = TRUE, limit = 10)

occ_search(scientificName = "Acinonyx jubatus", hasCoordinate = TRUE, limit = 10)

data(wrld_simpl)
plot(wrld_simpl, xlim = c(26, 37), ylim = c(-19, -8)) # Extent for Malawi
plot(wrld_simpl, xlim = c(-12, 40), ylim = c(-30, 25)) # Extent for Africa

points(gbif_cheetah$decimalLongitude, gbif_cheetah$decimalLatitude, 
       col = 'red',  pch = 19)
cheetah_location <- sf::st_as_sf(gbif_cheetah,
                                 coords = c("decimalLongitude", "decimalLatitude"),
                                 crs = 4326, agr = "constant")

mapview::mapview(cheetah_location)



# Some coordinates appearing in North America. let's clean them
library(CoordinateCleaner)

# We use only those data entries with coordinate information - Note that you 
# don't need this if you have used the hasCoordinate=T in the occ_search() function:
gbif_cheetah <- subset(gbif_cheetah, !is.na(decimalLatitude))

# We now clean the coordinates and check for outliers - see ?clean_coordinates for more options
cl_gbif_cheetah <- clean_coordinates(gbif_cheetah, lon = "decimalLongitude", 
                                     lat = "decimalLatitude", countries = "countryCode", 
                                     tests = c("centroids", "outliers", "duplicates", "institutions"), 
                                     inst_rad = 10000)

points(gbif_cheetah$decimalLongitude, gbif_cheetah$decimalLatitude, 
       col = 'red', pch = 19)

points(gbif_cheetah$decimalLongitude[cl_gbif_cheetah$.summary], 
       gbif_cheetah$decimalLatitude[cl_gbif_cheetah$.summary], 
       col = 'blue', pch = 18)



gbif_cheetah <- occ_search(scientificName = "Acinonyx jubatus", 
                           hasCoordinate = TRUE, limit = 500) 

# We are just interested in the data frame containing the records
gbif_cheetah <- gbif_cheetah$data


# Create new column 
gbif_cheetah$presence <- 1

# Drop other columns
gbif_cheetah_trim <- gbif_cheetah[,c('decimalLongitude', 'decimalLatitude', 'presence')]

# Convert to spatial object
coordinates(gbif_cheetah_trim) <- ~decimalLongitude + decimalLatitude

mapview::mapview(gbif_cheetah_trim)

# Get predictor variables
# https://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/bio_10m_bil.zip
bio <- raster::getData('worldclim', var = 'bio', res = 10)

# Test collinear variables
v1 <- vifstep(bio)

v2 <- vifcor(bio, th = 0.9)

biom <- exclude(bio, v2)

# Map
plot(biom[[1]])
points(gbif_cheetah, cex=.5)

# Assign projection
# Raster () function without any argument creates an empty raster with a 
# geographic coordinate system
proj4string(gbif_cheetah_trim) <- projection(raster())
mapview::mapview(gbif_cheetah_trim)

presence_data <- sdmData(presence~., gbif_cheetah_trim, predictors = biom)
presence_data

nrow(gbif_cheetah_trim) # number of presence points

# > presence_data
# class                                 : sdmdata 
# =========================================================== 
# number of species                     :  1 
# species names                         :  presence 
# number of features                    :  12 
# feature names                         :  bio2, bio3, bio4, ... 
# type                                  :  Presence-Only 
# has independet test data?             :  FALSE 
# number of records                     :  495 
# has Coordinates?                      :  TRUE 

# Generate pseudo absence
presence_absence_data <- sdmData(presence~., gbif_cheetah_trim, predictors = biom,
                                 bg = list(n= 1000), method = 'gRandom') # 1000 random locations

presence_absence_data <- sdmData(presence~., gbif_cheetah_trim, predictors = biom,
                                 bg = list(n= 1000))

model <- sdm(presence~.,presence_absence_data, methods = c('glm','gbm','svm','rf', 'brt', 'mars'),
             replication = c('boot'), n = 2)

# getmethodNames()

predicted <- predict(model, biom, 'predictions.img')

gui(model)
plot(predicted[[c(1, 3, 5, 7, 9)]])
plot(predicted[[c(1)]])

ensembled <- ensemble(model, biom,  'ensemble.img',          # current
                      setting = list(methods = 'weights', 
                                      stat = "AUC",
                                      opt = 2))

ensembled2 <- ensemble(model, biom,  'ensemble.img', 
                       setting = list(methods = 'weights', 
                                      stat = "TSS",
                                      opt = 2))

biof <- raster::getData('CMIP5', var = 'bio', res = 10, rcp = 85, year = 70, model = 'AC')
biof2 <- biof[[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)]]

names(biof2)

names(biof2) <- names(biom) # match names

names(biof2)

pf <- predict(model, biof2, 'predictionsf.img')

ensembled_future <- calc(pf, mean)

ensembled_future_weighted <- ensemble(model, biof2, 'ensembled future.img',
                setting = list(method = 'weights', stat = "TSS", opt = 2))


plot(stack(ensembled_future, ensembled2))

mapview::mapview(stack(ensembled_future, ensembled2))

# Calculate change numerically
change <- ensembled_future - ensembled2

plot(change)

mapview::mapview(change)

getEvaluation(model, stat = c('AUC', 'TSS', 'threshold'), opt = 2)