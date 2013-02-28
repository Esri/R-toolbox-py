normalize = function(x){
    normX = (x - mean(x)) / sqrt(var(x))
    }

#### Get/Parse Arguments ####
Args = commandArgs()
inputFC = sub(".shp", "", Args[5], ignore.case = TRUE)
outputFC = sub(".shp", "", Args[6], ignore.case = TRUE)
numClusters = as.integer(Args[7])
clusterMethod = Args[8]
fields = Args[9]
useLocation = as.integer(Args[10])
useFields = FALSE

#### Create Field Names is Applicable ####
if (fields != "NA"){
    varNames = strsplit(fields, ";")
    varNames = c(unlist(varNames))
    useFields = TRUE
    }

#### Import the Cluster Library ####
print("Loading Libraries....")
require(clustTool)

#### Using Maptools For Shapefiles ####
library(maptools)  	

print("Begin Calculations....")
shp = readShapePoints(inputFC)
centroids = coordinates(shp)

#### Create Boolean For Using Location ####
#### Set to True: if useLocation requested or no fields provided ####
noFields = (useFields == FALSE)
addLocation = (useLocation | noFields)

#### Add Locations ####
if (addLocation){
    XCOO = normalize(centroids[,1])
    YCOO = normalize(centroids[,2])
    allVars = cbind(XCOO, YCOO)
    } else {
    allVars = NULL
    }

#### Add Fields ####
if (useFields){
    for (i in 1:length(varNames)){
        allVars = cbind(allVars, normalize(shp[[varNames[i]]]))
        }
    } 

#### Create Data Frame for Analysis ####
new = data.frame(allVars)

#### Run Algorithm ####
newCL = clust(new, k = numClusters, method = clusterMethod)

#### Write Output to New FC ####
shp$CLUSTER = newCL$cluster
writeSpatialShape(shp, outputFC)

#### Print Output Summary ####
print(newCL)
delim = paste(rep('-', 55), collapse = "")
print(delim)
print(names(newCL))
print(delim)
print(delim)
print("Centers")
print(newCL$centers)
print(delim)
print(delim)
print(newCL$valMeasures)
print(delim)

print("Calculations Complete...")
