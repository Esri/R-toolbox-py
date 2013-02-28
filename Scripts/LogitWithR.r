#### Load Libraries ####
print("Loading Libraries....")
library(maptools)  		
library(sm)

#### Check Version for Required "lrm" Package ####
checkRVersion = function(checkMajor, checkMinor){
    majorBool = as.numeric(R.version$major) >= checkMajor
    minorBool = as.numeric(R.version$minor) >= checkMinor
    majorBool & minorBool
    }

versionBool = checkRVersion(2, 14)
if (versionBool){
    library(rms)
}else{
    require(Design)
}

#### Get Arguments ####
Args = commandArgs()
inputFC = sub(".shp", "", Args[5], ignore.case = TRUE)
outputFC = sub(".shp", "", Args[6], ignore.case = TRUE)
dependentVar = Args[7]
independentVarString = Args[8]
usePenalty = as.integer(Args[9])
usePenalty = usePenalty == 1
coefTable = sub(".dbf", "", Args[10], ignore.case = TRUE)
diagTable = sub(".dbf", "", Args[11], ignore.case = TRUE)
print(paste(commandArgs(), collapse=" "))

#### Get Ind Var Names ####
independentVars = strsplit(independentVarString, ";")
independentVars = c(unlist(independentVars))

### Make Formula ####
form = as.formula(paste(dependentVar, paste(independentVars, collapse='+'), sep='~'))

print("Begin Calculations....")
### Using Maptools ####
shp = readShapeSpatial(inputFC)

### Do Logit ####
print("Logit....")
#fit = lrm(form, shp) 
fit = lrm(form, shp, x = TRUE, y = TRUE) 

print("Adjustment....")
### AIC Alternative Measure of Model Performance ####
bf = pentrace(fit, seq(.2,1,by=.05))
if (usePenalty) {
    pen = bf$penalty
} else {
    pen = 0.0
}

allPens = bf$results.all[,1]
allAICs = bf$results.all[,3]
for (i in 1:length(allPens)){ 
    penValue = allPens[i]
    if (penValue == pen){
        aic = allAICs[i]
        }
    }

if (usePenalty){
    fit = update(fit, penalty = bf$penalty)
    }

### Residuals ####
res = residuals.lrm(fit)
resOut = c(res)
resSTD = (resOut - mean(resOut)) / sqrt(var(resOut))

### Create Output Shape File ####
print("Writing Output....")
shp$Residual = resOut
shp$StdResid = resSTD
writeSpatialShape(shp, outputFC)

### Write Coefficient DBF Table ####
allIndVars = c("Intercept")
allIndVars = append(allIndVars, independentVars)
k = length(allIndVars)
d = matrix(0, k, 4)
d[,1] = fit$coefficients
d[,2] = sqrt(diag(fit$var))
d[,3] = d[,1] / d[,2]
d[,4] = pnorm(abs(d[,3]), lower.tail = FALSE) * 2.0
coefList = list("Variable" = allIndVars, "Coef" = d[,1], 
               "StdError" = d[,2], "Wald" = d[,3], 
               "Prob" = d[,4])
coefFrame = data.frame(coefList)
write.dbf(coefFrame, coefTable)

### Write Diagnostic DBF Table ####
diagNames = names(fit$stats)
allStats = c(as.vector(fit$stats), pen, aic)
diagValues = matrix(allStats, length(diagNames), 1)
diagList = list("Diag_Name" = diagNames, "Diag_Value" = diagValues)
diagFrame = data.frame(diagList)
write.dbf(diagFrame, diagTable)

print("Calculations Complete...")
