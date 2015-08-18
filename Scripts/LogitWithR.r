#### Load Libraries ####
print("Loading Libraries....")
library(maptools)  		
library(sm)
library(foreign)

#### Check Version for Required "lrm" Package ####
checkRVersion = function(checkMajor, checkMinor){
    minorBool = as.numeric(R.version$minor) >= checkMinor	#will work until v3.14
    majorBool = as.numeric(R.version$major) == checkMajor
	if (!minorBool)	{
		majorBool = as.numeric(R.version$major) > checkMajor
		if (majorBool)	{
			minorBool = TRUE
			}
		}
    if (majorBool & minorBool)	{
		library(rms)
		}
	else	{
		library(Design)
		}
    }

versionBool = checkRVersion(2, 14)

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
res = resid(fit)
resOut = c(res)
resSTD = (resOut - mean(resOut)) / sqrt(var(resOut))

### Create Output Shape File ####
print("Writing Output....")
shp$Residual = resOut
shp$StdResid = resSTD
writeSpatialShape(shp, outputFC)

### Write Coefficient DBF Table ####
coefFrame <- data.frame(c("Intercept",independentVars))
names(coefFrame)[1] <- "Coefficients"
coefFrame <- cbind(coefFrame,coef(summary(fit)))
write.dbf(coefFrame, coefTable)

### Write Diagnostic DBF Table ####
diagNames = names(fit$stats)
allStats = c(as.vector(fit$stats), pen, aic)
diagValues = matrix(allStats, length(diagNames), 1)
diagList = list("Diag_Name" = diagNames, "Diag_Value" = diagValues)
diagFrame = data.frame(diagList)
write.dbf(diagFrame, diagTable)

print("Calculations Complete...")
