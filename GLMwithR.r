#### Load Libraries ####

library(maptools)  		
library(sm)
require(foreign)
require(MASS)
require(GGally)
require(texreg)

#### Get Arguments ####
Args = commandArgs()
inputFC = sub(".shp", "", Args[5], ignore.case = TRUE)
outputFC = sub(".shp", "", Args[6], ignore.case = TRUE)
dependentVar = Args[7]
independentVarString = Args[8]
coefTable = sub(".dbf", "", Args[9], ignore.case = TRUE)
diagTable = sub(".dbf", "", Args[10], ignore.case = TRUE)
print(paste(commandArgs(), collapse=" "))

#### Get Ind Var Names ####
independentVars = strsplit(independentVarString, ";")
independentVars = c(unlist(independentVars))

### Make Formula ###
form = as.formula(paste(dependentVar, paste(independentVars, collapse='+'), sep='~'))

print("Begin Calculations....")
### Using Maptools ###
shp = readShapeSpatial(inputFC)

### Check Linear Relationships Between Dependent and Independent Variables ###
shp_tab <- as.data.frame(shp)
var_tab <- data.frame(shp_tab[dependentVar])
for (i in 1:length(independentVars))	{
	var_tab <- cbind(var_tab,shp_tab[independentVars[i]])
	}
png("nb_cov_test.png")
ggpairs(var_tab)
dev.off()

### Fit Negative-Binomial GLM ###
var_tab[,1] <- as.integer(var_tab[,1])
nb = glm.nb(form, var_tab) 

### Print Table of Results for Coefficents ###
fit_tab <- extract.glm(nb)
png('nb_coefs.png')
plotreg(nb, custom.model.names="Neg. Binomial", custom.note="Darker bars: 1 standard error.  Lighter bars: 95% CI.")
dev.off()

### Check Model Assumptions of Normality and Independence ###
shp$res <- resid(nb)

png("nb_norm_test.png")
qqnorm(shp$res)
qqline(shp$res)
dev.off()
png("nb_indep_test.png")
qplot(shp_tab[,dependentVar],shp$res)
dev.off()

### Residuals ####
res = resid(nb)
resOut = c(res)
resSTD = (resOut - mean(resOut)) / sqrt(var(resOut))
nb$var <- var(resOut)

### Create Output Shape File ####
shp$Residual = resOut
shp$StdResid = resSTD
writeSpatialShape(shp, outputFC)

### Write Coefficient DBF Table ####

coefFrame <- data.frame(c("Intercept",independentVars))
names(coefFrame)[1] <- "Coefficients"
coefFrame <- cbind(coefFrame,coef(summary(nb)))
write.dbf(coefFrame, coefTable)

### Calculate AIC and BIC for Candidate Models ###

### Make list of all unique combinations using ridiculous iterative algorithm ###

indepvars <- c(2:(length(independentVars)+1))	# column number for each indep var in var_tab
N <- 0													#initialize variables for use in function
cb <- list()
mb <- list()
cbi <- list()
check <- list()
save <- list()
flag <- 0
unique <- 0
vmb <- 0
vcb <- 0

comb <- function(	indepvars=indepvars){
	N <- length(indepvars)
	cb <- list()

for (k in 1:N)	{										#size of combo
	if (k == 1)	{
		for (i in 1:N)	{								#for each member
			mb <- append(mb,indepvars[i])				#add member to combo
			cb <- append(cb,mb)							#new combo
			save <- append(save,mb)
			mb <- list()									#reset member list for next combo
			}
		}
	if (k > 1)	{
		check <- save
		save <- list()
		for (i in 1:length(check))	{
			for (j in 1:N)	{
				unique <- 1								#assume new combo is unique
				if (!(indepvars[j] %in% unlist(check[i])))	{	#new member not already in combo
					mb <- append(mb, check[i])			#new combo includes old combo
					mb <- append(mb, indepvars[j])		#adds new member to combo
					if (length(save) != 0)	{
						for (m in 1:length(save))	{			#for each combo already in list
							vmb <- sort(unlist(mb))			#sort new combo and old combo for comparison
							vsave <- sort(unlist(save[m]))
							if (length(vmb)==length(vsave))	{	#for combos of equal length
								flag <- 0
								for (n in 1:length(vmb))	{	#elementwise comparison of members
									if (vmb[n]!=vsave[n])	{		#flag if candidate combo different than combo in list
										flag <- flag + 1
										}
									}
								if (flag==0)	{				#if no flags, label new combo "not unique"
									unique <- 0
									}
								}
							}
						}
					} else	{
					unique <- 0
					}
				if (unique==1)	{
					cb <- append(cb, list(mb))			#add unique combo to list
					save <- append(save, list(mb))
					}
				mb <- list()
				}
			}
		}
	k <- k + 1
	}
	return(cb)		#combo list contains all unique combinations
	}

cb <- comb(indepvars)

### Build Formulas for AIC and BIC Tests ###
v_list <- list()
varlist <- list()
for (i in 1:length(cb))	{
	for (j in 1:length(unlist(cb[i])))	{
		v_list <- append(v_list, paste("var_tab[ ,", unlist(cb[i])[j], "]", sep=''))
		}
	varlist <- append(varlist, list(v_list))
	v_list <- list()
	}

fm <- list()
formlist <- list()
for (i in 1:length(cb))	{
	fm <-	paste("var_tab[,1]","~",
			paste(c(unlist(varlist[i])), collapse = '+'),
			sep=" ")
	formlist <- append(formlist, fm)
	}

### Build Matrix of AIC and BIC Tests ###
diag_tab <- matrix(0, ncol=3, nrow=length(formlist), dimnames=list(list(),
				c("Model","AIC","BIC")))

model <- 0

for (i in 1:length(formlist))	{
	diag_tab[i,1] <- paste(independentVars[c((unlist(cb[i]))-1)], collapse = '+')
	model <- glm.nb(as.formula(paste(formlist[i])))
	diag_tab[i,2] <- AIC(model)
	diag_tab[i,3] <- BIC(model)
	}

### Write Diagnostic DBF Table ####
write.dbf(diag_tab, diagTable)

print("Calculations Complete...")

