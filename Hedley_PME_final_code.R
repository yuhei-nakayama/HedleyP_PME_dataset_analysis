#packages required
library(psych)
library(lavaan)
library(nFactors)
library(Hmisc)
library(ggplot2)
library(viridis)
library(rnaturalearth)
library(raster)
library(sp)

##Resin subset
#Main dataset used for analyses: entire dataset filtered to include only observations with Resin Pi data
#Exploratory factor analysis > Confirmatory factor analysis > structural equation modeling
#Factor analysis conducted for P fractions data only, and then PME activity data included in structural equation modeling

###Data input and description
#Data with resin (subset)
Pfrac.PME.resin <- read.csv("Pfract_resin_final.csv", header = T)
#data without PME (subsetting P fractions only)
Pfrac.PME.resin.2 <- Pfrac.PME.resin[,-7]

#Bivariate spearman correlation between variables
rcorr(as.matrix(Pfrac.PME.resin), type = "pearson")
rcorr(as.matrix(Pfrac.PME.resin), type = "spearman")


###Exploratory and confirmatory factor analysis (P fractions)
#Factor Analysis on Pfractions only

##Exploratory Factor Analysis
#extract correlation matrix from P fractions data
Pfrac.PME.resin.cor <- cor((Pfrac.PME.resin.2), use = "complete")

#Choose number of factors
ev <- eigen(Pfrac.PME.resin.cor, symmetric = T)
ap <- parallel(subject=nrow(Pfrac.PME.resin),var=6, rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)
#3 factors suggested based on Scree Plot

#factor analysis: maximum likelihood method

#with varimax rotation
FAresin <- fa(Pfrac.PME.resin.cor, nfactors = 3, fm = 'ml', SMC = T ,rotate = "varimax", scores = "regression")
FAresin
round(FAresin$residual,2)
FAresin$values

#tried oblimin rotation check if obique transformation changes the loadings, but the effect was small and therefore we proceeded with varimax rotation
FAresin.2 <- fa(Pfrac.PME.resin.cor, nfactors = 3, fm = 'ml', SMC = T ,rotate = "oblimin", scores = "regression")
FAresin.2
round(FAresin.2$residual,2)
FAresin.2$values

#initial model based on exploratory factor analysis output
#based on the factor analysis result (>=0.5 loading cutoff)
#f1 =~ Resin.Pi + Bic.Pi + acid.Pi
#f2 =~ Bic.Po + OH.Po
#f3 =~ OH.Pi

#based on the factor analysis result (>=0.4 loading cutoff)
#f1 =~ Resin.Pi + Bic.Pi + acid.Pi
#f2 =~ Bic.Po + OH.Po
#f3 =~ Bic.Pi + OH.Pi

##Confirmatory Factor Analysis
#original model with orthogonal factors (given that Varimax rotation was used), based on >=0.5 loading cutoff
mod.resin0 <- '
f1 =~ Resin.Pi + Bic.Pi + acid.Pi
f2 =~ Bic.Po + OH.Po
f3 =~ OH.Pi
f1 ~~ 0*f2
f2 ~~ 0*f3
f1 ~~ 0*f3'
cfa.resin0 <- cfa(mod.resin0, data = Pfrac.PME.resin.2, std.lv = T, std.ov = T)
#error, std error could not be computed

#original model with orthogonal factors (given that Varimax rotation was used), f2 loadings fixed to avoid error
mod.resin1 <- '
f1 =~ Resin.Pi + Bic.Pi + acid.Pi
f2 =~ a*Bic.Po + a*OH.Po
f3 =~ OH.Pi
f1 ~~ 0*f2
f2 ~~ 0*f3
f1 ~~ 0*f3'
cfa.resin1 <- cfa(mod.resin1, data = Pfrac.PME.resin.2, std.lv = T, std.ov = T)
#no error

# testing for correlated factors
mod.resin1.1 <- '
f1 =~ Resin.Pi + Bic.Pi + acid.Pi
f2 =~ a*Bic.Po + a*OH.Po
f3 =~ OH.Pi
f1 ~~ f2
f2 ~~ f3
f1 ~~ f3'
cfa.resin1.1 <- cfa(mod.resin1.1, data = Pfrac.PME.resin.2, std.lv = T, std.ov = T)
summary(cfa.resin1.1, standardized=T, fit.measure=T, rsq=T)
##no error
#covariance between f1 ~~ f2 and f2 ~~ f3 not significant

#compare uncorrelated factors (based on Varimax) vs. correlated factors
anova(cfa.resin1, cfa.resin1.1)
#model 1.1 is better than model 1 (orthogonal factors), with smaller AIC and p-value = 0.000858
#covariance between all factors supported by likelihood ratio test, but covariance between f1 ~~ f2 and f2 ~~ f3 not statistically significant

#loading of f2 fixed, only f1 ~~ f3 with f1 ~~ 0*f2 and f2 ~~ 0*f3 to remove covariances that were not statistically significant
mod.resin1.2 <- '
f1 =~ Resin.Pi + Bic.Pi + acid.Pi
f2 =~ a*Bic.Po + a*OH.Po
f3 =~ OH.Pi
f1 ~~ 0*f2
f2 ~~ 0*f3
f1 ~~ f3'
cfa.resin1.2 <- cfa(mod.resin1.2, data = Pfrac.PME.resin.2, std.lv = T, std.ov = T)
## no error warning

#compare uncorrelated factors (based on Varimax) vs. correlated factors
anova(cfa.resin1.1, cfa.resin1.2)
#cfa.resin1.2 has lower AIC, larger df, and p-value > 0.4741, so cfa.resin1.2 is supported, with only f1 ~~ f3 covarying.

#model 1.2 summary
summary(cfa.resin1.2, standardized=T, fit.measure=T, rsq=T)
resid(cfa.resin1.2)
#fit is unacceptable, CFI = 0.834, SRMR = 0.085 (cut off: CFI > 0.95, SRMR < 0.08)
#given that exploratory factor analysis shows the loading of Bic Pi onto 3rd factor is higher than the liberal cutoff of 0.40, Bic Pi included in the third factor to test the model fit.
#also, resid(cfa.resin1.2) shows large residual covariance of Bic.Pi and OH Pi (0.346).

# adding Bic.Pi to f3 to model 1.2
mod.resin1.3 <- '
f1 =~ Resin.Pi + Bic.Pi + acid.Pi
f2 =~ a*Bic.Po + a*OH.Po
f3 =~ OH.Pi + Bic.Pi
f1 ~~ 0*f2
f2 ~~ 0*f3
f1 ~~ f3'
cfa.resin1.3 <- cfa(mod.resin1.3, data = Pfrac.PME.resin.2, std.lv = T, std.ov = T)
summary(cfa.resin1.3, standardized=T, fit.measure=T, rsq=T)
parameterEstimates(cfa.resin1.3, ci=TRUE, level=0.95, boot.ci.type = "bca.simple", standardized=TRUE)
standardizedSolution(cfa.resin1.3, type="std.all")
#model fit significantly improved (CFI = 0.992, SRMR = 0.035)

#compare models (addition of Bic Pi on f3)
anova(cfa.resin1.2, cfa.resin1.3)
#model1.3 has lower AIC, and p-value < 2.2e-16, thus addition of Bic Pi to f3 supported.


###Structural Equation modeling
#Structural equation modeling including PME activity as a response variable
#model 1, adapting model 1.3 from CFA
mod.resin.sem1 <- '
f1 =~ Resin.Pi + Bic.Pi + acid.Pi
f2 =~ a*Bic.Po + a*OH.Po
f3 =~ OH.Pi + Bic.Pi
f1 ~~ 0*f2
f2 ~~ 0*f3
f1 ~~ f3
PME ~ f1 + f2 + f3'
sem.resin1 <- sem(mod.resin.sem1, data = Pfrac.PME.resin, std.ov = T, std.lv = T, orthogonal = T)
summary(sem.resin1, standardized=T, fit.measure=T, rsq=T)
#bias-corrected and accelerated bootstrap for estimating confidence intervals
parameterEstimates(sem.resin1, ci=TRUE, level=0.95, boot.ci.type = "bca.simple", standardized=TRUE)
#obtain standardized solutions
standardizedSolution(sem.resin1, type="std.all")
#check residuals
resid(sem.resin1, "cor")

#all path coefficients are significant, and no error in the model

#tried unfixing loading on f2 to allow individual estimation of Bic.Po and OH.Po --> produced negative variance (for OH.Po), so not done

#testing feedback loop: does PME increase inorganic P? or decrease organic P?
#model 2 (adding f1 ~ PME)
mod.resin.sem2 <- '
f1 =~ Resin.Pi + Bic.Pi + acid.Pi
f2 =~ a*Bic.Po + a*OH.Po
f3 =~ OH.Pi + Bic.Pi
f1 ~~ 0*f2
f2 ~~ 0*f3
f1 ~~ f3
PME ~ f1 + f2 + f3
f1 ~ PME'
sem.resin2 <- sem(mod.resin.sem2, data = Pfrac.PME.resin, std.ov = T, std.lv = T, orthogonal = T)
summary(sem.resin2, standardized=T, fit.measure=T, rsq=T)

#f1 ~ PME is not significant (p-value = 0.174)
anova(sem.resin1, sem.resin2)
#feedback loop not justified, no decrease in AIC and p-value = 0.19

#same is true for f2 ~ PME
#model 3 (adding f2 ~ PME)
mod.resin.sem3 <- '
f1 =~ Resin.Pi + Bic.Pi + acid.Pi
f2 =~ a*Bic.Po + a*OH.Po
f3 =~ OH.Pi + Bic.Pi
f1 ~~ 0*f2
f2 ~~ 0*f3
f1 ~~ f3
PME ~ f1 + f2 + f3
f2 ~ PME'
sem.resin3 <- sem(mod.resin.sem3, data = Pfrac.PME.resin, std.ov = T, std.lv = T, orthogonal = T)
summary(sem.resin3, standardized=T, fit.measure=T, rsq=T)

#f2 ~ PME is not significant (p-value = 0.195)
anova(sem.resin1, sem.resin3)
#feedback loop not justified, increase in AIC and p-value = 0.2074


##Resin subset (resin.Pi + Bic.Pi = Labile.Pi)
#what if we sum resin.Pi and Bic.Pi together?
  
###Data organization
#add column of labile Pi as sum of resin.Pi and bic.Pi, then remove resin.Pi and bic.Pi
Pfrac.PME.resin.labile <- Pfrac.PME.resin
Pfrac.PME.resin.labile$Labile.Pi <- Pfrac.PME.resin.labile$Resin.Pi + Pfrac.PME.resin.labile$Bic.Pi
Pfrac.PME.resin.labile.2 <- Pfrac.PME.resin.labile[,-c(1,2,7)]


###Exploratory factor analysis
#Exploratory Factor Analysis
#correlation matrix of untransformed data without PME
Pfrac.PME.resin.labile.cor.2 <- cor(Pfrac.PME.resin.labile.2, use = "complete")

# Choose number of factors #
ev <- eigen(Pfrac.PME.resin.labile.cor.2, symmetric = T)
ap <- parallel(subject=nrow(Pfrac.PME.resin.labile),var=5, rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)
#use of 2 factors suggested

#maximum likelihood method
FAresin.labile2 <- fa(Pfrac.PME.resin.labile.cor.2, nfactors = 2, fm = 'ml', SMC = T, rotate = "varimax", scores = "regression")
FAresin.labile2


##Resin subset (Bic.Pi + OH.Pi = Bic.OH.Pi)
#what if we sum bic.Pi and OH.Pi together?

  
###Data organization
#add column of Bic.OH Pi as sum of OH.Pi and Bic.Pi, then remove Bic.Pi and OH.Pi
Pfrac.PME.resin.Bic.OH <- Pfrac.PME.resin
Pfrac.PME.resin.Bic.OH$Bic.OH.Pi <- Pfrac.PME.resin.Bic.OH$Bic.Pi + Pfrac.PME.resin.Bic.OH$OH.Pi
Pfrac.PME.resin.Bic.OH.2 <- Pfrac.PME.resin.Bic.OH[,-c(2,4,7)]

###Exploratory factor analysis
#Exploratory Factor Analysis
#correlation matrix of untransformed data without PME
Pfrac.PME.resin.Bic.OH.cor.2 <- cor(Pfrac.PME.resin.Bic.OH.2, use = "complete")

# Choose number of factors #
ev <- eigen(Pfrac.PME.resin.Bic.OH.cor.2, symmetric = T)
ap <- parallel(subject=nrow(Pfrac.PME.resin.Bic.OH),var=5, rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)
#use of 2 factors suggested

#maximum likelihood method
FAresin.Bic.OH <- fa(Pfrac.PME.resin.Bic.OH.cor.2, nfactors = 2, fm = 'ml', SMC = T, rotate = "varimax", scores = "regression")
FAresin.Bic.OH


##Entire dataset with Labile Pi and Po
#full, entire data including all observations with labile Pi and labile Po

###Data input
#full, entire data including all observations with labile Pi and labile Po
Pfrac.PME.labile <- read.csv("Pfract_labile_final.csv", header = T)
#data without PME
Pfrac.PME.labile.2 <- Pfrac.PME.labile[,-6]

###Exploratory Factor Analysis
#Exploratory Factor Analysis
#correlation matrix
Pfrac.PME.labile.cor <- cor(Pfrac.PME.labile.2, use = "complete")

#Choose number of factors #
ev <- eigen(Pfrac.PME.labile.cor, symmetric = T)
ap <- parallel(subject=nrow(Pfrac.PME.labile),var=5, rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)
# use 2/1 factors for data without PME

#maximum likelihood method
FAlabile <- fa(Pfrac.PME.labile.cor, nfactors = 2, fm = 'ml', SMC = T, rotate = "varimax", scores = "regression")
FAlabile

##Resin subset other data description
#Resin subset dataset, including other variables (TC, pH, etc.)

###Data input
###Resin subset data including other potentially relevant variables (TC, pH, etc.)
Pfrac.PME.resin.descrip <- read.csv("Pfrac_resin_descrip_final.csv", header = T)

###Descriptive statistics
##descriptive stats of P fractions, PME activity, and other variables
psych::describe(Pfrac.PME.resin.descrip[,c(3,4,11:16,20,21,25,27,30)])

##descriptive stats of acid Pi and pH by weathering groping
psych::describeBy(Pfrac.PME.resin.descrip[,c(16,20)], Pfrac.PME.resin.descrip$Degree.Weathering)

###Bivariate correlations between variables
##pH subset
Pfrac.PME.resin.pH <- Pfrac.PME.resin.descrip[,c(11:16,20,27)] #subset pH, Pfractions, PME
Pfrac.PME.resin.pH.omit <- na.omit(Pfrac.PME.resin.pH) #omit missing observations
rcorr(as.matrix(Pfrac.PME.resin.pH.omit), type = "spearman") #Spearman correlation

##TC subset
Pfrac.PME.resin.TC <- Pfrac.PME.resin.descrip[,c(11:16,21,27)] #subset TC, Pfractions, PME
Pfrac.PME.resin.TC.omit <- na.omit(Pfrac.PME.resin.TC) #omit missing observations
rcorr(as.matrix(Pfrac.PME.resin.TC.omit), type = "spearman") #Spearman correlation

###C:Po ratio
##TC subset with TC:Po
Pfrac.PME.resin.TC.Po <- Pfrac.PME.resin.TC.omit #assign another dataframe
Pfrac.PME.resin.TC.Po$Po <- Pfrac.PME.resin.TC.Po$Bic.Po + Pfrac.PME.resin.TC.Po$OH.Po #add column of Po = Bic.Po + OH.Po
Pfrac.PME.resin.TC.Po$TC.Po <- 10000*(Pfrac.PME.resin.TC.Po$TC)/Pfrac.PME.resin.TC.Po$Po
rcorr(as.matrix(Pfrac.PME.resin.TC.Po), type = "spearman") #Spearman correlation

#descriptive statistics
psych::describe(Pfrac.PME.resin.TC.Po$TC.Po)

#test correlations for C:Po <100, <200, (100-200), 200-300, and >300 separately to see if threshold applies. Some threshold also tested with extreme values removed
a <- subset(Pfrac.PME.resin.TC.Po, TC.Po < 100)
b <- subset(Pfrac.PME.resin.TC.Po, TC.Po < 200)
c <- subset(Pfrac.PME.resin.TC.Po, TC.Po > 100 & TC.Po < 200)
d <- subset(Pfrac.PME.resin.TC.Po, TC.Po > 200 & TC.Po < 300)
e <- subset(Pfrac.PME.resin.TC.Po, TC.Po > 300)
f <- subset(Pfrac.PME.resin.TC.Po, TC.Po > 200)
g <- subset(Pfrac.PME.resin.TC.Po, TC.Po > 200 & TC.Po < 20000)
h <- subset(Pfrac.PME.resin.TC.Po, TC.Po > 300 & TC.Po < 20000)
i <- subset(Pfrac.PME.resin.TC.Po, TC.Po < 20000)

cor.test(a$PME,a$TC.Po,method="spearman")
cor.test(b$PME,b$TC.Po,method="spearman")
cor.test(c$PME,c$TC.Po,method="spearman")
cor.test(d$PME,d$TC.Po,method="spearman")
cor.test(e$PME,e$TC.Po,method="spearman")
cor.test(f$PME,f$TC.Po,method="spearman")
cor.test(g$PME,g$TC.Po,method="spearman")
cor.test(h$PME,h$TC.Po,method="spearman")
cor.test(i$PME,i$TC.Po,method="spearman")

#descriptive statistics with extreme values removed
psych::describe(i)

###Bivariate correlations between variables and extracted factor scores
##correlation with factor scores (from the 1st section: exploratory factor analysis of resin subset)
FAresin_scores<-as.matrix(factor.scores(Pfrac.PME.resin.2, FAresin, method = "regression")$scores) #extract factor scores

##combine factor scores with initial dataset
Pfrac.PME.resin.descrip.fascores <- cbind(FAresin_scores,Pfrac.PME.resin.descrip)

###Correlation of factor scores with other soil properties
##pH subset
Pfrac.PME.resin.fascores.pH <- Pfrac.PME.resin.descrip.fascores[,c(1:3,23)] #subset pH, Pfractions, PME
Pfrac.PME.resin.fascores.pH.omit <- na.omit(Pfrac.PME.resin.fascores.pH) #omit missing observations
rcorr(as.matrix(Pfrac.PME.resin.fascores.pH.omit), type = "spearman") #Spearman correlation

##TC subset
Pfrac.PME.resin.fascores.TC <- Pfrac.PME.resin.descrip.fascores[,c(1:3,24)] #subset TC, Pfractions, PME
Pfrac.PME.resin.fascores.TC.omit <- na.omit(Pfrac.PME.resin.fascores.TC) #omit missing observations
rcorr(as.matrix(Pfrac.PME.resin.fascores.TC.omit), type = "spearman") #Spearman correlation

##Figures (Global Map Climate Distribution)
###Plot of Mean annual temperature and precipitation distribution (resin subset)
##subset and organize data for plot making
Pfrac.PME.resin.clim <- Pfrac.PME.resin.descrip[,c(3,4)]
#add frequency of each MAT + MAP to the climate data
Pfrac.PME.resin.clim.freq <- as.data.frame(table(Pfrac.PME.resin.clim))
Pfrac.PME.resin.clim.freq2 <- subset(Pfrac.PME.resin.clim.freq, Freq > 0)
#convert data to numeric
Pfrac.PME.resin.clim.freq2$MAT <- as.numeric(as.character(Pfrac.PME.resin.clim.freq2$MAT))
Pfrac.PME.resin.clim.freq2$MAP <- as.numeric(as.character(Pfrac.PME.resin.clim.freq2$MAP))
str(Pfrac.PME.resin.clim.freq2)

#plot mean annual temperature and precipitation 
p2<-ggplot(data = Pfrac.PME.resin.clim.freq2, aes(x=MAT, y=MAP)) +
  geom_point(aes(size = Freq, fill = Freq, alpha = 0.5), shape = 21, stroke = FALSE) +
  scale_fill_viridis(name = "Number of observations", option = "plasma", limits = c(1,80), breaks = seq(0,80,by=20)) +
  scale_size_continuous(name = "Number of observations", range = c(1, 25), limits = c(1,80), breaks = seq(0,80,by=20)) +
  theme_light() +
  guides(fill = guide_legend(), size = guide_legend()) +
  labs(x = "Mean Annual Temperature (\u00B0C)", y = "Mean Annual Precipitation (mm)", size = "Number of observations")
p2

###Creating map to show global distribution of observations (resin subset)
#map setup
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

##subset organize data for plot making
Pfrac.PME.resin.coord <- Pfrac.PME.resin.descrip[,c(1,2)]
#add frequency of each Lat+Long to the climate data
Pfrac.PME.resin.coord.freq <- as.data.frame(table(Pfrac.PME.resin.coord))
Pfrac.PME.resin.coord.freq2 <- subset(Pfrac.PME.resin.coord.freq, Freq > 0)
#convert data to numeric
Pfrac.PME.resin.coord.freq2$Longitude <- as.numeric(as.character(Pfrac.PME.resin.coord.freq2$Longitude))
Pfrac.PME.resin.coord.freq2$Latitude <- as.numeric(as.character(Pfrac.PME.resin.coord.freq2$Latitude))
str(Pfrac.PME.resin.coord.freq2)

#map drawing
map2 <- ggplot(data=world) + 
  geom_sf(fill = "white") +
  theme_bw() +
  geom_point(shape = 21, data = Pfrac.PME.resin.coord.freq2, aes(x=Longitude, y=Latitude, size = Freq, fill = Freq), alpha = 0.4) +
  scale_size_continuous(name = "Number of observations", range = c(1, 10), limits = c(1,80), breaks = seq(0,80,by=20)) +
  scale_fill_viridis(name = "Number of observations", option = "plasma", limits = c(1,80), breaks = seq(0,80,by=20)) +
  guides(fill = guide_legend(), size = guide_legend())

map2

##Code for Extracting climate data from WorldClim (remove # to run)
##data input
#coords <- read.csv("coordinates.csv", header = T)
#r <- getData("worldclim", var="bio", res=10)
#r <- r[[c(1,12)]]
#names(r) <- c("Temp","Prec")

#points <- SpatialPoints(coords, proj4string = r@crs)
#values <- extract(r,points)

#map_mat_data <- cbind.data.frame(coordinates(points),values)
#map_mat_data