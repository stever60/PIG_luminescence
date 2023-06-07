setwd("/Users/Steve/Dropbox/BAS/Data/R/Luminescence")
library(Luminescence)
.pardefault <- par()
par(.pardefault)
#clear previous console if needed
remove (list = ls())

##PIG7 (EXP OR LIN) adn OT (EXP+LIN, EXP OR LIN)
#clear plot window if needed
dev.off()
#resets parameters to default - useful for rescaling if it goes off
.pardefault <- par()
par(.pardefault)

#plot_GrowthCurve - Fits and plot a growth curve for luminescence data (Lx/Tx against dose, with Test dose response Tx/Tn in third column)
sample <- read.csv ("inputGC_OT.csv")
plot_GrowthCurve(sample, mode = "interpolation", fit.method = "EXP+LIN", 
                 fit.force_through_origin = TRUE,fit.bounds = TRUE, 
                 fit.weights = TRUE, NumberIterations.MC = 1000)



#No fit to error weighting - plot_GrowthCurve - Fits and plot a growth curve for luminescence data (Lx/Tx against dose, with Test dose response Tx/Tn in third column)
sample <- read.csv ("inputGC_OT.csv")
plot_GrowthCurve(sample, mode = "interpolation", fit.method = "EXP", 
                 fit.force_through_origin = TRUE,fit.bounds = TRUE, 
                 fit.weights = FALSE, NumberIterations.MC = 1000)

#sample is a data frame with three columns for x=Dose,y=LxTx,z=LxTx.Error,
#y1=TnTx. The column for the test dose response is optional, but requires 'TnTx'
#as column name if used. For exponential fits at least three dose points including the natural should be provided.
#A dose response curve is produced for luminescence measurements using a regenerative or additive
#protocol. The function supports interpolation and extrapolation to calculate the equivalent dose.

#print to PDF with default settings
pdf()

dev.off()
#other options   
#na.rm = TRUE, fit.includingRepeatedRegPoints = FALSE, fit.NumberRegPoints = 8,fit.NumberRegPointsReal = 8, fit.force_through_origin = TRUE, 
#plot_GrowthCurve, fit.bounds = TRUE,
#NumberIterations.MC = 100, output.plot = TRUE, output.plotExtended = TRUE, fit.NumberRegPoints = NULL, output.plotExtended.single = FALSE,
#cex.global = 1, txtProgressBar = TRUE, verbose = TRUE)

#clear plot window if needed
dev.off()
#resets parameters to default - useful for rescaling if it goes off
.pardefault <- par()
par(.pardefault)



## PIG 7 DATA INPUT - GRAPH PLOTTING

#Temp Range cps regression De data for 250-340, 280-350, 290-330 C intervals 
data.DeValuesTR <- read.csv ("inputKDE_TempRanges.csv")

#data.1-data.3 - 24 lines for 3 TempRanges x 8 discs
data.TR1 <- data.DeValuesTR[1:8,]
data.TR2 <- data.DeValuesTR[9:16,]
data.TR3 <- data.DeValuesTR[17:24,]

#discs 1-8 - PIG7 Temp Ranges cps 250-340C cps integrals
data.TempRanges.1 <- data.DeValuesTR[1:8,]
#discs 1-8 - PIG7 Temp Ranges cps 280-350C cps integrals 
data.TempRanges.2 <- data.DeValuesTR[9:16,]
#discs 1-8 - Temp Ranges cps 290-330C cps integrals 
data.TempRanges.3 <- data.DeValuesTR[17:24,]

#merge data again into various combinations - different regression analysis, different plots - exploratory only
#discs 1-8 - Temp Ranges 250-340, 280-340, 300-340C cps integrals undertaken before 
#regression plotted separately as different colour on the same plot
data.TempRanges.4 <- list(data.TR1, data.TR2, data.TR3)
#discs 1-8 - Temp Ranges 250-340C, 280-340C cps integrals
data.TempRanges.5 <- list(data.TR1, data.TR2)



#Temp Range cps regression 8DM De data for 250-340, 280-350, 290-330 C intervals 
data.DeValuesTR8DM <- read.csv ("inputKDE_8DM_10C.csv")

#data.1-data.3 - 24 lines for 3 TempRanges x 8 discs
data.TR1 <- data.DeValuesTR8DM[1:9,]
data.TR2 <- data.DeValuesTR8DM[4:10,]
data.TR3 <- data.DeValuesTR8DM[5:8,]



## Input 10C bands regression De data then combined to 250-340, 280-350, 290-330 C temp ranges 
data.DeValues <- read.csv ("inputKDE_10C.csv") #import all the data including extrapolated De values
## OR import interpolated De values only ##
data.DeValues_inter <- read.csv ("inputKDE_10C_inter.csv") #see below for ranges

#data.1-data.9 - 72 lines for 6, 10C integrals (250-340, 280-350, 290-330 C) x 8 discs
data.1 <- data.DeValues[1:8,] #250-260 C
data.2 <- data.DeValues[9:16,] #260-270 C
data.3 <- data.DeValues[17:24,] #270-280 C
data.4 <- data.DeValues[25:32,] #280-290 C
data.5 <- data.DeValues[33:40,] #290-300 C
data.6 <- data.DeValues[41:48,] #300-310 C
data.7 <- data.DeValues[49:56,] #310-320 C
data.8 <- data.DeValues[57:64,] #320-330 C
data.9 <- data.DeValues[65:72,] #330-340 C
data.10 <- data.DeValues[73:80,] #340-350 C


#10C interval regression 
data.10C.1 <- data.DeValues[1:72,] #250-340C all data - single plot
data.10C.2 <- data.DeValues[25:80,] #280-350C all data - single plot
data.10C.3 <- data.DeValues[33:64,] #290-330C all data - single plot

#10C interval regression - discs 1-8 - 250-340C as separate 10C integral x 8 disc plots
data.10C.4 <- list(data.1, data.2, data.3, data.4, data.5, data.6, data.7, data.8, data.9)
#10C interval regression - discs 1-8 - 280-340C as separate 10C integral x 8 disc plots
data.10C.5 <- list(data.4, data.5, data.6,data.7, data.8, data.9)
#10C interval regression - discs 1-8 - 300-340C as separate 10C integral x 8 disc plots
data.10C.6 <- list(data.6, data.7, data.8, data.9)


## OR import interpolated 10C band integral data De values only ##
data.DeValues_inter <- read.csv ("inputKDE_10C_inter.csv")

#data.1-data.9 - 72 lines for 6, 10C integrals (250-340, 280-350, 290-330 C) x 8 discs
data.1 <- data.DeValues_inter[1:4,] #250-260 C
data.2 <- data.DeValues_inter[5:9,] #260-270 C
data.3 <- data.DeValues_inter[10:15,] #270-280 C
data.4 <- data.DeValues_inter[16:22,] #280-290 C
data.5 <- data.DeValues_inter[23:29,] #290-300 C
data.6 <- data.DeValues_inter[30:36,] #300-310 C
data.7 <- data.DeValues_inter[37:44,] #310-320 C
data.8 <- data.DeValues_inter[45:51,] #320-330 C
data.9 <- data.DeValues_inter[52:59,] #330-340 C
data.10 <- data.DeValues_inter[60:65,] #340-350 C

data.1

#10C interval regression 
data.10C.0 <- data.DeValues_inter[1:29,] #250-300C all data - single plot
data.10C.1 <- data.DeValues_inter[1:59,] #250-340C all data - single plot
data.10C.2 <- data.DeValues_inter[16:65,] #280-350C all data - single plot
data.10C.3 <- data.DeValues_inter[23:51,] #290-330C all data - single plot

data.10C.1

#10C interval regression - discs 1-8 - 250-340C as separate 10C integral x 8 disc plots
data.10C.4 <- list(data.1, data.2, data.3, data.4, data.5, data.6, data.7, data.8, data.9)
#10C interval regression - discs 1-8 - 280-350C as separate 10C integral x 8 disc plots
data.10C.5 <- list(data.4, data.5, data.6,data.7, data.8, data.9, data.10)
#10C interval regression - discs 1-8 - 290-330C as separate 10C integral x 8 disc plots
data.10C.6 <- list(data.5, data.6, data.7, data.8)

data.10C.7 <- list(data.10C.4, data.10C.5, data.10C.6)


#Cosmo ages - Fildes peninsula for KDE
data.cosmo <- read.csv ("inputKDE_cosmo.csv")

#data.1-data.3 - 24 lines for 3 TempRanges x 8 discs
data.R1 <- data.cosmo[1:3,]
data.R2 <- data.cosmo[4:6,]
data.R3 <- data.cosmo[7:9,]
data.R4 <- data.cosmo[1:9,]


## DATA PROCESSING - GRAPH PLOTTING

## Using Kernel Density Estimate and Abanico Plots
## Use for comprehensive presentation of data precision and its dispersion
# around a central value as well as illustration of a kernel density estimate, 
# histogram and/or dot plot of the dose values

#Kernel Density Estimate (KDE)plot 
#clear plot window if needed
dev.off()
plot_KDE(data.R4,
         range = "290-330 C", #remember to change this 
         stats = c("mean", "sd", "se"),
         summary = c("n", "mean", "sd.abs", "sd.rel", "se.abs", "median"), #"in.2s",
         summary.pos = c("sub"),
         pch = c(16, 16, 16, 15, 15, 15, 17, 17, 17), 
         col = c("black", "blue","red", "black", "blue","red", "black", "blue","red"))

#plot_AbanicoPlot Function to create a radial plot and weighted mean De value
dev.off()
plot_AbanicoPlot(data.R4, hist = FALSE, boxplot = TRUE, dispersion = c("qr"),
                 pch = c(16, 16, 16, 15, 15, 15, 17, 17, 17), 
                 col = c("black", "blue","red", "black", "blue","red", "black", "blue","red"),
                 z.0 = "mean.weighted", #error.bars = FALSE, #bar.col = c("grey"),
                 summary.method = "MCM", #"MCM" #"weighted" #"MCM"
                 summary = c("n", "mean", "sd.abs", "sd.rel","median","in.2s"),
                 summary.pos = "topleft",
                 stats = c("min", "max", "median", "mean"),
                 legend.pos = "topleft")

## Apply a 3-parameter Minimum Dose Model MAM and plot as summary Abanico plot
## Minimum age model applies the (un-)logged minimum age model (MAM) 
## after Galbraith et al. (1999) to a given De distribution

mam <- calc_MinDose(data = data.R4,
                    sigmab = 0.2,
                    plot = FALSE,
                    bootstrap = TRUE)
# Show summary table that contains the most relevant results
resmin <- get_RLum(mam, "summary")
resmin

# Plot the log likelihood profiles retroactively, because before we set plot = FALSE
# plot_RLum(mam) 


# Plot the dose distribution in an abanico plot and draw a line at the minimum dose estimate
#clear plot window
dev.off()
title = "Fildes ALL"
plot_AbanicoPlot(data = data.R4,
                 main = paste ("3-par. MAM: ", title),
                 line = mam,polygon.col = "none",
                 hist = FALSE,
                 boxplot = TRUE,
                 rug = TRUE,
                 dispersion = c("qr"),
                 pch = c(16, 16, 16, 15, 15, 15, 17, 17, 17), 
                 col = c("black", "blue","red", "darkblue", "cyan4","brown", 
                         "darkslategray", "deepskyblue","darkred"),
                 z.0 = "mean.weighted", #error.bars = FALSE, #bar.col = c("grey"),
                 summary.method = "MCM", #"MCM" #"weighted" #"MCM"
                 summary = c("n", "mean", "mean.weighted", "sd.abs", "sd.rel", "se.abs", "median", "in.2s"),
                 centrality = resmin$de,
                 stats = c("min", "max"),
                 line.col = "red",
                 grid.col = "none",
                 line.label = paste0(round(resmin$de, 0), "\U00B1",
                                     round(resmin$de_err, 0), " Gy"),
                 bw = 0.1,
                 ylim = c(-25, 18),
                 summary.pos = "topleft",
                 mtext = bquote("Parameters: " ~
                                  sigma[b] == .(get_RLum(mam, "args")$sigmab) ~ ", " ~
                                  gamma == .(round(log(resmin$de), 1)) ~ ", " ~
                                  sigma == .(round(resmin$sig, 1)) ~ ", " ~
                                  rho == .(round(resmin$p0, 2))))


##ADDITIONAL AGE MODELS 
## MAXAM
## Applies Function to fit the maximum age model to De data. This is a wrapper function that calls calc_MinDose()
# and applies a similiar approach as described in Olley et al. (2006). Data transformation
#To estimate the maximum dose population and its standard error, the three parameter minimum age
#model of Galbraith et al. (1999) is adapted. The measured De values are transformed as follows:
# 1. convert De values to natural logs
# 2. multiply the logged data to creat a mirror image of the De distribution
# 3. shift De values along x-axis by the smallest x-value found to obtain only positive values
# 4. combine in quadrature the measurement error associated with each De value with a relative error specified by sigmab
# 5. apply the MAM to these data
# When all calculations are done the results are then converted as follows
# 1. subtract the x-offset
# 2. multiply the natural logs by -1
# 3. take the exponent to obtain the maximum dose estimate in Gy

# clear plot window
dev.off()
# Run the MAM model,  save results to a variable and turn
# plotting of the log-likelihood profiles off.

maxam <- calc_MaxDose(data = data.TempRanges.1,
                      sigmab = 20, par = 3,
                      plot = FALSE,
                      bootstrap = FALSE)
# Show summary table that contains the most relevant results
resmax <- get_RLum(maxam, "summary")
resmax

# Use calc_FiniteMixture  to apply a finite mixture model (FMM) after Galbraith (2005) to a
# given De distribution as a check - the FMM is envisioned only for single grain data 
# (Galbraith and Green, 1990).

#This function fits a k-component mixture to a De distribution with differing known standard errors.
#Parameters (doses and mixing proportions) are estimated by maximum likelihood assuming that the
#log dose estimates are from a mixture of normal distributions.

#clear plot window
dev.off()
## show structure of the results
FMM
## (1) apply the finite mixture model
## avery small sigmab is used for in case the finite mixture model is not suitable
calc_FiniteMixture(data.TempRanges.1,
                   sigma  = 0.2, n.components = 2,
                   grain.probability = TRUE)
## (2) repeat the finite mixture model for 2, 3 and 4 maximum number of fitted
## components and save results
FMM<- calc_FiniteMixture(data.TempRanges.1,
                         sigmab = 0.2, n.components = c(2:4),
                         pdf.weight = TRUE, dose.scale = c(0, 300))
## show the results on equivalent dose, standard error and proportion of fitted components
res <- get_RLum(object = FMM, data.object = "components")
res
res1 <- get_RLum(FMM, "summary")
res1






#resets parameters to default - useful for rescaling if it goes off
.pardefault <- par()
par(.pardefault)

#plot_RadialPlot Function to create a Radial Plot

#Galbraith’s radial plot is produced on a logarithmic or a linear scale.

## split data set into sub-groups, can be for manipulation, and merge again
data.DeValues <- read.csv ("inputKDE.csv")

data.1 <- data.DeValues[1:8,]
data.2 <- data.DeValues[9:16,]
data.3 <- data.DeValues[17:24,]
#merge data again
data.4 <- list(data.1, data.2, data.3)

plot_RadialPlot(data.4, centrality = "mean.weighted", pch=c(19, 19, 19), col = c("blue", "black","red"), 
                summary = c("n", "in.2s", "mean.weighted", "se.weighted","median"), 
                summary.pos = "sub",
                stats = c("min", "max", "median"))

#with a brief statistical summary plot_RadialPlot(data = ExampleData.DeValues,
summary = c("n", "in.2s"))
#with another statistical summary as subheader plot_RadialPlot(data = ExampleData.DeValues,
summary = c("n", "in.2s", "mean.weighted", "median"), summary.pos = "sub")
#with central value defined 
central.value = 200,

#options: na.rm = TRUE, log.z = TRUE, central.value, centrality = "mean.weighted", mtext,
 # summary, summary.pos,plot_RadialPlot 249, legend, legend.pos, stats,
 # rug = FALSE, plot.ratio, bar.col, y.ticks = TRUE, grid.col,line,line.col, 
 # line.label, output = FALSE)






#resets parameters to default - useful for rescaling if it goes off
.pardefault <- par()
par(.pardefault)

## split data set into sub-groups, can be for manipulation, and merge again
data.DeValues <- read.csv ("inputKDE.csv")

data.1 <- data.DeValues[1:8,]
data.2 <- data.DeValues[9:16,]
data.3 <- data.DeValues[17:24,]
#merge data again
data.4 <- list(data.1, data.2, data.3)








## now with legend, colour, different points and smaller scale plot_RadialPlot(data = ExampleData.DeValues,
legend.text = "Sample 1", col = "tomato4",
plot_RadialPlot
253
bar.col = "peachpuff", pch = "R",
cex = 0.8)
## now without 2-sigma bar, y-axis, grid lines and central value line plot_RadialPlot(data = ExampleData.DeValues,
bar.col = "none", grid.col = "none", y.ticks = FALSE, lwd = 0)
## now with user-defined axes labels plot_RadialPlot(data = ExampleData.DeValues,
xlab = c("Data error (%)",
         "Data precision"),
ylab = "Scatter",
zlab = "Equivalent dose [Gy]")
## now with minimum, maximum and median value indicated plot_RadialPlot(data = ExampleData.DeValues,
central.value = 150,
stats = c("min", "max", "median"))
## now with a brief statistical summary plot_RadialPlot(data = ExampleData.DeValues,
summary = c("n", "in.2s"))
## now with another statistical summary as subheader plot_RadialPlot(data = ExampleData.DeValues,
summary = c("mean.weighted", "median"), summary.pos = "sub")
## now the data set is split into sub-groups, one is manipulated data.1 <- ExampleData.DeValues[1:15,]
data.2 <- ExampleData.DeValues[16:25,] * 1.3
## now a common dataset is created from the two subgroups data.3 <- list(data.1, data.2)
## now the two data sets are plotted in one plot plot_RadialPlot(data = data.3)
## now with some graphical modification plot_RadialPlot(data = data.3,
col = c("darkblue", "darkgreen"), bar.col = c("lightblue", "lightgreen"), pch = c(2, 6),
summary = c("n", "in.2s"),
summary.pos = "sub",
legend = c("Sample 1", "Sample 2"))