# TABLE OF CONTENTS
# 0. Prepare environment
# 1. Model Specification and inference
# 2. Analyze Results
                 


# 0. Prepare environment ----
rm(list=ls())
list.of.packages <- c("mgcv", "data.table","itsadug","parallel","plotfunctions","xtable")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

#set seed
set.seed(17)

#load packages
require(mgcv)
require(data.table)
require(itsadug)
require(plotfunctions)
require(xtable)

#sessioninfo
print("-------------------Session Info-----------------")
print(sessionInfo())

print("-------------------Data pre-processing-----------------")

print("Reading data")
X.e <- fread("/data/p275030/elektraGAMdata.csv", sep = ";", dec = ",")
X.e <- X.e[Datum %between% c("2015-01-01", "2016-12-31")]
X.e <- droplevels(X.e)
mylevels <- unique(X.e$EAN)
X.e$EAN <- factor(X.e$EAN, levels=mylevels)
X.e$weekend <- as.factor(X.e$weekend)


print("number of unique households")
print(length(levels(X.e$EAN)))

print("subsetting data")
X.pv <- subset(X.e, PV==TRUE)
X.nopv <- subset(X.e, PV==FALSE)
X.pv <- droplevels(X.pv)
X.nopv <- droplevels(X.nopv)


print("Number of households with PV")
print(length(levels(X.pv$EAN)))

print("Number of households without PV")
print(length(levels(X.nopv$EAN)))

#drop unused factors
#droplevels(X.pv$EAN)
#droplevels(X.nopv$EAN)

print("sampling EAN")
#pv owners
set.pv <- sample(levels(X.pv$EAN), 350)
X.pv.600 <- X.pv[X.pv$EAN %in% set.pv, ]
X.pv.600 <- droplevels(X.pv.600)
#non-pv
set.nopv <- sample(levels(X.nopv$EAN), 350)
X.nopv.600 <- X.nopv[X.nopv$EAN %in% set.nopv, ]
X.nopv.600 <- droplevels(X.nopv.600)


print("Final dataset PV households")
print(length(levels(X.pv.600$EAN)))

print("Final dataset non-PV households")
print(length(levels(X.nopv.600$EAN)))

print("joining files")
X.sample <- rbind(X.pv.600,X.nopv.600)
X.e <- droplevels(X.sample)
rm(list=setdiff(ls(), "X.e"))

#set contrasts
# change factor to ordered factor:
X.e$PV <- as.ordered(X.e$PV)
# change contrast to treatment coding (difference curves)
contrasts(X.e$PV) <- 'contr.treatment'

# 0. Model Specification, estimation, and inference ----
print("-------------------MODELS-----------------")

# m.0.0 : base model with only main effects ------------------------------------
print("Running mod m.0.0")
system.time(m.0.0 <- bam(
  net ~   PV + load24 +
    # hour of day effects
    s(hour_of_day, bs = "cc", k=23) +
    s(hour_of_day, by = PV, bs = "cc", k=23) +
    #month of year effects
    s(month, bs = "cc", k=11) +
    s(month, by=PV, bs = "cc", k=11),
  data = X.e,
  drop.unused.levels=FALSE,
  na.action="na.omit",
  discrete=TRUE,
  nthreads=8
))


# m.0.1 : base model with main effects and interactions ------------------------
print("Running mod m.0.1")
system.time(m.0.1 <- bam(
  net ~   PV + load24 +
    # hour of day effects
    s(hour_of_day, bs = "cc", k=23) +
    s(hour_of_day, by = PV, bs = "cc", k=23) +
    #month of year effects
    s(month, bs = "cc", k=11) +
    s(month, by=PV, bs = "cc", k=11)+
    #reference surface
    ti(hour_of_day, month, k = c(23,11)) +
    #difference surface
    ti(hour_of_day, month, by = PV, k = c(23,11)),
  data = X.e,
  drop.unused.levels=FALSE,
  na.action="na.omit",
  discrete=TRUE,
  nthreads=8
))

# m.0.2 : base model with main effects, interactions, and random smooths-------
print("Running mod m.0.2")
system.time(m.0.2 <- bam(
  net ~   PV + load24 +
    # hour of day effects
    s(hour_of_day, bs = "cc", k=23) +
    s(hour_of_day, by = PV, bs = "cc", k=23) +
    #month of year effects
    s(month, bs = "cc", k=11) +
    s(month, by=PV, bs = "cc", k=11)+
    #reference surface
    ti(hour_of_day, month, k = c(23,11)) +
    #difference surface
    ti(hour_of_day, month, by = PV, k = c(23,11))+
    #random smooth for each house
    s(hour_of_day,  EAN, bs = 'fs', k = 23,  m = 1),
  data = X.e,
  drop.unused.levels=FALSE,
  na.action="na.omit",
  discrete=TRUE,
  nthreads=22
))



# m.1 final model: main effects, interactions, random effects, and res --------
#include AR term
X.e <- start_event(X.e, column = "time", event = c("EAN"))
#autocorrelation parameter
valRho <- start_value_rho(m.0.2, plot = TRUE)

print("Running mod m.1")
system.time(m.1 <- bam(
  net ~   PV + load24 +
    # hour of day effects
    s(hour_of_day, bs = "cc", k=23) +
    s(hour_of_day, by = PV, bs = "cc", k=23) +
    #month of year effects
    s(month, bs = "cc", k=11) +
    s(month, by=PV, bs = "cc", k=11)+
    #reference surface
    ti(hour_of_day, month, k = c(23,11)) +
    #difference surface
    ti(hour_of_day, month, by = PV, k = c(23,11))+
    #random smooth for each house
    s(hour_of_day,  EAN, bs = 'fs', k = 23,  m = 1),
  data = X.e,
  drop.unused.levels=FALSE,
  na.action="na.omit",
  discrete=TRUE,
  AR.start = X.e$start.event,
  rho = valRho,
  nthreads=22 
))

print("saving models")
saveRDS(m.0.0, '/data/p275030/GAM/models/final/m00.rds')
saveRDS(m.0.1, '/data/p275030/GAM/models/final/m01.rds')
saveRDS(m.0.2, '/data/p275030/GAM/models/final/m02.rds')
saveRDS(m.1, '/data/p275030/GAM/models/final/m1.rds')

# 2. Model Comparison ----
print("-------------------Model Comparison-----------------")

# Model comparision -------------------------------------------------------
mcompAIC <- AIC(m.0.0, m.0.1, m.0.2)
print("Model comparison (base models with no residual autocorrelation)")
print(xtable(mcompAIC))

print("Model comparision 2 (with residual autocorrelation)")
mcompML <- compareML(m.0.2, m.1, suggest.report = TRUE)
print(xtable(mcompML$table))

rm(m.0.0, m.0.1)


#Generate plots
# 2. Results: Plots ----
print("-------------------Results: Graphs-----------------")


#residual plots
pdf("results/images/M02acf.pdf",width=6,height=4,paper='special') 
acf_resid(m.0.2, split_pred = list(X.e$EAN), main="Model m0 residual plot: Averaged ACF") #averaged residuals 
dev.off()

#Model with AR term
pdf("results/images/M1acf.pdf",width=6,height=4,paper='special') 
acf_resid(m.1, split_pred = list(X.e$EAN), main="Model m1 residual plot: Averaged ACF") #averaged residuals 
dev.off()


#Plot smooths


#smooth per month
for (d in c(1:12)) {
  pdf(paste("results/images/M1smooth",d,".pdf", sep="_"),width=6,height=4,paper='special') 
  par(mfrow = c(1, 1), cex = 0.7)
  #visualise smooths
  plot_smooth(
    m.1,
    view = "hour_of_day",
    plot_all = "PV",
    cond = list(month=d),
    rm.ranef = TRUE,
    xlab = "hour of day",
    ylab = "energy consumption (kWh)",
    hide.label = TRUE,
    print.summary=FALSE,
    col = c("#E64B35FF","#4DBBD5FF"),
    legend_plot_all = list(x=45,y=25)
  )
  # Add legend:
  legend('bottomright', 
         legend=c("Prosumers","Consumers"),
         lwd=2,
         col = c("#E64B35FF","#4DBBD5FF"),
         bty='n')
  
  dev.off()
}



#plot diff smooths
diffDF <- list()
for (d in c(1:12)) {
  pdf(paste("results/images/M1dailySmoothdiff",d,".pdf", sep="_"),width=6,height=4,paper='special') 
  par(mfrow = c(1, 1), cex = 0.7)
            out <-   plot_diff(m.1, view = c("hour_of_day"),
            comp = list(PV = c(TRUE, FALSE)),
            xlab = "hour of day",
            cond = list(month=d),
            ylab = "energy consumption (kWh)",
            main="",
            rm.ranef = TRUE,
            hide.label = TRUE,
            plot=TRUE,
            print.summary = FALSE
  )
  dev.off()
  x <- find_difference(out$est, out$CI, f=1, xVals=out$hour_of_day)
  diffDF[[d]] <- x
}

print("saving difference smooth values")
saveRDS(diffDF, '/data/p275030/GAM/models/final/diffDF.rds')

#plot interactions
#difference surfaces
pdf("results/images/M1dailyMonthlySmooths.pdf",width=6,height=4,paper='special') 
par(mfrow=c(1,2), cex = 0.5)
# Note: specify zlim when comparing two plots
fvisgam(m.1, view = c("hour_of_day","month"), 
        cond=list(PV = TRUE),
        plot.type='contour', 
        color='topo', 
        main='Prosumers',
        xlab = "hour of day",
        ylab = "month",
        zlim = c(-0.35,0.3))
fvisgam(m.1, 
        view = c("hour_of_day","month"),
        cond=list(PV = FALSE),
        plot.type='contour',
        color='topo', 
        main='Consumers',
        xlab = "hour of day",
        ylab = "month",
        zlim = c(-0.35,0.3))
dev.off()

print("-------------------Results: Model Check-----------------")

png("results/images/m1ModelCheck.png") 
par(mfrow = c(2,2))
m.1.check <- gam.check(m.1,  type="deviance")
dev.off()

print("saving model check")
saveRDS(m.1.check, '/data/p275030/GAM/models/final/m1ModelCheck.rds')

#Extract tables
print("-------------------Results: Model Tables-----------------")
m.1.tab <- gamtabs(m.1)
saveRDS(m.1.tab, 'results/tables/final/m1tab.rds')


print("I'm glad this worked out")
