# ------------------------------------------------------
# Program: TimeSeries2.R 
#  Author: Steven M. Boker
#    Date: Sun May 14 11:22:12 EDT 2017
#
#  Timeseries, State Space and Vector Plots Part 2
#
# ------------------------------------------------------

# ------------------------------------------------
# Load required libraries. 

library(lattice)
library(stats)

# ------------------------------------------------------
#  Create an error-free autoregressive time series and 
#  synchronous cross-regression relation
#  starting far from equilibrium

tX0 <- rep(NA, 201)
tX0[1] <- 100
for (i in 2:201) {
   tX0[i] <- .9 * tX0[i-1]
}

tY0 <- 5 + .5 * tX0 + rnorm(201, mean=0, sd=2)
tX0 <- tX0 + rnorm(201, mean=0, sd=2)

# ------------------------------------------------------
#  Create an error-free autoregressive time series and 
#  synchronous cross-regression relation 
#  starting at equilibrium

tX1 <- rep(NA, 201)
tX1[1] <- 1
for (i in 2:201) {
   tX1[i] <- .9 * tX1[i-1]
}

tY1 <- 5 + .5 * tX1 + rnorm(201, mean=0, sd=2)
tX1 <- tX1 + rnorm(201, mean=0, sd=2)

# ------------------------------------------------------
#  Create a timeseries of tX1
 
tTime <- seq(0,1000, by=5)/60
tData1a <- data.frame(x=tX1, y=tY1, time=tTime)

pdf("TimeSeries2TSPlot1aX.pdf", height=5, width=6)
plot(c(0,15), c(-15,15), 
     xlab="Time",
     ylab="X",
     type='n',
     pch=1)
lines(tData1a$time, tData1a$x, type='l')
dev.off()

pdf("TimeSeries2TSPlot1aY.pdf", height=5, width=6)
plot(c(0,15), c(-15,15), 
     xlab="Time",
     ylab="Y",
     type='n',
     pch=1)
lines(tData1a$time, tData1a$y, type='l')
dev.off()


#####################
#   MARKER A
#####################

# ------------------------------------------------------
#  Create an error dependent autoregressive time series and 
#  synchronous cross-regression relation

tX2 <- rep(NA, 201)
tX2[1] <- 1
for (i in 2:201) {
   tX2[i] <- .9 * tX2[i-1] + rnorm(1, mean=0, sd=2)
}

tY2 <- 5 + .5 * tX2 + rnorm(201, mean=0, sd=2)

tTime <- seq(0,1000, by=5)/60
tData2 <- data.frame(x=tX2, y=tY2, time=tTime)
summary(tData2)

pdf("TimeSeries2TSPlot2X.pdf", height=5, width=6)
plot(c(0,15), c(-15,15), 
     xlab="Time",
     ylab="X",
     type='n',
     pch=1)
lines(tData2$time, tData2$x, type='l')
dev.off()

pdf("TimeSeries2TSPlot2Y.pdf", height=5, width=6)
plot(c(0,15), c(-15,15), 
     xlab="Time",
     ylab="Y",
     type='n',
     pch=1)
lines(tData2$time, tData2$y, type='l')
dev.off()


#####################
#   MARKER B
#####################

# ------------------------------------------------------
#  Create a scrambled timeseries of tX1
 
tTimeRand1 <- order(runif(201,0,1))
tData1b <- data.frame(x=tX1[tTimeRand1], time=tTime)


pdf("TimeSeries2TSPlot1Xb.pdf", height=5, width=6)
plot(c(0,15), c(-15,15), 
     xlab="Time",
     ylab="X",
     type='n',
     pch=1)
lines(tData1b$time, tData1b$x, type='l')
dev.off()


#####################
#   MARKER C
#####################

# ------------------------------------------------------
#  Create a scrambled timeseries of tX2
 
tTimeRand2 <- order(runif(201,0,1))
tData2b <- data.frame(x=tX2[tTimeRand2], time=tTime)


pdf("TimeSeries2TSPlot2Xb.pdf", height=5, width=6)
plot(c(0,15), c(-15,15), 
     xlab="Time",
     ylab="X",
     type='n',
     pch=1)
lines(tData2b$time, tData2b$x, type='l')
dev.off()

#####################
#   MARKER D
#####################

# ------------------------------------------------------
#  Create a time-lagged timeseries of tX2
 
tau <- 3
tLen <- length(tX2)
tDataLag2 <- data.frame(x1=tX2[1:(tLen-tau)], x2=tX2[(1+tau):tLen])


pdf("TimeSeries2SSPlot1.pdf", height=5, width=5)
plot(c(-15,15), c(-15,15), 
     xlab="X(t)",
     ylab="X(t+3)",
     type='n',
     pch=1)
lines(tDataLag2$x1, tDataLag2$x2, type='l')
dev.off()

tau <- 15
tLen <- length(tX2)
tDataLag2 <- data.frame(x1=tX2[1:(tLen-tau)], x2=tX2[(1+tau):tLen])

pdf("TimeSeries2SSPlot1-10.pdf", height=5, width=5)
plot(c(-15,15), c(-15,15), 
     xlab="X(t)",
     ylab="X(t+15)",
     type='n',
     pch=1)
lines(tDataLag2$x1, tDataLag2$x2, type='l')
dev.off()

tau <- 3
tLen <- length(tX2)
tDataLag2Scram <- data.frame(x1=tData2b$x[1:(tLen-tau)], 
                        x2=tData2b$x[(1+tau):tLen])


pdf("TimeSeries2SSPlot2.pdf", height=5, width=5)
plot(c(-15,15), c(-15,15), 
     xlab="X(t)",
     ylab="X(t+3)",
     type='n',
     pch=1)
lines(tDataLag2Scram$x1, tDataLag2Scram$x2, type='l')
dev.off()


# ------------------------------------------------------
#  Create an autocorrelation function plot of tX0
#  the timeseries far from equilibrium
 
pdf("TimeSeries2TSPlotX0.pdf", height=5, width=6)
plot(c(0,15), c(-10,100), 
     xlab="Time",
     ylab="X",
     type='n',
     pch=1)
lines(tTime, tX0, type='l')
dev.off()

#####################
#   MARKER E
#####################

maxTau <- 100
tLen <- length(tX0)
tACF1 <- rep(NA, maxTau+1)

for (tau in (0:maxTau)) {
   t1 <- tX0[1:(tLen-tau)]
   t2 <- tX0[(1+tau):tLen]
   tSelect <- !is.na(t1) & !is.na(t2)
   if (length(t1[tSelect]) < 5) next
   tACF1[tau+1] <- cor(t1[tSelect], t2[tSelect])
}

backward <- seq(maxTau+1, 2, by= -1)
pdf("TimeSeries2ACF0.pdf", height=5, width=6)
plot(c(-maxTau,maxTau), c(-1,1), type="n",
     xlab="Lag",
     ylab="Autocorrelation")
lines(c(-maxTau:maxTau), c(tACF1[backward],tACF1), type="l", lty=1)
lines(c(-maxTau,maxTau), c(0,0), type="l", lty=2)
dev.off()


# ------------------------------------------------------
#  Create an autocorrelation function plot of tX1
 
maxTau <- 100
tLen <- length(tX1)
tACF1 <- rep(NA, maxTau+1)

for (tau in (0:maxTau)) {
   t1 <- tX1[1:(tLen-tau)]
   t2 <- tX1[(1+tau):tLen]
   tSelect <- !is.na(t1) & !is.na(t2)
   if (length(t1[tSelect]) < 5) next
   tACF1[tau+1] <- cor(t1[tSelect], t2[tSelect])
}

backward <- seq(maxTau+1, 2, by= -1)
pdf("TimeSeries2ACF1.pdf", height=5, width=6)
plot(c(-maxTau,maxTau), c(-1,1), type="n",
     xlab="Lag",
     ylab="Autocorrelation")
lines(c(-maxTau:maxTau), c(tACF1[backward],tACF1), type="l", lty=1)
lines(c(-maxTau,maxTau), c(0,0), type="l", lty=2)
dev.off()

# ------------------------------------------------------
#  Create an autocorrelation function plot of tX2
 
maxTau <- 100
tLen <- length(tX2)
tACF2 <- rep(NA, maxTau+1)

for (tau in (0:maxTau)) {
   t1 <- tX2[1:(tLen-tau)]
   t2 <- tX2[(1+tau):tLen]
   tSelect <- !is.na(t1) & !is.na(t2)
   if (length(t1[tSelect]) < 5) next
   tACF2[tau+1] <- cor(t1[tSelect], t2[tSelect])
}

backward <- seq(maxTau+1, 2, by= -1)
pdf("TimeSeries2ACF2.pdf", height=5, width=6)
plot(c(-maxTau,maxTau), c(-1,1), type="n",
     xlab="Lag",
     ylab="Autocorrelation")
lines(c(-maxTau:maxTau), c(tACF2[backward],tACF2), type="l", lty=1)
lines(c(-maxTau,maxTau), c(0,0), type="l", lty=2)
dev.off()

#####################
#   MARKER F
#####################

# ------------------------------------------------------
#  Using the acf() function from the library(stats)
 
pdf("TimeSeries2ACF0b.pdf", height=5, width=6)
tACF = acf(tX0, 
           lag.max=maxTau, 
           type="correlation", 
           plot=FALSE, 
           na.action=na.fail)
plot(tACF, 
     main="",
     ci=.95, 
     ci.col="blue", 
     ci.type="white",
     type="l", 
     cex=1.25)
dev.off()

pdf("TimeSeries2ACF1b.pdf", height=5, width=6)
tACF = acf(tX1, 
           lag.max=maxTau, 
           type="correlation", 
           plot=FALSE, 
           na.action=na.fail)
plot(tACF, 
     main="",
     ci=.95, 
     ci.col="blue", 
     ci.type="white",
     type="l", 
     cex=1.25)
dev.off()

pdf("TimeSeries2ACF2b.pdf", height=5, width=6)
tACF = acf(tX2, 
           lag.max=maxTau, 
           type="correlation", 
           plot=FALSE, 
           na.action=na.fail)
plot(tACF, 
     main="",
     ci=.95, 
     ci.col="blue", 
     ci.type="white",
     type="l", 
     cex=1.25)
dev.off()

pdf("TimeSeries2ACF2bMA.pdf", height=5, width=6)
tACF = acf(tX2, 
           lag.max=maxTau, 
           type="correlation", 
           plot=FALSE, 
           na.action=na.fail)
plot(tACF, 
     main="",
     ci=.95, 
     ci.col="blue", 
     ci.type="ma",
     type="l", 
     cex=1.25)
dev.off()

#####################
#   MARKER G
#####################

# ------------------------------------------------------
#  Create a cyclic time series and 
#  synchronous cross-regression relation

tX3 <- sin(seq(0,20, by=.2)) * .5
tY3 <- 5 + .5 * tX3 + rnorm(101,mean=0, sd=.5)
tX3 <- tX3 + rnorm(101,mean=0,sd=.5)

tData3 <- data.frame(x=tX3,y=tY3)
summary(tData3)

# ------------------------------------------------------
#  Create a timeseries plot of tX3

tLen <- length(tX3)
pdf("TimeSeries2TSPlot4.pdf", height=5, width=6)
plot(c(0,tLen), c(-2,2),
     xlab="Time",
     ylab="X",
     type='n',
     pch=1)
lines(c(1:tLen), tX3, type="l", lty=1)
lines(c(0,tLen), c(0,0), type="l", lty=2)
dev.off()


# ------------------------------------------------------
#  Create an autocorrelation function plot of tX3
 
maxTau <- 100
tLen <- length(tX3)
tACF3 <- rep(NA, maxTau+1)

for (tau in (0:maxTau)) {
   t1 <- tX3[1:(tLen-tau)]
   t2 <- tX3[(1+tau):tLen]
   tSelect <- !is.na(t1) & !is.na(t2)
   if (length(t1[tSelect]) < 5) next
   tACF3[tau+1] <- cor(t1[tSelect], t2[tSelect])
}

backward <- seq(maxTau+1, 2, by= -1)
pdf("TimeSeries2ACF3.pdf", height=5, width=6)
plot(c(-maxTau,maxTau), c(-1,1), type="n",
     xlab="Lag",
     ylab="Autocorrelation")
lines(c(-maxTau:maxTau), c(tACF3[backward],tACF3), type="l", lty=1)
lines(c(-maxTau,maxTau), c(0,0), type="l", lty=2)
dev.off()

pdf("TimeSeries2ACF3b.pdf", height=5, width=6)
tACF = acf(tX3, 
           lag.max=maxTau, 
           type="correlation", 
           plot=FALSE, 
           na.action=na.fail)
plot(tACF, 
     main="",
     ci=.95, 
     ci.col="blue", 
     ci.type="white",
     type="l", 
     cex=1.25)
dev.off()


# ------------------------------------------------------
#  Create a series of State Space Plots of tX3
 
tLen <- length(tX3)

for(tTau in seq(1, 30, by=2)) {
    tFrame <- data.frame(x1=tX3[1:(tLen-tTau)], x2=tX3[(1+tTau):tLen])
    tFile <- paste("TimeSeries2SSTau", tTau, ".pdf", sep="")
    pdf(tFile, height=5, width=5)
    plot(c(-3, 3), c(-3, 3), 
         xlab="X(t)",
         ylab=paste("X(t+", tTau, ")", sep=""),
         main="",
         col="blue",
         type='n')
    lines(tFrame$x1, tFrame$x2, pch=1, type='p')
    abline(lm(x2 ~ x1, data = tFrame))
    dev.off()
}

# ------------------------------------------------------
#  Create a nearly noise-free cyclic time series and 
#  synchronous cross-regression relation

tX4 <- sin(seq(0,20, by=.2)) * .8
tY4 <- 5 + .5 * tX3 + rnorm(101,mean=0, sd=.05)
tX4 <- tX4 + rnorm(101,mean=0,sd=.05)

tData4 <- data.frame(x=tX4,y=tY4)
summary(tData4)

tLen <- length(tX4)
pdf("TimeSeries2TSPlotX4.pdf", height=5, width=6)
plot(c(0,tLen), c(-1,1),
     xlab="Time",
     ylab="X",
     type='n',
     pch=1)
lines(c(1:tLen), tX4, type="l", lty=1)
lines(c(0,tLen), c(0,0), type="l", lty=2)
dev.off()

# ------------------------------------------------------
#  Create  autocorrelation function plots of tX4
 
maxTau <- 100
tLen <- length(tX4)
tACF4 <- rep(NA, maxTau+1)

for (tau in (0:maxTau)) {
   t1 <- tX4[1:(tLen-tau)]
   t2 <- tX4[(1+tau):tLen]
   tSelect <- !is.na(t1) & !is.na(t2)
   if (length(t1[tSelect]) < 5) next
   tACF4[tau+1] <- cor(t1[tSelect], t2[tSelect])
}

backward <- seq(maxTau+1, 2, by= -1)
pdf("TimeSeries2ACF4.pdf", height=5, width=6)
plot(c(-maxTau,maxTau), c(-1,1), type="n",
     xlab="Lag",
     ylab="Autocorrelation")
lines(c(-maxTau:maxTau), c(tACF4[backward],tACF4), type="l", lty=1)
lines(c(-maxTau,maxTau), c(0,0), type="l", lty=2)
dev.off()

pdf("TimeSeries2ACF4b.pdf", height=5, width=6)
tACF = acf(tX4, 
           lag.max=maxTau, 
           type="correlation", 
           plot=FALSE, 
           na.action=na.fail)
plot(tACF, 
     main="",
     ci=.95, 
     ci.col="blue", 
     ci.type="white",
     type="l", 
     cex=1.25)
dev.off()

pdf("TimeSeries2ACF4c.pdf", height=5, width=6)
tACF = ccf(tX4, tX4, 
           lag.max=maxTau, 
           type="correlation", 
           plot=FALSE, 
           na.action=na.fail)

plot(tACF, 
     main="",
     ci=.95, 
     ci.col="blue", 
     ci.type="white",
     type="l", 
     cex=1.25)

dev.off()


# ------------------------------------------------------
#  Create a series of State Space Plots of tX4
 
tLen <- length(tX4)

for(tTau in seq(1, 30, by=2)) {
    tFrame <- data.frame(x1=tX4[1:(tLen-tTau)], x2=tX4[(1+tTau):tLen])
    tFile <- paste("TimeSeries2SS4Tau", tTau, ".pdf", sep="")
    pdf(tFile, height=5, width=5)
    plot(c(-1, 1), c(-1, 1), 
         xlab="X(t)",
         ylab=paste("X(t+", tTau, ")", sep=""),
         main="",
         col="blue",
         type='n')
    lines(tFrame$x1, tFrame$x2, pch=1, type='p')
    abline(lm(x2 ~ x1, data = tFrame))
    dev.off()
}

# ------------------------------------------------------
#  Write the timeseries matrix to a textfile for input to
#  another program.

write.csv(tData3, "tSineData1.csv", row.names=FALSE)


# ------------------------------------------------------
# End of script






