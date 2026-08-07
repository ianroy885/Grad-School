# ------------------------------------------------------
# Program: TimeSeries1.R
#  Author: Steven M. Boker
#    Date: Sun May 14 11:17:59 EDT 2017
#
#  Timeseries Plots Part 1
#
# ------------------------------------------------------


# ------------------------------------------------
# Load required libraries. 

library(lattice)

# ------------------------------------
#  Read the gllaEmbedding function from disk

source("GLLAfunctions.R")


#####################
#   MARKER A
#####################

# ------------------------------------------------------
#  Create an error-free first order autoregressive time series and 
#  a synchronous cross-regression relation

tX1 <- rep(NA, 201)
tX1[1] <- 90
for (i in 2:201) {
   tX1[i] <- .975 * tX1[i-1]
}

tY1 <- 10 + .5 * tX1 + rnorm(201, mean=0, sd=5)
tX1 <- 10 + tX1 + rnorm(201, mean=0, sd=5)
tX1[tX1 <= 0] <- 0.01

tData1 <- data.frame(x=tX1, y=tY1)
summary(tData1)

pdf("TimeSeries1Hist1.pdf", height=5, width=6)
hist(tData1$x, nclass=10,
   main = "", 
   xlab = "X", 
   ylab = "Frequency",
   col="lightblue")
dev.off()

pdf("TimeSeries1QQNorm1.pdf", height=5, width=5)
print(qqmath(~ x, 
        distribution=qnorm, data=tData1,
        prepanel = prepanel.qqmathline, 
            panel = function(x, ...) {
            panel.qqmathline(x, ...)
            panel.qqmath(x, ...)
        },
        aspect=1,
        xlab = list(label="Normal Distribution",cex=1.25),
        ylab=list(label="Autoregressive Time Series",cex=1.25),
        scales=list(cex=1.25),
        main=list(label="",cex=1.25)
        )
    )
dev.off()

# ------------------------------------------------------
#  Explore Transforms of tX1
 
pdf("TimeSeries1QQNorm2.pdf", height=5, width=5)
print(qqmath(~ log(x), 
        distribution=qnorm, data=tData1,
        prepanel = prepanel.qqmathline, 
        panel = function(x, ...) {
            panel.qqmathline(x, ...)
            panel.qqmath(x, ...)
        },
        aspect=1,
        xlab = list(label="Normal Distribution",cex=1.25),
        ylab=list(label="Log of Autoregressive Time Series",cex=1.25),
        scales=list(cex=1.25),
        main=list(label="",cex=1.25)
        )
    )
dev.off()


#####################
#   MARKER B
#####################

# ------------------------------------------------------
#  Create a timeseries of tX1
 
tTime <- seq(0,1000, by=5)/60
tData1a <- data.frame(x=tX1, y=tY1, time=tTime)

pdf("TimeSeries1TSPlot1.pdf", height=5, width=6)
plot(tData1a$time, tData1a$x, 
     xlab="Time",
     ylab="X",
     type='p',
     pch=1,
     col="blue")
dev.off()

pdf("TimeSeries1TSPlot2.pdf", height=5, width=6)
plot(tData1a$time, tData1a$x, 
     xlab="Time",
     ylab="X",
     type='l',
     pch=1,
     col="blue")
dev.off()

pdf("TimeSeries1TSPlot3.pdf", height=5, width=6)
plot(tData1a$time, log(tData1a$x), 
     xlab="Time",
     ylab="log(X)",
     type='l',
     pch=1,
     col="blue")
dev.off()


#####################
#   MARKER C
#####################

# ------------------------------------------------------
#  Create a scrambled timeseries of tX1
 
tTimeRand <- order(runif(201,0,1))
tData1b <- data.frame(x=tX1, y=tY1, time=tTimeRand)

pdf("TimeSeries1QQNorm3.pdf", height=5, width=5)
print(qqmath(~ x, 
        distribution=qnorm, data=tData1b,
        prepanel = prepanel.qqmathline, 
        panel = function(x, ...) {
            panel.qqmathline(x, ...)
            panel.qqmath(x, ...)
        },
        aspect=1,
        xlab = list(label="Normal Distribution",cex=1.25),
        ylab=list(label="Shuffled Time Series",cex=1.25),
        scales=list(cex=1.25),
        main=list(label="",cex=1.25)
        )
    )
dev.off()

pdf("TimeSeries1TSPlot4.pdf", height=5, width=6)
plot(tData1b$time, tData1b$x, 
     xlab="Time",
     ylab="X",
     type='p',
     pch=1,
     col="blue")
dev.off()

#####################
#   MARKER D
#####################

# ------------------------------------------------------
#  Create a time-lagged timeseries of tX1 and plot state space
 
tau <- 5
tLen <- length(tX1)
tDataLag1 <- data.frame(x1=tX1[1:(tLen-tau)], 
                        x2=tX1[(1+tau):tLen])

pdf("TimeSeries1SSTau5.pdf", height=5, width=5)
plot(tDataLag1$x1, tDataLag1$x2, 
     xlab="X(t)",
     ylab="X(t+5)",
     main="",
     type='p',
     pch=1,
     col="blue")
dev.off()

tau <- 10
tLen <- length(tX1)
tDataLag10 <- data.frame(x1=tX1[1:(tLen-tau)], 
                         x2=tX1[(1+tau):tLen])

pdf("TimeSeries1SSTau10.pdf", height=5, width=5)
plot(tDataLag10$x1, tDataLag10$x2, 
     xlab="X(t)",
     ylab="X(t+10)",
     main="",
     type='p',
     pch=1,
     col="blue")
dev.off()

tau <- 15
tLen <- length(tX1)
tDataLag10 <- data.frame(x1=tX1[1:(tLen-tau)], 
                         x2=tX1[(1+tau):tLen])

pdf("TimeSeries1SSTau15.pdf", height=5, width=5)
plot(tDataLag10$x1, tDataLag10$x2, 
     xlab="X(t)",
     ylab="X(t+15)",
     main="",
     type='p',
     pch=1,
     col="blue")
dev.off()

tau <- 20
tLen <- length(tX1)
tDataLag20 <- data.frame(x1=tX1[1:(tLen-tau)], 
                         x2=tX1[(1+tau):tLen])

pdf("TimeSeries1SSTau20.pdf", height=5, width=5)
plot(tDataLag20$x1, tDataLag20$x2, 
     xlab="X(t)",
     ylab="X(t+20)",
     main="",
     type='p',
     pch=1,
     col="blue")
dev.off()

# ------------------------------------------------------
#  Plot state space with ordering expressed by lines
 

pdf("TimeSeries1SSTau20Lines.pdf", height=5, width=5)
plot(tDataLag20$x1, tDataLag20$x2, 
     xlab="X(t)",
     ylab="X(t+20)",
     main="",
     type='l',
     pch=1,
     col="blue")
dev.off()

#####################
#   MARKER E
#####################

# ------------------------------------------------------
#  Plot time shuffled state space
 
tau <- 20
tLen <- length(tX1)
tShuffledX <- tX1[tTimeRand]
tDataLagShuf <- data.frame(x1=tShuffledX[1:(tLen-tau)], 
                           x2=tShuffledX[(1+tau):tLen])

pdf("TimeSeries1SSShuffled.pdf", height=5, width=5)
plot(tDataLagShuf$x1, tDataLagShuf$x2, 
     xlab="X(t)",
     ylab="X(t+20) (shuffled order)",
     main="",
     type='p',
     pch=1,
     col="blue")
dev.off()

pdf("TimeSeries1SSShuffledLines.pdf", height=5, width=5)
plot(tDataLagShuf$x1, tDataLagShuf$x2, 
     xlab="X(t)",
     ylab="X(t+20) (shuffled order)",
     main="",
     type='l',
     pch=1,
     col="blue")
dev.off()

# ------------------------------------
#  Simulate a longitudinal style data set

totalSubjects <- 10
totalObservations <- 20

longMatrix <- matrix(NA, totalSubjects, totalObservations+1)

for (tSubject in 1:totalSubjects) {
    longMatrix[tSubject,1] <- tSubject
    tData <- rep(NA, totalObservations)
    tData[1] <- 80 + rnorm(1, mean=0, sd=5)
    tBeta <- .80 + rnorm(1, mean=0, sd=.05)
    for (i in 2:totalObservations) {
       tData[i] <- tBeta * tData[i-1]
    }
    tData <- 100 + rnorm(1, mean=0, sd=10) - tData
    tData  <- tData + rnorm(totalObservations, 0, 5)
    longMatrix[tSubject,2:(totalObservations+1)] <- tData
}

summary(longMatrix)

# ------------------------------------
#  Plot multiple time series from the longitudinal data set

pdf("TimeSeries1Long1.pdf", height=5, width=6)
plot(c(0,totalObservations), c(0,130), 
     xlab="Time",
     ylab="X",
     main="Multiple Time Series",
     type='n')
for (i in 1:totalSubjects){
    lines(longMatrix[i,2:(totalObservations+1)], type='l')
}
dev.off()

#####################
#   MARKER F
#####################

# ------------------------------------
#  Convert the longitudinal style data set to an mixed effects (HLM) style

mixedMatrix <- matrix(NA, totalSubjects*totalObservations, 3)

for (i in 1:totalSubjects) {
    tRows <- c(1:totalObservations) + ((i - 1) * totalObservations)
    mixedMatrix[tRows, 1] <- rep(longMatrix[i,1], totalObservations)
    mixedMatrix[tRows, 2] <- c(1:totalObservations)
    mixedMatrix[tRows, 3] <- longMatrix[i,2:(totalObservations+1)]
}

summary(mixedMatrix)

# ------------------------------------
#  Plot multiple time series from the mixed effects data set

theSubjectIDs <- longMatrix[,1]

pdf("TimeSeries1Long2.pdf", height=5, width=6)
plot(c(0,totalObservations), c(0,130), 
     xlab="Time",
     ylab="X",
     main="Multiple Time Series",
     type='n')
for (tSubject in theSubjectIDs){
    tSelect <- mixedMatrix[,1] == tSubject
    lines(mixedMatrix[tSelect,2], mixedMatrix[tSelect,3], type='l')
}
dev.off()

mixedFrame <- data.frame(id=mixedMatrix[,1],
                         occasion=mixedMatrix[,2],
                         x=mixedMatrix[,3])

summary(mixedFrame)


# ------------------------------------
#  Save the mixed effects data set to disk as a text file

write.csv(mixedFrame, file="mixedFrame.csv", row.names=FALSE)



# ------------------------------------
#  Time delay embed the mixed effects data set.

theTau <- 1

embeddedMatrix <- gllaEmbed(mixedFrame$x, embed=2, tau=theTau, groupby=mixedFrame$id, idColumn=TRUE)

# ------------------------------------
#  Plot multiple embeddings from the mixed effects data set

theSubjectIDs <- longMatrix[,1]

pdf("TimeSeries1LongSS1.pdf", height=5, width=5)
plot(c(0,150), c(0,150), 
     xlab="X(t)",
     ylab= paste("X(t+", theTau, ")", sep=""),
     main="Multiple State Spaces",
     type='n')
tColor <- 1
for (tSubject in theSubjectIDs){
    tSelect <- embeddedMatrix[,1] == tSubject
    lines(embeddedMatrix[tSelect,2], embeddedMatrix[tSelect,3], type='l', col=tColor)
    tColor <- tColor + 1
}
dev.off()


# ------------------------------------------------------
# End of script



