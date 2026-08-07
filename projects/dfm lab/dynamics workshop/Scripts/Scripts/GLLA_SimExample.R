# ---------------------------------------------------------------------
# Program: GLLA_SimExample170522.R  
#  Author: Steven M. Boker
#    Date: date
#
# This program runs a simulation of GLLA on a damped linear oscillator
# 
#
# ---------------------------------------------------------------------


# ---------------------------------------------
# Set up the constants and load libraries

library(deSolve)

source("GLLAfunctions.R")

# ----------------------------------
# Define the damped linear oscillator function.

DLOmodel <- function(t, prevState, parms) {
    x <- prevState[1] # x[t]
    y <- prevState[2] # dx[t]

    with(as.list(parms), 
        {
            dx <- y
            dy <- parms[1]*x + parms[2]*y + parms[3]*rnorm(1,mean=0,sd=1)
            res<-c(dx,dy)
            list(res)
        }
    )
}

# ---------------------------------------------
# Parameters for simulation.

N <- 300
eta <- -0.01
zeta <- -0.01
embedD <- 5
theOrder <- 2
dynError <- 0
tSNR <- 2.333
theTau <- 1

totalSamples <- N + 100
deltaT <- 1
totalInterval <- totalSamples * deltaT

theTimes  <- seq(0, totalInterval, length=totalSamples)  # the measurement occasions

# ----------------------------------
# Simulate a damped linear oscillator.

parms <- c(eta, zeta, dynError)
tOffsets <- c(1:N)
xstart <- c(x = 1, y = 0)
out1 <- as.data.frame(lsoda(xstart, theTimes, DLOmodel, parms))[tOffsets,]

# ----------------------------------
# Scale error for a chosen signal to noise ratio.

tSD <- sqrt(var(c(out1$x)))
tESD <- 1 / tSNR
tOscDataX <- (out1$x/tSD + rnorm(N, mean=0, sd=tESD))

# ----------------------------------
# Plot the true displacement and first derivative plus the noisy signal.

pdf("DLOExample0Error.pdf", height=5, width=6)
plot(c(min(theTimes), max(theTimes[1:N])), c(-3, 3),
     xlab="Time",
     ylab="Score",
     type='n')
lines(out1$time, out1$x/tSD, type='l', lwd=2, col='blue')
lines(out1$time, out1$y/tSD, type='l', lwd=2, col='green')
lines(out1$time, tOscDataX, type='l', lwd=2, col='red')
lines(c(min(theTimes), max(theTimes[1:N])), c(-0, 0), type='l', lty=2, col=1)
dev.off()

# ----------------------------------
# Calculate and print the W matrix.

wMatrix <- gllaWMatrix(embed=embedD, tau=theTau, deltaT=deltaT, order=theOrder)

wMatrix

# ----------------------------------
# Embed the data and calculate the derivatives.

tEmbedded <- gllaEmbed(tOscDataX, embed=embedD, tau=theTau, idColumn=FALSE)
GLLAEstimates <- tEmbedded %*% wMatrix
dimnames(GLLAEstimates) <- list(NULL, c("x", "dx/dt", "dx^2/dt^2"))

midSample <- ((embedD - 1) / 2) * theTau

# ----------------------------------
# Plot the true and estimated displacements and first derivatives.

tSelect <- 1:dim(GLLAEstimates)[1]
pdf("DLOExample0GLLAEstimates.pdf", height=5, width=6)
plot(c(min(theTimes), max(theTimes[tSelect])), c(-3, 3),
     xlab="Time",
     ylab="Score",
     type='n')
lines(out1$time[tSelect], out1$x[tSelect+midSample]/tSD, type='l', lwd=2, col='blue')
lines(out1$time[tSelect], out1$y[tSelect+midSample]/tSD, type='l', lwd=2, col='green')
lines(out1$time[tSelect], GLLAEstimates[,1], type='l', lwd=1, col='red')
lines(out1$time[tSelect], GLLAEstimates[,2], type='l', lwd=1, col='red')
lines(c(min(theTimes), max(theTimes[tSelect])), c(-0, 0), type='l', lty=2, col=1)
dev.off()
    
# ----------------------------------
# Pairs plot of the estimated derivatives.

pdf("DLOExample0GLLAPairs.pdf", height=5, width=6)
pairs(GLLAEstimates, cex=.5)
dev.off()
    
# ----------------------------------
# fit a linear model to the GLLA estimates.

cat("True eta=", eta, "  True zeta=", zeta, "\n\n", sep="")

tLM <- lm(GLLAEstimates[,3] ~ GLLAEstimates[,1] + GLLAEstimates[,2] - 1)
summary(tLM)


