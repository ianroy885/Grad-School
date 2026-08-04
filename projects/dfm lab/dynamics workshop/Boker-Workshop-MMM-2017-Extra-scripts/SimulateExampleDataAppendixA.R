# ---------------------------------------------------------------------
# Program: SimulateExampleDataAppendixA.R
#  Author: Steve Boker
#    Date: Wed Jan 27 13:58:17 EST 2016
#
#
# ---------------------------------------------------------------------


# ----------------------------------
# Read libraries and set options.

library(deSolve)
library(OpenMx)
library(psych)

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
eta <- -0.05
zeta <- -0.02
embedD <- 9
theOrder <- 2
dynError <- 0
tSNR <- 1
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
# Write a univariate data vector.

write.csv(tOscDataX, file="ExampleDataAppendixA.csv", row.names=FALSE)

