# ---------------------------------------------------------------------
# Program: SimulateExampleDataAppendixB.R
#  Author: Steve Boker
#    Date: Wed Jan 27 13:58:24 EST 2016
#
# This program simulates individual differences in equilibrium intercept. 
#
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# Variables 
# ---------------------------------------------------------------------
# ID -- The subject ID
# Occasion -- The occasion of measurement
# X  -- The variable of interest
# Z1 -- The first demographic covariate
# Z2 -- The second demographic covariate
# Z3 -- The third demographic covariate
# Z4 -- The fourth demographic covariate
# V1 -- A time-varying covariate
# ---------------------------------------------------------------------

# ----------------------------------
# Read libraries and set options.

options(width=110)

library(deSolve)
library(psych)
library(OpenMx)

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
# SIMULATE RANDOM ETA AND ZETA

# ---------------------------------------------
# Parameters for simulation.

P <- 50
N <- 50
meanEta <- -0.2
meanZeta <- -0.1
theOrder <- 2
dynError <- 0
tSNR <- 1
totalSamples <- P + 100
deltaT <- 1
totalInterval <- totalSamples * deltaT

theTimes  <- seq(0, totalInterval, length=totalSamples)  # the measurement occasions

sdInt <- 2
sdSlope <- .3

# ----------------------------------
# Create the data matrix in which to store the simulated data

simDataMatrix <- matrix(NA, N*P, 10)
dimnames(simDataMatrix) <- list(NULL, c("ID", "Occasion", "X", "Z1", "Z2", "Z3", "TrueInt", "TrueSlope", "TrueEta", "TrueZeta"))

# ----------------------------------
# Create the covariates

z4 <- rnorm(N, mean=0, sd=1)
z2 <- rnorm(N, mean=0, sd=1)
z3 <- rnorm(N, mean=0, sd=1)
z1 <- .5 * z4 + rnorm(N, mean=0, sd=1)
z1 <- z1 / sqrt(var(z4))

# ----------------------------------
# Loop to create multiple subjects' data

i <- 1

for(tID in 1:N) {
    
    # ----------------------------------
    # Simulate a damped linear oscillator.

    eta <- meanEta + rnorm(1, mean=0, sd=.1)
    zeta <- meanZeta + rnorm(1, mean=0, sd=.1)
    randInt <-  rnorm(1, mean=0, sd=sdInt)
    randSlope <-  0

    parms <- c(eta, zeta, dynError)
    randOffset <- runif(1, min=0, max=41)
    tOffsets <- c(1:P) + randOffset
    xstart <- c(x = 1, y = 0)
    out1 <- as.data.frame(lsoda(xstart, theTimes, DLOmodel, parms))[tOffsets,]

    # ----------------------------------
    # Scale error for a chosen signal to noise ratio.

    tSD <- sqrt(var(c(out1$x)))
    tESD <- 1 / tSNR
    tOscDataX <- (out1$x/tSD + rnorm(P, mean=0, sd=tESD))

    # ----------------------------------
    # Save the data into the matrix.

	v1 <- rnorm(length(tOscDataX))
    simDataMatrix[i:(i+P-1),1] <- tID
    simDataMatrix[i:(i+P-1),2] <- theTimes[1:P]
    simDataMatrix[i:(i+P-1),3] <- tOscDataX  + randInt  + randSlope*(0:(P-1)) 
    simDataMatrix[i:(i+P-1),4] <- z1[tID] 
    simDataMatrix[i:(i+P-1),5] <- z2[tID] 
    simDataMatrix[i:(i+P-1),6] <- z3[tID] 
    simDataMatrix[i:(i+P-1),7] <- randInt 
    simDataMatrix[i:(i+P-1),8] <- randSlope
    simDataMatrix[i:(i+P-1),9] <- eta 
    simDataMatrix[i:(i+P-1),10] <- zeta 

    i <- i + P
}

simDataFrame <- data.frame(simDataMatrix)
describe(simDataFrame)

# ----------------------------------
# Plot the noisy signal plus the random slopes and intercepts.

pdf("ExampleDataAppendixB.pdf", height=5, width=6)
plot(c(1, P), c(-10, 10),
     xlab="Sample Index",
     ylab="Score",
     type='n')
for(tID in unique(simDataFrame$ID)) {
    tSel <- simDataFrame$ID==tID
    lines(c(1:P), simDataFrame$X[tSel], type='l', lwd=2, col='blue')
}
lines(c(1,P), c(-0, 0), type='l', lty=1, col=1)
dev.off()

# ----------------------------------
# Write out the data frame.

write.csv(simDataFrame, file="ExampleDataAppendixB.csv", row.names=FALSE)

