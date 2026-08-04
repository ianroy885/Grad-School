# ---------------------------------------------------------------------
# Program: LDEMultivariateFIMLExample.R
#  Author: Steve Boker
#    Date: Sun May 14 14:03:15 EDT 2017
#
#
# ---------------------------------------------------------------------


# ----------------------------------
# Read libraries and set options.

options(width=110)

library(OpenMx)
library(psych)

source("GLLAfunctions.R")

mxOption(NULL, 'Number of Threads', parallel::detectCores())

# ----------------------------------
# Read the data.

tData <- as.matrix(read.table("LDEMultivariateExample.dat"))

describe(tData)

# ----------------------------------
# Set the embedding constants

embedD <- 11
theTau <- 1
deltaT <- 0.1

# ----------------------------------
# Plot the time series.

pdf("MultivariateTimeSeries.pdf", height=5, width=6)
plot(c(1, dim(tData)[1])*deltaT, c(-2, 2),
     xlab="Time",
     ylab="Score",
     type='n')
lines(c(1:dim(tData)[1])*deltaT, tData[,1], type='p', lwd=2, col='red')
lines(c(1:dim(tData)[1])*deltaT, tData[,2], type='p', lwd=2, col='green')
lines(c(1:dim(tData)[1])*deltaT, tData[,3], type='p', lwd=2, col='blue')
lines(c(1:dim(tData)[1])*deltaT, tData[,4], type='p', lwd=2, col='purple')
lines(c(1, dim(tData)[1])*deltaT, c(-0, 0), type='l', lty=1, col=1)
dev.off()

# ----------------------------------
# Time-delay embed the data.

numIndicators <- 4

tEmbedded <- cbind(gllaEmbed(tData[,1], embed=embedD, tau=theTau, label="x", idColumn=FALSE),
                   gllaEmbed(tData[,2], embed=embedD, tau=theTau, label="y", idColumn=FALSE),
                   gllaEmbed(tData[,3], embed=embedD, tau=theTau, label="z", idColumn=FALSE),
                   gllaEmbed(tData[,4], embed=embedD, tau=theTau, label="w", idColumn=FALSE))

# ----------------------------------
# Create the fixed LDE loading matrix.

L1 <- rep(1,embedD)
L2 <- c(1:embedD)*theTau*deltaT-mean(c(1:embedD)*theTau*deltaT)
L3 <-  (L2^2)/2

LMatrix <- cbind(L1,L2,L3)

# ----------------------------------
# Create a 2nd order Multivariate LDE model.

manifestVars <- dimnames(tEmbedded)[[2]]

ldeModelMulti1 <- mxModel("LDE_Model_1",
    mxMatrix("Full", numIndicators, 1,
        values=c(1, .2, .2, .2), 
        free=c(FALSE, TRUE, TRUE, TRUE), 
        labels=c(NA, "a", "b", "c"),
        name="LFree", 
        byrow=TRUE
    ),
    mxMatrix("Full",  
        values=LMatrix, 
        free=FALSE, 
        name="LFixed", 
        byrow=TRUE
    ),
    mxAlgebra(LFree %x% LFixed, name="L"),
    mxMatrix("Full", 3, 3, 
        values=c(0,  0, 0,
                 0,  0, 0,
                 -.2,-.2, 0), 
        labels=c(NA,  NA, NA,
                 NA,  NA, NA,
                 "eta","zeta", NA), 
        free=c(FALSE,FALSE,FALSE,
               FALSE,FALSE,FALSE,
                TRUE, TRUE,FALSE), 
        name="A", 
        byrow=TRUE
    ),
    mxMatrix("Symm", 3, 3,
        values=c(  .8,
                    0, .8,
                    0,  0, .8), 
        free=c( TRUE,
               FALSE,  TRUE,
               FALSE, FALSE, TRUE), 
        labels=c("Vx",
                 NA, "Vdx",
                 NA, NA, "Vd2x"), 
        name="S", 
        byrow=TRUE,
        lbound=c(0.00000001,
                 NA, 0.00000001,
                 NA, NA, 0.00000001)
    ),
    mxMatrix("Diag", embedD*numIndicators, embedD*numIndicators, 
        values=.8, 
        free=TRUE, 
        labels=c(rep("uX", embedD), rep("uY", embedD), rep("uZ", embedD), rep("uW", embedD)), 
        name="U",
        lbound=0.000001
    ),
    mxMatrix("Iden", 3, name="I"),
    mxAlgebra(L %*% solve(I-A) %*% S %*% t(solve(I-A)) %*% t(L) + U, 
        name="R", 
        dimnames = list(manifestVars, manifestVars)
    ),
    mxMatrix("Full", nrow=1, ncol=embedD*numIndicators,
        values=1,
        free=TRUE,
        labels=c(rep("mX", embedD), rep("mY", embedD), rep("mZ", embedD), rep("mW", embedD)),
        dimnames=list(NULL, manifestVars),
        name="M"
    ),
    mxExpectationNormal(covariance="R", means="M"),
    mxFitFunctionML(),
    mxData(tEmbedded, 
        type="raw"
    )
)

# ----------------------------------
# Fit the LDE model and examine the summary results.

ldeModel1MultiFit <- mxRun(ldeModelMulti1)

summary(ldeModel1MultiFit)


# ----------------------------------
# Use mxRefModels to calculate the saturated log likelihood for chi square and RMSEA.

tRefModel <- mxRefModels(ldeModel1MultiFit, run=TRUE) 

summary(ldeModel1MultiFit, refModels=tRefModel)


