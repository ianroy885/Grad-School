# ---------------------------------------------------------------------
# Program: CoupledLDEUnivariateExample.R
#  Author: Steve Boker
#    Date: Sun May 14 14:16:34 EDT 2017
#
#
# ---------------------------------------------------------------------


# ----------------------------------
# Read libraries and set options.

options(width=110)

library(OpenMx)
library(psych)

source("GLLAfunctions.R")
mxOption(NULL, 'Default optimizer', 'NPSOL')

# ----------------------------------
# Read the data.

tData <- read.csv("CoupledLDE.csv", header=TRUE)
describe(tData)

# ----------------------------------
# Plot the time series.

pdf("CoupledTimeSeries.pdf", height=5, width=6)
plot(c(1, dim(tData)[1]), c(-10, 10),
     xlab="Time",
     ylab="Score",
     type='n')
lines(c(1:dim(tData)[1]), tData[,1], type='p', lwd=2, col='red')
lines(c(1:dim(tData)[1]), tData[,2], type='p', lwd=2, col='blue')
lines(c(1, dim(tData)[1]), c(-0, 0), type='l', lty=1, col=1)
dev.off()



# ----------------------------------
# Time-delay embed the data.

embedD <- 7
theTau <- 1
deltaT <- 1

numIndicators <- 2

# ----------------------------------
# Time-delay embed the data.

tEmbedded <- cbind(gllaEmbed(tData[,1], embed=embedD, tau=theTau, label="x", idColumn=FALSE),
                   gllaEmbed(tData[,2], embed=embedD, tau=theTau, label="y", idColumn=FALSE))

# ----------------------------------
# Create the fixed LDE loading matrix.

L1 <- rep(1,embedD)
L2 <- c(1:embedD)*theTau*deltaT-mean(c(1:embedD)*theTau*deltaT)
L3 <-  (L2^2)/2
LMatrix <- cbind(L1,L2,L3)

# ----------------------------------
# Create a 2nd order Multivariate LDE model.

manifestVars <- dimnames(tEmbedded)[[2]]

ldeModelCoupled1 <- mxModel("LDE_Coupled_Model_1",
    mxMatrix("Iden", 2, name="I2"),
    mxMatrix("Full",  
        values=LMatrix, 
        free=FALSE, 
        name="LFixed", 
        byrow=TRUE
    ),
    mxAlgebra(I2 %x% LFixed, name="L"),
    mxMatrix("Full", 6, 6, 
        values=c(  0,  0,  0,  0,  0,  0,
                   0,  0,  0,  0,  0,  0,
                 -.2,-.2,  0, -.1, .1,  0,
                   0,  0,  0,  0,  0,  0,
                   0,  0,  0,  0,  0,  0,
                  -.1, .1,  0,-.2,-.2,  0), 
        labels=c(    NA,     NA,     NA,     NA,     NA,     NA,
                     NA,     NA,     NA,     NA,     NA,     NA,
                 "etaX","zetaX",     NA,"gammaX1","gammaX2", NA,
                     NA,     NA,     NA,     NA,     NA,     NA,
                     NA,     NA,     NA,     NA,     NA,     NA,
             "gammaY1","gammaY2",    NA, "etaY","zetaY",     NA), 
        free=c( F,F,F,F,F,F,
                F,F,F,F,F,F,
                T,T,F,T,T,F,
                F,F,F,F,F,F,
                F,F,F,F,F,F,
                T,T,F,T,T,F), 
        ubound=c(    NA,     NA,     NA,     NA,     NA,     NA,
                     NA,     NA,     NA,     NA,     NA,     NA,
                  -.001,     NA,     NA,     .4,     .4,     NA,
                     NA,     NA,     NA,     NA,     NA,     NA,
                     NA,     NA,     NA,     NA,     NA,     NA,
                     .4,     .4,     NA,  -.001,     NA,     NA), 
        lbound=c(    NA,     NA,     NA,     NA,     NA,     NA,
                     NA,     NA,     NA,     NA,     NA,     NA,
                    -.5,     NA,     NA,    -.4,    -.4,     NA,
                     NA,     NA,     NA,     NA,     NA,     NA,
                     NA,     NA,     NA,     NA,     NA,     NA,
                     -.4,   -.4,     NA,    -.5,     NA,     NA), 
        name="A", 
        byrow=TRUE
    ),
    mxMatrix("Symm", 6, 6,
        values=c(  .8,
                    0, .8,
                    0, 0, .8,
                    -.1, -.1, 0, .8,
                    -.1, -.1, 0, 0, .8,
                    0, 0, 0, 0, 0, .8), 
        free=c( T,
                F, T,
                F, F, T,
                T, T, F, T,
                T, T, F, F, T,
                F, F, F, F, F, T), 
        labels=c("VX",
                 NA, "VdX",
                 NA, NA, "Vd2X",
                 NA, NA, NA, "VY",
                 NA, NA, NA, NA, "VdY",
                 NA, NA, NA, NA, NA, "Vd2Y"), 
        name="S", 
        byrow=TRUE,
        lbound=c(0.00000001,
                 NA, 0.00000001,
                 NA, NA, 0.00000001,
                 NA, NA, NA, 0.00000001,
                 NA, NA, NA, NA, 0.00000001,
                 NA, NA, NA, NA, NA, 0.00000001)
    ),
    mxMatrix("Diag", embedD*numIndicators, embedD*numIndicators, 
        values=.8, 
        free=TRUE, 
        labels=c(rep("uX", embedD), rep("uY", embedD)), 
        name="U",
        lbound=0.000001
    ),
    mxMatrix("Iden", 6, name="I"),
    mxAlgebra(L %*% solve(I-A) %*% S %*% t(solve(I-A)) %*% t(L) + U, 
        name="R", 
        dimnames = list(manifestVars, manifestVars)
    ),
    mxExpectationNormal(covariance="R"),
    mxFitFunctionML(),
    mxData(cov(tEmbedded), 
        type="cov", 
        numObs=dim(tEmbedded)[1]
    )
)

# ----------------------------------
# Fit the LDE model and examine the summary results.

ldeModel1CoupledFit <- mxRun(ldeModelCoupled1)

summary(ldeModel1CoupledFit)

