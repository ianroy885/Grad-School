# ---------------------------------------------------------------------
# Program: LDEUnivariateExample.R
#  Author: Steve Boker
#    Date: Sun May 14 14:00:20 EDT 2017
#
#
# ---------------------------------------------------------------------


# ----------------------------------
# Read libraries and set options.

options(width=80)

require(OpenMx)
require(psych)

source("GLLAfunctions.R")


# ----------------------------------
# Read the data.

tData <- as.matrix(read.table("LDEUnivariateExample.dat"))

describe(tData)

# ----------------------------------
# Set the embedding constants

embedD <- 11
theTau <- 1
deltaT <- .1


# ----------------------------------
# Plot the time series.

pdf("UnivariateTimeSeries.pdf", height=5, width=6)
plot(c(1, dim(tData)[1])*deltaT, c(-5, 5),
     xlab="Time",
     ylab="Score",
     type='n')
lines(c(1:dim(tData)[1])*deltaT, tData[,1], type='p', lwd=2, col='red')
lines(c(1, dim(tData)[1])*deltaT, c(-0, 0), type='l', lty=1, col=1)
dev.off()

# ----------------------------------
# Time-delay embed the data.

tEmbedded <- gllaEmbed(tData[,1], embed=embedD, tau=theTau, idColumn=FALSE)

dim(tEmbedded)

# ----------------------------------
# Create the fixed LDE loading matrix.

L1 <- rep(1,embedD)
L2 <- c(1:embedD)*theTau*deltaT-mean(c(1:embedD)*theTau*deltaT)
L3 <-  (L2^2)/2
LMatrix <- cbind(L1,L2,L3)

# ----------------------------------
# Create 2nd order LDE model.

manifestVars <- dimnames(tEmbedded)[[2]]

ldeModel1 <- mxModel("LDE_Model_1",
    mxMatrix("Full",  
        values=LMatrix, 
        free=FALSE, 
        name="L", 
        byrow=TRUE
    ),
    mxMatrix("Full", 3, 3, 
        values=c(  0,  0, 0,
                   0,  0, 0,
                 -.2,-.2, 0), 
        labels=c(  NA,  NA, NA,
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
    mxMatrix("Diag", embedD, embedD, 
        values=.8, 
        free=TRUE, 
        labels="uX", 
        name="U",
        lbound=0.000001
    ),
    mxMatrix("Iden", 3, name="I"),
    mxAlgebra(L %*% solve(I-A) %*% S %*% t(solve(I-A)) %*% t(L) + U, 
        name="R", 
        dimnames = list(manifestVars, manifestVars)
    ),
    mxMatrix("Full", nrow=1, ncol=length(manifestVars),
        values=1,
        free=TRUE,
        labels="mX",
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

ldeModel1Fit <- mxRun(ldeModel1)

summary(ldeModel1Fit)



# ----------------------------------
# Use mxRefModels to calculate the saturated log likelihood for chi square and RMSEA.

tRefModel <- mxRefModels(ldeModel1Fit, run=TRUE) 

summary(ldeModel1Fit, refModels=tRefModel)

