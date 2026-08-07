# ---------------------------------------------------------------------
# Program: ChangesInDynamics_AppendixA.R
#  Author: Steve Boker
#    Date: Wed Jan 27 13:48:57 EST 2016
#
#   Fit a univariate latent differential equation model to ExampleDataAppendixA.csv
#
# ---------------------------------------------------------------------


# ----------------------------------
# Read libraries and set options.

library(psych)              # The psych library can be found on CRAN.
library(OpenMx)             # OpenMx can be downloaded from http://openmx.psyc.virginia.edu or CRAN
source("GLLAfunctions.R")   # This has the time delay embedding function we use.

# ----------------------------------
# Read the data from a csv file (which could be output from, e.g., SPSS or Excel).

tDataFrame <- read.csv("ExampleDataAppendixA.csv")

describe(tDataFrame)

# ----------------------------------
# Plot the time series.

pdf("TimeSeriesAppendixA.pdf", height=5, width=6)
plot(c(1, dim(tDataFrame)[1]), c(-5, 5),
     xlab="Time",
     ylab="Score",
     type='n')
lines(c(1:dim(tDataFrame)[1]), tDataFrame[,1], type='l', lwd=2, col='blue')
lines(c(1, dim(tDataFrame)[1]), c(-0, 0), type='l', lty=1, col=1)
dev.off()

# ----------------------------------
# Time-delay embed the data.

embedD <- 5    # Create 5 time delay embedding columns (5 dimensional embedding)
theTau <- 1    # Use a lag of 1 when embedding
deltaT <- 1.0  # A lag of 1 is equal an elapsed time of 1.0

tEmbedded <- gllaEmbed(tDataFrame[,1], embed=embedD, tau=theTau, idColumn=FALSE)

describe(tEmbedded)

# ----------------------------------
# Create the fixed LDE loading matrix (the matrix L in the article).

L1 <- rep(1,embedD)   # The first column of L estimates the displacement
L2 <- c(1:embedD)*theTau*deltaT-mean(c(1:embedD)*theTau*deltaT)  # the second column is the 1st derivative
L3 <-  (L2^2)/2  # the third column estimates the 2nd derivative
LMatrix <- cbind(L1,L2,L3)  # bind the three columns together to make L

# ----------------------------------
# Create the 2nd order LDE model using OpenMx.

manifestVars <- dimnames(tEmbedded)[[2]]  # read the names of the manifest variable columns from tEmbedded

ldeModel1 <- mxModel("LDE_Model_1",
    # Specify the L matrix by reading in the previously defined matrix
    mxMatrix("Full",     
        values=LMatrix, 
        free=FALSE, 
        name="L", 
        byrow=TRUE
    ),
    # Specify the A matrix to have two free parameters, eta and zeta.
    mxMatrix("Full", 3, 3, 
        values=c(  0,  0, 0,
                   0,  0, 0,
                 -.2,-.2, 0), 
        labels=c(  NA,  NA,    NA,
                   NA,  NA,    NA,
                 "eta","zeta", NA), 
        free=c(FALSE,FALSE,FALSE,
               FALSE,FALSE,FALSE,
                TRUE, TRUE,FALSE), 
        name="A", 
        byrow=TRUE
    ),
    # Specify the S matrix to contain the free variances and bound them to be positive.
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
    # Specify the U matrix to be diagonal and bound to be positive.
    mxMatrix("Diag", embedD, embedD, 
        values=.8, 
        free=TRUE, 
        labels="varX", 
        name="U",
        lbound=0.000001
    ),
    # Create an appropriately sized identity matrix
    mxMatrix("Iden", 3, name="I"),
    # Specify the model-expected covariance.
    mxAlgebra(L %*% solve(I-A) %*% S %*% t(solve(I-A)) %*% t(L) + U, 
        name="R", 
        dimnames = list(manifestVars, manifestVars)
    ),
    # specify the means matrix to be fixed at 0.
    mxMatrix("Full", nrow=1, ncol=length(manifestVars),
        values=0,
        free=FALSE,
        labels="meanX",
        dimnames=list(NULL, manifestVars),
        name="M"
    ),
    # Specify the expectation and maximum likelihood fit function
    mxExpectationNormal(covariance="R", means="M"),
    mxFitFunctionML(),
    # Specify which data to use when fitting the model
    mxData(tEmbedded, type="raw")
)

# ----------------------------------
# Fit the LDE model and examine the summary results.

ldeModel1Fit <- mxRun(ldeModel1)

summary(ldeModel1Fit)


# ----------------------------------
# Create an LDE model saturated at the latents.

manifestVars <- dimnames(tEmbedded)[[2]]

ldeModelNull <- mxModel("LDE_Model_null",
    mxMatrix("Full",  
        values=LMatrix, 
        free=FALSE, 
        name="L", 
        byrow=TRUE
    ),
    mxMatrix("Full", 3, 3, 
        values=c(  0,  0, 0,
                   0,  0, 0,
                   0,  0, 0), 
        free=c(FALSE,FALSE,FALSE,
               FALSE,FALSE,FALSE,
                FALSE, FALSE,FALSE), 
        name="A", 
        byrow=TRUE
    ),
    mxMatrix("Symm", 3, 3,
        values=c(  .8,
                    -.01, .8,
                    -.01, -.01, .8), 
        free=c( TRUE,
               TRUE,  TRUE,
               TRUE, TRUE, TRUE), 
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
    # Specify the expectation and maximum likelihood fit function
    mxExpectationNormal(covariance="R", means="M"),
    mxFitFunctionML(),
    # Specify which data to use when fitting the model
    mxData(tEmbedded, type="raw")
)

# ----------------------------------
# Fit the null LDE model and examine the summary results.

ldeModelNullFit <- mxRun(ldeModelNull)

summary(ldeModelNullFit)

# ----------------------------------
# Calculate a pseudo-rsquare for the second derivative.

1 - (mxEval(S[3,3], ldeModel1Fit) / mxEval(S[3,3], ldeModelNullFit))

