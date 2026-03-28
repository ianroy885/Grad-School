# ---------------------------------------------------------------------
# Program: ChangesInDynamics_AppendixD.R
#  Author: Steve Boker
#    Date: Wed Jan 27 13:49:34 EST 2016
#
#   Fit a univariate latent differential equation model with
#   individual differences in equilibria intercept and slope,
#   and individual differences in dynamics to ExampleDataAppendixC.csv
#
# ---------------------------------------------------------------------


# ----------------------------------
# Read libraries and set options.

library(psych)              # The psych library can be found on CRAN.
library(OpenMx)             # OpenMx can be downloaded from http://openmx.psyc.virginia.edu or CRAN
source("GLLAfunctions.R")   # This has the time delay embedding function we use.
mxOption(NULL, 'Number of Threads', parallel::detectCores())
mxOption(NULL, 'Default optimizer', 'NPSOL')

# ----------------------------------
# Read the data from a csv file (which could be output from, e.g., SPSS or Excel).

tDataFrame <- read.csv(file="ExampleDataAppendixD.csv", header=TRUE)

describe(tDataFrame) # We will only use the ID, X, and Occasion variables.

numSubjects <- length(unique(tDataFrame$ID))

# We will assume that all individuals are measured the same number of occasions as the 
# first individual in the data frame.
numOccasions <- length(tDataFrame$ID[tDataFrame$ID==unique(tDataFrame$ID)[1]])

# We need to convert subject IDs into an index from 1:numSubjects.
# The example csv file has this characteristic, but your data may have another ID scheme.
indexID <- as.numeric(as.character(as.factor(tDataFrame$ID)))

# ----------------------------------
# Plot the time series for the first 10 individuals in the data.

pdf("TimeSeriesAppendixD.pdf", height=5, width=6)
plot(c(1, numOccasions), c(-15, 15),
     xlab="Sample Index",
     ylab="Score",
     type='n')
for(tID in unique(tDataFrame$ID)[1:10]) {
    tSel <- tDataFrame$ID==tID
    lines(c(1:numOccasions), tDataFrame$X[tSel], type='l', lwd=2, col=rainbow(10)[tID])
}
lines(c(1,numOccasions), c(-0, 0), type='l', lty=1, col=1)
dev.off()


# ----------------------------------
# Set constants for embedding and create the embedded data matrix.

embedD <- 5    # Create 5 time delay embedding columns (5 dimensional embedding)
theTau <- 1    # Use a lag of 1 when embedding
deltaT <- 1.0  # A lag of 1 is equal an elapsed time of 1.0
centerEmbed <- floor((embedD + 1) / 2)  # which is the middle column of the embedded matrix

# We need to keep track of which data belongs to which person and thus we use the
#   indexID we created earlier to group the data that belongs to each individual.
# This embedding function assumes that each individual's Occasions are presorted
#   to be ordered from earliest to most recent.

tEmbeddedData <- gllaEmbed(tDataFrame$X, embed=embedD, tau=theTau, groupby=indexID, label="X", idColumn=TRUE)

# ----------------------------------
# Create dataframe for LDE Estimation.

# create a matrix to hold the embedded data and the occasion
ldeData <- cbind(tEmbeddedData, matrix(NA, dim(tEmbeddedData)[1], 1))
dimnames(ldeData) <- list(NULL, c(dimnames(tEmbeddedData)[[2]], "Occasion"))

embedCols <- dim(tEmbeddedData)[2]

# For each individual, we store the Occasion number of the middle column of each
#   row into the Occasion variable for that row.
# This again assumes that Occasions are presorted into ascending order for each individual.
for(tID in unique(tDataFrame$ID)) {
    covIndex <- c(centerEmbed:(centerEmbed-1+length(ldeData[ldeData[,1]==tID,1])))
    ldeData[ldeData[,1]==tID, embedCols+1] <- tDataFrame$Occasion[tDataFrame$ID==tID][covIndex]
}

ldeFrame <- data.frame(ldeData)
describe(ldeFrame)

# ----------------------------------
# Create the fixed LDE loading matrix.

L1 <- rep(1,embedD)
L2 <- c(1:embedD)*theTau*deltaT-mean(c(1:embedD)*theTau*deltaT)
L3 <-  (L2^2)/2
LMatrix <- cbind(L1,L2,L3)

# ----------------------------------
# Create the fixed LGC loading matrix to estimate both intercept and slope.

# A row of 1's to estimate the intercept.
LGC1 <- rep(1,embedD)
# A basis function increasing by (theTau*deltaT) to estimate the slope in the selected time units.
LGC2 <- (c(1:embedD)*theTau*deltaT-mean(c(1:embedD)*theTau*deltaT))
# Create the LGC matrix for later use in the mxModel.
LGCMatrix <- rbind(LGC1,LGC2)
# Create a matrix to expand the time-varying occasion to add to the LGCMatrix.
LGCMatrix2 <- rbind(rep(0,embedD), rep(deltaT,embedD))

# ----------------------------------
# Create a random intercepts and slopes model.

embeddedVars <- dimnames(tEmbeddedData)[[2]][2:embedCols]  # The names of the embedded variables
latentVars <- c("X", "dX", "d2X")  # The names of the latent variables
manifestVars <- c(embeddedVars)  # The names of the manifest variables
numManifestVars <- length(manifestVars)
numLatentVars <- length(latentVars)
numTotalVars <- numManifestVars + numLatentVars
# We recreate numSubjects in case some subject(s) might not have sufficient non-missing data.
numSubjects <- length(unique(ldeFrame$ID))  

multilevelLDE3 <- mxModel("multilevelLDE3", 
    # The Rand matrix holds one row for each indexID
    # The first column holds each individual's equilibrium intercept.
    # The second column holds each individual's equilibrium slope.
    # We add a little noise to the starting values to help the estimation begin.
    mxMatrix("Full", 
        values=matrix(c(0, 0)+rnorm(numSubjects*2,sd=0.1), nrow=numSubjects, ncol=4),
        free=TRUE, 
        byrow=TRUE,
        dimnames=list(NULL, c("randInt", "randSlope", "randEta", "randZeta")),
        name="Rand",
    ),
    # Create an mxAlgebra that will select an individual's row in the Rand matrix
    mxAlgebra(Rand[data.ID,], name="RandRow"),    
    # Specify the A matrix to have two free parameters, eta and zeta.
    mxMatrix("Full", nrow=numLatentVars, ncol=numLatentVars, 
        values=c(  0,  0, 0,
                   0,  0, 0,
                 -.2,-.2, 0), 
        labels=c(NA, NA, NA,
                 NA, NA, NA,
                "RandRow[1,3]", "RandRow[1,4]", NA),
        free=FALSE,
        byrow=TRUE,
        name="A"
    ),
    # Specify the S matrix to contain the free variances and bound them to be positive.
    mxMatrix("Symm", nrow=numLatentVars, ncol=numLatentVars,
        values=c(.8,
                  0,.8,
                  0, 0,.8),
        free=c(TRUE,
               FALSE,TRUE,
               FALSE,FALSE,TRUE),
        labels=c("VarX",
                 NA, "VarDX",
                 NA,NA, "VarD2X"), 
        lbound=c(0.00000001,
                 NA, 0.00000001,
                 NA, NA, 0.00000001),
        byrow=TRUE,
        name="S" 
    ),
    # Specify the U matrix to be diagonal and bound to be positive.
    mxMatrix("Diag", nrow=numManifestVars, ncol=numManifestVars, 
        values=.8, 
        free=TRUE, 
        labels="uX", 
        name="U",
        lbound=0.000001
    ),
    # Create an appropriately sized identity matrix
    mxMatrix("Iden", nrow=numLatentVars, name="I"),
    # Specify the L matrix by reading in the previously defined matrix
    mxMatrix("Full",  
        values=LMatrix, 
        free=FALSE, 
        byrow=TRUE,
        name="L" 
    ),
    # Specify the model-expected covariance.
    mxAlgebra(L %*% solve(I-A) %*% S %*% t(solve(I-A)) %*% t(L) + U, 
        dimnames = list(manifestVars, manifestVars),
        name="R" 
    ),
    # Specify the C matrix by reading in the previously defined fixed LGC matrix
    mxMatrix("Full", values=LGCMatrix, name="C"),
    # Specify the H matrix by reading in the previously defined fixed LGC filter matrix
    mxMatrix("Full", values=LGCMatrix2, name="H"),
    # Specify the J matrix to 
    mxMatrix("Full", nrow=2, ncol=2, 
            free=FALSE,
            labels=c(NA, NA,
                     NA, "data.Occasion"), 
            byrow=TRUE,
            name="J"),
    mxMatrix("Full", nrow=1, ncol=2,
        values=0, free=FALSE, labels=c("RandRow[1,1]", "RandRow[1,2]"), name="Mu"
    ),
    mxAlgebra(Mu %*% ((J %*% H) + C), 
        dimnames=list(NULL, embeddedVars), name="M"
    ),
    # Specify the expectation and maximum likelihood fit function
    mxExpectationNormal(covariance="R", means="M"),
    mxFitFunctionML(),
    # Specify which data to use when fitting the model
    mxData(ldeFrame, type="raw")
)

# ----------------------------------
# Fit the  model and examine the summary results.

multilevelLDE3Fit <- mxRun(multilevelLDE3)

summary(multilevelLDE3Fit)

# ----------------------------------
# Examine the individual differences results.

RandOut <- mxEval(Rand, model=multilevelLDE3Fit)

var(RandOut)
apply(RandOut, 2, mean)

# Since this is a simulation, we can check how well we estimated eta and zeta.
tSel <- tDataFrame$Occasion==0
round(mean(tDataFrame$TrueEta[tSel]),4)
round(mean(tDataFrame$TrueZeta[tSel]),4)

round(var(cbind(tDataFrame$TrueEta[tSel], tDataFrame$TrueZeta[tSel])), 8)

# ----------------------------------
# Plot the time series and estimated equilibria for the first 10 individuals in the data.

pdf("TimeSeriesAppendixD-fitted.pdf", height=5, width=6)
plot(c(1, numOccasions), c(-15, 15),
     xlab="Sample Index",
     ylab="Score",
     type='n')
for(tID in unique(tDataFrame$ID)[1:10]) {
    tSel <- tDataFrame$ID==tID
    lines(c(1:numOccasions), tDataFrame$X[tSel], type='l', lwd=2, col=rainbow(10)[tID])
    lines(c(0,numOccasions), c(RandOut[tID,1], RandOut[tID,1]+numOccasions*RandOut[tID,2]), type='l', lwd=2, col=rainbow(10)[tID])
}
lines(c(1,numOccasions), c(-0, 0), type='l', lty=1, col=1)
dev.off()

