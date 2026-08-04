# ---------------------------------------------------------------------
# Program: VectorFields1.R
#  Author: Steven M. Boker
#    Date: Sun May 14 12:10:55 EDT 2017
#
#  Vector Field plotting
#
# ---------------------------------------------------------------------

#####################
#   MARKER A
#####################

# ------------------------------------------------------
#  Create and plot two autoregressive time series 

tX1 <- rep(NA, 201)
tX1[1] <- 90
for (i in 2:201) {
   tX1[i] <- .98 * tX1[i-1]
}

tX1 <- 100 - tX1 + rnorm(201, mean=0, sd=10)

# Uncomment the next line if you want to include a ceiling effect.
#tX1[tX1 > 100] <- 100

tX2 <- rep(NA, 201)
tX2[1] <- 90
for (i in 2:201) {
   tX2[i] <- .99 * tX2[i-1]
}

tX2 <- 100 - tX2 + rnorm(201, mean=0, sd=10)
#tX2[tX2 > 100] <- 100

tX3 <- rep(NA, 201)
tX3[1] <- 90
for (i in 2:201) {
   tX3[i] <- .995 * tX3[i-1]
}

tX3 <- 100 - tX3 + rnorm(201, mean=0, sd=10)
#tX3[tX3 > 100] <- 100

pdf("VF1ARLongPlot1.pdf", height=5, width=6)
plot(c(0,70), c(0,120), type="n" ,
     xlab="Age",
     ylab="Total Score")
lines(c(0:200)/3, tX1, type="l", lty=1, col=1)
lines(c(0:200)/3, tX2, type="l", lty=1, col=2)
lines(c(0:200)/3, tX3, type="l", lty=1, col=4)
dev.off()


# ------------------------------------------------------
#  Plot state space plots of the Autoregressive Data
#  with two waves of measurement spaced 15 days apart


for (theTau in seq(1, 56, by=5)) {
    tfile <- paste("VF1AR1StateSpaceTau", theTau, ".pdf", sep="")
    tseq <- c(1:(length(tX1)-(theTau+1)))
    tEmbedFrame <- data.frame(x1=tX1[tseq],
                              x1tau=tX1[tseq+theTau],
                              x2=tX2[tseq],
                              x2tau=tX2[tseq+theTau],
                              x3=tX3[tseq],
                              x3tau=tX3[tseq+theTau]
                              )
    pdf(tfile, height=5, width=5)
    plot(c(0,120), c(0,120), type="n" ,
         xlab="Test Score (t)",
         ylab=paste("Test Score (t+", theTau, ")", sep=""))
    lines(tEmbedFrame$x1, tEmbedFrame$x1tau, type="p", col=1)
    abline(lm(x1tau ~ x1, data = tEmbedFrame), col=1)
    lines(tEmbedFrame$x2, tEmbedFrame$x2tau, type="p", col=2)
    abline(lm(x2tau ~ x2, data = tEmbedFrame), col=2)
    lines(tEmbedFrame$x3, tEmbedFrame$x3tau, type="p", col=4)
    abline(lm(x3tau ~ x3, data = tEmbedFrame), col=4)
    dev.off()
}

#####################
#   MARKER B
#####################

# ------------------------------------------------------
#  Create a timelagged dataframe (3-D Embedding)

tau <- 15
tX <- c(tX1, rep(NA,2*tau+1), tX2, rep(NA,2*tau+1), tX3)
tGroup <- c(rep(1, length(tX1)), rep(NA,2*tau+1), 
            rep(2, length(tX2)), rep(NA,2*tau+1),
            rep(4, length(tX3)), rep(NA,2*tau+1))
tAge <- c(1:length(tX1), rep(NA,2*tau+1), 
           1:length(tX2), rep(NA,2*tau+1), 
           1:length(tX3))/3
tLen <- length(tX)
ageData1 <- data.frame(x1=tX[1:(tLen-(2*tau))], 
                          x2=tX[(1+tau):(tLen-tau)],
                          x3=tX[(1+(2*tau)):tLen],
                          age=tAge[1:(tLen-(2*tau))],
                          group=tGroup[1:(tLen-(2*tau))])


# ------------------------------------------------------
#  Plot the first two columns as if they were a panel study with two
#  waves of measurement spaced 15 days apart

pdf("VF1ARLagged1.pdf", height=5, width=6)
plot(c(0,220)/3, c(0,120), type="n" ,
     xlab="Age",
     ylab="Total Score")
for (i in seq(1,length(ageData1$x1),by=2)) {
   if (is.na(ageData1$age[i]) | 
       ageData1$age[i] > (78-tau)) {
          next
   }
   lines(c(ageData1$age[i],ageData1$age[i]+tau), 
         c(ageData1$x1[i],ageData1$x2[i]), type="l", lty=1)
}
dev.off()


#####################
#   MARKER C
#####################

# ------------------------------------------------------
#  Plot the three columns as if they were a panel study with three
#  waves of measurement spaced 15 days apart

pdf("VF1ARLong2.pdf", height=5, width=6)
plot(c(0,220)/3, c(0,120), type="n" ,
     xlab="Age",
     ylab="Total Score")
for (i in seq(1,length(ageData1$x1),by=2)) {
   if (is.na(ageData1$age[i]) | 
       ageData1$age[i] > (78-(2*tau))) {
          next
   }
   lines(c(ageData1$age[i],
           ageData1$age[i]+tau,
           ageData1$age[i]+(2*tau)), 
         c(ageData1$x1[i],
           ageData1$x2[i],
           ageData1$x3[i]), type="l", lty=1)
}
dev.off()



# ------------------------------------------------------
#  Calculate the slope and midpoint for each pair from the
#  first two columns of ageData1

ageMid <- ageData1$age + (tau/2)
scoreMid <- (ageData1$x1 + ageData1$x2) / 2
slopeMid <- (ageData1$x2 - ageData1$x1) / tau
tSelect <- !is.na(ageMid) & !is.na(scoreMid) & !is.na(slopeMid)
ageSlopeFrame <- data.frame(
                 age = ageMid[tSelect],
                 group = ageData1$group[tSelect],
                 score = scoreMid[tSelect],
                 slope = slopeMid[tSelect])


#####################
#   MARKER D
#####################

# ------------------------------------------------------
#  Calculate a loess smoothed slope surface

slopeLoess <- loess(slope ~ age * score, data=ageSlopeFrame)

ageSlopeFrame2 <- data.frame(age=ageSlopeFrame$age,
                   score=ageSlopeFrame$score,
                   slope=ageSlopeFrame$slope,
                   fitted=predict(slopeLoess))



#####################
#   MARKER E
#####################

# ------------------------------------------------------
#  Calculate a slope field of the fitted slopes

smoothX <- seq(0, 75, by=6)
smoothY <- seq(0, 120, by=8)
smoothP <- matrix(NA, length(smoothX), length(smoothY))
smoothZ <- matrix(NA, length(smoothX), length(smoothY))
h <- 10
for (i in 1:length(smoothX)) {
   x <- smoothX[i]
   t1x <- abs(x - ageSlopeFrame2$age)/h
   for (j in 1:length(smoothY)) {
      y <- smoothY[j]
      t1y <- abs(y - ageSlopeFrame2$score)/h
      t2 <- rep(0, length(t1x))
      t2[t1x < 1 & t1y < 1] <- 1/2
      smoothP[i,j] <- sum(t2)/(length(ageSlopeFrame2$age)*h) 
      smoothZ[i,j] <- mean(ageSlopeFrame2$fitted[t1x < 1 & t1y < 1])
   }
}
smoothP <- 100 * (smoothP/sum(smoothP))

smoothSlopeFrame <- data.frame(
        age=rep(smoothX,length(smoothY)),
        score=rep(smoothY,each=length(smoothX)),
        density=c(smoothP),
        slope=c(smoothZ))



#####################
#   MARKER F
#####################

# ------------------------------------------------------
#  Plot a slope field of the fitted slopes

pdf("VF1ARSlope1.pdf", height=5, width=6)
plot(c(0,220)/3, c(0,120), type="n" ,
     xlab="Age",
     ylab="Total Score")
scale <- 2.5
for (i in 1:length(smoothSlopeFrame$age)) {
   tx <- smoothSlopeFrame$age[i]
   ty <- smoothSlopeFrame$score[i]
   ta <- smoothSlopeFrame$slope[i]
   td <- smoothSlopeFrame$density[i]
   txinc <- (scale * td) / sqrt( 1 + (ta^2))
   tyinc <- ta * txinc
   lines(c(tx+txinc, tx-txinc), c(ty+tyinc, ty-tyinc))
}
dev.off()


# ------------------------------------------------------
#  Plot a vector field of the fitted slopes

pdf("VF2ARVector1.pdf", height=5, width=6)
plot(c(0,220)/3, c(0,101), type="n",
     xlab="Age",
     ylab="Total Score")
scale <- 2.5
for (i in 1:length(smoothSlopeFrame$age)) {
   tx <- smoothSlopeFrame$age[i]
   ty <- smoothSlopeFrame$score[i]
   ta <- smoothSlopeFrame$slope[i]
   td <- smoothSlopeFrame$density[i]
   txinc <- (scale * td) / sqrt( 1 + (ta^2))
   tyinc <- ta * txinc
   tarrowlen <- sqrt(txinc^2 + tyinc^2) * .05
   if (is.na(tarrowlen) | tarrowlen < .001) next
   arrows(tx-txinc, ty-tyinc, tx+txinc, ty+tyinc, length = tarrowlen)
}
dev.off()


# ------------------------------------------------------
#  Plot two wave data line segments including groups as colors

pdf("VF2ARLagged1b.pdf", height=5, width=6)
plot(c(0,220)/3, c(0,120), type="n" ,
     xlab="Age",
     ylab="Total Score")
for (i in seq(1,length(ageData1$x1),by=2)) {
   if (is.na(ageData1$age[i]) | 
       ageData1$age[i] > (78-tau)) {
          next
   }
   lines(c(ageData1$age[i],ageData1$age[i]+tau), 
         c(ageData1$x1[i],ageData1$x2[i]), type="l", lty=1, 
         col=ageData1$group[i])
}
dev.off()

# ------------------------------------------------------
#  Plot two wave data line segments including groups as colors

pdf("VF2ARLong2b.pdf", height=5, width=6)
plot(c(0,220)/3, c(0,120), type="n" ,
     xlab="Age",
     ylab="Total Score")
for (i in seq(1,length(ageData1$x1),by=2)) {
   if (is.na(ageData1$age[i]) | 
       ageData1$age[i] > (78-(2*tau))) {
          next
   }
   lines(c(ageData1$age[i],
           ageData1$age[i]+tau,
           ageData1$age[i]+(2*tau)), 
         c(ageData1$x1[i],
           ageData1$x2[i],
           ageData1$x3[i]), type="l", lty=1, 
         col=ageData1$group[i])
}
dev.off()


# ------------------------------------------------------
#  Calculate a loess smoothed slope surfaces by Group

tSel <- ageSlopeFrame$group==1
tFrame <- data.frame(age=ageSlopeFrame$age[tSel],
                     score=ageSlopeFrame$score[tSel],
                     slope=ageSlopeFrame$slope[tSel])
tSlopeLoess <- loess(slope ~ age * score, data=tFrame)

ageSlopeFrame2G1 <- data.frame(age=tFrame$age,
                   score=tFrame$score,
                   slope=tFrame$slope,
                   fitted=predict(tSlopeLoess))

tSel <- ageSlopeFrame$group==2
tFrame <- data.frame(age=ageSlopeFrame$age[tSel],
                     score=ageSlopeFrame$score[tSel],
                     slope=ageSlopeFrame$slope[tSel])
tSlopeLoess <- loess(slope ~ age * score, data=tFrame)

ageSlopeFrame2G2 <- data.frame(age=tFrame$age,
                   score=tFrame$score,
                   slope=tFrame$slope,
                   fitted=predict(tSlopeLoess))

tSel <- ageSlopeFrame$group==4
tFrame <- data.frame(age=ageSlopeFrame$age[tSel],
                     score=ageSlopeFrame$score[tSel],
                     slope=ageSlopeFrame$slope[tSel])
tSlopeLoess <- loess(slope ~ age * score, data=tFrame)

ageSlopeFrame2G4 <- data.frame(age=tFrame$age,
                   score=tFrame$score,
                   slope=tFrame$slope,
                   fitted=predict(tSlopeLoess))



# ------------------------------------------------------
#  Calculate a slope field of the fitted slopes

smoothX <- seq(0, 75, by=6)
smoothY <- seq(0, 120, by=8)
smoothP <- matrix(NA, length(smoothX), length(smoothY))
smoothZ <- matrix(NA, length(smoothX), length(smoothY))
h <- 10
for (i in 1:length(smoothX)) {
   x <- smoothX[i]
   t1x <- abs(x - ageSlopeFrame2G1$age)/h
   for (j in 1:length(smoothY)) {
      y <- smoothY[j]
      t1y <- abs(y - ageSlopeFrame2G1$score)/h
      t2 <- rep(0, length(t1x))
      t2[t1x < 1 & t1y < 1] <- 1/2
      smoothP[i,j] <- sum(t2)/(length(ageSlopeFrame2G1$age)*h) 
      smoothZ[i,j] <- mean(ageSlopeFrame2G1$fitted[t1x < 1 & t1y < 1])
   }
}
smoothP <- 100 * (smoothP/sum(smoothP))

smoothSlopeFrameG1 <- data.frame(
        age=rep(smoothX,length(smoothY)),
        score=rep(smoothY,each=length(smoothX)),
        density=c(smoothP),
        slope=c(smoothZ))


smoothX <- seq(0, 75, by=6)
smoothY <- seq(0, 120, by=8)
smoothP <- matrix(NA, length(smoothX), length(smoothY))
smoothZ <- matrix(NA, length(smoothX), length(smoothY))
h <- 10
for (i in 1:length(smoothX)) {
   x <- smoothX[i]
   t1x <- abs(x - ageSlopeFrame2G2$age)/h
   for (j in 1:length(smoothY)) {
      y <- smoothY[j]
      t1y <- abs(y - ageSlopeFrame2G2$score)/h
      t2 <- rep(0, length(t1x))
      t2[t1x < 1 & t1y < 1] <- 1/2
      smoothP[i,j] <- sum(t2)/(length(ageSlopeFrame2G2$age)*h) 
      smoothZ[i,j] <- mean(ageSlopeFrame2G2$fitted[t1x < 1 & t1y < 1])
   }
}
smoothP <- 100 * (smoothP/sum(smoothP))

smoothSlopeFrameG2 <- data.frame(
        age=rep(smoothX,length(smoothY)),
        score=rep(smoothY,each=length(smoothX)),
        density=c(smoothP),
        slope=c(smoothZ))


smoothX <- seq(0, 75, by=6)
smoothY <- seq(0, 120, by=8)
smoothP <- matrix(NA, length(smoothX), length(smoothY))
smoothZ <- matrix(NA, length(smoothX), length(smoothY))
h <- 10
for (i in 1:length(smoothX)) {
   x <- smoothX[i]
   t1x <- abs(x - ageSlopeFrame2G4$age)/h
   for (j in 1:length(smoothY)) {
      y <- smoothY[j]
      t1y <- abs(y - ageSlopeFrame2G4$score)/h
      t2 <- rep(0, length(t1x))
      t2[t1x < 1 & t1y < 1] <- 1/2
      smoothP[i,j] <- sum(t2)/(length(ageSlopeFrame2G4$age)*h) 
      smoothZ[i,j] <- mean(ageSlopeFrame2G4$fitted[t1x < 1 & t1y < 1])
   }
}
smoothP <- 100 * (smoothP/sum(smoothP))

smoothSlopeFrameG4 <- data.frame(
        age=rep(smoothX,length(smoothY)),
        score=rep(smoothY,each=length(smoothX)),
        density=c(smoothP),
        slope=c(smoothZ))

# ------------------------------------------------------
#  Plot a slope field of the fitted slopes grouped by color

pdf("VF1ARSlope1Groups.pdf", height=5, width=6)
plot(c(0,220)/3, c(0,120), type="n" ,
     xlab="Age",
     ylab="Total Score")
scale <- 3
for (i in 1:length(smoothSlopeFrameG1$age)) {
   tx <- smoothSlopeFrameG1$age[i]
   ty <- smoothSlopeFrameG1$score[i]
   ta <- smoothSlopeFrameG1$slope[i]
   td <- smoothSlopeFrameG1$density[i]
   txinc <- (scale * td) / sqrt( 1 + (ta^2))
   tyinc <- ta * txinc
   lines(c(tx+txinc, tx-txinc), c(ty+tyinc, ty-tyinc), col=1)
}
for (i in 1:length(smoothSlopeFrameG2$age)) {
   tx <- smoothSlopeFrameG2$age[i]
   ty <- smoothSlopeFrameG2$score[i]
   ta <- smoothSlopeFrameG2$slope[i]
   td <- smoothSlopeFrameG2$density[i]
   txinc <- (scale * td) / sqrt( 1 + (ta^2))
   tyinc <- ta * txinc
   lines(c(tx+txinc, tx-txinc), c(ty+tyinc, ty-tyinc), col=2)
}
for (i in 1:length(smoothSlopeFrameG4$age)) {
   tx <- smoothSlopeFrameG4$age[i]
   ty <- smoothSlopeFrameG4$score[i]
   ta <- smoothSlopeFrameG4$slope[i]
   td <- smoothSlopeFrameG4$density[i]
   txinc <- (scale * td) / sqrt( 1 + (ta^2))
   tyinc <- ta * txinc
   lines(c(tx+txinc, tx-txinc), c(ty+tyinc, ty-tyinc), col=4)
}
dev.off()

# ------------------------------------------------------
#  Quit the program here

