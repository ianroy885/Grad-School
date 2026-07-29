library("mirt")

## mirt() example
data <- expand.table(LSAT7)
(mod1 <- mirt(data, 1))
coef(mod1)
summary(mod1)
residuals(mod1)
residuals(mod1, restype = "exp")

(mod2 <- mirt(data, 2))
summary(mod2, rotate = "oblimin", suppress = 0.25)
anova(mod1, mod2)

## bfactor() example
data("SAT12")
data <- key2binary(SAT12, key = c(1, 4, 5, 2, 3, 1, 2, 1, 3,
	1, 2, 4, 2, 1, 5, 3, 4, 4, 1, 4, 3, 3, 4, 1, 3, 5, 1, 3,
	1, 5, 4, 5))
specific <- c(2, 3, 2, 3, 3, 2, 1, 2, 1, 1, 1, 3, 1, 3, 1, 2,
	1, 1, 3, 3, 1, 1, 3, 1, 3, 3, 1, 3, 2, 3, 1, 2)
guess <- rep(0.1, 32)
b_mod1 <- bfactor(data, specific, guess)
coef(b_mod1)

b_mod2 <- bfactor(data, specific, guess, par.prior = list(int = c(0,
	4), int.items = 32))
coef(b_mod2)

## polymirt() example
#simulated data parameters from ?simdata
a <- matrix(c(
 .7471, .0250, .1428,
 .4595, .0097, .0692,
 .8613, .0067, .4040,
1.0141, .0080, .0470,
 .5521, .0204, .1482,
1.3547, .0064, .5362,
1.3761, .0861, .4676,
 .8525, .0383, .2574,
1.0113, .0055, .2024,
 .9212, .0119, .3044,
 .0026, .0119, .8036,
 .0008, .1905,1.1945,
 .0575, .0853, .7077,
 .0182, .3307,2.1414,
 .0256, .0478, .8551,
 .0246, .1496, .9348,
 .0262, .2872,1.3561,
 .0038, .2229, .8993,
 .0039, .4720, .7318,
 .0068, .0949, .6416,
 .3073, .9704, .0031,
 .1819, .4980, .0020,
 .4115,1.1136, .2008,
 .1536,1.7251, .0345,
 .1530, .6688, .0020,
 .2890,1.2419, .0220,
 .1341,1.4882, .0050,
 .0524, .4754, .0012,
 .2139, .4612, .0063,
 .1761,1.1200, .0870),30,3,byrow=TRUE)

d <- matrix(c(.1826,-.1924,-.4656,-.4336,-.4428,-.5845,-1.0403,
  .6431,.0122,.0912,.8082,-.1867,.4533,-1.8398,.4139,
  -.3004,-.1824,.5125,1.1342,.0230,.6172,-.1955,-.3668,
  -1.7590,-.2434,.4925,-.3410,.2896,.006,.0329),ncol=1)

simdata1 <- simdata(a, d, 2000)

p_mod <- polymirt(simdata1, 3)
p_mod
coef(p_mod)

## confmirt() example
#data from Appendix B
a <- matrix(c(1.5, NA, 0.5, NA, 1, NA, 1, 0.5, NA, 1.5, NA, 0.5,
	NA, 1, NA, 1), ncol = 2, byrow = TRUE)
d <- matrix(c(-1, NA, NA, -1.5, NA, NA, 1.5, NA, NA, 0, NA, NA,
	3, 2, -0.5, 2.5, 1, -1, 2, 0, NA, 1, NA, NA), ncol = 3, byrow = TRUE)
sigma <- diag(2)
sigma[1, 2] <- sigma[2, 1] <- 0.4
simdata2 <- simdata(a, d, 2000, sigma)

model.1 <- confmirt.model()
	F1 = 1-4
	F2 = 4-8
	COV = F1*F2
	
c_mod <- confmirt(simdata2, model.1, printcycles = FALSE)
c_mod
coef(c_mod)

model.2 <- confmirt.model()
	F1 = 1-4
	F2 = 4-8
	COV = F1*F2
	SLOPE = F1@1 eq 1.346, F1@3 eq F1@4, F2@7 eq F2@8

c_mod2 <- confmirt(simdata2, model.2, printcycles = FALSE)
coef(c_mod2)
anova(c_mod2, c_mod)

## Code from Appendix A
plot(mod1)
plot(mod2)
itemplot(mod1, combine = 5, auto.key=list(space="right"))
fscores(mod1)
dataWithFS <- fscores(mod1, full.scores = TRUE)

## Code from Appendix B
a <- matrix(c(1.5, NA, 0.5, NA, 1, NA, 1, 0.5, NA, 1.5, NA, 0.5,
	NA, 1, NA, 1), ncol = 2, byrow = TRUE)
d <- matrix(c(-1, NA, NA, -1.5, NA, NA, 1.5, NA, NA, 0, NA, NA,
	3, 2, -0.5, 2.5, 1, -1, 2, 0, NA, 1, NA, NA), ncol = 3, byrow = TRUE)
sigma <- diag(2)
sigma[1, 2] <- sigma[2, 1] <- 0.4
simdata2 <- simdata(a, d, 2000, sigma)
