#
# BAM_ANES2008.r
#
library(foreign)
library(rjags)
library(coda)
library(basicspace)
library(lattice)
#
#
ANES2008 <- read.dta("c:/Dropbox/BAM/ANES/ANES_2008/ANES2008.dta", convert.factors=FALSE)
attach(ANES2008)
#
#     
libcon.self <- V083069
libcon.self[libcon.self < 0] <- NA
#
libcon.Obama <- V083070a
libcon.Obama[libcon.Obama < 0] <- NA
#
libcon.McCain <- V083070b
libcon.McCain[libcon.McCain < 0] <- NA
#
libcon.DemParty <- V083071a
libcon.DemParty[libcon.DemParty < 0] <- NA
#
libcon.RepParty <- V083071b
libcon.RepParty[libcon.RepParty < 0] <- NA
#
#
T <- cbind(libcon.self, libcon.Obama, libcon.McCain, libcon.DemParty, libcon.RepParty)
#
# This Flips the Polarity of a Seven-Point Scale (put conservative
# responses on the right):
#
#T <- abs(T - 8)
#
T <- T - 4
#
libcon.self <- T[,1]
TT <- T[,-1]
#
cutoff <- 3
libcon_self <- libcon.self[rowSums(!is.na(TT)) >= cutoff]
TTT <- TT[rowSums(!is.na(TT)) >= cutoff,]
#
#
N <- nrow(TTT)
q <- ncol(TTT)
z <- TTT
#
alpha.starts <- runif(N,-2,2)
beta.starts <- runif(N,-2,2)
#
zhat.starts <- runif(q,-1,1)
zhat.starts[1] <- -1 + runif(1,-0.099,0.099)
zhat.starts[2] <- 1 + runif(1,-0.099,0.099)
#
inits <- function() {list (zhat=zhat.starts, a=alpha.starts, b=beta.starts)}
#
BAM.sim <- jags.model('c:/Dropbox/BAM/ANES/ANES_2008/BAM_JAGScode_ANES2008.bug',
    data = list('z' = z, 'q' = q, 'N' = N),
    inits = inits, n.chains = 2, n.adapt = 5000)
#
update(BAM.sim, n.iter=45000)
#
params <- c("zhat", "a", "b", "tauj", "sigmaj.sq")
samps <- coda.samples(BAM.sim, params, n.iter=10000, thin=10)
#
#
# Write results to disk:
save(samps,file="c:/Dropbox/BAM/ANES/ANES_2008/ANES2008_samps_libcon.Rda")
#
#
#   MCMC DIAGNOSTICS
#
# Summary:
summary(samps[,"a[1]"])
summary(samps[,"b[1]"])
summary(samps[,"zhat[1]"])
#
# Trace and Density Plots:
plot(samps[,"a[1]"])
plot(samps[,"b[1]"])
plot(samps[,"zhat[1]"])
#
# Autocorrelation Plots:
autocorr.plot(samps[,"a[1]"])
autocorr.plot(samps[,"b[1]"])
autocorr.plot(samps[,"zhat[1]"])
#
# Gelman-Rubin Diagnostic (need multiple chains):
gelman.diag(samps[,"a[1]"])
gelman.diag(samps[,"b[1]"])
gelman.diag(samps[,"zhat[1]"])
#
# Geweke Diagnostic:
geweke.diag(samps[,"a[1]"])
geweke.diag(samps[,"b[1]"])
geweke.diag(samps[,"zhat[1]"])
#