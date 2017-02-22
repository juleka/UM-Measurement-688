#
# BAM_ANES_libcon_scales.r
#
library(foreign)
library(rjags)
library(coda)
library(basicspace)
library(lattice)
library(latticeExtra)
library(overlap)
library(apsrtable)
library(lmtest)
library(DAMisc)
#
#
# %%%%%%%%%%%%%%%   2012   %%%%%%%%%%%%%%%
#
ANES2012 <- read.dta("c:/Dropbox/BAM/ANES/ANES_2012/ANES2012.dta", convert.factors=FALSE)
attach(ANES2012)
#
# =============================================================
#       LIBERAL-CONSERVATIVE SCALE RESULTS:
#
libcon.self.2012 <- libcpre_self
libcon.self.2012[libcon.self.2012 < 0] <- NA
#
libcon.Obama.2012 <- libcpre_dpc
libcon.Obama.2012[libcon.Obama.2012 < 0] <- NA
#
libcon.Romney.2012 <- libcpre_rpc
libcon.Romney.2012[libcon.Romney.2012 < 0] <- NA
#
libcon.DemParty.2012 <- libcpre_ptyd
libcon.DemParty.2012[libcon.DemParty.2012 < 0] <- NA
#
libcon.RepParty.2012 <- libcpre_ptyr
libcon.RepParty.2012[libcon.RepParty.2012 < 0] <- NA
#
#
libcon_T.2012 <- cbind(libcon.self.2012, libcon.Obama.2012, libcon.Romney.2012, libcon.DemParty.2012, libcon.RepParty.2012)
#
libcon_T.2012 <- libcon_T.2012 - 4
#
libcon.self.2012 <- libcon_T.2012[,1]
libcon_TT.2012 <- libcon_T.2012[,-1]
#
cutoff <- 3
libcon_self.2012 <- libcon.self.2012[rowSums(!is.na(libcon_TT.2012)) >= cutoff]
libcon_TTT.2012 <- libcon_TT.2012[rowSums(!is.na(libcon_TT.2012)) >= cutoff,]
#
# PARTY ID: 1 = STRONG DEM, 7 = STRONG REP
libcon_party.2012 <- pid_x[rowSums(!is.na(libcon_TT.2012)) >= cutoff]
libcon_party.2012[libcon_party.2012 < 0] <- NA
#
libcon_three_party.2012 <- libcon_party.2012
libcon_three_party.2012[libcon_three_party.2012==1 | libcon_three_party.2012==2 | libcon_three_party.2012==3] <- 1
libcon_three_party.2012[libcon_three_party.2012==4] <- 2
libcon_three_party.2012[libcon_three_party.2012==5 | libcon_three_party.2012==6 | libcon_three_party.2012==7] <- 3
#
#
# PRESIDENTIAL VOTE: 0 = OBAMA, 1 = ROMNEY
libcon_presvote.2012 <- presvote2012_x[rowSums(!is.na(libcon_TT.2012)) >= cutoff]
libcon_presvote.2012[libcon_presvote.2012==5] <- NA
libcon_presvote.2012[libcon_presvote.2012 < 0] <- NA
libcon_presvote.2012 <- libcon_presvote.2012 - 1 
#
load("c:/Dropbox/BAM/ANES/ANES_2012/ANES2012_samps_libcon.Rda")
libcon_samps.2012 <- samps
#
libcon_zhat.2012 <- summary(libcon_samps.2012[,grep("zhat", colnames(libcon_samps.2012[[1]]))])
libcon_sigmaj.sq.2012 <- summary(libcon_samps.2012[,grep("sigmaj.sq", colnames(libcon_samps.2012[[1]]))])
#
libcon_a.2012 <- libcon_samps.2012[,grep(pattern="a[", colnames(libcon_samps.2012[[1]]), fixed=TRUE),]
libcon_b.2012 <- libcon_samps.2012[,grep(pattern="b[", colnames(libcon_samps.2012[[1]]), fixed=TRUE),]
libcon_a_point_est.2012 <- summary(libcon_a.2012)[[1]][,1]
libcon_b_point_est.2012 <- summary(libcon_b.2012)[[1]][,1]
libcon_all.a.2012 <- do.call("rbind", libcon_a.2012)
libcon_all.b.2012 <- do.call("rbind", libcon_b.2012)
#
libcon_nsamp.2012 <- nrow(libcon_all.a.2012)
libcon_nresp.2012 <- ncol(libcon_all.a.2012)
#
libcon_idealpt.2012 <- rep(NA, libcon_nsamp.2012*libcon_nresp.2012)
dim(libcon_idealpt.2012) <- c(libcon_nsamp.2012,libcon_nresp.2012)
for (i in 1:libcon_nresp.2012){
for (j in 1:libcon_nsamp.2012){
libcon_idealpt.2012[j,i] <- ((libcon_self.2012[i] - libcon_all.a.2012[j,i]) / libcon_all.b.2012[j,i])
}}
#
libcon_idealpt.percentiles.2012 <- rep(NA,libcon_nresp.2012*3)
dim(libcon_idealpt.percentiles.2012) <- c(libcon_nresp.2012,3)
colnames(libcon_idealpt.percentiles.2012) <- c("5%","50%","95%")
for (i in 1:libcon_nresp.2012){
libcon_idealpt.percentiles.2012[i,] <- quantile(libcon_idealpt.2012[,i],probs=c(0.05,0.5,0.95),na.rm=TRUE)
}
libcon_idealpt.estimates.2012 <- libcon_idealpt.percentiles.2012[,2]
#
#
#
# Probit Model of Vote Choice for All Voters
res1.2012 <- glm(formula = libcon_presvote.2012 ~ libcon_self.2012 + libcon_party.2012, family = binomial(link = "probit"))
res2.2012 <- glm(formula = libcon_presvote.2012 ~ libcon_idealpt.estimates.2012 + libcon_party.2012, family = binomial(link = "probit"))
apsrtable(res1.2012,res2.2012)
#
# Probit Model of Vote Choice for Self-Identified Moderates and Slightly Liberal/Conservative
res1.mods.2012 <- glm(formula = libcon_presvote.2012[libcon_self.2012 >= -1 & libcon_self.2012 <= 1] ~ libcon_self.2012[libcon_self.2012 >= -1 & libcon_self.2012 <= 1] + libcon_party.2012[libcon_self.2012 >= -1 & libcon_self.2012 <= 1], family = binomial(link = "probit"))
res2.mods.2012 <- glm(formula = libcon_presvote.2012[libcon_self.2012 >= -1 & libcon_self.2012 <= 1] ~ libcon_idealpt.estimates.2012[libcon_self.2012 >= -1 & libcon_self.2012 <= 1] + libcon_party.2012[libcon_self.2012 >= -1 & libcon_self.2012 <= 1], family = binomial(link = "probit"))
apsrtable(res1.mods.2012,res2.mods.2012)
#
# Probit Model of Vote Choice for Only Self-Identified Moderates 
res1.modsonly.2012 <- glm(formula = libcon_presvote.2012[libcon_self.2012 == 0] ~ libcon_party.2012[libcon_self.2012 == 0], family = binomial(link = "probit"))
res2.modsonly.2012 <- glm(formula = libcon_presvote.2012[libcon_self.2012 == 0] ~ libcon_idealpt.estimates.2012[libcon_self.2012 == 0] + libcon_party.2012[libcon_self.2012 == 0], family = binomial(link = "probit"))
apsrtable(res1.modsonly.2012,res2.modsonly.2012)
#
#
# Overlap Estimates for Presidential Vote
Obama_idealpt.estimates.2012 <- na.omit(libcon_idealpt.estimates.2012[libcon_presvote.2012==0])
Romney_idealpt.estimates.2012 <- na.omit(libcon_idealpt.estimates.2012[libcon_presvote.2012==1])
#
Obama_libcon_self.2012 <- na.omit(libcon_self.2012[libcon_presvote.2012==0])
Romney_libcon_self.2012 <- na.omit(libcon_self.2012[libcon_presvote.2012==1])
#
#  Rescale between 0 and 2*pi
range02pi <- function(x){((x-min(x))/(max(x)-min(x)))*(2*pi)}
tmp.2012.correct <- range02pi(c(Obama_idealpt.estimates.2012,Romney_idealpt.estimates.2012))
tmp.2012.raw <- range02pi(c(Obama_libcon_self.2012,Romney_libcon_self.2012))
Obama_idealpt.estimates.2012.rescaled <- tmp.2012.correct[1:length(Obama_idealpt.estimates.2012)]
Romney_idealpt.estimates.2012.rescaled <- tmp.2012.correct[(length(Obama_idealpt.estimates.2012)+1):length(tmp.2012.correct)]
Obama_libcon_self.2012.rescaled <- tmp.2012.raw[1:length(Obama_libcon_self.2012)]
Romney_libcon_self.2012.rescaled <- tmp.2012.raw[(length(Obama_libcon_self.2012)+1):length(tmp.2012.raw)]
#
#  Bootstrapped standard errors for Overlap Estimates
Ntrials <- 1000
overlap.2012_idealpt.estimates <- vector("numeric", Ntrials)
for(i in 1:Ntrials){
overlap.2012_idealpt.estimates[[i]] <- 
	overlapEst(Obama_idealpt.estimates.2012.rescaled[sample(1:length(Obama_idealpt.estimates.2012.rescaled),length(Obama_idealpt.estimates.2012.rescaled),replace=TRUE)],
	Romney_idealpt.estimates.2012.rescaled[sample(1:length(Romney_idealpt.estimates.2012.rescaled),length(Romney_idealpt.estimates.2012.rescaled),replace=TRUE)])[2]
}
point.overlap.2012_idealpt.estimates <- overlapEst(Obama_idealpt.estimates.2012.rescaled,Romney_idealpt.estimates.2012.rescaled)[2]
summary.overlap.2012_idealpt.estimates <- c(point.overlap.2012_idealpt.estimates, 
	sd(overlap.2012_idealpt.estimates),
	(point.overlap.2012_idealpt.estimates-sd(overlap.2012_idealpt.estimates*1.96)),
	(point.overlap.2012_idealpt.estimates+sd(overlap.2012_idealpt.estimates*1.96)))
names(summary.overlap.2012_idealpt.estimates) <- c("point","sd","lower","upper")
#
#
overlap.2012_libcon.self <- vector("numeric", Ntrials)
for(i in 1:Ntrials){
overlap.2012_libcon.self[[i]] <- 
	overlapEst(Obama_libcon_self.2012.rescaled[sample(1:length(Obama_libcon_self.2012.rescaled),length(Obama_libcon_self.2012.rescaled),replace=TRUE)],
	Romney_libcon_self.2012.rescaled[sample(1:length(Romney_libcon_self.2012.rescaled),length(Romney_libcon_self.2012.rescaled),replace=TRUE)])[2]
}
point.overlap.2012_libcon.self <- overlapEst(Obama_libcon_self.2012.rescaled,Romney_libcon_self.2012.rescaled)[2]
summary.overlap.2012_libcon.self <- c(point.overlap.2012_libcon.self, 
	sd(overlap.2012_libcon.self),
	(point.overlap.2012_libcon.self-sd(overlap.2012_libcon.self)*1.96),
	(point.overlap.2012_libcon.self+sd(overlap.2012_libcon.self)*1.96))
names(summary.overlap.2012_libcon.self) <- c("point","sd","lower","upper")
#
print(summary.overlap.2012_idealpt.estimates)
print(summary.overlap.2012_libcon.self)
#
#  USE UNCERTAINTY IN IDEAL POINT ESTIMATES TO DEVELOP
#  BAYESIAN CREDIBLE INTERVALS FOR THE OVERLAP POINT ESTIMATES
  tmpjunk <- matrix(NA, nrow=Ntrials, ncol=libcon_nresp.2012)
	for (i in 1:Ntrials){
	for (j in 1:libcon_nresp.2012){
	tmpjunk[i,j] <- sample(libcon_idealpt.2012[,j], 1)
	}}
  #
  TEST.Obama_libcon_self.2012 <- list(NA)
  for (i in 1:Ntrials){
  TEST.Obama_libcon_self.2012[[i]] <- na.omit(tmpjunk[i,libcon_presvote.2012==0])
  }
  TEST.Romney_libcon_self.2012 <- list(NA)
  for (i in 1:Ntrials){
  TEST.Romney_libcon_self.2012[[i]] <- na.omit(tmpjunk[i,libcon_presvote.2012==1])
  }
  for (i in 1:Ntrials){
  TEST.Obama_libcon_self.2012[[i]] <- TEST.Obama_libcon_self.2012[[i]][TEST.Obama_libcon_self.2012[[i]] < 10 & TEST.Obama_libcon_self.2012[[i]] > -10]
  TEST.Romney_libcon_self.2012[[i]] <- TEST.Romney_libcon_self.2012[[i]][TEST.Romney_libcon_self.2012[[i]] < 10 & TEST.Romney_libcon_self.2012[[i]] > -10]
  }
  #
  TEST.tmp.2012.correct <- list(NA)
  for (i in 1:Ntrials){
  TEST.tmp.2012.correct[[i]] <- range02pi(c(TEST.Obama_libcon_self.2012[[i]],TEST.Romney_libcon_self.2012[[i]]))
  }
  TEST.Obama_idealpt.estimates.2012.rescaled <- list(NA)
  TEST.Romney_idealpt.estimates.2012.rescaled <- list(NA)
  for (i in 1:Ntrials){
  TEST.Obama_idealpt.estimates.2012.rescaled[[i]] <- TEST.tmp.2012.correct[[i]][1:length(TEST.Obama_libcon_self.2012[[i]])]
  TEST.Romney_idealpt.estimates.2012.rescaled[[i]] <- TEST.tmp.2012.correct[[i]][(length(TEST.Obama_libcon_self.2012[[i]])+1):length(TEST.tmp.2012.correct[[i]])]
  }
  #
  TEST.overlap.2012_idealpt.estimates <- vector("numeric", Ntrials)
  for(i in 1:Ntrials){
  TEST.overlap.2012_idealpt.estimates[[i]] <- 
	overlapEst(TEST.Obama_idealpt.estimates.2012.rescaled[[i]],TEST.Romney_idealpt.estimates.2012.rescaled[[i]])[2]
  }
  #
  TEST.summary.overlap.2012_idealpt.estimates <- quantile(TEST.overlap.2012_idealpt.estimates, probs = c(0.025,0.5,0.975))
  #
#
#
#
vote.summary.overlap.2012_idealpt.estimates <- c(mean(TEST.summary.overlap.2012_idealpt.estimates), sd(TEST.summary.overlap.2012_idealpt.estimates), TEST.summary.overlap.2012_idealpt.estimates[1], TEST.summary.overlap.2012_idealpt.estimates[3])
names(vote.summary.overlap.2012_idealpt.estimates) <- c("point","sd","lower","upper")
vote.summary.overlap.2012_libcon.self <- summary.overlap.2012_libcon.self
#
# =============================================================
#
# Overlap Estimates for Party Identification
Democrats_idealpt.estimates.2012 <- na.omit(libcon_idealpt.estimates.2012[libcon_three_party.2012==1])
Republicans_idealpt.estimates.2012 <- na.omit(libcon_idealpt.estimates.2012[libcon_three_party.2012==3])
#
#
#  Rescale between 0 and 2*pi
range02pi <- function(x){((x-min(x))/(max(x)-min(x)))*(2*pi)}
tmp.2012.correct <- range02pi(c(Democrats_idealpt.estimates.2012,Republicans_idealpt.estimates.2012))
tmp.2012.raw <- range02pi(c(Democrats_libcon_self.2012,Republicans_libcon_self.2012))
Democrats_idealpt.estimates.2012.rescaled <- tmp.2012.correct[1:length(Democrats_idealpt.estimates.2012)]
Republicans_idealpt.estimates.2012.rescaled <- tmp.2012.correct[(length(Democrats_idealpt.estimates.2012)+1):length(tmp.2012.correct)]
Democrats_libcon_self.2012.rescaled <- tmp.2012.raw[1:length(Democrats_libcon_self.2012)]
Republicans_libcon_self.2012.rescaled <- tmp.2012.raw[(length(Democrats_libcon_self.2012)+1):length(tmp.2012.raw)]
#
#  Bootstrapped standard errors for Overlap Estimates
Ntrials <- 1000
overlap.2012_idealpt.estimates <- vector("numeric", Ntrials)
for(i in 1:Ntrials){
overlap.2012_idealpt.estimates[[i]] <- 
	overlapEst(Democrats_idealpt.estimates.2012.rescaled[sample(1:length(Democrats_idealpt.estimates.2012.rescaled),length(Democrats_idealpt.estimates.2012.rescaled),replace=TRUE)],
	Republicans_idealpt.estimates.2012.rescaled[sample(1:length(Republicans_idealpt.estimates.2012.rescaled),length(Republicans_idealpt.estimates.2012.rescaled),replace=TRUE)])[2]
}
point.overlap.2012_idealpt.estimates <- overlapEst(Democrats_idealpt.estimates.2012.rescaled,Republicans_idealpt.estimates.2012.rescaled)[2]
summary.overlap.2012_idealpt.estimates <- c(point.overlap.2012_idealpt.estimates, 
	sd(overlap.2012_idealpt.estimates),
	(point.overlap.2012_idealpt.estimates-sd(overlap.2012_idealpt.estimates*1.96)),
	(point.overlap.2012_idealpt.estimates+sd(overlap.2012_idealpt.estimates*1.96)))
names(summary.overlap.2012_idealpt.estimates) <- c("point","sd","lower","upper")
#
#
overlap.2012_libcon.self <- vector("numeric", Ntrials)
for(i in 1:Ntrials){
overlap.2012_libcon.self[[i]] <- 
	overlapEst(Democrats_libcon_self.2012.rescaled[sample(1:length(Democrats_libcon_self.2012.rescaled),length(Democrats_libcon_self.2012.rescaled),replace=TRUE)],
	Republicans_libcon_self.2012.rescaled[sample(1:length(Republicans_libcon_self.2012.rescaled),length(Republicans_libcon_self.2012.rescaled),replace=TRUE)])[2]
}
point.overlap.2012_libcon.self <- overlapEst(Democrats_libcon_self.2012.rescaled,Republicans_libcon_self.2012.rescaled)[2]
summary.overlap.2012_libcon.self <- c(point.overlap.2012_libcon.self, 
	sd(overlap.2012_libcon.self),
	(point.overlap.2012_libcon.self-sd(overlap.2012_libcon.self)*1.96),
	(point.overlap.2012_libcon.self+sd(overlap.2012_libcon.self)*1.96))
names(summary.overlap.2012_libcon.self) <- c("point","sd","lower","upper")
#
print(summary.overlap.2012_idealpt.estimates)
print(summary.overlap.2012_libcon.self)
#
#  USE UNCERTAINTY IN IDEAL POINT ESTIMATES TO DEVELOP
#  BAYESIAN CREDIBLE INTERVALS FOR THE OVERLAP POINT ESTIMATES
  tmpjunk <- matrix(NA, nrow=Ntrials, ncol=libcon_nresp.2012)
	for (i in 1:Ntrials){
	for (j in 1:libcon_nresp.2012){
	tmpjunk[i,j] <- sample(libcon_idealpt.2012[,j], 1)
	}}
  #
  TEST.Democrats_libcon_self.2012 <- list(NA)
  for (i in 1:Ntrials){
  TEST.Democrats_libcon_self.2012[[i]] <- na.omit(tmpjunk[i,libcon_three_party.2012==1])
  }
  TEST.Republicans_libcon_self.2012 <- list(NA)
  for (i in 1:Ntrials){
  TEST.Republicans_libcon_self.2012[[i]] <- na.omit(tmpjunk[i,libcon_three_party.2012==3])
  }
  for (i in 1:Ntrials){
  TEST.Democrats_libcon_self.2012[[i]] <- TEST.Democrats_libcon_self.2012[[i]][TEST.Democrats_libcon_self.2012[[i]] < 10 & TEST.Democrats_libcon_self.2012[[i]] > -10]
  TEST.Republicans_libcon_self.2012[[i]] <- TEST.Republicans_libcon_self.2012[[i]][TEST.Republicans_libcon_self.2012[[i]] < 10 & TEST.Republicans_libcon_self.2012[[i]] > -10]
  }
  #
  TEST.tmp.2012.correct <- list(NA)
  for (i in 1:Ntrials){
  TEST.tmp.2012.correct[[i]] <- range02pi(c(TEST.Democrats_libcon_self.2012[[i]],TEST.Republicans_libcon_self.2012[[i]]))
  }
  TEST.Democrats_idealpt.estimates.2012.rescaled <- list(NA)
  TEST.Republicans_idealpt.estimates.2012.rescaled <- list(NA)
  for (i in 1:Ntrials){
  TEST.Democrats_idealpt.estimates.2012.rescaled[[i]] <- TEST.tmp.2012.correct[[i]][1:length(TEST.Democrats_libcon_self.2012[[i]])]
  TEST.Republicans_idealpt.estimates.2012.rescaled[[i]] <- TEST.tmp.2012.correct[[i]][(length(TEST.Democrats_libcon_self.2012[[i]])+1):length(TEST.tmp.2012.correct[[i]])]
  }
  #
  TEST.overlap.2012_idealpt.estimates <- vector("numeric", Ntrials)
  for(i in 1:Ntrials){
  TEST.overlap.2012_idealpt.estimates[[i]] <- 
	overlapEst(TEST.Democrats_idealpt.estimates.2012.rescaled[[i]],TEST.Republicans_idealpt.estimates.2012.rescaled[[i]])[2]
  }
  #
  TEST.summary.overlap.2012_idealpt.estimates <- quantile(TEST.overlap.2012_idealpt.estimates, probs = c(0.025,0.5,0.975))
  #
#
party.summary.overlap.2012_idealpt.estimates <- c(mean(TEST.summary.overlap.2012_idealpt.estimates), sd(TEST.summary.overlap.2012_idealpt.estimates), TEST.summary.overlap.2012_idealpt.estimates[1], TEST.summary.overlap.2012_idealpt.estimates[3])
names(party.summary.overlap.2012_idealpt.estimates) <- c("point","sd","lower","upper")
#party.summary.overlap.2012_idealpt.estimates <- summary.overlap.2012_idealpt.estimates
party.summary.overlap.2012_libcon.self <- summary.overlap.2012_libcon.self
#
#
# =============================================================
#
#
#
# %%%%%%%%%%%%%%%   2008   %%%%%%%%%%%%%%%
#
detach(ANES2012)
ANES2008 <- read.dta("c:/Dropbox/BAM/ANES/ANES_2008/ANES2008.dta", convert.factors=FALSE)
attach(ANES2008)
#
# =============================================================
#       LIBERAL-CONSERVATIVE SCALE RESULTS:
#
libcon.self.2008 <- V083069
libcon.self.2008[libcon.self.2008 < 0] <- NA
#
libcon.Obama.2008 <- V083070a
libcon.Obama.2008[libcon.Obama.2008 < 0] <- NA
#
libcon.McCain.2008 <- V083070b
libcon.McCain.2008[libcon.McCain.2008 < 0] <- NA
#
libcon.DemParty.2008 <- V083071a
libcon.DemParty.2008[libcon.DemParty.2008 < 0] <- NA
#
libcon.RepParty.2008 <- V083071b
libcon.RepParty.2008[libcon.RepParty.2008 < 0] <- NA
#
#
libcon_T.2008 <- cbind(libcon.self.2008, libcon.Obama.2008, libcon.McCain.2008, libcon.DemParty.2008, libcon.RepParty.2008)
#
libcon_T.2008 <- libcon_T.2008 - 4
#
libcon.self.2008 <- libcon_T.2008[,1]
libcon_TT.2008 <- libcon_T.2008[,-1]
#
cutoff <- 3
libcon_self.2008 <- libcon.self.2008[rowSums(!is.na(libcon_TT.2008)) >= cutoff]
libcon_TTT.2008 <- libcon_TT.2008[rowSums(!is.na(libcon_TT.2008)) >= cutoff,]
#
# PARTY ID: 1 = STRONG DEM, 7 = STRONG REP
libcon_party.2008 <- V083098x[rowSums(!is.na(libcon_TT.2008)) >= cutoff]
libcon_party.2008 <- libcon_party.2008 + 1
libcon_party.2008[libcon_party.2008==0] <- NA
#
libcon_three_party.2008 <- libcon_party.2008
libcon_three_party.2008[libcon_three_party.2008==1 | libcon_three_party.2008==2 | libcon_three_party.2008==3] <- 1
libcon_three_party.2008[libcon_three_party.2008==4] <- 2
libcon_three_party.2008[libcon_three_party.2008==5 | libcon_three_party.2008==6 | libcon_three_party.2008==7] <- 3
#
#
# PRESIDENTIAL VOTE: 0 = OBAMA, 1 = MCCAIN
libcon_presvote.2008 <- V083169a
libcon_presvote.2008[libcon_presvote.2008 < 0] <- NA
libcon_presvote.2008[libcon_presvote.2008 > 4] <- NA
libcon_presvote.2008[V085044a==1] <- 1
libcon_presvote.2008[V085044a==3] <- 2
libcon_presvote.2008 <- libcon_presvote.2008[rowSums(!is.na(libcon_TT.2008)) >= cutoff]
libcon_presvote.2008 <- libcon_presvote.2008 - 1
#
#
rm(samps)
load("c:/Dropbox/BAM/ANES/ANES_2008/ANES2008_samps_libcon.Rda")
libcon_samps.2008 <- samps
#
libcon_zhat.2008 <- summary(libcon_samps.2008[,grep("zhat", colnames(libcon_samps.2008[[1]]))])
libcon_sigmaj.sq.2008 <- summary(libcon_samps.2008[,grep("sigmaj.sq", colnames(libcon_samps.2008[[1]]))])
#
libcon_a.2008 <- libcon_samps.2008[,grep(pattern="a[", colnames(libcon_samps.2008[[1]]), fixed=TRUE),]
libcon_b.2008 <- libcon_samps.2008[,grep(pattern="b[", colnames(libcon_samps.2008[[1]]), fixed=TRUE),]
libcon_a_point_est.2008 <- summary(libcon_a.2008)[[1]][,1]
libcon_b_point_est.2008 <- summary(libcon_b.2008)[[1]][,1]
libcon_all.a.2008 <- do.call("rbind", libcon_a.2008)
libcon_all.b.2008 <- do.call("rbind", libcon_b.2008)
#
libcon_nsamp.2008 <- nrow(libcon_all.a.2008)
libcon_nresp.2008 <- ncol(libcon_all.a.2008)
#
libcon_idealpt.2008 <- rep(NA, libcon_nsamp.2008*libcon_nresp.2008)
dim(libcon_idealpt.2008) <- c(libcon_nsamp.2008,libcon_nresp.2008)
for (i in 1:libcon_nresp.2008){
for (j in 1:libcon_nsamp.2008){
libcon_idealpt.2008[j,i] <- ((libcon_self.2008[i] - libcon_all.a.2008[j,i]) / libcon_all.b.2008[j,i])
}}
#
libcon_idealpt.percentiles.2008 <- rep(NA,libcon_nresp.2008*3)
dim(libcon_idealpt.percentiles.2008) <- c(libcon_nresp.2008,3)
colnames(libcon_idealpt.percentiles.2008) <- c("5%","50%","95%")
for (i in 1:libcon_nresp.2008){
libcon_idealpt.percentiles.2008[i,] <- quantile(libcon_idealpt.2008[,i],probs=c(0.05,0.5,0.95),na.rm=TRUE)
}
libcon_idealpt.estimates.2008 <- libcon_idealpt.percentiles.2008[,2]
#
#
#
# Probit Model of Vote Choie for All Voters
res1.2008 <- glm(formula = libcon_presvote.2008 ~ libcon_self.2008 + libcon_party.2008, family = binomial(link = "probit"))
res2.2008 <- glm(formula = libcon_presvote.2008 ~ libcon_idealpt.estimates.2008 + libcon_party.2008, family = binomial(link = "probit"))
apsrtable(res1.2008,res2.2008)
#
# Probit Model of Vote Choice for Self-Identified Moderates and Slightly Liberal/Conservative
res1.mods.2008 <- glm(formula = libcon_presvote.2008[libcon_self.2008 >= -1 & libcon_self.2008 <= 1] ~ libcon_self.2008[libcon_self.2008 >= -1 & libcon_self.2008 <= 1] + libcon_party.2008[libcon_self.2008 >= -1 & libcon_self.2008 <= 1], family = binomial(link = "probit"))
res2.mods.2008 <- glm(formula = libcon_presvote.2008[libcon_self.2008 >= -1 & libcon_self.2008 <= 1] ~ libcon_idealpt.estimates.2008[libcon_self.2008 >= -1 & libcon_self.2008 <= 1] + libcon_party.2008[libcon_self.2008 >= -1 & libcon_self.2008 <= 1], family = binomial(link = "probit"))
apsrtable(res1.mods.2008,res2.mods.2008)
#
# Probit Model of Vote Choice for Only Self-Identified Moderates 
res1.modsonly.2008 <- glm(formula = libcon_presvote.2008[libcon_self.2008 == 0] ~ libcon_party.2008[libcon_self.2008 == 0], family = binomial(link = "probit"))
res2.modsonly.2008 <- glm(formula = libcon_presvote.2008[libcon_self.2008 == 0] ~ libcon_idealpt.estimates.2008[libcon_self.2008 == 0] + libcon_party.2008[libcon_self.2008 == 0], family = binomial(link = "probit"))
apsrtable(res1.modsonly.2008,res2.modsonly.2008)
#
#
# Overlap Estimates for Presidential Vote
Obama_idealpt.estimates.2008 <- na.omit(libcon_idealpt.estimates.2008[libcon_presvote.2008==0])
McCain_idealpt.estimates.2008 <- na.omit(libcon_idealpt.estimates.2008[libcon_presvote.2008==1])
#
Obama_libcon_self.2008 <- na.omit(libcon_self.2008[libcon_presvote.2008==0])
McCain_libcon_self.2008 <- na.omit(libcon_self.2008[libcon_presvote.2008==1])
#
#
#  Rescale between 0 and 2*pi
range02pi <- function(x){((x-min(x))/(max(x)-min(x)))*(2*pi)}
tmp.2008.correct <- range02pi(c(Obama_idealpt.estimates.2008,McCain_idealpt.estimates.2008))
tmp.2008.raw <- range02pi(c(Obama_libcon_self.2008,McCain_libcon_self.2008))
Obama_idealpt.estimates.2008.rescaled <- tmp.2008.correct[1:length(Obama_idealpt.estimates.2008)]
McCain_idealpt.estimates.2008.rescaled <- tmp.2008.correct[(length(Obama_idealpt.estimates.2008)+1):length(tmp.2008.correct)]
Obama_libcon_self.2008.rescaled <- tmp.2008.raw[1:length(Obama_libcon_self.2008)]
McCain_libcon_self.2008.rescaled <- tmp.2008.raw[(length(Obama_libcon_self.2008)+1):length(tmp.2008.raw)]
#
#  Bootstrapped standard errors for Overlap Estimates
Ntrials <- 1000
overlap.2008_idealpt.estimates <- vector("numeric", Ntrials)
for(i in 1:Ntrials){
overlap.2008_idealpt.estimates[[i]] <- 
	overlapEst(Obama_idealpt.estimates.2008.rescaled[sample(1:length(Obama_idealpt.estimates.2008.rescaled),length(Obama_idealpt.estimates.2008.rescaled),replace=TRUE)],
	McCain_idealpt.estimates.2008.rescaled[sample(1:length(McCain_idealpt.estimates.2008.rescaled),length(McCain_idealpt.estimates.2008.rescaled),replace=TRUE)])[2]
}
point.overlap.2008_idealpt.estimates <- overlapEst(Obama_idealpt.estimates.2008.rescaled,McCain_idealpt.estimates.2008.rescaled)[2]
summary.overlap.2008_idealpt.estimates <- c(point.overlap.2008_idealpt.estimates, 
	sd(overlap.2008_idealpt.estimates),
	(point.overlap.2008_idealpt.estimates-sd(overlap.2008_idealpt.estimates*1.96)),
	(point.overlap.2008_idealpt.estimates+sd(overlap.2008_idealpt.estimates*1.96)))
names(summary.overlap.2008_idealpt.estimates) <- c("point","sd","lower","upper")
#
#
overlap.2008_libcon.self <- vector("numeric", Ntrials)
for(i in 1:Ntrials){
overlap.2008_libcon.self[[i]] <- 
	overlapEst(Obama_libcon_self.2008.rescaled[sample(1:length(Obama_libcon_self.2008.rescaled),length(Obama_libcon_self.2008.rescaled),replace=TRUE)],
	McCain_libcon_self.2008.rescaled[sample(1:length(McCain_libcon_self.2008.rescaled),length(McCain_libcon_self.2008.rescaled),replace=TRUE)])[2]
}
point.overlap.2008_libcon.self <- overlapEst(Obama_libcon_self.2008.rescaled,McCain_libcon_self.2008.rescaled)[2]
summary.overlap.2008_libcon.self <- c(point.overlap.2008_libcon.self, 
	sd(overlap.2008_libcon.self),
	(point.overlap.2008_libcon.self-sd(overlap.2008_libcon.self)*1.96),
	(point.overlap.2008_libcon.self+sd(overlap.2008_libcon.self)*1.96))
names(summary.overlap.2008_libcon.self) <- c("point","sd","lower","upper")
#
print(summary.overlap.2008_idealpt.estimates)
print(summary.overlap.2008_libcon.self)
#
#  USE UNCERTAINTY IN IDEAL POINT ESTIMATES TO DEVELOP
#  BAYESIAN CREDIBLE INTERVALS FOR THE OVERLAP POINT ESTIMATES
  tmpjunk <- matrix(NA, nrow=Ntrials, ncol=libcon_nresp.2008)
	for (i in 1:Ntrials){
	for (j in 1:libcon_nresp.2008){
	tmpjunk[i,j] <- sample(libcon_idealpt.2008[,j], 1)
	}}
  #
  TEST.Obama_libcon_self.2008 <- list(NA)
  for (i in 1:Ntrials){
  TEST.Obama_libcon_self.2008[[i]] <- na.omit(tmpjunk[i,libcon_presvote.2008==0])
  }
  TEST.McCain_libcon_self.2008 <- list(NA)
  for (i in 1:Ntrials){
  TEST.McCain_libcon_self.2008[[i]] <- na.omit(tmpjunk[i,libcon_presvote.2008==1])
  }
  for (i in 1:Ntrials){
  TEST.Obama_libcon_self.2008[[i]] <- TEST.Obama_libcon_self.2008[[i]][TEST.Obama_libcon_self.2008[[i]] < 10 & TEST.Obama_libcon_self.2008[[i]] > -10]
  TEST.McCain_libcon_self.2008[[i]] <- TEST.McCain_libcon_self.2008[[i]][TEST.McCain_libcon_self.2008[[i]] < 10 & TEST.McCain_libcon_self.2008[[i]] > -10]
  }
  #
  TEST.tmp.2008.correct <- list(NA)
  for (i in 1:Ntrials){
  TEST.tmp.2008.correct[[i]] <- range02pi(c(TEST.Obama_libcon_self.2008[[i]],TEST.McCain_libcon_self.2008[[i]]))
  }
  TEST.Obama_idealpt.estimates.2008.rescaled <- list(NA)
  TEST.McCain_idealpt.estimates.2008.rescaled <- list(NA)
  for (i in 1:Ntrials){
  TEST.Obama_idealpt.estimates.2008.rescaled[[i]] <- TEST.tmp.2008.correct[[i]][1:length(TEST.Obama_libcon_self.2008[[i]])]
  TEST.McCain_idealpt.estimates.2008.rescaled[[i]] <- TEST.tmp.2008.correct[[i]][(length(TEST.Obama_libcon_self.2008[[i]])+1):length(TEST.tmp.2008.correct[[i]])]
  }
  #
  TEST.overlap.2008_idealpt.estimates <- vector("numeric", Ntrials)
  for(i in 1:Ntrials){
  TEST.overlap.2008_idealpt.estimates[[i]] <- 
	overlapEst(TEST.Obama_idealpt.estimates.2008.rescaled[[i]],TEST.McCain_idealpt.estimates.2008.rescaled[[i]])[2]
  }
  #
  TEST.summary.overlap.2008_idealpt.estimates <- quantile(TEST.overlap.2008_idealpt.estimates, probs = c(0.025,0.5,0.975))
  #
#
#
vote.summary.overlap.2008_idealpt.estimates <- c(mean(TEST.summary.overlap.2008_idealpt.estimates), sd(TEST.summary.overlap.2008_idealpt.estimates), TEST.summary.overlap.2008_idealpt.estimates[1], TEST.summary.overlap.2008_idealpt.estimates[3])
names(vote.summary.overlap.2008_idealpt.estimates) <- c("point","sd","lower","upper")
#vote.summary.overlap.2008_idealpt.estimates <- summary.overlap.2008_idealpt.estimates
vote.summary.overlap.2008_libcon.self <- summary.overlap.2008_libcon.self
#
# =============================================================
#
# Overlap Estimates for Party Identification
Democrats_idealpt.estimates.2008 <- na.omit(libcon_idealpt.estimates.2008[libcon_three_party.2008==1])
Republicans_idealpt.estimates.2008 <- na.omit(libcon_idealpt.estimates.2008[libcon_three_party.2008==3])
#
Democrats_libcon_self.2008 <- na.omit(libcon_self.2008[libcon_three_party.2008==1])
Republicans_libcon_self.2008 <- na.omit(libcon_self.2008[libcon_three_party.2008==3])
#
#
#  Rescale between 0 and 2*pi
range02pi <- function(x){((x-min(x))/(max(x)-min(x)))*(2*pi)}
tmp.2008.correct <- range02pi(c(Democrats_idealpt.estimates.2008,Republicans_idealpt.estimates.2008))
tmp.2008.raw <- range02pi(c(Democrats_libcon_self.2008,Republicans_libcon_self.2008))
Democrats_idealpt.estimates.2008.rescaled <- tmp.2008.correct[1:length(Democrats_idealpt.estimates.2008)]
Republicans_idealpt.estimates.2008.rescaled <- tmp.2008.correct[(length(Democrats_idealpt.estimates.2008)+1):length(tmp.2008.correct)]
Democrats_libcon_self.2008.rescaled <- tmp.2008.raw[1:length(Democrats_libcon_self.2008)]
Republicans_libcon_self.2008.rescaled <- tmp.2008.raw[(length(Democrats_libcon_self.2008)+1):length(tmp.2008.raw)]
#
#  Bootstrapped standard errors for Overlap Estimates
Ntrials <- 1000
overlap.2008_idealpt.estimates <- vector("numeric", Ntrials)
for(i in 1:Ntrials){
overlap.2008_idealpt.estimates[[i]] <- 
	overlapEst(Democrats_idealpt.estimates.2008.rescaled[sample(1:length(Democrats_idealpt.estimates.2008.rescaled),length(Democrats_idealpt.estimates.2008.rescaled),replace=TRUE)],
	Republicans_idealpt.estimates.2008.rescaled[sample(1:length(Republicans_idealpt.estimates.2008.rescaled),length(Republicans_idealpt.estimates.2008.rescaled),replace=TRUE)])[2]
}
point.overlap.2008_idealpt.estimates <- overlapEst(Democrats_idealpt.estimates.2008.rescaled,Republicans_idealpt.estimates.2008.rescaled)[2]
summary.overlap.2008_idealpt.estimates <- c(point.overlap.2008_idealpt.estimates, 
	sd(overlap.2008_idealpt.estimates),
	(point.overlap.2008_idealpt.estimates-sd(overlap.2008_idealpt.estimates*1.96)),
	(point.overlap.2008_idealpt.estimates+sd(overlap.2008_idealpt.estimates*1.96)))
names(summary.overlap.2008_idealpt.estimates) <- c("point","sd","lower","upper")
#
#
overlap.2008_libcon.self <- vector("numeric", Ntrials)
for(i in 1:Ntrials){
overlap.2008_libcon.self[[i]] <- 
	overlapEst(Democrats_libcon_self.2008.rescaled[sample(1:length(Democrats_libcon_self.2008.rescaled),length(Democrats_libcon_self.2008.rescaled),replace=TRUE)],
	Republicans_libcon_self.2008.rescaled[sample(1:length(Republicans_libcon_self.2008.rescaled),length(Republicans_libcon_self.2008.rescaled),replace=TRUE)])[2]
}
point.overlap.2008_libcon.self <- overlapEst(Democrats_libcon_self.2008.rescaled,Republicans_libcon_self.2008.rescaled)[2]
summary.overlap.2008_libcon.self <- c(point.overlap.2008_libcon.self, 
	sd(overlap.2008_libcon.self),
	(point.overlap.2008_libcon.self-sd(overlap.2008_libcon.self)*1.96),
	(point.overlap.2008_libcon.self+sd(overlap.2008_libcon.self)*1.96))
names(summary.overlap.2008_libcon.self) <- c("point","sd","lower","upper")
#
print(summary.overlap.2008_idealpt.estimates)
print(summary.overlap.2008_libcon.self)
#
#
#  USE UNCERTAINTY IN IDEAL POINT ESTIMATES TO DEVELOP
#  BAYESIAN CREDIBLE INTERVALS FOR THE OVERLAP POINT ESTIMATES
  tmpjunk <- matrix(NA, nrow=Ntrials, ncol=libcon_nresp.2008)
	for (i in 1:Ntrials){
	for (j in 1:libcon_nresp.2008){
	tmpjunk[i,j] <- sample(libcon_idealpt.2008[,j], 1)
	}}
  #
  TEST.Democrats_libcon_self.2008 <- list(NA)
  for (i in 1:Ntrials){
  TEST.Democrats_libcon_self.2008[[i]] <- na.omit(tmpjunk[i,libcon_three_party.2008==1])
  }
  TEST.Republicans_libcon_self.2008 <- list(NA)
  for (i in 1:Ntrials){
  TEST.Republicans_libcon_self.2008[[i]] <- na.omit(tmpjunk[i,libcon_three_party.2008==3])
  }
  for (i in 1:Ntrials){
  TEST.Democrats_libcon_self.2008[[i]] <- TEST.Democrats_libcon_self.2008[[i]][TEST.Democrats_libcon_self.2008[[i]] < 10 & TEST.Democrats_libcon_self.2008[[i]] > -10]
  TEST.Republicans_libcon_self.2008[[i]] <- TEST.Republicans_libcon_self.2008[[i]][TEST.Republicans_libcon_self.2008[[i]] < 10 & TEST.Republicans_libcon_self.2008[[i]] > -10]
  }
  #
  TEST.tmp.2008.correct <- list(NA)
  for (i in 1:Ntrials){
  TEST.tmp.2008.correct[[i]] <- range02pi(c(TEST.Democrats_libcon_self.2008[[i]],TEST.Republicans_libcon_self.2008[[i]]))
  }
  TEST.Democrats_idealpt.estimates.2008.rescaled <- list(NA)
  TEST.Republicans_idealpt.estimates.2008.rescaled <- list(NA)
  for (i in 1:Ntrials){
  TEST.Democrats_idealpt.estimates.2008.rescaled[[i]] <- TEST.tmp.2008.correct[[i]][1:length(TEST.Democrats_libcon_self.2008[[i]])]
  TEST.Republicans_idealpt.estimates.2008.rescaled[[i]] <- TEST.tmp.2008.correct[[i]][(length(TEST.Democrats_libcon_self.2008[[i]])+1):length(TEST.tmp.2008.correct[[i]])]
  }
  #
  TEST.overlap.2008_idealpt.estimates <- vector("numeric", Ntrials)
  for(i in 1:Ntrials){
  TEST.overlap.2008_idealpt.estimates[[i]] <- 
	overlapEst(TEST.Democrats_idealpt.estimates.2008.rescaled[[i]],TEST.Republicans_idealpt.estimates.2008.rescaled[[i]])[2]
  }
  #
  TEST.summary.overlap.2008_idealpt.estimates <- quantile(TEST.overlap.2008_idealpt.estimates, probs = c(0.025,0.5,0.975))
  #
#
party.summary.overlap.2008_idealpt.estimates <- c(mean(TEST.summary.overlap.2008_idealpt.estimates), sd(TEST.summary.overlap.2008_idealpt.estimates), TEST.summary.overlap.2008_idealpt.estimates[1], TEST.summary.overlap.2008_idealpt.estimates[3])
names(party.summary.overlap.2008_idealpt.estimates) <- c("point","sd","lower","upper")
#party.summary.overlap.2008_idealpt.estimates <- summary.overlap.2008_idealpt.estimates
party.summary.overlap.2008_libcon.self <- summary.overlap.2008_libcon.self
#
#
#
#
# =============================================================
# 
# 
#
#
#
# %%%%%%%%%%%%%%%   2004   %%%%%%%%%%%%%%%
#
detach(ANES2008)
ANES2004 <- read.dta("c:/Dropbox/BAM/ANES/ANES_2004/ANES2004.dta", convert.factors=FALSE)
attach(ANES2004)
#
# =============================================================
#       LIBERAL-CONSERVATIVE SCALE RESULTS:
#
libcon.self.2004 <- V043085
libcon.self.2004[libcon.self.2004 >= 80] <- NA
#
libcon.Bush.2004 <- V043087
libcon.Bush.2004[libcon.Bush.2004 >= 8] <- NA
#
libcon.Kerry.2004 <- V043088
libcon.Kerry.2004[libcon.Kerry.2004 >= 8] <- NA
#
libcon.Nader.2004 <- V043089
libcon.Nader.2004[libcon.Nader.2004 >= 8] <- NA
#
libcon.DemParty.2004 <- V043090
libcon.DemParty.2004[libcon.DemParty.2004 >= 8] <- NA
#
libcon.RepParty.2004 <- V043091
libcon.RepParty.2004[libcon.RepParty.2004 >= 8] <- NA
#
#
libcon_T.2004 <- cbind(libcon.self.2004, libcon.Bush.2004, libcon.Kerry.2004, libcon.Nader.2004, libcon.DemParty.2004, libcon.RepParty.2004)
#
libcon_T.2004 <- libcon_T.2004 - 4
#
libcon.self.2004 <- libcon_T.2004[,1]
libcon_TT.2004 <- libcon_T.2004[,-1]
#
cutoff <- 3
libcon_self.2004 <- libcon.self.2004[rowSums(!is.na(libcon_TT.2004)) >= cutoff]
libcon_TTT.2004 <- libcon_TT.2004[rowSums(!is.na(libcon_TT.2004)) >= cutoff,]
#
# PARTY ID: 1 = STRONG DEM, 7 = STRONG REP
libcon_party.2004 <- V043116[rowSums(!is.na(libcon_TT.2004)) >= cutoff]
libcon_party.2004 <- libcon_party.2004 + 1
libcon_party.2004[libcon_party.2004 >= 8] <- NA
#
libcon_three_party.2004 <- libcon_party.2004
libcon_three_party.2004[libcon_three_party.2004==1 | libcon_three_party.2004==2 | libcon_three_party.2004==3] <- 1
libcon_three_party.2004[libcon_three_party.2004==4] <- 2
libcon_three_party.2004[libcon_three_party.2004==5 | libcon_three_party.2004==6 | libcon_three_party.2004==7] <- 3
#
#
# PRESIDENTIAL VOTE: 0 = KERRY, 1 = BUSH
libcon_presvote.2004 <- V043203
libcon_presvote.2004[libcon_presvote.2004 >= 3] <- NA
libcon_presvote.2004 <- libcon_presvote.2004 - 1
libcon_presvote.2004[V045026==3] <- 1
libcon_presvote.2004[V045026==1] <- 0
libcon_presvote.2004 <- libcon_presvote.2004[rowSums(!is.na(libcon_TT.2004)) >= cutoff]
#
#
rm(samps)
load("c:/Dropbox/BAM/ANES/ANES_2004/ANES2004_samps_libcon.Rda")
libcon_samps.2004 <- samps
#
libcon_zhat.2004 <- summary(libcon_samps.2004[,grep("zhat", colnames(libcon_samps.2004[[1]]))])
libcon_sigmaj.sq.2004 <- summary(libcon_samps.2004[,grep("sigmaj.sq", colnames(libcon_samps.2004[[1]]))])
#
libcon_a.2004 <- libcon_samps.2004[,grep(pattern="a[", colnames(libcon_samps.2004[[1]]), fixed=TRUE),]
libcon_b.2004 <- libcon_samps.2004[,grep(pattern="b[", colnames(libcon_samps.2004[[1]]), fixed=TRUE),]
libcon_a_point_est.2004 <- summary(libcon_a.2004)[[1]][,1]
libcon_b_point_est.2004 <- summary(libcon_b.2004)[[1]][,1]
libcon_all.a.2004 <- do.call("rbind", libcon_a.2004)
libcon_all.b.2004 <- do.call("rbind", libcon_b.2004)
#
libcon_nsamp.2004 <- nrow(libcon_all.a.2004)
libcon_nresp.2004 <- ncol(libcon_all.a.2004)
#
libcon_idealpt.2004 <- rep(NA, libcon_nsamp.2004*libcon_nresp.2004)
dim(libcon_idealpt.2004) <- c(libcon_nsamp.2004,libcon_nresp.2004)
for (i in 1:libcon_nresp.2004){
for (j in 1:libcon_nsamp.2004){
libcon_idealpt.2004[j,i] <- ((libcon_self.2004[i] - libcon_all.a.2004[j,i]) / libcon_all.b.2004[j,i])
}}
#
libcon_idealpt.percentiles.2004 <- rep(NA,libcon_nresp.2004*3)
dim(libcon_idealpt.percentiles.2004) <- c(libcon_nresp.2004,3)
colnames(libcon_idealpt.percentiles.2004) <- c("5%","50%","95%")
for (i in 1:libcon_nresp.2004){
libcon_idealpt.percentiles.2004[i,] <- quantile(libcon_idealpt.2004[,i],probs=c(0.05,0.5,0.95),na.rm=TRUE)
}
libcon_idealpt.estimates.2004 <- libcon_idealpt.percentiles.2004[,2]
#
#
#
# Probit Model of Vote Choie for All Voters
res1.2004 <- glm(formula = libcon_presvote.2004 ~ libcon_self.2004 + libcon_party.2004, family = binomial(link = "probit"))
res2.2004 <- glm(formula = libcon_presvote.2004 ~ libcon_idealpt.estimates.2004 + libcon_party.2004, family = binomial(link = "probit"))
apsrtable(res1.2004,res2.2004)
#
# Probit Model of Vote Choice for Self-Identified Moderates and Slightly Liberal/Conservative
res1.mods.2004 <- glm(formula = libcon_presvote.2004[libcon_self.2004 >= -1 & libcon_self.2004 <= 1] ~ libcon_self.2004[libcon_self.2004 >= -1 & libcon_self.2004 <= 1] + libcon_party.2004[libcon_self.2004 >= -1 & libcon_self.2004 <= 1], family = binomial(link = "probit"))
res2.mods.2004 <- glm(formula = libcon_presvote.2004[libcon_self.2004 >= -1 & libcon_self.2004 <= 1] ~ libcon_idealpt.estimates.2004[libcon_self.2004 >= -1 & libcon_self.2004 <= 1] + libcon_party.2004[libcon_self.2004 >= -1 & libcon_self.2004 <= 1], family = binomial(link = "probit"))
apsrtable(res1.mods.2004,res2.mods.2004)
#
# Probit Model of Vote Choice for Only Self-Identified Moderates 
res1.modsonly.2004 <- glm(formula = libcon_presvote.2004[libcon_self.2004 == 0] ~ libcon_party.2004[libcon_self.2004 == 0], family = binomial(link = "probit"))
res2.modsonly.2004 <- glm(formula = libcon_presvote.2004[libcon_self.2004 == 0] ~ libcon_idealpt.estimates.2004[libcon_self.2004 == 0] + libcon_party.2004[libcon_self.2004 == 0], family = binomial(link = "probit"))
apsrtable(res1.modsonly.2004,res2.modsonly.2004)
#
#
# Overlap Estimates for Presidential Vote
Obama_idealpt.estimates.2004 <- na.omit(libcon_idealpt.estimates.2004[libcon_presvote.2004==0])
McCain_idealpt.estimates.2004 <- na.omit(libcon_idealpt.estimates.2004[libcon_presvote.2004==1])
#
Obama_libcon_self.2004 <- na.omit(libcon_self.2004[libcon_presvote.2004==0])
McCain_libcon_self.2004 <- na.omit(libcon_self.2004[libcon_presvote.2004==1])
#
#
#  Rescale between 0 and 2*pi
range02pi <- function(x){((x-min(x))/(max(x)-min(x)))*(2*pi)}
tmp.2004.correct <- range02pi(c(Obama_idealpt.estimates.2004,McCain_idealpt.estimates.2004))
tmp.2004.raw <- range02pi(c(Obama_libcon_self.2004,McCain_libcon_self.2004))
Obama_idealpt.estimates.2004.rescaled <- tmp.2004.correct[1:length(Obama_idealpt.estimates.2004)]
McCain_idealpt.estimates.2004.rescaled <- tmp.2004.correct[(length(Obama_idealpt.estimates.2004)+1):length(tmp.2004.correct)]
Obama_libcon_self.2004.rescaled <- tmp.2004.raw[1:length(Obama_libcon_self.2004)]
McCain_libcon_self.2004.rescaled <- tmp.2004.raw[(length(Obama_libcon_self.2004)+1):length(tmp.2004.raw)]
#
#  Bootstrapped standard errors for Overlap Estimates
Ntrials <- 1000
overlap.2004_idealpt.estimates <- vector("numeric", Ntrials)
for(i in 1:Ntrials){
overlap.2004_idealpt.estimates[[i]] <- 
	overlapEst(Obama_idealpt.estimates.2004.rescaled[sample(1:length(Obama_idealpt.estimates.2004.rescaled),length(Obama_idealpt.estimates.2004.rescaled),replace=TRUE)],
	McCain_idealpt.estimates.2004.rescaled[sample(1:length(McCain_idealpt.estimates.2004.rescaled),length(McCain_idealpt.estimates.2004.rescaled),replace=TRUE)])[2]
}
point.overlap.2004_idealpt.estimates <- overlapEst(Obama_idealpt.estimates.2004.rescaled,McCain_idealpt.estimates.2004.rescaled)[2]
summary.overlap.2004_idealpt.estimates <- c(point.overlap.2004_idealpt.estimates, 
	sd(overlap.2004_idealpt.estimates),
	(point.overlap.2004_idealpt.estimates-sd(overlap.2004_idealpt.estimates*1.96)),
	(point.overlap.2004_idealpt.estimates+sd(overlap.2004_idealpt.estimates*1.96)))
names(summary.overlap.2004_idealpt.estimates) <- c("point","sd","lower","upper")
#
#
overlap.2004_libcon.self <- vector("numeric", Ntrials)
for(i in 1:Ntrials){
overlap.2004_libcon.self[[i]] <- 
	overlapEst(Obama_libcon_self.2004.rescaled[sample(1:length(Obama_libcon_self.2004.rescaled),length(Obama_libcon_self.2004.rescaled),replace=TRUE)],
	McCain_libcon_self.2004.rescaled[sample(1:length(McCain_libcon_self.2004.rescaled),length(McCain_libcon_self.2004.rescaled),replace=TRUE)])[2]
}
point.overlap.2004_libcon.self <- overlapEst(Obama_libcon_self.2004.rescaled,McCain_libcon_self.2004.rescaled)[2]
summary.overlap.2004_libcon.self <- c(point.overlap.2004_libcon.self, 
	sd(overlap.2004_libcon.self),
	(point.overlap.2004_libcon.self-sd(overlap.2004_libcon.self)*1.96),
	(point.overlap.2004_libcon.self+sd(overlap.2004_libcon.self)*1.96))
names(summary.overlap.2004_libcon.self) <- c("point","sd","lower","upper")
#
print(summary.overlap.2004_idealpt.estimates)
print(summary.overlap.2004_libcon.self)
#
#  USE UNCERTAINTY IN IDEAL POINT ESTIMATES TO DEVELOP
#  BAYESIAN CREDIBLE INTERVALS FOR THE OVERLAP POINT ESTIMATES
  tmpjunk <- matrix(NA, nrow=Ntrials, ncol=libcon_nresp.2004)
	for (i in 1:Ntrials){
	for (j in 1:libcon_nresp.2004){
	tmpjunk[i,j] <- sample(libcon_idealpt.2004[,j], 1)
	}}
  #
  TEST.Kerry_libcon_self.2004 <- list(NA)
  for (i in 1:Ntrials){
  TEST.Kerry_libcon_self.2004[[i]] <- na.omit(tmpjunk[i,libcon_presvote.2004==0])
  }
  TEST.Bush_libcon_self.2004 <- list(NA)
  for (i in 1:Ntrials){
  TEST.Bush_libcon_self.2004[[i]] <- na.omit(tmpjunk[i,libcon_presvote.2004==1])
  }
  for (i in 1:Ntrials){
  TEST.Kerry_libcon_self.2004[[i]] <- TEST.Kerry_libcon_self.2004[[i]][TEST.Kerry_libcon_self.2004[[i]] < 10 & TEST.Kerry_libcon_self.2004[[i]] > -10]
  TEST.Bush_libcon_self.2004[[i]] <- TEST.Bush_libcon_self.2004[[i]][TEST.Bush_libcon_self.2004[[i]] < 10 & TEST.Bush_libcon_self.2004[[i]] > -10]
  }
  #
  TEST.tmp.2004.correct <- list(NA)
  for (i in 1:Ntrials){
  TEST.tmp.2004.correct[[i]] <- range02pi(c(TEST.Kerry_libcon_self.2004[[i]],TEST.Bush_libcon_self.2004[[i]]))
  }
  TEST.Kerry_idealpt.estimates.2004.rescaled <- list(NA)
  TEST.Bush_idealpt.estimates.2004.rescaled <- list(NA)
  for (i in 1:Ntrials){
  TEST.Kerry_idealpt.estimates.2004.rescaled[[i]] <- TEST.tmp.2004.correct[[i]][1:length(TEST.Kerry_libcon_self.2004[[i]])]
  TEST.Bush_idealpt.estimates.2004.rescaled[[i]] <- TEST.tmp.2004.correct[[i]][(length(TEST.Kerry_libcon_self.2004[[i]])+1):length(TEST.tmp.2004.correct[[i]])]
  }
  #
  TEST.overlap.2004_idealpt.estimates <- vector("numeric", Ntrials)
  for(i in 1:Ntrials){
  TEST.overlap.2004_idealpt.estimates[[i]] <- 
	overlapEst(TEST.Kerry_idealpt.estimates.2004.rescaled[[i]],TEST.Bush_idealpt.estimates.2004.rescaled[[i]])[2]
  }
  #
  TEST.summary.overlap.2004_idealpt.estimates <- quantile(TEST.overlap.2004_idealpt.estimates, probs = c(0.025,0.5,0.975))
  #
#
#
vote.summary.overlap.2004_idealpt.estimates <- c(mean(TEST.summary.overlap.2004_idealpt.estimates), sd(TEST.summary.overlap.2004_idealpt.estimates), TEST.summary.overlap.2004_idealpt.estimates[1], TEST.summary.overlap.2004_idealpt.estimates[3])
names(vote.summary.overlap.2004_idealpt.estimates) <- c("point","sd","lower","upper")
#vote.summary.overlap.2004_idealpt.estimates <- summary.overlap.2004_idealpt.estimates
vote.summary.overlap.2004_libcon.self <- summary.overlap.2004_libcon.self
#
# =============================================================
#
# Overlap Estimates for Party Identification
Democrats_idealpt.estimates.2004 <- na.omit(libcon_idealpt.estimates.2004[libcon_three_party.2004==1])
Republicans_idealpt.estimates.2004 <- na.omit(libcon_idealpt.estimates.2004[libcon_three_party.2004==3])
#
Democrats_libcon_self.2004 <- na.omit(libcon_self.2004[libcon_three_party.2004==1])
Republicans_libcon_self.2004 <- na.omit(libcon_self.2004[libcon_three_party.2004==3])
#
#
#  Rescale between 0 and 2*pi
range02pi <- function(x){((x-min(x))/(max(x)-min(x)))*(2*pi)}
tmp.2004.correct <- range02pi(c(Democrats_idealpt.estimates.2004,Republicans_idealpt.estimates.2004))
tmp.2004.raw <- range02pi(c(Democrats_libcon_self.2004,Republicans_libcon_self.2004))
Democrats_idealpt.estimates.2004.rescaled <- tmp.2004.correct[1:length(Democrats_idealpt.estimates.2004)]
Republicans_idealpt.estimates.2004.rescaled <- tmp.2004.correct[(length(Democrats_idealpt.estimates.2004)+1):length(tmp.2004.correct)]
Democrats_libcon_self.2004.rescaled <- tmp.2004.raw[1:length(Democrats_libcon_self.2004)]
Republicans_libcon_self.2004.rescaled <- tmp.2004.raw[(length(Democrats_libcon_self.2004)+1):length(tmp.2004.raw)]
#
#  Bootstrapped standard errors for Overlap Estimates
Ntrials <- 1000
overlap.2004_idealpt.estimates <- vector("numeric", Ntrials)
for(i in 1:Ntrials){
overlap.2004_idealpt.estimates[[i]] <- 
	overlapEst(Democrats_idealpt.estimates.2004.rescaled[sample(1:length(Democrats_idealpt.estimates.2004.rescaled),length(Democrats_idealpt.estimates.2004.rescaled),replace=TRUE)],
	Republicans_idealpt.estimates.2004.rescaled[sample(1:length(Republicans_idealpt.estimates.2004.rescaled),length(Republicans_idealpt.estimates.2004.rescaled),replace=TRUE)])[2]
}
point.overlap.2004_idealpt.estimates <- overlapEst(Democrats_idealpt.estimates.2004.rescaled,Republicans_idealpt.estimates.2004.rescaled)[2]
summary.overlap.2004_idealpt.estimates <- c(point.overlap.2004_idealpt.estimates, 
	sd(overlap.2004_idealpt.estimates),
	(point.overlap.2004_idealpt.estimates-sd(overlap.2004_idealpt.estimates*1.96)),
	(point.overlap.2004_idealpt.estimates+sd(overlap.2004_idealpt.estimates*1.96)))
names(summary.overlap.2004_idealpt.estimates) <- c("point","sd","lower","upper")
#
#
overlap.2004_libcon.self <- vector("numeric", Ntrials)
for(i in 1:Ntrials){
overlap.2004_libcon.self[[i]] <- 
	overlapEst(Democrats_libcon_self.2004.rescaled[sample(1:length(Democrats_libcon_self.2004.rescaled),length(Democrats_libcon_self.2004.rescaled),replace=TRUE)],
	Republicans_libcon_self.2004.rescaled[sample(1:length(Republicans_libcon_self.2004.rescaled),length(Republicans_libcon_self.2004.rescaled),replace=TRUE)])[2]
}
point.overlap.2004_libcon.self <- overlapEst(Democrats_libcon_self.2004.rescaled,Republicans_libcon_self.2004.rescaled)[2]
summary.overlap.2004_libcon.self <- c(point.overlap.2004_libcon.self, 
	sd(overlap.2004_libcon.self),
	(point.overlap.2004_libcon.self-sd(overlap.2004_libcon.self)*1.96),
	(point.overlap.2004_libcon.self+sd(overlap.2004_libcon.self)*1.96))
names(summary.overlap.2004_libcon.self) <- c("point","sd","lower","upper")
#
print(summary.overlap.2004_idealpt.estimates)
print(summary.overlap.2004_libcon.self)
#
#  USE UNCERTAINTY IN IDEAL POINT ESTIMATES TO DEVELOP
#  BAYESIAN CREDIBLE INTERVALS FOR THE OVERLAP POINT ESTIMATES
  tmpjunk <- matrix(NA, nrow=Ntrials, ncol=libcon_nresp.2004)
	for (i in 1:Ntrials){
	for (j in 1:libcon_nresp.2004){
	tmpjunk[i,j] <- sample(libcon_idealpt.2004[,j], 1)
	}}
  #
  TEST.Democrats_libcon_self.2004 <- list(NA)
  for (i in 1:Ntrials){
  TEST.Democrats_libcon_self.2004[[i]] <- na.omit(tmpjunk[i,libcon_three_party.2004==1])
  }
  TEST.Republicans_libcon_self.2004 <- list(NA)
  for (i in 1:Ntrials){
  TEST.Republicans_libcon_self.2004[[i]] <- na.omit(tmpjunk[i,libcon_three_party.2004==3])
  }
  for (i in 1:Ntrials){
  TEST.Democrats_libcon_self.2004[[i]] <- TEST.Democrats_libcon_self.2004[[i]][TEST.Democrats_libcon_self.2004[[i]] < 10 & TEST.Democrats_libcon_self.2004[[i]] > -10]
  TEST.Republicans_libcon_self.2004[[i]] <- TEST.Republicans_libcon_self.2004[[i]][TEST.Republicans_libcon_self.2004[[i]] < 10 & TEST.Republicans_libcon_self.2004[[i]] > -10]
  }
  #
  TEST.tmp.2004.correct <- list(NA)
  for (i in 1:Ntrials){
  TEST.tmp.2004.correct[[i]] <- range02pi(c(TEST.Democrats_libcon_self.2004[[i]],TEST.Republicans_libcon_self.2004[[i]]))
  }
  TEST.Democrats_idealpt.estimates.2004.rescaled <- list(NA)
  TEST.Republicans_idealpt.estimates.2004.rescaled <- list(NA)
  for (i in 1:Ntrials){
  TEST.Democrats_idealpt.estimates.2004.rescaled[[i]] <- TEST.tmp.2004.correct[[i]][1:length(TEST.Democrats_libcon_self.2004[[i]])]
  TEST.Republicans_idealpt.estimates.2004.rescaled[[i]] <- TEST.tmp.2004.correct[[i]][(length(TEST.Democrats_libcon_self.2004[[i]])+1):length(TEST.tmp.2004.correct[[i]])]
  }
  #
  TEST.overlap.2004_idealpt.estimates <- vector("numeric", Ntrials)
  for(i in 1:Ntrials){
  TEST.overlap.2004_idealpt.estimates[[i]] <- 
	overlapEst(TEST.Democrats_idealpt.estimates.2004.rescaled[[i]],TEST.Republicans_idealpt.estimates.2004.rescaled[[i]])[2]
  }
  #
  TEST.summary.overlap.2004_idealpt.estimates <- quantile(TEST.overlap.2004_idealpt.estimates, probs = c(0.025,0.5,0.975))
  #
#
party.summary.overlap.2004_idealpt.estimates <- c(mean(TEST.summary.overlap.2004_idealpt.estimates), sd(TEST.summary.overlap.2004_idealpt.estimates), TEST.summary.overlap.2004_idealpt.estimates[1], TEST.summary.overlap.2004_idealpt.estimates[3])
names(party.summary.overlap.2004_idealpt.estimates) <- c("point","sd","lower","upper")
party.summary.overlap.2004_libcon.self <- summary.overlap.2004_libcon.self
#
#
#
#
#
#
# =============================================================
#
#
# ALL VOTERS
apsrtable(res1.2004,res2.2004,res1.2008,res2.2008,res1.2012,res2.2012)
#
# SELF-IDENTIFIED MODERATES AND SLIGHTLY LIBERAL/CONSERVATIVE VOTERS
apsrtable(res1.mods.2004,res2.mods.2004,res1.mods.2008,res2.mods.2008,res1.mods.2012,res2.mods.2012)
#
pre(res1.mods.2004, mod2=NULL)
pre(res2.mods.2004, mod2=NULL)
pre(res1.mods.2008, mod2=NULL)
pre(res2.mods.2008, mod2=NULL)
pre(res1.mods.2012, mod2=NULL)
pre(res2.mods.2012, mod2=NULL)
#  
#
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# PLOTS ***** PLOTS ***** PLOTS ***** PLOTS ***** PLOTS ***** PLOTS 
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
#  DOT PLOT OF !!! PARTY !!! OVERLAP ESTIMATES
#
A1.04 <- rbind(party.summary.overlap.2004_libcon.self, party.summary.overlap.2004_idealpt.estimates)
A1.04 <- cbind(2004,A1.04)
#
A1.08 <- rbind(party.summary.overlap.2008_libcon.self, party.summary.overlap.2008_idealpt.estimates)
A1.08 <- cbind(2008,A1.08)
#
A1.12 <- rbind(party.summary.overlap.2012_libcon.self, party.summary.overlap.2012_idealpt.estimates)
A1.12 <- cbind(2012,A1.12)
#
A2 <- rbind(A1.04, A1.08, A1.12)
#
year <- A2[,1]
overlap <- A2[,2]
dd <- A2[,4]
ff <- A2[,5]
groups1 <- rep(c("Raw Placements","BAM Estimates"),3)
groups2 <- factor(groups1, levels=c("Raw Placements","BAM Estimates"))
#
#
#
dotplot(overlap ~ groups2 | as.factor(year), horizontal = FALSE,
	par.strip.text=list(cex=1.5),
	scales=list(cex=1.4),
	main=list(
		label="Democratic and Republican Identifiers (with Leaners)\nOverlap Point Estimates with 95% Confidence/Credible Intervals",
		cex=1.8),
	ylab=list(
		label="Overlap",
		cex=1.7),
	lower=dd,upper=ff,
	panel = function(x, y, subscripts, lower, upper, ...){
		panel.abline(h=seq(0.1,1,0.1), lwd=1, lty=1, col="light grey")
		panel.segments(x, lower[subscripts],
		x, upper[subscripts], col="black",lwd=3)
		panel.xyplot(x,y,pch=16,col="black",cex=1.6)
	},
	layout=c(3,1), ylim=c(0,0.65))
#
#
#  DOT PLOT OF !!! VOTE !!! OVERLAP ESTIMATES
#
A1.04 <- rbind(vote.summary.overlap.2004_libcon.self, vote.summary.overlap.2004_idealpt.estimates)
A1.04 <- cbind(2004,A1.04)
#
A1.08 <- rbind(vote.summary.overlap.2008_libcon.self, vote.summary.overlap.2008_idealpt.estimates)
A1.08 <- cbind(2008,A1.08)
#
A1.12 <- rbind(vote.summary.overlap.2012_libcon.self, vote.summary.overlap.2012_idealpt.estimates)
A1.12 <- cbind(2012,A1.12)
#
A2 <- rbind(A1.04, A1.08, A1.12)
#
year <- A2[,1]
overlap <- A2[,2]
dd <- A2[,4]
ff <- A2[,5]
groups1 <- rep(c("Raw Placements","BAM Estimates"),3)
groups2 <- factor(groups1, levels=c("Raw Placements","BAM Estimates"))
#
#
dotplot(overlap ~ groups2 | as.factor(year), horizontal = FALSE,
	par.strip.text=list(cex=1.5),
	scales=list(cex=1.4),
	main=list(
		label="Democratic and Republican Presidential Vote Choice\nOverlap Point Estimates with 95% Confidence/Credible Intervals",
		cex=1.8),
	ylab=list(
		label="Overlap",
		cex=1.7),
	lower=dd,upper=ff,
	panel = function(x, y, subscripts, lower, upper, ...){
		panel.abline(h=seq(0.1,1,0.1), lwd=1, lty=1, col="light grey")
		panel.segments(x, lower[subscripts],
		x, upper[subscripts], col="black",lwd=3)
		panel.xyplot(x,y,pch=16,col="black",cex=1.6)
	},
	layout=c(3,1), ylim=c(0,0.65))
#
#
# ===================================================
#
#  DENSITY PLOT OF THE ORIGINAL AND RE-SCALED IDEAL POINT ESTIMATES
#
#   %%2012%%
libcon.all.2012 <- na.omit(libcon_idealpt.estimates.2012)
libcon.positive.weights.2012 <- na.omit(libcon_idealpt.estimates.2012[libcon_b_point_est.2012 > 0])
#
libcon.all.Share.2012 <- length(libcon.all.2012)/(length(libcon.all.2012)+length(libcon.positive.weights.2012))
libcon.positive.weights.Share.2012 <- length(libcon.positive.weights.2012)/(length(libcon.all.2012)+length(libcon.positive.weights.2012))
#
libcon.all.dens.2012 <- density(libcon.all.2012)
libcon.all.dens.2012$y <- libcon.all.dens.2012$y*libcon.all.Share.2012
#
libcon.positive.weights.dens.2012 <- density(libcon.positive.weights.2012)
libcon.positive.weights.dens.2012$y <- libcon.positive.weights.dens.2012$y*libcon.positive.weights.Share.2012
#
#
pdf('c:/Dropbox/BAM/ANES/ANES2012_raw_corrected_libcon.pdf', width=8, height=4)
#
par(mfrow=c(1,2))
#
tmp <- hist(libcon_self.2012, breaks=-3:(max(libcon_self.2012, na.rm=TRUE)+1), 
	xaxt="n", yaxt="n", main="", xlab="", ylab="",
	right=FALSE, freq=FALSE, col="gray")
axis(1, at=tmp$mids, labels=1:7)
axis(2, las=1, at=seq(0,0.3,0.1))
# Main title
mtext("Raw Self-Placements",side=3,line=1.50,cex=1.2,font=2)
# x-axis title
mtext("Liberal-Conservative",side=1,line=3.25,cex=1.1)
# y-axis title
mtext("Density",side=2,line=2.75,cex=1.1)
#
#
plot(libcon.all.dens.2012, col="black", lwd=2, xlim=c(-3,3),
	yaxt="n", main="", xlab="", ylab="")
axis(2, las=1, at=seq(0,0.2,0.1))
# Main title
mtext("BAM Ideal Point Estimates",side=3,line=1.50,cex=1.2,font=2)
# x-axis title
mtext("Liberal-Conservative\nBlack = All, Gray = Positive Weights",side=1,line=3.75,cex=1.1)
# y-axis title
mtext("Density",side=2,line=2.75,cex=1.1)
lines(libcon.positive.weights.dens.2012, col="gray", lwd=2)
text(libcon_zhat.2012[[1]][1,1],0.01,"O",col="blue",font=2,cex=1.4)
text(libcon_zhat.2012[[1]][2,1],0.01,"R",col="red",font=2,cex=1.4)
#
dev.off()
#
#
#  BOX-AND-WHISKER PLOTS OF ALPHA VALUES BY IDEOLOGICAL
#  SELF-PLACEMENTS AND PARTY IDENTIFICATION
#
avalues <- summary(libcon_a.2012)[[1]][,1]
#
alpha.libcon.plot <- bwplot(avalues ~ libcon_self.2012, horizontal=FALSE,
	main=list(label="Ideological Self-Placement and Alpha (Shift) Term",cex=1.4),
	xlab=list(label="Self-Placement (Very Liberal to Very Conservative)",cex=1.3),
	ylab=list(label="Alpha Posterior Mean",cex=1.3),
	scales=list(cex=1.2)
	)
#
update(alpha.libcon.plot, panel = function(...) {
            panel.abline(h=0, lwd=2, lty=1, col="light grey")
            panel.bwplot(...)
        })
#
alpha.party.plot <- bwplot(avalues ~ libcon_party.2012, horizontal=FALSE,
	main=list(label="Party Identification and Alpha (Shift) Term",cex=1.4),
	xlab=list(label="Party ID (Strong Democrat to Strong Republican)",cex=1.3),
	ylab=list(label="Alpha Posterior Mean",cex=1.3),
	scales=list(cex=1.2)
	)
#
update(alpha.party.plot, panel = function(...) {
            panel.abline(h=0, lwd=2, lty=1, col="light grey")
            panel.bwplot(...)
        })
#
#
#  OLS REGRESSION OF ALPHA VALUES ON IDEOLOGY
#
library(apsrtable)
#
alpha.ideology.ols.1 <- lm(avalues ~ libcon_self.2012)
apsrtable(alpha.ideology.ols.1)
#
#
#
#  BOX-AND-WHISKER PLOTS OF BETA VALUES BY IDEOLOGICAL
#  SELF-PLACEMENTS AND PARTY IDENTIFICATION
#
bvalues <- summary(libcon_b.2012)[[1]][,1]
#
beta.libcon.plot <- bwplot(bvalues ~ libcon_self.2012, horizontal=FALSE,
	main=list(label="Ideological Self-Placement and Beta (Stretch) Term",cex=1.4),
	xlab=list(label="Self-Placement (Very Liberal to Very Conservative)",cex=1.3),
	ylab=list(label="Beta Posterior Mean",cex=1.3),
	scales=list(cex=1.2)
	)
#
update(beta.libcon.plot, panel = function(...) {
            panel.abline(h=0, lwd=2, lty=1, col="light grey")
            panel.bwplot(...)
        })
#
beta.party.plot <- bwplot(bvalues ~ libcon_party.2012, horizontal=FALSE,
	main=list(label="Party Identification and Beta (Stretch) Term",cex=1.4),
	xlab=list(label="Party ID (Strong Democrat to Strong Republican)",cex=1.3),
	ylab=list(label="Beta Posterior Mean",cex=1.3),
	scales=list(cex=1.2)
	)
#
update(beta.party.plot, panel = function(...) {
            panel.abline(h=0, lwd=2, lty=1, col="light grey")
            panel.bwplot(...)
        })
#
#
#
#
#  POLITICAL KNOWLEDGE BY NEGATIVE/POSITIVE BETA VALUES
#
bvalues.pos <- bvalues
bvalues.pos[bvalues.pos < 0] <- 0
bvalues.pos[bvalues.pos > 0] <- 1
#
knowledge.housemaj <- knowl_housemaj[rowSums(!is.na(libcon_TT.2012)) >= cutoff]
knowledge.housemaj[knowledge.housemaj == -6] <- NA
knowledge.housemaj[knowledge.housemaj == -7] <- NA
knowledge.housemaj[knowledge.housemaj == 1] <- 0
knowledge.housemaj[knowledge.housemaj == -9] <- 0
knowledge.housemaj[knowledge.housemaj == -8] <- 0
knowledge.housemaj[knowledge.housemaj == 2] <- 1
#
knowledge.senmaj <- knowl_senmaj[rowSums(!is.na(libcon_TT.2012)) >= cutoff]
knowledge.senmaj[knowledge.senmaj == -6] <- NA
knowledge.senmaj[knowledge.senmaj == -7] <- NA
knowledge.senmaj[knowledge.senmaj == 1] <- 1
knowledge.senmaj[knowledge.senmaj == -9] <- 0
knowledge.senmaj[knowledge.senmaj == -8] <- 0
knowledge.senmaj[knowledge.senmaj == 2] <- 0
#
table(bvalues.pos,knowledge.housemaj)
table(bvalues.pos,knowledge.senmaj)
#
#
#
#
#  ***ALTERNATE*** DENSITY PLOT OF THE ORIGINAL AND RE-SCALED IDEAL POINT ESTIMATES
#
#   %%2012%%
libcon.all.2012 <- na.omit(libcon_idealpt.estimates.2012)
libcon.positive.weights.2012 <- na.omit(libcon_idealpt.estimates.2012[libcon_b_point_est.2012 > 0])
#
libcon.all.Share.2012 <- length(libcon.all.2012)/(length(libcon.all.2012)+length(libcon.positive.weights.2012))
libcon.positive.weights.Share.2012 <- length(libcon.positive.weights.2012)/(length(libcon.all.2012)+length(libcon.positive.weights.2012))
#
libcon.all.dens.2012 <- density(libcon.all.2012)
libcon.all.dens.2012$y <- libcon.all.dens.2012$y*libcon.all.Share.2012
#
libcon.positive.weights.dens.2012 <- density(libcon.positive.weights.2012)
libcon.positive.weights.dens.2012$y <- libcon.positive.weights.dens.2012$y*libcon.positive.weights.Share.2012
#
#
pdf('c:/Dropbox/BAM/ANES/ANES2012_raw_corrected_libcon_alternate.pdf', width=8, height=4)
#
par(mfrow=c(1,2))
#
tmp <- hist(libcon_self.2012, breaks=-3:(max(libcon_self.2012, na.rm=TRUE)+1), 
	xaxt="n", yaxt="n", main="", xlab="", ylab="",
	right=FALSE, freq=TRUE, col="gray", ylim=c(-30,1700))
axis(1, at=tmp$mids, labels=1:7)
axis(2, las=1, at=seq(0,1750,250), labels=c("0","","500","","1000","","1500",""))
# Main title
mtext("Raw Self-Placements",side=3,line=1.50,cex=1.2,font=2)
# x-axis title
mtext("Liberal-Conservative",side=1,line=3.25,cex=1.1)
# y-axis title
mtext("Frequency",side=2,line=2.75,cex=1.1)
#
#
libcon.all.2012.HIST <- libcon.all.2012
libcon.all.2012.HIST[libcon.all.2012.HIST < -3] <- -3
libcon.all.2012.HIST[libcon.all.2012.HIST > 3] <- 3
#
libcon.pos.2012.HIST <- libcon.positive.weights.2012
libcon.pos.2012.HIST[libcon.pos.2012.HIST < -3] <- -3
libcon.pos.2012.HIST[libcon.pos.2012.HIST > 3] <- 3
#
tmp.2 <- hist(libcon.all.2012.HIST, breaks=seq(-3,3,by=0.25),
	xlab="",ylab="",main="",
	right=FALSE, freq=TRUE, col="gray75", ylim=c(-30,600))
tmp.3 <- hist(libcon.pos.2012.HIST, breaks=seq(-3,3,by=0.25),
	right=FALSE,freq=TRUE, col="gray50", add=TRUE)
text(libcon_zhat.2012[[1]][1,1],-27,"O",col="blue",font=2,cex=1.2)
text(libcon_zhat.2012[[1]][2,1],-27,"R",col="red",font=2,cex=1.2)
#
# Main title
mtext("BAM Ideal Point Estimates",side=3,line=1.50,cex=1.2,font=2)
# x-axis title
mtext("Liberal-Conservative",side=1,line=3.25,cex=1.1)
# y-axis title
mtext("Frequency",side=2,line=2.75,cex=1.1)
#
#
dev.off()
#
#
#
#
# Correlations between zhat and respondent placements
#
stimuli.est <- libcon_zhat.2012[[1]][,1]
political.info <- apply(z, 1, function(x){cor(x,stimuli.est,use="complete")})
hist(political.info,breaks=20,col="blue",ylim=c(0,4500),xlab="Correlation",main="Histogram of Correlations between Respondent Placements and\nStimuli Estimates in the 2012 ANES")
#
knowledge.housemaj <- knowl_housemaj[rowSums(!is.na(TT)) >= cutoff]
knowledge.housemaj[knowledge.housemaj == -6] <- NA
knowledge.housemaj[knowledge.housemaj == -7] <- NA
knowledge.housemaj[knowledge.housemaj == 1] <- 0
knowledge.housemaj[knowledge.housemaj == -9] <- 0
knowledge.housemaj[knowledge.housemaj == -8] <- 0
knowledge.housemaj[knowledge.housemaj == 2] <- 1
#
knowledge.senmaj <- knowl_senmaj[rowSums(!is.na(TT)) >= cutoff]
knowledge.senmaj[knowledge.senmaj == -6] <- NA
knowledge.senmaj[knowledge.senmaj == -7] <- NA
knowledge.senmaj[knowledge.senmaj == 1] <- 1
knowledge.senmaj[knowledge.senmaj == -9] <- 0
knowledge.senmaj[knowledge.senmaj == -8] <- 0
knowledge.senmaj[knowledge.senmaj == 2] <- 0
#
table(knowledge.housemaj[political.info <= -0.9])
table(knowledge.senmaj[political.info <= -0.9])
#
table(knowledge.housemaj[political.info > -0.9 & political.info < 0])
table(knowledge.senmaj[political.info > -0.9 & political.info < 0])
#
table(knowledge.housemaj[political.info >= 0])
table(knowledge.senmaj[political.info >= 0])
#
prop.test(rbind(table(knowledge.housemaj[political.info <= -0.9]), table(knowledge.housemaj[political.info > -0.9 & political.info < 0])))
prop.test(rbind(table(knowledge.senmaj[political.info <= -0.9]), table(knowledge.senmaj[political.info > -0.9 & political.info < 0])))
#