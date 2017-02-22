#
# BAM_CCES2010_Senate.r
#
library(foreign)
library(rjags)
library(coda)
library(basicspace)
library(ggplot2)
library(lattice)
library(maps)
library(maptools)
library(rgdal)
#
set.seed(1985)
#
CCES2010 <- read.dta("c:/Dropbox/BAM/CCES_2010/CCES2010.dta", convert.factors=FALSE)
#
attach(CCES2010)
#
libcon.Obama <- CC334C
libcon.Obama[libcon.Obama==8] <- NA
libcon.DemParty <- CC334D
libcon.DemParty[libcon.DemParty==8] <- NA
libcon.RepParty <- CC334E
libcon.RepParty[libcon.RepParty==8] <- NA
libcon.TeaParty <- CC334M
libcon.TeaParty[libcon.TeaParty==8] <- NA
#
cutoff <- 3   # This section removes respondents who placed less than 3 stimuli prior to sampling.
#
AAA <- cbind(libcon.Obama,libcon.DemParty,libcon.RepParty,libcon.TeaParty)
CCES2010.2 <- CCES2010[rowSums(!is.na(AAA)) >= cutoff,]
#
detach(CCES2010)
attach(CCES2010.2)
#
#
state <- V206
congdist <- V276
#
stateabbrev <- state
stateabbrev[stateabbrev==1] <- "AL"
stateabbrev[stateabbrev==2] <- "AK"
stateabbrev[stateabbrev==4] <- "AZ"
stateabbrev[stateabbrev==5] <- "AR"
stateabbrev[stateabbrev==6] <- "CA"
stateabbrev[stateabbrev==8] <- "CO"
stateabbrev[stateabbrev==9] <- "CT"
stateabbrev[stateabbrev==10] <- "DE"
stateabbrev[stateabbrev==12] <- "FL"
stateabbrev[stateabbrev==13] <- "GA"
stateabbrev[stateabbrev==15] <- "HI"
stateabbrev[stateabbrev==16] <- "ID"
stateabbrev[stateabbrev==17] <- "IL"
stateabbrev[stateabbrev==18] <- "IN"
stateabbrev[stateabbrev==19] <- "IA"
stateabbrev[stateabbrev==20] <- "KS"
stateabbrev[stateabbrev==21] <- "KY"
stateabbrev[stateabbrev==22] <- "LA"
stateabbrev[stateabbrev==23] <- "ME"
stateabbrev[stateabbrev==24] <- "MD"
stateabbrev[stateabbrev==25] <- "MA"
stateabbrev[stateabbrev==26] <- "MI"
stateabbrev[stateabbrev==27] <- "MN"
stateabbrev[stateabbrev==28] <- "MS"
stateabbrev[stateabbrev==29] <- "MO"
stateabbrev[stateabbrev==30] <- "MT"
stateabbrev[stateabbrev==31] <- "NE"
stateabbrev[stateabbrev==32] <- "NV"
stateabbrev[stateabbrev==33] <- "NH"
stateabbrev[stateabbrev==34] <- "NJ"
stateabbrev[stateabbrev==35] <- "NM"
stateabbrev[stateabbrev==36] <- "NY"
stateabbrev[stateabbrev==37] <- "NC"
stateabbrev[stateabbrev==38] <- "ND"
stateabbrev[stateabbrev==39] <- "OH"
stateabbrev[stateabbrev==40] <- "OK"
stateabbrev[stateabbrev==41] <- "OR"
stateabbrev[stateabbrev==42] <- "PA"
stateabbrev[stateabbrev==44] <- "RI"
stateabbrev[stateabbrev==45] <- "SC"
stateabbrev[stateabbrev==46] <- "SD"
stateabbrev[stateabbrev==47] <- "TN"
stateabbrev[stateabbrev==48] <- "TX"
stateabbrev[stateabbrev==49] <- "UT"
stateabbrev[stateabbrev==50] <- "VT"
stateabbrev[stateabbrev==51] <- "VA"
stateabbrev[stateabbrev==53] <- "WA"
stateabbrev[stateabbrev==54] <- "WV"
stateabbrev[stateabbrev==55] <- "WI"
stateabbrev[stateabbrev==56] <- "WY"
#
state.abbreviations <- c("AL","AK","AZ","AR","CA","CO","CT","DE","FL","GA","HI",
	"ID","IL","IN","IA","KS","KY","LA","ME","MD","MA","MI","MN","MS","MO","MT","NE",
	"NV","NH","NJ","NM","NY","NC","ND","OH","OK","OR","PA","RI","SC","SD","TN","TX",
	"UT","VT","VA","WA","WV","WI","WY")
#
CCES2010.3 <- NULL
for (i in state.abbreviations){
temp <- CCES2010.2[stateabbrev==i,]
temp <- temp[(sample(nrow(temp), 200, replace = TRUE)),]
CCES2010.3 <- rbind(CCES2010.3,temp)
}
#
detach(CCES2010.2)
attach(CCES2010.3)
#
#
state <- V206
congdist <- V276
#
stateabbrev <- state
stateabbrev[stateabbrev==1] <- "AL"
stateabbrev[stateabbrev==2] <- "AK"
stateabbrev[stateabbrev==4] <- "AZ"
stateabbrev[stateabbrev==5] <- "AR"
stateabbrev[stateabbrev==6] <- "CA"
stateabbrev[stateabbrev==8] <- "CO"
stateabbrev[stateabbrev==9] <- "CT"
stateabbrev[stateabbrev==10] <- "DE"
stateabbrev[stateabbrev==12] <- "FL"
stateabbrev[stateabbrev==13] <- "GA"
stateabbrev[stateabbrev==15] <- "HI"
stateabbrev[stateabbrev==16] <- "ID"
stateabbrev[stateabbrev==17] <- "IL"
stateabbrev[stateabbrev==18] <- "IN"
stateabbrev[stateabbrev==19] <- "IA"
stateabbrev[stateabbrev==20] <- "KS"
stateabbrev[stateabbrev==21] <- "KY"
stateabbrev[stateabbrev==22] <- "LA"
stateabbrev[stateabbrev==23] <- "ME"
stateabbrev[stateabbrev==24] <- "MD"
stateabbrev[stateabbrev==25] <- "MA"
stateabbrev[stateabbrev==26] <- "MI"
stateabbrev[stateabbrev==27] <- "MN"
stateabbrev[stateabbrev==28] <- "MS"
stateabbrev[stateabbrev==29] <- "MO"
stateabbrev[stateabbrev==30] <- "MT"
stateabbrev[stateabbrev==31] <- "NE"
stateabbrev[stateabbrev==32] <- "NV"
stateabbrev[stateabbrev==33] <- "NH"
stateabbrev[stateabbrev==34] <- "NJ"
stateabbrev[stateabbrev==35] <- "NM"
stateabbrev[stateabbrev==36] <- "NY"
stateabbrev[stateabbrev==37] <- "NC"
stateabbrev[stateabbrev==38] <- "ND"
stateabbrev[stateabbrev==39] <- "OH"
stateabbrev[stateabbrev==40] <- "OK"
stateabbrev[stateabbrev==41] <- "OR"
stateabbrev[stateabbrev==42] <- "PA"
stateabbrev[stateabbrev==44] <- "RI"
stateabbrev[stateabbrev==45] <- "SC"
stateabbrev[stateabbrev==46] <- "SD"
stateabbrev[stateabbrev==47] <- "TN"
stateabbrev[stateabbrev==48] <- "TX"
stateabbrev[stateabbrev==49] <- "UT"
stateabbrev[stateabbrev==50] <- "VT"
stateabbrev[stateabbrev==51] <- "VA"
stateabbrev[stateabbrev==53] <- "WA"
stateabbrev[stateabbrev==54] <- "WV"
stateabbrev[stateabbrev==55] <- "WI"
stateabbrev[stateabbrev==56] <- "WY"
# 
libcon.Dem.senatecandidate <- CC334H
libcon.Dem.senatecandidate[libcon.Dem.senatecandidate==8] <- NA
libcon.Rep.senatecandidate <- CC334I
libcon.Rep.senatecandidate[libcon.Rep.senatecandidate==8] <- NA
libcon.senateincumbent1 <- CC334F
libcon.senateincumbent1[libcon.senateincumbent1==8] <- NA
libcon.senateincumbent2 <- CC334G
libcon.senateincumbent2[libcon.senateincumbent2==8] <- NA
libcon.Dem.senatecandidateNY2 <- CC334Hb
libcon.Dem.senatecandidateNY2[libcon.Dem.senatecandidateNY2==8] <- NA
libcon.Rep.senatecandidateNY2 <- CC334Ib
libcon.Rep.senatecandidateNY2[libcon.Rep.senatecandidateNY2==8] <- NA
libcon.Governorraw <- CC334B
libcon.Governorraw[libcon.Governorraw==8] <- NA
#
xx <- model.matrix(~stateabbrev - 1)
libcon.Dem.senate <- xx*libcon.Dem.senatecandidate
libcon.Dem.senate[libcon.Dem.senate==0] <- NA
libcon.Dem.senate <- libcon.Dem.senate[,colSums(is.na(libcon.Dem.senate)) < nrow(libcon.Dem.senate)]
colnames(libcon.Dem.senate) <- paste(colnames(libcon.Dem.senate), "D", sep=".")
colnames(libcon.Dem.senate) <- gsub("stateabbrev", "", colnames(libcon.Dem.senate))
#
libcon.Rep.senate <- xx*libcon.Rep.senatecandidate
libcon.Rep.senate[libcon.Rep.senate==0] <- NA
libcon.Rep.senate <- libcon.Rep.senate[,colSums(is.na(libcon.Rep.senate)) < nrow(libcon.Rep.senate)]
colnames(libcon.Rep.senate) <- paste(colnames(libcon.Rep.senate), "R", sep=".")
colnames(libcon.Rep.senate) <- gsub("stateabbrev", "", colnames(libcon.Rep.senate))
#
libcon.inc1.senate <- xx*libcon.senateincumbent1
libcon.inc1.senate[libcon.inc1.senate==0] <- NA
libcon.inc1.senate <- libcon.inc1.senate[,colSums(is.na(libcon.inc1.senate)) < nrow(libcon.inc1.senate)]
colnames(libcon.inc1.senate) <- paste(colnames(libcon.inc1.senate), "inc1", sep=".")
colnames(libcon.inc1.senate) <- gsub("stateabbrev", "", colnames(libcon.inc1.senate))
#
libcon.inc2.senate <- xx*libcon.senateincumbent2
libcon.inc2.senate[libcon.inc2.senate==0] <- NA
libcon.inc2.senate <- libcon.inc2.senate[,colSums(is.na(libcon.inc2.senate)) < nrow(libcon.inc2.senate)]
colnames(libcon.inc2.senate) <- paste(colnames(libcon.inc2.senate), "inc2", sep=".")
colnames(libcon.inc2.senate) <- gsub("stateabbrev", "", colnames(libcon.inc2.senate))
#
libcon.Dem.senateNY2 <- xx*libcon.Dem.senatecandidateNY2
libcon.Dem.senateNY2[libcon.Dem.senateNY2==0] <- NA
libcon.Dem.senateNY2 <- libcon.Dem.senateNY2[,colSums(is.na(libcon.Dem.senateNY2)) < nrow(libcon.Dem.senateNY2)]
libcon.Dem.senate <- cbind(libcon.Dem.senate, libcon.Dem.senateNY2)
colnames(libcon.Dem.senate)[ncol(libcon.Dem.senate)] <- "NY2.D"
#
libcon.Rep.senateNY2 <- xx*libcon.Rep.senatecandidateNY2
libcon.Rep.senateNY2[libcon.Rep.senateNY2==0] <- NA
libcon.Rep.senateNY2 <- libcon.Rep.senateNY2[,colSums(is.na(libcon.Rep.senateNY2)) < nrow(libcon.Rep.senateNY2)]
libcon.Rep.senate <- cbind(libcon.Rep.senate, libcon.Rep.senateNY2)
colnames(libcon.Rep.senate)[ncol(libcon.Rep.senate)] <- "NY2.R"
#
libcon.Governor <- xx*libcon.Governorraw
libcon.Governor[libcon.Governor==0] <- NA
libcon.Governor <- libcon.Governor[,colSums(is.na(libcon.Governor)) < nrow(libcon.Governor)]
colnames(libcon.Governor) <- paste(colnames(libcon.Governor), "gov", sep=".")
colnames(libcon.Governor) <- gsub("stateabbrev", "", colnames(libcon.Governor))
#
#
Crist.I.FL <- libcon.Governor[,"FL.gov"]
#
libcon.Senate <- cbind(libcon.Dem.senate, libcon.Rep.senate, libcon.inc1.senate, libcon.inc2.senate,Crist.I.FL)
#
#
libcon.self <- CC334A
libcon.self[libcon.self==8] <- NA
libcon.Obama <- CC334C
libcon.Obama[libcon.Obama==8] <- NA
libcon.DemParty <- CC334D
libcon.DemParty[libcon.DemParty==8] <- NA
libcon.RepParty <- CC334E
libcon.RepParty[libcon.RepParty==8] <- NA
libcon.TeaParty <- CC334M
libcon.TeaParty[libcon.TeaParty==8] <- NA
#
#
libcon.common <- cbind(libcon.self, libcon.Obama, libcon.DemParty, libcon.RepParty, libcon.TeaParty)
#
partyid <- V212d
general_voter <- vote_gen10  ### Validated
primary_voter <- vote_pri10  ### Validated
#
#
T <- cbind(libcon.common, libcon.Senate)
#
T <- T - 4
#
self <- T[,1]
TT <- T[,-1]
#
N <- nrow(TT)
q <- ncol(TT)
z <- TT
#
alpha.starts <- runif(N,-2,2)
beta.starts <- runif(N,-2,2)
#
zhat.starts <- runif(q,-1,1)
zhat.starts[1] <- -1 + runif(1,-0.099,0.099)
zhat.starts[2] <- -1 + runif(1,-0.099,0.099)
zhat.starts[3] <- 1 + runif(1,-0.099,0.099)
zhat.starts[4] <- 1 + runif(1,-0.099,0.099)
#
inits <- function() {list (zhat=zhat.starts, a=alpha.starts, b=beta.starts)}
#
BAM.sim <- jags.model('c:/Dropbox/BAM/CCES_2010/BAM_JAGScode_CCES2010.bug',
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
save(samps,file="c:/Dropbox/BAM/CCES_2010/CCES2010_Senate_samps.Rda")
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
#
# Load results from disk:
load("c:/Dropbox/BAM/CCES_2010/CCES2010_Senate_samps.Rda") 
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
#
#   STIMULI ESTIMATES
#
libcon_zhat <- summary(samps[,grep("zhat", colnames(samps[[1]]))])
libcon_sigmaj.sq <- summary(samps[,grep("sigmaj.sq", colnames(samps[[1]]))])
#
#
#   Catterpillar plot of 2010 Senate Candidates with uncertainty bounds:
#
bayes.out <- cbind(libcon_zhat[[1]][,1],libcon_zhat[[2]][,1],libcon_zhat[[2]][,5])
colnames(bayes.out) <- c("Mean","LB","UB")
senate.results <- read.csv("c:/Dropbox/BAM/CCES_2010/Senate_results.csv", header=TRUE)
rownames(bayes.out) <- senate.results$LABEL2
# 
bayes.out2 <- bayes.out[order(bayes.out[,1]), ]
#
pdf('c:/Dropbox/BAM/CCES_2010/CCES2010_caterpillar.pdf', width=10, height=10)
rx1 <- dotplot(bayes.out2[1:70,1],cex=0.0, xlim=c(-1.5, 1.5),
	xlab=list(label = "Liberal - Conservative", cex=1.1),
	main="Stimuli Estimates with 95% Credible Intervals",
	scales=list(cex=0.65, font=2),
	points = list(pch = c(16, 4), col = "black", cex = 1),
	panel = function(x,y,subscripts,lower, upper){
	panel.dotplot(x, y,col = "white", lty = 1)
	panel.points(x,y, pch=16, col="black", cex=.7)
 	panel.segments(bayes.out2[1:70,2], y, bayes.out2[1:70,3], y, lwd=2)})

#
rx2 <- dotplot(bayes.out2[71:141,1],cex=0.0, xlim=c(-1.5, 1.5),
	xlab=list(label = "Liberal - Conservative", cex=1.1),
	main="Stimuli Estimates with 95% Credible Intervals",
	scales=list(cex=0.65, font=2),
	points = list(pch = c(16, 4), col = "black", cex = 1),
	panel = function(x,y,subscripts,lower, upper){
	panel.dotplot(x, y,col = "white", lty = 1)
	panel.points(x,y, pch=16, col="black", cex=.7)
 	panel.segments(bayes.out2[71:141,2], y, bayes.out2[71:141,3], y, lwd=2)})
#
#
print(rx1,position=c(0.0,0.0,0.5,1),more=TRUE)
print(rx2,position=c(0.5,0.0,1,1),more=TRUE)
#
dev.off()
#
#
#   IDEAL POINTS
#
libcon_self <- self
#
libcon_a <- samps[,grep(pattern="a[", colnames(samps[[1]]), fixed=TRUE),]
libcon_b <- samps[,grep(pattern="b[", colnames(samps[[1]]), fixed=TRUE),]
libcon_all.a <- do.call("rbind", libcon_a)
libcon_all.b <- do.call("rbind", libcon_b)
#
libcon_nsamp <- nrow(libcon_all.a)
libcon_nresp <- ncol(libcon_all.a)
#
libcon_idealpt <- rep(NA, libcon_nsamp*libcon_nresp)
dim(libcon_idealpt) <- c(libcon_nsamp,libcon_nresp)
for (i in 1:libcon_nresp){
for (j in 1:libcon_nsamp){
libcon_idealpt[j,i] <- ((libcon_self[i] - libcon_all.a[j,i]) / libcon_all.b[j,i])
}}
#
libcon_idealpt.percentiles <- rep(NA,libcon_nresp*3)
dim(libcon_idealpt.percentiles) <- c(libcon_nresp,3)
colnames(libcon_idealpt.percentiles) <- c("5%","50%","95%")
for (i in 1:libcon_nresp){
libcon_idealpt.percentiles[i,] <- quantile(libcon_idealpt[,i],probs=c(0.05,0.5,0.95),na.rm=TRUE)
}
libcon_idealpt.estimates <- libcon_idealpt.percentiles[,2]
#
#
#
# tapply(avalues, stateabbrev, mean)
#
avalues <- summary(libcon_a)[[1]][,1]
#
alpha.libcon.plot <- bwplot(avalues ~ stateabbrev, horizontal=FALSE,
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
# =================
#
library(apsrtable)
library(QuantPsyc)
#
alpha.ideology.ols.1 <- lm(avalues ~ libcon_self)
alpha.ideology.ols.2 <- lm(avalues ~ average)
apsrtable(alpha.ideology.ols.1,alpha.ideology.ols.2)
#
lm.beta(alpha.ideology.ols.1)
lm.beta(alpha.ideology.ols.2)
#
# =================
#
alpha.libcon.plot2 <- bwplot(avalues ~ average, horizontal=FALSE,
	main=list(label="Ideological Self-Placement and Alpha (Shift) Term",cex=1.4),
	xlab=list(label="Self-Placement (Very Liberal to Very Conservative)",cex=1.3),
	ylab=list(label="Alpha Posterior Mean",cex=1.3),
	scales=list(cex=1.2)
	)
#
update(alpha.libcon.plot2, panel = function(...) {
            panel.abline(h=0, lwd=2, lty=1, col="light grey")
            panel.bwplot(...)
        })
#
#
#  BETA VALUES
#
bvalues <- summary(libcon_b)[[1]][,1]
#
bvalues.pos <- bvalues
bvalues.pos[bvalues.pos < 0] <- 0
bvalues.pos[bvalues.pos > 0] <- 1
#
table(bvalues.pos,vote_gen10)
table(bvalues.pos,vote_pri10)
#
#
#
#
library(apsrtable)
#
senate.results <- read.csv("c:/Dropbox/BAM/CCES_2010/Senate_results.csv", header=TRUE)
#
result.all <- lm(senate.results$MEAN ~ senate.results$DWNOM.COMMON + senate.results$CFSCORE)
#
party.code <- senate.results$PARTY
#
result.democrats <- lm(senate.results$MEAN[party.code=="D" | party.code=="I"] ~ senate.results$DWNOM.COMMON[party.code=="D" | party.code=="I"] + 
	senate.results$CFSCORE[party.code=="D" | party.code=="I"])
#
result.republicans <- lm(senate.results$MEAN[party.code=="R"] ~ senate.results$DWNOM.COMMON[party.code=="R"] + 
	senate.results$CFSCORE[party.code=="R"])
#
apsrtable(result.all,result.democrats,result.republicans)
#
#
#
#
library(lattice)
#
test <- data.frame(cbind(senate.results$DWNOM.COMMON,senate.results$CFSCORE,senate.results$MEAN,senate.results$NAIVEMEAN))
colnames(test) <- c("DW-NOMINATE","CF","BAM","MEAN")
test <- test[senate.results$INCORCOMPRACE==1,]
test.parties <- senate.results$PARTY
test.parties <- test.parties[senate.results$INCORCOMPRACE==1]
#
pdf('c:/Dropbox/BAM/CCES_2010/CCES2010_correlations_bw.pdf', width=8, height=8)
panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))

  #  Overall Correlations
  r <- cor(x, y, use="complete")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r = ", txt, sep = "")
  text(0.5, 0.7, txt, col="black", cex=1.7, font=2)
  #  Party Correlations
  r.dem <- cor(x[test.parties!="R"], y[test.parties!="R"], use="complete")
  txt2 <- format(c(r.dem, 0.123456789), digits = digits)[1]
  txt2 <- paste("r (Dem) = ", txt2, sep = "")
  text(0.5, 0.5, txt2, col="gray60", cex=1.7, font=2)
  r.rep <- cor(x[test.parties=="R"], y[test.parties=="R"], use="complete")
  txt3 <- format(c(r.rep, 0.123456789), digits = digits)[1]
  txt3 <- paste("r (Rep) = ", txt3, sep = "")
  text(0.5, 0.3, txt3, col="gray20", cex=1.7, font=2)

}
pairs(test, cex.labels=1.8, labels=c("DW-NOMINATE\nScore","CF Score","BAM Score","Mean\nPlacement"), upper.panel = panel.cor,
lower.panel = function(x,y){
abline(v=seq(-2,2,0.5), lty=1, col='gray80')
abline(h=seq(-2,2,0.5), lty=1, col='gray80')
points(x[test.parties!="R"],y[test.parties!="R"], cex=1.2, pch="D", col="gray60") 
points(x[test.parties=="R"],y[test.parties=="R"], cex=1.2, pch="R", col="gray20")
lines(lowess(y[!is.na(x) & !is.na(y)] ~ x[!is.na(x) & !is.na(y)]), col="black", lwd=2)
}) 
#
dev.off()
#
#
#
#
#  MAPS
#
BAM <- read.csv("c:/Dropbox/BAM/CCES_2010/state_alpha_values.csv", header=TRUE)
head(BAM)
#
fixup <- function(usa,alaskaFix,hawaiiFix){
  alaska=usa[usa$STATE_NAME=="Alaska",]
  alaska = fix1(alaska,alaskaFix)
  proj4string(alaska) <- proj4string(usa)
  hawaii = usa[usa$STATE_NAME=="Hawaii",]
  hawaii = fix1(hawaii,hawaiiFix)
  proj4string(hawaii) <- proj4string(usa)
  usa = usa[! usa$STATE_NAME %in% c("Alaska","Hawaii"),]
  usa = rbind(usa,alaska,hawaii)
  return(usa)
}
#
fix1 <- function(object,params){
  r=params[1];scale=params[2];shift=params[3:4]
  object = elide(object,rotate=r)
  size = max(apply(bbox(object),1,diff))/scale
  object = elide(object,scale=size)
  object = elide(object,shift=shift)
  object
}
us <- readOGR(dsn = "states_21basic",layer="states")
usAEA <- spTransform(us,CRS("+init=epsg:2163"))
state <- fixup(usAEA,c(-35,1.5,-2800000,-2600000),c(-35,1,6800000,-1600000))
#
state@data <- cbind(state@data, BAM[match(state@data$STATE_ABBR, BAM[,1]), 3:5])
#
rwramp <- colorRampPalette(c("red", "white"), space="Lab")
wbramp <- colorRampPalette(c("white", "blue"), space="Lab")
#
seq.range <- function(x, n){
	x <- na.omit(x)
	seq(min(x), max(x), length=n)
}
#
ats <- c(seq(min(state@data$alphamean, na.rm=TRUE), 0, length=10), seq(0, max(state@data$alphamean, na.rm=TRUE), length=3)[-1])
ats <- c(-.5, -.3, -.25, -.2, -.15, -.1, -.05, 0, .05, .1, .15, .2, .35)+.025
png("alpha_map1.png", width=4, height=2.5, units="in", res=600)
spplot(state, "alphamean", col.regions=c(rwramp(7), wbramp(6)[-1]), at=ats)
dev.off()
#
ramp <- colorRampPalette(c("red", "white", "blue"), space="Lab")
ats <- quantile(state@data$alphamean, probs=seq(0, 1, by=.05), na.rm=TRUE) 
ats <- ats + c(-.1, rep(0, length(ats)-2), .1)
png("alpha_map2.png", width=4, height=2.5, units="in", res=600)
spplot(state, "alphamean", col.regions=ramp(length(ats)-1), at=ats)
dev.off()
#
ats <- seq.range(state@data$alphamean, 21)
png("alpha_map3.png", width=4, height=2.5, units="in", res=600)
spplot(state, "alphamean", col.regions=ramp(length(ats)-1), at=ats)
dev.off()
#
#
state@data$preds <- predict(lm(alphamean ~ meanobamaplace, data=state@data, na.action=na.exclude))
ats <- c(-.5, -.3, -.25, -.2, -.15, -.1, -.05, 0, .05, .1, .15, .2, .35)+.025
png("obamaplace_map_a1.png", width=4, height=2.5, units="in", res=600)
spplot(state, "preds", col.regions=ramp(length(ats)-1), at=ats)
dev.off()
#
ats <- quantile(state@data$alphamean, probs=seq(0, 1, by=.05), na.rm=TRUE)
png("obamaplace_map_a2.png", width=4, height=2.5, units="in", res=600)
spplot(state, "preds", col.regions=ramp(length(ats)-1), at=ats)
dev.off()
#
ats <- seq.range(state@data$alphamean, 21)
png("obamaplace_map_a3.png", width=4, height=2.5, units="in", res=600)
spplot(state, "preds", col.regions=ramp(length(ats)-1), at=ats)
dev.off()
#
#
state@data$preds08 <- predict(lm(alphamean ~ obama08, data=state@data, na.action=na.exclude))
#
ats <- c(-.5, -.3, -.25, -.2, -.15, -.1, -.05, 0, .05, .1, .15, .2, .35)+.025
png("obama08_map_a1.png", width=4, height=2.5, units="in", res=600)
spplot(state, "preds08", col.regions=ramp(length(ats)-1), at=ats)
dev.off()
#
ats <- quantile(state@data$alphamean, probs=seq(0, 1, by=.05), na.rm=TRUE)
png("obama08_map_a2.png", width=4, height=2.5, units="in", res=600)
spplot(state, "preds08", col.regions=ramp(length(ats)-1), at=ats)
dev.off()
#
ats <- seq.range(state@data$alphamean, 21)
png("obama08_map_a3.png", width=4, height=2.5, units="in", res=600)
spplot(state, "preds08", col.regions=ramp(length(ats)-1), at=ats)
dev.off()
#
#
state@data$meanobamaplace4 <- state@data$meanobamaplace+4
ats <- quantile(state@data$meanobamaplace4, probs=seq(0, 1, by=.05), na.rm=TRUE)
ats <- ats + c(-.1, rep(0, length(ats)-2), .1)
png("obamaplace_map.png", width=4, height=2.5, units="in", res=600)
spplot(state, "meanobamaplace4", col.regions=ramp(length(ats)-1), at=ats)
dev.off()
#
ats <- quantile(state@data$obama08, probs=seq(0, 1, by=.05), na.rm=TRUE)
png("obama08_map.png", width=4, height=2.5, units="in", res=600)
spplot(state, "obama08", col.regions=ramp(length(ats)-1), at=ats)
dev.off()
#
ats <- seq.range(state@data$meanobamaplace, 21)
png("obamaplace_map_even.png", width=4, height=2.5, units="in", res=600)
spplot(state, "meanobamaplace", col.regions=ramp(length(ats)-1), at=ats)
dev.off()
#
ats <- seq.range(state@data$obama08, 21)
png("obama08_map_even.png", width=4, height=2.5, units="in", res=600)
spplot(state, "obama08", col.regions=ramp(length(ats)-1), at=ats)
dev.off()
#
#
#
#
all_states <- map_data("state")
all_states
head(all_states)
BAM$region <- BAM$statename
Total <- merge(all_states, BAM, by="region")
head(Total)
#
Total$meanobamaplace <- Total$meanobamaplace + 4
#
p <- ggplot()
p <- p + geom_polygon(data=Total, aes(x=long, y=lat, group = group, fill=Total$meanobamaplace), color="white") + 
scale_fill_continuous(low = "red", high = "blue", guide="colorbar")
P1 <- p + theme_bw() + labs(fill = "Mean Obama Placement\n on Liberal-Conservative Scale", title = "Mean Obama Placement, 2010", x="", y="")
P1 + scale_y_continuous(breaks=c()) + scale_x_continuous(breaks=c()) + theme(panel.border = element_blank())
#
p <- ggplot()
p <- p + geom_polygon(data=Total, aes(x=long, y=lat, group = group, fill=Total$alphamean), color="white") + 
scale_fill_continuous(low = "red", high = "blue", guide="colorbar")
P1 <- p + theme_bw() + labs(fill = "Mean Alpha Value" ,title = "Mean Alpha Value", x="", y="")
P1 + scale_y_continuous(breaks=c()) + scale_x_continuous(breaks=c()) + theme(panel.border = element_blank())
#