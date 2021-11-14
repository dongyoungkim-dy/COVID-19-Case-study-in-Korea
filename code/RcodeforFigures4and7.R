---
  title: "Technical Vignette to Accompany `Towards reduction in bias in epidemic curves...'"
author: "Burstyn I, Goldstein ND, Gustafson P"
date: "03/04/2020"
---
  
  * There are five public files: This .Rmd document, two versions of its pdf output (Alberta and Philadelphia), and two csv files (three-columns: date, positive, negative) for Alberta and Philadelphia. respectively.

* Note the "switch" below to produce either AB or KR output.

* Note the "switch" below to echo the R commands or not in the pdf output.

* Note the "switch" below to additionally produce separate pdf files for manuscript figures. 

 

setwd("D:/01. 대학원지원/2021 캐나다 대학원 지원/UBC/corona19")


## which jurisdiction

JRSDCT <- "KR"


require(rjags)
require(MCMCvis)

### bespoke function to add CIs to plot
add.ci <- function(x, est, se, lvl=.95, trnc.lo=0) {
  for (i in 1:length(x)) {
    points(rep(x[i],2),
           pmax(est[i]+qnorm(1-(1-lvl)/2)*c(-1,1)*se[i],rep(trnc.lo,2)),
           type="l")
  }
}


if (JRSDCT=="KR") {
  dta <- read.csv("south korea march 1 brief.csv", header=T)
}

T.end <- dim(dta)[1]

ystr <- dta$positive[1:T.end]
n <- dta$negative[1:T.end] + dta$positive[1:T.end]
q.hat <- ystr/n
se <- sqrt(q.hat*(1-q.hat)/n)


### first knot must be at day 1
### last knot must be at day T.end

### choose three-equally spaced interior knots
knts <- c(1, round(T.end*(1:3)/4), T.end)
num.kn <- length(knts)
spc.kn <- knts[-1] - knts[-num.kn]


Proportion positive ($Y^{*}/n$ time-series), with 95% confidence intervals, and knot indicators


plot(1:T.end, q.hat,
     xlab=paste("Day (",dta$date[1]," - ",dta$date[T.end],")",sep=""),
     ylab="Proportion Positive", main="Korea",ylim=c(0, max(q.hat+2.1*se)))
add.ci((1:T.end), q.hat, se)
points(knts, rep(0,num.kn),pch=17,col="red")


### JAGS specification
genmod.JAGS <- "model{
  
  ### prior distribution
  
  ### how much weight on the linear component for sens
  ### (earlier version corresponds to no weight)
  
  ### sn.wt <- 0 no smoothness (as per earlier versions)
  ### sn.wt <- 1 perfectly linear
  sn.wt ~ dunif(0.5,0.9)  ### some unspecified amount of smoothing
  
  ### linear component endpoints
  sn.LL ~ dunif(sn.wt*sn.lo[1], sn.wt*sn.hi[1])
  sn.LR ~ dunif(sn.wt*sn.lo[num.kn], sn.wt*sn.hi[num.kn])
    
  ### prev, sn piecewise linear
  ### parameterized by value at knots
  for (i in 1:num.kn) {
    r.kn[i] ~ dunif(0, r.hi[i])
    
    ### linear component
    sn.L[i] <- ((knts[i]-knts[1])*sn.LR+(knts[num.kn]-knts[i])*sn.LL)/
               (knts[num.kn]-knts[1])
    
    ### jumpy component
    sn.J[i] ~ dunif((1-sn.wt)*sn.lo[i],(1-sn.wt)*sn.hi[i]) 
    
    ### two together
    sn.kn[i] ~ dsum(sn.L[i], sn.J[i])
  }
  
  sp ~ dunif(sp.lo,1)
  
  ### these imply the daily values
  for (i in 1:(num.kn-1)) {
    for (j in 0:(spc.kn[i]-1)) {
      r[knts[i]+j] <- ((spc.kn[i]-j)*r.kn[i]+j*r.kn[i+1])/(spc.kn[i])
    }  
  }    
  r[knts[num.kn]] <- r.kn[num.kn]
  
  ### these imply the daily values
  for (i in 1:(num.kn-1)) {
    for (j in 0:(spc.kn[i]-1)) {
      sn[knts[i]+j] <- ((spc.kn[i]-j)*sn.kn[i]+j*sn.kn[i+1])/(spc.kn[i])
    }  
  }    
  sn[knts[num.kn]] <- sn.kn[num.kn]
  
  for (i in 1:(knts[num.kn])) {
    y[i] ~ dbinom(r[i], n[i])            ### true positives
    ystr1[i] ~ dbinom(sn[i], y[i])       ### correct positives
    ystr2[i] ~ dbinom(1-sp, n[i]-y[i])   ### false positives
    ystr[i] ~ sum(ystr1[i], ystr2[i])
  }

}"


if (JRSDCT=="KR") {
  
  ### specify hyperparameters: 
  ### important this is inside this (cached) code chunk, 
  ### so a change in hyperparams triggers MCMC re-run)
  ### also one such chunk per jurisdiction, for same reason
  
  sn.lo.KR <- rep(0.4, num.kn)  # same bound at each knot
  sn.hi.KR <- rep(0.9, num.kn)  # same upper-bound at each knot
  r.hi.KR <- rep(0.15, num.kn)   # same upper-bound at each knot
  sp.lo.KR <- 0.95
  
  ### generative model, data go in
  
  mod <- jags.model(textConnection(genmod.JAGS),
                    data=list(knts=knts, num.kn=num.kn, spc.kn=spc.kn,
                              sp.lo=sp.lo.KR, sn.lo=sn.lo.KR, 
                              sn.hi=sn.hi.KR, r.hi=r.hi.KR,
                              ystr=ystr, n=n),
                    inits=list(y=round(1.2*ystr), ystr1=ystr, 
                               ystr2=rep(0,length(ystr))),
                    n.chains=4)          
  
  ###  MC output comes out
  
  opt.JAGS.KR <- coda.samples(mod, 
                              variable.names=c("sp","sn.wt","sn.kn","r.kn","y"),
                              n.iter=400000, n.thin=200)




if (JRSDCT=="KR") {
  opt.JAGS <- opt.JAGS.KR
  sn.lo <- sn.lo.KR; sn.hi <- sn.hi.KR
  r.hi <- r.hi.KR; sp.lo <- sp.lo.KR
}


###Hyperparameter settings


r.hi
sn.lo
sn.hi
sp.lo



### raw posterior draws for future use
mc.opt <-as.matrix(opt.JAGS)
mc.qnt <- summary(opt.JAGS)$quantiles
n.draws <- dim(mc.opt)[1]



###Adjusted daily positive tests


ndx <- (1:(dim(mc.qnt)[1]))[row.names(mc.qnt)=="y[1]"]
plot(1:T.end, mc.qnt[ndx:(ndx+T.end-1),"50%"],
     pch=19,
     main="Korea",
     ylim=c(0, max(mc.qnt[ndx:(ndx+T.end-1),"97.5%"])),
     ylab="Positive Tests",
     xlab=paste("Day (",dta$date[1]," - ",dta$date[T.end],")",sep="")
)
for (i in 1:T.end) {
  points(rep(i,2), mc.qnt[ndx+i-1,c("2.5%","97.5%")], type="l")
}

points(1:(knts[num.kn]), ystr, col="blue", pch=5)


##Specificity in testing pool (posterior distribution)


mc.qnt["sp",]


hist(mc.opt[,"sp"], breaks=seq(from=sp.lo,to=1,by=.0025), prob=T,
     xlab="Specificity",ylab="Density",main="Korea") 




###Sensitivity in testing pool (posterior mean and draws)

ndx <- sample(1:n.draws, size=75, replace=T)

plot(0,0, type="n", xlim=c(0,T.end), ylim=c(0.35,1),
     xlab=paste("Day (",dta$date[1]," - ",dta$date[T.end],")",sep=""), 
     ylab="Sensitivity", main="Korea",)

points(knts, sn.lo, type="l", lty=3)
points(knts, sn.hi, type="l", lty=3)

col.start.sn <- (1:(dim(mc.opt)[2]))[colnames(mc.opt)=="sn.kn[1]"]
for (i in ndx) {
  points(knts, mc.opt[i,col.start.sn:(col.start.sn+num.kn-1)],
         type="l", col="grey")
}         

points(knts,
       apply(mc.opt[,col.start.sn:(col.start.sn+num.kn-1)], 2, mean),
       lwd=2, col="red",type="l")


#Prevalence in testing pool (posterior mean and draws)


ndx <- sample(1:n.draws, size=50, replace=T)

plot(0,0, type="n", xlim=c(0,T.end), ylim=c(0,max(r.hi)),
     xlab=paste("Day (",dta$date[1]," - ",dta$date[T.end],")",sep=""), 
     ylab="Prevalence", main="Korea",)
points(knts, r.hi, type="l", lty=3)

col.start.r <- (1:(dim(mc.opt)[2]))[colnames(mc.opt)=="r.kn[1]"]
for (i in ndx) {
  points(knts, mc.opt[i,col.start.r:(col.start.r+num.kn-1)],
         type="l", col="grey")
}         

points(knts,
       apply(mc.opt[,col.start.r:(col.start.r+num.kn-1)], 2, mean),
       lwd=2, col="red",type="l")
```


$E(Y^{*}/n)=r\times sens + (1-r)*(1-spec)$ (posterior mean and draws)

r.kn <- mc.opt[,col.start.r:(col.start.r+num.kn-1)]
r.full <- matrix(NA,dim(mc.opt)[1],T.end)

for (i in 1:(num.kn-1)) {
  for (j in 0:(spc.kn[i]-1)) {
    r.full[,knts[i]+j]<- (1/spc.kn[i])*
      ( (spc.kn[i]-j)*r.kn[,i]+j*r.kn[,i+1] )
  }  
}    
r.full[,knts[num.kn]] <- r.kn[,num.kn]

sn.kn <- mc.opt[,col.start.sn:(col.start.sn+num.kn-1)]
sn.full <- matrix(NA,dim(mc.opt)[1],T.end)

for (i in 1:(num.kn-1)) {
  for (j in 0:(spc.kn[i]-1)) {
    sn.full[,knts[i]+j]<- (1/spc.kn[i])*
      ( (spc.kn[i]-j)*sn.kn[,i]+j*sn.kn[,i+1] )
  }  
}    
sn.full[,knts[num.kn]] <- sn.kn[,num.kn]

ndx <- sample(1:n.draws, size=50, replace=T)

plot(0,0, type="n", xlim=c(0,T.end), ylim=c(0,max(q.hat+2.1*se)),
     xlab=paste("Day (",dta$date[1]," - ",dta$date[T.end],")",sep=""), 
     ylab="Expected Count", main="Korea",)

for (i in ndx) {
  points(1:T.end, r.full[i,]*sn.full[i,]+ (1-r.full[i,])*(1-mc.opt[i,"sp"]),
         type="l", col="grey")
}         

points(1:T.end,
       apply(r.full*sn.full + (1-r.full)*(1-mc.opt[,"sp"]), 2, mean),
       lwd=2, col="red",type="l")

points(1:T.end, q.hat, ylim=c(0,1))
add.ci((1:T.end), q.hat, se)

##Some clues about MCMC numerical performance

```{r}
MCMCtrace(opt.JAGS, params="sp", pdf=F)

MCMCtrace(opt.JAGS, params="sn.wt", pdf=F)

MCMCsummary(opt.JAGS, params="sp")

MCMCsummary(opt.JAGS, params="sn.kn")

MCMCsummary(opt.JAGS, params="r.kn")

MCMCsummary(opt.JAGS, params="y")
```






