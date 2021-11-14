library(plotrix)
#Number of cases observed C*
cstr<-c(0,0,0,0,2,0,2,34,16,74,190,210,207,130,253,449,427,909,595)

T<-length(cstr)
####################################################
#time-invariant SN (the left-hands side of Figure E)
####################################################
#observed epi curve
days<-rep(1:T)
#95% confidence interval for observed counts under Poisson distribution
#http://ms.mcmaster.ca/peter/s743/poissonalpha.html
li<-qchisq(0.025, 2*cstr)/2
ui<-qchisq(0.975, 2*(cstr+1))/2
plotCI(x=days,y=cstr,li=li,ui=ui,
       col="blue", lwd=2, ylim=c(0,2700), xlim=c(1,19),
       ylab="predicted true cases (C)",
       xlab="presumed date of onset (t: Feburary 12-March 1, 2020)")
#used https://www.desmos.com/calculator/kx83qio7yl calculator to select param of beta dist for Sn
#adjusted epi curves MC sensitivity analysis
sim=10 #number of MC simulations
#SN=0.95
for (i in c(1:sim)){
  sn=rbeta(T, 18.05, 0.95)
  y=cstr/sn
  lines(days, y, col="black", lwd=0.01)
  rm(y, sn)}
#SN=0.8
for (i in c(1:sim)){
  sn=rbeta(T, 51.2, 12.8)
  y=cstr/sn
  lines(days, y, col="lightgrey", lwd=0.01)
  rm(y, sn)}
#SN=0.6
for (i in c(1:sim)){
  sn=rbeta(T, 57.6, 38.4)
  y=cstr/sn
  lines(days, y, col="red", lwd=0.01)
  rm(y, sn)}
#SN=0.4
#for (i in c(1:sim)){
#  sn=rbeta(T, 38.4, 57.6)
#  y=cstr/sn
#  lines(days, y, col="orange", lwd=0.01)
#  rm(y, sn)}

for (i in c(1:sim)){
  sn=rbeta(T, 38.4, 57.6)
  
  y=cstr/sn
  print(y)
  lines(days, y, col="orange", lwd=0.01)
}
max(y)
max(a)
abline(h=seq(1,3300,100),col="gray", lty=1)
legend(cex=0.5, x=1, y=2300, legend=c("Observed (C*) and 95%CI", "Sn=95%", "Sn=80%", "Sn=60%", "Sn=40%"),
       col=c("blue", "lightgrey", "red", "black", "orange"), lty=c(1,1,1,1,1), lwd=c(2,1,1,1,1), bg="ivory", bty = "grey", )

####################################################
#time-varying SN (the right-hands side of Figure E)
####################################################
#observed epi curve
days<-rep(1:T)
sim=20 #number of simulation realizations plotted
#95% confidence interval for observed counts under Poisson distribution
#http://ms.mcmaster.ca/peter/s743/poissonalpha.html
li<-qchisq(0.025, 2*cstr)/2
ui<-qchisq(0.975, 2*(cstr+1))/2
plotCI(x=days,y=cstr,li=li,ui=ui,
       col="blue", lwd=2, ylim=c(0,3100), xlim=c(1,19),
       ylab="predicted true cases (C)",
       xlab="presumed date of onset (t: Feburary 12-March 1, 2020)")
#adjusted epi curves MC sensitivity analysis
#SN increases 0.40 to 0.6 to .95
for (i in c(1:sim)){
  sn1=rbeta(1+T/3, 38.4, 57.6)
  sn2=rbeta(T/3, 57.6, 38.4)
  sn3=rbeta(T/3, 18.05, 0.95)
  sn=c(sn1, sn2, sn3)
  y=cstr/sn
  lines(days, y, col="green")
  rm(y, sn)}
#SN decreases 0.95 to 0.6 to 0.4
for (i in c(1:sim)){
  sn1=rbeta(1+T/3, 18.05, 0.95)
  sn2=rbeta(T/3, 57.6, 38.4)
  sn3=rbeta(T/3, 38.4, 57.6)
  sn=c(sn1, sn2, sn3)
  y=cstr/sn
  print(y)
  lines(days, y, col="brown")
}
max(y)
abline(h=seq(1,3300,100),col="gray", lty=1)
abline(v=c(1+T/3, 2*T/3),col="black", lty=1)

legend(cex=0.5, x=1, y=2700, legend=c("Observed (C*) and 95%CI", "Sn increases from 40, 60, 95% at break-points",
                            "Sn decreases from 95, 60, 40% at break-points", "break-points"),
       col=c("blue", "green", "brown", "black"), lty=c(1,1,1,1), lwd=c(2,1,1,1), bg="ivory", bty = "grey")

#########end of code for Korea##############

