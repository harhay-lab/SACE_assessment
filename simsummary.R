library(doParallel)
source("simfun.R")
load("ssim.Rdata")
load("fobj.Rdata")
load("bobj.RData")

nsim=length(bobj$bio)
tsace=ssim[[3]] 

cc=detectCores()
cl=makeCluster(cc)
registerDoParallel(cl)

mc=foreach(m=1:nsim,.combine='cbind',.multicombine=T,.maxcombine=nsim,
.errorhandling ="remove") %dopar%{
bmoi=bobj$bmo[[m]]
bsummary(bmoi,tsace)}

fc=foreach(m=1:nsim,.combine='cbind',.multicombine=T,.maxcombine=nsim,
.errorhandling ="remove") %dopar%{
bfoi=bobj$bfo[[m]]
bsummary(bfoi,tsace)}

ic=foreach(m=1:nsim,.combine='cbind',.multicombine=T,.maxcombine=nsim,
.errorhandling ="remove") %dopar%{
bioi=bobj$bio[[m]]
bsummary(bioi,tsace)}

qc=foreach(m=1:nsim,.combine='cbind',.multicombine=T,.maxcombine=nsim,
.errorhandling ="remove") %dopar%{
qoi=bobj$bqo[[m]]
bsummary(qoi,tsace)}

stopCluster(cl) 

rmc=rowMeans(mc)
rfc=rowMeans(fc)
ric=rowMeans(ic)
rqc=rowMeans(qc)

coverage=c(rmc[2],rfc[2],ric[2],rqc[2])
names(coverage)=c("ME","FE","IPW","WRV")
power=c(1-c(rmc[1],rfc[1],ric[1],rqc[1]),mean(fobj$wo<0.05))
names(power)=c("ME","FE","IPW","WRV","W/R")

#### computing time comparison

load("time.RData")

it=timeo[[3]][3]
ratio=c(timeo[[1]][3],timeo[[2]][3],timeo[[4]][3],timeo[[5]][3])/it
names(ratio)=c("ME","FE","WRV","W/R")
