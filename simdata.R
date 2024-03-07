library(doParallel)
source("simfun.R")

nc =40 # number of clusters
tt = 2 # the total variance
trho =0.05 # ICC
taut=tt*trho # variance of the random intercept
nsim = 5 # number of simulations

# simulate the random intercepts

eta1=matrix(rnorm(nc*nsim,mean=0,sd=sqrt(taut)),nc,nsim)
eta0=matrix(rnorm(nc*nsim,mean=0,sd=sqrt(taut)),nc,nsim)

#### simulate predictors

fsize = 30
shapep=4
scalep=15/2
tass=c(1,2,1)
tasn=c(-0.5,-1.5,-1)

cs="fixed"

if(cs=="fixed") {rn = matrix(fsize,nc,nsim)} # fixed cluster size

if(cs=="variable") { rn= simn(shapep, scalep, nc, nsim)} # variable cluster size

cc=detectCores()
cl=makeCluster(cc)
registerDoParallel(cl)

sx1=foreach(kk=1:nsim,.combine='rbind') %:%
foreach(ii=1:nc,.combine='rbind') %dopar%{
			simx(rn,tass,tasn,ii,kk)}

sx0=foreach(kk=1:nsim,.combine='rbind') %:%
foreach(ii=1:nc,.combine='rbind') %dopar%{
            simx(rn,tass,tasn,ii,kk)}

stopCluster(cl)
			
## simulate outcome variable

tbsn=c(-0.3,0.8,1.3)
tbss = c(-0.5,1, 1.5)
p=length(tbss)
gc=0.25
dsace=0.1#c(0,0.1)
tbss0=tbss+dsace
tsat=tt*(1-trho)

psmodel="logit"

cc=detectCores()
cl=makeCluster(cc)
registerDoParallel(cl)

if(psmodel=="logit"){
sda1=foreach(kk=1:nsim,.combine='rbind') %:%
foreach(ii=1:nc,.combine='rbind') %dopar%{

ceta1=eta1[ii,kk]
xk=subset(sx1,simid==kk)
xki=xk[which(xk$cid==ii),]
sxp=xlp(xki,gc*ceta1)
simy1(sxp,ceta1,tbss,tbsn,tsat,ii,kk)
}

sda0=foreach(kk=1:nsim,.combine='rbind') %:%
foreach(ii=1:nc,.combine='rbind') %dopar%{

ceta0=eta0[ii,kk]
xk0=subset(sx0,simid==kk)
xki0=xk0[which(xk0$cid==ii),]
sxp0=xlp(xki0,gc*ceta0)
simy0(sxp0,ceta0,tbss0,tsat,ii,kk)
}}

if(psmodel=="probit"){

sda1=foreach(kk=1:nsim,.combine='rbind') %:%
foreach(ii=1:nc,.combine='rbind') %dopar%{

ceta1=eta1[ii,kk]
xk=subset(sx1,simid==kk)
xki=xk[which(xk$cid==ii),]
sxp=xpp(xki,gc*ceta1)
simy1(sxp,ceta1,tbss,tbsn,tsat,ii,kk)
}

sda0=foreach(kk=1:nsim,.combine='rbind') %:%
foreach(ii=1:nc,.combine='rbind') %dopar%{

ceta0=eta0[ii,kk]
xk0=subset(sx0,simid==kk)
xki0=xk0[which(xk0$cid==ii),]
sxp0=xpp(xki0,gc*ceta0)
simy0(sxp0,ceta0,tbss0,tsat,ii,kk)
}}
	
stopCluster(cl)

simo = list(sda1=sda1,sda0=sda0)

save(simo,file="simdata.RData")

sm0=colMeans(sda0[,c("yss","ss","sn","nn")])
sm1=colMeans(sda1[,c("yss","ss","sn","nn")])

tsace=-dsace
ssim=list(sm0,sm1,tsace)
save(ssim,file="ssim.RData")
