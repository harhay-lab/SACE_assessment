library(doParallel)
source("fecode.R")
source("mecode.R")
source("simfun.R")

load("simdata.RData")

sda0=simo$sda0
sda1=simo$sda1

nsim=max(sda0$simid)
nc=max(sda0$cid)

cc=detectCores()
cl=makeCluster(cc)
registerDoParallel(cl)

ptime=proc.time()
# mixed effects approach
mo=foreach(m=1:nsim,.combine='c',.multicombine=T,.maxcombine=nsim,
.errorhandling ="remove") %dopar%{
source("mecode.R")
da0=subset(sda0,simid==m)
da1=subset(sda1,simid==m)
fit(da1,da0)}

mt=proc.time()-ptime

ptime=proc.time()
# fixed effects approach
fo=foreach(m=1:nsim,.combine='c',.multicombine=T,.maxcombine=nsim,
.errorhandling ="remove") %dopar%{
da0=subset(sda0,simid==m)
da1=subset(sda1,simid==m)
ffit(da1,da0)}

ft=proc.time()-ptime

ptime=proc.time()
# inverse probability weighting
io=foreach(m=1:nsim,.combine='c',.multicombine=T,.maxcombine=nsim,
.errorhandling ="remove") %dopar%{
da0=subset(sda0,simid==m)
da1=subset(sda1,simid==m)
ipe(da1,da0)}

it=proc.time()-ptime

ptime=proc.time()
# worst ranked value
qo=foreach(m=1:nsim,.combine='c',.multicombine=T,.maxcombine=nsim,
.errorhandling ="remove") %dopar%{
da0=subset(sda0,simid==m)
da1=subset(sda1,simid==m)
qreg(da1,da0)}

qrt=proc.time()-ptime

ptime=proc.time()
# win/ratio 
wo=foreach(m=1:nsim,.combine='c',.multicombine=T,.maxcombine=nsim,
.errorhandling ="remove") %dopar%{
da0=subset(sda0,simid==m)
da1=subset(sda1,simid==m)
wrf(da1,da0)}

wt=proc.time()-ptime

stopCluster(cl) 

timeo=list(mt,ft,it,qrt,wt)
save(timeo,file="time.RData")

fobj=list(mo=mo,fo=fo,io=io,qo=qo,wo=wo)

save(fobj,file="fobj.RData")

#####

nbp=5 # number of bootstrap samples

cl=makeCluster(cc)
registerDoParallel(cl)

### mixed effects approach

bmo=foreach(m=1:nsim) %:%
foreach(icount(nbp),.combine='c',.multicombine=T,.maxcombine=nbp,.errorhandling ="pass") %dopar% {
source("mecode.R")
da0=subset(sda0,simid==m)
da1=subset(sda1,simid==m)

s1=sample(nc,nc,replace=T)
s0=sample(nc,nc,replace=T)

bda1=bda0=NULL

for(jj in 1:nc){
da1j=subset(da1,cid==s1[jj])
da1j$cid=jj
bda1=rbind(bda1,da1j)}

for(kk in 1:nc){
da0k=subset(da0,cid==s0[kk])
da0k$cid=kk
bda0=rbind(bda0,da0k)}

t1=table(bda1$cid,bda1$yind)[,2]
if(any(t1==0)) stop
t0=table(bda0$cid,bda0$yind)[,2]
if(any(t0==0)) stop

fit(bda1,bda0)
}
stopCluster(cl) 

######## fixed effects approach

cl=makeCluster(cc)
registerDoParallel(cl)

bfo=foreach(m=1:nsim) %:% 
foreach(icount(nbp),.combine='c',.multicombine=T,.maxcombine=nbp,
.errorhandling ="pass") %dopar% {

da0=subset(sda0,simid==m)
da1=subset(sda1,simid==m)

m1=nrow(da1)
m0=nrow(da0)

s1=sample(m1,m1,replace=T)
s0=sample(m0,m0,replace=T)

bda1=da1[s1,]
bda0=da0[s0,]

t1=table(bda1$yind)[2]
if(any(t1==0)) stop
t0=table(bda0$yind)[2]
if(any(t0==0)) stop

ffit(bda1,bda0) 
}
stopCluster(cl)

######## inverse probability weighting

cl=makeCluster(cc)
registerDoParallel(cl)

bio=foreach(m=1:nsim)%:%
foreach(icount(nbp),.combine='c',.multicombine=T,.maxcombine=nbp,
.errorhandling ="pass") %dopar% {

da0=subset(sda0,simid==m)
da1=subset(sda1,simid==m)

m1=nrow(da1)
m0=nrow(da0)

s1=sample(m1,m1,replace=T)
s0=sample(m0,m0,replace=T)

bda1=da1[s1,]
bda0=da0[s0,]

t1=table(bda1$yind)[2]
if(t1==0) stop
t0=table(bda0$yind)[2]
if(t0==0) stop

ipe(bda1,bda0) 
}
stopCluster(cl)

######## worst ranked value

cl=makeCluster(cc)
registerDoParallel(cl)

bqo=foreach(m=1:nsim)%:%
foreach(icount(nbp),.combine='c',.multicombine=T,.maxcombine=nbp,
.errorhandling ="pass") %dopar% {

da0=subset(sda0,simid==m)
da1=subset(sda1,simid==m)

m1=nrow(da1)
m0=nrow(da0)

s1=sample(m1,m1,replace=T)
s0=sample(m0,m0,replace=T)

bda1=da1[s1,]
bda0=da0[s0,]

t1=table(bda1$yind)[2]
if(t1==0) stop
t0=table(bda0$yind)[2]
if(t0==0) stop

qreg(bda1,bda0) 
}
stopCluster(cl)

bobj=list(bmo=bmo,bfo=bfo,bio=bio,bqo=bqo)

save(bobj,file="bobj.RData")

