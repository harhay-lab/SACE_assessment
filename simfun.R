library(MASS)
library(BuyseTest)

# simulate cluster size from a Gamma distribution
# shp: shape parameter
# scp: scale paramter
# nc: number of clusters
# nsim: number of simulations
simn=function(shp,scp,nc,nsim){
	rn=matrix(round(rgamma(nc*nsim,shape=shp,scale=scp)),nc,nsim)
cond=(rn<5)
    if(any(cond)) rn[cond]=5
	return(rn)}

# simulate predictors
# rn: matrix of cluster size
# tass: alpha_ss
# tasn: alpha_sn
# ii: cluster id
# kk: simulation index

simx=function(rn,tass,tasn,ii,kk){
		
px=1/2
sdx=1

ni=rn[ii,kk]
	
x1=rbinom(ni,1,px)
x2=rnorm(ni,0,sdx)

xm=cbind(1,x1,x2)
	lss=xm%*%tass
        lsn=xm%*%tasn
			
xobj=data.frame(x0=1,x1=x1,x2=x2,id=1:ni,cid=ii,simid=kk,lss=lss,lsn=lsn)	
return(xobj)}

# simulate principal stratum membership probability from a multinomial 
# logit model
# gi: the strength of unobserved cluster-level principal membership confounding 

xlp=function(xobj,gi){

        ess=exp(xobj$lss+gi)
        esn=exp(xobj$lsn+gi)
        xobj$pnn=pnn=1/(1+ess+esn)
        xobj$pss=ess*pnn
        xobj$psn=esn*pnn
        
        return(xobj)
        }

# simulate principal stratum membership probability from a Probit model
# gi: the strength of unobserved cluster-level principal membership confounding 

xpp=function(xobj,gi){

        xobj$pss=pnorm(xobj$lss+gi)
        xobj$psn=(1-xobj$pss)*pnorm(xobj$lsn+gi)
        xobj$pnn=1-xobj$pss-xobj$psn
        
        return(xobj)
        }

# simulate the outcome variable in the treatment group
# eta1: random intercept
# tbss: beta_ss
# tbsn: beta_sn
# tsat: variance parameter
# ii: cluster id
# kk: simulation index
simy1=function(sx1,eta1,tbss,tbsn,tsat,ii,kk){

probi=as.matrix(sx1[,c("pss","psn","pnn")])
	xi=as.matrix(sx1[,c("x0","x1","x2")])
	ni=nrow(xi)
				
    sta=matrix(0,ni,3)
    for(j in 1:ni){
    sta[j,]=rmultinom(n=1,size=1,prob=probi[j,])
	}
			
    stsat=sqrt(tsat)	
	myiss=xi%*%tbss
	yiss=rnorm(ni,myiss,sd=stsat)+eta1		
	myisn=xi%*%tbsn
	yisn=rnorm(ni,myisn,sd=stsat)+eta1

	yi=rowSums(cbind(yiss,yisn)*sta[,1:2])
	yind=rowSums(sta[,1:2])
	dai=data.frame(y=yi,yss=yiss,ysn=yisn,pss=probi[,1],psn=probi[,2],pnn=probi[,3],
	ss=sta[,1],sn=sta[,2],nn=sta[,3],id=1:ni,cid=ii,simid=kk,x0=1,
	x1=xi[,"x1"],x2=xi[,"x2"],yind=yind)		

return(dai)}

# simulate the outcome variable in the control group
# eta0: random intercept
# tbss0: beta_ss0
# tsat: variance parameter
# ii: cluster id
# kk: simulation index
simy0=function(sx0,eta0,tbss0,tsat,ii,kk){

probi=as.matrix(sx0[,c("pss","psn","pnn")])
	xi=as.matrix(sx0[,c("x0","x1","x2")])
	ni=nrow(xi)
        
        sta=matrix(0,ni,3)
        for(j in 1:ni){
        sta[j,]=rmultinom(n=1,size=1,prob=probi[j,])
		}
	
        stsat=sqrt(tsat)
						
	myi=xi%*%tbss0
	yi=rnorm(ni,myi,sd=stsat)+eta0			

	dai=data.frame(y=yi*sta[,1],yss=yi,pss=probi[,1],psn=probi[,2],pnn=probi[,3],
	ss=sta[,1],sn=sta[,2],nn=sta[,3],id=1:ni,cid=ii,simid=kk,x0=1,
	x1=xi[,"x1"],x2=xi[,"x2"],yind=sta[,1])		

return(dai)}

# SACE estimation by the worst ranked value approach 

qreg=function(da1,da0){

	da0$treat=0
	da1$treat=1
	
	da0$ny=da0$y
	ind0=which(da0$yind==0)
	da0$ny[ind0]=quantile(da0$y[ind0],prob=0.95)

	da1$ny=da1$y
	ind1=which(da1$yind==0)
	da1$ny[ind1]=quantile(da1$y[ind1],prob=0.95)
		
	cn=c("treat","ny","y","yind","x1","x2")	
	da=rbind(da0[,cn],da1[,cn])
	
	reg=lm(ny~x1+x2+treat,data=da)
	nda1=nda0=da
	nda1$treat=1
	nda0$treat=0
	
	reg1=predict(reg,nda1)
	reg0=predict(reg,nda0)
	
return(mean(reg1)-mean(reg0))
}

# SACE estimation by the inverse probability weighting approach 
ipe=function(da1,da0){

	da0$treat=0
	da1$treat=1
	cn=c("treat","y","yind","x1","x2")	
	da=rbind(da0[,cn],da1[,cn])

	sda=subset(da,yind==1)
	go=glm(treat~x1+x2,data=sda,family="binomial")
	ps=go$fit
	
	ipwe=mean(sda$y*sda$treat/ps)-mean(sda$y*(1-sda$treat)/(1-ps))
	return(ipwe)
}

# SACE estimation by the win/ratio approach 
wrf=function(da1,da0){
	da0$treat=0
	da1$treat=1
	cn=c("treat","y","yind")	
	da=rbind(da0[,cn],da1[,cn])
        bobj=BuyseTest(endpoint=c("yind","y"),type=c("b","c"),
             treatment="treat",operator=c(">0","<0"),data=da,trace=0)
        ss=confint(bobj,statistic="WinRatio")$p.value[2]
        return(ss)
}

#### calculate the frequency of rejecting the null hypothesis SACE=0 
#### and the frequency of covering the true sace value
### tsace: true SACE value
bsummary=function(bobj,tsace){

pb=c(0.025,0.975)
if(is.list(bobj)) { bobj=bobj[which(names(bobj)=="")]
                    bobj=c(unlist(bobj))
                       }
ind=which(is.na(bobj))
if(length(ind)>0) bobj=bobj[-ind]

bq=quantile(bobj,probs=pb)

c0=((bq[1]<=0)&(bq[2]>=0)) # rejection
cds=((bq[1]<=tsace)&(bq[2]>=tsace)) # coverage
return(c(c0,cds))}	