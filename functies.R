OC<-function(N=NULL,n,a,p.max=1,type="both",step=0.01,rect=F,add=F,lty=2) {
	res<-T
	p.s<-seq(0,p.max,step)
	i<-1
	beta.lot<-p.s
	beta.pop<-p.s
	if((type=="both")||(type=="lot")) {
		n.samples.tot<-choose(N,n)
		for(p in p.s) {
			n.defect.lot<-round(p*N)
			n.samples.except<-0
			for(j in 0:a) {
				n.samples.except<-n.samples.except+choose(N-n.defect.lot,n-j)*choose(n.defect.lot,j)
			}
			beta.lot[i]<-n.samples.except/n.samples.tot
			i<-i+1
		}
		if(!add) plot(p.s,beta.lot,ylim=c(0,1),xlab="p",ylab="beta",type="l")
		if(add) lines(p.s,beta.lot,lty=lty)
	}
	beta.pop<-pbinom(a,size=n,prob=p.s)
	if(type=="pop") {
		if(!add) plot(p.s,beta.pop,ylim=c(0,1),xlab="p",ylab="beta",type="l")
            if(add) lines(p.s,beta.pop,lty=lty)
	}
	if(type=="both") {
		if(!add) lines(p.s,beta.pop,lty=2)
		if(add) lines(p.s,beta.pop,lty=(lty+1))
	}
	if(rect) {
		tmp<-100*beta.pop*p.s
		tmp.s<-seq(0,1,0.2)
		tmp1.s<-round(tmp.s*max(tmp),2)
		axis(4,at=tmp.s,labels=paste(as.character(tmp1.s),"%"),cex=0.5)
		if(!add) lines(p.s,(tmp/max(tmp)),lty=3)
		if(add) lines(p.s,(tmp/max(tmp)),lty=lty+2)
		mtext(side=4,"AOQ")
	}
	if(res) return(data.frame(p.s,beta.pop,beta.lot))
}

OC.plot<-function(p,beta,add=F,lty=1) {
	if(!add) plot(p,beta,ylim=c(0,1),xlab="p",ylab="beta",type="l",lty=lty)
	else lines(p,beta,lty=lty)
}

plot.OC<-function(OCobj,add=F,lty=1,type="pop") {
	if(is.null(OCobj$beta.pop)) {
		OCobj$beta.pop<-OCobj$beta # voor het variabele plan
	}
	if(type=="lot") OCobj$beta.pop<-OCobj$beta.lot
	p<-OCobj$p.s
	beta<-OCobj$beta.pop
	if(!add) plot(p,beta,ylim=c(0,1),xlab="p",ylab="beta",type="l",lty=lty)
	else lines(p,beta,lty=lty)
}



OC2<-function(n1,a1,r1,n2,a2,p.max=1,step=0.01,add=F,lty=2) {
	res<-T
	p.s<-seq(0,p.max,step)
	beta.1<-pbinom(a1,size=n1,prob=p.s)
	beta.2<-rep(0,length(p.s))
	for(m in (a1+1):(r1-1)) {
		beta.2<-beta.2+(pbinom((a2-m),n2,p.s)*dbinom(m,n1,p.s))
	}
	beta<-beta.1+beta.2
	if(!add) plot(p.s,beta,ylim=c(0,1),xlab="p",ylab="beta",type="l")
	if(add) lines(p.s,beta,lty=lty)
	if(res) return(data.frame(p.s,beta))
}

OCv<-function(n,Qc,p.max=1,step=0.01,add=F,lty=2) {
	res<-T
	p.s<-seq(0,p.max,step)
	q.s<-qnorm(1-p.s)
	beta<-pnorm(sqrt(n)*(q.s-Qc))
	if(!add) plot(p.s,beta,ylim=c(0,1),xlab="p",ylab="beta",type="l")
	if(add) lines(p.s,beta,lty=lty)
	if(res) return(data.frame(p.s,beta))
}


beta.pop.1<-function(n,a,p) {
	beta<-pbinom(a,size=n,prob=p)
	return(beta)
}

beta.lot.1<-function(n,a,p,N) {
	n.samples.tot<-choose(N,n)
	n.defect.lot<-round(p*N)
	n.samples.except<-0
	for(j in 0:a) {
		n.samples.except<-n.samples.except+choose(N-n.defect.lot,n-j)*choose(n.defect.lot,j)
	}
	beta<-n.samples.except/n.samples.tot
	return(beta)
}

beta.pop.2<-function(n1,a1,r1,n2,a2,p) {
	beta.1<-pbinom(a1,size=n1,prob=p)
	beta.2<-0
	for(m in (a1+1):(r1-1)) {
		beta.2<-beta.2+(pbinom((a2-m),n2,p)*dbinom(m,n1,p))
	}
	beta<-beta.1+beta.2
	return(beta)
}

beta.lot.2<-function(n1,a1,r1,n2,a2,p,N) {
	n.samples.tot.1<-choose(N,n1)
	n.defect.lot<-round(p*N)
	beta.1<-beta.lot.1(n1,a1,p=p,N=N)
	beta.2<-0
	for(m in (a1+1):(r1-1)) {
		beta.2<-beta.2+beta.lot.1(n2,(a2-m),p=p,N=N)*choose(n.defect.lot,m)*choose((N-n.defect.lot),(n1-m))/n.samples.tot.1
	}
	beta<-beta.1+beta.2
	return(beta)
}

OC2.lot<-function(N,n1,a1,r1,n2,a2,p.max=1,step=0.01,res=F,add=F,lty=2) {
	p.s<-seq(0,p.max,step)
	beta<-beta.lot.2(N=N,n1=n1,a1=a1,r1=r1,n2=n2,a2=a2,p=p.s)
	if(!add) plot(p.s,beta,ylim=c(0,1),xlab="p",ylab="beta",type="l")
	if(add) lines(p.s,beta,lty=lty)
	if(res) return(data.frame(p.s,beta))
}



AQL<-function(OCobj,beta=0.95,type="pop") {
	if(is.null(OCobj$beta.pop)) OCobj$beta.pop<-OCobj$beta # voor het variabele plan
	if(type=="lot") OCobj$beta.pop<-OCobj$beta.lot
	p<-OCobj$p.s[abs(OCobj$beta.pop-beta)<0.03]
	diff<-abs(OCobj$beta.pop[abs(OCobj$beta.pop-beta)<0.03]-beta)
	p<-p[diff==min(diff)]
	return(p)
}

LTPD<-function(OCobj,beta=0.10) {
	if(is.null(OCobj$beta.pop)) OCobj$beta.pop<-OCobj$beta # voor het variabele plan
	p<-OCobj$p.s[abs(OCobj$beta.pop-beta)<0.03]
	diff<-abs(OCobj$beta.pop[abs(OCobj$beta.pop-beta)<0.03]-beta)
	p<-p[diff==min(diff)]
	return(p)
}

AOQL<-function(OCobj) {
	if(is.null(OCobj$beta.pop)) OCobj$beta.pop<-OCobj$beta # voor het variabele plan
	pb<-OCobj$p.s*OCobj$beta.pop
	aoql<-max(pb)
	return(aoql)
}

read.plot<-function(OCobj,p=NULL,beta=NULL,rect=F) {
	if(is.null(OCobj$beta.pop)) OCobj$beta.pop<-OCobj$beta # voor het variabele plan
	if(rect) OCobj$beta.pop<-OCobj$beta.pop*OCobj$p.s
	if(is.null(p)) {
		p2<-OCobj$p.s[abs(OCobj$beta.pop-beta)<0.03]
	   diff<-abs(OCobj$beta.pop[abs(OCobj$beta.pop-beta)<0.03]-beta)
	   p2<-p2[diff==min(diff)]
	   print(c(p2,beta))
	}
	if(is.null(beta)) {
		beta2<-OCobj$beta.pop[abs(OCobj$p.s-p)<0.01]
	   diff<-abs(OCobj$p.s[abs(OCobj$p.s-p)<0.01]-p)
	   beta2<-beta2[diff==min(diff)]
		print(c(p,beta2))
	}
	return()
}

############################################################################################
####################    STATISTICAL PROCESS CONTROL   ######################################
############################################################################################

MakeControlLimit <- function(statistics.vector, lowerBound = 0.005, upperBound=0.995){
  #returns the LCL qnd UCL for a vector of statistics of simulated samples and the lower and upper bound
  statistics.vector <- sort(as.vector(statistics.vector))
  number.simulations <- length(statistics.vector)
  LCL.estimate <- statistics.vector[round(lowerBound*number.simulations)]
  UCL.estimate <- statistics.vector[round(upperBound*number.simulations)]
  return(c(LCL.estimate,UCL.estimate))
}