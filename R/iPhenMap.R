## iPhenMap.R
# Wrapper function to load data, preprocess, and run iPhenMAP analysis

#' Run iPhenMAP Analysis
#'
#' This function reads clinical, CNA, proteomics, metabolomics, and RNA-seq data,
#' preprocesses and matches samples, and then runs the iPhenMAP model.
#'
#' @param data_dir Path to the parent data directory containing subfolders: Clinical, CNA, Proteomics, Metabolomics, RNAseq.
#' @param q Number of latent variables (default 10)
#' @param mcmc Number of MCMC iterations (default 100000)
#' @param ... Additional parameters passed to \code{iphen} (g, c, d, nu0, Sigo, s2o, nu0k, noiseMod, scaleType, iniType, time, status)
#' @return Saves results in "iPhenMAPresults" and returns the final MCMC output list.
#' @export
#’ @name iPhenMAP
#’ @title Fit the iPhenMAP Model

iphen <- function(D, X,
                  q = 10,
                  mcmc = 10000,
                  g = nrow(X),
                  c = 1,
                  d = 1,
                  nu0 = .01,
                  Sigo = .01,
                  s2o = .01,
                  nu0k = .01,
                  noiseMod = "unique",
                  scaleType = "unit",
                  iniType = "EVD",
                  time = NULL,
                  status = NULL) {
  # … your existing implementation …
}









iphen<-function(D,X,q=10,mcmc=10000,g=1000,c=1,d=1,nu0=.01,Sigo=.01,s2o=.01,nu0k=.01,noiseMod='unique',scaleType='unit',iniType='EVD',time=NULL,status=NULL){

  ####################################
  ## installing missing R packages
  installed<-installed.packages()[,1] 
  required<-c("MCMCglmm","cluster","fpc","dplyr","MetabolAnalyze","ggplot2","RColorBrewer","ggrepel","circlize","MASS","ConsensusClusterPlus",
              "Rcpp","scales","digest","gridExtra","reshape2","plyr","lattice","cowplot","ComplexHeatmap","MCMCpack","siggenes",
              "pamr","mixtools","colorspace","GetoptLong","Hmisc","plotrix","survminer","AnnotationDbi","org.Hs.eg.db","ReactomePA",
              "psych","mclust","survival","devtools","githubinstall","NMF")
  toinstall<-required[!(required %in% installed)]
  if(length(toinstall) != 0){
    if(as.numeric(paste0((unlist(getRversion()))[1:2],collapse = '.'))>3.4){
      BiocManager::install(toinstall)}else{
    source("https://bioconductor.org/biocLite.R")
    biocLite(toinstall)
  }
  } 
  lapply(required, require, character.only = TRUE)
  n.cores=detectCores()
  ## theme for ggplot
  theme=theme(strip.text.y = element_text(),#rotate strip text to horizontal
              axis.text = element_text(colour = "black", family="Helvetica", size=18),
              axis.title.x = element_text(colour = "black", family="Helvetica",face="bold", size=20),
              axis.title.y = element_text(colour = "black", family="Helvetica",face="bold", size=20),
              axis.text.x = element_text(colour = "black", family="Helvetica", size=18),#remove y axis labels
              legend.background = element_rect(fill = "white"),
              panel.grid.major = element_line(colour = "white"),
              legend.title = element_text(colour = "black", family="Helvetica",face="bold", size=20),
              legend.text = element_text(colour="black", family="Helvetica", size = 18, face = "bold"),
              strip.text.x = element_text(colour = "black", family="Helvetica", size = 18),
              plot.title = element_text(size=20, family="Helvetica", face="bold"),
              panel.border=element_rect(colour="black",size=2)) #format title
  
  ########################
  ## required functions 
  ########################
  plotSc<-function(sc, l1, l2, tau, col, grps=grps, cex=1.2, pchs=16){
    xtemp = range(sc[,c(l1)])
    xlims = c(xtemp[1]-.01, xtemp[2]+.01)
    ytemp = range(sc[,c(l2)])
    ylims = c(ytemp[1]-.01,ytemp[2]+.01)
    plot(1, type="n", xlab=paste("Dimension ", colnames(sc)[l1], sep=""), 
         ylab=paste("Dimension ", colnames(sc)[l2], sep=""), 
         main='MV-clustering',
         xlim=xlims, 
         ylim=ylims)
    points(sc[,l1],sc[,l2], cex=cex, col=scales::alpha(col[tau],.8), pch=pchs)
    co=c(min(c(l1,l2)),max(c(l1,l2)))
    ntau=table(tau)
    for(i in 1:length(grps)){
      if(ntau[i]>5) mixtools:::ellipse(mu=colMeans(sc[tau==i,co]),sigma=cov(sc[tau==i,co]),alpha=.25,npoints=nrow(sc[tau==i,]),col=col[i],lwd=2,lty=2)
    }
    abline(h=0, col=adjustcolor("grey50",.2), lty=2,lwd=3)
    abline(v=0, col=adjustcolor("grey50",.2), lty=2,lwd=3)
    legend('bottomleft',grps,col=col,pch=16,bty='n',pt.cex=2,ncol=1)
    box(lwd=4)
  }
  
  phenplot <- function(cest, l1, l2, dela=0.01, nlv,labs){
    x <- cest[,l1]
    y <- cest[,l2]
    ind1<-rep(.9,nrow(cest)) 
    ind1[abs(x)>=1.96]<-1.15
    ind1[abs(y)>=1.96]<-1.15
    ind2<-rep("grey50",nrow(cest)) 
    ind2[abs(x)>=1.96]<-"blue4"
    ind2[abs(y)>=1.96]<-"blue4"
    rg1=c(min(x)-dela,max(x)+dela); rg2=c(min(y)-dela,max(y)+dela)
    par(mgp=c(1.5,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
    plot(x, y, xlab=paste0("scaled coeff MV-",l1), ylab=paste0("scaled coeff MV-",l2),xlim=rg1,ylim=rg2,pch=16,col=scales::alpha(ind2,.2),cex=4,main="Phenotypes",cex.lab=1.5)
    text(x, y, labs, col=scales::alpha(ind2,.75),cex=ind1) 
    abline(h=c(-1.96,1.96), col="red", lty=2,lwd=1.5)
    abline(v=c(-1.96,1.96), col="red", lty=2,lwd=1.5)
    box(lwd=4)
  }

  plotNice <- function(oload, l1, l2, dela=0.001, nlv,labs){
    x <- oload[,l1]
    y <- oload[,l2]
    ind1<-rep(.7,nrow(oload)) 
    ind1[abs(x)>=3]<-.9
    ind1[abs(y)>=3]<-.9
    ind2<-rep("grey50",nrow(oload)) 
    ind2[abs(x)>=3]<-"blue4"
    ind2[abs(y)>=3]<-"blue4"
    rg1=c(min(x)-dela,max(x)+dela); rg2=c(min(y)-dela,max(y)+dela)
    par(mgp=c(1.5,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
    plot(x, y, xlab=paste0("scaled coeff MV-",l1), ylab=paste0("scaled coeff MV-",l2),xlim=rg1,ylim=rg2,pch=16,col=scales::alpha(ind2,.2),main="",cex.lab=1.5)
    text(x, y, labs, col=adjustcolor(ind2,.6),cex=ind1) 
    abline(h=c(-3,3), col="red", lty=2,lwd=1.5)
    abline(v=c(-3,3), col="red", lty=2,lwd=1.5)
    box(lwd=4)
  }
    
  loadingplot=function(dat,onames=NULL,ndim=1,threshold=.2,fac=NULL,lab=NULL,cex.lab=0.5,lab.jitter=0,srt=0,adj=NULL){
    myLoadings=(dat[,ndim])^2
    sgn=sign(dat[,ndim])
    lsign=as.numeric(sgn<0)+1
    myLoad=myLoadings*sgn
    names(myLoad)=rownames(dat)
    x=myLoad
    hthresh=quantile(x,1-threshold/2)
    lthresh=quantile(x,threshold/2)
    
    ## some checks
    if(is.data.frame(x)|is.matrix(x)){
      if(is.null(lab)) {lab<-rownames(x)}
      x=x[,ndim]
    }else{ if(is.null(lab)){lab=names(x)} }
    
    lab <- names(x)  # x is a mtrix 
    if(!is.null(onames)){ lab=onames}
    
    if(!is.numeric(x)) stop("x is not numeric")
    if(any(is.na(x))) stop("NA entries in x")
    
    ## handle lab
    if(is.null(lab)) {lab=1:length(x)}
    
    ## handle fac
    if(!is.null(fac)){
      fac <- factor(fac, levels=unique(fac))
      grp.idx <- cumsum(table(fac)) + 0.5
      grp.lab.idx <- tapply(1:length(x), fac, mean) # average of each group
      grp.lab <- names(grp.idx)
      grp.idx <- grp.idx[-length(grp.idx)]
    }
    
    ## preliminary computations
    y.min <- min(min(x),0)
    y.max <- max(max(x),0)
    y.offset <- (y.max)*0.04
    x.offset <- abs(y.min)*0.04
    
    ## start the plot
    at=1:length(x)
    dat=cbind(at, x)  
    plot(dat, type='h', xlab=expression(bold("features")), ylab=paste0('MV',ndim,'-relevence'),cex=1,lwd=2,
         main='', xaxt="n",xlim=c(-15,max(at)), ylim=c(y.min*1.1,y.max*1.1),
         col=c("orange","darkgreen")[lsign])
    
    ## add groups of variables (optional)
    if(!is.null(fac)){
      abline(v=grp.idx,col='black',lty=2,lwd=2) # split groups of variables
      text(x=grp.lab.idx,y=y.max*1.12, labels=grp.lab, cex=1,font=4) # annotate groups
    }
    
    ## annotate variables that are above the threshold
    if(sum((x<lthresh)|(x>hthresh))>0){
      x.ann= at[(x<lthresh)|(x>hthresh)]
      x.ann= jitter(x.ann,factor=lab.jitter)
      yout=x[(x<lthresh)|(x>hthresh)]
      y.ann=yout
      y.ann[yout>hthresh]<-(y.ann[yout>hthresh]+ y.offset)
      y.ann[yout<lthresh]<-(y.ann[yout<lthresh]- x.offset) 
      y.ann= jitter(y.ann,factor=lab.jitter)
      txt.ann= lab[(x<lthresh)|(x>hthresh)]
      text(x=x.ann, y=y.ann, label=txt.ann, col=(c("orange","darkgreen")[lsign])[(x<lthresh)|(x>hthresh)],cex=cex.lab, srt=srt, adj=adj)
      #abline(h=c(lthresh,hthresh), col="grey50",lty=2,lwd=1.5)
      box(lwd=4)
      
      ## build the result
      res<-list(threshold=threshold,
                var.names=txt.ann,
                var.idx=which((x<lthresh)|(x>hthresh)),
                var.values=x[(x<lthresh)|(x>hthresh)])
      return(invisible(res))
    }
    return(NULL) # if no point above threshold
  } # end loadingplot

  unlistfn<-function(x,g,byrow) matrix(unlist(x),ncol=g,byrow=byrow)
  unlistarray<-function(x,p,q,g) array(unlist(x), c(p,q,g))
  rowVar<-function(x){
    rowSums((x-rowMeans(x))^2)/(dim(x)[2]-1)
  }
  
  ## fast extraction of diagonal entries
  dag<-function(x,p){x[1+0:(p-1)*(p+1)]}
  
  ## fast extraction of diagonal entries and summing
  mdiag<-function(x,p){sum(x[1+0:(p-1)*(p+1)])}
  
  fnx=function(p,q){ 
    nq = c(0,rep(p,q))
    lapply(1:q,function(i){ (sum(nq[1:i])+1):(sum(nq[1:i])+nq[i+1]) })
  }
  chol.inv<-function(x){chol2inv(chol(x))}
  
  scaling=function(tD, scaleType,ovSD){
    if(scaleType=="pareto"){ dt=t(tD/(rowVar(tD)^(.25))) }
    if(scaleType=="unit"){ dt=t(tD/(rowVar(tD)^(.5))) }
    if(scaleType=='view'){ dt=t(tD/sd(as.vector(tD))) } # scaling each data matrix Xi by some norm ||Xi||. Then each data type will contribute equally to the integrative solution. 
    if(scaleType=='overall'){ dt=t(tD/ovSD) }
    if(scaleType=="none"){ dt=t(tD) }
    return(dt)}
  
  iniMod=function(dat,N,S,p,q,Iq,iniType){
    if(q>N) stop("Shut down ... q cannot be greater than N.")
    if((iniType=="EVD")|(iniType=="SVD")){
      if(iniType=="EVD"){
        tryCatch(if(p>10000){
          tp=10000;tmp=eigs_sym(S,k=tp,which="LM");tp=length(!is.na(tmp$val))
        }else{ tp=p;tmp=eigen(S) },error=function(e) NULL)
        sig=abs(sum(tmp$val[(q+1):length(tmp$val)])/(length(tmp$val)-q))
        w=tmp$vec[,1:q,drop=F] }
      if(iniType=="SVD"){
        svdY=svd(dat)
        v=svdY$v[,1:q,drop=F]
        if(p>N){ sig=abs(sum(svdY$d[(q+1):N])/(N-q))
        }else{ sig=abs(sum(svdY$d[(q+1):p])/(p-q)) }
        w=v%*%(svdY$d[1:q]*Iq)  }
      u=(dat%*%(w%*%chol.inv(crossprod(w)+sig*Iq)))
    }else{
      u=scale(matrix(rnorm(N*q),N,q)) # random init LVs
      sig=1
      decomM=eigen( crossprod(u)/sig+Iq )
      tmp=decomM$vectors*outer(rep(1,q),1/sqrt(decomM$values))
      w=crossprod(dat,u)%*%(tcrossprod(tmp))/sig + matrix(rnorm(p*q),p,q)%*%tmp
    }
    list(u=u,sig=sig,w=w)}
  ###############################################################################
  
  ########################
  ## fit the model
  itime=proc.time()     
  ## constants 
  n=nrow(D[[1]]); nT=length(D); Q = q;  repn=rep(1,n)
  desT=1:nT; L=ncol(X)+1; rL=rep(0,L); Il=diag(L) 
  p_t = unlist(lapply(D,ncol)); p = sum(p_t)     
  nq = c(0,rep(p,q)); np = c(0,p_t) 
  qndx=lapply(1:q,function(i){ (sum(nq[1:i])+1):(sum(nq[1:i])+nq[i+1]) })
  fin=rep(1:nT,p_t) 
  ########
  
  ## preprocessing the data
  Y=NULL; D_t<-list()
  indx<-id<-vector("list",length=nT) 
  ovSD=sd(as.vector(unlist(D)))
  for(t in 1:nT){
    astart=(sum(np[1:t])+1)
    aend=(sum(np[1:t])+np[t+1])
    indx[[t]]=astart:aend
    id[[t]]=rep(1,p_t[t])
    D_t[[t]]=t(t(D[[t]])-colMeans(D[[t]]))   # mean center
    D_t[[t]]=scaling(t(D_t[[t]]),scaleType,ovSD)  # scale data
    Y=cbind(Y,D_t[[t]])} 
  tY=t(Y)     
  ## normalize the covariates
  tcov=cbind(repn,standardize(as.matrix(X))) 
  covars=t(tcov)
  
  ########
  set.seed(123) 
  ## initialize parameters 
  St=(tY%*%Y)/n      # data covariance matrix
  iniUWS=iniMod(Y,n,St,p,q,diag(q),iniType)
  isig=rep(1/iniUWS$sig,p)
  W=iniUWS$w; U=iniUWS$u
  theta=matrix(1,p,q)
  z=matrix(1,nT,q)
  ppi=rep(0.5,q)
  
  ## some proirs and constants
  #nu0=.5; sig0=.001; Sigo=nu0*sig0
  XtX <- crossprod(tcov)
  cholCv=chol( XtX + (1/g)*Il )
  IR=backsolve(cholCv,Il)    
  v.b=tcrossprod(IR)
  Beta = v.b%*%crossprod(tcov,U) 
  tbeta = t(Beta)
  v.B=tcrossprod(v.b,tcov)
  residue <- U - tcov%*%Beta
  s2 <- dag(crossprod(residue),q)/n
  is2 = 1/s2  
  Hg = diag(n) - tcrossprod(tcov%*%v.b,tcov)
  
  ###########
  ## running the mcmc sampler
  S=mcmc; nthin=1000; burn=250; si=0
  thinstep=mcmc/nthin; iniS=1
  thin=seq(1,S,by=thinstep)
  burnPeriod=thinstep*burn
  tsamp=(S-burnPeriod)/thinstep
  sig.store=matrix(NA,p,tsamp)
  
  nthin=1250; burn=500; si=0
  thinstep=mcmc/nthin; iniS=1
  thin=seq(1,mcmc,by=thinstep)
  t0=round((thinstep*burn)/2.667,0)
  t1=round((thinstep*burn)/16,0)
  altburn=rep(c(FALSE,TRUE),11)
  bPeriods=rep(altburn,c(t0-1,1,rep(c(t1-1,1),10)))
  burnPeriod=length(bPeriods)
  sel=rep(TRUE,Q)
  fsel=rep(TRUE,p)
  tsamp=(mcmc-burnPeriod)/thinstep
  for(s in 1:mcmc){  
    ## update loadings
    UY=crossprod(U,Y); iUU=crossprod(U)
    diagU=diag(iUU); diag(iUU)=0
    wnorm=rnorm(p*q)
    for(k in 1:q){
      v.l=diagU[k]*isig+theta[,k]
      m.l=isig*c(UY[k,]-W%*%iUU[k,])*(v.l^(-1))      # Note here if lapply columns of W are not updated
      W[,k]=(v.l^(-0.5))*wnorm[qndx[[k]]]+m.l
      lpy=0.5*( log(theta[,k])-log(v.l)+v.l*m.l^2 )  # NOTE:here lpy is view spcific
      if(s>iniS){ 
        for(t in 1:nT){  # wiki logit
          Rt=sum(lpy[indx[[t]]])+log(ppi[k])-log(1-ppi[k]) 
          z[t,k]=rbinom(1,1,(1+exp(-Rt))^(-1))
        }
      } # NOTE: runif(1) generate fractions btwn 0 & 1 equally likelly  
      W[,k]=W[,k]*rep(z[,k],p_t) }#k
    
    ## remove redundant LVs
    if(s <= burnPeriod){ 
      if(bPeriods[s]){
        sel=(colSums(z)!=0)
        if(sum(sel)==0){ stop("Shut down ... no LVs found.") } 
        q=sum(sel)
        z=z[,sel,drop=FALSE]
        W=W[,sel,drop=FALSE]
        theta=theta[,sel,drop=FALSE]
        U=U[,sel,drop=FALSE]
        Beta=Beta[,sel,drop=FALSE]
        is2=is2[sel]
        ppi=ppi[sel,drop=FALSE]
      } 
    }
    if(s %% 1000 == 0){ cat("Iteration ",s,"; ",p," features; ",q," components. \r\r\r") }
    
    ## update LVs
    Iq=diag(q)
    iS2=is2*Iq
    Wisig=W*isig
    uR=backsolve(chol(iS2+crossprod(Wisig,W)),Iq)
    U=((Y%*%Wisig)+tcov%*%(Beta%*%iS2))%*%tcrossprod(uR)+tcrossprod(matrix(rnorm(n*q),n,q),uR)
    
    # draw sig
    resid=Y-tcrossprod(U,W) 
    psqr=colSums(resid^2)  
    if(noiseMod=='unique'){ isig=rgamma(p,shape=.5*n+nu0,rate=Sigo+.5*psqr) }else{  
      isig=do.call('c',lapply(desT,function(t){ rgamma(1,shape=.5*n*p_t[[t]]+nu0[t],rate=Sigo[t]+.5*sum(psqr[indx[[t]]]))*id[[t]] }))
    }
    
    # draw s2
    is2=rgamma(q,.5*n+nu0k,.5*dag(crossprod(U,Hg)%*%U,q)+s2o)
    
    # draw Beta
    Beta = v.B%*%U + tcrossprod(IR,matrix(rnorm(q*L),q,L)*(is2^(-.5)))
    
    ## update ppi 
    if(s>iniS){ ppi=(colSums(z)+1)/(nT+2) } 
    
    ## update theta
    W2=W^2
    theta=unlistfn(lapply(1:p,function(j){ rgamma(q,shape=.5*z[fin[j],]+c,rate=.5*W2[j,]+d) }),q,byrow=TRUE)
    
    ## Storing parameter estimates and latent variables
    if(s > burnPeriod){
      if(any(s==thin)){
        si=si+1
        if(si==1){
          target=W
          loadings.store <- oload.store <- array(NA, c(p, q, tsamp))
          scores.store <- array(NA, c(n, q, tsamp))
          b.store <- ob.store <- array(NA, c(L, q, tsamp))
          scores.store[,,si]<-U
          b.store[,,si] <- ob.store[,,si] <- Beta
          loadings.store[,,si] <- oload.store[,,si] <- target 
          is2.store=matrix(NA,q,tsamp)
          sig.store=matrix(NA,p,tsamp)
        }else{
          ob.store[,,si] <- Beta
          oload.store[,,si] <- W
          # correct rotational ambiquity
          rs=MCMCpack:::procrustes(W, target, translation=FALSE, dilation=FALSE)
          ld=rs$X.new
          for(j in 1:nT){ ld[indx[[j]],z[j,]==0]=0 }
          scores.store[,,si]<-U%*%rs$R
          b.store[,,si]<-Beta%*%rs$R
          loadings.store[,,si]<-ld 
          is2.store[,si]<-is2
          sig.store[,si]<-isig
        } 
      }
    }
  }#s  
  
  ## Storing parameter estimates and latent variables
  d1 <- unlist(lapply(indx,function(x){x[1]})) # Find the first index of each view.
  if(TRUE) theta <- theta[d1,]  # grouped
  
  ## creating folder for storing results and setting the working directory
  dir.create("iPhenMAPresults", showWarnings = FALSE)
  dir.create("iPhenMAPresults/TracePlots", showWarnings = FALSE)
  dir.create("iPhenMAPresults/PoV", showWarnings = FALSE)
  dir.create("iPhenMAPresults/Phenotypes", showWarnings = FALSE)
  dir.create("iPhenMAPresults/Phenotypes/Pairwise", showWarnings = FALSE)
  dir.create("iPhenMAPresults/Prognosis", showWarnings = FALSE)
  dir.create("iPhenMAPresults/MetaGenes", showWarnings = FALSE)
  dir.create("iPhenMAPresults/MetaGenes/Pairwise", showWarnings = FALSE)
  dir.create("iPhenMAPresults/Clustering", showWarnings = FALSE)
  dir.create("iPhenMAPresults/Clustering/ModelBased", showWarnings = FALSE)
  dir.create("iPhenMAPresults/Clustering/ModelBased/DFA", showWarnings = FALSE)
  dir.create("iPhenMAPresults/Clustering/ModelBased/Pairwise", showWarnings = FALSE)
  dir.create("iPhenMAPresults/Clustering/Consensus", showWarnings = FALSE)
  dir.create("iPhenMAPresults/Clustering/Consensus/DFA", showWarnings = FALSE)
  dir.create("iPhenMAPresults/Clustering/Consensus/DFA/MVs", showWarnings = FALSE)
  dir.create("iPhenMAPresults/Clustering/Consensus/Pairwise", showWarnings = FALSE)
  dir.create("iPhenMAPresults/MCMCsamples", showWarnings = FALSE)
  dirp<-getwd()
  setwd(paste0(dirp,"/iPhenMAPresults"))
  
  ## Diagnostics plots: trace plots
  pdf("TracePlots/tracePlots.pdf")
  par(mfrow=c(2,2))
  sq=sample(1:q,1)
  sp=sample(1:p,ceiling(.2*p))
  sp= if(p<10) 1:p else sample(1:p,ceiling(.2*p))
  matplot(t(scores.store[sample(1:n,ceiling(.25*n)),sq,]), type="l",main="Factors trace plot",
          xlab="thinned iteration number",ylab=expression(bold(U)))
  matplot(t(loadings.store[sp,sq,]), type="l",main="Loadings trace plot",
          xlab="thinned iteration number",ylab=expression(bold(W)))
  matplot(t(b.store[,sq,]), type="l",main="Regression coefficients trace plot",
          xlab="thinned iteration number",ylab=expression(bold(beta)))
  matplot(t(sig.store[sp,]), type="l",main="Covariance trace plot",
          xlab="thinned iteration number",ylab=expression(bold(sigma^2)))
  dev.off()
  
  ## summarize 
  keep=colSums(z)>0
  scores=data.frame(apply(scores.store,c(1,2), median)[,keep]); dim(scores)
  wstore=apply(loadings.store,c(1,2), median)[,keep]; dim(wstore)
  bstore=apply(b.store,c(1,2), median)[,keep]; dim(bstore)
  z=z[,keep,drop=FALSE]
  rownames(z) = if(is.null(names(D))) paste0('plat:',1:nrow(z)) else names(D)
  colnames(z)<-colnames(scores)<-paste0('MV',1:ncol(scores))
  fnames=unlist(lapply(D,colnames))
    
  ## PoV
  varfil=matrix(NA,nrow(z),ncol(scores))
  for(s in 1:nrow(z)){
    Vtot<-sum(diag(crossprod(D_t[[s]]))) 
    for(k in 1:ncol(scores)){
      if(z[s,k]==1){   
        mfit=wstore[indx[[s]],k,drop=FALSE]%*%t(scores[,k,drop=FALSE]) 
      }else{ mfit=0 }
      resid = t(D_t[[s]])-mfit
      rss=sum(diag(tcrossprod(resid)))
      vprop = 1-rss/Vtot
      varfil[s,k]=if(vprop<0) 0 else vprop
    }
  }
  colnames(varfil)<-colnames(scores)
  vfil=NULL
  for(j in 1:nT){ vfil=rbind(vfil,cumsum(varfil[j,])) }
  rownames(vfil)=rownames(z)
  df2=reshape2::melt(t(vfil))
  colnames(df2)=c('MV','Type','PoV')
  write.table(df2,"PoV/PoV.txt",quote=TRUE,sep="\t",row.names=TRUE)
  p1<-ggplot(df2, aes(x=MV, y=PoV, group=factor(Type))) +ylim(c(0,1))+
    geom_line(aes(color=factor(Type)), size=1.5, alpha = 0.7)+
    geom_point(aes(color=factor(Type)), size=3)+
    labs(title="",x=" ", y="PoV")+theme_bw()+theme+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(filename="PoV/PoV.pdf", device=cairo_pdf,plot = p1, width = 20, height = 20, units = "cm")
  
  pdf(file='PoV/ActivatedMVs.pdf')
  print(Heatmap(
    z,name="Active",
    cluster_columns = FALSE,
    cluster_rows = FALSE,rect_gp = gpar(col= "black",lwd=3.5),
    row_names_side = "left", row_names_gp = gpar(fontsize = 17),
    column_names_gp = gpar(fontsize = 17),
    show_row_dend = FALSE,  show_column_dend = FALSE,show_heatmap_legend=FALSE,
    col = colorRamp2(c(0, 1), c("white", "grey50"))
  ))
  dev.off()
  
  ## Regression coefficient plots
  sbet=b.store[-1,keep,]; dim(sbet)
  if(L<3){
    beta<-t(as.matrix(apply(sbet,1, median)))
    bload<-t(as.matrix(apply(sbet,1,mean)))
    bsd<-t(as.matrix(apply(sbet,1,sd)))
  }else{
    beta<-apply(sbet,c(1,2), median)
    bload<-apply(sbet,c(1,2),mean)
    bsd<-apply(sbet,c(1,2),sd)
  }
  cest<-bload/bsd
  rownames(cest)<-rownames(beta)<-colnames(X)
  colnames(cest)<-colnames(beta)<-colnames(scores)
  eload<-data.frame(wstore)
  SD=apply(loadings.store[,keep,],c(1,2),sd)
  sdload<-data.frame(apply(loadings.store[,keep,],c(1,2),mean)/SD)
  colnames(eload)<-colnames(sdload)<-colnames(scores)
  
  # pheno heatmap
  hDat=cest
  hDat[abs(hDat)<1.96]<-0
  pdf('Phenotypes/phenotypesHeatmap.pdf')
  print(Heatmap(
    hDat,name="Effect",
    cluster_columns = FALSE,
    cluster_rows = FALSE,rect_gp = gpar(col= "black",lwd=.5),
    row_names_side = "left", row_names_gp = gpar(fontsize = 8.5),
    column_names_gp = gpar(fontsize = 11),
    show_row_dend = TRUE,  show_column_dend = FALSE,
    col = colorRamp2(c(-(round(max(hDat),0)-1), 0, (round(max(hDat),0)-1)), c("darkgreen", "white", "orange"))
  ))
  dev.off()
  pdf(paste0("Phenotypes/scaled-Regression-CoeffPlot.pdf"))
  par(mgp=c(1.5,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
  x <- cest[,1]
  y <- cest[,2]
  ind1<-rep(.9,nrow(cest)) 
  ind1[abs(x)>=1.96]<-1.15
  ind1[abs(y)>=1.96]<-1.15
  ind2<-rep("grey50",nrow(cest)) 
  ind2[abs(x)>=1.96]<-"blue4"
  ind2[abs(y)>=1.96]<-"blue4"
  rg1=range(cest[,1]); rg2=range(cest[,2])
  plot(x,y,xlim=c(min(x)-.5,max(x)+.5),ylim=c(min(y)-.5,max(y)+.5), 
       xlab=expression(bold("scaled coeff MV:1")), ylab=expression(bold("scaled coeff MV:2")),pch=16,col=scales::alpha(ind2,.2),main="Phenotypes")
  text(x,y, colnames(X),col=scales::alpha(ind2,.8), cex=ind1)
  abline(h=c(-1.96,1.96), col="red", lty=2,lwd=1.5)
  abline(v=c(-1.96,1.96), col="red", lty=2,lwd=1.5)
  box(lwd=4)
  dev.off()
  write.table(beta,"Phenotypes/regressionCoefficients.txt",quote=TRUE,sep="\t",row.names=TRUE)
  write.table(cest,"Phenotypes/scaledRegressionCoefficients.txt",quote=TRUE,sep="\t",row.names=TRUE)
  
  ## MV prognosis
  if(!is.null(time)){
    sfit=Surv(time,status)
    td=data.frame(sfit,scores)
    formula=as.formula(paste0('sfit~',paste(paste0(colnames(scores), collapse = '+')))) 
    cox=coxph(formula,td)
    rcox=summary(cox)
    Cancer=colnames(scores)
    df_data=data.frame(Cancer,rcox$coefficients[,1],confint(cox)); dim(df_data)
    rownames(df_data)=NULL; colnames(df_data)=c('Cancer','HR','HR_lower','HR_upper')
    df_data$Cancer<-factor(df_data$Cancer, levels=df_data$Cancer[q:1])
    p <- ggplot(df_data, aes(x=Cancer, y=HR, ymin=HR_lower, ymax=HR_upper)) + 
      geom_linerange(size=3,col='darkgreen') +
      geom_hline(aes(yintercept=0), lty=4,col=2,size=.8) +
      geom_point(size=1, shape=21, fill=4,colour = "white", stroke = 1) +
      coord_flip() +   labs(title="",x=" ", y="Log10 (HR)")+
      theme(axis.text.y = element_text( hjust = 1,family="Helvetica",face="bold", size=15))
    ggtitle(" ") + theme#
    ggsave(filename="Prognosis/MetaVariableProgronosis.pdf", device=cairo_pdf,plot = p, width = 20, height = 20, units = "cm")
  }
  
  ## Loadings plots
  snams = if(sum(table(colnames(Y))>1)>0) paste0('p',1:ncol(Y)) else colnames(Y)
  rownames(eload)<-rownames(sdload)<-snams
  par(mgp=c(1.5,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
  pdf(paste0('MetaGenes/LoadingsMatplot.pdf'),width=8)
  matplot((eload),type="l",main='',xlab='features',ylab='loadings'); box(lwd=3)
  abline(h=c(-2,2), col="red", lty=2)
  dev.off()
  nact=colSums(z)
  tload=sdload
  tload[is.na(tload)]<-0
  for(k in 1:ncol(scores)){
    ifeat<-colnames(Y)[(abs(tload[,k])>3)]
    tld=data.frame(eload)
    lnames=colnames(Y); lnames[is.na(match(lnames,ifeat))]<-' '
    dType=if(is.null(names(D))){
      factor(do.call('c',lapply(1:nT,function(x) rep(paste0('Platform-',1:nT)[x],p_t[x]))))
    }else{ factor(do.call('c',lapply(1:nT,function(x) rep(names(D)[x],p_t[x])))) }
    pdf(paste0('MetaGenes/MV-',k,'-scaledLoadings.pdf'),width=8)
    par(mgp=c(1.5,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
    loadingplot(tload,onames=lnames,ndim=k,threshold=1,fac=dType) 
    abline(h=c(-3^2,3^2),lty=2,col=2,lwd=2)
    box(lwd=5)
    dev.off()
    # pairwise loadings
    for(j in 1:ncol(scores)){
      if((k!=j)&(k<j)){
        pdf(paste0('MetaGenes/Pairwise/',colnames(scores)[k],'-vs-',colnames(scores)[j],'-loadings.pdf'),width=8)
        plotNice(oload=tload,l1=k,l2=j,dela=0.5,nlv=colnames(scores),labs=colnames(Y))
        dev.off()
      }
    }
    # pairwise phenotypes
    for(j in 1:ncol(scores)){
      if((k!=j)&(k<j)){
        pdf(paste0('Phenotypes/Pairwise/',colnames(scores)[k],'-vs-',colnames(scores)[j],'-phenotypes.pdf'),width=8)
        phenplot(cest=cest,l1=k,l2=j,dela=0.5,nlv=colnames(scores),labs=colnames(X))
        dev.off()
      }
    }
    
    # Divergent plot
    ord=order(tld[,k],decreasing = TRUE)
    otl=tld[ord,]
    tla=data.frame(otl[c(1:5,(nrow(tld)-4):nrow(tld)),k]); dim(tla) 
    colnames(tla)='ld'
    mtc=tla
    mtc$`gename`=rownames(otl)[c(1:5,(nrow(tld)-4):nrow(tld))]  
    mtc$ld=round(tla$ld, 2)  
    mtc$mpg_type <- ifelse(mtc$ld < 0, "below", "above")  
    mtc <- mtc[order(mtc$ld), ]  # sort
    mtc$`gename` <- factor(mtc$`gename`, levels = mtc$`gename`)  
    # Diverging Barcharts
    gn1=ggplot(mtc, aes(x=`gename`, y=ld, label=ld)) + 
      geom_bar(stat='identity', aes(fill=mpg_type), width=.5)  +
      scale_fill_manual(name="Mileage", 
                        labels = c("Above Average", "Below Average"), 
                        values = c("above"="indian Red", "below"="#C7E9C0")) + 
      labs(subtitle="", title= "", x='',y=paste0('MV-',k)) +
      coord_flip()  + theme_bw()+ theme(legend.position = 'none')
    ggsave(filename =paste0('MetaGenes/Top10-features-MV-',k,'.pdf'), device=cairo_pdf,plot=gn1, width = 20, height = 20, units = "cm")
  }
  fil1=paste0("MetaGenes/loadingsCoefficients.txt")
  fil2=paste0("MetaGenes/scaledLoadingsCoefficients.txt")
  write.table(eload,fil1,quote=TRUE,sep="\t",row.names=TRUE)
  write.table(sdload,fil2,quote=TRUE,sep="\t",row.names=TRUE) 
  
  ## Clustering MVs 
  if(n<60) gmax=6 else gmax=8 
  # consensus clustering  
  kmccp=ConsensusClusterPlus(t(scores), maxK=gmax, reps=1000, clusterAlg='km', title='km', seed=123) 
  consensusMatrices=lapply(kmccp[2:gmax], function(x){matrix(unlist(x['consensusMatrix']), n, n)}) 
  consensuscls=lapply(kmccp[2:gmax], function(x){unlist(x['consensusClass'])}) 
  sil=function(clss, cns){ meanWidth=mean(silhouette(clss, dmatrix=1-cns)[,3]) }
  silWidth=mapply(sil, consensuscls, consensusMatrices)
  coph=sapply(consensusMatrices, cophcor)
  subr=2:gmax
  tdat3=as.data.frame(cbind(coph,silWidth,subr)) 
  colnames(tdat3)=c("coph","silh","subtypes")
  cop=tdat3$coph[-1]
  gopt=(3:gmax)[which(cop==max(cop)[length(max(cop))])]
  mm3=reshape2::melt(subset(tdat3, select=c(subtypes,coph, silh)), id.var="subtypes")
  fg3=ggplot(mm3, aes(x = subtypes, y = value)) + geom_line(aes(color = variable),size=2,linetype = 4) +theme_bw()+theme+
    scale_colour_manual("variable",labels=c("coph","silh"),values=brewer.pal(n=5,name="Set1")[1:2])+
    geom_vline(aes(xintercept=gopt), lty=4,col=2,size=1) +
    facet_grid(variable ~ ., scales = "free_y") + theme(legend.position = "none");fg3
  ggsave(filename="Clustering/Consensus/Cophenetic_coeff.pdf", device=cairo_pdf,plot=fg3)
  dtau=kmccp[[gopt]]$consensusClass
  nfns=length(names(table(dtau)))
  grps<-paste("iSub",1:nfns,sep="")
  COLS=brewer.pal(8, 'Set1')
  cols=COLS[1:gopt]
  rg1=range(scores[,1]); rg2=range(scores[,2])
  
  ## 3D plot
  gdata=data.frame(scores,Cluster=cols[dtau])
  colnames(gdata)<-c(paste0('MV',1:ncol(scores)),'Cluster')
  if(ncol(scores)>2){
    pdf(file='Clustering/Consensus/3DLVplot.pdf',pointsize=15)
    pp=(cloud(MV1 ~ MV2*MV3 ,data=gdata,
          groups=Cluster,
          col=adjustcolor(cols,.8),
          screen = list(x = -90, y = 70), 
          zoom = .9,
          pch=16,cex=2,
          scales=list(lwd=2, col = "blue"), ## formating the arrows
          par.box=list(lwd=4, col="black"), ## formating the axis lines
          par.settings=list(axis.line = list(lwd = 0)) ## removing the outer border
    ))
    print(pp)
    dev.off()
  }

  ddir0="Clustering/Consensus/DFA/MVs/DensityPlots"
  dir.create(ddir0, showWarnings = FALSE)
  ddir001="Clustering/Consensus/DFA/MVs/BoxPlots"
  dir.create(ddir001, showWarnings = FALSE)
  scpval=rep(NA,ncol(scores))
  for(i in 1:ncol(scores)){
    d12<-data.frame(x=scores[,i],y=grps[dtau])
    colnames(d12)<-c("lv", "iSubs")
    scpval[i]<-kruskal.test(d12$lv~factor(d12$iSubs))$p.value
    (p12<-ggplot(d12, aes(factor(iSubs), lv,fill = factor(iSubs))) + geom_boxplot()+
        scale_fill_manual("iSubs",labels=grps,values=cols)+theme_bw()+theme+theme(legend.position="none")+
        labs(title=paste0('p = ',round(scpval[i],4)),x=" ",y=paste0("MV",i) ))
    ggsave(filename =paste0(ddir001,'/MV',i,'_by_iSubs_boxplot.pdf'), device=cairo_pdf,plot = p12)
    wd = as_data_frame(data.frame(biofac = factor(grps[dtau]), weight = scores[,i]))
    l <- density(wd$weight)
    den2<-ggplot(wd, aes(x = weight)) +
      geom_density(aes(fill = biofac), alpha = 0.8) + 
      xlim(range(l$x))+
      labs(title="",x=paste0("MV",i)) +theme_bw()+theme+
      scale_fill_manual("subtypes",labels=grps,values=cols)+
      guides(colour = guide_legend(override.aes = list(linetype=0)))
    pdf(paste0(ddir0,'/MV',i,'_by_iSubs_density_plot.pdf'))
    print(den2)
    dev.off()
  }
   
  # sam-pam
  dd=NULL
  for(i in 1:nT){ dd=cbind(dd,D[[i]])}
  for(i in 1:nT){
    ddir=paste0("Clustering/Consensus/DFA/",names(D)[i])
    dir.create(ddir, showWarnings = FALSE)
    set.seed(123)
    d1=t(Y[,indx[[i]]]); dim(d1) 
    d1u=t(dd[,indx[[i]]]); dim(d1u) 
    tryCatch({
      sam.out=NULL
      sam.out<-sam(d1, dtau, rand=123)},error=function(ex) {print(ex$message)})
    if(!is.null(sam.out)){
      samTab=data.frame(attributes(summary(sam.out))$mat.fdr)
      write.table(samTab,paste0(ddir,"/SAM_table.txt"),quote=TRUE,sep="\t",row.names=FALSE)
      if(any(samTab$FDR<.05)){
        if(any(samTab$FDR>0.04999)){
          rdel=findDelta(sam.out, fdr=0.05)
          del = if(is.matrix(rdel)) rdel[1,1] else rdel[1]
          pdf(paste0(ddir,"/SAM_plot_delta_",del,".pdf"))
          plot(sam.out,del)
          box(lwd=3)
          dev.off()
          featsel=names(summary(sam.out,del,entrez=FALSE)@row.sig.genes)
          write.table(matrix(featsel),paste0(ddir,"/SAM_genes_delta_",del,".txt"),quote=TRUE,sep="\t",row.names=FALSE)
        }else{ 
          del=min(samTab[samTab$False<.15,1])
          pdf(paste0(ddir,"/SAM_plot_delta_",del,".pdf"))
          plot(sam.out,del)
          box(lwd=3)
          dev.off()
          featsel=names(summary(sam.out,del,entrez=FALSE)@row.sig.genes)
          write.table(matrix(featsel),paste0(ddir,"/SAM_genes_delta_",del,".txt"),quote=TRUE,sep="\t",row.names=FALSE)
        }
        if(length(featsel)>2){
          m0=match(featsel,rownames(d1))
          w0=which(!is.na(m0)); tD=d1[m0[w0],]; dim(tD)
          Ypcen<-tD; dim(Ypcen)
          tDu=d1u[m0[w0],]; Ypcenu<-tDu
          mydata=list(x=apply(Ypcen,2,as.numeric),y=dtau, geneid=rownames(Ypcen),genenames=rownames(Ypcen))
          mydatau=list(x=apply(Ypcenu,2,as.numeric),y=dtau, geneid=rownames(Ypcenu),genenames=rownames(Ypcenu))
          mytrain=pamr.train(mydata)
          mytrainu=pamr.train(mydatau)
          mycv=pamr.cv(mytrain,mydata, nfold=5)
          mycvu=pamr.cv(mytrainu,mydatau, nfold=5)
          thresh=mycv$threshold[which(mycv$error==min(mycv$error))][1]
          threshu=mycvu$threshold[which(mycvu$error==min(mycvu$error))][1]
          pamr.confusion(mycv, threshold=thresh)
          pamr.confusion(mycvu, threshold=threshu)
          pamcen<-pamr.listgenes(mytrain, mydata, threshold=thresh)
          pamcenu<-pamr.listgenes(mytrainu, mydatau, threshold=threshu)
          colnames(pamcen)<-c("genes",colnames(pamcen)[-1]); dim(pamcen) 
          colnames(pamcenu)<-c("genes",colnames(pamcenu)[-1]); dim(pamcenu) 
          pdf(paste0(ddir,"/PAM_plot_threshold_",thresh,".pdf"))
          pamr.plotcv(mycv)
          dev.off()
          write.table(pamcen,paste0(ddir,"/PAM_centroid_threshold_",thresh,".txt"),quote=TRUE,sep="\t",row.names=FALSE)
          write.table(pamcenu,paste0(ddir,"/PAM_centroid_unscaled_threshold_",threshu,".txt"),quote=TRUE,sep="\t",row.names=FALSE)
          
          # heatmap
          m4=match(pamcen[,1],rownames(Ypcen)); 
          w4=which(!is.na(m4)); 
          G0=Ypcen[m4[w4],]; dim(G0)
          ord<-order(dtau)
          df<-data.frame(dtau[ord]) 
          colnames(df)<-c('Subtypes')
          dlist=list(Subtypes=c('1'=COLS[1],'2'=COLS[2],'3'=COLS[3],'4'=COLS[4],
                                '5'=COLS[5],'6'=COLS[6],'7'=COLS[7],'8'=COLS[8])[1:gopt])      
          ha = HeatmapAnnotation(df = df,col=dlist, which = "col")
          odat=data.frame(G0[,ord]); dim(odat)
          #fscaled = sweep(odat,1,apply(odat,1,median),'-'); dim(fscaled)
          fscaled = odat
          pdf(paste0(ddir,"/Heatmap.pdf"))
          print(Heatmap( 
            fscaled,name="X",
            cluster_columns = FALSE,
            cluster_rows = TRUE,
            show_heatmap_legend=TRUE,
            #top_annotation=ha,top_annotation_height = unit(8, "mm"),
            row_names_gp = gpar(fontsize = 7),
            show_column_names = FALSE,
            show_row_dend = TRUE,  show_column_dend = FALSE,
            col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
            heatmap_legend_param=list(color_bar = "continuous",at=seq(-3,3, by=3))
          ))
          dev.off() 
        }
      }
    } 
  }
 
  # model based 
  set.seed(123)
  yBIC2<-mclustBIC(scores,G=2:gmax,modelNames=c('EEI','EII','EEE')) 
  rs2<-summary(yBIC2,scores)
  tau<-rs2$classification  
  ngrps<-length(names(table(tau)))
  mgrps<-paste("iSub",1:ngrps,sep="")
  mcols=COLS[1:ngrps]
  d4 <- data.frame(x=2:gmax, y=yBIC2[,1])
  colnames(d4)<-c("models","BIC")
  rownames(d4)<-paste("g",2:gmax,sep="_")
  p4 <-ggplot(data=d4, aes(models, BIC)) + geom_point(size=3)+geom_line()+labs(title="", x="Number of subtypes", 
         y="BIC")+theme_bw() 
  ggsave(filename ='Clustering/ModelBased/clustering_BICplot.pdf', device=cairo_pdf,plot = p4, width = 20, height = 20, units = "cm")
  
  # subtypes prognosis
  if(!is.null(time)){
    sfit= Surv(time,status)
    ccd=cbind(time,status,dtau);
	ccd=as.data.frame(ccd)
    fitccc<- survminer::surv_fit(sfit ~ factor(dtau),data=ccd)
    resccc<-ggsurvplot(fitccc, main = "Survival curve",
                       pval = TRUE, # Add p-value
                       palette =cols,
                       conf.int = FALSE,
                       risk.table = TRUE,
                       risk.table.col = "strata",
                       risk.table.y.text.col = TRUE
    ); resccc
    ggsave(filename="Clustering/Consensus/prognosis-consensus-clustering.pdf", device=cairo_pdf,plot=print(resccc), width = 20, height = 20, units = "cm")
    mbc=data.frame(time,status,tau)
    fitmbc<- survminer::surv_fit(sfit ~ factor(tau),data=mbc)
    resmbc<-ggsurvplot(fitmbc, main = "Survival curve",
                       pval = TRUE, # Add p-value
                       palette =mcols,
                       conf.int = FALSE,
                       risk.table = TRUE,
                       risk.table.col = "strata",
                       risk.table.y.text.col = TRUE
    ); resmbc
    ggsave(filename="Clustering/ModelBased/prognosis-model-based-clustering.pdf", device=cairo_pdf,plot=print(resmbc), width = 20, height = 20, units = "cm")
  }
  
  # pairwise scores plot, pathway analysis, and heatmaps
  ddir0001="Clustering/MV_heatmaps"
  dir.create(ddir0001, showWarnings = FALSE)
  ddir2="PathwayAnalysis"
  dir.create(ddir2, showWarnings = FALSE)
  for(i in 1:ncol(scores)){
    for(j in 1:ncol(scores)){
      if((i!=j)&(i<j)){
        pdf(paste0('Clustering/Consensus/Pairwise/',colnames(scores)[i],'-vs-',colnames(scores)[j],'.pdf'),width=8)
        par(mgp=c(1.5,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
        plotSc(scores,l1=i,l2=j,tau=dtau, col=cols,grps=grps, cex=2, pchs=16)
        box(lwd=4)
        dev.off()
        pdf(paste0('Clustering/ModelBased/Pairwise/',colnames(scores)[i],'-vs-',colnames(scores)[j],'.pdf'),width=8)
        par(mgp=c(1.5,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
        plotSc(scores,l1=i,l2=j,tau=tau, col=mcols,grps=mgrps, cex=2, pchs=16)
        box(lwd=4)
        dev.off()
      }
    }
    # pathway analysis
    expgenesup <- expgenesdown <- progenesup <- progenesdown <- NULL
    expgenesup = if(z[1,i]==1) rownames(tload)[tload[indx[[1]],i] > 2]
    expgenesdown = if(z[1,i]==1) rownames(tload)[tload[indx[[1]],i] < (-2)]
    if(z[2,i]==1){
      tfeatup = rownames(tload)[ tload[indx[[2]],i] > 2 ]
      tfeatdown = rownames(tload)[ tload[indx[[2]],i] < (-2) ]
      if(length(tfeatup)!=0) {
        progenesup=unlist(lapply(1:length(tfeatup), function(x) strsplit(tfeatup[x], "_")[[1]][2]))
        progenesup=progenesup[!is.na(progenesup)]
      }
      if(length(tfeatdown)!=0) {
        progenesdown=unlist(lapply(1:length(tfeatdown), function(x) strsplit(tfeatdown[x], "_")[[1]][2]))
        progenesdown=progenesdown[!is.na(progenesdown)]
      }
    } 
    pfeatup=c(expgenesup,progenesup)
    pfeatdown=c(expgenesdown,progenesdown)
    if(!is.null(pfeatup)){
      if(length(pfeatup)!=0){
        deup=AnnotationDbi::select(org.Hs.eg.db, keys = pfeatup,columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")
        deup=deup[!is.na(deup[,2]),2]; length(deup)
        xup <- enrichPathway(gene=deup,pvalueCutoff=0.1, readable=T)
        head(as.data.frame(xup))
        pdf(paste0(ddir2,"/pathway_MV",i,"_up_regulated_genes.pdf"))
        print(dotplot(xup, showCategory=15))
        dev.off()
      }
    }
    if(!is.null(pfeatdown)){
      if(length(pfeatdown)!=0){
        dedown=AnnotationDbi::select(org.Hs.eg.db, keys = pfeatdown,columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")
        dedown=dedown[!is.na(dedown[,2]),2]; length(dedown)
        xdown <- enrichPathway(gene=dedown,pvalueCutoff=0.1, readable=T)
        head(as.data.frame(xdown))
        pdf(paste0(ddir2,"/pathway_MV",i,"_down_regulated_genes.pdf"))
        print(dotplot(xdown, showCategory=15))
        dev.off()
      }
    }
    afys<-datsplit<-NULL
    for(t in 1:nT){
      if(z[t,i]==1){
        ytj=Y[,indx[[t]]]
        xo=order(abs(tload[indx[[t]],i]) ,decreasing=TRUE)
        ys=ytj[,xo]
        optm=5
        fys=ys[,1:optm]; dim(fys)
        afys=cbind(afys,fys)
        datsplit=c(datsplit,rep(names(D)[t],optm))
      }
    }
    g0=t(afys); dim(g0)
    ord<-order(scores[,i])
    df<-data.frame(scores[ord,i],dtau[ord],tau[ord]) 
    colnames(df)<-c(colnames(scores)[i],'consensus','modelBased')
    csub=c('1'=COLS[1],'2'=COLS[2],'3'=COLS[3],'4'=COLS[4],
      '5'=COLS[5],'6'=COLS[6],'7'=COLS[7],'8'=COLS[8])[1:gopt]
    sci=colorRamp2(c(range(df[,1])[1],0,range(df[,1])[2]), c("blue",'white',"red"))
    dlist=list(mv=sci,consensus=csub) 
    names(dlist)[1]<-colnames(df)[1]
    ha = HeatmapAnnotation(df = df,col=dlist, which = "col")
    odat=data.frame(g0[,ord]); dim(odat)
    fscaled = sweep(odat,1,apply(odat,1,median),'-'); dim(fscaled)
    #fscaled = odat
    pdf(paste0(ddir0001,"/heatmap_MV",i,".pdf"))
      print(Heatmap( 
        fscaled,name="expression",
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_heatmap_legend=TRUE,
        #top_annotation=ha,top_annotation_height = unit(30, "mm"),
        row_names_gp = gpar(fontsize = 7),
        show_column_names = FALSE,
        split=  datsplit,
        show_row_dend = TRUE,  show_column_dend = FALSE,
        col = colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red")),
        heatmap_legend_param=list(color_bar = "continuous",at=seq(-2.5,2.5, by=2.5))
      ))
    dev.off()
  }

  if(is.null(rownames(Y))) rownames(Y)=1:nrow(Y)
  sc=data.frame(rownames(Y),scores,tau,dtau)
  colnames(sc)=c('samples',paste0('MV',1:q),'MBC-subtypes','ConsensusClusters')
  fil5=paste0("Clustering/MV_clustering.txt")
  write.table(sc,fil5,quote=TRUE,sep="\t")
  itime <- proc.time()[1] - itime
  
  # sam-pam
  for(i in 1:nT){
    ddirmb=paste0("Clustering/ModelBased/DFA/",names(D)[i])
    dir.create(ddirmb, showWarnings = FALSE)
    set.seed(123)
    d1=t(Y[,indx[[i]]]); dim(d1) 
    d1u=t(dd[,indx[[i]]]); dim(d1u) 
    tryCatch({
      sam.out=NULL
      sam.out<-sam(d1, tau, rand=123)},error=function(ex) {print(ex$message)})
    if(!is.null(sam.out)){
      samTab=data.frame(attributes(summary(sam.out))$mat.fdr)
      write.table(samTab,paste0(ddirmb,"/SAM_table.txt"),quote=TRUE,sep="\t",row.names=FALSE)
      if(any(samTab$FDR<.05)){
        if(any(samTab$FDR>0.04999)){
          rdel=findDelta(sam.out, fdr=0.05)
          del = if(is.matrix(rdel)) rdel[1,1] else rdel[1]
          pdf(paste0(ddirmb,"/SAM_plot_delta_",del,".pdf"))
          plot(sam.out,del)
          box(lwd=3)
          dev.off()
          featsel=names(summary(sam.out,del,entrez=FALSE)@row.sig.genes)
          write.table(matrix(featsel),paste0(ddirmb,"/SAM_genes_delta_",del,".txt"),quote=TRUE,sep="\t",row.names=FALSE)
        }else{ 
          del=min(samTab[samTab$False<.15,1])
          pdf(paste0(ddirmb,"/SAM_plot_delta_",del,".pdf"))
          plot(sam.out,del)
          box(lwd=3)
          dev.off()
          featsel=names(summary(sam.out,del,entrez=FALSE)@row.sig.genes)
          write.table(matrix(featsel),paste0(ddirmb,"/SAM_genes_delta_",del,".txt"),quote=TRUE,sep="\t",row.names=FALSE)
        }
        if(length(featsel)>2){
          m0=match(featsel,rownames(d1))
          w0=which(!is.na(m0)); tD=d1[m0[w0],]; dim(tD)
          Ypcen<-tD; dim(Ypcen)
          tDu=d1u[m0[w0],]; Ypcenu<-tDu
          mydata=list(x=apply(Ypcen,2,as.numeric),y=tau, geneid=rownames(Ypcen),genenames=rownames(Ypcen))
          mydatau=list(x=apply(Ypcenu,2,as.numeric),y=tau, geneid=rownames(Ypcenu),genenames=rownames(Ypcenu))
          mytrain=pamr.train(mydata)
          mytrainu=pamr.train(mydatau)
          mycv=pamr.cv(mytrain,mydata, nfold=5)
          mycvu=pamr.cv(mytrainu,mydatau, nfold=5)
          thresh=mycv$threshold[which(mycv$error==min(mycv$error))][1]
          threshu=mycvu$threshold[which(mycvu$error==min(mycvu$error))][1]
          pamr.confusion(mycv, threshold=thresh)
          pamr.confusion(mycvu, threshold=threshu)
          pamcen<-pamr.listgenes(mytrain, mydata, threshold=thresh)
          pamcenu<-pamr.listgenes(mytrainu, mydatau, threshold=threshu)
          colnames(pamcen)<-c("genes",colnames(pamcen)[-1]); dim(pamcen) 
          colnames(pamcenu)<-c("genes",colnames(pamcenu)[-1]); dim(pamcenu) 
          pdf(paste0(ddirmb,"/PAM_plot_threshold_",thresh,".pdf"))
          pamr.plotcv(mycv)
          dev.off()
          write.table(pamcen,paste0(ddirmb,"/PAM_centroid_threshold_",thresh,".txt"),quote=TRUE,sep="\t",row.names=FALSE)
          write.table(pamcenu,paste0(ddirmb,"/PAM_centroid_unscaled_threshold_",threshu,".txt"),quote=TRUE,sep="\t",row.names=FALSE)
          
          # heatmap
          m4=match(pamcen[,1],rownames(Ypcen)); 
          w4=which(!is.na(m4)); 
          G0=Ypcen[m4[w4],]; dim(G0)
          ord<-order(tau)
          df<-data.frame(tau[ord]) 
          colnames(df)<-c('Subtypes')
          dlist=list(Subtypes=c('1'=COLS[1],'2'=COLS[2],'3'=COLS[3],'4'=COLS[4],
                                '5'=COLS[5],'6'=COLS[6],'7'=COLS[7],'8'=COLS[8])[1:ngrps])      
          ha = HeatmapAnnotation(df = df,col=dlist, which = "col")
          odat=data.frame(G0[,ord]); dim(odat)
          #fscaled = sweep(odat,1,apply(odat,1,median),'-'); dim(fscaled)
          fscaled = odat
          pdf(paste0(ddirmb,"/Heatmap.pdf"))
          print(Heatmap( 
            fscaled,name="X",
            cluster_columns = FALSE,
            cluster_rows = TRUE,
            show_heatmap_legend=TRUE,
            #top_annotation=ha,top_annotation_height = unit(8, "mm"),
            row_names_gp = gpar(fontsize = 7),
            show_column_names = FALSE,
            show_row_dend = TRUE,  show_column_dend = FALSE,
            col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
            heatmap_legend_param=list(color_bar = "continuous",at=seq(-3,3, by=3))
          ))
          dev.off() 
        }
      }
    } 
  }
  
  
  ## Output
  out=list(u.store=scores.store,
           b.store=b.store, 
           w.store=loadings.store,
           s.store=sig.store,
           target=target,
           is2.store=is2.store,
           ob.store=ob.store,
           oload.store=oload.store,
           z=z, 
           itime=itime)
  save(out,file="MCMCsamples/iPhenMAP-mcmcSamples.RData")
}

