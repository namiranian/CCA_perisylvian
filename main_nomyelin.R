################sparsity selection################
if(TRUE){
  library(PMA)
  library(foreach)
  library(iterators)
  library(parallel)
  library(doParallel)
  numCores <- detectCores()
  registerDoParallel(numCores-1)  # use multicore, set to the number of our cores
  perisylvians=c(11,12,13,14,29,30,63,64,65,66,79,80,81,82,85,86)
  maindir="C:/SF_data/data_byregions/nomyelin/sparcityselection/"
  paths="C:/SF_data/morphologicalfeats.csv"
}

getfeatures <- function(path,S=FALSE){
  feat <- read.csv(path,header = FALSE)
  if (sum(feat[1,])==0){feat=feat[-1,]}
  shuffleind <- read.csv("C:/SF_data/suffleind.csv",header = TRUE)
  shuffleind=shuffleind[,2]
  if(S){
    initind=c(1,79,157,235,313,391)
    rmind=c(initind[1]:(initind[2]-1))
    feat <- feat[shuffleind,-rmind]
  }else{
    feat <- feat[shuffleind,]
  }
}

pmacca_jacknif=function(bestpenaltyx,bestpenaltyz,pathx,pathz,savepath, partname='',k=1){
  pacman::p_load(pacman,PMA)
  x <- getfeatures(path=pathx)
  z<- getfeatures(path=pathz)#,S=TRUE
  
  N=nrow(x)
  ccperms <- matrix(nrow = N, ncol= 1)
  us <- matrix(nrow = N, ncol= dim(x)[2])
  vs <- matrix(nrow = N, ncol= dim(z)[2])
  
  for (i in c(1:N)){
    X <- x[-i, ]
    Z <- z[-i, ]
    out <- PMA::CCA(X, Z, typex = "standard", typez = "standard",
               niter = 500, K = k, penaltyx = bestpenaltyx,
               penaltyz = bestpenaltyz,trace = FALSE)
    
    sdx <- apply(X, 2, sd)
    sdz <- apply(Z, 2, sd)
    X <- scale(X, TRUE, sdx)
    Z <- scale(Z, TRUE, sdz)
    
    if (mean(out$u == 0) != 1 && mean(out$v == 0) != 1) {
      ccperms[i,] <- cor(t(t(X)) %*% t(t(out$u)), t(t(Z)) %*% t(t(out$v)))
    }else {
      ccperms[i,] <- 0
    }
    us[i,]=out$u
    vs[i,]=out$v
    
  }
  
  savedir=paste(savepath,partname,sep="")
  write.csv(us, paste(savedir,"us.csv",sep=""))
  write.csv(vs, paste(savedir,"vs.csv",sep=""))
  write.csv(ccperms, paste(savedir,"ds.csv",sep=""))
}

pmacca_main=function(bestpenaltyx,bestpenaltyz,pathx,pathz,savepath,partname,permnums=500,maincal=FALSE,otherloadings=FALSE){
  pacman::p_load(pacman,PMA)
  x <- getfeatures(path=pathx)
  z<- getfeatures(path=pathz,S=TRUE)
  
  k= 1#min(dim(x)[2],dim(z)[2])
  if (maincal){
    mainout <- PMA::CCA(x, z, typex = "standard", typez = "standard",
                        niter = 500, K = k,
                        penaltyx = bestpenaltyx,penaltyz = bestpenaltyz,trace=FALSE)
    indmax=which.max(mainout$cors)
    mainout$cors=mainout$cors[indmax]
    mainout$u=mainout$u[,indmax]
    mainout$v=mainout$v[,indmax]
    
    savedir=paste(savepath,'/main',partname,sep="")
    write.csv(mainout$cors, paste(savedir,"cors.csv",sep=""))
    write.csv(mainout$u, paste(savedir,"u.csv",sep=""))
    write.csv(mainout$v, paste(savedir,"v.csv",sep=""))
    
    if(otherloadings){
      sdx <- apply(x, 2, sd)
      sdz <- apply(z, 2, sd)
      x <- scale(x, TRUE, sdx)
      z <- scale(z, TRUE, sdz)
      
      markerx = t(t(x)) %*% (t(mainout$u))
      markerz = t(t(z)) %*% (t(mainout$v))
      
      for (i in c(1:featx)){
        mainout$interloadingx <- cor(markerz, X[,i])
        mainout$intraloadingx <- cor(markerx, X[,i])
      }
      for (i in c(1:featz)){
        mainout$interloadingz <- cor(markerx, z[,i])
        mainout$intraloadingz <- cor(markerz, z[,i])
      }
      write.csv(mainout$interloadingx, paste(savedir,"interloadingsu.csv",sep=""))
      write.csv(mainout$intraloadingx, paste(savedir,"intraloadingsu.csv",sep=""))
      write.csv(mainout$interloadingz, paste(savedir,"interloadingsv.csv",sep=""))
      write.csv(mainout$intraloadingz, paste(savedir,"intraloadingsv.csv",sep=""))
    }
    
  }
  else{
  
  }
}

permcca2=function(bestpenaltyx,bestpenaltyz,pathx,pathz,k=1,permnums=500,maind=0){
  pacman::p_load(pacman,PMA)
  
  x <- read.csv(pathx,header = FALSE)
  z <- read.csv(pathz,header = FALSE)
  
  if (k==0){
    k=min(dim(x)[2],dim(z)[2])
  }
  
  if (sum(x[1,])==0){x=x[-1,]}
  if (sum(z[1,])==0){z=z[-1,]}
  shuffleind <- read.csv("C:/SF_data/suffleind.csv",header = TRUE)
  shuffleind=shuffleind[,2]
  x <- x[shuffleind,]
  z <- z[shuffleind,]
  
  N=nrow(z)
  ccperms <- matrix(nrow = permnums, ncol= 1)
  for (perm in c(1:permnums)){
    sampz <- sample(1:nrow(z))
    sampx <- sample(1:nrow(x))
    X<-x[sampx, ]
    Z<- z[sampz, ]
    out <- CCA(X, Z, typex = "standard", typez = "standard",
               niter = 30, K = k, penaltyx = bestpenaltyx,
               penaltyz = bestpenaltyz,trace = FALSE)
    indmax=which.max(out$cors)
    out$cors=out$cors[indmax]
    out$u=out$u[,indmax]
    out$v=out$v[,indmax]
    
    sdx <- apply(X, 2, sd)
    sdz <- apply(Z, 2, sd)
    X <- scale(X, TRUE, sdx)
    Z <- scale(Z, TRUE, sdz)
    
    if (mean(out$u == 0) != 1 && mean(out$v == 0) != 1) {
      ccperms[perm,] <- abs(cor(t(t(X)) %*% t(t(out$u)), t(t(Z)) %*% t(t(out$v))))
    }
    else {
      ccperms[perm,] <- 0
      print('here')
    }
  }
  pvals=length(which(ccperms>maind))/permnums
  # pvals <- apply(sweep(ccperms, 1, maind, "-") >= 0, 1, mean)
  return(pvals)
}

savedir=paste(maindir,"jacknife/",sep="")
temp=seq(0,1,by=0.1)
temp[0]=0.01
temp[11]=0.98
besttemp<- data.frame(expand.grid(bestpenalx = temp,bestpenalz = temp))

########jacknife of each sparsity
foreach (indxbst = seq_along(besttemp$bestpenalx))%dopar%{
  bestx=besttemp$bestpenalx[indxbst]
  bestz=besttemp$bestpenalz[indxbst]
  for(region in c(11,12)){
    pathf=paste("C:/SF_data/FCs/FC_",region,"_to_all.csv",sep="")
    savepath=paste(savedir,'region',region,sep="")
    if(!dir.exists(savepath)){dir.create(savepath)}
    pmacca_jacknif(bestpenaltyx=bestx,bestpenaltyz=bestz,pathx=pathf,pathz=paths,savepath,
                   partname=paste("/f",bestx,"_s",bestz,sep=""),k=1)
  }
}

######## main
savedir=paste(maindir,"perms/",sep="")
foreach (indxbst = seq_along(besttemp$bestpenalx))%dopar%{
  bestx=besttemp$bestpenalx[indxbst]
  bestz=besttemp$bestpenalz[indxbst]
  for(region in c(11,12)){
    pathf=paste("C:/SF_data/FCs/FC_",region,"_to_all.csv",sep="")
    savepath=paste(savedir,'region',region,sep="")
    if(!dir.exists(savepath)){dir.create(savepath)}
    pmacca_main(bestpenaltyx=bestx,bestpenaltyz=bestz,
                   pathx=pathf,pathz=paths,savepath,partname=paste("f",bestx,"_s",bestz,sep=""),
                   permnums=500,maincal=TRUE,otherloadings=FALSE)
  }
}

####### pvals
for(region in c(11,12)){
  # perms<- data.frame(expand.grid(bestpenalx = temp,bestpenalz = temp,pvals=1,CC=0))
  for (indxbst in seq_along(besttemp$bestpenalx)){
    pathf=paste("C:/SF_data/FCs/FC_",region,"_to_all.csv",sep="")
    bestx=besttemp$bestpenalx[indxbst]
    bestz=besttemp$bestpenalz[indxbst]
    partname=paste("f",bestx,"_s",bestz,sep="")
    
    savedir=paste(maindir,"jacknife/",sep="")
    savepath=paste(savedir,'region',region,'/',sep="")
    ds <-read.csv(paste(savepath,partname,"ds.csv",sep=""))[,2]
    print(paste("f",bestx,"_s",bestz,"meand:",mean(abs(ds))))
    
    savedir=paste(maindir,"perms/",sep="")
    savepath=paste(savedir,'region',region,sep="")
    maind <- read.csv(paste(savepath,'/main',partname,"cors.csv",sep=""))[,2]
    print(paste("f",bestx,"_s",bestz,"main:",maind))
    
    # besttemp$pvals[indxbst]=permcca2(bestpenaltyx=bestx,bestpenaltyz=bestz,pathx=pathf,pathz=paths,k=1,permnums=1000,maind=mean(abs(ds)))
    # print(paste("f",bestx,"_s",bestz,"***pval***:",mean(besttemp$pvals[indxbst])))
  }
  write.csv(besttemp,paste(savepath,"perms.csv",sep=""), row.names = FALSE)
}

region=12
bestx=0.7
bestz=0.8
pathf=paste("C:/SF_data/FCs/FC_",region,"_to_all.csv",sep="")
partname=paste("f",bestx,"_s",bestz,sep="")

savedir=paste(maindir,"jacknife/",sep="")
savepath=paste(savedir,'region',region,'/',sep="")
ds <-read.csv(paste(savepath,partname,"ds.csv",sep=""))[,2]
print(paste("f",bestx,"_s",bestz,"meand:",mean(abs(ds))))

savedir=paste(maindir,"perms/",sep="")
savepath=paste(savedir,'region',region,sep="")
maind <- read.csv(paste(savepath,'/main',partname,"cors.csv",sep=""))[,2]
print(paste("f",bestx,"_s",bestz,"main:",maind))



################bootstaraping to validate weights################
if(TRUE){
  library(PMA)
  library(foreach)
  library(iterators)
  library(parallel)
  library(doParallel)
  numCores <- detectCores()
  registerDoParallel(numCores-1)  # use multicore, set to the number of our cores
  perisylvians=c(11,12,13,14,29,30,63,64,65,66,79,80,81,82,85,86)
  maindir="C:/SF_data/data_byregions/nomyelin/permsdata/"
  paths="C:/SF_data/morphologicalfeats.csv"
  penalties <- read.csv("C:/SF_data/data_byregions/nomyelin/sparcityselection/penalties.csv",header = TRUE)
}

getfeatures <- function(path,S=FALSE){
  feat <- read.csv(path,header = FALSE)
  if (sum(feat[1,])==0){feat=feat[-1,]}
  # shuffleind <- read.csv("C:/SF_data/suffleind.csv",header = TRUE)
  # shuffleind=shuffleind[,2]
  if(S){
    initind=c(1,79,157,235,313,391)
    rmind=c(initind[1]:(initind[2]-1))
    feat <- feat[,-rmind]
  }
  # else{
  #   feat <- feat[shuffleind,]
  # }
  return(feat)
}

permcca=function(bestpenaltyx,bestpenaltyz,pathx,pathz,savepath,partname,permnums=500,otherloadings=FALSE){
  pacman::p_load(pacman,PMA)
  x <- getfeatures(path=pathx)
  z<- getfeatures(path=pathz,S=TRUE)
  
  k= 1
  N=nrow(z)
  featx=ncol(x)
  featz=ncol(z)
  corrs <- matrix(nrow = permnums, ncol= 1)
  us <- matrix(0,nrow = permnums, ncol= featx)
  vs <- matrix(0,nrow = permnums, ncol= featz)
  interloadingx <- matrix(nrow = permnums, ncol= featx)
  intraloadingx <- matrix(nrow = permnums, ncol= featx)
  interloadingz <- matrix(nrow = permnums, ncol= featz)
  intraloadingz <- matrix(nrow = permnums, ncol= featz)
  indperms <- read.csv("C:/SF_data/data_byregions/29-10/permsdata/permsind.csv",header = TRUE)
  
  for (perm in c(1:permnums)){
    ind <- sample(N,N,replace = TRUE)#unlist(indperms[perm,])
    Z <- z[ind,]
    X <- x[ind,]
    out <- CCA(X, Z, typex = "standard", typez = "standard",
               niter = 500, K = k, penaltyx = bestpenaltyx,
               penaltyz = bestpenaltyz,trace = FALSE)
    indmax=which.max(out$cors)
    corrs[perm,]=out$cors[indmax]
    us[perm,]=out$u[,indmax]
    vs[perm,]=out$v[,indmax]
    
    if(otherloadings){
      sdx <- apply(X, 2, sd)
      sdz <- apply(Z, 2, sd)
      X <- scale(X, TRUE, sdx)
      Z <- scale(Z, TRUE, sdz)
      
      markerx = t(t(X)) %*% (us[perm,])
      markerz = t(t(Z)) %*% (vs[perm,])
      
      for (i in c(1:featx)){
        interloadingx[perm,i] <- cor(markerz, X[,i])
        intraloadingx[perm,i] <- cor(markerx, X[,i])
      }
      for (i in c(1:featz)){
        interloadingz[perm,i] <- cor(markerx, Z[,i])
        intraloadingz[perm,i] <- cor(markerz, Z[,i])
      }
      
    }
    
  }
  
    savedir=paste(savepath,'/perm',perm,partname,sep="")
    write.csv(corrs, paste(savedir,"cors.csv",sep=""))
    write.csv(us, paste(savedir,"u.csv",sep=""))
    write.csv(vs, paste(savedir,"v.csv",sep=""))
    if(otherloadings){
      write.csv(interloadingx, paste(savedir,"interloadingsu.csv",sep=""))
      write.csv(intraloadingx, paste(savedir,"intraloadingsu.csv",sep=""))
      write.csv(interloadingz, paste(savedir,"interloadingsv.csv",sep=""))
      write.csv(intraloadingz, paste(savedir,"intraloadingsv.csv",sep=""))
    }
  
}

###################sampling
x <- read.csv(paths,header = TRUE)
N=dim(x)[1]
samplnum=5000
samples <- matrix(nrow = samplnum, ncol= N)
for (i in c(1:samplnum)){samples[i,]=sample(N,N,replace = TRUE)}
for (perm in c(1:samplnum)){
  for (ind in c(1:N)){
    if (!(ind %in% samples[perm,])){ print(paste(ind,'not in',perm))}
  }
}
write.csv(samples,paste(savedir,"permsind.csv",sep=""), row.names = FALSE)


#################apllying CCA_bootstrapping
foreach(region=c(11,12))%dopar%{
  bests=penalties[which(penalties[,1]==region),2:3]
  bestx=bests[1]
  bestz=bests[2]
  pathf=paste("C:/SF_data/FCs/FC_",region,"_to_all.csv",sep="")
  savepath=paste(maindir,region,'/',sep="")
  if(!dir.exists(savepath)){dir.create(savepath)}
  permcca(bestpenaltyx=bestx,bestpenaltyz=bestz,pathx=pathf,pathz=paths,
          savepath,partname='',permnums=1000,otherloadings=FALSE)
  
}

################visualization (data preparation)################
if(TRUE){
  library(ggplot2)
  library(circlize)
  # require(utils)
  # library(dplyr)
  # library(GenomicRanges)
  # library(ComplexHeatmap)
  perisylvians=c(11,12,13,14,29,30,63,64,65,66,79,80,81,82,85,86)
  maindir="C:/SF_data/data_byregions/nomyelin/permsdata/"
  paths="C:/SF_data/morphologicalfeats.csv"
  penalties <- read.csv("C:/SF_data/data_byregions/nomyelin/sparcityselection/penalties.csv",header = TRUE)
}

bootplotfunc <- function(org,boot,savename,region,fromplotcirc=FALSE,absactive=FALSE){
  sampnum <- dim(boot)[2]
  boot1 <- matrix(unlist(boot), nrow=dim(boot)[1], ncol=89)
  boot1 <- cbind(boot1[,1:region-1],matrix(0,nrow = dim(boot1)[1],ncol = 1),boot1[,region:sampnum])
  org1 <- rbind(matrix(org[1:region-1],nrow = region-1,ncol = 1),matrix(0,nrow = 1,ncol = 1),
                matrix(org[region:sampnum],nrow = sampnum-region+1,ncol = 1))
  if(fromplotcirc){
    out=list(org=org1,boot=boot1)
  }else{
    if (absactive){
      out=bootplot2(org=abs(org1),boot=abs(boot1),plotname=paste(savename,'abs',sep=""))
    }else{
      out=bootplot2(org=(org1),boot=(boot1),plotname=savename)
    }
  }
}
bootplotstruct <- function(org,boot,savename,fromplotcirc=FALSE,absactive=FALSE){
  boot <- matrix(unlist(boot), nrow=dim(boot)[1], ncol=dim(boot)[2] )
  strucnames=c('thick','sulc','curv','surf')
  # strucfeatind <- read.csv("C:/SF_data/Sroisnum.csv",header = FALSE,check.names = FALSE)
  roidel=c(37, 38, 41, 42, 71, 72, 73, 74, 75, 76, 77, 78)
  roi <-c(1:90)
  roi <- roi[-roidel]
  initind=c(1,79,157,235,313)
  
  for(i in c(1:(length(initind)-1))){
    bootemp <- matrix(0,dim(boot)[1],90)
    bootemp[,roi] <- boot[,initind[i]:(initind[i+1]-1)]
    orgtemp <- matrix(0,90)
    orgtemp[roi]<- org[initind[i]:(initind[i+1]-1)]
    if(fromplotcirc){
      assign(strucnames[i],list(org=orgtemp,boot=bootemp))
    }else{
      if (absactive){
        assign(strucnames[i],bootplot2(org=abs(orgtemp),boot=abs(bootemp),plotname=paste(savename,strucnames[i],'abs',sep='')))
      }else{ assign(strucnames[i],bootplot2(org=(orgtemp),boot=(bootemp),plotname=paste(savename,strucnames[i])))}
    }
  }
  out<-list(thick=thick, sulc=sulc, curv=curv, surf=surf)
}
bootstats2 <- function(org,bootdata, uplim,lwlim){
  # boot.stats <- sapply(seq_along(1:length(org)), function(i) org[i] - quantile(bootdata[,i] - org[i],c(uplim,lwlim),na.rm =T))
  boot.stats <- sapply(seq_along(1:length(org)), function(i)  quantile(bootdata[,i],c(uplim,lwlim),na.rm =T))
  boot.stats <- as.data.frame(t(boot.stats))
  colnames(boot.stats)<-c("low","high")
  zerocoef=(bootdata==0)*1
  boot.stats$ci <- sapply(seq_along(1:length(org)),function(i) sum(zerocoef[,i])/length(zerocoef[,i]))#(boot.stats$high - boot.stats$low)
  boot.stats$load <- org
  boot.stats$fea <- 1:length(org)
  boot.stats
}
bootplot2 <- function(org, boot,plotname,uplim=0.975,lwlim=0.025,plotsave=FALSE){
  # org <- sapply(seq_along(1:dim(boot)[2]), function(i) mean(boot[,i]))
  btst <- bootstats2(org,bootdata=boot,uplim=0.975,lwlim=0.025)
  cv <- apply(boot,2, function(x) sd(x) / mean(x) * 100)
  ind=is.na(cv)
  btst$cv <-sapply(seq_along(ind), function(x) if(ind[x]){if(sd(boot[,x])==0){0}else{500}}else{cv[x]})
  
  if (plotsave){
    ggplot(btst,aes(fea,load))+
      geom_point(aes(colour = low * high > 0)) +
      geom_errorbar(aes(ymax = high, ymin = low),  width=0.25) +
      theme_classic()
    ggsave(filename = paste(plotname,".tiff"), width = 10, height = 5, device='tiff', dpi=400)
  }
  btst
}

inval <- function(vec,svgname,sigvec,plotflag=TRUE){ 
  #pl:properloads
  invalidperms <- c(which(vec$low * vec$high <=0))
  mincheck= min(abs(vec$low),abs(vec$high)) > abs(vec$load)
  maxcheck= max(abs(vec$low),abs(vec$high)) < abs(vec$load)
  cvcheck=abs(vec$cv)>20
  invalidloads <- c(which( mincheck| maxcheck | cvcheck))
  sigvec[is.na(sigvec)]=0
  invalqvalue=which(sigvec!=1)
  invalids=c(invalidperms,invalidloads,invalqvalue)
  
  vec$load_shape=matrix('significant_loads',length(vec$low))
  vec$load_shape[which(mincheck+maxcheck>0)]='loads out of CI'
  vec$load_shape[which(sigvec==0)]='insignificant_loads'
  
  vec$CI_color=matrix('proper CI',length(vec$low))
  vec$CI_color[which(cvcheck)]= 'cv > 20'
  vec$CI_color[invalidperms]='zero in CI'
  if(plotflag){
    image=ggplot(vec,aes(x = fea, y = load, col = CI_color, group= load_shape)) +
      geom_point(aes(shape = load_shape)) +
      geom_errorbar(aes(ymax = high, ymin = low),  width=0.25) +
      scale_shape_manual(values=c(4, 20,1))+
      # scale_color_manual(values = c("#ca7dcc","#18a558","#ff0000"))
      scale_color_manual(values = c("#ca7dcc","#FF0000","#18a558")) +
      theme(text = element_text(size = 20))  
    ggsave(file=svgname, plot=image, width=12, height=8)
  }
  return(unique(invalids))
}
ttestweights <- function(strucvar2,funcvar2){
  pval=matrix(100,nrow = 5,ncol = 90)
  for(i in c(1:90)){
    pval[1,i]=t.test(strucvar2$thick$boot[,i],mu=strucvar2$thick$org[i,1])$p.value
    pval[2,i]=t.test(strucvar2$sulc$boot[,i],mu=strucvar2$sulc$org[i,1])$p.value
    pval[3,i]=t.test(strucvar2$curv$boot[,i],mu=strucvar2$curv$org[i,1])$p.value
    pval[4,i]=t.test(strucvar2$surf$boot[,i],mu=strucvar2$surf$org[i,1])$p.value
    pval[5,i]=t.test(funcvar2$boot[,i],mu=funcvar2$org[i,1])$p.value
  }
  qvalueS=matrix(p.adjust(c(pval[1:4,]), method="BH"),nrow = 4,ncol=90)
  qvalueF=p.adjust(c(pval[5,]), method="BH")
  qvalues=matrix(p.adjust(c(pval), method="BH"),nrow = 5,ncol=90)
  sig_qvalueS=(qvalueS<0.05)*1
  sig_qvalueF=(qvalueF<0.05)*1
  sig_qvalue=(qvalues<0.05)*1
  return(sig_qvalue)
}


plotcirc2 <- function(orgvecs,u_boost,v_boost,savepath,region,absactive=FALSE){
  
  nameinfo="CI95CV20"
  if(absactive){
    svgname=paste(savepath,region,nameinfo,"_abs",sep="")
  }else{
    svgname=paste(savepath,region,nameinfo,sep="")}
  
  funcvar=bootplotfunc(org=orgvecs$u,boot=u_boost,savename=paste(savepath,"functions"),region = region,absactive = absactive)
  funcvar2=bootplotfunc(org=orgvecs$u,boot=u_boost,savename=paste(savepath,"functions"),region = region,fromplotcirc=TRUE,absactive = absactive)
  strucvar=bootplotstruct(org=orgvecs$v,boot=v_boost,savename=paste(savepath),absactive = absactive)
  strucvar2=bootplotstruct(org=orgvecs$v,boot=v_boost,savename=paste(savepath),fromplotcirc=TRUE,absactive = absactive)
  
  sig_qvalue=ttestweights(strucvar2,funcvar2)
  
  th <- strucvar$thick
  sul <- strucvar$sulc
  cur <- strucvar$curv
  sur <- strucvar$surf
  
  thind=inval(th,paste(svgname,'thickloads.svg',sep=""),sigvec = sig_qvalue[1,])
  sulind=inval(sul,paste(svgname,'sulcloads.svg',sep=""),sigvec = sig_qvalue[2,])
  curind=inval(cur,paste(svgname,'curvloads.svg',sep=""),sigvec = sig_qvalue[3,])
  surind=inval(sur,paste(svgname,'surfloads.svg',sep=""),sigvec = sig_qvalue[4,])
  
  thick=strucvar2$thick$boot
  thick[,thind]=matrix(0,dim(thick)[1],)
  sulc=strucvar2$sulc$boot
  sulc[,sulind]=matrix(0,dim(sulc)[1],)
  curv=strucvar2$curv$boot
  curv[,curind]=matrix(0,dim(curv)[1],)
  surf=strucvar2$surf$boot
  surf[,surind]=matrix(0,dim(surf)[1],)
  
  funcind=inval(funcvar,paste(svgname,'funcloads.svg',sep=""),sigvec = sig_qvalue[5,])
  funcin=funcvar2$boot
  funcin[,funcind]=matrix(0,dim(funcin)[1],)
  ind=which(abs(colMeans(funcin))!=0)
  
  
  sectorsvec=c(c(1:45)*2,c(44:0)*2+1)
  pos=matrix(0.5,dim(thick)[2],)
  bgcolors = rand_color(45)

  svg(paste(svgname,'S_F.svg',sep=""))
  
  circos.par("track.height" = 0.1,start.degree = 90)
  circos.initialize(sectorsvec, xlim = c(0, 0.5))
  
  circos.track(ylim = c(-0.2, 0.2), panel.fun = function(x, y) {
    pos=get.cell.meta.data("xcenter")
    l = c(1:90) == CELL_META$sector.index
    circos.boxplot(thick[,l],get.cell.meta.data("xcenter"),cex = 0.1)
  },bg.border='burlywood3')
  
  circos.track(ylim = c(-0.2, 0.2), panel.fun = function(x, y) {
    pos=get.cell.meta.data("xcenter")
    l = c(1:90) == CELL_META$sector.index
    circos.boxplot(sulc[,l],pos,cex = 0.2)
  },bg.border='cornsilk4')
  
  circos.track(ylim = c(-0.2, 0.2), panel.fun = function(x, y) {
    pos=get.cell.meta.data("xcenter")
    l = c(1:90) == CELL_META$sector.index
    circos.boxplot(curv[,l],pos,cex = 0.2)
  },bg.border='bisque3')
  
  circos.track(ylim = c(-0.2, 0.2), panel.fun = function(x, y) {
    pos=get.cell.meta.data("xcenter")
    l = c(1:90) == CELL_META$sector.index
    circos.boxplot(surf[,l],pos,cex = 0.2)
  },bg.border='azure3')
  
  if(length(funcind)<90){
    loads=funcvar$load[ind]
    loads=loads[loads!=0]
    
    minval=min(loads)
    maxval=max(loads)
    step=(maxval-minval)/4
    rangeload=seq(minval,maxval,by=step)
    # col_fun = colorRamp2(rangeload, c("yellow","green","chartreuse","cyan","blue","darkblue","darkorchid","deeppink","deeppink4","red"))
    
    if(maxval<0){
      col_fun = colorRamp2(rangeload, c("deeppink4","blue","cyan","green","yellow"))
    }else if(minval>0){
      col_fun = colorRamp2(rangeload, c("yellow","green","cyan","blue","deeppink4"))
    }else{
      col_fun = colorRamp2(rangeload, c("deeppink4","yellow","cyan","green","blue"))
    }
    
    for (i in ind){
      circos.link(region, 0.5, i, 0.5,
                  col = add_transparency(col_fun(funcvar$load[i])))
    }
    lgd_links = Legend(at = round(c(minval,minval+(maxval-minval)/2,maxval),digits = 3), col_fun = col_fun,
                       title_position = "topleft", title = "FCs")
    lgd_list_vertical = packLegend(lgd_links)
    draw(lgd_list_vertical, x = unit(7, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
  }
  circos.clear()
  dev.off()
}


mainvecs <-function(maindir){
  dsmain <-read.csv(paste(maindir,"cors.csv",sep=""))[,2]
  usmain <-read.csv(paste(maindir,"u.csv",sep=""))[,-1]
  vsmain <-read.csv(paste(maindir,"v.csv",sep=""))[,-1]
  out <-list(v=vsmain,u=usmain,d=dsmain)
  return(out)

}

mainweightsdir="C:/SF_data/data_byregions/nomyelin/sparcityselection/perms/"
roinds=matrix(0,nrow= 4,ncol=90)
roiloads=matrix(0,nrow=2,ncol=90) #
for(region in c(11,12,29, 30, 63, 79, 81,82)){
  bests=penalties[which(penalties[,1]==region),2:3]
  bestx=bests[1]
  bestz=bests[2]
  pathf=paste("C:/SF_data/FCs/FC_",region,"_to_all.csv",sep="")
  orgvecs=mainvecs(maindir=paste(mainweightsdir,'region',region,'/mainf',bestx,'_s',bestz,sep=""))
  
  u_boost <- read.csv(paste(maindir,region,"/perm1000u.csv",sep=""))[,-1]
  v_boost <- read.csv(paste(maindir,region,"/perm1000v.csv",sep=""))[,-1]
  savepath=paste(maindir,region,'/visig/',sep="")
  if(!dir.exists(savepath)){dir.create(savepath)}
  # plotcirc_basevalues (orgvecs,u_boost,v_boost,savepath=savepath,
  #                      region,absactive=FALSE,pathf,paths)
  plotcirc2 (orgvecs,u_boost,v_boost,savepath,region,absactive = FALSE)
  # plotcircnets(orgvecs,u_boost,v_boost,savepath=paste(savepath,partname,sep=''),region,absactive = FALSE)
  # roijointinfos=plotbrainregion_pre (orgvecs,u_boost,v_boost,savepath=paste(savepath,partname,sep=''),
  #                             region,prevalue = roinds,prevals = roiloads,absactive = FALSE,
  #                             corrplotflag=TRUE,pathx=pathf,pathz=paths)
  # roinds=roijointinfos$inds
  # roiloads=roijointinfos$vals
  # 
  # plotcorrsofresiduals (orgvecs,u_boost,v_boost,savepath=paste(savepath,partname,sep=''),
  #                       region,pathx=pathf,pathz=paths)
  
  # residulSFmap(u_boot,v_boot,region,savepath)

}

write.csv(roinds,"C:/SF_data/data_byregions/29-10/permsdatall/absjointroivalue.csv")
write.csv(roiloads,"C:/SF_data/data_byregions/29-10/permsdatall/absjointroiloads.csv")









sig_regions=c(11,12,13,14,29,30,63,64,65,66,80,79,81,82,85,86)

partname='first'
partname2=''

region=86

for(region in sig_regions){
  mainfile=paste(savedirmain,'/',region,sep="")
  maindir=mainfile
  pathf=paste("C:/SF_data/FCs/FC_",region,"_to_all.csv",sep="")
  savedirmain="C:/SF_data/data_byregions/29-10/permsdatall/main/"
  saveperm="C:/SF_data/data_byregions/29-10/permsdatall/"
  paths="C:/SF_data/morphologicalfeats.csv"
  dsmain <-read.csv(paste(maindir,"cors.csv",sep=""))[,2]
  usmain <-read.csv(paste(maindir,"u.csv",sep=""))[,-1]
  vsmain <-read.csv(paste(maindir,"v.csv",sep=""))[,-1]
  
  match=1
  usmain <- usmain[,match]
  vsmain <- vsmain[,match]
  dsmain <- dsmain[match]
  
  pathx=pathf
  pathz=paths
  
  X <- read.csv(pathx,header = FALSE)
  if (sum(X[1,])==0){X=X[-1,]}
  Z <- read.csv(pathz,header = FALSE) 
  if (sum(Z[1,])==0){Z=Z[-1,]}
  sdx <- apply(X, 2, sd)
  sdz <- apply(Z, 2, sd)
  X <- scale(X, TRUE, sdx)
  Z <- scale(Z, TRUE, sdz)
  
  functional_marker = t(t(X)) %*% (usmain)
  structural_marker = t(t(Z)) %*% (vsmain)
  corrs <- cor.test(structural_marker, functional_marker)
  
  plotname=paste(mainfile,'mainweights',partname,partname2,'.svg',sep = "")
  svg(plotname)
  plot(functional_marker, structural_marker, pch = 19, col = "lightblue")# Creating the plot
  abline(lm(structural_marker ~ functional_marker), col = "red", lwd = 3)# Regression line
  text(paste("Correlation:", round(corrs$estimate, 2)), x = 5, y = 4)# Pearson correlation
  text(paste("p_value:", round(corrs$p.value, 4)), x = 5, y = 0)# pvalue
  dev.off()
  
  
  
  mainfile=paste(savedirmain,'/',region,sep="")
  pathf=paste("C:/SF_data/FCs/FC_",region,"_to_all.csv",sep="")
  orgvecs=mainvecs(maindir=mainfile,selectype=partname,outbootype=partname2,pathx=pathf,pathz=paths,
                   plotorrs=TRUE,plotname=paste(mainfile,'mainweights',partname,partname2,'.png',sep = ""))
  
  savepath=paste(saveperm,region,'/visig',partname2,'/',sep="")
  if(!dir.exists(savepath)){dir.create(savepath)}
  
  u_boost <- read.csv(paste(saveperm,region,'/',partname,partname2,"u.csv",sep=""))[,-1]
  v_boost <- read.csv(paste(saveperm,region,'/',partname,partname2,"v.csv",sep=""))[,-1]
  
  nameinfo="CI95CV20"
  absactive=FALSE
  funcvar=bootplotfunc(org=orgvecs$u,boot=u_boost,savename=paste(savepath,"functions"),region = region,absactive = absactive)
  funcvar2=bootplotfunc(org=orgvecs$u,boot=u_boost,savename=paste(savepath,"functions"),region = region,fromplotcirc=TRUE,absactive = absactive)
  strucvar=bootplotstruct(org=orgvecs$v,boot=v_boost,savename=paste(savepath),absactive = absactive)
  strucvar2=bootplotstruct(org=orgvecs$v,boot=v_boost,savename=paste(savepath),fromplotcirc=TRUE,absactive = absactive)
  
  sig_qvalue=ttestweights(strucvar2,funcvar2)
  svgname=paste(savepath,region,nameinfo,sep="")
  my <- strucvar$myelin
  myind=inval(my,paste(svgname,'myelinloads.svg',sep=""),sigvec = sig_qvalue[1,])
  
  myelin=strucvar2$myelin$boot
  myelin[,myind]=matrix(0,dim(myelin)[1],)
  
  myelin_org=strucvar2$myelin$org
  myelin_org[myind]=0
  
  funcind=inval(funcvar,paste(svgname,'funcloads.svg',sep=""),sigvec = sig_qvalue[6,])
  funcin=funcvar2$boot
  funcin[,funcind]=matrix(0,dim(funcin)[1],)
  
  funcin_org=funcvar2$org
  funcin_org[funcind]=0
  
  pathx=pathf
  pathz=paths
  X <- read.csv(pathx,header = FALSE)
  Z <- read.csv(pathz,header = FALSE)  
  sdx <- apply(X, 2, sd)
  sdz <- apply(Z, 2, sd)
  X <- scale(X, TRUE, sdx)
  Z <- scale(Z, TRUE, sdz)
  
  roidel=c(37, 38, 41, 42, 71, 72, 73, 74, 75, 76, 77, 78)
  roi <-c(1:90)
  roi <- roi[-roidel]
  initind=c(1,79,157,235,313,391)
  
  functional_marker = t(t(X)) %*% (funcin_org[-region])
  structural_marker = t(t(Z[,initind[1]:(initind[2]-1)])) %*% (myelin_org[-roidel])
  corrs <- cor.test(structural_marker, functional_marker)
  
  plotname=paste(savepath,region,nameinfo,'residweights.svg',sep="")
  ypos=mean(structural_marker)
  svg(plotname)
  plot(functional_marker, structural_marker, pch = 19, col = "lightblue")# Creating the plot
  abline(lm(structural_marker ~ functional_marker), col = "red", lwd = 3)# Regression line
  text(paste("Correlation:", round(corrs$estimate, 2)), x = 0.5, y = 1)# Pearson correlation
  text(paste("p_value:", round(corrs$p.value, 4)), x = 0.5, y = -0.5)# Pearson correlation
  dev.off()
}

library("corrplot")
for(region in sig_regions){
  mainfile=paste(savedirmain,'/',region,sep="")
  maindir=mainfile
  pathf=paste("C:/SF_data/FCs/FC_",region,"_to_all.csv",sep="")
  savedirmain="C:/SF_data/data_byregions/29-10/permsdatall/main/"
  saveperm="C:/SF_data/data_byregions/29-10/permsdatall/"
  paths="C:/SF_data/morphologicalfeats.csv"
  
  pathx=pathf
  pathz=paths
  X <- read.csv(pathx,header = FALSE)
  if (sum(X[1,])==0){X=X[-1,]}
  Z <- read.csv(pathz,header = FALSE) 
  if (sum(Z[1,])==0){Z=Z[-1,]}
  sdx <- apply(X, 2, sd)
  sdz <- apply(Z, 2, sd)
  X <- scale(X, TRUE, sdx)
  Z <- scale(Z, TRUE, sdz)
  
  corval=matrix(0,nrow = dim(X)[2],ncol=78)
  pval=matrix(100,nrow = dim(X)[2],ncol=78)
  for (featx in c(1:dim(X)[2])){
    functional_marker = X[,featx]
    assign(paste('S_F',featx),matrix)
    for (featz in c(initind[5]:(initind[6]-1))){
      structural_marker = Z[,featz]
      corrs <- cor.test(structural_marker, functional_marker)
      corval[featx,featz-(initind[5]-1)]=corrs$estimate
      pval[featx,featz-(initind[5]-1)]=corrs$p.value
    }
  }
  colnames(corval)=c(1:90)[-roidel]
  rownames(corval)=c(1:90)[-region]
  # colnames(pval)=c(1:90)[-roidel]
  # rownames(pval)=c(1:90)[-region]
  nosigvals=which(pval>0.05)
  corval[nosigvals]=0
  
  # strucnames=c('myelin','thick','sulc','curv','surf')
  svgname=paste(mainfile,'solocorr_FSurf',partname,partname2,'.svg',sep = "")
  svg(svgname)
  corrplot(round(corval,2),method = "number",number.cex=0.2,tl.cex=0.3,is.corr=FALSE,col.lim = c(min(corval)-0.1,max(corval)+0.1))
  dev.off()
  
}



