# G:/SF9_2024/CCA/release2_com/upload/FALFF/funcs_FALFF_PSR.R

get_data <-function(xread,zread,featid,selectedrois){
  x <- matrix(unlist(xread),nrow=dim(xread)[1],ncol=dim(xread)[2])
  z <- matrix(unlist(zread),nrow=dim(zread)[1],ncol=dim(zread)[2])
  newx=x[,selectedrois]
  # if(i==6){z<-z}else{z <- z[,initind[i]:(initind[i+1]-1)]}
  initind=c(1,79,157,235,313,391)
  
  roidel=c(37, 38, 41, 42, 71, 72, 73, 74, 75, 76, 77, 78)
  roilist <-c(1:90)
  roilist <- roilist[-roidel]      
  selectedrois_instruc=c()
  for (j in selectedrois){selectedrois_instruc=c(selectedrois_instruc,which(roilist==j))}
  newz=c()
  for (ii in c(1:5)){
    tempz=z[,initind[ii]:(initind[ii+1]-1)]
    newz=cbind(newz,tempz[,selectedrois_instruc])
  }
  return(list(x=newx,z=newz))
}

###################feature_selection(SCCA)###############
########Permutation
permcca_parallel <- function(penalsx, penalsz, x, z, k = 1, permnums = 500) {
  pacman::p_load(pacman, PMA)
  ftrans <- function(x) { return(0.5 * log((1 + x) / (1 - x))) }
  
  if (k == 0) { k = min(dim(x)[2], dim(z)[2]) }
  N = nrow(z)
  pvals = c()
  zstats = c()
  cors = c()
  penalsnum = length(penalsx)
  
  # Set up parallel backend
  no_cores <- detectCores() - 2
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  set.seed(42)
  
  for (i in 1:penalsnum) {
    penalx = penalsx[i]
    penalz = penalsz[i]
    
    ccperms <- foreach(perm = 1:permnums, .combine = rbind, .packages = 'PMA') %dopar% {
      set.seed(perm+7)
      sampz <- sample(1:nrow(z))
      sampx <- sample(1:nrow(x))
      X <- x[sampx, ]
      Z <- z[sampz, ]
      out <- CCA(X, Z, typex = "standard", typez = "standard",
                 niter = 30, K = k, penaltyx = penalx,
                 penaltyz = penalz, trace = FALSE)
      
      if (mean(out$u == 0) != 1 && mean(out$v == 0) != 1) {
        return(out$cors)
      } else {
        return(0)
      }
    }
    
    out <- CCA(x, z, typex = "standard", typez = "standard",
               niter = 30, K = k, penaltyx = penalx,
               penaltyz = penalz, trace = FALSE)
    maincor = out$cors
    pvals = c(pvals, length(which(ccperms >= maincor)) / permnums)
    cc.norm <- ftrans(maincor)
    ccperm.norm <- ftrans(ccperms)
    zstats <- c(zstats, abs(cc.norm - mean(ccperm.norm)) / sd(ccperm.norm))
    cors = c(cors, maincor)
  }
  
  stopCluster(cl)
  
  indmax = which.max(zstats)
  
  return(list(pvals = pvals, zstats = zstats, pvalbestz = pvals[indmax],
              bestpenaltyx = penalsx[indmax], bestpenaltyz = penalsz[indmax], cors = cors))
}

########main
pmacca_main=function(bestpenaltyx,bestpenaltyz,x,z,savepath,partname){
  pacman::p_load(pacman,PMA)
  
  k= min(dim(x)[2],dim(z)[2])
  xcrossz_list=matrix(nrow=dim(x)[2],ncol=k)
  zcrossx_list=matrix(nrow=dim(z)[2],ncol=k)
  zloading_list=matrix(nrow=dim(z)[2],ncol=k)
  xloading_list=matrix(nrow=dim(x)[2],ncol=k)
  
  mainout <- PMA::CCA(x, z, typex = "standard", typez = "standard",niter = 500, K = k,
                      penaltyx = bestpenaltyx,penaltyz = bestpenaltyz,trace=FALSE)
  
  savedir=paste(savepath,'/main',partname,sep="")
  write.csv(mainout$cors, paste(savedir,"cors.csv",sep=""))
  write.csv(mainout$u, paste(savedir,"u.csv",sep=""))
  write.csv(mainout$v, paste(savedir,"v.csv",sep=""))
  sdx <- apply(x, 2, sd)
  sdz <- apply(z, 2, sd)
  X <- scale(x, TRUE, TRUE)
  Z <- scale(z, TRUE, TRUE)
  
  for (i in 1:k) {
    xcrossz_list[,i] <- cor(X, Z%*%mainout$v[, i])
    zcrossx_list[,i]<- cor(Z, X%*%mainout$u[, i])
    zloading_list[,i] <- cor(Z, Z%*%mainout$v[, i])
    xloading_list[,i] <- cor(X, X%*%mainout$u[, i])
  }
  write.csv(xcrossz_list, paste(savedir,"xcross.csv",sep=""))
  write.csv(zcrossx_list, paste(savedir,"zcrosS.csv",sep=""))
  write.csv(zloading_list, paste(savedir,"zloading.csv",sep=""))
  write.csv(xloading_list, paste(savedir,"Xloading.csv",sep=""))
  
}

#######significant canonical components (permutation)
permcca2_parallel2 <- function(penalsx, penalsz, x, z, k = 1, permnums = 500) {
  pacman::p_load(pacman, PMA)
  ftrans <- function(x) { return(0.5 * log((1 + x) / (1 - x))) }
  
  if (k == 0) { k = min(dim(x)[2], dim(z)[2]) }
  N = nrow(z)
  pvals = c()
  zstats = c()
  cors = c()
  penalsnum = length(penalsx)
  
  # Set up parallel backend
  no_cores <- detectCores() - 2
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  for (i in 1:penalsnum) {
    penalx = penalsx[i]
    penalz = penalsz[i]
    ccperms <- foreach(perm = 1:permnums, .combine = rbind, .packages = 'PMA') %dopar% {
      set.seed(perm+7)
      sampz <- sample(1:nrow(z))
      sampx <- sample(1:nrow(x))
      X <- x[sampx, ]
      Z <- z[sampz, ]
      out <- CCA(X, Z, typex = "standard", typez = "standard",
                 niter = 30, K = k, penaltyx = penalx,
                 penaltyz = penalz, trace = FALSE)
      
      if (mean(out$u == 0) != 1 && mean(out$v == 0) != 1) {
        return(out$cors)
      } else {
        return(0)
      }
    }
  }
  stopCluster(cl)
  
  return(ccperms)
}

#######bootstraping stable canonical weights
Bootscca <- function(bestpenaltyx, bestpenaltyz, x, z, k=0,savepath, partname, iternums = 500, otherloadings = FALSE) {
  pacman::p_load(pacman,PMA)
  if (k==0){k <- min(ncol(x),ncol(z))}
  N <- nrow(z)
  featx <- ncol(x)
  featz <- ncol(z)
  
  # Initialize lists to store matrices
  corrs_list <- vector("list", k)
  us_list <- vector("list", k)
  vs_list <- vector("list", k)
  set.seed(7)
  
  for (i in 1:k) {
    corrs_list[[i]] <- matrix(nrow = iternums, ncol = 1)
    us_list[[i]] <- matrix(nrow = iternums, ncol = featx)
    vs_list[[i]] <- matrix(nrow = iternums, ncol = featz)
  }
  scoresX=matrix(nrow=iternums,ncol = k)
  scoresZ=matrix(nrow=iternums,ncol = k)
  for (perm in 1:iternums) {
    ind <- sample(N, N, replace = TRUE)
    Z <- z[ind, ]
    X <- x[ind, ]
    out <- CCA(X, Z, typex = "standard", typez = "standard",
               niter = 500, K = k, penaltyx = bestpenaltyx,
               penaltyz = bestpenaltyz, trace = FALSE)
    
    sdx <- apply(X, 2, sd)
    sdz <- apply(Z, 2, sd)
    X <- scale(X, TRUE, sdx)
    Z <- scale(Z, TRUE, sdz)
    
    for (i in 1:k) {
      corrs_list[[i]][perm, ] <- out$cors[i]
      us_list[[i]][perm, ] <- out$u[, i]
      vs_list[[i]][perm, ] <- out$v[, i]
    }
  }
  
  for (i in 1:k) {
    savedir <- paste0(savepath, '/perm', i, partname)
    write.csv(corrs_list[[i]], paste0(savedir, "cors_comp", i, ".csv"))
    write.csv(us_list[[i]], paste0(savedir, "u_comp", i, ".csv"))
    write.csv(vs_list[[i]], paste0(savedir, "v_comp", i, ".csv"))
  }
  
}

###############visualization(selected features)
bootstats2 <- function(org,bootdata, uplim,lwlim){
    boot.stats <- sapply(seq_along(1:length(org)), function(i) quantile(bootdata[,i],c(uplim,lwlim),na.rm =T))
    boot.stats <- as.data.frame(t(boot.stats))
    colnames(boot.stats)<-c("low","high")
    boot.stats$ci <- abs(boot.stats$high - boot.stats$low)
    boot.stats$load <- sapply(seq_along(1:length(org)), function(i) mean(bootdata[,i]))
    boot.stats$fea <- 1:length(org)
    boot.stats$org <- org[,1]
    cv <- apply(bootdata,2, function(x) sd(x) / mean(x) * 100)
    ind=is.na(cv)
    boot.stats$cv <-sapply(seq_along(ind), function(x) if(ind[x]){if(sd(bootdata[,x])==0){0}else{500}}else{cv[x]})
    boot.stats
  }
  
bootplot2 <- function(org, boot,plotname,uplim=0.975,lwlim=0.025,plotsave=FALSE){
    btst <- bootstats2(org,bootdata=boot,uplim=0.975,lwlim=0.025)
    if (plotsave){
      ggplot(btst,aes(fea,load))+
        geom_point(aes(colour = low * high > 0)) +
        geom_errorbar(aes(ymax = high, ymin = low),  width=0.25) +
        theme_classic()
      ggsave(filename = paste(plotname,".tiff"), width = 10, height = 5, device='tiff', dpi=400)
    }
    btst
  }
  
bootdata4allregions<- function(org,boot,struc_data,info){
    boot=(sapply(seq_along(1:length(org)), function(i) org[i]-(boot[,i] - org[i])))
    if(struc_data){
      boot <- matrix(unlist(boot), nrow=dim(boot)[1], ncol=dim(boot)[2])
      strucnames=info
      roi <-c(11,12,13,14,29,30,63,64,65,66,79,80,81,82,85,86)
      
      if(length(org)>16){
        initind=c(1,17,2*16+1,3*16+1,4*16+1,5*16+1)
        for(i in c(1:(length(initind)-1))){
          bootemp <- matrix(0,dim(boot)[1],16)
          bootemp <- boot[,initind[i]:(initind[i+1]-1)]
          # orgtemp <- matrix(0,16)
          orgtemp<- matrix(unlist(org[initind[i]:(initind[i+1]-1)]), ncol=16)
          assign(strucnames[i],list(org=orgtemp,boot=bootemp))
        }
        out<-list(myelin=myelin, thick=thick, sulc=sulc, curv=curv, surf=surf)
        
      }else{
        bootemp <- matrix(0,dim(boot)[1],16)
        bootemp[,roi] <- boot
        # orgtemp <- matrix(0,16)
        orgtemp<- matrix(unlist(org), ncol=16)
        out=list(org=orgtemp,boot=bootemp)
      }
      
    }else{
      sampnum <- dim(boot)[2]
      boot1 <- matrix(unlist(boot), nrow=dim(boot)[1], ncol=16)
      org1 <- matrix(unlist(org), ncol=16)
      out=list(org=org1,boot=boot1)
    }
    return(out)
  }
  
plotcirc3 <- function(orgvecs,u_boost,v_boost,savepath,region=NULL){
    
    nameinfo="CI95_visall_"
    svgname=paste(savepath,nameinfo,sep="")
    
    struc.feat.info=c('myelin', 'thick', 'sulc', 'curv', 'surf')
    U_vec=bootdata4allregions(org=orgvecs$u,boot=u_boost,struc_data=FALSE,info=NULL)
    V_vec=bootdata4allregions(org=orgvecs$v,boot=v_boost,struc_data=TRUE,info=struc.feat.info)
    
    funcvar=bootplot2(org=abs(U_vec$org),boot=abs(U_vec$boot),plotname=paste(savepath,"functions"))
    
    curv=bootplot2(org=abs(V_vec$curv$org),boot=abs(V_vec$curv$boot),plotname=paste(savename,"curv"))
    thick=bootplot2(org=abs(V_vec$thick$org),boot=abs(V_vec$thick$boot),plotname=paste(savename,"thick"))
    myelin=bootplot2(org=abs(V_vec$myelin$org),boot=abs(V_vec$myelin$boot),plotname=paste(savename,"myelin"))
    sulc=bootplot2(org=abs(V_vec$sulc$org),boot=abs(V_vec$sulc$boot),plotname=paste(savename,"sulc"))
    surf=bootplot2(org=abs(V_vec$surf$org),boot=abs(V_vec$surf$boot),plotname=paste(savename,"surf"))
    
    mainloads=list(func=funcvar$load,myelin=myelin$load,thick=thick$load,curv=curv$load,sulc=sulc$load,surf=surf$load)
    mainpower=list(func=sum((funcvar$load)^2),myelin=sum((myelin$load)^2),thick=sum((thick$load)^2),curv=sum((curv$load)^2),sulc=sum((sulc$load)^2),surf=sum((surf$load)^2))
    strucloads=list(myelin=myelin$load,thick=thick$load,curv=curv$load,sulc=sulc$load,surf=surf$load)
    
    sig_qvalue=ttestweights(strucvar=V_vec,funcvar=U_vec)
    
    myind=inval(vec=myelin,paste(svgname,'myelinloads.svg',sep=""),sigvec = sig_qvalue[1,])
    myelin$load[myind] = matrix(0,length(myind))
    
    thind=inval(thick,paste(svgname,'thickloads.svg',sep=""),sigvec = sig_qvalue[2,])
    thick$load[thind] = matrix(0,length(thind))
    
    sulind=inval(sulc,paste(svgname,'sulcloads.svg',sep=""),sigvec = sig_qvalue[3,])
    sulc$load[sulind] = matrix(0,length(sulind))
    
    curind=inval(curv,paste(svgname,'curvloads.svg',sep=""),sigvec = sig_qvalue[4,])
    curv$load[curind] = matrix(0,length(curind))
    
    surind=inval(surf,paste(svgname,'surfloads.svg',sep=""),sigvec = sig_qvalue[5,])
    surf$load[surind] = matrix(0,length(surind))
    
    funcind=inval(funcvar,paste(svgname,'funcloads.svg',sep=""),sigvec = sig_qvalue[6,])
    funcvar$load[funcind]=matrix(0,length(funcind))
    
    sigloads=list(func=funcvar$load,myelin=myelin$load,thick=thick$load,curv=curv$load,sulc=sulc$load,surf=surf$load)
    sigpower=list(func=sum((funcvar$load)^2),myelin=sum((myelin$load)^2),thick=sum((thick$load)^2),curv=sum((curv$load)^2),sulc=sum((sulc$load)^2),surf=sum((surf$load)^2))
    
    strucloads=list(myelin=myelin$load,thick=thick$load,curv=curv$load,sulc=sulc$load,surf=surf$load)
    
    return(out=list(mainloads=mainloads,sigloads=sigloads,mainpower=mainpower,sigpower=sigpower))
  }
  
plotcirc4 <- function(orgvecs,u_boost,v_boost,savepath,region=NULL){
    nameinfo="CI95_vis"
    svgname=paste0(savepath,region,nameinfo)
    
    struc.feat.info=strsplit(savepath,"/")[[1]]
    struc.feat.info=struc.feat.info[length(struc.feat.info)-3]
    U_vec=bootdata4allregions(org=orgvecs$u,boot=u_boost,struc_data=FALSE,info=region)
    V_vec=bootdata4allregions(org=orgvecs$v,boot=v_boost,struc_data=TRUE,info=struc.feat.info)
    
    funcvar=bootplot2(org=U_vec$org,boot=U_vec$boot,plotname=paste(savepath,"functions"))
    strucvar=bootplot2(org=V_vec$org,boot=V_vec$boot,plotname=paste(savepath,struc.feat.info))
    
    mainloads=list(func=funcvar$load,struc=strucvar$load)
    mainpower=list(func=sum((funcvar$load)^2),struc=sum((strucvar$load)^2))
    
    cirplotfunc(svgname=paste0(svgname,'_',struc.feat.info,'_FC_main.svg'),strucloads=strucvar$load,funcloads=funcvar$load,region=region,plotmode='singlestruc')
    
    
    sig_qvalue=ttestweights(strucvar=V_vec,funcvar=U_vec,solofeat = TRUE)
    
    strucind=inval(vec=strucvar,paste0(svgname,struc.feat.info,'loads.svg'),sigvec = sig_qvalue[1,],plotflag = FALSE)
    strucvar$load[strucind] = matrix(0,length(strucind))
    
    funcind=inval(funcvar,paste0(svgname,'funcloads.svg'),sigvec = sig_qvalue[2,],plotflag = FALSE)
    funcvar$load[funcind]=matrix(0,length(funcind))
    
    cirplotfunc(svgname=paste0(svgname,'_',struc.feat.info,'_FC_sig.svg'),strucloads=strucvar$load,funcloads=funcvar$load,region=region,plotmode='singlestruc')
    
    sigpower=list(func=funcvar$load,struc=strucvar$load)
    sigloads=list(func=sum((funcvar$load)^2),struc=sum((strucvar$load)^2))
    
    return(out=list(mainloads=mainloads,sigloads=sigloads,mainpower=mainpower,sigpower=sigpower))
  }
  
inval <- function(vec,svgname,sigvec,plotflag=TRUE,cvth=30){ #pl:properloads
    cvth=3000000
    invalidperms <- c(which(vec$low * vec$high <=0)) #confidence intervals overlapped with zero
    mincheck= min(abs(vec$low),abs(vec$high)) > abs(vec$load) #loads out of CI
    maxcheck= max(abs(vec$low),abs(vec$high)) < abs(vec$load) #loads out of CI
    cvcheck=abs(vec$cv)>=cvth #cv more than threshold
    invalidloads <- c(which( mincheck| maxcheck | cvcheck))
    invalqvalue=which(sigvec!=1)
    invalids=c(invalidperms,invalidloads)#,invalqvalue)
    
    vec$load_shape=matrix('significant_loads',length(vec$low))
    vec$load_shape[which(mincheck+maxcheck>0)]='loads out of CI'
    vec$load_shape[which(sigvec==0)]='insignificant_loads'
    
    vec$CI_color=matrix('proper CI',length(vec$low))
    vec$CI_color[which(cvcheck)]= paste('cv >',cvth)
    vec$CI_color[invalidperms]='zero in CI'
    if(plotflag){
      image=ggplot(vec,aes(x = fea, y = org, col = CI_color, group= load_shape)) +
        geom_point(aes(shape = load_shape)) +
        geom_errorbar(aes(ymax = high, ymin = low),  width=0.25) +
        scale_shape_manual(values=c(4, 20,1))+
        # scale_color_manual(values = c("#ca7dcc","#18a558","#ff0000"))
        scale_color_manual(values = c("#ca7dcc","#18a558","#FF0000"))+
        theme(text = element_text(size = 20)) 
      ggsave(file=svgname, plot=image, width=8, height=5)
    }
    return(unique(invalids))
  }
  
ttestweights <- function(strucvar,funcvar,solofeat=FALSE){
    if (solofeat){
      pval=matrix(100,nrow = 2,ncol = 16)
      for(i in c(1:90)){
        pval[1,i]=t.test(strucvar$boot[,i],mu=strucvar$org[i,1])$p.value
        pval[2,i]=t.test(funcvar$boot[,i],mu=funcvar$org[i,1])$p.value
      }
      pval[is.nan(pval)] = 0
      qvalues=matrix(p.adjust(c(pval), method="BH"),nrow = 2,ncol=90)
      sig_qvalue=(qvalues<0.05)*1
    }else{
      pval=matrix(100,nrow = 6,ncol = 16)
      for(i in c(1:16)){
        pval[1,i]=t.test(strucvar$myelin$boot[,i],mu=strucvar$myelin$org[1,i])$p.value
        pval[2,i]=t.test(strucvar$thick$boot[,i],mu=strucvar$thick$org[1,i])$p.value
        pval[3,i]=t.test(strucvar$sulc$boot[,i],mu=strucvar$sulc$org[1,i])$p.value
        pval[4,i]=t.test(strucvar$curv$boot[,i],mu=strucvar$curv$org[1,i])$p.value
        pval[5,i]=t.test(strucvar$surf$boot[,i],mu=strucvar$surf$org[1,i])$p.value
        pval[6,i]=t.test(funcvar$boot[,i],mu=funcvar$org[1,i])$p.value
      }
      pval[is.nan(pval)] = 0
      qvalueS=matrix(p.adjust(c(pval[1:5,]), method="BH"),nrow = 5,ncol=16)
      qvalueF=p.adjust(c(pval[6,]), method="BH")
      qvalues=matrix(p.adjust(c(pval), method="BH"),nrow = 6,ncol=16)
      sig_qvalueS=(qvalueS<0.05)*1
      sig_qvalueF=(qvalueF<0.05)*1
      sig_qvalue=(qvalues<0.05)*1
    }
    return(sig_qvalue)
  }
  
mainvecs <-function(maindir,outbootype,plotorrs=FALSE,plotname=FALSE,match=NULL){
    dsmain <-read.csv(paste(maindir,"cors.csv",sep=""))[,2]
    usmain <-read.csv(paste(maindir,"u.csv",sep=""))[,-1]
    vsmain <-read.csv(paste(maindir,"v.csv",sep=""))[,-1]
    if(is.null(match)){match=1}
    if(match=="max"){match=which.max(dsmain)}
    usmain <- usmain[,match]
    vsmain <- vsmain[,match]
    dsmain <- dsmain[match]
    out <-list(v=vsmain,u=usmain,d=dsmain)
    return (out)
  }


###################CS-CF generation###############
########CCA for selected features
sig_corrs <- function(x,y,zeros=FALSE){
  # Calculate correlations between columns of X and columns of Y
  correlations <- c()
  for (i in 1:ncol(x)) {
    for (j in  1:ncol(y)) { #c(i,i+16,i+32,i+48,i+64)
      corinfo=cor.test(scale(x[, i]), scale(y[, j]))
      if (corinfo$p.value<0.05){correlations <- c(correlations,corinfo$estimate)}
      else{if(zeros){correlations <- c(correlations,0)}}
    }
  }
  return(correlations)
}

plotcorvalues<-function(x,y,CCC,plotname){
  # Calculate correlations between columns of X and columns of Y
  correlations <- c()
  correlations=sig_corrs(x,y)
  # Create a histogram with a vertical red line
  ggplot(data.frame(correlations), aes(x = correlations)) +
    geom_histogram(binwidth = 0.05, fill = "#ccc", color = "black") +
    geom_vline(aes(xintercept = CCC), color = "red", linetype = "dashed", size = 1) +
    annotate("text", x = CCC, y = 10, label = paste("CCC =", round(CCC,digits = 2)), angle = 90, vjust =1.5, color = "red", size = 4) +
    labs(title = " ", x = "Correlation Values", y = " ") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # White background
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(), 
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      plot.title = element_text(size = 20))
  ggsave(filename = paste(plotname,"hist_original_CCC.tiff"), width = 5, height = 5, device='tiff', dpi=400)
  
  correlation_matrix <- cor(x, y)
  
  roi_names=c('opr_l','opr_r','trian_l','trian_r','insula_l','insula_r','supr_l','supr_r','ang_l','ang_r','hsl_l','hsl_r','STG_l','STG_r','MTG_l','MTG_r')
  categories_S=c('myelin', 'thick', 'sulc', 'curv', 'surf')
  
  correlation_df <- melt(t(correlation_matrix))
  colnames(correlation_df) <- c("Feature_S", "Feature_F", "Correlation")
  correlation_df$Feature_F <- factor(correlation_df$Feature_F, levels = 1:16, labels = roi_names)
  correlation_df$Feature_S <- factor(correlation_df$Feature_S, levels = 1:80, labels = rep(roi_names,5))
  correlation_df$Category_S <- rep(categories_S, each = 16)
  
  ggplot(correlation_df, aes(x = Feature_S, y = Feature_F, fill = Correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(min(correlation_matrix), max(correlation_matrix)), name = "Correlation") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = "Correlation Heatmap", x = "Features in S", y = "Features in F") +
    facet_wrap(~ Category_S, scales = "free_x", nrow = 1)  # Divide into 5 categories
  ggsave(filename = paste(plotname,"corr_heatmap_original.tiff"), width = 10, height = 5, device='tiff', dpi=400)
  
  
}

cca_comps_main <- function(rois,x,z,savepath,partname=""){
  
  X <- scale(x[,rois[rois<=dim(x)[2]]], TRUE, TRUE)
  Z <- scale(z[,rois], TRUE, TRUE)
  k= min(dim(X)[2],dim(Z)[2])
  cca_result<-cca(X,Z)
  
  N = nrow(Z)
  featx=ncol(x)
  featz=ncol(z)
  # Set up parallel backend
  no_cores <- detectCores() - 2
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  permnums=1000
  
  #permutation test
  ccperms <- foreach(perm = 1:permnums, .combine = rbind, .packages = 'yacca') %dopar% {
    set.seed(perm+222)
    sampz <- sample(1:nrow(z))
    sampx <- sample(1:nrow(x))
    Xperm <- X[sampx, ]
    Zperm <- Z[sampz, ]
    out <- cca(Xperm, Zperm)
    return(out[["corr"]])
  }
  stopCluster(cl)
  
  maincor = cca_result$corr
  
  pvals=(sapply(seq_along(1:k), function(i) length(which(ccperms[,i] >= maincor[i])) / permnums))
  comp=which(pvals<0.05)[1] #select the first significant component
  if(length(comp)>0){
    comp=comp[1]
    print(paste("component",comp, "has significant results."))
    
    #bootstaping for 90%CI of CCC
    bootnums=5000
    uplim=0.95
    lwlim=0.05
    corrs_list<- matrix(nrow = bootnums, ncol = 1)
    us_list <- matrix(nrow = bootnums, ncol = ncol(X))
    vs_list <- matrix(nrow = bootnums, ncol = ncol(Z))
    
    for (perm in 1:bootnums) {
      set.seed(perm+222)
      ind <- sample(N, N, replace = TRUE)
      Zboot <- Z[ind, ]
      Xboot <- X[ind, ]
      out <- cca(Xboot, Zboot)
      
      Xboot <- scale(Xboot, TRUE, TRUE)
      Zboot <- scale(Zboot, TRUE, TRUE)
      
      corrs_list[perm] <- out$corr[comp]
      us_list[perm, ] <- out$xcoef[, comp]
      vs_list[perm, ] <- out$ycoef[, comp]
    }
    
    markerx=X%*%cca_result[["xcoef"]][, comp]
    markerz=Z%*%cca_result[["ycoef"]][, comp]
    #correlation between individual features and canonical features
    xcrossz_main=cor(x, markerz)
    xcrossz_main_sig=sig_corrs(x, markerz)
    zcrossx_main=cor(z,markerx)
    zcrossx_main_sig=sig_corrs(z,markerx)
    xloading_main=cor(x, markerx)
    xloading_main_sig=sig_corrs(x, markerx)
    zloading_main=cor(z,markerz)
    zloading_main_sig=sig_corrs(z,markerz)
    
    corsboot=abs(maincor[comp])-(abs(corrs_list)-abs(maincor[comp]))
    corrs.ci <-  quantile((corsboot),c(uplim,lwlim),na.rm =T)
    
    #plot histogram
    plotcorvalues(x=x,y=z,CCC=maincor[comp],plotname=paste0(savepath,"FALFF"))
    plotdata <- data.frame(x = markerx, y = markerz)
    
    model <- lm(y ~ x, data = plotdata)# Fit a linear model
    
    # Extract the slope and p-value using broom
    model_summary <- tidy(model)
    slope <- model_summary$estimate[2]
    p_value <- model_summary$p.value[2]
    Rsquared_adj <- summary(model)$adj.r.squared
    
    # Create scatter plot with regression line and slope
    ggplot(plotdata, aes(x = x, y = y)) +
      geom_point(color = "blue") +
      geom_smooth(method = "lm", color = "red", se = FALSE) +
      annotate("text", x = Inf, y = Inf, label = paste("Slope =", round(slope, 2)), hjust = 1.1, vjust = 2, color = "red", size = 5) +
      labs(title = paste("pvalue=",round(p_value, 4)," Adj_R",supsc(2),"=",round(Rsquared_adj, 4)), x = "CF feature", y = "CS feature") +
      theme_classic() +
      theme(axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            plot.title = element_text(size = 14))
    ggsave(filename = paste0(savepath,"scatter_cca_variates.png"), width = 5, height = 3, dpi = 600)
    
    bootinfo=list(corrs.ci=corrs.ci,corrs.mean=mean(corsboot))
    maininfo=list(corr=maincor[comp],xcoef=cca_result[["xcoef"]][, comp],ycoef=cca_result[["ycoef"]][, comp],
                  xcrossz_main=xcrossz_main,zcrossx_main=zcrossx_main,
                  xloading_main=xloading_main,zloading_main=zloading_main)
    siginfo=list(xcrossz_main_sig=xcrossz_main_sig,zcrossx_main_sig=zcrossx_main_sig,
                 xloading_main_sig=xloading_main_sig,zloading_main_sig=zloading_main_sig)
    
    return(list(boot=bootinfo,main=maininfo,sigmain=siginfo))
  }else{
    print("no significant CCA components found for data")
  }
}

save_cca_weights<-function(rois,funccoefs,struccoefs,savename){
  funcs=matrix(0,ncol = 16,nrow=1)
  funcs[1,rois[rois<17]]=funccoefs
  strucs=matrix(0,ncol=80,nrow=1)
  strucs[1,rois]=struccoefs
  
  write.csv(funcs,paste0(savename,"_func.csv"))
  write.csv(strucs,paste0(savename,"_struc.csv"))
}
