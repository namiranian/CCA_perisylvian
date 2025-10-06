get_data <-function(xread,zread,featid,selectedrois,seed){
  x <- matrix(unlist(xread),nrow=dim(xread)[1],ncol=dim(xread)[2])
  z <- matrix(unlist(zread),nrow=dim(zread)[1],ncol=dim(zread)[2])
  
  initind=c(1,79,157,235,313,391)
  
  roidel=c(37, 38, 41, 42, 71, 72, 73, 74, 75, 76, 77, 78)
  roilist <-c(1:90)
  roilist <- roilist[-roidel]   
  # roilist=selectedrois
  selectedrois_instruc=c()
  for (j in selectedrois){selectedrois_instruc=c(selectedrois_instruc,which(roilist==j))}
  newz=c()
  for (ii in featid){
    tempz=z[,initind[ii]:(initind[ii+1]-1)]
    newz=cbind(newz,tempz[,selectedrois_instruc])
  }
  
  selectedrois_func <- selectedrois[selectedrois != seed]
  selectedrois_func[which(selectedrois_func>seed)]=selectedrois_func[which(selectedrois_func>seed)]-1
  newx=x[,selectedrois_func]
  
  return(list(x=(newx),z=(newz)))
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
  write.csv(zcrossx_list, paste(savedir,"zcross.csv",sep=""))
  write.csv(zloading_list, paste(savedir,"zloading.csv",sep=""))
  write.csv(xloading_list, paste(savedir,"xloading.csv",sep=""))
  
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
  }
  stopCluster(cl)
  
  return(ccperms)
}

#######bootstraping stable canonical weights
Bootcca <- function(bestpenaltyx, bestpenaltyz, x, z, k,savepath, partname, iternums = 500, otherloadings = FALSE) {
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

bootstats2 <- function(org,bootdata,bias, uplim,lwlim){
  # boot.stats <- sapply(seq_along(1:length(org)), function(i) org[i] - quantile(bootdata[i,] - org[i],c(uplim,lwlim),na.rm =T))
  
  boot.stats <- sapply(seq_along(1:length(org)), function(i) quantile(bootdata[,i],c(uplim,lwlim),na.rm =T))
  boot.stats <- as.data.frame(t(boot.stats))
  colnames(boot.stats)<-c("low","high")
  boot.stats$ci <- abs(boot.stats$high - boot.stats$low)
  boot.stats$load <- sapply(seq_along(1:length(org)), function(i) mean(bootdata[,i]))
  boot.stats$fea <- 1:length(org)
  boot.stats$org <- matrix(unlist(org),nrow=length(org))
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

#https://github.com/QianLi9423/Project_SCZ_controllability/blob/main/Code_sCCA/ccafunction.R
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
    seed=info
    region_id=which(perisylvians==seed)
    sampnum <- dim(boot)[2]
    boot1 <- matrix(unlist(boot), nrow=dim(boot)[1], ncol=dim(boot)[2])
    org1 <- matrix(unlist(org), ncol=dim(boot)[2])
    
    if(region_id==1){
      boot1 <- cbind(matrix(0,nrow = dim(boot1)[1],ncol = 1),boot1)
      org1 <- rbind(matrix(0,nrow = 1,ncol = 1),matrix(org[region_id:sampnum],nrow = sampnum-region_id+1,ncol = 1))
    }else{
      boot1 <- cbind(boot1[,1:region_id-1],matrix(0,nrow = dim(boot1)[1],ncol = 1),boot1[,region_id:sampnum])
      org1 <- rbind(matrix(org[1:region_id-1],nrow = region_id-1,ncol = 1),matrix(0,nrow = 1,ncol = 1),
                    matrix(org[region_id:sampnum],nrow = sampnum-region_id+1,ncol = 1))
    }
    
    out=list(org=org1,boot=boot1)
  }
  return(out)
}

plotcirc3 <- function(orgvecs,u_boost,v_boost,cor_boost,savepath,seed){
  
  nameinfo="CI95_visall_"
  svgname=paste(savepath,nameinfo,sep="")
  
  struc.feat.info=c('myelin', 'thick', 'sulc', 'curv', 'surf')
  U_vec=bootdata4allregions(org=orgvecs$u,boot=u_boost,struc_data=FALSE,info=seed)
  V_vec=bootdata4allregions(org=orgvecs$v,boot=v_boost,struc_data=TRUE,info=struc.feat.info)
  
  funcvar=bootplot2(org=abs(U_vec$org),boot=abs(U_vec$boot),plotname=paste(savepath,"functions"))
  
  corr_bc=(sapply(seq_along(1:length(orgvecs$d)), function(i) orgvecs$d[i]-(mean(cor_boost) - orgvecs$d[i])))
  
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
  
  return(out=list(mainloads=mainloads,sigloads=sigloads,mainpower=mainpower,sigpower=sigpower,corr_bc=corr_bc))
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
      pval[6,i]=t.test(funcvar$boot[,i],mu=funcvar$org[i,1])$p.value
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

plotcorvalues<-function(x,y,CCC,seed,plotname){
  library(reshape2)
  
  # Calculate correlations between columns of X and columns of Y
  correlations=sig_corrs(x,y)
  
  # Create a histogram with a vertical red line
  ggplot(data.frame(correlations), aes(x = correlations)) +
    geom_histogram(binwidth = 0.05, fill = "#ccc", color = "black") +
    geom_vline(aes(xintercept = CCC), color = "red", linetype = "dashed", size = 1) +
    annotate("text", x = CCC, y = 10, label = paste("CCC =", round(CCC,digits = 2)), angle = 90, vjust =1.5, color = "red", size = 6) +
    labs(title = " ", x = "Correlation Values", y = " ") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # White background
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(), 
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 16),
      plot.title = element_text(size = 20))
  ggsave(filename = paste(plotname,"hist_original_CCC.tiff"), width = 5, height = 5, device='tiff', dpi=400)
  
  temp <- cor(x, y)
  correlation_matrix=matrix(0,nrow=16,ncol=80)
  correlation_matrix[1:(seed-1),]=temp[1:(seed-1),]
  correlation_matrix[(seed+1):16,]=temp[seed:15,]
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

cca_comps_main <- function(rois,seed,x,z,savepath,partname=""){
  
  perisylvians=c(11,12,13,14,29,30,63,64,65,66,79,80,81,82,85,86)
  seed=which(perisylvians==seed)
  Z <- scale(z[,rois], TRUE, TRUE)
  roif=rois[rois<=16]
  roiadd=rois[rois>16]%%16
  roif=unique(sort(c(roif,roiadd)))
  selectedrois_func <- roif[roif != seed]
  selectedrois_func[which(selectedrois_func>seed)]=selectedrois_func[which(selectedrois_func>seed)]-1
  X <- scale(x[,selectedrois_func], TRUE, TRUE)
  
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
  comp=which(pvals<0.05)
  if(length(comp)>0){
    comp=comp[1]
    print(paste("component",comp, "has significant results."))
    markerx=X%*%cca_result[["xcoef"]][, comp]
    markerz=Z%*%cca_result[["ycoef"]][, comp]
    xcrossz_main=cor(x, markerz)
    xcrossz_main_sig=sig_corrs(x, markerz)
    zcrossx_main=cor(z,markerx)
    zcrossx_main_sig=sig_corrs(z,markerx)
    xloading_main=cor(x, markerx)
    xloading_main_sig=sig_corrs(x, markerx)
    zloading_main=cor(z,markerz)
    zloading_main_sig=sig_corrs(z,markerz)
    
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
    }
    
    corsboot=abs(maincor[comp])-(abs(corrs_list)-abs(maincor[comp]))
    print(mean(corsboot))
    corrs.ci <-  quantile((corsboot),c(uplim,lwlim),na.rm =T)
    corrs.cv=mean(corsboot)/sd(corsboot)*100
    plotcorvalues(x=x,y=z,CCC=maincor[comp],seed=seed,plotname=paste0(savepath,"/FC_",seed))
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
    ggsave(filename = paste0(savepath,"/scatter_cca_variates.png"), width = 5, height = 3, dpi = 600)
    
    
    xweights=matrix(0,nrow=dim(x)[2])
    zweights=matrix(0,nrow=dim(z)[2])
    xweights[selectedrois_func]=cca_result[["xcoef"]][, comp]
    zweights[rois]=cca_result[["ycoef"]][, comp]
    
    bootinfo=list(corrs.ci=corrs.ci,corrs.mean=mean(corsboot))
    maininfo=list(corr=maincor[comp],xcoef=xweights,ycoef=zweights,
                  xcrossz_main=xcrossz_main,zcrossx_main=zcrossx_main,
                  xloading_main=xloading_main,zloading_main=zloading_main)
    siginfo=list(xcrossz_main_sig=xcrossz_main_sig,zcrossx_main_sig=zcrossx_main_sig,
                 xloading_main_sig=xloading_main_sig,zloading_main_sig=zloading_main_sig)
    return(list(boot=bootinfo,main=maininfo,sigmain=siginfo))
  }else{
    print("no significant CCA components found for data")
  }
}

############plot CCA results
cirplotfunc <- function(svgname,strucloads,funcloads,region_id,plotmode){
  print(svgname)
  if(plotmode=="allstrucfeat"){
    sectorsvec=c(c(1:8)*2,c(7:0)*2+1)
    region_id=which(sectorsvec==region_id)
    roiabbrev=c('operG','operG','trianG','trianG','insula','insula',
                'suprG','suprG','AG','AG','HG','HG','STG','STG','MTG','MTG')[sectorsvec]
    myelin=(strucloads[1:16])[sectorsvec]
    thick=(strucloads[17:32])[sectorsvec]
    curv=(strucloads[33:48])[sectorsvec]
    sulc=(strucloads[49:64])[sectorsvec]
    surf=(strucloads[65:80])[sectorsvec]
    minval=min(strucloads)
    maxval=max(strucloads)
    funcloads=funcloads[sectorsvec]
    if (maxval-minval!=0){
      sectorsvec=c(1:16)
      svg(svgname,width=10,height = 10)
      circos.par("track.height" = 0.5, start.degree = 90, cell.padding = c(0, 0, 0, 0),gap.degree=0)
      circos.initialize(sectorsvec, xlim = c(0, 1))
      
      sector_design = function(value,minval,maxval) {
        circos.rect(0.015,minval+0.01,0.25,maxval-0.01, col = "#6666FF",border = NA)#99FF99
        circos.rect(0.2,minval+0.01,0.4,maxval-0.01, col = "#66FFFF",border= NA)#AAFFFF
        circos.rect(0.4,minval+0.01,0.6,maxval-0.01, col = "#66FF33",border = NA)#9999FF
        circos.rect(0.6,minval+0.01,0.8,maxval-0.01, col = "#FFFF66",border= NA)#FFAAAA
        circos.rect(0.8,minval+0.01,0.985,maxval-0.01, col = "#FF6666",border= NA)#CCCCCC
        for (i in c(1:5)){
          circos.barplot(value[i],0.1+(i-1)*0.2,col="black",border = "white",bar_width=0.17)
        }
      }
      
      circos.track(ylim = c(minval, maxval), sectors = sectorsvec, bg.border = "black",bg.col ="black",
                   panel.fun = function(x, y) {
                     l = get.cell.meta.data("sector.numeric.index")
                     value= c(myelin[l], thick[l] ,sulc[l] ,curv[l] ,surf[l])
                     # value=c(0.2,-0.1,0.4,-0.05,0)
                     sector_design(value,minval=minval,maxval=maxval)
                     circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + mm_y(4),
                                 roiabbrev[l],niceFacing = TRUE,cex = 1.5)
                   })
      
      for (i in c(1,4,8,12)){
        circos.yaxis(side = "right", sector.index = i,labels.cex = 0.6,tick.length = 0.3)
      }
      
      
      ind=which(funcloads!=0)
      if(length(ind)>0){
        minval=min(funcloads[ind])
        maxval=max(funcloads[ind])
        step=(maxval-minval)/4
        rangeload=seq(minval,maxval,by=step)
        
        if(length(rangeload)>1){
          minval=min(funcloads[ind])
          maxval=max(funcloads[ind])
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
            circos.link(region_id, 0.5, i, 0.5,
                        col = add_transparency(col_fun(funcloads[i])))
          }
          lgd_links = Legend(at = round(c(minval,minval+(maxval-minval)/2,maxval),digits = 3), col_fun = col_fun,
                             title_position = "topleft", title = "FCs")
          
          textvalue=c("myelin","thickness","sulcal_depth","curvature","surface_area")
          lgd_lines = Legend(at = textvalue, type = "lines",
                             legend_gp = gpar(col = c("#6666FF","#66FFFF","#66FF33","#FFFF66","#FF6666"), lwd = 4), title_position = "topleft",
                             title = "structural\nfeatures")
          lgd_list_vertical = packLegend(lgd_links,lgd_lines)
          
          draw(lgd_list_vertical, x = unit(7, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
        }
      }
      circos.clear()
      dev.off()
      
    }else{
      print("No plot: range of loads is zero.")
    }
  }
}

cirplotfunc1 <- function(svgname,strucloads,funcloads,region_id,plotmode="halfrange"){
  print(svgname)
    sectorsvec=c(c(1:8)*2,c(7:0)*2+1)
    region_id=which(sectorsvec==region_id)
    roiabbrev=c('operG','operG','trianG','trianG','insula','insula',
                'suprG','suprG','AG','AG','HG','HG','STG','STG','MTG','MTG')[sectorsvec]
    myelin=(strucloads[1:16])[sectorsvec]
    thick=(strucloads[17:32])[sectorsvec]
    curv=(strucloads[33:48])[sectorsvec]
    sulc=(strucloads[49:64])[sectorsvec]
    surf=(strucloads[65:80])[sectorsvec]
    minval=min(strucloads)
    maxval=max(strucloads)
    funcloads=funcloads[sectorsvec]
    if (maxval-minval!=0){
      sectorsvec=c(1:16)
      svg(svgname,width=10,height = 10)
      circos.par("track.height" = 0.5, start.degree = 90, cell.padding = c(0, 0, 0, 0),gap.degree=0)
      circos.initialize(sectorsvec, xlim = c(0, 1))
      
      sector_design = function(value,minval,maxval) {
        circos.rect(0.015,minval+0.01,0.25,maxval-0.01, col = "#6666FF",border = NA)#99FF99
        circos.rect(0.2,minval+0.01,0.4,maxval-0.01, col = "#66FFFF",border= NA)#AAFFFF
        circos.rect(0.4,minval+0.01,0.6,maxval-0.01, col = "#66FF33",border = NA)#9999FF
        circos.rect(0.6,minval+0.01,0.8,maxval-0.01, col = "#FFFF66",border= NA)#FFAAAA
        circos.rect(0.8,minval+0.01,0.985,maxval-0.01, col = "#FF6666",border= NA)#CCCCCC
        for (i in c(1:5)){
          circos.barplot(value[i],0.1+(i-1)*0.2,col="black",border = "white",bar_width=0.17)
        }
      }
      
      circos.track(ylim = c(minval, maxval), sectors = sectorsvec, bg.border = "black",bg.col ="black",
                   panel.fun = function(x, y) {
                     l = get.cell.meta.data("sector.numeric.index")
                     value= c(myelin[l], thick[l] ,sulc[l] ,curv[l] ,surf[l])
                     # value=c(0.2,-0.1,0.4,-0.05,0)
                     sector_design(value,minval=minval,maxval=maxval)
                     circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + mm_y(4),
                                 roiabbrev[l],niceFacing = TRUE,cex = 1.5)
                   })
      
      for (i in c(1,4,8,12)){
        circos.yaxis(side = "right", sector.index = i,labels.cex = 0.6,tick.length = 0.3)
      }
      
      
      ind=which(funcloads!=0)
      if(length(ind)>0){
        if (plotmode=="fullrange"){
          rangeload=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1)
          col_fun = colorRamp2(rangeload, c("#120175","#1F02CA","#2E0CFC","#ADA0FE","white","#FFB3B3","#FF7575","#D60000","#760000"))
          showrange=c(-1,0,1)
        }else{
          rangeload=c(-0.6,-0.45,-0.3,-0.15,0,0.15,0.3,0.45,0.6)
          col_fun = colorRamp2(rangeload, c("#120175","#1F02CA","#2E0CFC","#ADA0FE","white","#FFB3B3","#FF7575","#D60000","#760000"))
          showrange=c(-0.6,0,0.6)
          
        }
        # rangeload=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1)
        # col_fun = colorRamp2(rangeload, c("#120175","#1F02CA","#2E0CFC","#ADA0FE","white","#FFB3B3","#FF7575","#D60000","#760000"))
        
        
          for (i in ind){
            circos.link(region_id, 0.5, i, 0.5,
                        col = add_transparency(col_fun(funcloads[i])))
          }
          lgd_links = Legend(at = showrange, col_fun = col_fun, title_position = "topleft", title = "FCs")
          
          textvalue=c("myelin","thickness","sulcal_depth","curvature","surface_area")
          lgd_lines = Legend(at = textvalue, type = "lines",
                             legend_gp = gpar(col = c("#6666FF","#66FFFF","#66FF33","#FFFF66","#FF6666"), lwd = 4), title_position = "topleft",
                             title = "structural\nfeatures")
          lgd_list_vertical = packLegend(lgd_links,lgd_lines)
          
          draw(lgd_list_vertical, x = unit(7, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
        
      }
      circos.clear()
      dev.off()
      
    }else{
      print("No plot: range of loads is zero.")
    }
  
}
