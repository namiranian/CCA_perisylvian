library(PMA)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
library(plyr)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
maindirecory="C:/Users/Namiranan/Downloads/project_root/"
source(paste0(maindirecory,"FC/funcs_FC.R"))

numCores <- detectCores()
registerDoParallel(numCores-1)  # use multicore, set to the number of our cores
perisylvians=c(11,12,13,14,29,30,63,64,65,66,79,80,81,82,85,86)
initind=c(1,79,157,235,313,391)
strucnames=c('myelin', 'thick', 'sulc', 'curv', 'surf','all')

#paths to save and load data
savedir= paste0(maindirecory,"FC/")
paths=paste0(maindirecory,"input_data/age_reg_out_metrics/morphologicalfeats.csv")
zread=read.csv(paths,header = TRUE)
pathf_dir=paste0(maindirecory,"input_data/age_reg_out_metrics/FC_")
savepath1=paste0(savedir,"feature_selection/")
if(!dir.exists(savepath1)){dir.create(savepath1)}
savepath2=paste0(savedir,"CS-CF/")
if(!dir.exists(savepath2)){dir.create(savepath2)}


#####################Permutation#############
temp=seq(0,1,by=0.1)
temp[1]=0.01
temp[11]=0.98
penals <- data.frame(expand.grid(penalx = temp,penalz = temp))
permsinfo=c()
for(seed in perisylvians){
  pathf=paste0(pathf_dir,seed,"_to_all.csv")
  xread <- read.csv(pathf,header = TRUE)

  data<-get_data(xread,zread,featid=c(1:5),selectedrois=perisylvians,seed=seed)
  perm.out <-permcca_parallel(penalsx=penals$penalx,penalsz=penals$penalz,
                              x=data[["x"]],z=data[["z"]],k=1,permnums=500)

inds=which(perm.out[["pvals"]]<0.05)
if(length(inds)!=0){
if(length(inds)>1){
  temp=perm.out[["cors"]][inds]
  maxcor=temp[which.max(temp)]
  indbest=inds[which(temp==maxcor)]
  if(length(indbest)>1){
    # temp=perm.out[["cors"]][indbest]
    temp.featnum <- penals$penalx[indbest]*89+penals$penalz[indbest]*390
    indbest=indbest[which.max(temp.featnum)]
  }
}else{
  indbest=which.min(perm.out[["pvals"]])
}
permsinfo=rbind(permsinfo,c("seeds"=seed,"pval"=perm.out[["pvals"]][indbest],
                            "penaltyx"=penals$penalx[indbest],
                            "penaltyz"=penals$penalz[indbest],
                            "cors"=perm.out[["cors"]][indbest]))

}

}
write.csv(permsinfo,paste0(savepath1,"permsinfo_seeds.csv"))


######## main########
penalties <- read.csv(paste0(savepath1,"permsinfo_seeds.csv"),header = TRUE)
sig_seeds=penalties$seeds#c(13,14,29,30,63,64,79,80,81,82)

savedir_temp=paste0(savepath1,'/main/')
if(!dir.exists(savedir_temp)){dir.create(savedir_temp)}

for(seed in sig_seeds){
  bests=penalties[which(penalties$seeds==seed),]
  bestx=bests$penaltyx
  bestz=bests$penaltyz
  pathf=paste0(pathf_dir,seed,"_to_all.csv")
  xread <- read.csv(pathf,header = TRUE)
  
  data<-get_data(xread,zread,featid=c(1:5),selectedrois=perisylvians,seed=seed)
  
  savepath=paste0(savedir_temp,'/FC_',seed)
  if(!dir.exists(savepath)){dir.create(savepath)}
  pmacca_main(bestpenaltyx=bestx,bestpenaltyz=bestz,data[["x"]],data[["z"]],
              savepath,partname=paste0("f",bestx,"_s",bestz,"main"))
}

#################checking significance of other canonical components (permutation) ###############
penalties <- read.csv(paste0(savepath1,"permsinfo_seeds.csv"),header = TRUE)
sig_seeds=penalties$seeds #c(13,14,29,30,63,64,79,80,81,82)
permnums=500
  for(seed in sig_seeds){
    bests=penalties[which(penalties$seeds==seed),]
    bestx=bests$penaltyx
    bestz=bests$penaltyz
    pathf=paste0(pathf_dir,seed,"_to_all.csv")
    xread <- read.csv(pathf,header = TRUE)
    
    data<-get_data(xread,zread,featid=c(1:5),selectedrois=perisylvians,seed=seed)
    
    perms_cor=permcca2_parallel2(penalsx=bestx,penalsz=bestz,data[["x"]],data[["z"]], k=0,permnums=permnums)
  
    savepath=paste0(savepath1,'/main/FC_',seed)
  maincor <- read.csv(paste0(savepath,"/mainf",bestx,"_s",bestz,"maincors.csv"),header = TRUE)[,-1]
  pvals <- sapply(seq_along(1:length(maincor)), function(i) length(which(perms_cor[,i] >= maincor[i]))) / permnums
  components=which(pvals<0.05)
  print(components)
  write.csv(components,paste0(savepath1,"sigcomponents_permutation_seed_",seed,".csv"))
}

#################apllying SCCA_bootstrapping ###############
penalties <- read.csv(paste0(savepath1,"permsinfo_seeds.csv"),header = TRUE)
sig_seeds=penalties$seeds
savepath=paste0(savepath1,'/boot/')
if(!dir.exists(savepath)){dir.create(savepath)}
numCores <- detectCores()
registerDoParallel(numCores-1) 
foreach(seed=sig_seeds)%dopar%{
  bests=penalties[which(penalties$seeds==seed),]
  bestx=bests$penaltyx
  bestz=bests$penaltyz
  pathf=paste0(pathf_dir,seed,"_to_all.csv")
  xread <- read.csv(pathf,header = TRUE)
  
  data<-get_data(xread,zread,featid=c(1:5),selectedrois=perisylvians,seed=seed)
  
  savepath=paste0(savepath1,'/boot/FC_',seed)
  if(!dir.exists(savepath)){dir.create(savepath)}
  components= read.csv(paste0(savepath1,"sigcomponents_permutation_seed_",seed,".csv"))[,-1]
    Bootcca(bestpenaltyx=bestx,bestpenaltyz=bestz,x=data[["x"]],z=data[["z"]],k=max(components),
          savepath,partname='',iternums=5000)
}


###############visualization (data preparation)##############
savepath=paste0(savepath1,'/vsig/')
if(!dir.exists(savepath)){dir.create(savepath)}
for(seed in sig_seeds){
  bests=penalties[which(penalties$seeds==seed),]
  bestx=bests$penaltyx
  bestz=bests$penaltyz
  
  savepath=paste0(savepath1,'/vsig/FC_',seed)
  if(!dir.exists(savepath)){dir.create(savepath)}
  
  comps <- read.csv(paste0(savepath1,"sigcomponents_permutation_seed_",seed,".csv"))[,-1]
  rois=c()
  sum.abs.feature=matrix(0,nrow=1,ncol=6)
  
  for (comp in comps){
    orgvecs<-mainvecs(maindir=paste0(savepath1,"/main/FC_",seed,"/mainf",bestx,"_s",bestz,"main"),
                      match=comp)
  u_boost <- read.csv(paste0(savepath1,"/boot/FC_",seed,"/perm",comp,"u_comp",comp,".csv"))[,-1]
  v_boost <- read.csv(paste0(savepath1,"/boot/FC_",seed,"/perm",comp,"v_comp",comp,".csv"))[,-1]
  cor_boost <- read.csv(paste0(savepath1,"/boot/FC_",seed,"/perm",comp,"cors_comp",comp,".csv"))[,-1]
  
  weights=plotcirc3 (orgvecs,u_boost,v_boost,cor_boost,savepath=paste(savepath,'/FC_',seed,"_1000_comp",comp,"_",sep=''),seed=seed)
  
  # write.csv(weights$mainloads,paste0(savepath,"_1000_comp",comp,"_mainloads.csv"))
  # write.csv(weights$sigloads,paste0(savepath,"_1000_comp",comp,"_sigloads.csv"))
  # write.csv(weights$mainpower,paste0(savepath,"_1000_comp",comp,"_mainpower.csv"))
  # write.csv(weights$sigpower,paste0(savepath,"_1000_comp",comp,"_sigpower.csv"))
  # write.csv(weights$corr_bc,paste0(savepath,"_1000_comp",comp,"_corr.csv"))
  
  sigload=matrix(unlist(weights$sigloads),nrow=16,ncol=6)
  loads=colSums(abs(sigload))
  sum.abs.feature=sum.abs.feature+loads
  for (iter in comps){
    rois=append(rois,c(which(abs(sigload[,1])>0)))
  }
  for (iter in comps){
    strucloads=as.vector((sigload[,2:6]))
    rois=append(rois,which(abs(strucloads)>0))
  }
  }
  rois=sort(unique(c(unlist(rois))))
  if(length(rois)>0){
    write.csv(rois,paste0(savepath,"/rois_using_all_components.csv"))
    print(rois)
    write.csv(sum.abs.feature,paste0(savepath,"/weights_using_all_components.csv"))
  }
  
}

#####CCA WITH SELECTED ROIS#####
library(yacca)
library(broom)
library(reshape2)
penalties <- read.csv(paste0(savepath1,"permsinfo_seeds.csv"),header = TRUE)
sig_seeds=penalties$seeds
for(seed in sig_seeds){
  print(paste("seed=",seed))
  pathf=paste0(pathf_dir,seed,"_to_all.csv")
  xread <- read.csv(pathf,header = TRUE)
  
  data<-get_data(xread,zread,featid=c(1:5),selectedrois=perisylvians,seed=seed)
  savepath=paste0(savepath2,'/FC',seed)
  if(!dir.exists(savepath)){dir.create(savepath)}

  roipath=paste0(savepath1,'/vsig/FC_',seed)
  rois <- read.csv(paste0(roipath,'/rois_using_all_components.csv'),header = TRUE)[,-1]
  # rois <- read.csv("G:/SF9_2024/CCA/release2_com/upload_pre/FALFF16/FC16/FC80/boot/visig/rois_using_all_components.csv",header = TRUE)[,-1]

  cca_boot=cca_comps_main(rois,seed,x=data[["x"]],z=data[["z"]],savepath)
  if (is.list(cca_boot)){
      print(paste("corr.ci=",cca_boot$boot$corrs.ci))
      if(cca_boot$boot$corrs.ci[2]>0.3){
        write.csv(cca_boot$boot$corrs.ci,paste0(savepath,"/5000_",seed,"_bootcorr_ci.csv"))
        
        write.csv(cca_boot$main[["xcoef"]],paste0(savepath,"/func_weights",seed,".csv"))
        write.csv(cca_boot$main[["ycoef"]],paste0(savepath,"/struc_weights",seed,".csv"))
        
        write.csv(cca_boot$main[["xcrossz_main"]],paste0(savepath,"/func_crossloading_CS",seed,".csv"))
        write.csv(cca_boot$main[["zcrossx_main"]],paste0(savepath,"/struc_crossloading_CF",seed,".csv"))

        write.csv(cca_boot$main[["xloading_main"]],paste0(savepath,"/func_loading",seed,".csv"))
        write.csv(cca_boot$main[["zloading_main"]],paste0(savepath,"/struc_loading",seed,".csv"))
        
        write.csv(cca_boot$main[["corr"]],paste0(savepath,"/maincorr_",seed,".csv"))
        
      }
  }
}
  
############plot CCA results

selected_seeds=c(30,79,80)
for(seed in selected_seeds){
  print(paste("seed=",seed))
  roipath=paste0(savepath1,'/vsig/FC_',seed)
  rois <- read.csv(paste0(roipath,'/rois_using_all_components.csv'),header = TRUE)[,-1]
  
  savepath=paste0(savepath2,'/FC',seed)
  #plot weights
  temp=read.csv(paste0(savepath,"/func_weights",seed,".csv"))[,-1]
  region_id=which(perisylvians==seed)
  funcloads <- c(temp[1:region_id-1],0,temp[region_id:length(temp)])
  strucloads=read.csv(paste0(savepath,"/struc_weights",seed,".csv"))[,-1]
  cirplotfunc1(svgname=paste0(savepath,"/weights_main.svg"),strucloads,funcloads,region_id=region_id,plotmode="fullrange")
  
  #plot crossloadings
  savepath=paste0(savepath2,'/FC',seed)
  temp=read.csv(paste0(savepath,"/func_crossloading_CS",seed,".csv"))[,-1]
  region_id=which(perisylvians==seed)
  funcloads <- c(temp[1:region_id-1],0,temp[region_id:length(temp)])
  strucloads=read.csv(paste0(savepath,"/struc_crossloading_CF",seed,".csv"))[,-1]
  cirplotfunc1(svgname=paste0(savepath,"/crossloading_main.svg"),strucloads,funcloads,region_id=region_id,plotmode="corrs")
  # Compute the sum of absolute values for each feature
  crossloading_feats <- sapply(1:5, function(i) {sum(abs(strucloads[((i - 1) * 16 + 1):(i * 16)]))})
  write.csv(crossloading_feats,paste0(savepath,"/crossloading_strucfeats_",seed,".csv"))
  
  #plot loadings
  savepath=paste0(savepath2,'/FC',seed)
  temp=read.csv(paste0(savepath,"/func_loading",seed,".csv"))[,-1]
  region_id=which(perisylvians==seed)
  funcloads <- c(temp[1:region_id-1],0,temp[region_id:length(temp)])
  strucloads=read.csv(paste0(savepath,"/struc_loading",seed,".csv"))[,-1]
  cirplotfunc1(svgname=paste0(savepath,"/loading_main.svg"),strucloads,funcloads,region_id=region_id,plotmode="fullrange")
  # Compute the sum of absolute values for each feature
  loading_feats <- sapply(1:5, function(i) {sum(abs(strucloads[((i - 1) * 16 + 1):(i * 16)]))})
  write.csv(loading_feats,paste0(savepath,"/loading_strucfeats_",seed,".csv"))
  
}
