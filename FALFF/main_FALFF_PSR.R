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
source(paste0(maindirecory,"FALFF/funcs_FALFF_PSR.R"))

 # use multicore, set to the number of our cores
numCores <- detectCores()
registerDoParallel(numCores-1) 
# #PSRs,structural features name and index
perisylvians=c(11,12,13,14,29,30,63,64,65,66,79,80,81,82,85,86)
initind=c(1,79,157,235,313,391)
strucnames=c('myelin', 'thick', 'sulc', 'curv', 'surf','all')

#paths to save and load data
savedir= paste0(maindirecory,"FALFF/")
paths=paste0(maindirecory,"input_data/age_reg_out_metrics/morphologicalfeats.csv")
pathf=paste0(maindirecory,"input_data/age_reg_out_metrics/Falff.csv")
xread=read.csv(pathf,header = TRUE)
zread=read.csv(paths,header = TRUE)
savepath1=paste0(savedir,"feature_selection/")
if(!dir.exists(savepath1)){dir.create(savepath1)}
savepath2=paste0(savedir,"CS-CF/")
if(!dir.exists(savepath2)){dir.create(savepath2)}
#selectedrois=perisylvians

#####################Permutation#############
temp=seq(0,1,by=0.1)
temp[1]=0.01
temp[11]=0.98
penals <- data.frame(expand.grid(penalx = temp,penalz = temp))
permsinfo=c()
data=get_data(xread,zread,featid=Null,perisylvians)
#permutation_test
perm.out <-permcca_parallel(penalsx=penals$penalx,penalsz=penals$penalz,
                            x=data[["x"]],z=data[["z"]],k=1,permnums=1000)
inds=which(perm.out[["pvals"]]<0.05)#significant sparsity
if(length(inds)>1){
  temp=perm.out[["cors"]][inds]
  maxcor=perm.out[["cors"]][which.max(temp)]
  indbest=inds[which(temp==maxcor)]
  if(length(indbest)>1){
    # temp=perm.out[["cors"]][indbest]
    temp.featnum <- penals$penalx[indbest]*89+penals$penalz[indbest]*390
    indbest=indbest[which.max(temp.featnum)]
  }
}else{
  indbest=which.min(perm.out[["pvals"]])
  print("There is no significant sparsity level")
}
i=6
permsinfo=rbind(permsinfo,c("struc_feat"=strucnames[i],
                            "zstat"=perm.out[["zstats"]][indbest],
                            "pval"=perm.out[["pvals"]][indbest],
                            "penaltyx"=penals$penalx[indbest],
                            "penaltyz"=penals$penalz[indbest],
                            "cors"=perm.out[["cors"]][indbest]))

write.csv(permsinfo,paste0(savepath1,"permsinfo_com_1.csv"))

######## main########
penalties <- read.csv(paste0(savepath1,"permsinfo_com_1.csv"),header = TRUE)[,-1]
xread=read.csv(pathf,header = TRUE)
zread=read.csv(paths,header = TRUE)

for (i in c(6)){
  print(strucnames[i])
  data=get_data(xread,zread,featid=Null,perisylvians)
  
  savepath=paste(savepath1,strucnames[i],"_main/",sep="")
  if(!dir.exists(savepath)){dir.create(savepath)}
  penalind=which(penalties$struc_feat==strucnames[i])
  bestx=penalties$penaltyx[penalind]
  bestz=penalties$penaltyz[penalind]
  pmacca_main(bestpenaltyx=bestx,bestpenaltyz=bestz,data[["x"]],data[["z"]],
              savepath,partname=paste("f",bestx,"_s",bestz,"main",sep=""))
}


#################checking significance of other canonical components (permutation) ###############
penalties <- read.csv(paste0(savepath1,"permsinfo_com_1.csv"),header = TRUE)[,-1]

xread=read.csv(pathf,header = TRUE)
zread=read.csv(paths,header = TRUE)
for (i in c(6)){
  print(strucnames[i])
  
  data=get_data(xread,zread,featid=Null,perisylvians)
  permnums=1000
  # if(penalties[["zstat"]][penalind]>1.96){
  # savepath=paste(savepath1,strucnames[i],"/",sep="")
  # if(!dir.exists(savepath)){dir.create(savepath)}
  penalind=which(penalties$struc_feat==strucnames[i])
  bestx=penalties$penaltyx[penalind]
  bestz=penalties$penaltyz[penalind]
  perms_cor=permcca2_parallel2(penalsx=bestx,penalsz=bestz,data[["x"]],data[["z"]], k=0,permnums=permnums)
  # }
}
maincor <- read.csv(paste0(savepath,"mainf",bestx,"_s",bestz,"maincors.csv"),header = TRUE)[,-1]
pvals <- sapply(seq_along(1:length(maincor)), function(i) length(which(perms_cor[,i] >= maincor[i]))) / permnums
components=which(pvals<0.05)
write.csv(components,paste0(savepath1,"sigcomponents_permutation.csv"))

#################apllying SCCA_bootstrapping ###############
penalties <- read.csv(paste0(savepath1,"permsinfo_com_1.csv"),header = TRUE)[,-1]
components<- read.csv(paste0(savepath1,"sigcomponents_permutation.csv"))[,-1]
data=get_data(xread,zread,featid=Null,perisylvians)

for (i in c(6)){
  print(strucnames[i])
  penalind=which(penalties$struc_feat==strucnames[i])
  bestx=penalties$penaltyx[penalind]
  bestz=penalties$penaltyz[penalind]
  if(penalties[["zstat"]][penalind]>1.96){
    savepath=paste(savepath1,strucnames[i],"_boot",sep="")
    if(!dir.exists(savepath)){dir.create(savepath)}
    Bootscca(bestpenaltyx=bestx,bestpenaltyz=bestz,x=data[["x"]],z=data[["z"]],k=16,
            savepath,partname='',iternums=5000)
  }
}
###############visualization (feature and roi selection)##############
penalties <- read.csv(paste0(savepath1,"permsinfo_com_1.csv"),header = TRUE)[,-1]
comps <- read.csv(paste0(savepath1,"sigcomponents_permutation.csv"),header = TRUE)[,-1]
for (i in c(6)){
  print(strucnames[i])
  
  penalind=which(penalties$struc_feat==strucnames[i])
  bestx=penalties$penaltyx[penalind]
  bestz=penalties$penaltyz[penalind]
  
  # comps <- read.csv(paste0(savepath1,strucnames[i],"/boot/corrfdr_comps.csv"),header = TRUE)[,-1]
  # comps=1#which(comps<0.05)
  savepath=paste(savepath1,strucnames[i],"_visig/",sep="")
  if(!dir.exists(savepath)){dir.create(savepath)}
  mainfile=paste(savepath1,strucnames[i],"_main",sep="")
  bootpath=paste(savepath1,strucnames[i],"_boot",sep="")
  rois=list()
  sum.abs.feature=matrix(0,nrow=1,ncol=6)
  for (comp in comps){
    orgvecs<-mainvecs(maindir=paste(mainfile,"/mainf",bestx,"_s",bestz,"main",sep = ""),
                      match=comp)
    u_boost <- read.csv(paste(bootpath,"/perm",comp,"u_comp",comp,".csv",sep=""))[,-1]
    v_boost <- read.csv(paste(bootpath,"/perm",comp,"v_comp",comp,".csv",sep=""))[,-1]
    if (i<6){
      weights=plotcirc4 (orgvecs,u_boost,v_boost,savepath=paste(savepath,"_1000_comp",comp,"_",sep=''))
    }else{
      weights=plotcirc3 (orgvecs,u_boost,v_boost,savepath=paste(savepath,"_1000_comp",comp,"_",sep=''))
    }
    # write.csv(weights$mainloads,paste0(savepath,"_1000_comp",comp,"_mainloads_",strucnames[i],".csv"))
    # write.csv(weights$sigloads,paste0(savepath,"_1000_comp",comp,"_sigloads_",strucnames[i],".csv"))
    # write.csv(weights$mainpower,paste0(savepath,"_1000_comp",comp,"_mainpower_",strucnames[i],".csv"))
    # write.csv(weights$sigpower,paste0(savepath,"_1000_comp",comp,"_sigpower_",strucnames[i],".csv"))
    
    sigload=matrix(unlist(weights$sigloads),nrow=16,ncol=6)
    loads=colSums(abs(sigload))
    sum.abs.feature=sum.abs.feature+loads
    rois=append(rois,c(which(abs(sigload[,1])>0)))
    strucloads=as.vector((sigload[,2:6]))
    rois=append(rois,which(abs(strucloads)>0))
  }
  rois=sort(unique(c(unlist(rois))))
  write.csv(rois,paste0(savepath,"selected_features.csv"))
  
}


#############CCA for selected rois from components####
library(yacca)
library(broom)
library(reshape2)
library(reporter)

i=6
savepath=paste(savepath1,strucnames[i],"_visig/",sep="")
rois <- read.csv(paste0(savepath,"selected_features.csv"),header = TRUE)[,-1]
data=get_data(xread,zread,featid=Null,perisylvians)

CSF_info=cca_comps_main(rois,x=data[["x"]],z=data[["z"]],savepath2)

if (is.list(CSF_info)){
  print(paste("corr.ci=",CSF_info$boot$corrs.ci))
  if(CSF_info$boot$corrs.ci[2]>0.3){
    
    save_cca_weights(rois,funccoefs=CSF_info$main[["xcoef"]],struccoefs=CSF_info$main[["ycoef"]],
                  savename=paste0(savepath2,"cca_weights"))
    
    write.csv(CSF_info$main[["xcrossz_main"]],paste0(savepath2,"func_crossloads.csv"))
    write.csv(CSF_info$main[["xloading_main"]],paste0(savepath2,"func_func_loads.csv"))
    
    write.csv(CSF_info$main[["zcrossx_main"]],paste0(savepath2,"struc_crossloads.csv"))
    write.csv(CSF_info$main[["zloading_main"]],paste0(savepath2,"struc_struc_loads.csv"))
    
    write.csv(CSF_info$sigmain[["xcrossz_main_sig"]],paste0(savepath2,"sig_func_crossloads.csv"))
    write.csv(CSF_info$sigmain[["xloading_main_sig"]],paste0(savepath2,"sig_func_loads.csv"))
    
    write.csv(CSF_info$sigmain[["zcrossx_main_sig"]],paste0(savepath2,"sig_struc_crossloads.csv"))
    write.csv(CSF_info$sigmain[["zloading_main_sig"]],paste0(savepath2,"sig_struc_loads.csv"))
    
    write.csv(CSF_info$main[["corr"]],paste0(savepath2,"/maincorr.csv"))
    
    
  }
  # write.csv(aa$boot$,paste0(savepath,"_1000_comp",comp,"_mainloads.csv"))
}


