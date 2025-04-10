## 040925
## Run all patients
##############################
library(R.matlab)
library(readxl)
library(parallel)
library(doSNOW)

###########################################
# ---- Analysis script --------------------------------------------------------

pts <- dipsaus::parse_svec("2,4,13")
pipeline_xls <- readxl::read_xlsx("/Volumes/bigbrain/ACL/RAVE_Projects/PipelineScripts/FragilityEEGDataset_pipeline_update_040825.xlsx")

pathdataBHI<-'/Volumes/bigbrain/ACL/RAVE_Projects/BHIDataPackage/'
pathres<-'/Volumes/bigbrain/ACL/RAVE_Projects/Fragility/ResEZFragility_125_75/'

# all epoch data has been resampled to this frequency
fs=500
  
for(i in pts){
  
  patname <- pipeline_xls$subject[i]
  print(patname)

  ## add channel names to the rows
  goodChannels <- dipsaus::parse_svec(pipeline_xls$good_electrodes[i])
  sozChannels  <- dipsaus::parse_svec(pipeline_xls$soz[i])
  channelfile<-paste(pathdataBHI,patname,"IctalEcoGChannels.xls",sep="")
  channelNames <- read_excel(channelfile)
  sozIndex<-which(goodChannels%in%sozChannels==TRUE)
  sozNames<-channelNames$name[sozChannels]
  
  
  ictal_runs <- dipsaus::parse_svec(pipeline_xls$ictal_runs[i])
  #ictal_runs<-1
  for(j in ictal_runs){
  
  #j=1
    
    print(paste("seizure",as.character(j)))
    
  #datafile<-paste(pathdataBHI,patname,"_seizure",as.character(j),"m30sp30s.mat",sep="")
    datafile<-paste(pathdataBHI,patname,"_seizure",as.character(j),"m10sp10s.mat",sep="")
    datamat<-readMat(datafile)
    #ptEpochRaw<-datamat$data[,10001:20001]
    ptEpochRaw<-datamat$data
    
    rownames(ptEpochRaw) <- channelNames$name[goodChannels]
    
    ## Add time stamps to the columns
    times <- seq(-10, 10, length.out=ncol(ptEpochRaw))
    times_with_sign <- ifelse(times >= 0, paste0("+", times), as.character(times))
    colnames(ptEpochRaw)<-times_with_sign 
    
    ptEcoGt<-ptEpochRaw
    attr(ptEcoGt, "sozIndex") <- sozIndex
    attr(ptEcoGt, "sozNames") <- sozNames
  
    epocht <- Epoch(ptEcoGt)
    window <- 125
    step <- 75
   
    cl <- parallel::makeCluster(4, type = "SOCK")
    doSNOW::registerDoSNOW(cl)
    
    title <- paste(patname,'seizure',as.character(j))
    fragtest<-calcAdjFrag(
      epoch = epocht, window = window,
      step = step, parallel = TRUE, progress = TRUE
    )
    
    
    ## stop the parallel backend
    parallel::stopCluster(cl)
    
    
    R2pt<-fragtest$R2
    lambdapt<-fragtest$lambdas
    fragpt<-fragtest$frag
    fragrankpt<-fragtest$frag_ranked
    startT<-fragtest$startTimes
    colnames(fragpt) <-startT
    colnames(fragrankpt) <-startT
    colnames(R2pt) <-startT
    
    write.csv(fragpt,paste(pathres,"fragmat_",patname,"sz",as.character(j),'_',as.character(window),'_',as.character(step),'.csv',sep=""))
    write.csv(fragrankpt,paste(pathres,"fragmatrank_",patname,"sz",as.character(j),'_',as.character(window),'_',as.character(step),'.csv',sep=""))
    write.csv(R2pt,paste(pathres,"R2_",patname,"sz",as.character(j),'_',as.character(window),'_',as.character(step),'.csv',sep=""))
    write.csv(lambdapt,paste(pathres,"lambda_",patname,"sz",as.character(j),'_',as.character(window),'_',as.character(step),'.csv',sep=""))
 
    ## Compute quantiles, mean and standard deviation for two electrodes group marked as soz non marked as soz
    ptfragstat <- fragStat(fragtest, sozIndex)
    quantilept<-ptfragstat$qmatrix
    colnames(quantilept) <-startT
    
    distribpt<-data.frame(ptfragstat$meanSOZ,ptfragstat$sdSOZ,ptfragstat$meanRef,ptfragstat$sdRef)
    rownames(distribpt)<-startT
    
    write.csv(quantilept,paste(pathres,"quantile_",patname,"sz",as.character(j),'_',as.character(window),'_',as.character(step),'.csv',sep=""))   
    write.csv(distribpt,paste(pathres,"distrib_",patname,"sz",as.character(j),'_',as.character(window),'_',as.character(step),'.csv',sep=""))   
    
    ## Result visualization
  
    #plotheatmap<-plotFragHeatmapranked(fragtest, sozIndex)
    plotheatmap<-plotFragHeatmap(fragtest, sozIndex)
    plotheatmap<-plotheatmap+ggplot2::ggtitle(paste("Fragility heatmap for ",title))
    # Add a vertical line at Seizure Onset
    #plotheatmap <- plotheatmap + ggplot2::geom_vline(xintercept = as.numeric(0), color = "black", linetype = "dashed", size = 1)
    plotheatmap
    ggplot2::ggsave(paste(pathres,'FragilityHeatmap_',patname,'sz',as.character(j),'_',as.character(window),'_',as.character(step),'.png',sep=""))
    
    plotDistribution<-plotFragDistribution(fragtest, sozIndex)
    plotDistribution<-plotDistribution+ggplot2::ggtitle((paste("Pooled fragility distribution for ",title)))
    plotDistribution
    ggplot2::ggsave(paste(pathres,'FragilityDistribution_',patname,'sz',as.character(j),'_',as.character(window),'_',as.character(step),'.png',sep=""))
  
    plotquantile<-plotFragQuantile(fragtest, sozIndex)
    plotquantile<-plotquantile+ggplot2::ggtitle(paste("Pooled fragility quantiles for ",title))
    plotquantile
    ggplot2::ggsave(paste(pathres,'FragilityQuantile_',patname,'sz',as.character(j),'_',as.character(window),'_',as.character(step),'.png',sep=""))
    
  
    
  }

}
