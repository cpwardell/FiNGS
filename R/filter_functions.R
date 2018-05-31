## This file should be sourced and contains all the functions necessary to 
## create/plot/filter data

## #rugs are commented out; it takes ~4 minutes to draw a GIAB plots with 4 million points...
## but only 5 seconds with no #rugs

## Function to read in filter data
filterreader=function(pathtotumor,pathtonormal){
  ## Column names for all filter output
  cnames=c("UID","CHR","POS","REF","ALT","refcount","altcount","varianttype","depth","vaf","raf","oaf",
           "medianbaseq","medianbaseqref","medianbaseqalt","medianmapq","medianmapqref","medianmapqalt",
           "zeros","zerospersite","softreadlengthsrefmean","softreadlengthsaltmean","goodoffsetproportion",
           "distancetoend1median","mad1","distancetoend2median","mad2","distancetoend1medianref","madref1",
           "distancetoend2medianref","madref2","distancetoend1medianalt","madalt1","distancetoend2medianalt",
           "madalt2","shortestdistancetoendmedian","madaltshort","sb","gsb","fishp","F1R2","F2R1","FoxoG","refld","altld",
           "refsecondprop","altsecondprop","refbadorientationprop","altbadorientationprop","refmatecontigcount",
           "altmatecontigcount")
  tumor=read.table(pathtotumor,header=F,stringsAsFactors=F)  
  normal=read.table(pathtonormal,header=F,stringsAsFactors=F)
  colnames(tumor)=paste0(cnames[1:ncol(tumor)],".tumor")
  colnames(normal)=paste0(cnames[1:ncol(normal)],".normal")  
  all=cbind(tumor,normal)
  ## Add in mutation type
  all$sixtypes=NA
  all[all$REF.tumor=="T"&all$ALT.tumor=="G","sixtypes"]="T>G/A>C"
  all[all$REF.tumor=="A"&all$ALT.tumor=="C","sixtypes"]="T>G/A>C"
  all[all$REF.tumor=="T"&all$ALT.tumor=="C","sixtypes"]="T>C/A>G"
  all[all$REF.tumor=="A"&all$ALT.tumor=="G","sixtypes"]="T>C/A>G"
  all[all$REF.tumor=="T"&all$ALT.tumor=="A","sixtypes"]="T>A/A>T"
  all[all$REF.tumor=="A"&all$ALT.tumor=="T","sixtypes"]="T>A/A>T"
  all[all$REF.tumor=="C"&all$ALT.tumor=="T","sixtypes"]="C>T/G>A"
  all[all$REF.tumor=="G"&all$ALT.tumor=="A","sixtypes"]="C>T/G>A"
  all[all$REF.tumor=="C"&all$ALT.tumor=="G","sixtypes"]="C>G/G>C"
  all[all$REF.tumor=="G"&all$ALT.tumor=="C","sixtypes"]="C>G/G>C"
  all[all$REF.tumor=="C"&all$ALT.tumor=="A","sixtypes"]="C>A/G>T"
  all[all$REF.tumor=="G"&all$ALT.tumor=="T","sixtypes"]="C>A/G>T"
  return(all)
}

############ Recall/precision and F1
#TP=true and PASS
#FP=not true and PASS
#TN=not true and FAIL
#FN=true and FAIL

#Sensitivity = TP/(TP+FN) # aka "Recall", the proportion of true positives recovered
#Specificity = TN/(TN+FP)
#Precision = TP/(TP+FP) # the proportion of positives that are true positives

senspec=function(TP,FP,TN,FN){
  sensitivity = round(TP/(TP+FN),3)
  precision = round(TP/(TP+FP),3)
  #specificity = round(TN/(TN+FP),3)
  fscore=round(2*TP/(2*TP+FP+FN),3)
  
  print(paste("TOTAL is:",TP+FP+TN+FN))
  print(paste("Positives is:",TP+FP))
  print(paste("Negatives is:",TN+FN))
  print(paste("TP is:",TP))
  print(paste("FP is:",FP))
  print(paste("TN is:",TN))
  print(paste("FN is:",FN))
  print(paste("Recall is:",sensitivity))
  print(paste("Precision is:",precision))
  print(paste("F score is:",fscore))
  print("")
  #print(paste("Specificity is:",specificity))
}

########## MINIMUM DEPTH ###################
#normal depth<=10
#tumor depth<=10
filterminimumdepth=function(exome,minimumdepth){
  exome$filter.depth.minimum.tumor=(exome$depth.tumor>=minimumdepth)
  exome$filter.depth.minimum.normal=(exome$depth.normal>=minimumdepth)
  
  ## Plots for this filter
  xlims=range(pretty(exome$depth.tumor,nclass.Sturges(exome$depth.tumor)))
  hist(exome$depth.tumor,breaks=100,main="Minimum depth (tumor)",xlim=xlims,xlab="Depth")
  par(new=T)
  try(plot(density(exome$depth.tumor,na.rm=T),col="black",lwd=2,
       main="",ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F,xlim=xlims),silent=TRUE)
  abline(v=minimumdepth,lwd=2,col="black",lty=2)

  xlims=range(pretty(exome$depth.normal,nclass.Sturges(exome$depth.normal)))
  hist(exome$depth.normal,breaks=100,main="Minimum depth (normal)",xlim=xlims,xlab="Depth")
  par(new=T)
  try(plot(density(exome$depth.normal,na.rm=T),col="black",lwd=2,
       main="",ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F,xlim=xlims),silent=TRUE)
  abline(v=minimumdepth,lwd=2,col="black",lty=2)

  return(exome)
}
########## MAXIMUM DEPTH ###################
#maximum depth=mean(x)+sd(x)
filtermaximumdepth=function(exome,maxdepth){
  #maxdepthn=(mean(exome$depth.normal)+sd(exome$depth.normal))
  #maxdeptht=(mean(exome$depth.tumor)+sd(exome$depth.tumor))
  maxdepthn=maxdepth
  maxdeptht=maxdepth
  exome$filter.depth.maximum.normal=(exome$depth.normal<maxdepthn)
  exome$filter.depth.maximum.tumor=(exome$depth.tumor<maxdeptht)
  
  ## Plots for this filter
  xlims=range(pretty(exome$depth.tumor,nclass.Sturges(exome$depth.tumor)))
  hist(exome$depth.tumor,breaks=100,main="Maximum depth (tumor)",xlim=xlims,xlab="Depth")
  par(new=T)
  try(plot(density(exome$depth.tumor,na.rm=T),col="black",lwd=2,
       main="",ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F,xlim=xlims),silent=TRUE)
  abline(v=maxdeptht,lwd=2,col="black",lty=2)

  xlims=range(pretty(exome$depth.normal,nclass.Sturges(exome$depth.normal)))
  hist(exome$depth.normal,breaks=100,main="Maximum depth (normal)",xlim=xlims,xlab="Depth")
  par(new=T)
  try(plot(density(exome$depth.normal,na.rm=T),col="black",lwd=2,
       main="",ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F,xlim=xlims),silent=TRUE)
  abline(v=maxdepthn,lwd=2,col="black",lty=2)

  return(exome)
}

########### Minimum number of variant supporting reads in tumor ################
## Must have at least 3 supporting reads at any depth
filterminaltcount=function(exome,minimumaltcount){
  exome$filter.altcount.minimum.tumor=(exome$altcount.tumor>=minimumaltcount)
  
  ## Plots for this filter
  xlims=range(pretty(exome$altcount.tumor,nclass.Sturges(exome$altcount.tumor)))
  hist(exome$altcount.tumor,breaks=100,main="Minimum ALT count (tumor)",xlim=xlims,xlab="ALT count")
  par(new=T)
  try(plot(density(exome$altcount.tumor,na.rm=T),col="black",lwd=2,
       main="",ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F,xlim=xlims),silent=TRUE)
  abline(v=minimumaltcount,lwd=2,col="black",lty=2)

  return(exome)
}

#################### Minimum base quality score ############
## medianbaseq must be greater than 30 for tumor ref and alt reads and normal ref reads
filterminbasequality=function(exome,minimumbaseq){
  exome$filter.medianbaseqalt.minimum.tumor=(exome$medianbaseqalt.tumor>=minimumbaseq)
  exome$filter.medianbaseqref.minimum.tumor=(exome$medianbaseqref.tumor>=minimumbaseq)
  exome$filter.medianbaseqref.minimum.normal=(exome$medianbaseqref.normal>=minimumbaseq)
  
  ## Pick up NA values
  ## NA ALT in tumor means no ALT reads - FAIL
  ## NA REF in tumor means no REF reads; could be homozygous or shallow - PASS
  ## NA REF in normal means no REF reads - FAIL
  exome$filter.medianbaseqalt.minimum.tumor[is.na(exome$filter.medianbaseqalt.minimum.tumor)]=FALSE
  exome$filter.medianbaseqref.minimum.tumor[is.na(exome$filter.medianbaseqref.minimum.tumor)]=TRUE
  exome$filter.medianbaseqref.minimum.normal[is.na(exome$filter.medianbaseqref.minimum.normal)]=FALSE
  
  ## Plots for this filter
  xlims=range(pretty(exome$medianbaseqalt.tumor,nclass.Sturges(exome$medianbaseqalt.tumor)))
  hist(exome$medianbaseqalt.tumor,breaks=100,main="Minimum base quality for ALT (tumor)",xlim=xlims,xlab="Base quality")
  par(new=T)
  try(plot(density(exome$medianbaseqalt.tumor,na.rm=T),col="black",lwd=2,
       main="",ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F,xlim=xlims),silent=TRUE)
  abline(v=minimumbaseq,lwd=2,col="black",lty=2)

  xlims=range(pretty(exome$medianbaseqref.tumor,nclass.Sturges(exome$medianbaseqref.tumor)))
  hist(exome$medianbaseqref.tumor,breaks=100,main="Minimum base quality for REF (tumor)",xlim=xlims,xlab="Base quality")
  par(new=T)
  try(plot(density(exome$medianbaseqref.tumor,na.rm=T),col="black",lwd=2,
       main="",ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F,xlim=xlims),silent=TRUE)
  abline(v=minimumbaseq,lwd=2,col="black",lty=2)

  xlims=range(pretty(exome$medianbaseqref.normal,nclass.Sturges(exome$medianbaseqref.normal)))
  hist(exome$medianbaseqref.normal,breaks=100,main="Minimum base quality for REF (normal)",xlim=xlims,xlab="Base quality")
  par(new=T)
  try(plot(density(exome$medianbaseqref.normal,na.rm=T),col="black",lwd=2,
       main="",ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F,xlim=xlims),silent=TRUE)
  abline(v=minimumbaseq,lwd=2,col="black",lty=2)

  return(exome)
}

######### Zero mapping quality reads filter #######################
# Proportion of zero mapping quality reads must be less than 10% 
filterzeroproportion=function(exome,zeroproportion){
  exome$filter.zerospersite.maximum.tumor=(exome$zerospersite.tumor<zeroproportion)
  exome$filter.zerospersite.maximum.normal=(exome$zerospersite.normal<zeroproportion)
  
  ## Pick up NA values.  These mean zero depth, so fail them
  exome$filter.zerospersite.maximum.tumor[is.na(exome$filter.zerospersite.maximum.tumor)]=FALSE
  exome$filter.zerospersite.maximum.normal[is.na(exome$filter.zerospersite.maximum.normal)]=FALSE

  ## Plots for this filter
  xlims=range(pretty(exome$zerospersite.tumor,nclass.Sturges(exome$zerospersite.tumor)))
  hist(exome$zerospersite.tumor,breaks=100,main="Proportion of zero mapq reads (tumor)",xlim=xlims,xlab="Proportion of zero mapq reads")
  par(new=T)
  try(plot(density(exome$zerospersite.tumor,na.rm=T),col="black",lwd=2,
       main="",ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F,xlim=xlims),silent=TRUE)
  abline(v=zeroproportion,lwd=2,col="black",lty=2)

  xlims=range(pretty(exome$zerospersite.normal,nclass.Sturges(exome$zerospersite.normal)))
  hist(exome$zerospersite.normal,breaks=100,main="Proportion of zero mapq reads (normal)",xlim=xlims,xlab="Proportion of zero mapq reads")
  par(new=T)
  try(plot(density(exome$zerospersite.normal,na.rm=T),col="black",lwd=2,
       main="",ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F,xlim=xlims),silent=TRUE)
  abline(v=zeroproportion,lwd=2,col="black",lty=2)

  return(exome)
}

######### Strand bias filter #######################
# Maximum strand bias = 1.25
# Top proportion = 10%
filterstrandbias=function(exome,maximumstrandbias){
  exome$filter.sb.tenpc.tumor=(exome$sb.tumor<quantile(exome$sb.tumor,prob=1-maximumstrandbias,na.rm=T))
  
  ## Pick up NA values.  Pass anything that didn't calculate
  exome$filter.sb.tenpc.tumor[is.na(exome$filter.sb.tenpc.tumor)]=TRUE
  
  ## Plots for this filter
  xlims=range(pretty(exome$sb.tumor,nclass.Sturges(exome$sb.tumor)))
  hist(exome$sb.tumor,breaks=100,main="Strand bias (top 10%) (tumor)",xlim=xlims,xlab="Strand bias")
  par(new=T)
  try(plot(density(exome$sb.tumor,na.rm=T),col="black",lwd=2,
       main="",ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F,xlim=xlims),silent=TRUE)
  abline(v=quantile(exome$sb.tumor,prob=0.9,na.rm=T),lwd=2,col="black",lty=2)

  return(exome)
}

#################### Minimum mapping quality score ############
## Median mapping quality of alt reads in the tumor must be greater than 50
filterminmapquality=function(exome,minimummapq){
  exome$filter.medianmapqalt.minimum.tumor=(exome$medianmapqalt.tumor>=minimummapq)
  
  ## Pick up NA values.  These mean zero depth, so fail them
  exome$filter.medianmapqalt.minimum.tumor[is.na(exome$filter.medianmapqalt.minimum.tumor)]=FALSE
  
  ## Plots for this filter
  xlims=range(pretty(exome$medianmapqalt.tumor,nclass.Sturges(exome$sb.tumor)))
  hist(exome$medianmapqalt.tumor,breaks=100,main="Minimum median mapping quality (tumor)",xlim=xlims,xlab="Mapping quality")
  par(new=T)
  try(plot(density(exome$medianmapqalt.tumor,na.rm=T),col="black",lwd=2,
       main="",ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F,xlim=xlims),silent=TRUE)
  abline(v=minimummapq,lwd=2,col="black",lty=2)

  return(exome)
}

#################### Maximum mapping quality difference ############
## Median mapping quality of alt reads in the tumor must be within 5 of
## median mapping quality of ref reads in the normal
filterminmapqualitydifference=function(exome,maximumdifference){
  exome$filter.medianmapq.difference=(abs(exome$medianmapqalt.tumor-exome$medianmapqref.normal)<maximumdifference)
  
  ## Pick up NA values.  These mean zero depth, so fail them
  exome$filter.medianmapq.difference[is.na(exome$filter.medianmapq.difference)]=FALSE
  
  ## Plots for this filter
  xlims=range(pretty(abs(exome$medianmapqalt.tumor-exome$medianmapqref.normal),nclass.Sturges(abs(exome$medianmapqalt.tumor-exome$medianmapqref.normal))))
  hist(abs(exome$medianmapqalt.tumor-exome$medianmapqref.normal),breaks=100,main="Maximum median mapping quality difference (tumor/normal)",xlim=xlims,xlab="Mapping quality difference")
  par(new=T)
  try(plot(density(abs(exome$medianmapqalt.tumor-exome$medianmapqref.normal),na.rm=T),col="black",lwd=2,
       main="",ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F,xlim=xlims),silent=TRUE)
  abline(v=maximumdifference,lwd=2,col="black",lty=2)

  return(exome)
}

#################### Maximum distance to alignment end ############
#The median shortest distance of the variant position within the read to either aligned end is less than 10 in the tumor
filterenddistance=function(exome,distance){
  exome$filter.medianaltenddistance.tumor=(exome$shortestdistancetoendmedian.tumor>distance)
  
  ## Pick up NA values.  These mean zero depth, so fail them
  exome$filter.medianaltenddistance.tumor[is.na(exome$filter.medianaltenddistance.tumor)]=FALSE
  
  ## Plots for this filter
  xlims=range(pretty(exome$shortestdistancetoendmedian.tumor,nclass.Sturges(exome$shortestdistancetoendmedian.tumor)))
  hist(exome$shortestdistancetoendmedian.tumor,breaks=100,main="Shortest median distance to end of read (tumor)",xlim=xlims,xlab="Median distance to end of read")
  par(new=T)
  try(plot(density(exome$shortestdistancetoendmedian.tumor,na.rm=T),col="black",lwd=2,
       main="",ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F,xlim=xlims),silent=TRUE)
  abline(v=distance,lwd=2,col="black",lty=2)

  return(exome)
}

#################### Maximum distance to alignment end ############
#The median shortest distance of the variant position within the read to either aligned end is less than 10 in the tumor
filterenddistancemad=function(exome,mad){
  exome$filter.medianaltenddistance.mad.tumor=(exome$madaltshort.tumor>mad)
  
  ## Pick up NA values.  These mean zero depth, so fail them
  exome$filter.medianaltenddistance.mad.tumor[is.na(exome$filter.medianaltenddistance.mad.tumor)]=FALSE
  
  ## Plots for this filter
  xlims=range(pretty(exome$madaltshort.tumor,nclass.Sturges(exome$madaltshort.tumor)))
  hist(exome$madaltshort.tumor,breaks=100,main="MAD of shortest median distance to end of read (tumor)",xlim=xlims,xlab="MAD of median distance to end of read")
  par(new=T)
  try(plot(density(exome$madaltshort.tumor,na.rm=T),col="black",lwd=2,
       main="",ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F,xlim=xlims),silent=TRUE)
  abline(v=mad,lwd=2,col="black",lty=2)

  return(exome)
}

#################### Maximum edit distance  ############
#The maximum median edit distance between tumor reads and the reference geneome must be <=4
filtereditdistance=function(exome,editdistance){
  exome$filter.editdistancealt.tumor=(exome$altld.tumor<=editdistance)
  exome$filter.editdistanceref.tumor=(exome$refld.tumor<=editdistance-1)
  
  ## Pick up NA values.  Fail zero alt records, pass zero ref records
  exome$filter.editdistancealt.tumor[is.na(exome$filter.editdistancealt.tumor)]=FALSE
  exome$filter.editdistanceref.tumor[is.na(exome$filter.editdistanceref.tumor)]=TRUE
  
  ## Plots for this filter
  xlims=range(pretty(exome$altld.tumor,nclass.Sturges(exome$altld.tumor)))
  hist(exome$altld.tumor,breaks=100,main="Maximum edit distance between ALT reads and ref genome (tumor)",xlim=xlims,xlab="Edit distance")
  par(new=T)
  try(plot(density(exome$altld.tumor,na.rm=T),col="black",lwd=2,
       main="",ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F,xlim=xlims),silent=TRUE)
  abline(v=editdistance,lwd=2,col="black",lty=2)

  xlims=range(pretty(exome$refld.tumor,nclass.Sturges(exome$refld.tumor)))
  hist(exome$refld.tumor,breaks=100,main="Maximum edit distance between REF reads and ref genome (tumor)",xlim=xlims,xlab="Edit distance")
  par(new=T)
  try(plot(density(exome$refld.tumor,na.rm=T),col="black",lwd=2,
       main="",ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F,xlim=xlims),silent=TRUE)
  abline(v=editdistance-1,lwd=2,col="black",lty=2)

  return(exome)
}

#################### Maximum VAF allowed in normal ############
#The maximum VAF in the normal is less than 0.05
filtermaxvafnormal=function(exome,vaf){
  exome$filter.maxvaf.normal=(exome$vaf.normal<vaf)
  
  ## Pick up NA values.  These mean zero depth, so fail them
  exome$filter.maxvaf.normal[is.na(exome$filter.maxvaf.normal)]=FALSE
  
  ## Plots for this filter
  xlims=range(pretty(exome$vaf.normal,nclass.Sturges(exome$vaf.normal)))
  hist(exome$vaf.normal,breaks=100,main="Maximum VAF (normal)",xlim=xlims,xlab="VAF")
  par(new=T)
  try(plot(density(exome$vaf.normal,na.rm=T),col="black",lwd=2,
       main="",ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F,xlim=xlims),silent=TRUE)
  abline(v=vaf,lwd=2,col="black",lty=2)

  return(exome)
}

#################### Minimum VAF in tumor ############
#The minimum VAF in the tumor is greater than 0.05
filterminvaftumor=function(exome,vaf){
  exome$filter.minvaf.tumor=(exome$vaf.tumor>=vaf)
  
  ## Pick up NA values.  These mean zero depth, so fail them
  exome$filter.minvaf.tumor[is.na(exome$filter.minvaf.tumor)]=FALSE
  
  ## Plots for this filter
  xlims=range(pretty(exome$vaf.tumor,nclass.Sturges(exome$vaf.tumor)))
  hist(exome$vaf.tumor,breaks=100,main="Minimum VAF (tumor)",xlim=xlims,xlab="VAF")
  par(new=T)
  try(plot(density(exome$vaf.tumor,na.rm=T),col="black",lwd=2,
       main="",ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F,xlim=xlims),silent=TRUE)
  abline(v=vaf,lwd=2,col="black",lty=2)

  return(exome)
}

#################### Maximum OAF allowed in tumor ############
#The maximum OAF in the normal is less than 0.05
filtermaxoaftumor=function(exome,vaf){
  exome$filter.maxoaf.tumor=(exome$oaf.tumor<vaf)
  
  ## Pick up NA values.  These mean zero depth, so fail them
  exome$filter.maxoaf.tumor[is.na(exome$filter.maxoaf.tumor)]=FALSE
  
  ## Plots for this filter
  xlims=range(pretty(exome$oaf.tumor,nclass.Sturges(exome$oaf.tumor)))
  hist(exome$oaf.tumor,breaks=100,main="Maximum OAF (tumor)",xlim=xlims,xlab="OAF")
  par(new=T)
  try(plot(density(exome$oaf.tumor,na.rm=T),col="black",lwd=2,
       main="",ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F,xlim=xlims),silent=TRUE)
    abline(v=vaf,lwd=2,col="black",lty=2)

  return(exome)
}

#################### Maximum secondary alignments allowed in ALT reads in tumor ############
#The maximum proportion of secondary alignments in ALT reads in the tumor is less than 0.05
filtermaxaltsecondtumor=function(exome,maxsecond){
  exome$filter.maxsecond.tumor=(exome$altsecondprop.tumor<maxsecond)
  
  ## Pick up NA values.  These mean zero depth, so fail them
  exome$filter.maxsecond.tumor[is.na(exome$filter.maxsecond.tumor)]=FALSE
  
  ## Plots for this filter
  xlims=range(pretty(exome$altsecondprop.tumor,nclass.Sturges(exome$altsecondprop.tumor)))
  hist(exome$altsecondprop.tumor,breaks=100,main="Maximum secondary alignments (tumor)",xlim=xlims,xlab="Proportion of secondary alignments")
  par(new=T)
  try(plot(density(exome$altsecondprop.tumor,na.rm=T),col="black",lwd=2,
       main="",ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F,xlim=xlims),silent=TRUE)
  abline(v=maxsecond,lwd=2,col="black",lty=2)

  return(exome)
}

#################### Maximum bad orientation allowed in REF reads in normal ############
#The maximum proportion of bad orientation in REF reads in the normal is 0.2
#This indicates low quality sequence/mapping errors
#It would exclude genuine mutations at sites of germline inversions
filtermaxrefbadorientnormal=function(exome,maxbadorient){
  exome$filter.maxbadorient.normal=(exome$refbadorientationprop.normal<maxbadorient)
  
  ## Pick up NA values.  These mean zero depth, so fail them
  exome$filter.maxbadorient.normal[is.na(exome$filter.maxbadorient.normal)]=FALSE
  
  ## Plots for this filter
  xlims=range(pretty(exome$refbadorientationprop.normal,nclass.Sturges(exome$refbadorientationprop.normal)))
  hist(exome$refbadorientationprop.normal,breaks=100,main="Maximum bad orientation REF reads (normal)",xlim=xlims,xlab="Proporotion of bad orientation reads")
  par(new=T)
  try(plot(density(exome$altsecondprop.tumor,na.rm=T),col="black",lwd=2,
       main="",ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F,xlim=xlims),silent=TRUE)
  abline(v=maxbadorient,lwd=2,col="black",lty=2)

  return(exome)
}


## Unused/developing filters

#################### Maximum contigs in REF reads in normal ############
#The maximum number of contigs in REF reads in the normal is 3
#Unmapped reads or a stray mismapped read are fine; 4 or more suggests
#that the region is poorly mappable or frequently mismapped
filtermaxrefcontignormal=function(exome,maxcontig){
  exome$filter.maxcontig.normal=(exome$refmatecontigcount.normal<maxcontig)
  
  ## Pick up NA values.  These mean zero depth, so fail them
  exome$filter.maxcontig.normal[is.na(exome$filter.maxcontig.normal)]=FALSE
  
  ## Plots for this filter
  xlims=range(pretty(exome$refmatecontigcount.normal,nclass.Sturges(exome$refmatecontigcount.normal)))
  hist(exome$refmatecontigcount.normal,breaks=100,main="Maximum contigs in REF reads (normal)",xlim=xlims,xlab="Number of contigs present")
  par(new=T)
  try(plot(density(exome$altsecondprop.tumor,na.rm=T),col="black",lwd=2,
       main="",ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F,xlim=xlims),silent=TRUE)
  abline(v=maxcontig,lwd=2,col="black",lty=2)

  return(exome)
}


#################### FoxoG in tumor ############
#Classify C>A|G>T variants as either artifact or non-artifact
#This is achieved by testing against 2 binomial distributions and selecting
#the lower pvalue. One is centered on 0.5, the other at 0.9 (artifact).
#The tested value is the proportion of reads in a specific orientation
filterfoxog=function(exome,foxog){
  ## Default is all values pass
  exome$filter.foxog=TRUE
  
  ## List of mutation types (sorted) 
  sixtypes=sort(unique(exome$sixtypes))

  ## Pvalues for a mutation belonging to 2 binomial distributions
  exome$phalf=NA
  exome$poxog=NA
  
  ## Separate calcuations for C>A and G>T
  ctype=which(exome$REF.tumor%in%c("C") & exome$sixtypes%in%sixtypes[1])
  for(i in 1:length(ctype)){
    try((exome$phalf[ctype[i]]=binom.test(exome$F2R1.tumor[ctype[i]],exome$F2R1.tumor[ctype[i]]+exome$F1R2.tumor[ctype[i]],p=foxog)$p.value),silent=TRUE)
    try((exome$poxog[ctype[i]]=binom.test(exome$F2R1.tumor[ctype[i]],exome$F2R1.tumor[ctype[i]]+exome$F1R2.tumor[ctype[i]],p=0.5)$p.value),silent=TRUE)
  }
  
  gtype=which(exome$REF.tumor%in%c("G") & exome$sixtypes%in%sixtypes[1])
  for(i in 1:length(gtype)){
    try((exome$phalf[gtype[i]]=binom.test(exome$F1R2.tumor[gtype[i]],exome$F1R2.tumor[gtype[i]]+exome$F2R1.tumor[gtype[i]],p=foxog)$p.value),silent=TRUE)
    try((exome$poxog[gtype[i]]=binom.test(exome$F1R2.tumor[gtype[i]],exome$F1R2.tumor[gtype[i]]+exome$F2R1.tumor[gtype[i]],p=0.5)$p.value),silent=TRUE)
  }
  
  ## Fail any row where the variant is more likely to belong to the OxoG binomial
  exome$filter.foxog[exome$poxog<exome$phalf]=FALSE
  
  ## Some plots
  ## pre-filtering profile
  try(plot(density(exome$FoxoG.tumor[exome$sixtypes%in%"C>A/G>T"],na.rm=T),
       main="FoxoG density pre-filtering",xlab="Artifact orientation proportion",lwd=2,xlim=c(0,1)),silent=TRUE)
  for(i in 2:length(sixtypes)){
    try(lines(density(exome$FoxoG.tumor[exome$sixtypes%in%sixtypes[i]],na.rm=T),col=rainbow(6)[i]),silent=TRUE)
  }
  legend("topleft",legend=sixtypes,col=c("black",rainbow(6)[2:6]),pch=15)
  
  ## post filtering profile
  try(plot(density(exome$FoxoG.tumor[exome$sixtypes%in%"C>A/G>T" & exome$filter.foxog],na.rm=T),
       main="FoxoG density post-filtering",xlab="Artifact orientation proportion",lwd=2,xlim=c(0,1)),silent=TRUE)
  for(i in 2:length(sixtypes)){
    try(lines(density(exome$FoxoG.tumor[exome$sixtypes%in%sixtypes[i] & exome$filter.foxog],na.rm=T),col=rainbow(6)[i]),silent=TRUE)
  }
  legend("topleft",legend=sixtypes,col=c("black",rainbow(6)[2:6]),pch=15)
  
  ## barplot
  try(barplot(table(exome$sixtypes),ylab="SNVs",main="Mutation type counts pre-FoxoG-filtering",col=c("black",rainbow(6)[2:6]),las=3),silent=TRUE)
  try(barplot(table(exome$sixtypes[exome$filter.foxog]),ylab="SNVs",main="Mutation type counts post-FoxoG-filtering",col=c("black",rainbow(6)[2:6]),las=3),silent=TRUE)
  
  return(exome)
}

## Filters out any variants that are not in the provided bedfile
filterontarget=function(exome,bedfile){
  ## Default is all values pass
  exome$filter.ontarget=TRUE
  
  ## If no bedfile is provided, just return
  #if(is.na(bedfile)){return(exome)}
  if(bedfile=="None"){return(exome)}
 
  ## Read in bedfile
  bed=read.table(bedfile,header=F,stringsAsFactors=F)
  
  ## Cast bedfile to genomic ranges object
  ## We add 1 to the positions because bed files have 0-based coordinates, VCFs are 1-based
  grbed=GRanges(seqnames=Rle(bed[,1]),ranges=IRanges(bed[,2]+1,end=bed[,3]+1))
  
  ## Cast variant coordinates to genomic ranges object
  grexome = GRanges(seqnames=Rle(exome$CHR.tumor),ranges=IRanges(exome$POS.tumor,end=exome$POS.tumor))

  ## Execute the filter - we suppress warnings abouts missing contigs in one or the other
  exome$filter.ontarget=suppressWarnings(grexome%within%grbed)
  
  return(exome)
}

## Filters out any variants that are not uniquely alignable
filteralign=function(exome,bedfile){
  ## Default is all values pass
  exome$filter.alignability=TRUE
  
  ## If no bedfile is provided, just return
  #if(is.na(bedfile)){return(exome)}
  if(bedfile=="None"){return(exome)}
  
  ## Read in bedfile
  bed=read.table(bedfile,header=F,stringsAsFactors=F)
  
  ## Cast bedfile to genomic ranges object
  ## We add 1 to the positions because bed files have 0-based coordinates, VCFs are 1-based
  grbed=GRanges(seqnames=Rle(bed[,1]),ranges=IRanges(bed[,2]+1,end=bed[,3]+1))
  
  ## Cast variant coordinates to genomic ranges object
  grexome = GRanges(seqnames=Rle(exome$CHR.tumor),ranges=IRanges(exome$POS.tumor,end=exome$POS.tumor))
  
  ## Execute the filter - we suppress warnings abouts missing contigs in one or the other
  exome$filter.alignability=suppressWarnings(grexome%within%grbed)
  
  return(exome)
}
