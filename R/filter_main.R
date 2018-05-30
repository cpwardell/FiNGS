## This is the main method for running filters

## Load packages
suppressMessages(library(VariantAnnotation))

## Suppress all warnings
options(warn=-1)

## Store time script was invoked for output names
thetime=strftime(Sys.time(),"%Y%m%d%H%M%S")

## Collect arguments
args = commandArgs(TRUE)
codedir=args[1]
parameters=args[2]
resultsdir=args[3]
pathtotumor=args[4]
pathtonormal=args[5]
pathtovcf=args[6]
pass=args[7]
passonly=args[8]

if(parameters=="default"){
  parameters=paste0(codedir,"/R/filter_parameters.txt")
}

## Source functions
source(paste0(codedir,"/R/filter_functions.R"))

## Read in filtering parameters
params=read.table(parameters,header=F,row.names=1)

# Read in mutect data for exome
pdata=filterreader(pathtotumor,pathtonormal)



### Perform filtering

## Write plots to a single PDF
pdf(paste0(resultsdir,"/filter.plots.",thetime,".pdf"), onefile=TRUE,paper="letter")

pdata=filterminimumdepth(pdata,params["minimumdepth",])
pdata=filtermaximumdepth(pdata,params["maximumdepth",])
pdata=filterminaltcount(pdata,params["minaltcount",])
pdata=filterminbasequality(pdata,params["minbasequality",])
pdata=filterzeroproportion(pdata,params["zeroproportion",])
pdata=filterstrandbias(pdata,params["strandbias",])
pdata=filterminmapquality(pdata,params["minmapquality",])
pdata=filterminmapqualitydifference(pdata,params["minmapqualitydifference",])
pdata=filterenddistance(pdata,params["enddistance",])
pdata=filterenddistancemad(pdata,params["enddistancemad",])
pdata=filtereditdistance(pdata,params["editdistance",])
pdata=filtermaxvafnormal(pdata,params["maxvafnormal",])
pdata=filterminvaftumor(pdata,params["minvaftumor",])
pdata=filtermaxoaftumor(pdata,params["maxoaftumor",])
pdata=filtermaxaltsecondtumor(pdata,params["maxsecondtumor",])
pdata=filtermaxrefbadorientnormal(pdata,params["maxbadorient",])
pdata=filterfoxog(pdata,params["foxog",])

## Turn off graphics device
graphics.off()

## Write the file containing everything...
write.table(pdata,file=paste0(resultsdir,"/filtered.results.",thetime,".txt"),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

## Read the original VCF, and write a filtered version
invcf=readVcf(pathtovcf,genome="hg")

## If --PASS is set, cut down the VCF to only PASSing variants
if(pass=="True"){
  invcf=invcf[fixed(invcf)$FILTER%in%"PASS"]
}

## Set up FILTER values for VCF
filtercols=grep("filter",colnames(pdata),value=T)
filtervals=rep("",nrow(pdata))
for(i in 1:length(filtercols)){
  filtervals[!pdata[,filtercols[i]]]=paste(filtervals[!pdata[,filtercols[i]]],filtercols[i],sep=";")
}
## Clean up any starting semicolons
filtervals=sub("^;","",filtervals)

## Set any row that passed all filters to PASS
filtervals[filtervals==""]="PASS"

## Edit VCF; update header and change filter values
fixed(invcf)$FILTER=filtervals

## If --PASSonly is set, strip VCF to only PASSing variants
if(passonly=="True"){
  invcf=invcf[fixed(invcf)$FILTER%in%"PASS"]
}


#1 define filter descriptions and names
filterheader=data.frame(Description=c(
  paste("Minimum depth in tumor:",params["minimumdepth",]),
  paste("Minimum depth in normal:",params["minimumdepth",]),
  paste("Maximum depth in normal:",params["maximumdepth",]),
  paste("Maximum depth in tumor:",params["maximumdepth",]),
  paste("Minimum number of ALT reads in tumor:",params["minaltcount",]),
  paste("Minimum median base quality of ALT reads in tumor:",params["minbasequality",]),
  paste("Minimum median base quality of REF reads in tumor:",params["minbasequality",]),
  paste("Minimum median base quality of REF reads in normal:",params["minbasequality",]),
  paste("Maximum proportion of zero mapping quality reads in tumor:",params["zeroproportion",]),
  paste("Maximum proportion of zero mapping quality reads in normal:",params["zeroproportion",]),
  paste("Strand bias exclusion proportion:",params["strandbias",]),
  paste("Minimum median mapping quality of ALT reads in tumor:",params["minmapquality",]),
  paste("Maximum difference between median mapping quality of ALT reads in tumor and REF reads in normal:",params["minmapqualitydifference",]),
  paste("Maximum median shortest distance to either aligned end in tumor:",params["enddistance",]),
  paste("Minimum MAD of ALT position in tumor:",params["enddistancemad",]),
  paste("Maximum edit distance of ALT reads in tumor:",params["editdistance",]),
  paste("Maximum edit distance of REF reads in tumor:",params["editdistance",]-1),
  paste("Maximum VAF in normal:",params["maxvafnormal",]),
  paste("Minimum VAF in tumor:",params["minvaftumor",]),
  paste("Maximum OAF in tumor:",params["maxoaftumor",]),
  paste("Maximum proportion of secondary alignments in tumor:",params["maxsecondtumor",]),
  paste("Maximum proportion of inversion orientation reads in normal:",params["maxbadorient",]),
  paste("FoxoG artefact proportion:",params["foxog",])))
rownames(filterheader)=filtercols

#2 Assign these to the relevant slot of the VCF
fixed(header(invcf))$FILTER=filterheader

## Write the finalized VCF
writeVcf(invcf,filename=paste0(resultsdir,"/filtered.results.",thetime,".vcf"))

