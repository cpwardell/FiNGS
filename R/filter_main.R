## This is the main method for running filters

## Load packages
library(VariantAnnotation)

## Collect arguments
args = commandArgs(TRUE)
codedir=args[1]
parameters=args[2]
resultsdir=args[3]
pathtotumor=args[4]
pathtonormal=args[5]
pathtovcf=args[6]

if(parameters=="default"){
  parameters=paste0(codedir,"/R/filter_parameters.txt")
}

## Source functions
source(paste0(codedir,"/R/filter_functions.R"))

## Read in filtering parameters
params=read.table(parameters,header=F,row.names=1)

# Read in mutect data for exome
rawdata=filterreader(pathtotumor,pathtonormal)

## Create a copy of the raw data to edit
pdata=rawdata

## Perform filtering

#source("//144.30.235.61/runs/chris_working_dir/project_snef/SNEF/R/filter_functions.R")
#params=read.table("//144.30.235.61/runs/chris_working_dir/project_snef/SNEF/R/filter_parameters.txt",header=F,row.names=1)

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
pdata=filtermaxrefcontignormal(pdata,params["maxcontig",])

## Write the file containing everything...
write.table(pdata,file=paste0(resultsdir,"/filtered.results.txt"),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

## Read the original VCF, and write a filtered version
invcf=readVcf(pathtovcf,genome="GRCh38")

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

#1 define filter descriptions and names
filterheader=data.frame(Description=filtercols) # descriptions
rownames(filterheader)=c("Minimum depth in tumor of 10x",
              "Minimum depth in normal of 10x",
              "Maximum depth in normal of 1000x",
              "Maximum depth in tumor of 1000x",
              "Minimum number of ALT reads in tumor of 3",
              "Minimum base quality of ALT reads in tumor of 30x",
              "Minimum base quality of REF reads in tumor of 30x",
              "Minimum base quality of REF reads in normal of 30x",
              "Proportion of zero mapping quality reads in tumor must be less than 10%",
              "Proportion of zero mapping quality reads in normal must be less than 10%",
              "Tumor sample strand bias top 10%",
              "Minimum median mapping quality of ALT reads in tumor greater than 50",
              "Difference between median mapping quality of ALT reads in tumor and REF reads in the normal greater than 5",
              "Median shortest distance to either aligned end in tumor less than 10",
              "MAD of ALT position in read less than 3",
              "Maximum edit distance of ALT reads in tumor less than 4",
              "Maximum edit distance of REF reads in tumor less than 3",
              "Maximum VAF in normal of 0.03",
              "Minimum VAF in tumor of 0.00",
              "Maximum OAF in tumor of 0.04")

#2 Assign these to the relevant slot of the VCF
fixed(header(invcf))$FILTER=filterheader

## Write the finalized VCF
writeVcf(invcf,filename=paste0(resultsdir,"/filtered.results.vcf"))
