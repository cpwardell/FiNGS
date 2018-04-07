#!/usr/bin/env python3

## Title: SNEF: Super NExtgen Filter
## Author: CP Wardell, UAMS
## Description: Takes a VCF file and BAM files used to generate it and generates various metrics that will be used to filter the VCF

import os
import vcf
import pysam
import logging
import copy
import argparse
import math
import gzip
import subprocess

from functions.shared_functions import *
from functions.vaf import *
from joblib import Parallel, delayed


#from functions.alignability import *
#from functions.ontarget import *
#from functions.oxog import *
#from functions.directioner import *
#from functions.altlevendist import *

## Gather command line args
parser = argparse.ArgumentParser()
parser.add_argument("-v", type=str, help="absolute path to vcf file",required=True)
parser.add_argument("-t", type=str, help="absolute path to tumor bam",required=True)
parser.add_argument("-n", type=str, help="absolute path to normal bam",required=True)
parser.add_argument("-a", type=str, help="absolute path to alignability track",required=False)
parser.add_argument("-b", type=str, help="absolute path to exome bed file",required=False)
parser.add_argument("-d", type=str, help="absolute path to results directory",required=False,default="snef")
parser.add_argument("-p", type=str, help="absolute path to filtering parameters",required=False,default="default")
parser.add_argument("-c", type=int, help="chunk size to process",required=False,default=1000)
parser.add_argument("-m", type=int, help="maximum read depth to process",required=False,default=1000)
parser.add_argument("-j", type=int, help="number of threads to use (set to -1 to use all threads)",required=False,default=1)
parser.add_argument("--debug", help="Turn debugging statements on",action="store_true")
args = parser.parse_args()

## Assign command line args to more friendly variable names
vcfpath=args.v
tbampath=args.t
nbampath=args.n
alignabilitytrack=args.a
bedfile=args.b
resultsdir=args.d
parameters=args.p
chunksize=args.c
maxdepth=args.m
njobs=args.j

## Create directory to put results in
## Trycatch prevents exception if directory already exists or location is unwritable 
try:
	os.mkdir(resultsdir)
except Exception as e:
	#logging.debug("WARNING: error creating results directory: "+str(e))
	pass

## Turn logging on if desired
if(args.debug):
  logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(levelname)s %(message)s',datefmt='%Y-%m-%d %H:%M:%S',filename=resultsdir+'/debug.log.txt',filemode='w')
  console = logging.StreamHandler()
  console.setLevel(logging.DEBUG)
  formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
  console.setFormatter(formatter)
  logging.getLogger('').addHandler(console)
  logging.debug("Debugging mode enabled")

## Write the command line parameters to the log file
logging.debug("VCF is: "+str(vcfpath))
logging.debug("Tumor bamfile is: "+str(tbampath))
logging.debug("Normal bamfile is: "+str(nbampath))
logging.debug("Alignability file is: "+str(alignabilitytrack))
logging.debug("Bedfile is: "+str(bedfile))
logging.debug("Chunksize is: "+str(chunksize)+" records per chunk")
logging.debug("Maximum read depth is: "+str(maxdepth))

## Set up some test variables
# Small VCF
#vcfpath="/illumina/runs/chris_working_dir/project_mgp/UK_data/GATK.Mutect2.snpeff/_EGAR00001321720_EGAS00001001147_D2CRHACXX_6_744.MarkDuplicates.mdup...annotated.vcf"
#tbampath="/illumina/runs/ExternalData/UK_dataset/EGAS00001001147/_EGAR00001321720_EGAS00001001147_D2CRHACXX_6_744.bam"
#nbampath="/illumina/runs/ExternalData/UK_dataset/EGAS00001001147/_EGAR00001321721_EGAS00001001147_D2BEPACXX_5_676.bam"

# Larger VCF
#vcfpath="/illumina/runs/chris_working_dir/project_mgp/UK_data/GATK.Mutect2.snpeff/_EGAR00001321710_EGAS00001001147_D2FR9ACXX_2_900.MarkDuplicates.mdup...annotated.vcf"
#tbampath="/illumina/runs/ExternalData/UK_dataset/EGAS00001001147/_EGAR00001321710_EGAS00001001147_D2FR9ACXX_2_900.bam"

# Example exome and alignability track
#alignabilitytrack="wgEncodeCrgMapabilityAlign100mer.bedGraph.gz"
#bedfile="agilent_sureselect_v5_plus_custom_baits.bed.gz"


# ## Define filename prefixes for output..
tbasename = os.path.basename(tbampath)
#tfilename = os.path.splitext(tbasename)[0]
tfilename = "tumor"
#nbasename = os.path.basename(nbampath) 
nbasename = "normal"
nfilename = os.path.splitext(nbasename)[0]

## Find size of VCF; vital for chunking and estimating runtime
logging.debug("VCF line counting begun")
nvcflines=vcfcount(vcfpath)
maxchunks=math.ceil(nvcflines/chunksize)
logging.debug("VCF contains "+str(nvcflines)+" records")

## Multithreaded code:
## VCF must be opened once per bamfile

## Tumor bam
vcf_reader = vcf.Reader(open(vcfpath, 'r'))
Parallel(n_jobs=njobs, backend="threading")(delayed(vaf)(vcfchunk,tbampath, resultsdir+"/"+tfilename+"."+str(chunknumber+1)+".txt",chunknumber,maxchunks,maxdepth) for chunknumber,vcfchunk in enumerate(vcfyield(vcf_reader,chunksize,nvcflines)))

## Normal bam
vcf_reader = vcf.Reader(open(vcfpath, 'r'))
Parallel(n_jobs=njobs, backend="threading")(delayed(vaf)(vcfchunk,nbampath, resultsdir+"/"+nfilename+"."+str(chunknumber+1)+".txt",chunknumber,maxchunks,maxdepth) for chunknumber,vcfchunk in enumerate(vcfyield(vcf_reader,chunksize,nvcflines)))
## End of multithreaded code:

## Test for incomplete results - would allow restarting
## Test length of all output files exist and make sure they match chunk size
## If fails, rerun those chunks only.  
## Last chunk will be an exception of nvcflines-(maxchunks*chunksize)
## We could add a call to a "chunkchecker" function that holds a list of chunks to run
## Initially this would be 0 to maxchunks, and on rerun would be the bad chunks only

## Singlethreaded code
#chunknumber=1
#vcf_reader = vcf.Reader(open(vcfpath, 'r')) # regular
#for vcfchunk in vcfyield(vcf_reader, chunksize,nvcflines):
#	logging.debug("Processing chunk "+str(chunknumber)+" of "+str(maxchunks)+" containing "+str(len(vcfchunk))+" records")
	## Clean up VCF by editing CHROM names.  We want "1", not "chr1"
#	vcfchunk=cleanchroms(vcfchunk)
	
	## Produce metrics from tumor and normal bams
	#logging.debug("Launching tumor bam metrics")
	#vaf(vcfchunk,tbampath,resultsdir+"/"+tfilename+"."+str(chunknumber)+".txt")
	#logging.debug("Completed tumor bam metrics")
	#logging.debug("Launching normal bam metrics")
	#vaf(vcfchunk,nbampath,resultsdir+"/"+nfilename+"."+str(chunknumber)+".txt")
	#logging.debug("Completed normal bam metrics")

#	chunknumber+=1

## End of singlethreaded code


## Run tests
## Write complete VCF with PASS or FAIL in column 
## or write trimmed VCF?  Make this an option?


## Combine all chunks into a single file per sample
## check that they're the same number of lines as the input
logging.debug("Concatenating tumor data")
tdata=resultsdir+"/"+tfilename+".combined.txt"
with open(tdata, 'w') as outfile:
	tdatacount=0
	for i in range(maxchunks):
		for line in open(resultsdir+"/"+tfilename+"."+str(i+1)+".txt"):
			tdatacount+=1
			outfile.write(line)

logging.debug("Concatenating normal data")
ndata=resultsdir+"/"+nfilename+".combined.txt"
with open(ndata, 'w') as outfile:
	ndatacount=0
	for i in range(maxchunks):
		for line in open(resultsdir+"/"+nfilename+"."+str(i+1)+".txt"):
			ndatacount+=1
			outfile.write(line)

## Compress all output
#logging.debug("Compressing all output")
#subprocess.call("gzip *txt", shell=True)

if(nvcflines==tdatacount==ndatacount):
	logging.debug("VCF, tumor data and normal data all contain "+str(nvcflines)+" lines")
else:
	logging.debug("WARNING: mismatch between number of lines in VCF and/or tumor and/or normal data")
	logging.debug("VCF contains "+str(nvcflines)+" lines")
	logging.debug("Tumor data contains "+str(tdatacount)+" lines")
	logging.debug("Normal data contains "+str(ndatacount)+" lines")


## Instead of re-writing tests in Python, we perform them using an R script, which also makes plotting easier

#exit()

## This tells us where the code is stored:
codedirectory = os.path.dirname(os.path.realpath(__file__))

logging.debug("Launching R script to perform filtering")
subprocess.call("Rscript --vanilla "+codedirectory+"/R/filter_main.R "+codedirectory+" "+parameters+" "+resultsdir+" "+tdata+" "+ndata+" "+vcfpath, shell=True)
logging.debug("R script complete")

exit()

# 
# ##Alignability 
# #if alignabilitytrack is not "skip":
# #  alignability(copy.deepcopy(myvcf),alignabilitytrack,tfilename+".alignability.txt")
# #
# ##On-target (i.e. within exome) - COMPLETE
# ### Skips this test if argument set to "skip"
# #if bedfile is not "skip":
# #  ontarget(copy.deepcopy(myvcf),bedfile,tfilename+".ontarget.txt")
# #
