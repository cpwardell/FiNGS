#!/usr/bin/env python3

## Title: FiNGS: Filters for Next Generation Sequencing
## Author: CP Wardell
## Description: Takes a VCF file from any variant caller
## and the tumor/normal BAM files of Illumina paired-end data.
## Multiple metrics are calculated and used to filter the VCF.
## Designed and tested for SNVs; may work on indels.

#### TO DO LIST #########################################

## Clean-up option; delete intermediate "chunks"
## Can we gzip the pre-R output? - do not do this; read.table sometimes crashes
## R should write/append to the log file
## R should output either a complete VCF or just PASSed results
## Option to only use "PASS" results from variant caller "trustvc" option
## Test for existing output and run R stage only
## Some way to pass parameter file in 
## Docker build file
## Remove number of normal contigs test?
## Write alignability filter
## Write ontarget filter
## Write foxog filter
## Delete deprecated files from functions dir
## Write a decent readme.md file for github

#CRITICAL	50
#ERROR		40
#WARNING	30
#INFO		20 # this should be standard turned on level, fairly minimal
#DEBUG		10 # this should be detailed messages
#NOTSET		0 # this should be optional default level (silent)

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

#### END OF TO DO LIST #########################################

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
from functions.primary import *
from joblib import Parallel, delayed

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
parser.add_argument("--logging", help="Set logging level (default is INFO, can be DEBUG for more detail or NOTSET for silent)",required=False,default="INFO")
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
## Trycatch prevents exception if directory alr exists or location is unwritable 
try:
	os.mkdir(resultsdir)
except Exception as e:
	#logging.debug("WARNING: error creating results directory: "+str(e))
	pass

## Set logging to desired level
if(args.logging=="INFO"):
  logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s',datefmt='%Y-%m-%d %H:%M:%S',filename=resultsdir+'/log.txt',filemode='w')
  console = logging.StreamHandler()
  console.setLevel(logging.INFO)
  formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
  console.setFormatter(formatter)
  logging.getLogger('').addHandler(console)

if(args.logging=="DEBUG"):
  logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s',datefmt='%Y-%m-%d %H:%M:%S',filename=resultsdir+'/log.txt',filemode='w')
  console = logging.StreamHandler()
  console.setLevel(logging.DEBUG)
  formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
  console.setFormatter(formatter)
  logging.getLogger('').addHandler(console)
  logging.debug("Debugging mode enabled")

## Write the command line parameters to the log file
logging.info("VCF is: "+str(vcfpath))
logging.info("Tumor bamfile is: "+str(tbampath))
logging.info("Normal bamfile is: "+str(nbampath))
logging.info("Alignability file is: "+str(alignabilitytrack))
logging.info("Bedfile is: "+str(bedfile))
logging.info("Chunksize is: "+str(chunksize)+" records per chunk")
logging.info("Maximum read depth is: "+str(maxdepth))

## Define filename prefixes for output..
tbasename = os.path.basename(tbampath)
#tfilename = os.path.splitext(tbasename)[0]
tfilename = "tumor"
#nbasename = os.path.basename(nbampath) 
nbasename = "normal"
nfilename = os.path.splitext(nbasename)[0]

## Find size of VCF; vital for chunking and estimating runtime
logging.info("VCF line counting begun")
nvcflines=vcfcount(vcfpath)
maxchunks=math.ceil(nvcflines/chunksize)
logging.info("VCF contains "+str(nvcflines)+" records")

## Multithreaded code:
## VCF must be opened once per bamfile

## Tumor bam
vcf_reader = vcf.Reader(open(vcfpath, 'r'))
Parallel(n_jobs=njobs, backend="threading")(delayed(primary)(vcfchunk,tbampath, resultsdir+"/"+tfilename+"."+str(chunknumber+1)+".txt",chunknumber,maxchunks,maxdepth) for chunknumber,vcfchunk in enumerate(vcfyield(vcf_reader,chunksize,nvcflines)))

## Normal bam
vcf_reader = vcf.Reader(open(vcfpath, 'r'))
Parallel(n_jobs=njobs, backend="threading")(delayed(primary)(vcfchunk,nbampath, resultsdir+"/"+nfilename+"."+str(chunknumber+1)+".txt",chunknumber,maxchunks,maxdepth) for chunknumber,vcfchunk in enumerate(vcfyield(vcf_reader,chunksize,nvcflines)))
## End of multithreaded code:

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
	logging.info("VCF, tumor data and normal data all contain "+str(nvcflines)+" lines")
else:
	logging.warning("WARNING: mismatch between number of lines in VCF and/or tumor and/or normal data")
	logging.debug("VCF contains "+str(nvcflines)+" lines")
	logging.debug("Tumor data contains "+str(tdatacount)+" lines")
	logging.debug("Normal data contains "+str(ndatacount)+" lines")


## Instead of re-writing tests in Python, we perform them using an R script, which also makes plotting easier

## This tells us where the code is stored:
codedirectory = os.path.dirname(os.path.realpath(__file__))

logging.info("Launching R script to perform filtering")
subprocess.call("Rscript --vanilla "+codedirectory+"/R/filter_main.R "+codedirectory+" "+parameters+" "+resultsdir+" "+tdata+" "+ndata+" "+vcfpath, shell=True)
logging.info("R script complete")

exit()

