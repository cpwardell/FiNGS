#!/usr/bin/env python3

# Copyright 2018 CP Wardell
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

## Title: FiNGS: Filters for Next Generation Sequencing
## Author: CP Wardell
## Description: Takes a VCF file from any variant caller
## and the aligned tumor/normal BAM files of Illumina paired-end data.
## Multiple metrics are calculated and used to filter the VCF.
## Designed and tested for SNVs; may work on indels.

#### TO DO LIST #########################################

## Do not gzip the pre-R output; read.table sometimes crashes in R

## Tumor-only (single BAM) mode
## Enable alignability filter
## Enable ontarget filter
## Enable foxog filter
## Complete Docker integration (build file)

#### END OF TO DO LIST #########################################

import os
import sys
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
## Create a new argparse class that will print the help message by default
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)
parser = MyParser(description="FiNGS: Filters for Next Generation Sequencing",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-v", type=str, help="absolute path to VCF file",required=True)
parser.add_argument("-t", type=str, help="absolute path to tumor BAM",required=True)
parser.add_argument("-n", type=str, help="absolute path to normal BAM",required=True)
#parser.add_argument("-a", type=str, help="absolute path to alignability track",required=False)
#parser.add_argument("-b", type=str, help="absolute path to BED file",required=False)
parser.add_argument("-d", type=str, help="absolute path to output directory",required=False,default="results")
parser.add_argument("-p", type=str, help="absolute path to filtering parameters (default is FiNGS/R/filter_parameters.txt",required=False,default="default")
parser.add_argument("-c", type=int, help="number of records to process per chunk",required=False,default=100)
parser.add_argument("-m", type=int, help="maximum read depth to process",required=False,default=1000)
parser.add_argument("-j", type=int, help="number of processors to use (default is -1, use all available resources)",required=False,default=-1)
parser.add_argument("--logging", help="Set logging level (default is INFO, can be DEBUG for more detail or NOTSET for silent)",required=False,default="INFO")
parser.add_argument("--donotcleanup", help="Keep intermediate files (not recommended, will keep intermediate chunks)",required=False,default=False,action='store_true')
parser.add_argument("--overwrite", help="Overwrite previous results if they exist?",required=False,default=False,action='store_true')
parser.add_argument("--PASSonlyin", help="Only use variants with that the original caller PASSed?",required=False,default=False,action='store_true')
parser.add_argument("--PASSonlyout", help="Only write PASS variants to the output VCF",required=False,default=False,action='store_true')
args = parser.parse_args()

## Turn arguments into a nice string for printing
printargs=str(sys.argv)
printargs=printargs.replace(",","")
printargs=printargs.replace("'","")
printargs=printargs.replace("[","")
printargs=printargs.replace("]","")

## Assign command line args to more friendly variable names
vcfpath=args.v
tbampath=args.t
nbampath=args.n
#alignabilitytrack=args.a
#bedfile=args.b
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
logging.info("FiNGS was invoked using this command: "+printargs)
logging.info("VCF is: "+str(vcfpath))
logging.info("Tumor BAM is: "+str(tbampath))
logging.info("Normal BAM is: "+str(nbampath))
logging.info("Output directory is: "+str(resultsdir))
logging.info("Filtering parameter file is: "+str(parameters))
#logging.info("Alignability file is: "+str(alignabilitytrack))
#logging.info("Bedfile is: "+str(bedfile))
logging.info("Number of records per chunk is: "+str(chunksize))
logging.info("Maximum read depth is: "+str(maxdepth))
logging.info("Processor threads used is: "+str(njobs))
logging.info("Logging level is: "+str(args.logging))
logging.info("Keep intermediate files?: "+str(args.donotcleanup))
logging.info("Overwrite existing output?: "+str(args.overwrite))
logging.info("Process only caller-PASSed variants?: "+str(args.PASSonlyin))
logging.info("Output only FiNGS-PASSed variants?: "+str(args.PASSonlyout))

## Define filename prefixes for output..
tbasename = os.path.basename(tbampath)
#tfilename = os.path.splitext(tbasename)[0]
tfilename = "tumor"
#nbasename = os.path.basename(nbampath) 
nbasename = "normal"
nfilename = os.path.splitext(nbasename)[0]

## Check for existing data; if they exist, skip straight to the R filtering phase
tdata=resultsdir+"/"+tfilename+".combined.txt"
ndata=resultsdir+"/"+nfilename+".combined.txt"

if(args.overwrite or (not os.path.exists(ndata) and not os.path.exists(tdata))):

	## Find size of VCF; vital for chunking and estimating runtime
	logging.info("VCF line counting begun")
	nvcflines=vcfcount(vcfpath)
	maxchunks=math.ceil(nvcflines/chunksize)
	logging.info("VCF contains "+str(nvcflines)+" records")

	## Multithreaded code:
	## VCF must be opened once per bamfile

	## Tumor bam
	vcf_reader = vcf.Reader(open(vcfpath, 'r'))
	Parallel(n_jobs=njobs, backend="threading")(delayed(primary)(vcfchunk,tbampath, resultsdir+"/"+tfilename+"."+str(chunknumber+1)+".txt",chunknumber,maxchunks,maxdepth,args.PASSonlyin) for chunknumber,vcfchunk in enumerate(vcfyield(vcf_reader,chunksize,nvcflines)))

	## Normal bam
	vcf_reader = vcf.Reader(open(vcfpath, 'r'))
	Parallel(n_jobs=njobs, backend="threading")(delayed(primary)(vcfchunk,nbampath, resultsdir+"/"+nfilename+"."+str(chunknumber+1)+".txt",chunknumber,maxchunks,maxdepth,args.PASSonlyin) for chunknumber,vcfchunk in enumerate(vcfyield(vcf_reader,chunksize,nvcflines)))
	## End of multithreaded code:

	## Combine all chunks into a single file per sample
	## check that they're the same number of lines as the input
	logging.debug("Concatenating tumor data")
	with open(tdata, 'w') as outfile:
		tdatacount=0
		for i in range(maxchunks):
			for line in open(resultsdir+"/"+tfilename+"."+str(i+1)+".txt"):
				tdatacount+=1
				outfile.write(line)

	logging.debug("Concatenating normal data")
	with open(ndata, 'w') as outfile:
		ndatacount=0
		for i in range(maxchunks):
			for line in open(resultsdir+"/"+nfilename+"."+str(i+1)+".txt"):
				ndatacount+=1
				outfile.write(line)

	## Delete intermediate chunks if donotcleanup is true
	if(not args.donotcleanup):
		logging.debug("Deleting intermediate chunks")
		for i in range(maxchunks):
			os.remove(resultsdir+"/"+nfilename+"."+str(i+1)+".txt")
			os.remove(resultsdir+"/"+tfilename+"."+str(i+1)+".txt")

	## Compress all output
	#logging.debug("Compressing all output")
	#subprocess.call("gzip *txt", shell=True)

	if(nvcflines==tdatacount==ndatacount):
		logging.debug("VCF, tumor data and normal data all contain "+str(nvcflines)+" lines")
	else:
		logging.debug("Potentional mismatch between number of lines in VCF and/or tumor and/or normal data; is the --PASSonlyin option being used?")
		logging.debug("VCF contains "+str(nvcflines)+" lines")
		logging.debug("Tumor data contains "+str(tdatacount)+" lines")
		logging.debug("Normal data contains "+str(ndatacount)+" lines")

else:
	logging.info("Previous results found; skipping ahead to filtering steps")

## Instead of re-writing tests in Python, we perform them using an R script, which also makes plotting easier

## This tells us where the code is stored:
codedirectory = os.path.dirname(os.path.realpath(__file__))

logging.info("Launching R script to perform filtering")
subprocess.call("Rscript --vanilla "+codedirectory+"/R/filter_main.R "+codedirectory+" "+parameters+" "+resultsdir+" "+tdata+" "+ndata+" "+vcfpath+" "+str(args.PASSonlyin)+" "+str(args.PASSonlyout), shell=True)
logging.info("R script complete")

exit()

