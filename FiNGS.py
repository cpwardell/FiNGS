#!/usr/bin/env python3

# Copyright 2018-2019 CP Wardell
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
## Designed and tested for SNVs; currently excludes indels

#### TO DO LIST #########################################

## Tumor-only (single BAM) mode
## Add indel support

#### END OF TO DO LIST #########################################
import os
import sys
import vcf
import logging
import argparse
import math
import gzip
import csv
import datetime

from functions.shared_functions import vcfcount, vcfyield, quantilewithnas
from functions.primary import primary
from functions.filter_functions import applyfilters
from joblib import Parallel, delayed


def main():

    ## Gather command line args
    ## Create a new argparse class that will print the help message by default
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)
    parser=MyParser(description="FiNGS: Filters for Next Generation Sequencing", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", type=str, help="absolute path to VCF file", required=True)
    parser.add_argument("-t", type=str, help="absolute path to tumor BAM", required=True)
    parser.add_argument("-n", type=str, help="absolute path to normal BAM", required=True)
    parser.add_argument("-r", type=str, help="absolute path to faidx indexed reference genome; required if using \'repeats\' filter", required=False, default="None")
    parser.add_argument("-d", type=str, help="absolute path to output directory", required=False,default="results")
    parser.add_argument("-p", type=str, help="absolute path to filtering parameters (default is filter_parameters.txt", required=False, default="filter_parameters.txt")
    parser.add_argument("-c", type=int, help="number of records to process per chunk",required=False,default=100)
    parser.add_argument("-m", type=int, help="maximum read depth to process",required=False,default=1000)
    parser.add_argument("-j", type=int, help="number of processors to use (default is -1, use all available resources)", required=False, default=-1)
    parser.add_argument("--ICGC", help="Use filters identical to those recommended by the ICGC (Alioto et al, 2015). Overrides \'-p\' flag",required=False, default=False,action='store_true')
    parser.add_argument("--logging", help="Set logging level (default is INFO, can be DEBUG for more detail or NOTSET for silent)", required=False, default="INFO")
    parser.add_argument("--overwrite", help="Overwrite previous results if they exist", required=False, default=False, action='store_true')
    parser.add_argument("--PASSonlyin", help="Only use variants with that the original caller PASSed?", required=False, default=False, action='store_true')
    parser.add_argument("--PASSonlyout", help="Only write PASS variants to the output VCF", required=False, default=False, action='store_true')
    args=parser.parse_args()

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
    referencegenome=args.r
    resultsdir=args.d
    parameters=args.p
    chunksize=args.c
    maxdepth=args.m
    njobs=args.j

    ## Use ICGC parameters if flag is set
    if(args.ICGC):
        parameters="icgc_filter_parameters.txt"

    ## Create directory to put results in
    ## Trycatch prevents exception if location is unwritable 
    try:
        os.makedirs(resultsdir,exist_ok=True)
    except Exception as e:
        print("CRITICAL ERROR: results directory could not be created: "+str(e))
        print("If using Docker, the current directory must be writeable by any user")
        sys.exit()

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
    starttime=datetime.datetime.now()
    logging.info("FiNGS was invoked using this command: "+printargs)
    logging.info("Start time is: "+str(starttime))
    logging.info("VCF is: "+str(vcfpath))
    logging.info("Tumor BAM is: "+str(tbampath))
    logging.info("Normal BAM is: "+str(nbampath))
    logging.info("Output directory is: "+str(resultsdir))
    logging.info("Filtering parameter file is: "+str(parameters))
    logging.info("Reference genome file is: "+str(referencegenome))
    logging.info("Number of records per chunk is: "+str(chunksize))
    logging.info("Maximum read depth is: "+str(maxdepth))
    logging.info("Processor threads used is: "+str(njobs))
    logging.info("Logging level is: "+str(args.logging))
    logging.info("Overwrite existing output?: "+str(args.overwrite))
    logging.info("Process only caller-PASSed variants?: "+str(args.PASSonlyin))
    logging.info("Output only FiNGS-PASSed variants?: "+str(args.PASSonlyout))

    
    ## If the parameter file is the default value or if using ICGC parameters, create an absolute path to it;
    if(parameters=="filter_parameters.txt" or args.ICGC):
        codedirectory = os.path.dirname(os.path.realpath(__file__))
        parameters = os.path.join(codedirectory,parameters)
    ## Read filtering parameters into a dictionary
    logging.info("Reading filtering parameters from this file: "+parameters)    
    pdict={}
    try:
        with open(parameters) as tsvin:
            tsvin = csv.reader(tsvin, delimiter='\t')
            for row in tsvin:
                pdict[row[0]]=float(row[1])
    except Exception as e:
        print("CRITICAL ERROR: unreadable filter parameter file: "+str(e))
        sys.exit()

    ## Check for existing data; if they exist, skip straight to the filtering phase
    tdata=resultsdir+"/tumor.combined.txt.gz"
    ndata=resultsdir+"/normal.combined.txt.gz"
    sdata=resultsdir+"/summarystats.txt.gz"

    if(args.overwrite or (not os.path.exists(ndata) and not os.path.exists(tdata))):
        ## Find size of VCF; vital for chunking and estimating runtime
        logging.info("VCF line counting begun")
        try:
            nvcflines=vcfcount(vcfpath)
        except Exception as e:
            print("CRITICAL ERROR: unreadable VCF file: "+str(e))
            sys.exit()
        maxchunks=math.ceil(nvcflines/chunksize)
        logging.info("VCF contains "+str(nvcflines)+" records")

        ## Multithreaded code:
        ## VCF must be opened once per bamfile

        ## Tumor bam
        vcf_reader = vcf.Reader(open(vcfpath, 'r'))
        tumvars=Parallel(n_jobs=njobs, backend="threading")(delayed(primary)(vcfchunk,tbampath,"tumor",chunknumber,maxchunks,maxdepth,args.PASSonlyin) for chunknumber,vcfchunk in enumerate(vcfyield(vcf_reader,chunksize,nvcflines)))

        ## Normal bam
        vcf_reader = vcf.Reader(open(vcfpath, 'r'))
        normvars=Parallel(n_jobs=njobs, backend="threading")(delayed(primary)(vcfchunk,nbampath,"normal",chunknumber,maxchunks,maxdepth,args.PASSonlyin) for chunknumber,vcfchunk in enumerate(vcfyield(vcf_reader,chunksize,nvcflines)))
        ## End of multithreaded code:

        ## Variables to store aggregate values for later filtering e.g. global averages
        strandbiastumor=list()

        ## Output all chunks into a single file per sample
        logging.info("Writing tumor data to temporary file: "+tdata)
        tout=gzip.open(tdata, 'wt')
        for LINE in tumvars:
            if(len(LINE) is not 0): # protects against no variants in chunk e.g. if none are PASS
                tcsv = csv.reader(LINE, delimiter='\t')
                strandbiastumor.append([item[37] for item in tcsv])  # hardcoded location of sb metric
                print(*LINE,sep="\n",file=tout)
        tout.close()

        logging.info("Writing normal data to temporary file: "+ndata)
        nout=gzip.open(ndata, 'wt')
        for LINE in normvars:
            if(len(LINE) is not 0): # protects against no variants in chunk e.g. if none are PASS
                print(*LINE,sep="\n",file=nout)
        nout.close()

        ## Convert list of lists into a single list
        strandbiastumor=[item for sublist in strandbiastumor for item in sublist]
        ## If strandbias cutoff proportion isn't a specified filter, set the value to 0 so it can still be written
        try:
            sbquantile=pdict["strandbiasprop"]
        except:
            sbquantile=0
        strandbiastumorq=quantilewithnas(strandbiastumor,sbquantile)
        ## Write any summary stats to a file
        logging.info("Writing summary data to temporary file: "+sdata)
        sout=gzip.open(sdata, 'wt')
        print("strandbiastumorq"+"\t"+str(strandbiastumorq),file=sout)
        sout.close()

    else:
        logging.info("Previous results found; skipping ahead to filtering steps")

    ## Filter data
    applyfilters(tdata,ndata,sdata,pdict,resultsdir,vcfpath,referencegenome,args.PASSonlyin,args.PASSonlyout)

    endtime=datetime.datetime.now()
    logging.info("End time is: "+str(endtime))
    logging.info("Total time taken: "+str(endtime-starttime))

    exit()


## Execute main method
if __name__ == '__main__':
    main()
