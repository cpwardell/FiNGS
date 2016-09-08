#!/usr/bin/env python2.7

## Title: SNEF: Super NExtgen Filter
## Author: CP Wardell, UAMS
## Description: Takes a VCF file and BAM files used to generate it and generates various metrics that will be used to filter the VCF


import os
import vcf
import pysam
import logging
import copy
import argparse

from functions.shared_functions import *
from functions.mapping_quality import *
from functions.depth import *
from functions.base_quality import *
from functions.vaf import *
from functions.alignability import *
from functions.ontarget import *
from functions.oxog import *
from functions.directioner import *


## Gather command line args
parser = argparse.ArgumentParser()
parser.add_argument("-v", type=str, help="full path to vcf file",required=True)
parser.add_argument("-t", type=str, help="full path to tumor bam",required=True)
parser.add_argument("-n", type=str, help="full path to normal bam",required=True)
parser.add_argument("-a", type=str, help="full path to alignability track",required=True)
parser.add_argument("-b", type=str, help="full path to exome bed file",required=True)
parser.add_argument("--debug", help="Turn debugging statements on",action="store_true")
args = parser.parse_args()

## Turn logging on if desired
if(args.debug):
  logging.basicConfig(level=logging.DEBUG)
  logging.debug("Debugging mode enabled")

## Assign command line args to more friendly variable names
vcfpath=args.v
tbampath=args.t
nbampath=args.n
alignabilitytrack=args.a
bedfile=args.b

print(vcfpath)
print(tbampath)
print(nbampath)
print(alignabilitytrack)
print(bedfile)

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

## Read in VCF...
myvcf=vcflist(vcfpath)
logging.debug("VCF ingested, contains "+str(len(myvcf))+" records")

## Clean up VCF by editing CHROM names.  We want "1", not "chr1"
myvcf=cleanchroms(myvcf)

## Define filename prefixes for output..
tbasename = os.path.basename(tbampath)
tfilename = os.path.splitext(tbasename)[0]
nbasename = os.path.basename(nbampath) 
nfilename = os.path.splitext(nbasename)[0]


################################################################
## Begin tests #################################################
################################################################

#Mapping score - COMPLETE
#mapping_quality(copy.deepcopy(myvcf),tbampath,tfilename+".mapping_quality.txt")

#Depth - COMPLETE
#depth(copy.deepcopy(myvcf),tbampath,tfilename+".depth.txt")
#depth(copy.deepcopy(myvcf),nbampath,nfilename+".depth.txt")

#Base quality score - COMPLETE
#base_quality(copy.deepcopy(myvcf),tbampath,tfilename+".base_quality.txt")

#VAF (tumor and normal) - COMPLETE
#vaf(copy.deepcopy(myvcf),tbampath,tfilename+".vaf.txt")
#vaf(copy.deepcopy(myvcf),nbampath,nfilename+".vaf.txt")

#Alignability - COMPLETE
#alignability(copy.deepcopy(myvcf),alignabilitytrack,tfilename+".alignability.txt")

#On-target (i.e. within exome) - COMPLETE
#ontarget(copy.deepcopy(myvcf),bedfile,tfilename+".ontarget.txt")

#OxoG artefacts - COMPLETE
#oxog(copy.deepcopy(myvcf),tbampath,tfilename+".oxog.txt")

#Read direction bias - errors
directioner(copy.deepcopy(myvcf),tbampath,tfilename+".read_direction.txt")

## Other filters to implement...
#Homopolymers ???
#Local quality/read metrics ???




