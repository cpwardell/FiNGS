# Shared functions

import vcf
import gzip
import os

## Function to count the number of records in a VCF file
def vcfcount(vcfpath):

  ## Check ending of vcfpath; if "vcf", ok, if "gz", open file
  ## Otherwise, quit with an error
  if(vcfpath.endswith("vcf")):
    vcf = open(vcfpath, 'r')
    comment="#" # character version
  if(vcfpath.endswith("gz")):
    vcf = gzip.open(vcfpath, 'r')
    comment=b'#' # byte version
  if not (vcfpath.endswith(("vcf","gz"))):
    logging.debug("ERROR: VCF doesn't end with \"vcf\" or \"gz\":"+vcfpath)
    exit("ERROR: VCF doesn't end with \"vcf\" or \"gz\":"+vcfpath)

  recordcount=0
  for line in vcf:
    if not line.startswith(comment):
      recordcount+=1
  return(recordcount)

def vcfyield(vcfbuffer,chunksize,finalrecord):
  recordcount=0
  vcfchunk=list()
  for record in vcfbuffer:
    vcfchunk.append(record)
    recordcount+=1
    if(len(vcfchunk)==chunksize or recordcount==finalrecord):
      yield(vcfchunk)
      vcfchunk=list()


## Function to read whole VCF files and return an object
def vcflist(vcfpath):
  myvcf=list()
  vcf_reader = vcf.Reader(open(vcfpath, 'r'))
  for record in vcf_reader:
    #print record
    myvcf.append(record)
  return(myvcf)
 
## Function for calculating the mean; we don't want to use numpy
def mean(numbers):
    try:
      x=float(sum(numbers))/len(numbers)
    except:
      x=0
    return(x)

## Function for removing "chr" from chromosome names
## Also rename "M" to "MT"
def cleanchroms(myvcf):
  for idx,variant in enumerate(myvcf):
    variant.CHROM=cleanchrom(variant.CHROM)
  return(myvcf)

## Function to clean a single chromosome name
def cleanchrom(CHROM):
  ## Remove "chr" string from chrom
  if not CHROM.find("chr"):
    CHROM=CHROM[3:]

  if CHROM=="M":
    CHROM="MT"

  return(CHROM)


   
