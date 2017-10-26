# Shared functions

import vcf

## Function to read VCF files and return an object
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


   
