# Shared functions
import vcf
import gzip
import logging
import numpy

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
#def cleanchroms(myvcf):
#    for idx,variant in enumerate(myvcf):
#        variant.CHROM=cleanchrom(variant.CHROM)
#    return(myvcf)

## Function to clean a single chromosome name
#def cleanchrom(CHROM):
#    ## Remove "chr" string from chrom
#    if not CHROM.find("chr"):
#        CHROM=CHROM[3:]
#    if CHROM=="M":
#        CHROM="MT"
#    return(CHROM)

## Add mutation type 
def sixtypes(ref,alt):
    if((ref=="T" and alt=="G") or (ref=="A" and alt=="C")):
        mtype="T>G/A>C"
    if((ref=="T" and alt=="C") or (ref=="A" and alt=="G")):
        mtype="T>C/A>G"
    if((ref=="T" and alt=="A") or (ref=="A" and alt=="T")):
        mtype="T>A/A>T"
    if((ref=="C" and alt=="T") or (ref=="G" and alt=="A")):
        mtype="C>T/G>A"
    if((ref=="C" and alt=="G") or (ref=="G" and alt=="C")):
        mtype="C>G/G>C"
    if((ref=="C" and alt=="A") or (ref=="G" and alt=="T")):
        mtype="C>A/G>T"
    if((ref not in ["A","C","G","T"]) or (alt not in ["A","C","G","T"])):
        mtype="INDEL"
    return(mtype)
    
## Calculates quantiles for lists containing non-numeric values
## Returns 0 for for empty vectors
def quantilewithnas(listwithnas,quantile):
    numerics=[]
    result=0
    for element in listwithnas:
        try:
            numerics.append(float(element))
        except:
            continue
        try:
            result=numpy.quantile(numerics,1-quantile)
        except:
            continue
    return(result)

