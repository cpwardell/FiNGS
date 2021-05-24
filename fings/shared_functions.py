# Shared functions
import gzip
import logging
import numpy


## Function to count the number of records in a VCF file
def vcfcount(vcfpath):
    ## Check ending of vcfpath; if "vcf", ok, if "gz", open file
    ## Otherwise, quit with an error
    if(vcfpath.endswith("vcf")):
        vcf = open(vcfpath, 'r')
        comment="#"  # character version
    if(vcfpath.endswith("gz")):
        vcf = gzip.open(vcfpath, 'r')
        comment=b'#'  # byte version
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


## Calculates quantiles for lists containing non-numeric values
## Returns 0 for for empty vectors
def quantilewithnas(listwithnas,quantile):
    numerics=[]
    for element in listwithnas:
        try:
            numerics.append(float(element))
        except:
            continue
    try:
        result=numpy.quantile(numerics,1-quantile)
    except:
        result=0
    return(result)

## Writes header for tumor and normal filter metrics
def writeheader(file):
    cnames=["UID","CHR","POS","REF","ALT","refcount","altcount","varianttype","depth","vaf","raf","oaf",
            "medianbaseq","medianbaseqref","medianbaseqalt","medianmapq","medianmapqref","medianmapqalt",
            "zeros","zerospersite","softreadlengthsrefmean","softreadlengthsaltmean","goodoffsetproportion",
            "distancetoend1median","mad1","distancetoend2median","mad2","distancetoend1medianref","madref1",
            "distancetoend2medianref","madref2","distancetoend1medianalt","madalt1","distancetoend2medianalt",
            "madalt2","shortestdistancetoendmedian","madaltshort","sb","gsb","fishp",
            "FR","FA","RR","RA","altsb","refsb","allsb","F1R2","F2R1","FoxoG","refld","altld",
            "refsecondprop","altsecondprop","refbadorientationprop","altbadorientationprop","refmatecontigcount",
            "altmatecontigcount","sixtypes"]
    file.write('\t'.join(cnames[0:]) + '\n')
