## Mapping quality calculator
## Also produces proportion of reads that are mapping quality zero

import logging
import pysam
import statistics
from functions.shared_functions import *

def mapping_quality(myvcf,bampath,filename):

  log=open(filename,"w")

  ## Open the bamfile.  We do it ONCE and ONCE ONLY for speed reasons
  logging.debug("Opening bam file")
  samfile=pysam.Samfile(bampath,"rb") # rb = "read bam"

  for idx,variant in enumerate(myvcf):
     
    mapq=[]

    ## Remove 1 to POS to convert from 1-based to 0-based
    ## VCFs are 1 based, SAM coordinates are 0 based
    variant.POS=variant.POS-1

    for alignedread in samfile.fetch(variant.CHROM,variant.POS,variant.POS+1):
      if alignedread.is_duplicate:
        continue
      else:
        mapq.append(alignedread.mapq)

    ## Count how many MAPQ==0 reads there are
    zeroes=0
    for i in range(0,len(mapq),1):
        if (mapq[i] == 0):
            zeroes += 1

    ## THIS IS WHERE WE WRITE OUTPUT
    variant.POS=variant.POS+1 # return variant.POS to original 1-based value
    resultstring=str(variant.CHROM)+":"+str(variant.POS)+"\t"+str(statistics.median(mapq))+"\t"+str(zeroes/len(mapq))
    print(resultstring)
    print(resultstring,file=log)
  
  ## Close the bamfile
  samfile.close()

  ## Close the logfile
  log.close()

