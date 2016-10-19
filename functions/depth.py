## Depth quality calculator

from __future__ import print_function
import logging
import pysam

from functions.shared_functions import *

def depth(myvcf,bampath,filename):

  log=open(filename,"w")

  ## Open the bamfile.  We do it ONCE and ONCE ONLY for speed reasons
  logging.debug("Opening bam file")
  samfile=pysam.Samfile(bampath,"rb") # rb = "read bam"

  for idx,variant in enumerate(myvcf):
    depth=0

    ## Remove 1 to POS to convert from 1-based to 0-based
    ## VCFs are 1 based, SAM coordinates are 0 based
    variant.POS=variant.POS-1

    for alignedread in samfile.fetch(variant.CHROM,variant.POS,variant.POS+1):
      if not alignedread.is_proper_pair:
        continue
      else:
        depth=depth+1

    ## THIS IS WHERE WE WRITE OUTPUT
    variant.POS=variant.POS+1 # return variant.POS to original 1-based value
    print(str(variant.CHROM)+":"+str(variant.POS)+"\t"+str(depth),file=log)
    print(str(variant.CHROM)+":"+str(variant.POS)+"\t"+str(depth))
  
  ## Close the bamfile
  samfile.close()
  ## Close the logfile
  log.close()

