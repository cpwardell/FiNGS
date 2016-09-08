## Read direction counter
## August 2016; doesn't work on indels.  Several types of read will be skipped, 
## including non-paired reads and low quality reads

from __future__ import print_function
import logging
import pysam

from functions.shared_functions import *

def directioner(myvcf,bampath,filename):

  log=open(filename,"w")

  ## Open the bamfile.  We do it ONCE and ONCE ONLY for speed reasons
  logging.debug("Opening bam file")
  samfile=pysam.Samfile(bampath,"rb") # rb = "read bam"

  for idx,variant in enumerate(myvcf):
   
    ## Determine if a variant is a SNV, insert or deletion
    if(len(variant.REF)==1 and len(str(variant.ALT[0]))==1):
      varianttype="SNV"
    else:
      varianttype="indel"

    ## Remove 1 to POS to convert from 1-based to 0-based
    ## VCFs are 1 based, SAM coordinates are 0 based
    variant.POS=variant.POS-1

    AF=0 # ALT FORWARD
    AR=0 # ALT REVERSE
    RF=0 # REF FORWARD
    RR=0 # REF REVERSE

    for alignedread in samfile.fetch(variant.CHROM,variant.POS,variant.POS+1):
      if not alignedread.is_proper_pair:
        continue
      else:
        ## Which base in the read is at the position we want?  Use the
        ## "aligned_pairs" list of tuples to determine this
        offset = [item for item in alignedread.aligned_pairs if item[1] == variant.POS][0][0]
        if(offset!=None and varianttype=="SNV" and alignedread.seq[offset]==str(variant.ALT[0])):
          if(alignedread.is_reverse):
            AR+=1
          if(alignedread.mate_is_reverse):
            AF+=1
        if(offset!=None and varianttype=="SNV" and alignedread.seq[offset]==variant.REF):
          if(alignedread.is_reverse):
            RR+=1
          if(alignedread.mate_is_reverse):
            RF+=1

    ## Write output
    print(variant.CHROM+":"+str(variant.POS+1)+"\t"+variant.REF+"\t"+str(variant.ALT[0])+"\t"+varianttype+"\t"+str(AF)+"\t"+str(AR)+"\t"+(str(RF)+"\t"+str(RR)))
    print(str(AF)+"\t"+str(AR)+"\t"+(str(RF)+"\t"+str(RR)),file=log)
  
  ## Close the bamfile
  samfile.close()

  ## Close the logfile
  log.close()


