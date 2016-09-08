## VAF (variant allele frequency) calculator

## Note that for indels, it will ONLY count reads that support the ENTIRE indel

from __future__ import print_function
import logging
import pysam

from functions.shared_functions import *

def vaf(myvcf,bampath,filename):

  log=open(filename,"w")

  ## Open the bamfile.  We do it ONCE and ONCE ONLY for speed reasons
  logging.debug("Opening bam file")
  samfile=pysam.Samfile(bampath,"rb") # rb = "read bam"

  for idx,variant in enumerate(myvcf):
    ## Determine if a variant is a SNV, insert or deletion
    if(len(variant.REF)==1 and len(str(variant.ALT[0]))==1):
      varianttype="SNV"
    elif(len(variant.REF)!=1):
      varianttype="del"
    elif(len(str(variant.ALT[0]))!=1):
      varianttype="ins"

    ## Exception for Strelka where some ALT alleles are encoded as "."
    if varianttype=="ins" and variant.ALT[0] is None:
      varianttype="SNV"

    ## Variables to count the number of REF and ALT reads
    refcount=0
    altcount=0

    ## Remove 1 to POS to convert from 1-based to 0-based
    ## VCFs are 1 based, SAM coordinates are 0 based
    variant.POS=variant.POS-1

    for alignedread in samfile.fetch(variant.CHROM,variant.POS,variant.POS+1):
      if not alignedread.is_proper_pair:
        continue
      else:
        ## Which base in the read is at the position we want?  Use the
        ## "aligned_pairs" list of tuples to determine this
        offset = [item for item in alignedread.aligned_pairs if item[1] == variant.POS][0][0]
        if(offset!=None):
          ## Counters for SNVs
          if(varianttype=="SNV"): 
            if(alignedread.seq[offset] == variant.REF[0]):
              refcount=refcount+1
            if(alignedread.seq[offset] == str(variant.ALT[0])):
              altcount=altcount+1
          ## Counters for insertions
          if(varianttype=="ins"):
            if(alignedread.seq[offset] == variant.REF[0] and
               alignedread.cigarstring==str(alignedread.rlen)+"M"):
              refcount=refcount+1
            ## We want reads that contain indels (i.e. the cigar string isn't readlengthM e.g. 74M)
            try:
              if(alignedread.seq[offset:offset+len(str(variant.ALT[0]))] == str(variant.ALT[0]) and
                 alignedread.cigarstring!=str(alignedread.rlen)+"M"):
                altcount=altcount+1
            except:
              pass
          ## Counters for deletions 
          if(varianttype=="del"):
          ## We want reads that DO NOT contain indels (i.e. the cigar string is readlengthM e.g. 74M)
            if(alignedread.seq[offset:offset+len(variant.REF[0])] == variant.REF[0] and
               alignedread.cigarstring==str(alignedread.rlen)+"M"):
              refcount=refcount+1
          try:
            if(alignedread.seq[offset] == str(variant.ALT[0]) and
               alignedread.cigarstring!=str(alignedread.rlen)+"M"):
              altcount=altcount+1
          except:
            pass
    ## Calculate VAF
    try:
      vaf=float(altcount)/(refcount+altcount)
    except:
      vaf=0
    
    ## THIS IS WHERE WE WRITE OUTPUT
    #print(vaf)
    print(variant.REF+"\t"+str(variant.ALT[0])+"\t"+str(refcount)+"\t"+str(altcount)+"\t"+varianttype+"\t"+str(vaf))
    print(variant.REF+"\t"+str(variant.ALT[0])+"\t"+str(refcount)+"\t"+str(altcount)+"\t"+varianttype+"\t"+str(vaf),file=log)
  
  ## Close the bamfile
  samfile.close()

  ## Close the logfile
  log.close()

