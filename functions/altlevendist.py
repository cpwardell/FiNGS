## VAF (variant allele frequency) calculator

## Note that for indels, it will ONLY count reads that support the ENTIRE indel

from __future__ import print_function
import logging
import pysam
import editdistance # to calculate distance between strings

from functions.shared_functions import *

## Define a median function:
def median(mylist):
  sorts = sorted(mylist)
  length = len(sorts)
  if not length % 2:
    return (sorts[length / 2] + sorts[length / 2 - 1]) / 2.0
  return sorts[length / 2]


def altlevendist(myvcf,bampath,filename):

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

    ## Variable to hold Levenschtein distance
    levdists=[]
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
              levdists.append(editdistance.eval(alignedread.get_reference_sequence(),alignedread.query_alignment_sequence))
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
                levdists.append(editdistance.eval(alignedread.get_reference_sequence(),alignedread.query_alignment_sequence))
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
              levdists.append(editdistance.eval(alignedread.get_reference_sequence(),alignedread.query_alignment_sequence))
          except:
            pass
    ## Calculate median Levenschtein distance for this variant
    try:
      #ld=float(sum(levdists))/len(levdists)
      ld=median(levdists)
    except:
      ld=0

    
    ## THIS IS WHERE WE WRITE OUTPUT
    variant.POS=variant.POS+1 # return variant.POS to original 1-based value
    print(str(variant.CHROM)+":"+str(variant.POS)+"\t"+variant.REF+"\t"+str(variant.ALT[0])+"\t"+str(refcount)+"\t"+str(altcount)+"\t"+varianttype+"\t"+str(ld))
    print(str(variant.CHROM)+":"+str(variant.POS)+"\t"+variant.REF+"\t"+str(variant.ALT[0])+"\t"+str(refcount)+"\t"+str(altcount)+"\t"+varianttype+"\t"+str(ld),file=log)
  
  ## Close the bamfile
  samfile.close()

  ## Close the logfile
  log.close()

