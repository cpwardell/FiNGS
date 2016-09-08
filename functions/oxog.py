## OxoG contamination calculator

from __future__ import print_function
import logging
import pysam

from functions.shared_functions import *

def oxog(myvcf,bampath,filename):

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

    FoxoG="NA"
    F2R1=float() # Orientation for FoxoG calc
    F1R2=float() # Orientation for FoxoG calc
    score="NA"

    for alignedread in samfile.fetch(variant.CHROM,variant.POS,variant.POS+1):
      if not alignedread.is_proper_pair:
        continue
      else:
        ## Which base in the read is at the position we want?  Use the
        ## "aligned_pairs" list of tuples to determine this
        offset = [item for item in alignedread.aligned_pairs if item[1] == variant.POS][0][0]
        if(offset!=None and varianttype=="SNV" and alignedread.seq[offset]==str(variant.ALT[0])):
          if(alignedread.is_read1 and alignedread.is_reverse): # 83/pPr1 F2R1
            F2R1+=1
          if(alignedread.is_read2 and alignedread.mate_is_reverse): # 163/pPR2 F2R1
            F2R1+=1
          if(alignedread.is_read1 and alignedread.mate_is_reverse): # 99/pPR1 F1R2
            F1R2+=1
          if(alignedread.is_read2 and alignedread.is_reverse): # 147/pPr1 F1R2
            F1R2+=1

        ## Calculation of FoxoG - only necessary if C>A|G>T SNV
        ## Equation is: ALT_F1R2/(ALT_F1R2 + ALT_F2R1) or ALT_F2R1/(ALT_F1R2 + ALT_F2R1)
        ## C>anything:  numerator is ALT_F2R1
        ## A>anything:  numerator is ALT_F2R1
        ## G>anything:  numerator is ALT_F1R2
        ## T>anything:  numerator is ALT_F1R2
        if((variant.REF=="C" and str(variant.ALT[0])=="A") or (variant.REF=="G" and str(variant.ALT[0])=="T")):
          ## If sum of F1R2 and F2R1 is zero, all reads have an indel in them, so it should be removed
          if((F1R2 + F2R1)!=0):
            if(variant.REF=="C"):
              FoxoG = F2R1/(F1R2 + F2R1)
            if(variant.REF=="G"):
              FoxoG = F1R2/(F1R2 + F2R1)
          ## If FoxoG is still "NA" at this point, it must be rubbish, so set it to 1
          if(FoxoG=="NA"):
            FoxoG=1
          #score=-10+(100/3)*FoxoG
          score=FoxoG


    ## THIS IS WHERE WE WRITE OUTPUT
    #print((variant.REF+"\t"+str(variant.ALT[0])=="A")+str(score))
    #print(variant.REF+"\t"+str(variant.ALT[0])+"\t"+str(score))
    print(score)
    print(score,file=log)
  
  ## Close the bamfile
  samfile.close()

  ## Close the logfile
  log.close()


