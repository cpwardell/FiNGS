## On target test

from __future__ import print_function
import logging
import pysam

from functions.shared_functions import *

def ontarget(myvcf,bedfile,filename):

  log=open(filename,"w")

  ## Read in alignability track
  bedtab=pysam.Tabixfile(bedfile)

  for idx,variant in enumerate(myvcf):

    ## Remove 1 to POS to convert from 1-based to 0-based
    ## VCFs are 1 based, SAM coordinates are 0 based
    variant.POS=variant.POS-1

    ## This tells us if it is present in the bed file (1=present, 0=absent)
    try: 
      print(len(list(bedtab.fetch(variant.CHROM, variant.POS, variant.POS+1))),file=log)
      print(len(list(bedtab.fetch(variant.CHROM, variant.POS, variant.POS+1))))
    except:
      print(0,file=log)
      print(0)


  ## Close the logfile
  log.close()

