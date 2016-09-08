## Alignability test
## Note that alignability track is created like so:
## bgzip wgEncodeCrgMapabilityAlign100mer.bedGraph ## compress using bgzip (from Samtools)
## tabix -p bed wgEncodeCrgMapabilityAlign100mer.bedGraph.gz ## index using tabix

from __future__ import print_function
import logging
import pysam

from functions.shared_functions import *

def alignability(myvcf,alignabilitytrack,filename):

  log=open(filename,"w")

  ## Read in alignability track
  alignfile=pysam.Tabixfile(alignabilitytrack)

  for idx,variant in enumerate(myvcf):

    ## Alignability track has "chr" contigs, so we need to add "chr" if
    ## not present in vcf records...
    if variant.CHROM.find("chr"):
      variant.CHROM="chr"+str(variant.CHROM)

    ## Remove 1 to POS to convert from 1-based to 0-based
    ## VCFs are 1 based, SAM coordinates are 0 based
    variant.POS=variant.POS-1

    #print(variant.CHROM+"\t"+str(variant.POS))

    ## Exceptions to catch non-standard CHROMs
    try:
      ## This statement catches variants that are not in the alignability file
      for record in alignfile.fetch(variant.CHROM, variant.POS, variant.POS+1):
        alignability=float(record.split("\t")[3])
        ## THIS IS WHERE WE WRITE OUTPUT
        print(alignability)
        print(alignability,file=log)

      if(len(list(alignfile.fetch(variant.CHROM, variant.POS, variant.POS+1)))==0):
        print(0)
        print(0,file=log)

    except:
      print(0)
      print(0,file=log)

  ## Close the logfile
  log.close()


