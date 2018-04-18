## Primary function.  There is one large function rather than many smaller ones
## to ensure the BAM file(s) are read once and once only
## Note that for indels, it will ONLY count reads that support the ENTIRE indel

import logging
import pysam
import statistics
import scipy.stats
import editdistance
import numpy
from functions.shared_functions import *

def primary(myvcf,bampath,filename,chunknumber,maxchunks,maxdepth,PASS):

  log=open(filename,"w")

  logging.info("Started processing chunk "+filename+" ("+str(chunknumber+1)+" of "+str(maxchunks)+") containing "+str(len(myvcf))+" records")

  ## Open the bamfile.  We do it ONCE and ONCE ONLY for speed reasons
  #logging.debug("Opening bam file")
  samfile=pysam.Samfile(bampath,"rb") # rb = "read bam"

  for idx,variant in enumerate(myvcf):
    ## Variables to hold depth, basequality, mapping quality
    depth=0
    baseq=[]
    baseqref=[]
    baseqalt=[]
    mapq=[]
    mapqref=[]
    mapqalt=[]

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

    ## Variables to count secondary alignments
    refsecond=0
    altsecond=0

    ## Variable for counting orientation
    refbadorientation=0
    altbadorientation=0

    ## Variable for counting mate contigs
    refmatecontigs=[]
    altmatecontigs=[]

    ## Variables to store direction of ALT/REF reads to calculate strand bias
    FR=0 # Forward strand, reference allele
    FA=0 # Forward strand, alternate allele
    RR=0 # Reverse strand, reference allele
    RA=0 # Reverse strand, alternate allele

    ## Variables for FoxoG calculations
    FoxoG="NA"
    F2R1=float() # Orientation for FoxoG calc
    F1R2=float() # Orientation for FoxoG calc

    ## Variable to store read lengths, variant offsets and distance to end of reads
    hardreadlengths=[] # actual length of read
    softreadlengths=[] # length of read after soft-clipping
    softreadlengthsref=[]
    softreadlengthsalt=[]
    offsets=[]
    distancetoend1=[]
    distancetoend2=[]
    distancetoend1ref=[]
    distancetoend2ref=[]
    distancetoend1alt=[]
    distancetoend2alt=[]

    ## Variables to hold Levenschtein distances
    reflevdists=[]
    altlevdists=[]

    ## If we are only using variants PASSed by the caller, skip the samfile pileup phase
    if(not PASS or (PASS and variant.FILTER==[])):
      
    ## Remember that VCFs are 1 based, SAM coordinates are 0 based ##
      for alignedread in samfile.fetch(variant.CHROM,variant.POS-1,variant.POS):
        if alignedread.is_duplicate:
          continue
        else:
          ## If depth is beyond maxdepth, break loop and calculate results
          ## Also, emit a debugging statement
          if(depth>=maxdepth):
            logging.info("WARNING: variant "+str(idx+1)+" of "+str(len(myvcf))+" in chunk "+str(chunknumber+1)+" of "+str(maxchunks)+" at position "+str(variant.CHROM)+":"+str(variant.POS)+" is deeper than maximum depth of "+str(maxdepth)+", calculations truncated at this depth")
            break

          ## Increment depth by 1
          depth+=1

          ## Append new mapping quality
          mapq.append(alignedread.mapq)

          ## Which base in the read is at the position we want?  Use the
          ## "aligned_pairs" list of tuples to determine this
          try:
            offset = [item for item in alignedread.aligned_pairs if item[1] == variant.POS-1][0][0]
          except:
            offset=None

          if(offset!=None):

            ## Append new base quality (we subtract 33 because SAM specification tells us to)
            baseq.append(ord(alignedread.qual[offset])-33)

            ## How long are the reads?  Store them in a list
            hardreadlengths.append(alignedread.infer_query_length()) # hard read lengths
            softreadlengths.append(alignedread.query_alignment_end) # soft-clipped read lengths
            offsets.append(offset)
            #distancetoend.append(abs(min(offset,offset-alignedread.infer_query_length()))) # hard read lengths
            de1=abs(offset-alignedread.query_alignment_start)
            de2=abs(offset-alignedread.query_alignment_end+1)
            distancetoend1.append(de1) # distance to lefthand soft-clipped read end
            distancetoend2.append(de2) # distance to righthand soft-clipped read end
            
            #print(min(offset,offset-alignedread.query_alignment_end))
            #print(str(alignedread.query_name)+"\t"+str(offset)+"\t"+str(abs(offset-alignedread.query_alignment_end))+"\t"+str(alignedread.query_alignment_start)+"\t"+str(alignedread.query_alignment_end-offset))

            #print(str(alignedread.query_name)+"\t"+str(offset)+"\t"+str(de1)+"\t"+str(alignedread.seq[offset])+"\t"+str(de2)+"\t"+str(alignedread.is_proper_pair)+"\t"+str(alignedread.is_secondary))
            
            #print(str(alignedread.query_name)+"\t"+str(alignedread.next_reference_id)+"\t"+str(alignedread.next_reference_id))
            #print(str(alignedread.reference_id)+"\t"+str(alignedread.next_reference_id))

            ## Counters for SNVs
            if(varianttype=="SNV"): 
              if(alignedread.seq[offset] == variant.REF[0]):
                refcount+=1
                if alignedread.is_secondary:
                  refsecond+=1
                if not alignedread.is_proper_pair:
                  refbadorientation+=1

                ## Append the contig identifier of the mate to an array
                refmatecontigs.append(alignedread.next_reference_id)

                ## Append new base quality (we subtract 33 because SAM specification tells us to)
                baseqref.append(ord(alignedread.qual[offset])-33)
                mapqref.append(alignedread.mapq)
                softreadlengthsref.append(alignedread.query_alignment_end) # soft-clipped read lengths
                #distancetoendref.append(abs(min(offset,offset-alignedread.infer_query_length()))) # hard read lengths
                distancetoend1ref.append(de1) # soft-clipped read length
                distancetoend2ref.append(de2) # soft-clipped read length
                try:
                  reflevdists.append(editdistance.eval(alignedread.get_reference_sequence(),alignedread.query_alignment_sequence))
                except:
                  reflevdists.append(1)
                if not (alignedread.is_reverse):
                  FR+=1
                else:
                  RR+=1
              if(alignedread.seq[offset] == str(variant.ALT[0])):
                altcount+=1
                if alignedread.is_secondary:
                  altsecond+=1
                if not alignedread.is_proper_pair:
                  altbadorientation+=1

                ## Append the contig identifier of the mate to an array
                altmatecontigs.append(alignedread.next_reference_id)

                baseqalt.append(ord(alignedread.qual[offset])-33)
                mapqalt.append(alignedread.mapq)
                softreadlengthsalt.append(alignedread.query_alignment_end) # soft-clipped read lengths
                #distancetoendalt.append(abs(min(offset,offset-alignedread.infer_query_length()))) # hard read lengths
                distancetoend1alt.append(de1) # soft-clipped read length
                distancetoend2alt.append(de2) # soft-clipped read length
                ## Calculate edit distance; get_reference_sequence() will fail without an MD tag in the BAM file, so if there's an exception set
                ## a default value of 1
                try:
                  altlevdists.append(editdistance.eval(alignedread.get_reference_sequence(),alignedread.query_alignment_sequence))
                except:
                  altlevdists.append(1)
                if not (alignedread.is_reverse):
                  FA+=1
                else:
                  RA+=1
                ## Counters for FoxoG calculations
                if(alignedread.is_read1 and alignedread.is_reverse): # 83/pPr1 F2R1
                  F2R1+=1
                if(alignedread.is_read2 and alignedread.mate_is_reverse): # 163/pPR2 F2R1
                  F2R1+=1
                if(alignedread.is_read1 and alignedread.mate_is_reverse): # 99/pPR1 F1R2
                  F1R2+=1
                if(alignedread.is_read2 and alignedread.is_reverse): # 147/pPr1 F1R2
                  F1R2+=1

            ## Counters for insertions
            if(varianttype=="ins"):
              if(alignedread.seq[offset] == variant.REF[0] and
                 alignedread.cigarstring==str(alignedread.rlen)+"M"):
                refcount+=1
                baseqref.append(ord(alignedread.qual[offset])-33)
                mapqref.append(alignedread.mapq)
                softreadlengthsref.append(alignedread.query_alignment_end) # soft-clipped read lengths
                #distancetoendref.append(abs(min(offset,offset-alignedread.infer_query_length()))) # hard read lengths
                distancetoend1ref.append(de1) # soft-clipped read length
                distancetoend2ref.append(de2) # soft-clipped read length
              ## We want reads that contain indels (i.e. the cigar string isn't readlengthM e.g. 74M)
              try:
                if(alignedread.seq[offset:offset+len(str(variant.ALT[0]))] == str(variant.ALT[0]) and
                   alignedread.cigarstring!=str(alignedread.rlen)+"M"):
                  altcount+=1
                  baseqalt.append(ord(alignedread.qual[offset])-33)
                  mapqalt.append(alignedread.mapq)
                  softreadlengthsalt.append(alignedread.query_alignment_end) # soft-clipped read lengths
                  #distancetoendalt.append(abs(min(offset,offset-alignedread.infer_query_length()))) # hard read lengths
                  distancetoend1alt.append(de1) # soft-clipped read length
                  distancetoend2alt.append(de2) # soft-clipped read length
              except:
                pass
            ## Counters for deletions 
            if(varianttype=="del"):
            ## We want reads that DO NOT contain indels (i.e. the cigar string is readlengthM e.g. 74M)
              if(alignedread.seq[offset:offset+len(variant.REF[0])] == variant.REF[0] and
                 alignedread.cigarstring==str(alignedread.rlen)+"M"):
                  refcount+=1
                  baseqref.append(ord(alignedread.qual[offset])-33)
                  mapqref.append(alignedread.mapq)
                  softreadlengthsref.append(alignedread.query_alignment_end) # soft-clipped read lengths
                  #distancetoendref.append(abs(min(offset,offset-alignedread.infer_query_length()))) # hard read lengths
                  distancetoend1ref.append(de1) # soft-clipped read length
                  distancetoend2ref.append(de2) # soft-clipped read length
              try:
                if(alignedread.seq[offset] == str(variant.ALT[0]) and
                   alignedread.cigarstring!=str(alignedread.rlen)+"M"):
                  altcount+=1
                  baseqalt.append(ord(alignedread.qual[offset])-33)
                  mapqalt.append(alignedread.mapq)
                  softreadlengthsalt.append(alignedread.query_alignment_end) # soft-clipped read lengths
                  #distancetoendalt.append(abs(min(offset,offset-alignedread.infer_query_length()))) # hard read lengths
                  distancetoend1alt.append(de1) # soft-clipped read length
                  distancetoend2alt.append(de2) # soft-clipped read length
              except:
                pass

      ## Calculation of FoxoG - only necessary if C>A|G>T SNV
      ## Equation is: ALT_F1R2/(ALT_F1R2 + ALT_F2R1) or ALT_F2R1/(ALT_F1R2 + ALT_F2R1)
      ## C>anything or A>anything:  numerator is ALT_F2R1
      ## G>anything or T>anything:  numerator is ALT_F1R2
      if((variant.REF=="C" and str(variant.ALT[0])=="A") or (variant.REF=="G" and str(variant.ALT[0])=="T")):
        ## If sum of F1R2 and F2R1 is zero, all reads have an indel in them, so it should be removed
        if((F1R2 + F2R1)!=0):
          if(variant.REF=="C"):
            FoxoG = F2R1/(F1R2 + F2R1)
          if(variant.REF=="G"):
            FoxoG = F1R2/(F1R2 + F2R1)
          FoxoG=round(FoxoG,3)
        ## If FoxoG is still "NA" at this point, it must be rubbish, so set it to 1
        #if(FoxoG=="NA"):
        #  FoxoG=1
        #score=-10+(100/3)*FoxoG
        
      ## Calculate VAF, RAF and OAF
      try:
        vaf=round(float(altcount)/(depth),3)
      except:
        vaf="NA"
      try:
        raf=round(float(refcount)/(depth),3)
      except:
        raf="NA"
      try:
        oaf=round(((depth)-float(altcount)-float(refcount))/(depth),3)
      except:
        oaf="NA"


      ## Calculate strand bias and fisher's exact test pvalue for it
      def strandbias(FR,FA,RR,RA):
        try:
          sb=abs( ( FA/(FR+FA) ) - ( RA/(RR+RA) ) ) /  ((FA+RA)/(FR+FA+RR+RA))
          sb=round(sb,3)
        except:
          sb="NA"
        return(sb)

      def gatkstrandbias(FR,FA,RR,RA):
        try:
          sba=(( FA/(FR+FA) ) * ( RR/(RR+RA) )) / ((FR+RR)/(FR+FA+RR+RA))
          sba=round(sba,3)
        except:
          sba="NA"
        try:
          sbb=(( RA/(RR+RA) ) * ( FR/(FR+FA) )) / ((FR+RR)/(FR+FA+RR+RA))
          sbb=round(sbb,3)
        except:
          sbb="NA"
        return(max(sba,sbb))

      sb=strandbias(FR,FA,RR,RA)
      gsb=gatkstrandbias(FR,FA,RR,RA)
      fishp=scipy.stats.fisher_exact([[FA,RA],[FR,RR]])[1]
      fishp="%.3g" % fishp

      ## Calculate read length and proportion of the variants occurring in the first 2/3 of the reads
      hardmeanlength=round(mean(hardreadlengths),3)
      softmeanlength=round(mean(softreadlengths),3)
      softreadlengthsrefmean=round(mean(softreadlengthsref),3) # soft-clipped read lengths
      softreadlengthsaltmean=round(mean(softreadlengthsalt),3) # soft-clipped read lengths
      ## Number of offsets that are less than 2/3rd of the mean read length
      ## This can fail if the offsets list is empty
      try:
        goodoffsets=len(list(filter(lambda x: x<(2/3*softmeanlength),offsets)))
      except:
        goodoffsets="NA"
      try:
        goodoffsetproportion=round(goodoffsets/len(offsets),3)
      except:
        goodoffsetproportion="NA"

      ## Median distance to alignment ends and MADs of these
      try:
        distancetoend1median=statistics.median(distancetoend1)
      except:
        distancetoend1median="NA"
      try:
        mad1=round(statistics.median([abs(x - distancetoend1median) for x in distancetoend1]),3)
      except:
        mad1="NA"
      try:
        distancetoend2median=statistics.median(distancetoend2)
      except:
        distancetoend2median="NA"
      try:
        mad2=round(statistics.median([abs(x - distancetoend2median) for x in distancetoend2]),3)
      except:
        mad2="NA"

      try:
        distancetoend1medianref=statistics.median(distancetoend1ref)
      except:
        distancetoend1medianref="NA"
      try:
        madref1=round(statistics.median([abs(x - distancetoend1medianref) for x in distancetoend1ref]),3)
      except:
        madref1="NA"
      try:
        distancetoend2medianref=statistics.median(distancetoend2ref)
      except:
        distancetoend2medianref="NA"
      try:
        madref2=round(statistics.median([abs(x - distancetoend2medianref) for x in distancetoend2ref]),3)
      except:
        madref2="NA"


      try:
        distancetoend1medianalt=statistics.median(distancetoend1alt)
      except:
        distancetoend1medianalt="NA"
      try:
        madalt1=round(statistics.median([abs(x - distancetoend1medianalt) for x in distancetoend1alt]),3)
      except:
        madalt1="NA"
      try:
        distancetoend2medianalt=statistics.median(distancetoend2alt)
      except:
        distancetoend2medianalt="NA"
      try:
        madalt2=round(statistics.median([abs(x - distancetoend2medianalt) for x in distancetoend2alt]),3)
      except:
        madalt2="NA"
      try:
        shortestdistancetoend=list(map(min,zip(distancetoend1alt,distancetoend2alt)))
      except:
        shortestdistancetoend="NA"
      try:
        shortestdistancetoendmedian=round(statistics.median(shortestdistancetoend),3)
      except:
        shortestdistancetoendmedian="NA"
      try:
        madaltshort=round(statistics.median([abs(x - shortestdistancetoendmedian) for x in shortestdistancetoend]),3)
      except:
        madaltshort="NA"


      ## Count how many MAPQ==0 reads there are and normalize for depth
      try:
        zeros=0
        for i in range(0,len(mapq),1):
            if (mapq[i] == 0):
                zeros += 1
        zerospersite=round(zeros/len(mapq),3)
      except:
        zeros=0
        zerospersite="NA"

      ## Calculate median base quality and map quality score
      try:
        medianbaseq=str(statistics.median(baseq))
      except:
        medianbaseq="NA"
      try:
        medianbaseqref=str(statistics.median(baseqref))
      except:
        medianbaseqref="NA"
      try:
        medianbaseqalt=str(statistics.median(baseqalt))
      except:
        medianbaseqalt="NA"
      try:
        medianmapq=str(statistics.median(mapq))
      except:
        medianmapq="NA"
      try:
        medianmapqref=str(statistics.median(mapqref))
      except:
        medianmapqref="NA"
      try:
        medianmapqalt=str(statistics.median(mapqalt))
      except:
        medianmapqalt="NA"

      ## Calculate median Levenschtein distances for this variant
      try:
        refld=statistics.median(reflevdists)
      except:
        refld="NA"
      try:
        altld=statistics.median(altlevdists)
      except:
        altld="NA"

      ## Calculate proportion of reads that are secondary alignments
      try:
        refsecondprop=round(refsecond/refcount,3)
      except:
        refsecondprop="NA"
      try:
        altsecondprop=round(altsecond/altcount,3)
      except:
        altsecondprop="NA"

      ## Calculate proportion of reads that have bad orientation (i.e. inversion orientation)
      try:
        refbadorientationprop=round(refbadorientation/refcount,3)
      except:
        refbadorientationprop="NA"
      try:
        altbadorientationprop=round(altbadorientation/altcount,3)
      except:
        altbadorientationprop="NA"


      ## How many contigs do mates map to?
      try:
        refmatecontigcount=len(numpy.unique((refmatecontigs)))
      except:
        refmatecontigcount="NA"
      try:
        altmatecontigcount=len(numpy.unique((altmatecontigs)))
      except:
        altmatecontigcount="NA"

      ## THIS IS WHERE WE WRITE OUTPUT
      outputstring=str(variant.CHROM)+":"+str(variant.POS)+":"+str(variant.REF)+":"+str(variant.ALT[0])+"\t"+\
      str(variant.CHROM)+"\t"+\
      str(variant.POS)+"\t"+\
      str(variant.REF)+"\t"+\
      str(variant.ALT[0])+"\t"+\
      str(refcount)+"\t"+\
      str(altcount)+"\t"+\
      str(varianttype)+"\t"+\
      str(depth)+"\t"+\
      str(vaf)+"\t"+\
      str(raf)+"\t"+\
      str(oaf)+"\t"+\
      str(medianbaseq)+"\t"+\
      str(medianbaseqref)+"\t"+\
      str(medianbaseqalt)+"\t"+\
      str(medianmapq)+"\t"+\
      str(medianmapqref)+"\t"+\
      str(medianmapqalt)+"\t"+\
      str(zeros)+"\t"+\
      str(zerospersite)+"\t"+\
      str(softreadlengthsrefmean)+"\t"+\
      str(softreadlengthsaltmean)+"\t"+\
      str(goodoffsetproportion)+"\t"+\
      str(distancetoend1median)+"\t"+\
      str(mad1)+"\t"+\
      str(distancetoend2median)+"\t"+\
      str(mad2)+"\t"+\
      str(distancetoend1medianref)+"\t"+\
      str(madref1)+"\t"+\
      str(distancetoend2medianref)+"\t"+\
      str(madref2)+"\t"+\
      str(distancetoend1medianalt)+"\t"+\
      str(madalt1)+"\t"+\
      str(distancetoend2medianalt)+"\t"+\
      str(madalt2)+"\t"+\
      str(shortestdistancetoendmedian)+"\t"+\
      str(madaltshort)+"\t"+\
      str(sb)+"\t"+\
      str(gsb)+"\t"+\
      str(fishp)+"\t"+\
      str(FoxoG)+"\t"+\
      str(refld)+"\t"+\
      str(altld)+"\t"+\
      str(refsecondprop)+"\t"+\
      str(altsecondprop)+"\t"+\
      str(refbadorientationprop)+"\t"+\
      str(altbadorientationprop)+"\t"+\
      str(refmatecontigcount)+"\t"+\
      str(altmatecontigcount)

      #print(outputstring) # print to screen
      print(outputstring,file=log) # print to log file
      
  
  ## Close the bamfile
  samfile.close()

  ## Close the logfile
  log.close()

  ## Announce completion
  logging.info("Finished processing chunk "+filename+" ("+str(chunknumber+1)+" of "+str(maxchunks)+") containing "+str(len(myvcf))+" records")

