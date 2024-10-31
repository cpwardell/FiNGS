## Primary function. This is called once for every chunk of SNVs
## To ensure speed, the BAM file(s) are read ONCE and ONCE only

import sys
import logging
import pysam
import statistics
import scipy.stats
import editdistance
import numpy
import vcfpy

def primary(myvcf,filename,chunknumber,maxchunks,args):

    logging.info("Started processing "+filename+" chunk "+str(chunknumber+1)+" of "+str(maxchunks)+" containing "+str(len(myvcf))+" records")

    ## Create a list to store the results from all variants in it
    allvariants=[]

    ## Open the bamfile. We do it ONCE and ONCE ONLY for speed reasons
    logging.debug("Opening bam file")
    samfile=pysam.Samfile(args.t,"rb")  # rb = "read bam"

    #################### START ITERATING THROUGH VARIANTS IN CHUNK ###########################
    for idx,variant in enumerate(myvcf):
        ## Determine if a variant is a SNV, insert or deletion
        varianttype=variant.ALT[0].type

        ## If the variant isn't a SNV, or is multiallelic, skip it
        if(varianttype!="SNV" or len(variant.REF)!=1): next

        ## Variable declarations:
        depth=0  # depth
        baseq,baseqref,baseqalt=[],[],[]  # basequality
        mapq,mapqref,mapqalt=[],[],[]  # mapping quality
        refcount,altcount=0,0  # number of REF and ALT reads
        refsecond,altsecond=0,0  # secondary alignment counts
        refbadorientation,altbadorientation=0,0  # orientation counts
        refmatecontigs,altmatecontigs=[],[]  # mate contig counts
        FR,FA,RR,RA=0,0,0,0  # 1st letter is strand: F forward, R reverse. 2nd letter is allele: R reference, A alternate
        F2R1,F1R2=float(),float()  # FoxoG orientation variables
        hardreadlengths,softreadlengths=[],[]  # read length and after soft-clipping
        softreadlengthsref,softreadlengthsalt=[],[]  # soft-clipped read lengths of REF/ALT containing reads
        offsets=[]  # this is the position in the read of the variant
        distancetoend1,distancetoend2=[],[]  # distance to ends from the variant position
        distancetoend1ref,distancetoend2ref=[],[]  # distance to ends from the variant position for REF reads
        distancetoend1alt,distancetoend2alt=[],[]  # distance to ends from the variant position for ALT reads
        reflevdists,altlevdists=[],[]  # Levenschtein distances

        ## If we are only using variants PASSed by the caller, skip the samfile pileup phase
        if(not args.PASSonlyin or (args.PASSonlyin and variant.FILTER[0]=="PASS")):
            #################### START ITERATING THROUGH READS AT VARIANT POSITION ###########################
            for alignedread in samfile.fetch(variant.CHROM,variant.POS-1,variant.POS):  # VCFs are 1 based, SAM coordinates are 0 based
                if alignedread.is_duplicate:  # ignore duplicates
                    continue
                else:
                    ## If depth is beyond maxdepth, break loop and calculate results and emit a warning to the user
                    if(depth>=args.m):
                        logging.info("WARNING: variant "+str(idx+1)+" of "+str(len(myvcf))+" in chunk "+str(chunknumber+1)+" of "+str(maxchunks)+" at position "+str(variant.CHROM)+":"+str(variant.POS)+" is deeper than maximum depth of "+str(args.m)+", calculations truncated at this depth")
                        break

                    ## Increment depth by 1 and append new mapping quality
                    depth+=1
                    mapq.append(alignedread.mapq)

                    ## Which base in the read is at the position we want?
                    offset=offsetget(alignedread.aligned_pairs,variant.POS)
                    if(offset is not None):
                        baseq.append(ord(alignedread.qual[offset])-33)  # append new base quality (we subtract 33 because SAM specification tells us to)
                        hardreadlengths.append(alignedread.infer_query_length())  # hard read lengths
                        softreadlengths.append(alignedread.query_alignment_end)  # soft-clipped read lengths
                        offsets.append(offset)
                        de1=abs(offset-alignedread.query_alignment_start)
                        de2=abs(offset-alignedread.query_alignment_end+1)
                        distancetoend1.append(de1)  # distance to lefthand soft-clipped read end
                        distancetoend2.append(de2)  # distance to righthand soft-clipped read end

                        ## Counters for REF allele reads
                        if(alignedread.seq[offset] == variant.REF[0]):
                            refcount+=1
                            if alignedread.is_secondary: refsecond+=1
                            if not alignedread.is_proper_pair: refbadorientation+=1
                            refmatecontigs.append(alignedread.next_reference_id)  # append the contig identifier of the mate to an array
                            baseqref.append(ord(alignedread.qual[offset])-33)  # append new base quality (we subtract 33 because SAM specification tells us to)
                            mapqref.append(alignedread.mapq)
                            softreadlengthsref.append(alignedread.query_alignment_end)  # soft-clipped read lengths
                            distancetoend1ref.append(de1)  # soft-clipped read length
                            distancetoend2ref.append(de2)  # soft-clipped read length
                            reflevdists.append(levcalc(alignedread))  # calculate edit distance
                            if not (alignedread.is_reverse):
                                FR+=1
                            else:
                                RR+=1

                        ## Counters for ALT alelle reads
                        if(alignedread.seq[offset] == str(variant.ALT[0].value)):
                            altcount+=1
                            if alignedread.is_secondary: altsecond+=1
                            if not alignedread.is_proper_pair: altbadorientation+=1
                            altmatecontigs.append(alignedread.next_reference_id)  # append the contig identifier of the mate to an array
                            baseqalt.append(ord(alignedread.qual[offset])-33)  # append new base quality (we subtract 33 because SAM specification tells us to)
                            mapqalt.append(alignedread.mapq)
                            softreadlengthsalt.append(alignedread.query_alignment_end)  # soft-clipped read lengths
                            distancetoend1alt.append(de1)  # soft-clipped read length
                            distancetoend2alt.append(de2)  # soft-clipped read length
                            altlevdists.append(levcalc(alignedread))  # calculate edit distance
                            if not (alignedread.is_reverse):
                                FA+=1
                            else:
                                RA+=1
                            ## Counters for FoxoG calculations
                            if(alignedread.is_read1 and alignedread.is_reverse): F2R1+=1  # 83/pPr1 F2R1
                            if(alignedread.is_read2 and alignedread.mate_is_reverse): F2R1+=1  # 163/pPR2 F2R1
                            if(alignedread.is_read1 and alignedread.mate_is_reverse): F1R2+=1  # 99/pPR1 F1R2
                            if(alignedread.is_read2 and alignedread.is_reverse): F1R2+=1  # 147/pPr1 F1R2

            ## Calculate metrics
            FoxoG=foxogcalc(variant.REF,F1R2,F2R1)  # FoxoG
            vaf=vafcalc(altcount,depth)  # VAF
            raf=vafcalc(refcount,depth)  # RAF
            oaf=oafcalc(refcount,altcount,depth)  # OAF
            altsb=oneallelestrandbias(FA,RA)
            refsb=oneallelestrandbias(FR,RR)
            allsb=simplestrandbias(FR,FA,RR,RA)
            sb=strandbias(FR,FA,RR,RA)  # Strand bias
            gsb=gatkstrandbias(FR,FA,RR,RA)  # GATK strand bias
            fishp=scipy.stats.fisher_exact([[FA,RA],[FR,RR]])[1]  # Fisher's exact test pvalue for strand bias
            fishp="%.3g" % fishp
            softmeanlength=roundmean(softreadlengths)  # soft-clipped read length mean
            softreadlengthsrefmean=roundmean(softreadlengthsref)  # mean length of REF soft-clipped reads
            softreadlengthsaltmean=roundmean(softreadlengthsalt)  # mean length of ALT soft-clipped reads
            goodoffsets=goodoffsetsfun(softmeanlength,offsets)
            goodoffsetproportion=goodoffsetproportionfun(goodoffsets,offsets)
            distancetoend1median=medianfun(distancetoend1)
            mad1=madfun(distancetoend1median,distancetoend1)
            distancetoend2median=medianfun(distancetoend2)
            mad2=madfun(distancetoend2median,distancetoend2)
            distancetoend1medianref=medianfun(distancetoend1ref)
            madref1=madfun(distancetoend1medianref,distancetoend1ref)
            distancetoend2medianref=medianfun(distancetoend2ref)
            madref2=madfun(distancetoend2medianref,distancetoend2ref)
            distancetoend1medianalt=medianfun(distancetoend1alt)
            madalt1=madfun(distancetoend1medianalt,distancetoend1alt)
            distancetoend2medianalt=medianfun(distancetoend2alt)
            madalt2=madfun(distancetoend2medianalt,distancetoend2alt)
            shortestdistancetoend=listmapminzip(distancetoend1alt,distancetoend2alt)
            shortestdistancetoendmedian=medianfun(shortestdistancetoend)
            madaltshort=madfun(shortestdistancetoendmedian,shortestdistancetoend)
            zeros,zerospersite=zerofun(mapq)
            medianbaseq=medianfun(baseq)
            medianbaseqref=medianfun(baseqref)
            medianbaseqalt=medianfun(baseqalt)
            medianmapq=medianfun(mapq)
            medianmapqref=medianfun(mapqref)
            medianmapqalt=medianfun(mapqalt)
            refld=medianfun(reflevdists)
            altld=medianfun(altlevdists)
            refsecondprop=roundpropfun(refsecond,refcount)
            altsecondprop=roundpropfun(altsecond,altcount)
            refbadorientationprop=roundpropfun(refbadorientation,refcount)
            altbadorientationprop=roundpropfun(altbadorientation,altcount)
            refmatecontigcount=lenunique(refmatecontigs)
            altmatecontigcount=lenunique(altmatecontigs)
            mtype=sixtypes(variant.REF,variant.ALT[0].value)  # what is the mutation type?

            ## Generate output string
            UID=str(variant.CHROM)+":"+str(variant.POS)+":"+str(variant.REF)+":"+str(variant.ALT[0].value)
            outputstring="\t".join(str(x) for x in [UID,
                                                    variant.CHROM,variant.POS,variant.REF,variant.ALT[0].value,
                                                    refcount,altcount,varianttype,depth,vaf,raf,oaf,
                                                    medianbaseq,medianbaseqref,medianbaseqalt,
                                                    medianmapq,medianmapqref,medianmapqalt,
                                                    zeros,zerospersite,
                                                    softreadlengthsrefmean,softreadlengthsaltmean,
                                                    goodoffsetproportion,
                                                    distancetoend1median,mad1,distancetoend2median,mad2,
                                                    distancetoend1medianref,madref1,distancetoend2medianref,madref2,
                                                    distancetoend1medianalt,madalt1,distancetoend2medianalt,madalt2,
                                                    shortestdistancetoendmedian,madaltshort,
                                                    sb,gsb,fishp,FR,FA,RR,RA,altsb,refsb,allsb,
                                                    F1R2,F2R1,FoxoG,
                                                    refld,altld,
                                                    refsecondprop,altsecondprop,
                                                    refbadorientationprop,altbadorientationprop,
                                                    refmatecontigcount,altmatecontigcount,
                                                    mtype])

            ## Append string to output list
            allvariants.append(outputstring)

    ## Close the bamfile
    samfile.close()

    ## Announce completion
    logging.info("Finished processing "+filename+" chunk "+str(chunknumber+1)+" of "+str(maxchunks)+" containing "+str(len(myvcf))+" records")

    ## Return variants for this chunk
    return(allvariants)


################ FUNCTIONS ######################

## Calculate extremely basic strand bias
def oneallelestrandbias(forward,reverse):
    try:
        sb=min(forward/(forward+reverse),reverse/(reverse+forward))
        sb=round(sb,3)
    except:
        sb="NA"
    return(sb)

def simplestrandbias(FR,FA,RR,RA):
    try:
        forward=FR+FA
        reverse=RR+RA
        sb=min(forward/(forward+reverse),reverse/(reverse+forward))
        sb=round(sb,3)
    except:
        sb="NA"
    return(sb)


## Calculate strand bias
def strandbias(FR,FA,RR,RA):
    try:
        sb=abs((FA/(FR+FA))-(RA/(RR+RA)))/((FA+RA)/(FR+FA+RR+RA))
        sb=round(sb,3)
    except:
        sb="NA"
    return(sb)


## Calculate strand bias as defined by GATK
def gatkstrandbias(FR,FA,RR,RA):
    try:
        sba=((FA/(FR+FA))*(RR/(RR+RA)))/((FR+RR)/(FR+FA+RR+RA))
        sba=round(sba,3)
    except:
        sba="NA"
    try:
        sbb=((RA/(RR+RA)) * (FR/(FR+FA)))/((FR+RR)/(FR+FA+RR+RA))
        sbb=round(sbb,3)
    except:
        sbb="NA"
    return(max(sba,sbb))


## Calculate allele frequency (VAF or RAF)
def vafcalc(allelecount,depth):
    try:
        vaf=round(float(allelecount)/(depth),3)
    except:
        vaf="NA"
    return(vaf)


## Calculate other allele frequency (OAF)
## "other" alleles are anything that isn't the REF or ALT allele in the VCF file
def oafcalc(refcount,altcount,depth):
    try:
        oaf=round(((depth)-float(altcount)-float(refcount))/(depth),3)
    except:
        oaf="NA"
    return(oaf)


## Calculation of FoxoG - only necessary if C>A|G>T SNV
## Equation is: ALT_F1R2/(ALT_F1R2 + ALT_F2R1) or ALT_F2R1/(ALT_F1R2 + ALT_F2R1)
## C>anything or A>anything:    numerator is ALT_F2R1
## G>anything or T>anything:    numerator is ALT_F1R2
## If sum of F1R2 and F2R1 is zero, all reads have an indel in them, so it should be removed
def foxogcalc(REF,F1R2,F2R1):
    FoxoG="NA"
    try:
        if((F1R2 + F2R1)!=0):
            if(REF=="C" or REF=="A"):
                FoxoG = F2R1/(F1R2 + F2R1)
            if(REF=="G" or REF=="T"):
                FoxoG = F1R2/(F1R2 + F2R1)
            FoxoG=round(FoxoG,3)
    except:
        FoxoG="NA"
    return(FoxoG)

## Determine the position of the variant within the read
def offsetget(aligned_pairs,POS):
    try:
        offset = [item for item in aligned_pairs if item[1] == POS-1][0][0]
    except:
        offset=None
    return(offset)


## Calculate Levenschtein distances (edit distances)
## get_reference_sequence() will fail without an MD tag in the BAM file, hence the exception
def levcalc(alignedread):
    try:
        distance=editdistance.eval(alignedread.get_reference_sequence(),alignedread.query_alignment_sequence)
    except:
        distance=1
    return(distance)


## Number of offsets that are less than 2/3rd of the mean read length
## This can fail if the offsets list is empty
def goodoffsetsfun(softmeanlength,offsets):
    try:
        goodoffsets=len(list(filter(lambda x: x<(2/3*softmeanlength),offsets)))
    except:
        goodoffsets="NA"
    return(goodoffsets)


## Proportion of offsets that are good (see goodoffsetsfun())
def goodoffsetproportionfun(goodoffsets,offsets):
    try:
        goodoffsetproportion=round(goodoffsets/len(offsets),3)
    except:
        goodoffsetproportion="NA"
    return(goodoffsetproportion)


## Median distance to alignment ends
def medianfun(alist):
    try:
        themedian=round(statistics.median(alist),3)
    except:
        themedian="NA"
    return(themedian)


## Median Absolute Deviation (MAD) of distance to alignment ends
def madfun(themedian,alist):
    try:
        mad=round(statistics.median([abs(x - themedian) for x in alist]),3)
    except:
        mad="NA"
    return(mad)


def listmapminzip(list1,list2):
    try:
        minlist=list(map(min,zip(list1,list2)))
    except:
        minlist="NA"
    return(minlist)


## Count how many MAPQ==0 reads there are and normalize for depth
def zerofun(mapq):
    try:
        zeros=0
        for i in range(0,len(mapq),1):
                if (mapq[i] == 0):
                        zeros += 1
        zerospersite=round(zeros/len(mapq),3)
    except:
        zeros=0
        zerospersite="NA"
    return([zeros,zerospersite])


def roundpropfun(num,denom):
    try:
        roundprop=round(num/denom,3)
    except:
        roundprop="NA"
    return(roundprop)


def lenunique(alist):
    try:
        lenuniquecount=len(numpy.unique((alist)))
    except:
        lenuniquecount="NA"
    return(lenuniquecount)


## Add mutation type
def sixtypes(ref,alt):
    ref=str(ref)
    alt=str(alt)
    if((ref=="T" and alt=="G") or (ref=="A" and alt=="C")):
        mtype="T>G/A>C"
    if((ref=="T" and alt=="C") or (ref=="A" and alt=="G")):
        mtype="T>C/A>G"
    if((ref=="T" and alt=="A") or (ref=="A" and alt=="T")):
        mtype="T>A/A>T"
    if((ref=="C" and alt=="T") or (ref=="G" and alt=="A")):
        mtype="C>T/G>A"
    if((ref=="C" and alt=="G") or (ref=="G" and alt=="C")):
        mtype="C>G/G>C"
    if((ref=="C" and alt=="A") or (ref=="G" and alt=="T")):
        mtype="C>A/G>T"
    if((ref not in ["A","C","G","T"]) or (alt not in ["A","C","G","T"])):
        mtype="INDEL"
    return(mtype)


## Function for calculating the mean; we don't want to use numpy
def roundmean(numbers):
        try:
            x=round(float(sum(numbers))/len(numbers),3)
        except:
            x=0
        return(x)
