## Functions to perform filtering

import logging
import sys
import os
import gzip
import csv
import scipy.stats
import vcf
import pysam
import numpy
import pandas
import seaborn as sns
import statsmodels.nonparametric.api as smnp
import matplotlib.backends.backend_pdf

from matplotlib import pyplot
from collections import OrderedDict
from itertools import repeat, product

## Define column names for tdata and ndata
cnames=["UID","CHR","POS","REF","ALT","refcount","altcount","varianttype","depth","vaf","raf","oaf",
        "medianbaseq","medianbaseqref","medianbaseqalt","medianmapq","medianmapqref","medianmapqalt",
        "zeros","zerospersite","softreadlengthsrefmean","softreadlengthsaltmean","goodoffsetproportion",
        "distancetoend1median","mad1","distancetoend2median","mad2","distancetoend1medianref","madref1",
        "distancetoend2medianref","madref2","distancetoend1medianalt","madalt1","distancetoend2medianalt",
        "madalt2","shortestdistancetoendmedian","madaltshort","sb","gsb","fishp",
        "FR","FA","RR","RA","altsb","refsb","allsb","F1R2","F2R1","FoxoG","refld","altld",
        "refsecondprop","altsecondprop","refbadorientationprop","altbadorientationprop","refmatecontigcount",
        "altmatecontigcount","sixtypes"]

## Define some colors for plotting
bloodred="#870202"
badred="#E41A1C"
goodgreen="#4DAF4A"

## Main filter function; a list of the filters that are applied
def applyfilters(tdata,ndata,sdata,pdict,resultsdir,vcfpath,referencegenome,PASSonlyin,PASSonlyout):

    ## Read in metrics and apply filters
    ## We read the data back from a temporary file in case to make it easier to recover from crashes
    ## and to speed up filtering; gathering the metrics is the rate limiting step

    ## Use input VCF as a template
    vcf_reader = vcf.Reader(open(vcfpath, 'r'))

    ## Read in summary data from file
    try:
        with gzip.open(sdata,"rt",newline="") as sd:
            scsv = csv.reader(sd, delimiter='\t')
            sdict={}
            for row in scsv:
                sdict[row[0]]=float(row[1])
    except Exception as e:
        print("CRITICAL ERROR: unreadable file: "+str(e))
        sys.exit()

    ## Santize list of filters; does the pdict object only contain real filter names?
    filtercheck(pdict)

    ## Generate list of all possible repeats to check for with the repeats filter
    try:
        repeatlist=generaterepeats(int(pdict["repeats"]))
    except:
        repeatlist=generaterepeats(12)  # if no value is supplied, use 12 as a default

    ## Create a list of filter results
    flist=[]
    for pair in metricyielder(tdata,ndata):
        ## Create a dictionary of filters and populate it
        fdict=filterbuild(pair,pdict,sdict,referencegenome,repeatlist,vcf_reader)
        flist.append(fdict)

    ## Now we have all mutations in an object, we can run filters that reference other mutations
    filtergroup(flist,pdict,vcf_reader)

    ## Write the filter results (a list of dictionaries) to a file
    fdata=resultsdir+"/filterresults.txt.gz"
    logging.info("Writing filter results to file: "+fdata)
    ordered_fieldnames = flist[0].keys()
    with gzip.open(fdata,'wt') as fout:
        dw=csv.DictWriter(fout, delimiter='\t', fieldnames=ordered_fieldnames)
        dw.writeheader()
        for fdict in flist:
            dw.writerow(fdict)
    fout.close()

    ####################################################
    ## Use the filter results to write VCF results
    vfile=resultsdir+"/"+os.path.splitext(os.path.basename(vcfpath))[0]+".filtered.vcf"
    logging.info("Writing filtered VCF to file: "+vfile)

    with open(vfile, "w") as out_vcf:
        writer = vcf.Writer(out_vcf, vcf_reader)
        passindex=0
        for idx,record in enumerate(vcf_reader):
            ## Working with all variants, emit all variants
            if(not PASSonlyin and not PASSonlyout):
                record.FILTER=judgement(flist[idx])
                writer.write_record(record)

            ## Working with all variants, emit only PASS variants
            if(not PASSonlyin and PASSonlyout):
                record.FILTER=judgement(flist[idx])
                if("PASS" in record.FILTER):
                    writer.write_record(record)

            ## Working with only PASS variants, emit all variants
            if(PASSonlyin and not PASSonlyout):
                if(len(record.FILTER) == 0):
                    record.FILTER=judgement(flist[passindex])
                    passindex+=1
                    writer.write_record(record)

            ## Working with only PASS variants, emit only PASS variants
            if(PASSonlyin and PASSonlyout):
                if(len(record.FILTER) == 0):
                    record.FILTER=judgement(flist[passindex])
                    passindex+=1
                    if("PASS" in record.FILTER):
                        writer.write_record(record)

    writer.close()
    logging.info("Writing VCF complete")

    ## Produce plots
    ## To produce plots, we need the whole dataset;
    ## read in tdata/ndata
    ## Loop through pdict and produce plots according to what's in there
    logging.info("Writing plots to file: "+resultsdir+"/plots.pdf")
    filterplotter(tdata,ndata,fdata,pdict,sdict,resultsdir)
    logging.info("Writing plots complete")

def filterplotter(tdata,ndata,fdata,pdict,sdict,resultsdir):
    ## Read in tdata and ndata
    try:
        tdf = pandas.read_csv(tdata, compression='gzip', sep="\t", header=None, names=cnames)
        ndf = pandas.read_csv(ndata, compression='gzip', sep="\t", header=None, names=cnames)

        ## Crate pdf output file and put plots in it; summary table first
        pdf = matplotlib.backends.backend_pdf.PdfPages(resultsdir+"/plots.pdf")
        summaryplots(fdata,pdict,pdf)
        for key in pdict.keys():
            if key == "minimumdepth":
                densityplot(tdf["depth"],pdict["minimumdepth"],"Minimum tumor depth: ","right",pdf)
                densityplot(ndf["depth"],pdict["minimumdepth"],"Minimum normal depth: ","right",pdf)
            if key == "maximumdepth":
                densityplot(tdf["depth"],pdict["maximumdepth"],"Maximum tumor depth: ","left",pdf)
                densityplot(ndf["depth"],pdict["maximumdepth"],"Maximum normal depth: ","left",pdf)
            if key == "minaltcount":
                densityplot(tdf["altcount"],pdict["minaltcount"],"Minimum number of ALT reads in tumor: ","right",pdf)
            if key == "minbasequality":
                densityplot(tdf["medianbaseqalt"],pdict["minbasequality"],"Minimum median base quality for ALT in tumor: ","right",pdf)
                densityplot(tdf["medianbaseqref"],pdict["minbasequality"],"Minimum median base quality for REF in tumor: ","right",pdf)
                densityplot(ndf["medianbaseqref"],pdict["minbasequality"],"Minimum median base quality for REF in normal: ","right",pdf)
            if key == "zeroproportion":
                densityplot(tdf["zerospersite"],pdict["zeroproportion"],"Maximum proportion of zero mapping quality reads in tumor: ","left",pdf)
                densityplot(ndf["zerospersite"],pdict["zeroproportion"],"Maximum proportion of zero mapping quality reads in normal: ","left",pdf)
            if key == "minmapquality":
                densityplot(tdf["medianmapqalt"],pdict["minmapquality"],"Minimum median mapping quality of ALT reads in tumor: ","right",pdf)
            if key == "enddistance":
                densityplot(tdf["shortestdistancetoendmedian"],pdict["enddistance"],"Minimum median shortest distance\nto either aligned end in tumor: ","right",pdf)
            if key == "enddistancemad":
                densityplot(tdf["madaltshort"],pdict["enddistancemad"],"Minimum MAD of median shortest distance\nto either aligned end in tumor: ","right",pdf)
            if key == "editdistance":
                densityplot(tdf["altld"],pdict["editdistance"],"Maximum edit distance of ALT reads in tumor: ","left",pdf)
                densityplot(tdf["refld"],pdict["editdistance"]-1,"Maximum edit distance of REF reads in tumor: ","left",pdf)
            if key == "minvaftumor":
                densityplot(tdf["vaf"],pdict["minvaftumor"],"Minimum VAF in tumor: ","right",pdf)
            if key == "maxoaftumor":
                densityplot(tdf["oaf"],pdict["maxoaftumor"],"Maximum other allele frequency (OAF) in tumor: ","left",pdf)
            if key == "maxsecondtumor":
                densityplot(tdf["altsecondprop"],pdict["maxsecondtumor"],"Maximum proportion of secondary alignments in tumor: ","left",pdf)
            if key == "strandbiasprop":
                densityplot(tdf["sb"],sdict["strandbiastumorq"],"Strand bias exclusion proportion: "+str(pdict["strandbiasprop"])+"\nThreshold: ","left",pdf)
            if key == "strandbias":
                densityplot(tdf["sb"],pdict["strandbias"],"Maximum strand bias in tumor: ","left",pdf)
            if key == "strandbiassimple":
                densityplot(tdf["altsb"],pdict["strandbiassimple"],"Maximum strand bias of ALT reads in tumor: ","right",pdf)
            if key == "maxaltcount":
                densityplot(ndf["altcount"],pdict["maxaltcount"],"Maximum number of ALT reads in normal: ","left",pdf)
            if key == "maxvafnormal":
                densityplot(ndf["vaf"],pdict["maxvafnormal"],"Maximum VAF in normal: ","left",pdf)
            if key == "maxbadorient":
                densityplot(ndf["refbadorientationprop"],pdict["maxbadorient"],"Maximum proportion of inversion orientation\nreads in normal: ","left",pdf)
            if key == "minmapqualitydifference":
                densityplot(abs(ndf["medianmapqref"]-tdf["medianmapqalt"]),pdict["minmapqualitydifference"],"Maximum difference between median mapping quality\nof ALT reads in tumor and REF reads in normal: ","left",pdf)
#            if key == "foxog":  # NEED AN X Y SCATTER PLOT
#                scatterplot(tdf["F1R2"],tdf["F2R1"],tdf["sixtypes"],"FoxoG plot",pdf)
#                fdict["foxogtumor"]=filterfoxog(pair[0]["REF"],pair[0]["sixtypes"],pair[0]["F1R2"],pair[0]["F2R1"],pdict["foxog"])  # TUMOR
#                vcf_reader.filters[18]=["foxogtumor","Maximum FoxoG artefact proportion: "+str(pdict["foxog"])]
        pdf.close()

    except Exception as e:
        print("CRITICAL ERROR: error during plotting: "+str(e))
        sys.exit()

## Summary plot function
def summaryplots(fdata,pdict,pdf):
    ## Read filter data back in
    fdf = pandas.read_csv(fdata, compression='gzip', sep="\t", header=0)
    
    ## Get rid of non-filter columns
    fdf=fdf[fdf.columns[range(4,fdf.shape[1])]]
    
    ## Create summary dataframe
    sumfdf=pandas.DataFrame(fdf,index=["PASS","Fail"]) 
    sumfdf["TOTAL"]=-1 # Assign value to total column so it is always sorted to top
    for column in fdf.columns:
        sumfdf[column]=[list(fdf[column]).count(True),list(fdf[column]).count(False)]
    sumfdf=sumfdf.transpose()
    sumfdf=sumfdf.sort_values(by="PASS")

    ## Assign real counts to total row
    ## If all elements sum to length of row, increment counter
    totalpass=0
    for idx, row in fdf.iterrows():
        if sum(row) == fdf.shape[1]:
            totalpass+=1
    sumfdf.loc["TOTAL"]=[totalpass,fdf.shape[0]-totalpass]

    ## Create axes for plot and set font size for labels
    ## [left, bottom, width, height]
    sns.set(font_scale=0.60,style="white") 
    ax=pyplot.axes([0.21,0.11,0.775,0.77]) 

    ## Produce plot and reset size of annotation text in boxes
    plt=sns.heatmap(sumfdf,annot=True, fmt="d", cbar=False, cmap="YlGnBu",ax=ax,annot_kws={"size": 10})#,ax=ax)
    #plt.get_figure().savefig("results/summary.png")
    pdf.savefig()
    plt.get_figure().clf() # clears figure to prevent multiple overlapping plots

## Scatterplot function for foxog data
#def scatterplot(F1,F2,sixtypes,foxocut):

def densityplot(pcol,vline,title,goodside,pdf):
    ## Discard any nan values
    pcol=pcol.dropna()
    ## Produces the initial seaborn plot
    ## We suppress warnings in case empty lists are passed in
    numpy.seterr(divide='ignore', invalid='ignore')
    plt=sns.distplot(pcol, color="black", kde_kws={"shade": False}, rug=True, rug_kws={"height": 0.02, "color": "black"}, hist=False)#,hist_kws={"color": "grey", "alpha": 0})
    plt.axvline(vline,linestyle="--",color=bloodred)
    plt.set_title(title+str(round(vline,3)))

    ## This block calculates the areas either side of the threshold line for shading
    ## We ensure y values are never less than 0
    kde = smnp.KDEUnivariate(pcol.astype(numpy.float64))
    kde.fit("gau", "scott", "gau", gridsize=100, cut=3, clip=(-numpy.inf, numpy.inf))
    x, y = kde.support, kde.density
    y=[max(x,0) for x in y] 
 
    ## Apply shading
    if(goodside=="left"):
        plt.fill_between(x,0,y,color=goodgreen, where = x < vline)
        plt.fill_between(x,0,y,color=badred, alpha=1, where = x > vline)
    if(goodside=="right"):
        plt.fill_between(x,0,y,color=badred, where = x < vline)
        plt.fill_between(x,0,y,color=goodgreen, alpha=1, where = x > vline)
    pdf.savefig()
    plt.get_figure().clf() # clears figure to prevent multiple overlapping plots



def filtergroup(flist,pdict,vcf_reader):
    ## Iterate through list of filters and add any as required
    for key in pdict.keys():
        if key == "snvcluster50":
            vcf_reader.filters[98]=["snvcluster50","Maximum number of mutations within 50 bp: "+str(int(pdict["snvcluster50"]))]
            snvcluster(flist,50,pdict["snvcluster50"])
        if key == "snvcluster100":
            vcf_reader.filters[99]=["snvcluster100","Maximum number of mutations within 100 bp: "+str(int(pdict["snvcluster100"]))]
            snvcluster(flist,100,pdict["snvcluster100"])


def snvcluster(flist,nlength,nsnvs):
    nlength=int(nlength)
    nsnvs=int(nsnvs)
    ## Iterate through all mutations
    for idx,fdict in enumerate(flist):
        ## We want to look at nsnv+1 slots before and after the current index
        slots=numpy.array(range(idx-(nsnvs+1),idx+nsnvs+2))
        slots=slots[slots>=0]  # exclude negative values as they are out of bounds
        slots=slots[slots<=len(flist)-1]  # exclude positive values that are out of out of bounds
        slots=set(slots)-set([idx])  # exclude the current index slot
        ## Iterate through the slots and count how many match the criteria
        ## 1) Same chromosome
        ## 2) Within nlength/2 bases
        ## 3) Minimum 5% VAF (mutations with no VAF fail this)
        ## 4) Minimum 2 supporting reads
        localmuts=0
        for slot in slots:
            strikes=0  # counter for matching criteria
            if(flist[slot]["CHR"]==fdict["CHR"]):
                strikes+=1
            if(abs(int(flist[slot]["POS"])-int(fdict["POS"]))<=round(nlength/2)):
                strikes+=1
            try:
                if(float(flist[slot]["vaf"])>=0.05):
                    strikes+=1
            except:
                strikes+=1
            if(int(flist[slot]["altcount"])>=2):
                strikes+=1
            ## If all conditions are met, iterate counter
            if(strikes==4):
                localmuts+=1

        ## If there are nsnvs or more, fail the filter
        if(localmuts>=nsnvs):
            fdict["snvcluster"+str(nlength)]=False
        else:
            fdict["snvcluster"+str(nlength)]=True


## Define which repeats we need to look for; repeat length must be a factor of sequence length
## e.g. we don't look for 4mers in sequence of length 13
def generaterepeats(nlength):
    ## Always create FOR 1mers
    allmers=[]
    allmers.append(''.join(list(repeat("A",nlength))))
    allmers.append(''.join(list(repeat("C",nlength))))
    allmers.append(''.join(list(repeat("G",nlength))))
    allmers.append(''.join(list(repeat("T",nlength))))

    ## All possible 2mers
    if((nlength % 2)==0):
        y=list(product("ACGT","ACGT"))
        z=["%s%s" % x for x in y]
        for mer in z:
            full=''.join(list(repeat(mer,int(nlength/2))))
            allmers.append(full)

    ## All possible 3mers
    if((nlength % 3)==0):
        y=list(product("ACGT","ACGT","ACGT"))
        z=["%s%s%s" % x for x in y]
        for mer in z:
            full=''.join(list(repeat(mer,int(nlength/3))))
            allmers.append(full)

    ## All possible 4mers
    if((nlength % 4)==0):
        y=list(product("ACGT","ACGT","ACGT","ACGT"))
        z=["%s%s%s%s" % x for x in y]
        for mer in z:
            full=''.join(list(repeat(mer,int(nlength/4))))
            allmers.append(full)

    ## Cast list to a set; this makes all elements unique
    allmers=set(allmers)
    return(allmers)


## Yields pairs of lists of metrics for filtering
def metricyielder(tdata,ndata):

    logging.info("Reading tumor data from this file: "+tdata)
    logging.info("Reading normal data from this file: "+ndata)
    try:
        with gzip.open(tdata,"rt",newline="") as td, gzip.open(ndata,"rt",newline="") as nd:
            tcsv = csv.reader(td, delimiter='\t')
            ncsv = csv.reader(nd, delimiter='\t')

            ## We loop through both the tumor and normal data together and yield pairs of metrics 
            for trow, nrow in zip(tcsv,ncsv):
                ## Convert list to dictionary
                tdict={}
                ndict={}
                for i in range(0,len(cnames)):
                    tdict[cnames[i]]=trow[i]
                    ndict[cnames[i]]=nrow[i]
                yield(tdict,ndict)
    except Exception as e:
        print("CRITICAL ERROR: unreadable file: "+str(e))
        sys.exit()


## Minimum depth
def filterminimumdepth(testvalue,minimumdepth):
    result=False
    if(float(testvalue) >= minimumdepth):
        result=True
    return(result)


## Maximum depth
def filtermaximumdepth(testvalue,maximumdepth):
    result=False
    if(float(testvalue) < maximumdepth):
        result=True
    return(result)


## Maximum depth
def filterminaltcount(testvalue,minaltcount):
    result=False
    if(float(testvalue) >= minaltcount):
        result=True
    return(result)


## Maximum number of ALT reads in normal
def filtermaxaltcount(testvalue,maxaltcount):
    result=False
    if(float(testvalue) <= maxaltcount):
        result=True
    return(result)


## Minimum base
def filterminbasequality(testvalue,minimumbaseq,naaction=False):
    ## Cast testvalue to float. If set to NA, follow instructions below
    ## NA ALT in tumor means no ALT reads - FAIL
    ## NA REF in tumor means no REF reads; could be homozygous ALT or shallow - PASS
    ## NA REF in normal means no REF reads - FAIL
    result=False
    isna=False
    try:
        testvalue=float(testvalue)
    except:
        testvalue=0
        isna=True
    if(testvalue >= minimumbaseq):
        result=True
    ## Set filter to pass if NA action is set to true
    if(isna):
        result=naaction
    return(result)


def filterzeroproportion(testvalue,zeroproportion):
    result=False
    try:
        testvalue=float(testvalue)
    except:
        testvalue=1
    if(testvalue < zeroproportion):
        result=True
    return(result)


def filterminmapquality(testvalue,minmapquality):
    result=False
    try:
        testvalue=float(testvalue)
    except:
        testvalue=0
    if(testvalue >= minmapquality):
        result=True
    return(result)


#################### Maximum mapping quality difference ############
## Median mapping quality of alt reads in the tumor must be within
## minmapqualitydifference of median mapping quality of ref reads in the normal
def filterminmapqualitydifference(testvaluet,testvaluen,minmapqualitydifference):
    result=False
    try:
        testvalue=abs(float(testvaluet)-float(testvaluen))
    except:
        testvalue=1000
    if(testvalue < minmapqualitydifference):
        result=True
    return(result)


def filterenddistance(testvalue,enddistance):
    result=False
    try:
        testvalue=float(testvalue)
    except:
        testvalue=0
    if(testvalue > enddistance):
        result=True
    return(result)


def filterenddistancemad(testvalue,enddistancemad):
    result=False
    try:
        testvalue=float(testvalue)
    except:
        testvalue=0
    if(testvalue > enddistancemad):
        result=True
    return(result)


## Fail zero alt records, pass zero ref records
def filtereditdistance(testvalue,editdistance,naaction=False):
    result=False
    isna=False
    try:
        testvalue=float(testvalue)
    except:
        testvalue=0
        isna=True
    if(testvalue <= editdistance):
        result=True
    ## Set filter to pass if NA action is set to true
    if(isna):
        result=naaction
    return(result)


def filtermaxvaf(testvalue,maxvaf):
    result=False
    try:
        testvalue=float(testvalue)
    except:
        testvalue=1
    if(testvalue < maxvaf):
        result=True
    return(result)


def filterminvaf(testvalue,minvaf):
    result=False
    try:
        testvalue=float(testvalue)
    except:
        testvalue=0
    if(testvalue >= minvaf):
        result=True
    return(result)


def filtermaxaltsecond(testvalue,maxsecond):
    result=False
    try:
        testvalue=float(testvalue)
    except:
        testvalue=1000
    if(testvalue < maxsecond):
        result=True
    return(result)


def filtermaxbadorient(testvalue,maxbadorient):
    result=False
    try:
        testvalue=float(testvalue)
    except:
        testvalue=1000
    if(testvalue < maxbadorient):
        result=True
    return(result)


def filterfoxog(ref,sixtypes,F1R2,F2R1,foxog):
    result=False

    ## Cast F1R2 and F2R1 to floats
    try:
        F1R2=float(F1R2)
        F2R1=float(F2R1)
    except:
        F1R2=1
        F2R1=1

    ## Separate tests for REF C and REF G
    if(sixtypes == "C>A/G>T" and ref == "C"):
        phalf=scipy.stats.binom_test(F2R1,F2R1+F1R2,0.5)
        pthreshold=scipy.stats.binom_test(F2R1,F2R1+F1R2,foxog)
        if(phalf > pthreshold):
            result=True

    if(sixtypes == "C>A/G>T" and ref == "G"):
        phalf=scipy.stats.binom_test(F1R2,F1R2+F2R1,0.5)
        pthreshold=scipy.stats.binom_test(F1R2,F1R2+F2R1,foxog)
        if(phalf > pthreshold):
            result=True

    ## Pass any non C>A/G>T mutations
    if(sixtypes != "C>A/G>T"):
        result=True
    return(result)


## Variants with no strand bias value FAIL by default
def filterstrandbiasprop(testvalue,strandbiastumorq):
    result=False
    try:
        testvalue=float(testvalue)
    except:
        testvalue=1000
    if(testvalue < strandbiastumorq):
        result=True
    return(result)


## Variants with no strand bias value FAIL by default
def filterstrandbias(testvalue,strandbias):
    result=False
    try:
        testvalue=float(testvalue)
    except:
        testvalue=0
    if(testvalue <= strandbias):
        result=True
    return(result)

## Variants with no strand bias value PASS by default
def filterstrandbiassimple(altsb,allsb,strandbias):
    result=False
    try:
        altsb=float(altsb)
        allsb=float(allsb)
    except:
        altsb=0
        allsb=0
    if((altsb > strandbias) or (altsb<=strandbias and allsb<=strandbias)):
        result=True
    return(result)

## Repeat filter
def filterrepeats(chr,pos,nlength,repeatlist,referencegenome):
    result=False
    ## Ensure that a reference genome has been supplied by the user:
    if(referencegenome=="None"):
        print("CRITICAL ERROR: if \'repeats\' filter is in use, you MUST supply the absolute path to a faidx indexed reference genome (-r flag)")
        sys.exit()

    ## Fetch the n bases upstream and downstream of the variant; note we reverse the upstream sequence
    refgenome=pysam.FastaFile(referencegenome)
    upstream=refgenome.fetch(str(chr),int(pos)-nlength-1,int(pos)-1)[::-1]
    downstream=refgenome.fetch(str(chr),int(pos),int(pos)+nlength)

    ## Only PASS if the upstream and downstream sequences are not in the list of repeats
    if((upstream not in repeatlist) and (downstream not in repeatlist)):
        result=True

    return(result)


## Ensure that list of filters received is ok
def filtercheck(pdict):
    logging.info("Checking filter names")
    truekeys=["minimumdepth",
              "maximumdepth",
              "minaltcount",
              "maxaltcount",
              "minbasequality",
              "zeroproportion",
              "strandbiasprop",
              "strandbias",
              "strandbiassimple",
              "repeats",
              "minmapquality",
              "minmapqualitydifference",
              "enddistance",
              "enddistancemad",
              "editdistance",
              "maxvafnormal",
              "minvaftumor",
              "maxoaftumor",
              "maxsecondtumor",
              "maxbadorient",
              "foxog",
              "snvcluster50",
              "snvcluster100"]
    for key in pdict.keys():
        if key in truekeys:
            continue
        else:
            print("CRITICAL ERROR: check filter parameter file (-p flag, default is \"filter_parameters.txt\"),\nthis filter is not in the acceptable list:\t"+str(key)+"\nAcceptable filters are:"+str(truekeys))
            sys.exit()
    logging.info("Checking filter names - complete")


## Returns a dictionary of filters
## and edits the header of the output VCF
def filterbuild(pair,pdict,sdict,referencegenome,repeatlist,vcf_reader):
    fdict=OrderedDict()
    fdict["CHR"]=pair[0]["CHR"]
    fdict["POS"]=pair[0]["POS"]
    fdict["altcount"]=pair[0]["altcount"]
    fdict["vaf"]=pair[0]["vaf"]
    vcf_reader.filters=OrderedDict()  # reset the VCF filter descriptions in the header
    for key in pdict.keys():
        if key == "minimumdepth":
            fdict["minimumdepthtumor"]=filterminimumdepth(pair[0]["depth"],pdict["minimumdepth"])  # TUMOR
            vcf_reader.filters[0]=["minimumdepthtumor","Minimum depth in tumor: "+str(int(pdict["minimumdepth"]))]
            fdict["minimumdepthnormal"]=filterminimumdepth(pair[1]["depth"],pdict["minimumdepth"])  # NORMAL
            vcf_reader.filters[1]=["minimumdepthnormal","Minimum depth in normal: "+str(int(pdict["minimumdepth"]))]
        if key == "maximumdepth":
            fdict["maximumdepthtumor"]=filtermaximumdepth(pair[0]["depth"],pdict["maximumdepth"])  # TUMOR
            vcf_reader.filters[2]=["maximumdepthtumor","Maximum depth in tumor: "+str(int(pdict["maximumdepth"]))]
            fdict["maximumdepthnormal"]=filtermaximumdepth(pair[1]["depth"],pdict["maximumdepth"])  # NORMAL
            vcf_reader.filters[3]=["maximumdepthnormal","Maximum depth in normal: "+str(int(pdict["maximumdepth"]))]
        if key == "minaltcount":
            fdict["minaltcounttumor"]=filterminaltcount(pair[0]["altcount"],pdict["minaltcount"])  # TUMOR
            vcf_reader.filters[4]=["minaltcounttumor","Minimum number of ALT reads in tumor: "+str(int(pdict["minaltcount"]))]
        if key == "minbasequality":
            fdict["minbasequalityalttumor"]=filterminbasequality(pair[0]["medianbaseqalt"],pdict["minbasequality"])  # TUMOR
            vcf_reader.filters[5]=["minbasequalityalttumor","Minimum median base quality of ALT reads in tumor: "+str(pdict["minbasequality"])]
            fdict["minbasequalityreftumor"]=filterminbasequality(pair[0]["medianbaseqref"],pdict["minbasequality"],True)  # TUMOR
            vcf_reader.filters[6]=["minbasequalityreftumor","Minimum median base quality of REF reads in tumor: "+str(pdict["minbasequality"])]
            fdict["minbasequalityrefnormal"]=filterminbasequality(pair[1]["medianbaseqref"],pdict["minbasequality"])  # NORMAL
            vcf_reader.filters[7]=["minbasequalityrefnormal","Minimum median base quality of REF reads in normal: "+str(pdict["minbasequality"])]
        if key == "zeroproportion":
            fdict["zeroproportiontumor"]=filterzeroproportion(pair[0]["zerospersite"],pdict["zeroproportion"])  # TUMOR
            vcf_reader.filters[8]=["zeroproportiontumor","Maximum proportion of zero mapping quality reads in tumor: "+str(pdict["zeroproportion"])]
            fdict["zeroproportionnormal"]=filterzeroproportion(pair[1]["zerospersite"],pdict["zeroproportion"])  # NORMAL
            vcf_reader.filters[9]=["zeroproportionnormal","Maximum proportion of zero mapping quality reads in normal: "+str(pdict["zeroproportion"])]
        if key == "minmapquality":
            fdict["minmapqualityalttumor"]=filterminmapquality(pair[0]["medianmapqalt"],pdict["minmapquality"])  # TUMOR
            vcf_reader.filters[10]=["minmapqualityalttumor","Minimum median mapping quality of ALT reads in tumor: "+str(pdict["minmapquality"])]
        if key == "enddistance":
            fdict["enddistancetumor"]=filterenddistance(pair[0]["shortestdistancetoendmedian"],pdict["enddistance"])  # TUMOR
            vcf_reader.filters[11]=["enddistancetumor","Minimum median shortest distance to either aligned end in tumor: "+str(pdict["enddistance"])]
        if key == "enddistancemad":
            fdict["enddistancemadtumor"]=filterenddistancemad(pair[0]["madaltshort"],pdict["enddistancemad"])  # TUMOR
            vcf_reader.filters[12]=["enddistancemadtumor","Minimum median absolute deviation (MAD) of median shortest distance to either aligned end in tumor: "+str(pdict["enddistancemad"])]
        if key == "editdistance":
            fdict["editdistancealttumor"]=filtereditdistance(pair[0]["altld"],pdict["editdistance"])  # TUMOR
            vcf_reader.filters[13]=["editdistancealttumor","Maximum edit distance of ALT reads in tumor: "+str(pdict["editdistance"])]
            fdict["editdistancereftumor"]=filtereditdistance(pair[0]["refld"],pdict["editdistance"]-1,True)  # TUMOR
            vcf_reader.filters[14]=["editdistancereftumor","Maximum edit distance of REF reads in tumor: "+str(pdict["editdistance"])]
        if key == "minvaftumor":
            fdict["minvaftumor"]=filterminvaf(pair[0]["vaf"],pdict["minvaftumor"])  # TUMOR
            vcf_reader.filters[15]=["minvaftumor","Minimum VAF in tumor: "+str(pdict["minvaftumor"])]
        if key == "maxoaftumor":
            fdict["maxoaftumor"]=filtermaxvaf(pair[0]["oaf"],pdict["maxoaftumor"])  # TUMOR
            vcf_reader.filters[16]=["maxoaftumor","Maximum other allele frequency (OAF) in tumor: "+str(pdict["maxoaftumor"])]
        if key == "maxsecondtumor":
            fdict["maxaltsecondtumor"]=filtermaxaltsecond(pair[0]["altsecondprop"],pdict["maxsecondtumor"])  # TUMOR
            vcf_reader.filters[17]=["maxaltsecondtumor","Maximum proportion of secondary alignments in tumor: "+str(pdict["maxsecondtumor"])]
        if key == "foxog":
            fdict["foxogtumor"]=filterfoxog(pair[0]["REF"],pair[0]["sixtypes"],pair[0]["F1R2"],pair[0]["F2R1"],pdict["foxog"])  # TUMOR
            vcf_reader.filters[18]=["foxogtumor","Maximum FoxoG artefact proportion: "+str(pdict["foxog"])]
        if key == "strandbiasprop":
            fdict["strandbiasprop"]=filterstrandbiasprop(pair[0]["sb"],sdict["strandbiastumorq"])  # TUMOR
            vcf_reader.filters[19]=["strandbiasprop","Strand bias exclusion proportion: "+str(pdict["strandbiasprop"])+" (>"+str(round(sdict["strandbiastumorq"],3))+")"]
        if key == "strandbias":
            fdict["strandbias"]=filterstrandbias(pair[0]["sb"],pdict["strandbias"])  # TUMOR
            vcf_reader.filters[20]=["strandbias","Maximum strand bias in tumor: "+str(pdict["strandbias"])]
        if key == "strandbiassimple":
            fdict["strandbiassimple"]=filterstrandbiassimple(pair[0]["altsb"],pair[0]["allsb"],pdict["strandbiassimple"])  # TUMOR
            vcf_reader.filters[20]=["strandbiassimple","Maximum strand bias in ALT reads in tumor: "+str(pdict["strandbiassimple"])]
        if key == "repeats":
            fdict["repeats"]=filterrepeats(pair[0]["CHR"],pair[0]["POS"],pdict["repeats"],repeatlist,referencegenome)  # TUMOR
            vcf_reader.filters[21]=["repeats","Maximum length of 1/2/3/4mer repeats around the variant position: "+str(int(pdict["repeats"]))]
        if key == "maxaltcount":
            fdict["maxaltcount"]=filtermaxaltcount(pair[1]["altcount"],pdict["maxaltcount"])  # NORMAL
            vcf_reader.filters[22]=["maxaltcount","Maximum number of ALT reads in normal: "+str(int(pdict["maxaltcount"]))]
        if key == "maxvafnormal":
            fdict["maxvafnormal"]=filtermaxvaf(pair[1]["vaf"],pdict["maxvafnormal"])  # NORMAL
            vcf_reader.filters[23]=["maxvafnormal","Maximum VAF in normal: "+str(pdict["maxvafnormal"])]
        if key == "maxbadorient":
            fdict["maxrefbadorientnormal"]=filtermaxbadorient(pair[1]["refbadorientationprop"],pdict["maxbadorient"])  # NORMAL
            vcf_reader.filters[24]=["maxrefbadorientnormal","Maximum proportion of inversion orientation reads in normal: "+str(pdict["maxbadorient"])]
        if key == "minmapqualitydifference":
            fdict["minmapqualitydifferencent"]=filterminmapqualitydifference(pair[0]["medianmapqalt"],pair[1]["medianmapqref"],pdict["minmapqualitydifference"])  # TUMOR AND NORMAL
            vcf_reader.filters[25]=["minmapqualitydifferencent","Maximum difference between median mapping quality of ALT reads in tumor and REF reads in normal: "+str(pdict["minmapqualitydifference"])]
    return(fdict)


## Determines filtering result
def judgement(fdict):
    result=[]
    for key in fdict.keys():
        if(not fdict[key]):
            result.append(key)
    if(len(result)==0):
        result=["PASS"]
    return(result)
