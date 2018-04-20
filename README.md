# FiNGS

## Quickstart guide

Full installation instructions and documentation are below, but the fastest way to get up and running is:
```
git clone https://github.com/cpwardell/FiNGS.git
cd FiNGS
chmod +x install_FiNGS_dependencies.sh
sudo ./install_FiNGS_dependencies.sh
python3 FiNGS.py -n /path/to/normal.bam -t /path/to/tumor.bam -v /path/to/somaticvariants.vcf
```

## Suggested usage
+ Use default filtering thresholds
+ Use every available processor
+ Only consider variants with a PASS filter value from the variant caller
+ Only emit variants that PASS all filters  

```
python3 FiNGS.py -n /path/to/normal.bam -t /path/to/tumor.bam -v /path/to/somaticvariants.vcf -j -1 --PASSonlyin --PASSonlyout
```

## Introduction


## Availability
FiNGS is released under the [Apache License, Version 2.0](https://www.apache.org/licenses/LICENSE-2.0.html). The latest source code is [freely available at GitHub](https://github.com/cpwardell/FiNGS).

## Installation
FiNGS depends on Python3, R, and several modules/packages for each.  The included script `FiNGS/install_FiNGS_dependencies.sh` should install all dependencies, but this can be manually performed as below.

On Debian-based systems such as Ubuntu, you can install Python3 and pip, the recommended tool for installing Python packages using the following commands:
```
apt-get -y install python3
apt-get -y install python3-pip
```

Now Python3 and pip are installed, you may install the required Python packages:
```
pip3 install PyVCF
pip3 install pysam
pip3 install editdistance
pip3 install scipy
pip3 install joblib
```

To install R and some components required by Bioconductor, use the following commands:
```
apt-get -y install r-base
apt-get -y install libcurl4-openssl-dev
apt-get -y install libxml2-dev
apt-get -y install zlib1g-dev
```

From within R, run the following commands to install the required R package:
```
source("https://bioconductor.org/biocLite.R")
biocLite("VariantAnnotation")
```

## Basic usage and outputs
In the simplest case, provide FiNGS with paths to a VCF and the BAM files used to produce it.
```
python3 FiNGS.py -n /path/to/normal.bam -t /path/to/tumor.bam -v /path/to/somaticvariants.vcf
```

FiNGS will create a directory called "results" containing the following files, where "date" is the time of execution expressed as a string e.g. "20180420155129" 2018/04/20 15:51:59 local time.

+ `filtered.results.date.vcf`  
A VCF containing the filtering results. Descriptions of the filters used and their threshold values are stored in the header, and PASS/fail status stored in the FILTER column of each record
+ `filtered.results.date.txt`  
A text file giving the raw data for every variant and TRUE/FALSE values for each filter. A dictionary explaining the contents of each column is below
+ `filter.plots.date.vcf`  
Plots for every filter applied in a single PDF. Data is plotted as a histogram with a kernel density plot overlaid, and a dashed vertical line demarcating the threshold used 
+ `log.txt`  
The log file for the run. Contains the command line arguments used and a complete log of the run
+ `normal.combined.txt`  
All intermediate chunks of the normal BAM file analysis are combined into this file. This and the matching combined tumor file are the input for the R filtering step
+ `tumor.combined.txt`  
All intermediate chunks of the tumor BAM file analysis are combined into this file.  This and the matching combined normal file are the input for the R filtering step

## Advanced usage

+ **-v** Required; path to the somatic variant VCF from any variant caller  
+ **-t** Required; path to the tumor BAM file used to produce the VCF  
+ **-n** Required; path to the normal BAM file used to produce the VCF  
+ **-d** Optional; path to output directory.  Default is to create a directory called "results" in the current working directory
+ **-p** Optional; path to file specifying filtering parameters.  Details on filters and default values is provided below
Default is a file located at `FiNGS/R/filter_parameters.txt`.  To alter these parameters, copy and edit this file, and provide a path to it
+ **-c** Optional; chunk size, the number of records to process per chunk. Default is 100
+ **-m** Optional; maximum read depth. Reads beyond this depth will be ignored and a warning emitted to the log file.Default is 1000 
+ **-j** Optional; number of processors to use. Use -1 to use all available. Default is 1 
+ **--logging** Optional; change logging level. Default is INFO, can be DEBUG for more detail or NOTSET for silent
+ **--donotcleanup** Optional; do not delete intermediate files  
+ **--overwrite** Optional; allow overwriting of existing results if rerunning  
+ **--PASSonlyin** Optional; only consider variants with a PASS in the filter field of the input VCF  
+ **--PASSonlyout** Optional; only write variants that PASS all filters to the output VCF

## Getting help
Run FiNGS with no additional arguments to get the help file.  If there's something not adddressed here, or if you need further help, raise an issue on GitHub or find me online.

## Citing FiNGS
A paper is being prepared for submission shortly and will be referenced here when available.

## Description of filters
FiNGS assesses variants using 22 filters.  Below are descriptions of each filter and the default values used.
 
1. **Minimum depth in tumor**.  Default value is 10x.  
This filter excludes variants that do not have adequate sequencing depth in the tumor sample.
 
2. **Minimum depth in normal**.  Default value is 10x.  
This filter excludes variants that do not have adequate sequencing depth in the normal sample.

3. **Maximum depth in normal**.  Default value is 1000x.  
This filter excludes variants that are unusually deep.  Some variants have apparent depths of thousands of reads, particularly around the centromeres, telomeres and other repeat-rich regions.

4. **Maximum depth in tumor**.  Default value is 1000x.  
This filter excludes variants that are unusually deep.  Some variants have apparent depths of thousands of reads, particularly around the centromeres, telomeres and other repeat-rich regions.

5. **Minimum number of ALT reads in tumor**.  Default value is 3.  
At low depths, a variant might have VAF above the default threshold of 5%, but very few reads supporting it.  This default value ensures that such a 5% variant would have a minimum of 60x depth.

6. **Minimum median base quality of ALT reads in tumor**.  Default value is 30.  
Excludes variants where the supporting reads in the tumor have low base quality.  Taking the median of all reads would allow non-variant reads to weight the median base quality.

7. **Minimum median base quality of REF reads in tumor**.  Default value is 30.  
Excludes variants where REF reads in the tumor have low base quality.  Variants with no REF reads (i.e. homozygous ALT) pass this filter.

8. **Minimum median base quality of REF reads in normal**.  Default value is 30.  
Excludes variants where REF reads in the normal have low base quality.

9. **Maximum proportion of zero mapping quality reads in tumor**.  Default value is 0.05.  
A high proportion of zero mapping quality reads suggests a region with poorly mappability.

10. **Maximum proportion of zero mapping quality reads in normal**.  Default value is 0.05.  
A high proportion of zero mapping quality reads suggests a region with poorly mappability.

11. **Strand bias exclusion proportion**.  Default value is 0.1.  
Strand bias is defined below [(Guo et al., 2012)](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-666) where a,c represent the forward and reverse strand allele counts of REF reads and b,d represent the forward and reverse strand allele counts of ALT reads.  The topmost proportion of biased variants is removed.
|b/(a+b)-d/(c+d)|/ ((b+d)/(a+b+c+d))

12. **Minimum median mapping quality of ALT reads in tumor**.  Default value is 50.  
Excludes variants where the supporting reads in the tumor have low mapping quality.  Taking the median of all reads would allow non-variant reads to weight the median mapping quality.

13. **Maximum difference between median mapping quality of ALT reads in tumor and REF reads in normal**.  Default value is 5.  
Excludes variants with a large difference in mapping quality between ALT reads in the tumor and REF reads in the matched normal sample.  Mismatches in the ALT reads will lower mapping quality, but this difference should be limited.

14. **Maximum median shortest distance to either aligned end in tumor**.  Default value is 10.  
Excludes variants that consistently occur only the ends of reads.

15. **Minimum MAD of ALT position in tumor**.  Default value is 3.  
Excludes variants that consistently occur in the same position within reads.

16. **Maximum edit distance of ALT reads in tumor**.  Default value is 4.  
Limits the number of mismatches possible within variant supporting reads in the tumor.  Default value allows for multiple mismatches, allowing SNPs and small indels to occur in addition to the variant of interest.

17. **Maximum edit distance of REF reads in tumor**.  Default value is 3.  
Limits the number of mismatches possible within reference reads in the normal.  Default value allows for multiple mismatches, allowing SNPs and small indels to occur in addition to the variant of interest.  Multiple mismatches in the reference reads in the normal would be indicative of a low quality region in the reference genome.

18. **Maximum VAF in normal**.  Default value is 0.03.  
Excludes variants with too many supporting reads in the normal sample.  Some support is allowed, which accounts for sequencing errors and tumor cross-contamination.

19. **Minimum VAF in tumor**.  Default value is 0.05.  
Excludes variants with too few supporting reads in the tumor sample.  It is difficult to differentiate between genuine low-frequency, subclonal variants, and sequencing errors or contamination.

20. **Maximum OAF in tumor**.  Default value is 0.04.  
Excludes variants where multiple alleles are present.  Other Allele Frequency (OAF) is the proportion of reads that support an allele other than the reference or alternate alleles.

21. **Maximum proportion of secondary alignments in tumor**.  Default value is 0.05.  
Excludes variants with a high proportion of secondary alignments.  Secondary alignments are indicative of poor mappability.

22. **Maximum proportion of inversion orientation reads in normal**.  Default value is 0.2.  
Excludes variants with a high proportion of reads in the normal sample that have inverted orientation.  This suggests a region of poor mappability.


## Dictionary of values reported in the metrics files
The `filtered.results.date.txt` output file is the `normal.combined.txt` and `tumor.combined.txt` files joined by column, with additional columns for each filter containing TRUE/FALSE for whether or not that variant passed that filter. Each row is a single variant.  
Column names are delimited by full stops (periods). The first half defines the data type and the second the source. For example, depth.tumor and depth.normal contain the depth of a variant for the tumor and normal respectively.


Column | Description | Example
--- | --- | ---
**UID** | Thingsly | blah
**CHR** | Thingsly | blah
**POS** | Thingsly | blah
**REF** | Thingsly | blah
**ALT** | Thingsly | blah


refcount
altcount
varianttype
depth
vaf
raf
oaf
medianbaseq
medianbaseqref
medianbaseqalt
medianmapq
medianmapqref
medianmapqalt
zeros
zerospersite
softreadlengthsrefmean
softreadlengthsaltmean
goodoffsetproportion
distancetoend1median
mad1
distancetoend2median
mad2
distancetoend1medianref
madref1
distancetoend2medianref
madref2
distancetoend1medianalt
madalt1
distancetoend2medianalt
madalt2
shortestdistancetoendmedian
madaltshort
sb
gsb
fishp
FoxoG
refld
altld
refsecondprop
altsecondprop
refbadorientationprop
altbadorientationprop
refmatecontigcount
altmatecontigcount




