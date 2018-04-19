# FiNGS

## Quickstart guide

Full installation instructions and documentation are below, but the fastest way to get up and running is:

	git clone https://github.com/cpwardell/FiNGS.git
	cd FiNGS
	python3 FiNGS.py -n /path/to/normal.bam -t /path/to/tumor.bam -v /path/to/somaticvariants.vcf

## Introduction


## Availability
FiNGS is released under the [Apache License, Version 2.0][1]. The latest source code is [freely available at github][2].

## Installation
FiNGS depends on Python3, R, and several modules/packages for each.

On Debian-based systems, you can install Python3 and pip, the recommended tool for installing Python packages using the following commands:

apt-get -y install python3
apt-get -y install python3-pip

Now Python3 and pip are installed, you may install the required packages:

pip3 install PyVCF
pip3 install pysam
pip3 install editdistance
pip3 install scipy
pip3 install joblib

To install R and some components required by Bioconductor, use the following commands:
apt-get -y install r-base
apt-get -y install libcurl4-openssl-dev
apt-get -y install libxml2-dev
apt-get -y install zlib1g-dev

From within R, run the following commands to install the required R package:
source("https://bioconductor.org/biocLite.R")
biocLite("VariantAnnotation")

## Docker usage...


## Description of filters
At release FiNGS assesses variants using 22 filters.  Below are descriptions of each filter, the default values used and a plot to illustrate how the chosen thresholds perform.  The GIAB WGS data is used in these plots.  Each is a histogram showing the distribution of the data, with a three kernel density plots fitted over the top; a black line for all the data, and green and red lines to show the distribution of true and false positives respectively.  The vertical dashed line represents the threshold value.


 
Filter 1: Minimum depth in tumor.  Default value is 10x
This filter excludes variants that do not have adequate sequencing depth in the tumor sample.
 
Filter 2: Minimum depth in normal.  Default value is 10x
This filter excludes variants that do not have adequate sequencing depth in the normal sample.

Filter 3: Maximum depth in normal.  Default value is 1000x
This filter excludes variants that are unusually deep.  Some variants have apparent depths of thousands of reads, particularly around the centromeres, telomeres and other repeat-rich regions.

Filter 4: Maximum depth in tumor.  Default value is 1000x
This filter excludes variants that are unusually deep.  Some variants have apparent depths of thousands of reads, particularly around the centromeres, telomeres and other repeat-rich regions.

Filter 5: Minimum number of ALT reads in tumor.  Default value is 3
At low depths, a variant might have VAF above the default threshold of 5%, but very few reads supporting it.  This default value ensures that such a 5% variant would have a minimum of 60x depth.

Filter 6: Minimum median base quality of ALT reads in tumor.  Default value is 30
Excludes variants where the supporting reads in the tumor have low base quality.  Taking the median of all reads would allow non-variant reads to weight the median base quality.

Filter 7: Minimum median base quality of REF reads in tumor.  Default value is 30
Excludes variants where REF reads in the tumor have low base quality.  Variants with no REF reads (i.e. homozygous ALT) pass this filter.

Filter 8: Minimum median base quality of REF reads in normal.  Default value is 30
Excludes variants where REF reads in the normal have low base quality.

Filter 9: Maximum proportion of zero mapping quality reads in tumor.  Default value is 0.05
A high proportion of zero mapping quality reads suggests a region with poorly mappability.

Filter 10: Maximum proportion of zero mapping quality reads in normal.  Default value is 0.05
A high proportion of zero mapping quality reads suggests a region with poorly mappability.

Filter 11: Strand bias exclusion proportion.  Default value is 0.1
Strand bias is defined below (Guo et al., 2012) where a,c represent the forward and reverse strand allele counts of REF reads and b,d represent the forward and reverse strand allele counts of ALT reads.  The topmost proportion of biased variants is removed.
|b/(a+b)-d/(c+d)|/ ((b+d)/(a+b+c+d))

Filter 12: Minimum median mapping quality of ALT reads in tumor.  Default value is 50
Excludes variants where the supporting reads in the tumor have low mapping quality.  Taking the median of all reads would allow non-variant reads to weight the median mapping quality.

Filter 13: Maximum difference between median mapping quality of ALT reads in tumor and REF reads in normal.  Default value is 5
Excludes variants with a large difference in mapping quality between ALT reads in the tumor and REF reads in the matched normal sample.  Mismatches in the ALT reads will lower mapping quality, but this difference should be limited.

Filter 14: Maximum median shortest distance to either aligned end in tumor.  Default value is 10
Excludes variants that consistently occur only the ends of reads.

Filter 15: Minimum MAD of ALT position in tumor.  Default value is 3
Excludes variants that consistently occur in the same position within reads.

Filter 16. Maximum edit distance of ALT reads in tumor.  Default value is 4
Limits the number of mismatches possible within variant supporting reads in the tumor.  Default value allows for multiple mismatches, allowing SNPs and small indels to occur in addition to the variant of interest.

Filter 17: Maximum edit distance of REF reads in tumor.  Default value is 3
Limits the number of mismatches possible within reference reads in the normal.  Default value allows for multiple mismatches, allowing SNPs and small indels to occur in addition to the variant of interest.  Multiple mismatches in the reference reads in the normal would be indicative of a low quality region in the reference genome.

Filter 18: Maximum VAF in normal.  Default value is 0.03
Excludes variants with too many supporting reads in the normal sample.  Some support is allowed, which accounts for sequencing errors and tumor cross-contamination.

Filter 19: Minimum VAF in tumor.  Default value is 0.05
Excludes variants with too few supporting reads in the tumor sample.  It is difficult to differentiate between genuine low-frequency, subclonal variants, and sequencing errors or contamination.

Filter 20: Maximum OAF in tumor.  Default value is 0.04
Excludes variants where multiple alleles are present.  Other Allele Frequency (OAF) is the proportion of reads that support an allele other than the reference or alternate alleles.

Filter 21: Maximum proportion of secondary alignments in tumor.  Default value is 0.05
Excludes variants with a high proportion of secondary alignments.  Secondary alignments are indicative of poor mappability.

Filter 22: Maximum proportion of inversion orientation reads in normal.  Default value is 0.2
Excludes variants with a high proportion of reads in the normal sample that have inverted orientation.  This suggests a region of poor mappability.




## Dictionary of values reported in the metrics files


[1]: https://www.apache.org/licenses/LICENSE-2.0.html
[2]: https://github.com/cpwardell/FiNGS