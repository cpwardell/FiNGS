<p align="center">
  <img src="/www/fings_logo_20190809_color_banner.png">
</p>

# <span style="color:#8a0707">F</span></span><span style="color:#b2b2b2">***i***</span>NGS: <span style="color:#8a0707">***F***</span><span style="color:#b2b2b2">***i***</span>lters for ***N***ext ***G***eneration ***S***equencing

## Key features
- **Filters SNVs from any variant caller to remove false positives**
- **Calculates metrics based on BAM files and provides filtering not possible with other tools**
- **Fully user-configurable filtering (including which filters to use and their thresholds)**
- **Option to use filters identical to ICGC recommendations**

## Introduction
Somatic variant callers compare matched pairs of tumor-normal samples to produce variant calls. The results can be extremely rich in false positives due to confounding factors such as the purity of the samples, artifacts introduced by sequencing chemistry, the alignment algorithm and the incomplete and repetitive nature of reference genomes.  
It has become common practice to attempt to ameliorate these effects using a variety of filtering techniques, including taking the intersect of results from multiple variant callers and employing some post-calling filtering. This ad-hoc filtering varies greatly between laboratories.  Attempts have been made to standardize the methodology of this filtering, with recommendations produced by the International Cancer Genome Consortium (ICGC) [(Alioto et al., 2015)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4682041/).  
We have developed Filters for Next Generation Sequencing (FiNGS), software written specifically to address these filtering issues. FiNGS can implement the ICGC filtering standards and has filters and thresholds that are fully configurable, which substantially increases the precision of results and provides high quality variants for further analysis. 

## Availability
FiNGS is open source and released under the [Apache License, Version 2.0](https://www.apache.org/licenses/LICENSE-2.0.html). The latest source code is [freely available at GitHub](https://github.com/cpwardell/FiNGS).

## Dependencies
```Python 3``` and these Python packages:  
```
pyvcf  
pysam  
numpy  
scipy  
pandas  
joblib  
seaborn  
statsmodels  
editdistance
```

## Quickstart guide, with example data and test
You have a number of options for installing and running FiNGS: 
1. Download using Bioconda (preferred)
2. Download directly from GitHub
3. Download using pip3
4. Download using Docker/Singularity

### Bioconda installation
FiNGS is on Bioconda here: https://bioconda.github.io/recipes/fings/README.html  
This method is preferred because it gives the cleanest environment and installs all dependencies, including system ones required to build dependencies like the pysam package.  You can follow the set-up instructions for Bioconda from [here](https://bioconda.github.io/user/install.html); assuming you already have a Conda installation, set up the Bioconda channels like so:

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Create a new conda environment for FiNGS and activate it
```
conda create -n fings python=3.7
conda activate fings
```

Now install the package (dependencies are installed automatically)
```
conda install fings
```

Locate the FiNGS package, navigate to the example data and run the script to confirm your installation is correct.
```
## Navigate to this directory /path/to/conda/envs/fings/lib/python3.7/site-packages/fings/exampledata/
./test.sh
```

### GitHub installation
The advantage of a GitHub installation is that you know exactly where your FiNGS code is; we still advise managing the dependencies listed above using Conda.
Clone the git repository and run the included example script in the `exampledata` directory.
```
git clone https://github.com/cpwardell/FiNGS.git
cd FiNGS/fings/exampledata
./test.sh
```

### python3-pip installation
Assuming you have Python3 and pip3, you can install using pip because FiNGS lives on PyPi here: https://pypi.org/project/fings/
```
pip3 install fings
```
Locate the FiNGS package, navigate to the example data and run the script to confirm your installation is correct.
```
## Use pip3 to tell you where FiNGS has been installed
pip3 show fings
cd /go/to/fings/exampledata/
./test.sh
```

### Docker installation and usage 
This guide assumes you have Docker installed and have some basic knowledge on how to use it. The Dockerfile builds an image based on the official Miniconda3 image and pulls the current Bioconda version of FiNGS, _not_ the current GitHub version.

You can either build your own image or pull it from Docker Hub.

#### Pulling from Docker Hub

FiNGS Docker Hub page is here: https://hub.docker.com/r/cpwardell/fings

```
docker pull cpwardell/fings
```
#### Building local image
You need to get a copy of the Dockerfile in this repository; below we use "wget" on Linux to download it, but you could just as easily 
copy and paste the link in your web browser and "right click/save as" the file. The Docker build command works identically in both Bash on Linux and PowerShell on Windows
and assumes that you're in the same directory as the dockerfile named "Dockerfile".

```
# Download the Dockerfile from this address:
wget https://raw.githubusercontent.com/cpwardell/FiNGS/master/Dockerfile
# Build the image and call it "fings" (lowercase)
docker build -t fings .
```

#### Suggested Docker usage
The default entrypoint for the image is launching FiNGS, so you can treat it much like the command-line version. However, you **MUST** ensure that the data you wish to use is somewhere the Docker container can access (by mounting directories using the -v argument of Docker) and that you either specify the working directory (-w argument of Docker) or specify the output directory (-d option of FiNGS).  If you fail to do this, the process will either fail or write to a directory internal to the Docker image. Note that the `-u` argument ensures that files created by the Docker container will be owned by the user invoking the process.
```
## Simple example assumes user has downloaded the included exampledata directory
docker run -it -v $PWD:/local -w /local -u $UID:$UID fings -n normal.bam -t tumor.bam -v s2.raw.vcf --PASSonlyin --PASSonlyout

## Longer example with more complex setup
docker run -v /path/to/tumorbamdir:/tumorbamdir -v /path/to/normalbamdir:/normalbamdir -v /path/to/vcfdir:/vcfdir -v $PWD:/local -w /local -u $UID:$UID fings -n /normalbamdir/normal.bam -t /tumorbamdir/tumor.bam -v /vcfdir/somatic.vcf --PASSonlyin --PASSonlyout
```

## Suggested usage
+ Use default filtering thresholds (either our filters or ICGC filters)
+ Use every available processor
+ Only consider variants with a PASS filter value from the variant caller
+ Only emit variants that PASS all filters  

```
python3 /path/to/FiNGS.py -n /path/to/normal.bam -t /path/to/tumor.bam -v /path/to/somaticvariants.vcf --PASSonlyin --PASSonlyout
```

+ ICGC mode:

```
python3 /path/to/FiNGS.py -n /path/to/normal.bam -t /path/to/tumor.bam -v /path/to/somaticvariants.vcf -r /path/to/reference/genome.fasta --PASSonlyin --PASSonlyout --ICGC
```

FiNGS will create a directory called "results" containing the following files:

+ `inputvcf.filtered.vcf`  
A VCF containing the filtered results. Descriptions of the filters used and their threshold values are stored in the header, and PASS/fail status stored in the FILTER column of each record
+ `plots.pdf`  
Plots for every filter applied in a single PDF. The first page shows a table of the PASS/Fail counts for each filter. Subsequent pages show kernel density plots of the data used, with a dashed vertical line demarcating the threshold used; the red region failed the filter, the green region passed
+ `log.txt`  
The log file for the run. Contains the command line arguments used and a complete log of the run
+ `filterresults.txt.gz`  
A gzipped text file giving the results of each filter for every variant 
+ `tumor.combined.txt.gz` and `normal.combined.txt.gz`  
All metrics collected and used for filtering are stored in these gzipped text files. A dictionary explaining the contents of each column is below
+ `summarystats.txt.gz`
A gzipped text file containing summary stats that may be used for filtering

## Further notes
The following arguments and flags are available:  
+ **-v** Required; path to the somatic variant VCF from any variant caller
+ **-t** Required; path to the tumor BAM file used to produce the VCF
+ **-n** Required; path to the normal BAM file used to produce the VCF
+ **-r** Optional; absolute path to faidx indexed reference genome; required if using 'repeats' filter 
+ **-d** Optional; path to output directory.  Default is to create a directory called "results" in the current working directory
+ **-p** Optional; path to file specifying filtering parameters. Details on filters and default values is provided below. Default is a file located at `FiNGS/filter_parameters.txt`
+ **-c** Optional; chunk size, the number of variants to process per chunk. Default is 100
+ **-m** Optional; maximum read depth. Reads beyond this depth will be ignored and a warning emitted to the log file. Default is 1000 
+ **-j** Optional; number of processors to use. -1 uses all available. Default is -1 
+ **--ICGC** Optional; use filters identical to those recommended by the ICGC (Alioto et al, 2015). File is located at `FiNGS/icgc_filter_parameters.txt`
+ **--logging** Optional; change logging level. Default is INFO, can be DEBUG for more detail or NOTSET for silent
+ **--overwrite** Optional; allow overwriting of existing results if rerunning  
+ **--PASSonlyin** Optional; only consider variants with a PASS in the filter field of the input VCF  
+ **--PASSonlyout** Optional; only write variants that PASS all filters to the output VCF

## Getting help
Run FiNGS with no additional arguments to get the help file. If there's something not adddressed here, or if you need further help, raise an issue on GitHub or find me online.

## Citing FiNGS
A paper is being prepared for submission shortly and will be referenced here when available.

## Description of filters
FiNGS assesses variants using any combination of these possible filters. Below is a table describing them, their default thresholds and ICGC thresholds.  NA values mean that the filter is not employed in eithe the default or ICGC mode.  **Users can create their own tab-delimited parameter text file using *any* combination of filters and thresholds, and pass it in using the -p argument.**

Filter name | Description | Default value | ICGC value
--- | --- | --- | ---
**minaltcount** | Minimum number of ALT reads in tumor | 3 | 4
**minbasequality** | Minimum median base quality (separate filters for ALT reads in tumor, REF reads in tumor and REF reads in normal) | 30 | 30
**minmapquality** | Minimum median mapping quality of ALT reads in tumor | 50 | 40
**minmapqualitydifference** | Maximum difference between median mapping quality of ALT reads in tumor and REF reads in normal | 5 | 5
**enddistance** | Maximum median shortest distance to either aligned end in tumor | 10 | 10
**enddistancemad** | Minimum MAD of ALT position in tumor | 3 | 3
**zeroproportion** | Maximum proportion of zero mapping quality reads in tumor and normal | 0.05 | 0.1
**minimumdepth** | Minimum depth in tumor and normal | 10 | NA
**maximumdepth** | Maximum depth in tumor and normal | 1000 | NA
**minvaftumor** | Minimum VAF in tumor | 0.05 | NA
**maxvafnormal** | Maximum VAF in normal | 0.03 | NA
**maxoaftumor** | Maximum OAF in tumor | 0.04 | NA
**foxog** | FoxoG artifact proportion (see note below) | 0.9 | NA
**editdistance** | Maximum edit distance of ALT reads in tumor (maximum edit distance of REF reads in tumor is 1 less than this value) | 4 | NA
**maxsecondtumor** | Maximum proportion of secondary alignments in tumor | 0.05 | NA
**maxbadorient** | Maximum proportion of inversion orientation reads in normal | 0.2 | NA
**strandbiasprop** | Strand bias exclusion proportion (see note below) | 0.1 | NA
**strandbiassimple** | Maximum strand bias (see note below) | NA | 0.02
**maxaltcount** | Maximum number of ALT reads in normal | NA | 1
**snvcluster50** | Maximum number of mutations within 50 bp (see note below) | NA | 2
**snvcluster100** | Maximum number of mutations within 100 bp (see note below) | NA | 4
**repeats** | Maximum length of 1/2/3/4mer repeats around the variant position (see note below) | NA | 12

+ **Note on _foxog_ filter**  
C>A|G>T variants can be oxidation artifacts introduced during library preparation [(Costello et al, 2013)](http://doi.org/10.1093/nar/gks1443).
These "OxoG" artifacts have a telltale read orientation, with the majority of ALT reads in the artifact orientation. All C>A|G>T variants are classified as being part of two binomial distributions, one centered at 0.5 (50% artifact orientation, 50% non-artifact orientation) and the other at the filter value (defautl is 0.9, which is 90% artifact orientation reads).  C>A|G>T variants classified as OxoG are removed.

+ **Note on _strandbiasprop_ filter**  
Strand bias is defined below [(Guo et al., 2012)](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-666) where a,c represent the forward and reverse strand allele counts of REF reads and b,d represent the forward and reverse strand allele counts of ALT reads.  The topmost proportion of biased variants is removed (e.g. suggested value is 0.1, leading to the top 10% variants with the highest strand bias being removed. Note that this is not the GATK strand bias, which is calculated differently. Note that this filter is available, but is *not* implemented in either of the default settings.  
Strand bias: |b/(a+b)-d/(c+d)|/ ((b+d)/(a+b+c+d))  
GATK strand bias: max(((b/(a+b))&ast;(c/(c+d)))/((a+c)/(a+b+c+d)),((d/(c+d))&ast;(a/(a+b)))/((a+c)/(a+b+c+d)))

+ **Note on _strandbiassimple_ filter**  
Maximum strand bias. This is defined very simply as the minimum proportion of reads in either direction. e.g. if there were 100 reads and only 1 were forward, strand bias would be 0.01.
Strand bias: min(forward/(forward+reverse),reverse/(reverse+forward)) or min((a+b)/(a+b+c+d),(c+d)/(a+b+c+d))

+ **Note on _SNVcluster50_ and _SNVcluster100_ filters**  
Maximum number of SNVs in the input VCF in a 50 or 100 bp window centered on the current SNV.  SNVs must have at least 2 reads supporting them and a VAF>=5%.

+ **Note on _repeats_ filter**  
Maximum length of 1/2/3/4mers surrounding the SNV. Lengths must be factors of repeats e.g. n=8 would only consider 1/2/4mers because 3 is not a factor of 8.  Repeated regions frequently result in false positive variants.  When using this filter, a reference genome *must* be supplied so FiNGS can find the flanking sequences.

## Dictionary of values reported in the metrics files
The gzipped `tumor.combined.txt.gz` and `normal.combined.txt.gz` output files contain all metrics calculated for every input variant. Each row is a single variant. Non-integer values are rounded to 3 decimal places. Not all of the values reported are used for filtering.

Column | Description | Example
--- | --- | ---
**UID** | Unique Identifier for variant | 1:931362:G:A
**CHR** | Chromosome | 1
**POS** | Position on chromsome | 931362
**REF** | Reference allele | G
**ALT** | Alternate allele | A
**refcount** | Count of REF alleles | 100
**altcount** | Count of ALT alleles | 19
**varianttype** | SNV or INDEL | SNV
**depth** | Depth of all reads | 122
**vaf** | Variant Allele Frequency (altcount/depth) | 0.156
**raf** | Reference Allele Frequency (refcount/depth) | 0.82
**oaf** | Other Allele Frequency ((depth-altcount-refcount)/depth) | 0.025
**medianbaseq** | Median base quality (all reads) | 32
**medianbaseqref** | Median base quality (REF reads only) | 34.5
**medianbaseqalt** | Median base quality (ALT reads only) | 32
**medianmapq** | Median mapping quality (all reads) | 60
**medianmapqref** | Median mapping quality (REF reads only) | 60
**medianmapqalt** | Median mapping quality (ALT reads only) | 60
**zeros** | Total number of zero mapping quality reads | 0
**zerospersite** | Proportion of reads that have zero mapping quality | 0
**softreadlengthsrefmean** | Mean length of REF reads after soft clipping | 147.68
**softreadlengthsaltmean** | Mean length of ALT reads after soft clipping | 151
**goodoffsetproportion** | Proportion of variants that occur within the first 2/3rds of the mean read length | 0.664
**distancetoend1median** | Median distance to lefthand soft-clipped read end (all reads) | 74.5
**mad1** | Median absolute deviation of distancetoend1median | 34
**distancetoend2median** | Median distance to righthand soft-clipped read end (all reads) | 70
**mad2** | Median absolute deviation of distancetoend2median | 34.5
**distancetoend1medianref** | Median distance to lefthand soft-clipped read end (REF reads only) | 76
**madref1** | Median absolute deviation of distancetoend1medianref | 31.5
**distancetoend2medianref** | Median distance to righthand soft-clipped read end (REF reads only) | 64
**madref2** | Median absolute deviation of distancetoend2medianref | 29
**distancetoend1medianalt** | Median distance to lefthand soft-clipped read end (ALT reads only) | 65
**madalt1** | Median absolute deviation of distancetoend1medianalt | 25
**distancetoend2medianalt** | Median distance to righthand soft-clipped read end (ALT reads only) | 85
**madalt2** | Median absolute deviation of distancetoend2medianalt | 25
**shortestdistancetoendmedian** | Median of shortest distance of distancetoend1alt and distancetoend2alt | 42
**madaltshort** | Median absolute deviation of shortestdistancetoendmedian | 19
**sb** | Strand bias, see definition above | 0.22
**gsb** | GATK strand bias, see definition above | 0.185
**fishp** | P value for Fisher's exact test for strand bias | 0.614
**FR** | Count of forward reads supporting REF allele | 36
**FA** | Count of forward reads supporting ALT allele | 8
**RR** | Count of reverse reads supporting REF allele | 64
**RA** | Count of reverse reads supporting ALT allele | 11
**altsb** | Simple strand bias of ALT reads | 0.421
**refsb** | Simple strand bias of REF reads | 0.36
**allsb** | Simple strand bias of all reads | 0.37
**F1R2** | Count of reads in FoxoG orientation 1 | 8 
**F2R1** | Count of reads in FoxoG orientation 2 | 11
**FoxoG** | Oxoguanine artifact orientation proportion, only relevant for for C>A or G>T mutations, defined in [Costello et al, 2013](https://academic.oup.com/nar/article/41/6/e67/2902364) | 0.421
**refld** | Edit distance for REF reads | 2
**altld** | Edit distance for ALT reads | 3
**refsecondprop** | Proportion of REF reads that have secondary alignments | 0
**altsecondprop** | Proportion of ALT reads that have secondary alignments | 0
**refbadorientationprop** | Proportion of REF reads with an inverted orientation | 0
**altbadorientationprop** | Proportion of ALT reads with an inverted orientation | 0
**refmatecontigcount** | Number of contigs seen in REF reads | 1
**altmatecontigcount** | Number of contigs seen in ALT reads | 1
**sixtypes** | Types of SNV (of the six possible types) | C>T/G>A

