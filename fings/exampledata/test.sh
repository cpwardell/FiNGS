#!/bin/bash

## Examples of how to run FiNGS on the example data, with number of PASSing variants for each

## Recommended settings
fings -n normal.bam -t tumor.bam -v s2.raw.vcf --PASSonlyin --PASSonlyout # 19 PASS, 19 total

## ICGC mode
#fings -n normal.bam -t tumor.bam -v s2.raw.vcf -r /path/to/reference/genome.fa --PASSonlyin --PASSonlyout --ICGC # 23 PASS

## Other options
#fings -n normal.bam -t tumor.bam -v s2.raw.vcf --PASSonlyin # 19 PASS, 38 total
#fings -n normal.bam -t tumor.bam -v s2.raw.vcf --PASSonlyout # 20 PASS, 20 total
#fings -n normal.bam -t tumor.bam -v s2.raw.vcf # 20 PASS, 1207 total


