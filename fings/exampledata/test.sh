#!/bin/bash

## Examples of how to run FiNGS on the example data, with number of PASSing variants for each

## Recommended settings
fings -n normal.bam -t tumor.bam -v s2.raw.vcf --PASSonlyin --PASSonlyout # 16 PASS

## ICGC mode
#fings -n normal.bam -t tumor.bam -v s2.raw.vcf -r /path/to/reference/genome.fa --PASSonlyin --PASSonlyout --ICGC # 23 PASS

## Other options
#fings -n normal.bam -t tumor.bam -v s2.raw.vcf --PASSonlyin # 16 PASS, 38 total
#fings -n normal.bam -t tumor.bam -v s2.raw.vcf --PASSonlyout # 20 PASS
#fings -n normal.bam -t tumor.bam -v s2.raw.vcf # 20 PASS, 1207 total


