#!/bin/bash

## Examples of how to run FiNGS on the example data, with number of PASSing variants for each

## Recommended settings
python3 ../FiNGS.py -n normal.bam -t tumor.bam -v s2.raw.vcf --PASSonlyin --PASSonlyout # 16 PASS

## ICGC mode
#python3 ../FiNGS.py -n normal.bam -t tumor.bam -v s2.raw.vcf -r /path/to/reference/genome.fa --PASSonlyin --PASSonlyout --ICGC # 23 PASS

## Other options
#python3 ../FiNGS.py -n normal.bam -t tumor.bam -v s2.raw.vcf --PASSonlyin # 16 PASS, 38 total
#python3 ../FiNGS.py -n normal.bam -t tumor.bam -v s2.raw.vcf --PASSonlyout # 20 PASS
#python3 ../FiNGS.py -n normal.bam -t tumor.bam -v s2.raw.vcf # 20 PASS, 1207 total


