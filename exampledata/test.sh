#!/bin/bash

python3 ../FiNGS.py -n normal.bam -t tumor.bam -v s2.raw.vcf --PASSonlyin --PASSonlyout # 16 PASS only
#python3 ../FiNGS.py -n normal.bam -t tumor.bam -v s2.raw.vcf --PASSonlyin # 16 PASS, 38 total
#python3 ../FiNGS.py -n normal.bam -t tumor.bam -v s2.raw.vcf --PASSonlyout # 20 PASS
#python3 ../FiNGS.py -n normal.bam -t tumor.bam -v s2.raw.vcf # 20 PASS, 1207 total


