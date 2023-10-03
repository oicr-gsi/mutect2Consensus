#!/bin/bash
set -euo pipefail

cd $1

module load vcftools/0.1.16 2>/dev/null

ls | sed 's/.*\.//' | sort | uniq -c

for i in *.vcf.gz; do vcf-query -l ${i}; done | sort

for i in *.vcf.gz; do zcat ${i} | grep -v ^# | cut -f 1 | uniq| md5sum; done | sort

for i in *.maf.gz; do zcat ${i} | md5sum; done | sort
