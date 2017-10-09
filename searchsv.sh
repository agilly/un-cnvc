#!/bin/bash
chr=$1
start=$2
end=$3

#########################################
## Warning: This script uses 2 threads ##
#########################################

## Edit below to generate the path for your VCFs
prefix=""
suffix=""

## By default, runs 10Mb chunks. Reduce if your machine is low in RAM or has a lousy processor.

bcftools view -r 'chr'$chr':'${start}'0000000-'${end}'0000000' ${prefix}${chr}${suffix}.vcf.gz  | bcftools query -f '%CHROM\t%POS\t[%DP\t]\n' > $chr.${start}-${end}M.alldp
./call_sv.R $chr.${start}-${end}M.alldp chr$chr.avg.depth.persample $chr.${start}-${end}M.sv
