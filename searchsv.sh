#!/bin/bash
chr=$1
start=$2
end=$3
prefix="/lustre/scratch115/realdata/mdt0/projects/t144_helic_15x/analysis/HA/release/chr"
suffix=""
bcftools view -r 'chr'$chr':'${start}'0000000-'${end}'0000000' ${prefix}${chr}${suffix}.vcf.gz  | bcftools query -f '%CHROM\t%POS\t[%DP\t]\n' > $chr.${start}-${end}M.alldp
./call_sv.R $chr.${start}-${end}M.alldp chr$chr.avg.depth.persample $chr.${start}-${end}M.sv
