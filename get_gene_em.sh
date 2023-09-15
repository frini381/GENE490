#!/bin/bash
#chsh -s /bin/bash

##variables assigned
id=$1 #4468
region=$2 #chr7:156983249-157021198
gene=$3 #MNX1
wkdir=/mnt/hcs/WCHP_Clinical_Genetics/Em/Genome_Analysis/${id}/Gene_Search
inputvcf=/mnt/hcs/WCHP_Clinical_Genetics/Em/Genome_Analysis/${id}/${id}_family.vcf.gz


##modules
module purge
module load BCFtools/1.16-GCC-11.2.0


##############
####code######
##############
mkdir -p $wkdir
if [[ ! -f ${inputvcf}.tbi ]] 
then
	bcftools index --tbi ${inputvcf}
fi	
bcftools view --regions $region --output ${wkdir}/${gene}.vcf.gz --output-type z ${inputvcf}
vcf_human_readable.sh ${wkdir}/${gene}.vcf.gz

exit 0