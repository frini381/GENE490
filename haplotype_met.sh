##variables assigned
id=2559_2562_2560_2561
region=chr7:116671196-116799386
gene=MET
wkdir=/mnt/hcs/WCHP_Clinical_Genetics/Niamh/${id}/Haplotypes
inputvcf=/mnt/hcs/WCHP_Clinical_Genetics/Niamh/${id}/${id}_20230707_hg38_split_ann.vcf.gz


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
/home/frini381/bin/vcf_human_readable.sh ${wkdir}/${gene}.vcf.gz

exit 0

#looking at met and region 10kb either side


#MAT
bcftools filter -i 'GT[0]="ref" && GT[1]="het" && GT[2]="het" && GT[3]="het"' /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2562_2560_2561/Haplotypes/MET.vcf.gz -Oz -o /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2562_2560_2561/Haplotypes/MET_mat.vcf.gz


#PAT
bcftools filter -i 'GT[0]="het" && GT[1]="ref" && GT[2]="het" && GT[3]="het"' /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2562_2560_2561/Haplotypes/MET.vcf.gz -Oz -o /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2562_2560_2561/Haplotypes/MET_pat.vcf.gz


bin/vcf_human_readable.sh /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2562_2560_2561/Haplotypes/MET_mat.vcf.gz

bin/vcf_human_readable.sh /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2562_2560_2561/Haplotypes/MET_pat.vcf.gz

