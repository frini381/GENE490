module load BCFtools
id=2559_2562_2560_2561
wkdir=/mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2562_2560_2561

######Het in all three - pat

input=${wkdir}/${id}_20230707_hg38_Split_ann_A0.01_HIGHMOD_DPGQ.vcf.gz
output=${wkdir}/${id}_20230707_hg38_Split_ann_A0.01_HIGHMOD_DPGQ_CompHet1.vcf.gz
if [[ ! -f ${output} ]] #if the output doesn't exist then...
then
    echo "${id} ARchet filter"
               bcftools view -i '(FORMAT/GT[0]=="het" && FORMAT/GT[1]=="ref" && FORMAT/GT[2]=="het" && FORMAT/GT[3]=="het")' --output ${output} --output-type z ${input}
##           bcftools view -i '(FORMAT/GT[0]=="het" & FORMAT/GT[1]=="het")' --output ${output} --output-type z ${input}
fi

/home/frini381/bin/vcf_human_readable.sh $output

input=${wkdir}/${id}_20230707_hg38_Split_ann_A0.01_HIGHMOD_DPGQ_CompHet1.vcf.gz
output=${wkdir}/${id}_20230707_hg38_Split_ann_A0.01_HIGHMOD_DPGQ_CompHet1_genelist.txt

bcftools +split-vep ${input} -f '%Gene\n' | sed 's/,/\n/g' | sort -u > $output #making the genelist

######Het in just sons - mat

input=${wkdir}/${id}_20230707_hg38_Split_ann_A0.01_HIGHMOD_DPGQ.vcf.gz
output=${wkdir}/${id}_20230707_hg38_Split_ann_A0.01_HIGHMOD_DPGQ_CompHet2.vcf.gz
if [[ ! -f ${output} ]]
then
    echo "${id} ARchet filter"
               bcftools view -i '(FORMAT/GT[0]=="ref" && FORMAT/GT[1]=="het" && FORMAT/GT[2]=="het" && FORMAT/GT[3]=="het")' --output ${output} --output-type z ${input}
##           bcftools view -i '(FORMAT/GT[0]=="het" & FORMAT/GT[1]=="het")' --output ${output} --output-type z ${input}
fi

/home/frini381/bin/vcf_human_readable.sh $output

input=${wkdir}/${id}_20230707_hg38_Split_ann_A0.01_HIGHMOD_DPGQ_CompHet2.vcf.gz
output=${wkdir}/${id}_20230707_hg38_Split_ann_A0.01_HIGHMOD_DPGQ_CompHet2_genelist.txt

bcftools +split-vep ${input} -f '%Gene\n' | sed 's/,/\n/g' | sort -u > $output #making the genelist

######What is in common between two

pat=${wkdir}/${id}_20230707_hg38_Split_ann_A0.01_HIGHMOD_DPGQ_CompHet1_genelist.txt
mat=${wkdir}/${id}_20230707_hg38_Split_ann_A0.01_HIGHMOD_DPGQ_CompHet2_genelist.txt
comm -1 -2 $pat $mat > ${wkdir}/${id}_20230106_hg38_split_ann_AF0.01_HIGHMOD_CompHet_shared_genelist.txt



/mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2562_2560_2561/2559_2562_2560_2561_20230707_hg38_Split_ann_A0.01_HIGHMOD_DPGQ.vcf.gz