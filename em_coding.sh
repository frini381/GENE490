###############################
###### AR chet model ##########
###############################

###########pat list#######

fam_id=2559_2560_2561
IDs=2559,2560,2561

input=${wkdir}/${id}_family_rarity.vcf.gz
output=${wkdir}/${id}_family_rarity_ARchet_pat.vcf.gz
if [[ ! -f ${output} ]]
then
    echo "${id} ARchet filter"
               bcftools view -i '(FORMAT/GT[0]=="het" && FORMAT/GT[1]=="het" && FORMAT/GT[2]=="ref")' --output ${output} --output-type z ${input}
##           bcftools view -i '(FORMAT/GT[0]=="het" & FORMAT/GT[1]=="het")' --output ${output} --output-type z ${input}
fi

 ##coding region only
input=${wkdir}/${id}_family_rarity_ARchet_pat.vcf.gz
output=${wkdir}/${id}_family_rarity_ARchet_pat_coding.vcf.gz

rtg vcffilter -i ${input} -o ${output} --javascript /home/jamem126/bin/filter-csq-damaging.js

##quality where depth is greater than or equal to 4
input=${wkdir}/${id}_family_rarity_ARchet_pat_coding.vcf.gz
output=${wkdir}/${id}_family_rarity_ARchet_pat_coding_qual.vcf.gz

bcftools view -i 'ALT!="*" && (FORMAT/DP[0]>=4 && FORMAT/DP[1]>=4 && FORMAT/DP[2]>=4)' --output ${output} --output-type z ${input} 
vcf_human_readable.sh $output

input=${wkdir}/${id}_family_rarity_ARchet_pat_coding_qual.vcf.gz
output=${wkdir}/${id}_family_rarity_ARchet_pat_coding_qual_genelist.txt

bcftools +split-vep ${input} -f '%Gene\n' | sed 's/,/\n/g' | sort -u > $output

################mat list#################
input=${wkdir}/${id}_family_rarity.vcf.gz
output=${wkdir}/${id}_family_rarity_ARchet_mat.vcf.gz
if [[ ! -f ${output} ]]
then
    echo "${id} ARchet filter"
               bcftools view -i '(FORMAT/GT[0]=="het" && FORMAT/GT[1]=="ref" && FORMAT/GT[2]=="het")' --output ${output} --output-type z ${input}
##           bcftools view -i '(FORMAT/GT[0]=="het" & FORMAT/GT[1]=="het")' --output ${output} --output-type z ${input}
fi

 ##coding region only
input=${wkdir}/${id}_family_rarity_ARchet_mat.vcf.gz
output=${wkdir}/${id}_family_rarity_ARchet_mat_coding.vcf.gz

rtg vcffilter -i ${input} -o ${output} --javascript /home/jamem126/bin/filter-csq-damaging.js

##quality where depth is greater than or equal to 4
input=${wkdir}/${id}_family_rarity_ARchet_mat_coding.vcf.gz
output=${wkdir}/${id}_family_rarity_ARchet_mat_coding_qual.vcf.gz

bcftools view -i 'ALT!="*" && (FORMAT/DP[0]>=4 && FORMAT/DP[1]>=4 && FORMAT/DP[2]>=4)' --output ${output} --output-type z ${input} 
vcf_human_readable.sh $output


input=${wkdir}/${id}_family_rarity_ARchet_mat_coding_qual.vcf.gz
output=${wkdir}/${id}_family_rarity_ARchet_mat_coding_qual_genelist.txt

bcftools +split-vep ${input} -f '%Gene\n' | sed 's/,/\n/g' | sort -u > $output


############compare gene lists#######################
pat=${wkdir}/${id}_family_rarity_ARchet_pat_coding_qual_genelist.txt
mat=${wkdir}/${id}_family_rarity_ARchet_mat_coding_qual_genelist.txt
comm -1 -2 $pat $mat > ${wkdir}/${id}_family_rarity_ARchet_shared_coding_qual_genelist.txt
