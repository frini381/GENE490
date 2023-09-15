dsmcfs0.uod.otago.ac.nz

module load BCFtools

#remember to load Xming first
module load IGV 

#Split out family vcf
srun -c 2 --mem 1G bin/SplitOut.sl /mnt/hcs/WCHP_Clinical_Genetics/SequenceData/GRCh37/VariantCalls/Genomic/20230417_Genomic_GATK_Whole_Cohort/SplitMultiallelics/20230417_Genomic_GATK_Whole_Cohort_Split_Annotated_1/20230417_Genomic_GATK_Whole_Cohort_Split_ann.vcf.gz /mnt/hcs/WCHP_Clinical_Genetics/Niamh/4584_4585_4586/4584_4585_4586_20230417_GRCh37_split_ann.vcf.gz


#STRUCTURAL VARIANTS


#################
###variables#####
#################
fam_id=4584_4585_4586
IDs=4584,4585,4586
types=DEL,DUP,INS,INV,BND


#################
###modules#######
#################
module purge
module load rtg-core
module load BCFtools/1.16-GCC-11.2.0



rm list_to_merge_*
echo $types | tr ',' '\n' | while read -r type
do
    echo $IDs | tr ',' '\n' | while read -r id
    do
        echo $id
        file=/mnt/hcs/WCHP_Clinical_Genetics/SequenceData/GRCh37/StructuralCalls/Genomic/Manta/${id}/Manta_Filtering/${id}_manta_sv_duphold_ann_${type}*Rarity_Coding.vcf.gz
        ls $file
        echo $file >> list_to_merge_${type}.txt
    done
    bcftools merge -l list_to_merge_${type}.txt -o /mnt/hcs/WCHP_Clinical_Genetics/Niamh/4584_4585_4586/Manta/${fam_id}_manta_sv_duphold_ann_${type}_Rarity_Coding_fam.vcf.gz   
    ls /mnt/hcs/WCHP_Clinical_Genetics/Niamh/4584_4585_4586/Manta/${fam_id}_manta_sv_duphold_ann*Rarity_Coding_fam.vcf.gz > vcf_list.txt
done
for sv in DEL DUP INS INV BND
   	do 
   	bin/vcf_human_readable.sh /mnt/hcs/WCHP_Clinical_Genetics/Niamh/4584_4585_4586/Manta/${fam_id}_manta_sv_duphold_ann_${sv}_Rarity_Coding_fam.vcf.gz
done



#Compound HET



module load BCFtools
id=4584_4585_4586
wkdir=/mnt/hcs/WCHP_Clinical_Genetics/Niamh/4584_4585_4586

######Het in 4585 and 4584, absent in 4586 - pat

input=/mnt/hcs/WCHP_Clinical_Genetics/Niamh/4584_4585_4586/4584_4585_4586_20230417_GRCh37_split_ann_A0.01_HIGHMOD_DPGQ.vcf.gz
output=${wkdir}/${id}_20230417_GRCh37_split_ann_AF0.01_HIGHMOD_CompHet1.vcf.gz
if [[ ! -f ${output} ]] #if the output doesn't exist then...
then
    echo "${id} ARchet filter"
               bcftools view -i '(FORMAT/GT[0]=="het" && FORMAT/GT[1]=="ref" && FORMAT/GT[2]=="het")' --output ${output} --output-type z ${input}
##           bcftools view -i '(FORMAT/GT[0]=="het" & FORMAT/GT[1]=="het")' --output ${output} --output-type z ${input}
fi

/home/frini381/bin/vcf_human_readable.sh $output

input=${wkdir}/${id}_20230417_GRCh37_split_ann_AF0.01_HIGHMOD_CompHet1.vcf.gz
output=${wkdir}/${id}_20230417_GRCh37_split_ann_AF0.01_HIGHMOD_CompHet1_genelist.txt

bcftools +split-vep ${input} -f '%Gene\n' | sed 's/,/\n/g' | sort -u > $output #making the genelist

######Het in 4586 & 4584, absent in 4585 - mat

input=/mnt/hcs/WCHP_Clinical_Genetics/Niamh/4584_4585_4586/4584_4585_4586_20230417_GRCh37_split_ann_A0.01_HIGHMOD_DPGQ.vcf.gz
output=${wkdir}/${id}_20230417_GRCh37_split_ann_AF0.01_HIGHMOD_CompHet2.vcf.gz
if [[ ! -f ${output} ]]
then
    echo "${id} ARchet filter"
               bcftools view -i '(FORMAT/GT[0]=="ref" && FORMAT/GT[1]=="het" && FORMAT/GT[2]=="het")' --output ${output} --output-type z ${input}
##           bcftools view -i '(FORMAT/GT[0]=="het" & FORMAT/GT[1]=="het")' --output ${output} --output-type z ${input}
fi

/home/frini381/bin/vcf_human_readable.sh $output

input=${wkdir}/${id}_20230417_GRCh37_split_ann_AF0.01_HIGHMOD_CompHet2.vcf.gz
output=${wkdir}/${id}_20230417_GRCh37_split_ann_AF0.01_HIGHMOD_CompHet2_genelist.txt

bcftools +split-vep ${input} -f '%Gene\n' | sed 's/,/\n/g' | sort -u > $output #making the genelist

######What is in common between two

pat=${wkdir}/${id}_20230417_GRCh37_split_ann_AF0.01_HIGHMOD_CompHet1_genelist.txt
mat=${wkdir}/${id}_20230417_GRCh37_split_ann_AF0.01_HIGHMOD_CompHet2_genelist.txt
comm -1 -2 $pat $mat > ${wkdir}/${id}_20230417_GRCh37_split_ann_AF0.01_HIGHMOD_CompHet_shared_genelist.txt



### Filtering on quality

bcftools filter -i '(FORMAT/DP[0]>=8 && FORMAT/DP[1]>=8 && FORMAT/DP[2]>=8)' /mnt/hcs/WCHP_Clinical_Genetics/Niamh/4584_4585_4586/4584_4585_4586_20230417_GRCh37_split_ann_A0.01_HIGHMOD.vcf.gz -Oz -o /mnt/hcs/WCHP_Clinical_Genetics/Niamh/4584_4585_4586/4584_4585_4586_20230417_GRCh37_split_ann_A0.01_HIGHMOD_DP.vcf.gz


bcftools filter -i '(FORMAT/GQ[0]>=10 && FORMAT/GQ[1]>=10 && FORMAT/GQ[2]>=10)' /mnt/hcs/WCHP_Clinical_Genetics/Niamh/4584_4585_4586/4584_4585_4586_20230417_GRCh37_split_ann_A0.01_HIGHMOD_DP.vcf.gz -Oz -o /mnt/hcs/WCHP_Clinical_Genetics/Niamh/4584_4585_4586/4584_4585_4586_20230417_GRCh37_split_ann_A0.01_HIGHMOD_DPGQ.vcf.gz
