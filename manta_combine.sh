#################
###variables#####
#################
fam_id=2559_2562_2560_2561
IDs=2259,2562,2560,2561
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
        file=/mnt/hcs/WCHP_Clinical_Genetics/SequenceData/hg38/StructuralCalls/Genomic/Manta/${id}/Manta_Filtering/${id}_manta_sv_duphold_ann_${type}_DHFFC_Masked_Rarity_Coding.vcf.gz
        ls $file
        echo $file >> list_to_merge_${type}.txt
    done
    bcftools merge -l list_to_merge_${type}.txt -o /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2562_2560_2561/Manta/${fam_id}_manta_sv_duphold_ann_${type}_Rarity_Coding_fam.vcf.gz   
    ls /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2562_2560_2561/Manta/${fam_id}_manta_sv_duphold_ann*Rarity_Coding_fam.vcf.gz > vcf_list.txt
done
for sv in DEL DUP INS INV BND
   	do 
   	bin/vcf_human_readable.sh /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2562_2560_2561/Manta/${fam_id}_manta_sv_duphold_ann_${sv}_Rarity_Coding_fam.vcf.gz
done

