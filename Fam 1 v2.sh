module load BCFtools

#remember to load Xming first
module load IGV 

/bin/vcf_human_readable.sh 

#if WCHP doesn't show up
ls /mnt/hcs/

srun -c 2 --mem 1G bin/SplitOut.sl /mnt/hcs/WCHP_Clinical_Genetics/SequenceData/hg38/VariantCalls/Genomic/20230707_Genomic_GATK_Genomics_Med_Cohort/SplitMultiallelics/20230707_Genomic_GATK_Genomics_Med_Cohort_Split_Annotated_1/20230707_Genomic_GATK_Genomics_Med_Cohort_Split_ann.vcf.gz /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2562_2560_2561/2559_2562_2560_2561_20230707_hg38_Split_ann.vcf.gz 2559,2562,2560,2561

bcftools filter -i '(gnomADv3_AF[*] <= '0.01' | gnomADv3_AF[*] = ".")' /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2562_2560_2561/2559_2562_2560_2561_20230707_hg38_Split_ann.vcf.gz -Oz -o /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2562_2560_2561/2559_2562_2560_2561_20230707_hg38_Split_ann_A0.01.vcf.gz


bcftools filter -i 'CSQ[*] ~ "|HIGH|" | CSQ[*] ~ "|MODERATE|"' /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2562_2560_2561/2559_2562_2560_2561_20230707_hg38_Split_ann_A0.01.vcf.gz -Oz -o /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2562_2560_2561/2559_2562_2560_2561_20230707_hg38_Split_ann_A0.01_HIGHMOD.vcf.gz 