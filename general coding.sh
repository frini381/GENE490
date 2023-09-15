ssh -Y frini381@dsmcfs0.uod.otago.ac.nz

dsmcfs0.uod.otago.ac.nz

module load BCFtools

#remember to load Xming first
module load IGV 

/home/frini381/bin/vcf_human_readable.sh

#if WCHP doesn't show up
ls /mnt/hcs/

#Gene name extraction
bcftools filter -i 'CSQ[*] ~ "|ARF6|" | CSQ[*] ~ "|CBL|" | CSQ[*] ~ "|COL11A1|" | CSQ[*] ~ "|COLL11A2|" | CSQ[*] ~ "|COL1A1|" | CSQ[*] ~ "|COL1A2|" | CSQ[*] ~ "|COL24A1|" | CSQ[*] ~ "|COL27A1|" | CSQ[*] ~ "|COL2A1|" | CSQ[*] ~ "|COL3A1|" | CSQ[*] ~ "|COL5A1|" | CSQ[*] ~ "|COL5A2|" | CSQ[*] ~ "|COL5A3|" | CSQ[*] ~ "|CRK|" | CSQ[*] ~ "|CRKL|" | CSQ[*] ~ "|DOCK7|" | CSQ[*] ~ "|EPS15|" | CSQ[*] ~ "|FN1|" | CSQ[*] ~ "|GAB1|" | CSQ[*] ~ "|GGA3|" | CSQ[*] ~ "|GRB2|" | CSQ[*] ~ "|HGF|" | CSQ[*] ~ "|HGFAC|" | CSQ[*] ~ "|HGS|" | CSQ[*] ~ "|HPN|" | CSQ[*] ~ "|HRAS|" | CSQ[*] ~ "|ITGA2|" | CSQ[*] ~ "|ITGA3|" | CSQ[*] ~ "|ITGB1|" | CSQ[*] ~ "|KRAS|" | CSQ[*] ~ "|LAMA1|" | CSQ[*] ~ "|LAMA2|" | CSQ[*] ~ "|LAMA3|" | CSQ[*] ~ "|LAMA4|" | CSQ[*] ~ "|LAMA5|" | CSQ[*] ~ "|LAMB1|" | CSQ[*] ~ "|LAMB2|" | CSQ[*] ~ "|LAMB3|" | CSQ[*] ~ "|LAMC1|" | CSQ[*] ~ "|LAMC2|" | CSQ[*] ~ "|LAMC3|" | CSQ[*] ~ "|LRIG1|" | CSQ[*] ~ "|MET|" | CSQ[*] ~ "|MUC20|" | CSQ[*] ~ "|NRAS|" | CSQ[*] ~ "|PIK3CA|" | CSQ[*] ~ "|PIKR1|" | CSQ[*] ~ "|PTK2|" | CSQ[*] ~ "|PTPN1|" | CSQ[*] ~ "|PTPN11|" | CSQ[*] ~ "|PTPN2|" | CSQ[*] ~ "|PTPRJ|" | CSQ[*] ~ "|RAB4A|" | CSQ[*] ~ "|RAB4B|" | CSQ[*] ~ "|RAC1|" | CSQ[*] ~ "|RANBP10|" | CSQ[*] ~ "|RANBP9|" | CSQ[*] ~ "|RAP1A|" | CSQ[*] ~ "|RAP1B|" | CSQ[*] ~ "|RAPGEF1|" | CSQ[*] ~ "|RPS27A|" | CSQ[*] ~ "|SH3GL1|" | CSQ[*] ~ "|SH3GL2|" | CSQ[*] ~ "|SH3GL3|" | CSQ[*] ~ "|SH3KBP1|" | CSQ[*] ~ "|SHC1|" | CSQ[*] ~ "|SOS1|" | CSQ[*] ~ "|SPINT1|" | CSQ[*] ~ "|SPINT2|" | CSQ[*] ~ "|SRC|" | CSQ[*] ~ "|STAM|" | CSQ[*] ~ "|STAM2|" | CSQ[*] ~ "|STAT3|" | CSQ[*] ~ "|TNS3|" | CSQ[*] ~ "|TNS4|" | CSQ[*] ~ "|UBA52|" | CSQ[*] ~ "|UBB|" | CSQ[*] ~ "|UBC|" | CSQ[*] ~ "|USP8|"' /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2560_2561/2559_2560_2561_20230106_hg38_split_ann.vcf.gz -Oz -o /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2560_2561/2559_2560_2561_20230106_hg38_split_ann_METGenes.vcf.gz

#Split vcf for a single pedigree
srun -c 2 --mem 1G bin/SplitOut.sl /mnt/hcs/WCHP_Clinical_Genetics/SequenceData/hg38/VariantCalls/Genomic/20230106_Genomic_GATK_Genomics_Med_Cohort/SplitMultiallelics/20230106_Genomic_GATK_Genomics_Med_Cohort_Split_Annotated_1/20230106_Genomic_GATK_Genomics_Med_Cohort_Split_ann.vcf.gz /mnt/hcs/WCHP_Clinical_Genetics/Niamh/4430/4427_4428_4429_4430_20230106_hg38_Split_ann.vcf.gz 4427,4428,4429,4430


bcftools filter -i '(gnomADv3_AF[*] <= '0.01' | gnomADv3_AF[*] = ".")' /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2560_2561/2559_2560_2561_20230106_hg38_Split_ann_ChrX.vcf.gz -Oz -o /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2560_2561/2559_2560_2561_20230106_hg38_Split_ann_ChrX_A0.01.vcf.gz 

bcftools filter -i 'CSQ[*] ~ "|HIGH|" | CSQ[*] ~ "|MODERATE|"' /mnt/hcs/WCHP_Clinical_Genetics/Niamh/4430/4427_4428_4429_4430_20230106_hg38_Split_ann_A0.01.vcf.gz -Oz -o /mnt/hcs/WCHP_Clinical_Genetics/Niamh/4430/4427_4428_4429_4430_20230106_hg38_Split_ann_A0.01_HIGHMOD.vcf.gz 

bcftools filter -i 'ClinVar_CLNSIG ~ "Pathogenic\|Likely_pathogenic\|Affects"' /mnt/hcs/WCHP_Clinical_Genetics/Niamh/4430/4427_4428_4429_4430_20230106_hg38_Split_ann.vcf.gz -Oz -o /mnt/hcs/WCHP_Clinical_Genetics/Niamh/4430/4427_4428_4429_4430_20230106_hg38_Split_ann_VarEff.vcf.gz 

bin/vcf_human_readable.sh /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2560_2561/Gene_Search/MET_A0.01.vcf.gz

#Chromosome X extraction
bcftools filter -i 'CHROM ~ "chrX"' /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2560_2561/2559_2560_2561_20230106_hg38_split_ann.vcf.gz -Oz -o /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2560_2561/2559_2560_2561_20230106_hg38_split_ann_ChrX.vcf.gz

#Inheritance Filtering
bcftools filter -i 'GT[0]="ref" && GT[1]="het" && GT[2]="het"' /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2560_2561/Gene_Search/MET.vcf.gz -Oz -o /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2560_2561/Gene_Search/MET_A0.01_pat.vcf.gz

#Rare MET variants
bcftools filter -i '(gnomADv3_AF[*] <= '0.01' | gnomADv3_AF[*] = ".")' /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2560_2561/Gene_Search/MET.vcf.gz -Oz -o /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2560_2561/Gene_Search/MET_A0.01.vcf.gz

#Rare MET variants all three have 
bcftools filter -i 'GT[0]="het" && GT[1]="het" && GT[2]="het"' /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2560_2561/Gene_Search/MET_A0.01.vcf.gz -Oz -o /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2560_2561/Gene_Search/MET_A0.01_pat.vcf.gz