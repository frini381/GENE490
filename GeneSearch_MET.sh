#!/bin/bash 
#Ben Halliday 
#Created 11/22


###  Modules  ###
module purge
module load rtg-core
module load BCFtools
module load snpEff

input=/mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2562_2560_2561/2559_2562_2560_2561_20230707_hg38_Split_ann.vcf.gz
id=2560
regionsfile=/mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2562_2560_2561/bed_file_regions_of_met_genes.tsv
wkdir=/mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2562_2560_2561/GeneSearch/2559_2562_2560_2561_20230707_hg38_Split_ann

mkdir -p ${wkdir}

    vcf_human_readable () {
        input=$1
        temp=$(echo ${input} | sed "s/\.gz//g")
        output=$(echo ${input} | sed "s/\.vcf.gz/_out\.tsv/g")
        bgzip -d $input -c > ${temp}
        variables=; for var in `less ${temp} | grep "^##FORMAT" | cut -d\, -f1 | sed "s/##FORMAT=<ID=//g"`; do if [ -z "$variables" ]; then variables="GEN["'$i'"].$var"; else variables="$variables GEN["'$i'"].$var"; fi; done; 
        number=$(cat ${temp} | grep '^#CHROM' | cut -f10-  | tr '\t' "\n" | wc -l) && number="$(($number-1))"
        sample= && for i in `seq 0 1 $number`; do if [ -z "$sample" ]; then sample=$(eval echo $variables); else sample="${sample} "$(eval echo $variables)""; fi; done;
        varall=; for var in `less ${temp} | grep "##INFO" | cut -d\, -f1 | sed "s/##INFO=<ID=//g" | grep -v "-"`; do if [ -z "$varall" ]; then varall="$var"; else varall="$varall $var"; fi; done;
        /resource/easybuild/software/Java/1.8.0_121/bin/java -Xmx8g -jar $EBROOTSNPEFF/SnpSift.jar extractFields ${temp}  "CHROM" "POS" "REF" "ALT"  $varall $sample > ${output}
        rm ${temp}
    }

## Delete firsyt row then retain first and fifth column
sed 1d $regionsfile | cut -f1-4 | while read -r gene chr start end 
do
echo "${chr}:${start}-${end}"
rtg vcffilter -i ${input} -o ${wkdir}/${id}_family_${gene}.vcf.gz --region "${chr}:${start}-${end}" --Xforce
vcf_human_readable ${wkdir}/${id}_family_${gene}.vcf.gz
done



exit
