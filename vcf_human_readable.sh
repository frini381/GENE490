#!/bin/bash
# Written by Ben Halliday

###  Modules  ###
module purge
module load snpEff &> /dev/null
module load HTSlib/1.11-GCC-10.2.0 &> /dev/null

###  Inputs   ###
vcf=$1

vcf_human_readable () {
    input=$1
    temp=$(echo ${input} | sed "s/\.gz//g")
    output=$(echo ${input} | sed "s/\.vcf.gz/_out\.tsv/g")
    bgzip -d $input -c > ${temp}
    variables=; for var in `less ${temp} | grep "^##FORMAT" | cut -d\, -f1 | sed "s/##FORMAT=<ID=//g"`; do if [ -z "$variables" ]; then variables="GEN["'$i'"].$var"; else variables="$variables GEN["'$i'"].$var"; fi; done; 
    number=$(cat ${temp} | grep '^#CHROM' | cut -f10-  | tr '\t' "\n" | wc -l) && number="$(($number-1))"
    sample=;for i in `seq 0 1 $number`; do if [ -z "$sample" ]; then sample=$(eval echo $variables); else sample="${sample} "$(eval echo $variables)""; fi; done;
    varall=;for var in `less ${temp} | grep "##INFO" | cut -d\, -f1 | sed "s/##INFO=<ID=//g" | grep -v "-"`; do if [ -z "$varall" ]; then varall="$var"; else varall="$varall $var"; fi; done;
    /resource/easybuild/software/Java/1.8.0_121/bin/java -Xmx3g -jar $EBROOTSNPEFF/SnpSift.jar extractFields ${temp}  "CHROM" "POS" "REF" "ALT" "QUAL" "FILTER" $varall $sample > ${output}
    rm ${temp}
}

if [[ ! -f "$vcf" ]]
then
    echo "Input VCF file ${vcf} doesn't exist, exiting..."
else
    vcf_human_readable ${vcf}
    wait
fi
  
exit 0
