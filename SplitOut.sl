#!/bin/bash
#SBATCH --job-name SplitSamples
#SBATCH --time	12:00:00
#SBATCH --mem	1G
#SBATCH --cpus-per-task	2
#SBATCH --error SplitSamples_%j.out
#SBATCH --output SplitSamples_%j.out

# SpliitOut.sl
# Written by Mini on 19_07_2019
# This script is intended to take an input vcf and extract samples of interest from the entire cohort that have been variant called together as a batch
# USAGE: sbatch SpliitOut.sl input.vcf.gz output.vcf.gz sample1,sample2,sample3...

samples=$3
input=$1
output=$2
set -o pipefail
if [ -f ${samples} ];then
	samples="-S ${samples}"
else
	samples="-s ${samples}"
fi
module purge
module load BCFtools

$(which bcftools) view -a ${samples} ${input} --force-samples -Ou | $(which bcftools) view -e 'ALT="."' -Oz -o ${output} || exit 1
$(which bcftools) index -t ${output} || exit 1
echo "The samples ${samples} have been successfully extracted from ${input} and wriiten to ${output}"


