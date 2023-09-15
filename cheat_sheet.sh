ssh -Y frini381@dsmcfs0.uod.otago.ac.nz

dsmcfs0.uod.otago.ac.nz

module load BCFtools

module load IGV

/bin/vcf_human_readable.sh 

#Annotate
/resource/pipelines/Variant_Annotation/Annotate_my_VCF.sh -d /resource/pipelines/Variant_Annotation/defaults.txt -i /mnt/hcs/WCHP_Clinical_Genetics/SequenceData/GRCh37/VariantCalls/Genomic/Genomic_R_220802_TIMMOR1_DNA_M001-R_220707_TIMMOR1_DNA_M001/SplitMultiallelics/Genomic_R_220802_TIMMOR1_DNA_M001-R_220707_TIMMOR1_DNA_M001_Split.vcf.gz

#Use SplitVCF
srun -c 2 --mem 1G bin/SplitOut.sl /mnt/hcs/WCHP_Clinical_Genetics/SequenceData/GRCh37_old/VariantCalls/Genomic/Genomic_R_220802_TIMMOR1_DNA_M001-R_220707_TIMMOR1_DNA_M001/SplitMultiallelics/Genomic_R_220802_TIMMOR1_DNA_M001-R_220707_TIMMOR1_DNA_M001_Split_Annotated_1/Genomic_R_220802_TIMMOR1_DNA_M001-R_220707_TIMMOR1_DNA_M001_Split_ann.vcf.gz 4524_Genomic_R_220802_TIMMOR1_DNA_M001-R_220707_TIMMOR1_DNA_M001_Split_ann.vcf.gz 4524

#Copy vcf to relevant DHB
cp 4524_Genomic_R_220802_TIMMOR1_DNA_M001-R_220707_TIMMOR1_DNA_M001_Split_ann.vcf.gz /mnt/hcs/WCHP_Clinical_Genetics/DHB_Genomes/Christchurch/4524/

#Allele freq filter
bcftools filter -i '(gnomADv3_AF[*] <= '0.01' | gnomADv3_AF[*] = ".")' 4524_Genomic_R_220802_TIMMOR1_DNA_M001-R_220707_TIMMOR1_DNA_M001_Split_ann.vcf.gz -Oz -o 4524_Genomic_R_220802_TIMMOR1_DNA_M001-R_220707_TIMMOR1_DNA_M001_Split_ann_AF0.01.vcf.gz

#Remove * from REF/ALT columns
bcftools filter -e '(REF = "*" | ALT = "*")' 4524_Genomic_R_220802_TIMMOR1_DNA_M001-R_220707_TIMMOR1_DNA_M001_Split_ann_AF0.01.vcf.gz -Oz -o 4524_Genomic_R_220802_TIMMOR1_DNA_M001-R_220707_TIMMOR1_DNA_M001_Split_ann_AF0.01_nostar.vcf.gz

#Impact
 bcftools filter -i 'CSQ[*] ~ "|HIGH|" | CSQ[*] ~ "|MODERATE|"' 4524_Genomic_R_220802_TIMMOR1_DNA_M001-R_220707_TIMMOR1_DNA_M001_Split_ann_AF0.01_nostar.vcf.gz -Oz -o 4524_Genomic_R_220802_TIMMOR1_DNA_M001-R_220707_TIMMOR1_DNA_M001_Split_ann_AF0.01_nostar_HIGHMOD.vcf.gz

 #ClinVar
  bcftools filter -i 'ClinVar_CLNSIG ~ "Pathogenic\|Likely_pathogenic\|Affects"' 4524_Genomic_R_220802_TIMMOR1_DNA_M001-R_220707_TIMMOR1_DNA_M001_Split_ann_AF0.01_nostar_HIGHMOD.vcf.gz -Oz -o 4524_Genomic_R_220802_TIMMOR1_DNA_M001-R_220707_TIMMOR1_DNA_M001_Split_ann_AF0.01_nostar_HIGHMOD_ClinVar.vcf.gz

  #Text output
  bin/vcf_human_readable.sh 4458fam_Genomic_R_220802_TIMMOR1_DNA_M001-R_220707_TIMMOR1_DNA_M001_Split_ann_AF0.01_nostar_HIGHMOD_denovo.vcf.gz

  #Gene name extraction
  bcftools filter -i 'CSQ[*] ~ "||" | CSQ[*] ~ "||"' or 'ANN[*] ~ "|" | ANN[*] ~"|"'    4524_Genomic_R_220802_TIMMOR1_DNA_M001-R_220707_TIMMOR1_DNA_M001_Split_ann.vcf.gz -Oz -o 4524_Genomic_R_220802_TIMMOR1_DNA_M001-R_220707_TIMMOR1_DNA_M001_Split_ann_EIF2B.vcf.gz

  #SpliceAI filter
  module load VEP
  filter_vep --vcf_info_field SpliceAI --filter "DS_AG >= 0.5 or DS_AL >= 0.5 or DS_DG >= 0.5 or DS_DL >= 0.5" -i 4524_Genomic_R_220802_TIMMOR1_DNA_M001-R_220707_TIMMOR1_DNA_M001_Split_ann_AF0.01_nostar.vcf.gz --gz | bgzip > 4524_Genomic_R_220802_TIMMOR1_DNA_M001-R_220707_TIMMOR1_DNA_M001_Split_ann_AF0.01_nostar_SpliceAI.vcf.gz

#Our NeSI project
/nesi/project/uoo00032/

#To write the command file for fastqs (to start alingments on NeSI)
 /mnt/hcs/WCHP_Clinical_Genetics/SequenceData/Utilities/make_Fastq2VCF_commands.sh
  4100,4384,4406,4407,4421,4422,4481,4482,4483,4564,4566,4569,4570,4573

  #change permissions of my stuff on NeSI 4 = read 2 = write 1 = execute
chmod 770 /nesi/nobackup/uoo00032/fastqfiles/Genomic/

#look for any job submit files in this directory
ls -alh /nesi/nobackup/uoo00032/fastqfiles/Genomic/*/*_jobSubmit*

#set up alignment on NeSI - start with GRCh37
/nesi/nobackup/uoo00032/fastqfiles/Genomic/4384/4384_jobSubmit_GRCh37.sh

#look at squeue in a nice way
squeue -u $USER -o "%.10a %.9u %.70j %12i %10A %.8T %.10M %.17R %.4P %.2C %.10l %.30E"
squeue -u $USER -o "%.10a %.13u %.70j %12i %10A %.8T %.10M %.17R %.4P %.2C %.10l %.30E" | grep "RUNNING\|USER"

#quick way to look at reports of jobs running/completed/failed/uh oh sex problem! - run this at the end of jobs to help with what to do next
/nesi/project/uoo00032/Resources/bin/GRCh37/Fastq2VCFplus/CheckAlignments.sh

#De novo filter - genotype filtering
bcftools filter -i 'GT[0]="alt" && GT[1]="ref" && GT[2]="ref"' 4458fam_Genomic_R_220802_TIMMOR1_DNA_M001-R_220707_TIMMOR1_DNA_M001_Split_ann_AF0.001_nostar_HIGHMOD.vcf.gz -Oz -o 4458fam_Genomic_R_220802_TIMMOR1_DNA_M001-R_220707_TIMMOR1_DNA_M001_Split_ann_AF0.001_nostar_HIGHMOD_denovo.vcf.gz

#Compound Het Filtering
bcftools filter -i 'GT[0]="het" && GT[1]="het" && GT[2]="het"' /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2560_2561_20230106_hg38_split_ann.vcf.gz  FILE 1

bcftools filter -i 'GT[0]="ref" && GT[1]="het" && GT[2]="het"' /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2560_2561_20230106_hg38_split_ann.vcf.gz

bcftools +split-vep /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2560_2561_20230106_hg38_split_ann_CompHet1.vcf.gz -f '%CHROM:%POS %Gene\n' -Oz -o /mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2560_2561_20230106_hg38_split_ann_CompHet1_Genelist.txt

bcftools filter -i '(gnomADgenomes_AF[*] <= '0.001' | gnomADgenomes_AF[*] = ".") & (gnomADexomes_AF[*] <= '0.001' | gnomADexomes_AF[*] = ".")' 4458fam_Genomic_R_220802_TIMMOR1_DNA_M001-R_220707_TIMMOR1_DNA_M001_Split_ann.vcf.gz -Oz -o 4458fam_Genomic_R_220802_TIMMOR1_DNA_M001-R_220707_TIMMOR1_DNA_M001_Split_ann_AF0.001.vcf.gz

#for loop to make sample map >> amend 
for file in $(ls -d /mnt/hcs/WCHP_Clinical_Genetics/SequenceData/GRCh37/Alignments/Genomic/*/*HC.g.vcf.gz); do id=$(basename $(dirname $file)); echo -e "$id\t$file" >> SampleMap_1.txt; done

#steps after aligning on NeSI
I should be in tomorrow morning, but you can make a start if you want. I checked the GenomicsDBWorkspaces for GRCh37 and it looks like we added all samples from the previous batch to a database, so we can now just add the current samples, before we do a joint call.
So you need a "samplemap" which is a file that specifies the sample ID and then the full path to the haplotypecaller g.vcf.gz file that contains the call for that sample. It is a tab-separated file, with one sample per line. You can write that manually using a text editor, or here's some code that does something a bit more useful - it finds all of the available *HC.g.vcf.gz files in any specified alignment directory, compares this to the samples already in the database, and only writes entries to your samplemap file if they are not already in the database. This way you get only the samples you want, and ALL of the samples you want without too much thought. Just run the following commands (in your home directory on the cluster)- the first three are just specifying directories/filename etc, the real business is in the last one (and the # lines are just comments - you don't need to run them)
# here I am just giving a name to the samplemap file - I using the batch name of (most of) the samples that was provided by Garvan
samplemap=R_221021_TIMMOR1_DNA_M001.txt
# here I am just specifying the directory that you have your alignment directories in 
aligndir=/mnt/hcs/WCHP_Clinical_Genetics/SequenceData/GRCh37/Alignments/Genomic
# here I am specifying the GenomicsDBWorkspace that contains the previously entered 
workspace=/mnt/hcs/WCHP_Clinical_Genetics/SequenceData/GRCh37/GenomicsDBWorkspaces/Genomic_ungapped_contigs
# then run this command which picks up the values from the variables above. Make sure you run it as a single line (even if the email breaks it over more than one line.
#and next time you need to make a samplemap you can just specify the values you want and do it again
for i in $aligndir/*/*_HC.g.vcf.gz; do sampledir=$(dirname ${i});if ! grep "^$(basename ${sampledir})" $workspace/allsamples.txt > /dev/null 2>&1; then echo -e "$(basename ${sampledir})\t ${sampledir}/$(basename ${sampledir})_HC.g.vcf.gz" >> $samplemap;fi; done

That should produce a file in your home directory called R_221021_TIMMOR1_DNA_M001.txt, which is your samplemap.
You then need to run the script that does the import into the database:
/resource/pipelines/Import_to_GenomicDB/Import.sh
It will ask you a couple of questions:
- GenomicsDBWorkspace (which is the workspace shown above - but you do need to enter it)
- samplemap - ie the full path to your samplemap file which will be something like /home/wadem43p/R_221021_TIMMOR1_DNA_M001.txt
- Assembly which will be GRCh37
