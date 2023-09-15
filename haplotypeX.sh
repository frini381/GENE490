module purge; module load R/3.6.0-foss-2019a; R

 

### Start up ####

if(!require(tidyverse)){

    install.packages("tidyverse")

    library(tidyverse)

}

if(!require(data.table)){

    install.packages("data.table")

    library(data.table)

}

if(!require(viridis)){

    install.packages("viridis")

    library(viridis)

}

if(!require(reshape2)){

    install.packages("reshape2")

    library(reshape2)

}

 

####################################################################################################################
###                                                     Setup                                                    ###
####################################################################################################################

 

# gnomAD Constraint
file.copy("/mnt/hcs/WCHP_Clinical_Genetics/Ben/Resources/gnomAD/gnomad.v2.1.1.lof_metrics.by_transcript.txt","data/gnomAD_Transcript_Constraint.txt")

gnomAD <-read.table( file = "/mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2562_2560_2561/R_Analysis/data/gnomAD_Transcript_Constraint.txt",header = TRUE,sep = "\t") 


dir <- paste0("/mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2562_2560_2561") #set woring directory

analysisDir <- paste0(dir,"/R_Analysis") 

Dirs <- c("data","graphs")

ifelse(!dir.exists(file.path(analysisDir)), dir.create(file.path(analysisDir)),TRUE)

setwd(analysisDir)

for(subDir in Dirs) {

    ifelse(!dir.exists(file.path(analysisDir, subDir)), dir.create(file.path(analysisDir, subDir)), TRUE)

}

print(paste0("Directory Set to ", analysisDir))


#############################################################################################



    not_shared_file_1=paste0(dir,"/R_Analysis/data/2559_2562_2560_2561_20230707_hg38_Split_ann_ChrX_NShared1_out.tsv")

   

        NS1_data <- fread(file = not_shared_file_1, sep = '\t', header = TRUE,colClasses = c("character")) #set variants as a variable in R


  




    not_shared_file_2=paste0(dir,"/R_Analysis/data/2559_2562_2560_2561_20230707_hg38_Split_ann_ChrX_NShared2_out.tsv")
   

        NS2_data <- fread(file = not_shared_file_2, sep = '\t', header = TRUE,colClasses = c("character")) #set variants as a variable in R

 



    shared_file=paste0(dir,"/2559_2562_2560_2561_20230707_hg38_Split_ann_ChrX_Shared_out.tsv")
   

        S_data <- fread(file = not_shared_file_2, sep = '\t', header = TRUE,colClasses = c("character")) #


data_NS <- rbind(NS1_data,NS2_data) ## combining variants into single object

    # Add gnomAD columns 

column_names <- c("Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature","BIOTYPE","EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position","Protein_position","Amino_acids","Codons","Existing_variation","DISTANCE","STRAND","FLAGS","SYMBOL_SOURCE","HGNC_ID","CCDS","ENSP","REFSEQ_MATCH","SOURCE","REFSEQ_OFFSET","GIVEN_REF","USED_REF","BAM_EDIT","GENE_PHENO","HGVS_OFFSET","HGVSg","CLIN_SIG","SOMATIC","PHENO","MOTIF_NAME","MOTIF_POS","HIGH_INF_POS","MOTIF_SCORE_CHANGE","TRANSCRIPTION_FACTORS","existing_InFrame_oORFs","existing_OutOfFrame_oORFs","existing_uORFs","five_prime_UTR_variant_annotation","five_prime_UTR_variant_consequence") #CSQ variable names for each separate column

    data_NS_Split <- data_NS %>%  separate_rows(CSQ,sep=",") %>% tidyr::separate(col=CSQ, sep ='\\|',into=column_names, convert = TRUE, remove = FALSE, fill = "right") %>% select(-CSQ)

    data_S_Split <- S_data %>%  separate_rows(CSQ,sep=",") %>% tidyr::separate(col=CSQ, sep ='\\|',into=column_names, convert = TRUE, remove = FALSE, fill = "right") %>% select(-CSQ)


    data_NS_Split$gnomAD_Transcript_Canonical <- gnomAD$canonical[match(data_NS_Split$Feature,gnomAD$transcript)] 

    data_S_Split$gnomAD_Transcript_Canonical <- gnomAD$canonical[match(data_S_Split$Feature,gnomAD$transcript)] 


   data_NS_Split_gnomAD <- data_NS_Split %>% filter(gnomAD_Transcript_Canonical == "true") #only canonical transcripts 

   data_S_Split_gnomAD <- data_S_Split %>% filter(gnomAD_Transcript_Canonical == "true") #only canonical transcripts 


# Split SpliceAI Field

 column_names2 <- sprintf("SpliceAI_%s", c("ALLELE","SYMBOL","DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL","DP_DG" ,"DP_DL"))

    
    data_NS_Split_gnomAD_SpliceAI <- data_NS_Split_gnomAD %>% separate_rows(SpliceAI,sep=",") %>% tidyr::separate(col=SpliceAI, sep ='\\|',into=column_names2, convert = TRUE, remove = FALSE, fill = "right") %>% select(-SpliceAI) 


    data_S_Split_gnomAD_SpliceAI <- data_S_Split_gnomAD %>% separate_rows(SpliceAI,sep=",") %>% tidyr::separate(col=SpliceAI, sep ='\\|',into=column_names2, convert = TRUE, remove = FALSE, fill = "right") %>% select(-SpliceAI)  



write.table(data_NS_Split_gnomAD_SpliceAI, file=paste0(dir, "/Haplotypes/2559_2562_2560_2561_20230707_hg38_Split_ann_ChrX_Split_NotShared_out.tsv"), quote=FALSE, sep='\t', row.names=F)


write.table(data_S_Split_gnomAD_SpliceAI, file=paste0(dir, "/Haplotypes/2559_2562_2560_2561_20230707_hg38_Split_ann_ChrX_Split_Shared_out.tsv"), quote=FALSE, sep='\t', row.names=F)

#tsv of variants that are in shared genes between mat and pat


