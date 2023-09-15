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

4584_4585_4586/4584_4585_4586_20230417_GRCh37_split_ann_AF0.01_HIGHMOD_DPGQ_out.tsv


dir <- paste0("/mnt/hcs/WCHP_Clinical_Genetics/Niamh/4584_4585_4586") #set woring directory

analysisDir <- paste0(dir,"/R_Analysis") 

Dirs <- c("data","graphs")

ifelse(!dir.exists(file.path(analysisDir)), dir.create(file.path(analysisDir)),TRUE)

setwd(analysisDir)

for(subDir in Dirs) {

    ifelse(!dir.exists(file.path(analysisDir, subDir)), dir.create(file.path(analysisDir, subDir)), TRUE)

}

print(paste0("Directory Set to ", analysisDir))


# gnomAD Constraint

file.copy("/mnt/hcs/WCHP_Clinical_Genetics/Ben/Resources/gnomAD/gnomad.v2.1.1.lof_metrics.by_transcript.txt","data/gnomAD_Transcript_Constraint.txt")

gnomAD <-read.table( file = "data/gnomAD_Transcript_Constraint.txt",header = TRUE,sep = "\t") 

#############################################################################################



    pat_file=paste0(dir,"/4584_4585_4586_20230417_GRCh37_split_ann_AF0.01_HIGHMOD_CompHet1_out.tsv")

   

        pat_data <- fread(file = pat_file, sep = '\t', header = TRUE,colClasses = c("character")) #set paternal variants as a variable in R


  
 


    # Split Arrays

    column_names <- c("Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature","BIOTYPE","EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position","Protein_position","Amino_acids","Codons","Existing_variation","DISTANCE","STRAND","FLAGS","SYMBOL_SOURCE","HGNC_ID","CCDS","ENSP","REFSEQ_MATCH","SOURCE","REFSEQ_OFFSET","GIVEN_REF","USED_REF","BAM_EDIT","GENE_PHENO","HGVS_OFFSET","HGVSg","CLIN_SIG","SOMATIC","PHENO","MOTIF_NAME","MOTIF_POS","HIGH_INF_POS","MOTIF_SCORE_CHANGE","TRANSCRIPTION_FACTORS","existing_InFrame_oORFs","existing_OutOfFrame_oORFs","existing_uORFs","five_prime_UTR_variant_annotation","five_prime_UTR_variant_consequence") #CSQ variable names for each separate column

    pat_data_Split <- pat_data %>%  separate_rows(CSQ,sep=",") %>% tidyr::separate(col=CSQ, sep ='\\|',into=column_names, convert = TRUE, remove = FALSE, fill = "right") %>% select(-CSQ) %>% filter(SOURCE == "Ensembl") #separate each CSQ data element into each column and remove CSQ and ANN columns and filter so only Emsemble gene ids

     write.table(pat_data_Split, file=paste0(dir, "/R_Analysis/data/4584_4585_4586_20230417_GRCh37_split_ann_AF0.01_HIGHMOD_CompHet1_DS_out.tsv"), quote=FALSE, sep='\t', row.names=F) #convert split paternal data into a tsv



select()



mat_file=paste0(dir,"/4584_4585_4586_20230417_GRCh37_split_ann_AF0.01_HIGHMOD_CompHet2_out.tsv") 
   

        mat_data <- fread(file = mat_file, sep = '\t', header = TRUE,colClasses = c("character")) #materal variants as variable in R

 
    # Split Arrays

   mat_data_Split <- mat_data %>%  separate_rows(CSQ,sep=",") %>% tidyr::separate(col=CSQ, sep ='\\|',into=column_names, convert = TRUE, remove = FALSE, fill = "right") %>% select(-CSQ) %>% filter(SOURCE == "Ensembl")  #separate each CSQ data element into each column and remove CSQ and ANN columns and filter so only Emsemble gene ids


     write.table(mat_data_Split, file=paste0(dir, "/R_Analysis/data/4584_4585_4586_20230417_GRCh37_split_ann_AF0.01_HIGHMOD_CompHet2_DS_out.tsv"), quote=FALSE, sep='\t', row.names=F) #convert split maternal data into a tsv




comphetgenes <- intersect(mat_data_Split$Gene,pat_data_Split$Gene) #overlapped gene ids between both maternal and paternal

mat_data_Split_comphetgenes <- mat_data_Split %>% filter(Gene %in% comphetgenes) ##remove the non-shared genes in each mat/pat object
pat_data_Split_comphetgenes <- pat_data_Split %>% filter(Gene %in% comphetgenes) 

data_Split_comphetgenes <- rbind(mat_data_Split,pat_data_Split) %>% arrange(Gene) %>% filter(BIOTYPE=="protein_coding") %>% filter(IMPACT!="MODIFIER") ## filtering for protein coding and removing modifier variants and combining mat/pat variants into single object

    # Add gnomAD columns 

    data_Split_comphetgenes$gnomAD_Transcript_Canonical <- gnomAD$canonical[match(data_Split_comphetgenes$Feature,gnomAD$transcript)] 

   data_Split_comphetgenes_gnomAD <- data_Split_comphetgenes %>% filter(gnomAD_Transcript_Canonical == "true") #only canonical transcripts 

    column_names2 <- sprintf("SpliceAI_%s", c("ALLELE","SYMBOL","DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL","DP_DG" ,"DP_DL"))

    data_Split_gnomAD_SpliceAI <- data_Split_comphetgenes_gnomAD %>% separate_rows(SpliceAI,sep=",") %>% tidyr::separate(col=SpliceAI, sep ='\\|',into=column_names2, convert = TRUE, remove = FALSE, fill = "right") %>% select(-SpliceAI)  

write.table(data_Split_gnomAD_SpliceAI, file=paste0(dir, "/R_Analysis/data/4584_4585_4586_20230417_GRCh37_split_ann_AF0.01_HIGHMOD_DS_AD_out.tsv"), quote=FALSE, sep='\t', row.names=F)

#tsv of variants that are in shared genes between mat and pat