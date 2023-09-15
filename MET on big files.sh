MET on big files


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

print(paste0("Modules loaded"))
 

####################################################################################################################
###                                                     Setup                                                    ###
####################################################################################################################

dir <- paste0("/mnt/hcs/WCHP_Clinical_Genetics/Niamh/4584_4585_4586") #set woring directory

analysisDir <- paste0(dir,"/R_Analysis") 

Dirs <- c("data","graphs")

ifelse(!dir.exists(file.path(analysisDir)), dir.create(file.path(analysisDir)),TRUE)

setwd(analysisDir)

for(subDir in Dirs) {

    ifelse(!dir.exists(file.path(analysisDir, subDir)), dir.create(file.path(analysisDir, subDir)), TRUE)

}

print(paste0("Directory Set to ", analysisDir))


file.copy("/mnt/hcs/WCHP_Clinical_Genetics/Ben/Resources/gnomAD/gnomad.v2.1.1.lof_metrics.by_transcript.txt","data/gnomAD_Transcript_Constraint.txt")

gnomAD <-read.table( file = "data/gnomAD_Transcript_Constraint.txt",header = TRUE,sep = "\t") 


###


trio_file=paste0(dir,"/4584_4585_4586_20230417_GRCh37_split_ann_out.tsv")

		trio_data <- fread(file = trio_file, sep = '\t', header = TRUE,colClasses = c("character"))

    
	column_names <- c("Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature","BIOTYPE","EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position","Protein_position","Amino_acids","Codons","Existing_variation","DISTANCE","STRAND","FLAGS","SYMBOL_SOURCE","HGNC_ID","CCDS","ENSP","REFSEQ_MATCH","SOURCE","REFSEQ_OFFSET","GIVEN_REF","USED_REF","BAM_EDIT","GENE_PHENO","HGVS_OFFSET","HGVSg","CLIN_SIG","SOMATIC","PHENO","MOTIF_NAME","MOTIF_POS","HIGH_INF_POS","MOTIF_SCORE_CHANGE","TRANSCRIPTION_FACTORS","existing_InFrame_oORFs","existing_OutOfFrame_oORFs","existing_uORFs","five_prime_UTR_variant_annotation","five_prime_UTR_variant_consequence")

	trio_data_split <- trio_data %>% separate_rows(CSQ,sep=",") %>% tidyr::separate(col=CSQ, sep ='\\|',into=column_names, convert = TRUE, remove = FALSE, fill = "right") %>% select(-ANN,-CSQ) %>% filter(SOURCE == "Ensembl")

	write.table(trio_data_split, file=paste0(dir, "/mnt/hcs/WCHP_Clinical_Genetics/Niamh/MET/4584_4585_4586_20230417_GRCh37_split_ann_DS_out.tsv"))


    Node_file=paste0(dir, "/string_protein_annotations.tsv")
		
		node_data <- fread(file = Node_file, sep = '\t', header = TRUE,colClasses = c("character")) #set paternal variants as a variable in R


	MET_genes <- intersect(trio_data_split$SYMBOL,node_data$node) #overlapped gene ids between both maternal and paternal

	data_MET_genes1 <- trio_data_split %>% filter(Gene %in% node_data)

  	write.table(data_MET_genes1, file=paste0(dir, "/mnt/hcs/WCHP_Clinical_Genetics/Niamh/MET/4584_4585_4586_20230417_GRCh37_split_ann_DS_MET_out.tsv"), quote=FALSE, sep='\t', row.names=F)



####################################################################################################################
###                                                     Setup                                                    ###
####################################################################################################################


dir <- paste0("/mnt/hcs/WCHP_Clinical_Genetics/Niamh/2559_2562_2560_2561") #set woring directory

analysisDir <- paste0(dir,"/R_Analysis") 

Dirs <- c("data","graphs")

ifelse(!dir.exists(file.path(analysisDir)), dir.create(file.path(analysisDir)),TRUE)

setwd(analysisDir)

for(subDir in Dirs) {

    ifelse(!dir.exists(file.path(analysisDir, subDir)), dir.create(file.path(analysisDir, subDir)), TRUE)

}

print(paste0("Directory Set to ", analysisDir))


file.copy("/mnt/hcs/WCHP_Clinical_Genetics/Ben/Resources/gnomAD/gnomad.v2.1.1.lof_metrics.by_transcript.txt","data/gnomAD_Transcript_Constraint.txt")

gnomAD <-read.table( file = "data/gnomAD_Transcript_Constraint.txt",header = TRUE,sep = "\t") 


###

	quart_file=paste0(dir,"/2559_2562_2560_2561_20230707_hg38_Split_ann_out.tsv")

		quart_data <- fread(file = quart_file, sep = '\t', header = TRUE,colClasses = c("character"))

    
	column_names <- c("Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature","BIOTYPE","EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position","Protein_position","Amino_acids","Codons","Existing_variation","DISTANCE","STRAND","FLAGS","SYMBOL_SOURCE","HGNC_ID","CCDS","ENSP","REFSEQ_MATCH","SOURCE","REFSEQ_OFFSET","GIVEN_REF","USED_REF","BAM_EDIT","GENE_PHENO","HGVS_OFFSET","HGVSg","CLIN_SIG","SOMATIC","PHENO","MOTIF_NAME","MOTIF_POS","HIGH_INF_POS","MOTIF_SCORE_CHANGE","TRANSCRIPTION_FACTORS","existing_InFrame_oORFs","existing_OutOfFrame_oORFs","existing_uORFs","five_prime_UTR_variant_annotation","five_prime_UTR_variant_consequence")

	quart_data_split <- quart_datat %>% separate_rows(CSQ,sep=",") %>% tidyr::separate(col=CSQ, sep ='\\|',into=column_names, convert = TRUE, remove = FALSE, fill = "right") %>% select(-CSQ) %>% filter(SOURCE == "Ensembl")

	write.table(quart_data_split, file=paste0(dir, "/mnt/hcs/WCHP_Clinical_Genetics/Niamh/MET/2559_2562_2560_2561_20230707_hg38_Split_DS_out.tsv"))


    Node_file=paste0(dir, "/string_protein_annotations.tsv")
		
		node_data <- fread(file = Node_file, sep = '\t', header = TRUE,colClasses = c("character")) #set paternal variants as a variable in R


	MET_genes <- intersect(quart_data_split$SYMBOL,node_data$node) #overlapped gene ids between both maternal and paternal

	data_MET_genes2 <- quart_data_split %>% filter(Gene %in% node_data)

  	write.table(data_MET_genes2, file=paste0(dir, "/mnt/hcs/WCHP_Clinical_Genetics/Niamh/MET/2559_2562_2560_2561_20230707_hg38_Split_DS_MET_out.tsv"), quote=FALSE, sep='\t', row.names=F)


