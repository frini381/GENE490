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

dir <- paste0("/mnt/hcs/WCHP_Clinical_Genetics/Niamh/MET") #set woring directory

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





#################################################################################################################

    t_file=paste0(dir,"/4584_4585_4586_20230417_GRCh37_split_ann_A0.01_HIGHMOD_DPGQ_DS_out.tsv")

   

        data <- fread(file = t_file, sep = '\t', header = TRUE,colClasses = c("character")) #set paternal variants as a variable in R


    Node_file=paste0(dir, "/string_protein_annotations.tsv")
		
		node_data <- fread(file = Node_file, sep = '\t', header = TRUE,colClasses = c("character")) #set paternal variants as a variable in R


MET_genes <- intersect(data$SYMBOL,node_data$node) #overlapped gene ids between both maternal and paternal

data_MET_genes <- data %>% filter(Gene %in% node_data)

  write.table(data_MET_genes, file=paste0(dir, "/4584_4585_4586_20230417_GRCh37_split_ann_AF0.01_HIGHMOD_DPGQ_DS_MET_out.tsv"), quote=FALSE, sep='\t', row.names=F)

#################################################################################################################


    t2_file=paste0(dir,"/2559_2562_2560_2561_20230707_hg38_Split_ann_A0.01_HIGHMOD_DPGQ_DS_out.tsv")

   

        data2 <- fread(file = t2_file, sep = '\t', header = TRUE,colClasses = c("character")) #set paternal variants as a variable in R


    Node2_file=paste0(dir, "/string_protein_annotations.tsv")
		
		node2_data <- fread(file = Node_file, sep = '\t', header = TRUE,colClasses = c("character")) #set paternal variants as a variable in R


MET2_genes <- intersect(data$SYMBOL,node_data$node) #overlapped gene ids between both maternal and paternal

data_MET2_genes <- data2 %>% filter(Gene %in% MET2_genes)

 write.table(data_MET2_genes, file=paste0(dir, "/2559_2562_2560_2561_20230707_hg38_Split_ann_A0.01_HIGHMOD_DPGQ_DS_MET_out.tsv"), quote=FALSE, sep='\t', row.names=F)