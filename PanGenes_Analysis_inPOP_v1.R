###########################################################################

# Script Comparing genes (Core-Variable) between population
# IM 24.11.17

#############################################################################




library(ggplot2)
require(scales)
library(RColorBrewer)
library(UpSetR)

#################################################################################################################
#

#########################################################################
######### CALCULATION OF Corge genes and variable genes per POP #########
#########################################################################


######### CALIFORNIA ############
#################################

setwd("~/Documents/EPFL/Vibrio_Chromosomal_rearrangements/PopGenomics_V2/Pop_California/CoreGenes/PanGenes/")
filesToProcess <- dir(pattern = "*\\.txt")  #files to process.
listOfFiles <- lapply(filesToProcess, function(x) tryCatch(read.table(x, header = F),colClasses ="NULL",
                                                           error= function (e) cbind.data.frame(V1="NA")))
names(listOfFiles)<-gsub(".txt","",filesToProcess)

listOfFiles<-lapply(listOfFiles,function(x) x[,1])

head(listOfFiles[[1]])
names(listOfFiles)

# Extract Pan genes for each sample
NonConserved_Cali<-list()
for  (i in 1: length(unique(gsub("Core_","",gsub("ALL_","",names(listOfFiles)) )))    ) {
  k<-unique(gsub("Core_","",gsub("ALL_","",names(listOfFiles)) ))[i]
NonConserved_Cali[[i]]<-listOfFiles[grep(paste(k,"(?=$)",sep = ""),names(listOfFiles), perl=T)][[1]]  [!listOfFiles[grep(paste(k,"(?=$)",sep = ""),names(listOfFiles), perl=T)][[1]]  %in%  listOfFiles[grep(paste(k,"(?=$)",sep = ""),names(listOfFiles), perl=T)][[2]]    ]
names(NonConserved_Cali)[i]<- paste("NonConserved",k,sep = "_")
}

#########################################################################
#### PLOTING CORE AND PAN GENES
CoreGenes<-listOfFiles[grep("Core",names(listOfFiles))]
Core_Cali<-CoreGenes
CoreGenes<-lapply(CoreGenes,function (x) length(x))
NonConserved<-lapply(NonConserved_Cali,function (x) length(x))
NonConserved<-data.frame(matrix(unlist(NonConserved), ncol = length(names(NonConserved)), byrow=T))
CoreGenes<-data.frame(matrix(unlist(CoreGenes), ncol = length(names(CoreGenes)), byrow=T))

df2<-as.data.frame(t(rbind.data.frame(rep(gsub("Core_","",names(listOfFiles)[grep("Core",names(listOfFiles))]),2),
                     as.numeric(cbind(NonConserved,CoreGenes)),c(rep("NonConserved",length(names(NonConserved))),rep("Core",length(names(NonConserved))) ) ))  )
df2$V2<-as.numeric(levels(df2$V2))[df2$V2]
df2$V3<-factor(df2$V3, levels=c("NonConserved","Core"))

pdf("Core_Pan_genes_CALIFORNIA.pdf",height = 8, width = 10)
ggplot(data=df2, aes(x=V1, y=V2, fill=V3)) + labs(fill = "Genes", y="Number of genes")+ scale_y_continuous(limits = c(0, 4300)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),axis.title.x=element_blank(),
                                    axis.text.y = element_text(size=12),axis.title.y= element_text(size=14)
                                    ) + scale_fill_grey() 
dev.off()
#########################################################################

######### MEXICO ############
#################################

setwd("~/Documents/EPFL/Vibrio_Chromosomal_rearrangements/PopGenomics_V2/Pop_Mexico/CoreGenes/PanGenes/")
filesToProcess <- dir(pattern = "*\\.txt")  #files to process.
listOfFiles <- lapply(filesToProcess, function(x) tryCatch(read.table(x, header = F),colClasses ="NULL",
                                                           error= function (e) cbind.data.frame(V1="NA")))
names(listOfFiles)<-gsub(".txt","",filesToProcess)

listOfFiles<-lapply(listOfFiles,function(x) x[,1])

head(listOfFiles[[1]])
names(listOfFiles)

# Extract Pan genes for each sample
NonConserved_Mex<-list()
for  (i in 1: length(unique(gsub("Core_","",gsub("ALL_","",names(listOfFiles)) )))    ) {
  k<-unique(gsub("Core_","",gsub("ALL_","",names(listOfFiles)) ))[i]
  NonConserved_Mex[[i]]<-listOfFiles[grep(paste(k,"(?=$)",sep = ""),names(listOfFiles), perl=T)][[1]]  [!listOfFiles[grep(paste(k,"(?=$)",sep = ""),names(listOfFiles), perl=T)][[1]]  %in%  listOfFiles[grep(paste(k,"(?=$)",sep = ""),names(listOfFiles), perl=T)][[2]]    ]
  names(NonConserved_Mex)[i]<- paste("NonConserved",k,sep = "_")
}

#########################################################################
#### PLOTING CORE AND PAN GENES
CoreGenes<-listOfFiles[grep("Core",names(listOfFiles))]
Core_Mex<-CoreGenes

CoreGenes<-lapply(CoreGenes,function (x) length(x))
NonConserved<-lapply(NonConserved_Mex,function (x) length(x))
NonConserved<-data.frame(matrix(unlist(NonConserved), ncol = length(names(NonConserved)), byrow=T))
CoreGenes<-data.frame(matrix(unlist(CoreGenes), ncol = length(names(CoreGenes)), byrow=T))

df2<-as.data.frame(t(rbind.data.frame(rep(gsub("Core_","",names(listOfFiles)[grep("Core",names(listOfFiles))]),2),
                                      as.numeric(cbind(NonConserved,CoreGenes)),c(rep("NonConserved",length(names(NonConserved))),rep("Core",length(names(NonConserved))) ) ))  )
df2$V2<-as.numeric(levels(df2$V2))[df2$V2]
df2$V3<-factor(df2$V3, levels=c("NonConserved","Core"))

pdf("Core_Pan_genes_MEXICO.pdf",height = 8, width = 10)
ggplot(data=df2, aes(x=V1, y=V2, fill=V3)) + labs(fill = "Genes", y="Number of genes")+ scale_y_continuous(limits = c(0, 4300)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),axis.title.x=element_blank(),
                                    axis.text.y = element_text(size=12),axis.title.y= element_text(size=14)
  ) + scale_fill_grey() 
dev.off()
#########################################################################

######### Massachusetts ############
#################################

setwd("~/Documents/EPFL/Vibrio_Chromosomal_rearrangements/PopGenomics_V2/Pop_Massachusetts/CoreGenes/PanGenes/")
filesToProcess <- dir(pattern = "*\\.txt")  #files to process.
listOfFiles <- lapply(filesToProcess, function(x) tryCatch(read.table(x, header = F),colClasses ="NULL",
                                                           error= function (e) cbind.data.frame(V1="NA")))
names(listOfFiles)<-gsub(".txt","",filesToProcess)

listOfFiles<-lapply(listOfFiles,function(x) x[,1])

head(listOfFiles[[1]])
names(listOfFiles)

# Extract Pan genes for each sample
NonConserved_Masa<-list()
for  (i in 1: length(unique(gsub("Core_","",gsub("ALL_","",names(listOfFiles)) )))    ) {
  k<-unique(gsub("Core_","",gsub("ALL_","",names(listOfFiles)) ))[i]
  NonConserved_Masa[[i]]<-listOfFiles[grep(paste(k,"(?=$)",sep = ""),names(listOfFiles), perl=T)][[1]]  [!listOfFiles[grep(paste(k,"(?=$)",sep = ""),names(listOfFiles), perl=T)][[1]]  %in%  listOfFiles[grep(paste(k,"(?=$)",sep = ""),names(listOfFiles), perl=T)][[2]]    ]
  names(NonConserved_Masa)[i]<- paste("NonConserved",k,sep = "_")
}

#########################################################################
#### PLOTING CORE AND PAN GENES
CoreGenes<-listOfFiles[grep("Core",names(listOfFiles))]
Core_Masa<-CoreGenes

CoreGenes<-lapply(CoreGenes,function (x) length(x))
NonConserved<-lapply(NonConserved_Masa,function (x) length(x))
NonConserved<-data.frame(matrix(unlist(NonConserved), ncol = length(names(NonConserved)), byrow=T))
CoreGenes<-data.frame(matrix(unlist(CoreGenes), ncol = length(names(CoreGenes)), byrow=T))

df2<-as.data.frame(t(rbind.data.frame(rep(gsub("Core_","",names(listOfFiles)[grep("Core",names(listOfFiles))]),2),
                                      as.numeric(cbind(NonConserved,CoreGenes)),c(rep("NonConserved",length(names(NonConserved))),rep("Core",length(names(NonConserved))) ) ))  )
df2$V2<-as.numeric(levels(df2$V2))[df2$V2]
df2$V3<-factor(df2$V3, levels=c("NonConserved","Core"))

pdf("Core_Pan_genes_MASSACHUSETTS.pdf",height = 8, width = 10)
ggplot(data=df2, aes(x=V1, y=V2, fill=V3)) + labs(fill = "Genes", y="Number of genes")+ scale_y_continuous(limits = c(0, 4300)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),axis.title.x=element_blank(),
                                    axis.text.y = element_text(size=12),axis.title.y= element_text(size=14)
  ) + scale_fill_grey() 
dev.off()
#########################################################################


######### Thailand ############
#################################

setwd("~/Documents/EPFL/Vibrio_Chromosomal_rearrangements/PopGenomics_V2/Pop_Thailand/CoreGenes/PanGenes/")
filesToProcess <- dir(pattern = "*\\.txt")  #files to process.
listOfFiles <- lapply(filesToProcess, function(x) tryCatch(read.table(x, header = F),colClasses ="NULL",
                                                           error= function (e) cbind.data.frame(V1="NA")))
names(listOfFiles)<-gsub(".txt","",filesToProcess)

listOfFiles<-lapply(listOfFiles,function(x) x[,1])

head(listOfFiles[[1]])
names(listOfFiles)

# Extract Pan genes for each sample
NonConserved_Thai<-list()
for  (i in 1: length(unique(gsub("Core_","",gsub("ALL_","",names(listOfFiles)) )))    ) {
  k<-unique(gsub("Core_","",gsub("ALL_","",names(listOfFiles)) ))[i]
  NonConserved_Thai[[i]]<-listOfFiles[grep(paste(k,"(?=$)",sep = ""),names(listOfFiles), perl=T)][[1]]  [!listOfFiles[grep(paste(k,"(?=$)",sep = ""),names(listOfFiles), perl=T)][[1]]  %in%  listOfFiles[grep(paste(k,"(?=$)",sep = ""),names(listOfFiles), perl=T)][[2]]    ]
  names(NonConserved_Thai)[i]<- paste("NonConserved",k,sep = "_")
}

#########################################################################
#### PLOTING CORE AND PAN GENES
CoreGenes<-listOfFiles[grep("Core",names(listOfFiles))]
Core_Thai<-CoreGenes

CoreGenes<-lapply(CoreGenes,function (x) length(x))
NonConserved<-lapply(NonConserved_Thai,function (x) length(x))
NonConserved<-data.frame(matrix(unlist(NonConserved), ncol = length(names(NonConserved)), byrow=T))
CoreGenes<-data.frame(matrix(unlist(CoreGenes), ncol = length(names(CoreGenes)), byrow=T))

df2<-as.data.frame(t(rbind.data.frame(rep(gsub("Core_","",names(listOfFiles)[grep("Core",names(listOfFiles))]),2),
                                      as.numeric(cbind(NonConserved,CoreGenes)),c(rep("NonConserved",length(names(NonConserved))),rep("Core",length(names(NonConserved))) ) ))  )
df2$V2<-as.numeric(levels(df2$V2))[df2$V2]
df2$V3<-factor(df2$V3, levels=c("NonConserved","Core"))

pdf("Core_Pan_genes_THAILAND.pdf",height = 8, width = 10)
ggplot(data=df2, aes(x=V1, y=V2, fill=V3)) + labs(fill = "Genes", y="Number of genes")+ scale_y_continuous(limits = c(0, 4300)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),axis.title.x=element_blank(),
                                    axis.text.y = element_text(size=12),axis.title.y= element_text(size=14)
  ) + scale_fill_grey() 
dev.off()
#########################################################################

######### HAITI ############
#################################

setwd("~/Documents/EPFL/Vibrio_Chromosomal_rearrangements/PopGenomics_V2/Pop_Haiti/CoreGenes/PanGenes/")
filesToProcess <- dir(pattern = "*\\.txt")  #files to process.
listOfFiles <- lapply(filesToProcess, function(x) tryCatch(read.table(x, header = F),colClasses ="NULL",
                                                           error= function (e) cbind.data.frame(V1="NA")))
names(listOfFiles)<-gsub(".txt","",filesToProcess)

listOfFiles<-lapply(listOfFiles,function(x) x[,1])

head(listOfFiles[[1]])
names(listOfFiles)

# Extract Pan genes for each sample
NonConserved_Haiti<-list()
for  (i in 1: length(unique(gsub("Core_","",gsub("ALL_","",names(listOfFiles)) )))    ) {
  k<-unique(gsub("Core_","",gsub("ALL_","",names(listOfFiles)) ))[i]
  NonConserved_Haiti[[i]]<-listOfFiles[grep(paste(k,"(?=$)",sep = ""),names(listOfFiles), perl=T)][[1]]  [!listOfFiles[grep(paste(k,"(?=$)",sep = ""),names(listOfFiles), perl=T)][[1]]  %in%  listOfFiles[grep(paste(k,"(?=$)",sep = ""),names(listOfFiles), perl=T)][[2]]    ]
  names(NonConserved_Haiti)[i]<- paste("NonConserved",k,sep = "_")
}

#########################################################################
#### PLOTING CORE AND PAN GENES
CoreGenes<-listOfFiles[grep("Core",names(listOfFiles))]
Core_Haiti<-CoreGenes

CoreGenes<-lapply(CoreGenes,function (x) length(x))
NonConserved<-lapply(NonConserved_Haiti,function (x) length(x))
NonConserved<-data.frame(matrix(unlist(NonConserved), ncol = length(names(NonConserved)), byrow=T))
CoreGenes<-data.frame(matrix(unlist(CoreGenes), ncol = length(names(CoreGenes)), byrow=T))

df2<-as.data.frame(t(rbind.data.frame(rep(gsub("Core_","",names(listOfFiles)[grep("Core",names(listOfFiles))]),2),
                                      as.numeric(cbind(NonConserved,CoreGenes)),c(rep("NonConserved",length(names(NonConserved))),rep("Core",length(names(NonConserved))) ) ))  )
df2$V2<-as.numeric(levels(df2$V2))[df2$V2]
df2$V3<-factor(df2$V3, levels=c("NonConserved","Core"))

pdf("Core_Pan_genes_HAITI.pdf",height = 8, width = 10)
ggplot(data=df2, aes(x=V1, y=V2, fill=V3)) + labs(fill = "Genes", y="Number of genes")+ scale_y_continuous(limits = c(0, 4300)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),axis.title.x=element_blank(),
                                    axis.text.y = element_text(size=12),axis.title.y= element_text(size=14)
  ) + scale_fill_grey() 
dev.off()
#########################################################################

######### PANDEMIC ############
#################################

setwd("~/Documents/EPFL/Vibrio_Chromosomal_rearrangements/PopGenomics_V2/Pandemic_PseudoPOP/CoreGenes/PanGenes/")
filesToProcess <- dir(pattern = "*\\.txt")  #files to process.
listOfFiles <- lapply(filesToProcess, function(x) tryCatch(read.table(x, header = F),colClasses ="NULL",
                                                           error= function (e) cbind.data.frame(V1="NA")))
names(listOfFiles)<-gsub(".txt","",filesToProcess)

listOfFiles<-lapply(listOfFiles,function(x) x[,1])

head(listOfFiles[[1]])
names(listOfFiles)

# Extract Pan genes for each sample
NonConserved_Pandemic<-list()
for  (i in 1: length(unique(gsub("Core_","",gsub("ALL_","",names(listOfFiles)) )))    ) {
  k<-unique(gsub("Core_","",gsub("ALL_","",names(listOfFiles)) ))[i]
  NonConserved_Pandemic[[i]]<-listOfFiles[grep(paste(k,"(?=$)",sep = ""),names(listOfFiles), perl=T)][[1]]  [!listOfFiles[grep(paste(k,"(?=$)",sep = ""),names(listOfFiles), perl=T)][[1]]  %in%  listOfFiles[grep(paste(k,"(?=$)",sep = ""),names(listOfFiles), perl=T)][[2]]    ]
  names(NonConserved_Pandemic)[i]<- paste("NonConserved",k,sep = "_")
}

#########################################################################
#### PLOTING CORE AND PAN GENES
CoreGenes<-listOfFiles[grep("Core",names(listOfFiles))]
Core_Pan<-CoreGenes

CoreGenes<-lapply(CoreGenes,function (x) length(x))
NonConserved<-lapply(NonConserved_Pandemic,function (x) length(x))
NonConserved<-data.frame(matrix(unlist(NonConserved), ncol = length(names(NonConserved)), byrow=T))
CoreGenes<-data.frame(matrix(unlist(CoreGenes), ncol = length(names(CoreGenes)), byrow=T))

df2<-as.data.frame(t(rbind.data.frame(rep(gsub("Core_","",names(listOfFiles)[grep("Core",names(listOfFiles))]),2),
                                      as.numeric(cbind(NonConserved,CoreGenes)),c(rep("NonConserved",length(names(NonConserved))),rep("Core",length(names(NonConserved))) ) ))  )
df2$V2<-as.numeric(levels(df2$V2))[df2$V2]
df2$V3<-factor(df2$V3, levels=c("NonConserved","Core"))

pdf("Core_Pan_genes_Pandemic.pdf",height = 8, width = 10)
ggplot(data=df2, aes(x=V1, y=V2, fill=V3)) + labs(fill = "Genes", y="Number of genes")+ scale_y_continuous(limits = c(0, 4300)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),axis.title.x=element_blank(),
                                    axis.text.y = element_text(size=12),axis.title.y= element_text(size=14)
  ) + scale_fill_grey() 
dev.off()
#########################################################################

######### ALL POPULATIONS ############
#################################
setwd("~/Documents/EPFL/Vibrio_Chromosomal_rearrangements/PopGenomics_V2/ALL_POPS/BetweenPOPS_VAR/PanGenes/")

filesToProcess <- dir(pattern = "*\\.txt")  #files to process.
listOfFiles <- lapply(filesToProcess, function(x) tryCatch(read.table(x, header = F),colClasses ="NULL",
                                                           error= function (e) cbind.data.frame(V1="NA")))
names(listOfFiles)<-gsub(".txt","",filesToProcess)

listOfFiles<-lapply(listOfFiles,function(x) x[,1])

head(listOfFiles[[1]])
names(listOfFiles)

# Extract Pan genes for each sample
NonConserved_ALL6POP<-list()
for  (i in 1: length(unique(gsub("Core_","",gsub("ALL_","",names(listOfFiles)) )))    ) {
  k<-unique(gsub("Core_","",gsub("ALL_","",names(listOfFiles)) ))[i]
  NonConserved_ALL6POP[[i]]<-listOfFiles[grep(paste(k,"(?=$)",sep = ""),names(listOfFiles), perl=T)][[1]]  [!listOfFiles[grep(paste(k,"(?=$)",sep = ""),names(listOfFiles), perl=T)][[1]]  %in%  listOfFiles[grep(paste(k,"(?=$)",sep = ""),names(listOfFiles), perl=T)][[2]]    ]
  names(NonConserved_ALL6POP)[i]<- paste("NonConserved",k,sep = "_")
}

#########################################################################
#### PLOTING CORE AND PAN GENES
CoreGenes<-listOfFiles[grep("Core",names(listOfFiles))]
#Core_Pan<-CoreGenes

CoreGenes<-lapply(CoreGenes,function (x) length(x))
NonConserved<-lapply(NonConserved_ALL6POP,function (x) length(x))
NonConserved<-data.frame(matrix(unlist(NonConserved), ncol = length(names(NonConserved)), byrow=T))
CoreGenes<-data.frame(matrix(unlist(CoreGenes), ncol = length(names(CoreGenes)), byrow=T))

df2<-as.data.frame(t(rbind.data.frame(rep(gsub("Core_","",names(listOfFiles)[grep("Core",names(listOfFiles))]),2),
                                      as.numeric(cbind(NonConserved,CoreGenes)),c(rep("NonConserved",length(names(NonConserved))),rep("Core",length(names(NonConserved))) ) ))  )
df2$V2<-as.numeric(levels(df2$V2))[df2$V2]
df2$V3<-factor(df2$V3, levels=c("NonConserved","Core"))
NonConserved_Btw_POPS_Nb <-df2[df2$V3=="NonConserved",-3] 

# order phylo
Populations<-read.csv("~/Documents/EPFL/Vibrio_Chromosomal_rearrangements/PopGenomics_V2/ALL_POPS/Phylogenies_ALL/Strains_Pop_vf.csv", h=T)
PopulationNew<-rep(Populations$strain,2)
df2<-df2[match(PopulationNew, df2$V1),]

df2$V1 <- factor(df2$V1, levels = unique(PopulationNew))

pdf("Core_Pan_genes_ALL.pdf",height = 8, width = 20)
ggplot(data=NULL, aes(x=df2$V1, y=df2$V2, fill=df2$V3,)) + labs(fill = "Genes", y="Number of genes")+ scale_y_continuous(limits = c(0, 4300)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),axis.title.x=element_blank(),
                                    axis.text.y = element_text(size=12),axis.title.y= element_text(size=14)
  ) 
dev.off()


###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################





#########################################################################
#########################################################################
## Correlation NON conserved and Phylogenetic distance
#########################################################################



HOMO<- list()
for ( i in 1:length(NonConserved_ALL6POP)){
  Homolog_info[[1]]
    tmp_homo<-Homolog_info[grep(paste(unlist(lapply(strsplit(as.character(NonConserved_ALL6POP[[i]][1]),split="_"), function (x) x[2])),"$",sep=""), names(Homolog_info), perl=T)]
  tmp_Gene<-as.data.frame(NonConserved_ALL6POP[[i]])
  names(tmp_Gene)<-c("Gene")
  HOMO[[i]]<-merge(tmp_Gene,tmp_homo[[1]],by="Gene", all.x=F)
}
names(NonConserved_ALL6POP)
######### METHOD TO MERGE WITH LESS MEMORY
# CHANGE TO DATA TABLE
library(data.table)

HOMO<-lapply(HOMO, function (x) x[!duplicated(x),] )

dt<-lapply(HOMO, function (x) data.table(x) ) # transform to data.table
dt<-lapply(dt,  function (x) setNames(x,nm=rev(colnames(x)))  )
dt<-lapply(dt,  function (x) as.data.frame(x) ) # go back to data.frame

N16961HOMO<-cbind.data.frame(Gene=unique(unlist(lapply(dt, function (x) x[,2]))))

N16961HOMO_dt<-data.table(N16961HOMO)  # go back to data.table
dt<-lapply(dt, function (x) data.table(x) )

for (i in 1:length(dt)) {
  N16961HOMO_dt<-N16961HOMO_dt[dt[[i]],on="Gene"]
}

N16961HOMO<-as.data.frame(N16961HOMO_dt) # FINALLY RETURN TO DATA.FRAME

for (i in 1:dim(N16961HOMO)[2]) {
  N16961HOMO[,i]<-as.character(N16961HOMO[,i])
}


colnames(N16961HOMO)<-c("REF",unlist(lapply(HOMO, function (x) colnames(x)[2])))
N16961HOMO[is.na(N16961HOMO)]<-0
N16961HOMO[N16961HOMO!=0]<-1
N16961HOMO<-N16961HOMO[,-1]

for (i in 1:dim(N16961HOMO)[2]){
  N16961HOMO[,i]<-as.numeric(N16961HOMO[,i])
}

colnames(N16961HOMO)

VariableGenes_Dist<-as.data.frame(t(N16961HOMO))


rownames(VariableGenes_Dist)

write.table()

plot(PatristicDistMatrix)


plot(VariableGenes_Dist[1:2])
VariableGenes_Dist<-vegan::vegdist(VariableGenes_Dist,binary = T,method = "bray")
VariableGenes_Dist<-dist(N16961HOMO,) 
colnames(VariableGenes_Dist)<-paste("A",1:dim(N16961HOMO)[2],sep="_")
# Phylogenetic distance
library(ape)

Populations<-read.csv("~/Documents/EPFL/Vibrio_Chromosomal_rearrangements/PopGenomics_V2/ALL_POPS/Phylogenies_ALL/Strains_Pop_vf.csv", h=T)
Populations$strain<-gsub("-","_",Populations$strain)
tree <- read.tree('~/Documents/EPFL/Vibrio_Chromosomal_rearrangements/PopGenomics_V2/ALL_POPS/Phylogenies_ALL/aLRT_WAG/ALL6POP_OUT_Vibrio_CoreGenes_MafftAlign.phy_phyml_tree.txt')
tree$tip.label<-gsub("_faa","", gsub("Core_","", tree$tip.label)   )
tree  <- root(tree, outgroup = "Vf_ATCC33809",resolve.root = T)
tree <- drop.tip(tree, "Vf_ATCC33809")

tree$node.label
PatristicDistMatrix<-cophenetic(tree)



############################################################################################################################################################################################################################
###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################






#########################################################################
######################################################################
################### ACP ##############################################
# then MErge and create an ACP to see wich samples are similar per POP


HOMO<- list()
for ( i in 1:length(NonConserved_Cali)){
Homolog_info[[1]]
tmp_homo<-Homolog_info[grep(unlist(lapply(strsplit(as.character(NonConserved_Cali[[i]][1]),split="_"), function (x) x[2])), names(Homolog_info))]
tmp_Gene<-as.data.frame(NonConserved_Cali[[i]])
names(tmp_Gene)<-c("Gene")
HOMO[[i]]<-merge(tmp_Gene,tmp_homo[[1]],by="Gene", all.x=F)
}



N16961HOMO<-cbind.data.frame(GENE=unique(unlist(lapply(HOMO, function (x) x[,2]))))

for (i in 1:length(HOMO)) {
N16961HOMO<-merge(N16961HOMO,HOMO[[i]],by.x="GENE",by.y=colnames(HOMO[[i]])[2], all.x=T, all.y=T)
N16961HOMO[,i]<-as.character(N16961HOMO[,i])
}
N16961HOMO[,i+1]<-as.character(N16961HOMO[,i+1])
summary(N16961HOMO)
colnames(N16961HOMO)<-c("REF",unlist(lapply(HOMO, function (x) colnames(x)[2])))
#rownames(N16961HOMO)<-N16961HOMO[,1]

N16961HOMO[is.na(N16961HOMO)]<-0
N16961HOMO[N16961HOMO!=0]<-1
N16961HOMO<-N16961HOMO[,-1]

for (i in 1:dim(N16961HOMO)[2]){
  N16961HOMO[,i]<-as.numeric(N16961HOMO[,i])
}



######## PCA
PC<-prcomp(t(N16961HOMO))
PC$sdev/ sum(PC$sdev)
PCi<-data.frame(PC$x,Strain=colnames(N16961HOMO))
ggplot(PCi,aes(x=PC1,y=PC2,col=Strain, label=colnames(N16961HOMO))) + geom_text(aes(fontface=2)) +
  theme_classic()+ theme(legend.position="none") + labs(x = paste("PC1",round((PC$sdev/ sum(PC$sdev))[1],digits=4)*100," % of variation"),y= paste("PC2",round((PC$sdev/ sum(PC$sdev))[2],digits=4)*100," % of variation")) + scale_x_continuous(labels = comma)
# Find Biological process of pandemic in total and per sample. 


#################################################################################################################
# PCA ALL SAMPLES ALL POPULATIONS NON CONSERVED GENES
ALL_POPS<-c(NonConserved_Cali,NonConserved_Haiti,NonConserved_Masa,NonConserved_Mex,NonConserved_Thai,NonConserved_Pandemic)

HOMO<- list()
for ( i in 1:length(ALL_POPS)){
  Homolog_info[[70]]
  tmp_homo<-Homolog_info[grep(unlist(lapply(strsplit(as.character(ALL_POPS[[i]][1]),split="_"), function (x) x[2])), names(Homolog_info))]
  tmp_Gene<-as.data.frame(ALL_POPS[[i]])
  names(tmp_Gene)<-c("Gene")
  HOMO[[i]]<-merge(tmp_Gene,tmp_homo[[1]],by="Gene", all.x=F)
}



N16961HOMO<-cbind.data.frame(GENE=unique(unlist(lapply(HOMO, function (x) x[,2]))))

for (i in 1:length(HOMO)) {
  N16961HOMO<-merge(N16961HOMO,HOMO[[i]],by.x="GENE",by.y=colnames(HOMO[[i]])[2], all.x=T, all.y=T)
  N16961HOMO[,i]<-as.character(N16961HOMO[,i])
}
N16961HOMO[,i+1]<-as.character(N16961HOMO[,i+1])
head(N16961HOMO)
colnames(N16961HOMO)<-c("REF",unlist(lapply(HOMO, function (x) colnames(x)[2])))
#rownames(N16961HOMO)<-N16961HOMO[,1]


N16961HOMO[is.na(N16961HOMO)]<-0
N16961HOMO[N16961HOMO!=0]<-1
N16961HOMO<-N16961HOMO[,-1]

for (i in 1:dim(N16961HOMO)[2]){
  N16961HOMO[,i]<-as.numeric(N16961HOMO[,i])
}

poblaciones<-c(rep("California",length(NonConserved_Cali)),rep("Haiti",length(NonConserved_Haiti)),rep("Massachusetts",length(NonConserved_Masa)),rep("Mexico",length(NonConserved_Mex)),
               rep("Thailand",length(NonConserved_Thai)),rep("Pandemic",length(NonConserved_Pandemic)))

######## PCA
PC<-prcomp(t(N16961HOMO))
PC$sdev/ sum(PC$sdev)
PCi<-data.frame(PC$x,Strain=colnames(N16961HOMO),Pop=poblaciones)
ggplot(PCi,aes(x=PC1,y=PC2,col=poblaciones, label=colnames(N16961HOMO))) + geom_text(aes(fontface=2)) +
  theme_classic()+ theme(legend.position="none") + labs(x = paste("PC1",round((PC$sdev/ sum(PC$sdev))[1],digits=4)*100," % of variation"),y= paste("PC2",round((PC$sdev/ sum(PC$sdev))[2],digits=4)*100," % of variation")) + scale_x_continuous(labels = comma)


N16961HOMO





###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################


# PAN GENES PROBLEM OF ID SIMILAR TO N16961 BIAS.
# Sets PLOT NON CONSERVED GENES

ALL_POPS<-list(California=unlist(NonConserved_Cali),Haiti=unlist(NonConserved_Haiti),Massachusetts=unlist(NonConserved_Masa),Mexico=unlist(NonConserved_Mex),Thailand=unlist(NonConserved_Thai),Pandemic=unlist(NonConserved_Pandemic))

HOMO<- list()
for ( i in 1:length(ALL_POPS)){
  tmp_homo<-Homolog_info[grep(unlist(lapply(strsplit(as.character(ALL_POPS[[i]][1]),split="_"), function (x) x[2])), names(Homolog_info))]
  tmp_Gene<-as.data.frame(ALL_POPS[[i]])
  names(tmp_Gene)<-c("Gene")
  HOMO[[i]]<-merge(tmp_Gene,tmp_homo[[1]],by="Gene", all.x=F)
}



N16961HOMO<-cbind.data.frame(GENE=unique(unlist(lapply(HOMO, function (x) x[,2]))))

for (i in 1:length(HOMO)) {
  N16961HOMO<-merge(N16961HOMO,HOMO[[i]],by.x="GENE",by.y=colnames(HOMO[[i]])[2], all.x=T, all.y=T)
  N16961HOMO[,i]<-as.character(N16961HOMO[,i])
}
N16961HOMO[,i+1]<-as.character(N16961HOMO[,i+1])
head(N16961HOMO)
colnames(N16961HOMO)<-c("BLA",names(ALL_POPS))
#rownames(N16961HOMO)<-N16961HOMO[,1]


N16961HOMO[is.na(N16961HOMO)]<-0
N16961HOMO[N16961HOMO!=0]<-1
N16961HOMO<-N16961HOMO[,-1]

for (i in 1:dim(N16961HOMO)[2]){
  N16961HOMO[,i]<-as.numeric(N16961HOMO[,i])
}
N16961HOMO
pdf("Sets_plot_unique_variable_genes_POP.pdf",height = 8,width = 20)
upset(N16961HOMO,  number.angles = 0, point.size = 3.5, line.size = 2, keep.order = T, sets = c("Pandemic","Mexico","California","Massachusetts","Haiti","Thailand"),
      mainbar.y.label = "Gene intersections \n between populations", sets.x.label = "Unique non-conserved \ngenes in Pop", order.by = "freq",
      text.scale = c(2, 4,2, 2, 2, 3))
dev.off()

# extract names 


#############################################################

upset(N16961HOMO,  number.angles = 0, point.size = 3.5, line.size = 2, keep.order = T, sets = c("Pandemic","Mexico","California","Massachusetts","Haiti","Thailand"),
      mainbar.y.label = "Gene intersections \n between populations", sets.x.label = "Unique non-conserved \ngenes in Pop", order.by = "freq",
      text.scale = c(2, 4,2, 2, 2, 3))

NAMES
#################################################################################################################
# PCA ALL SAMPLES ALL POPULATIONS ALL GENES

ALL_Cali<-list()
for (i in 1:length(Core_Cali)){
  ALL_Cali[[i]]<-c(as.character(Core_Cali[[i]]),as.character(NonConserved_Cali[[i]]))
  names(ALL_Cali)[i]<-names(Core_Cali)
}
ALL_Haiti<-list()
for (i in 1:length(Core_Haiti)){
  ALL_Haiti[[i]]<-c(as.character(Core_Haiti[[i]]),as.character(NonConserved_Haiti[[i]]))
  names(ALL_Haiti)[i]<-names(Core_Haiti)
}
ALL_Masa<-list()
for (i in 1:length(Core_Masa)){
  ALL_Masa[[i]]<-c(as.character(Core_Masa[[i]]),as.character(NonConserved_Masa[[i]]))
  names(ALL_Masa)[i]<-names(Core_Masa)
}
ALL_Mex<-list()
for (i in 1:length(Core_Mex)){
  ALL_Mex[[i]]<-c(as.character(Core_Mex[[i]]),as.character(NonConserved_Mex[[i]]))
  names(ALL_Mex)[i]<-names(Core_Mex)
}
ALL_Thai<-list()
for (i in 1:length(Core_Thai)){
  ALL_Thai[[i]]<-c(as.character(Core_Thai[[i]]),as.character(NonConserved_Thai[[i]]))
  names(ALL_Thai)[i]<-names(Core_Thai)
}
ALL_Pan<-list()
for (i in 1:length(Core_Pan)){
  ALL_Pan[[i]]<-c(as.character(Core_Pan[[i]]),as.character(NonConserved_Pandemic[[i]]))
  names(ALL_Pan)[i]<-names(Core_Pan)
}

ALL_POPS<-c(ALL_Cali,ALL_Haiti,ALL_Masa,ALL_Mex,ALL_Thai,ALL_Pan)
ALL_POPS
HOMO<- list()
for ( i in 1:length(ALL_POPS)){
  Homolog_info[[70]]
  tmp_homo<-Homolog_info[grep(unlist(lapply(strsplit(as.character(ALL_POPS[[i]][1]),split="_"), function (x) x[2])), names(Homolog_info))]
  tmp_Gene<-as.data.frame(ALL_POPS[[i]])
  names(tmp_Gene)<-c("Gene")
  HOMO[[i]]<-merge(tmp_Gene,tmp_homo[[1]],by="Gene", all.x=F)
}

#################################################################################################################

## CORE GENES COMPARISON IN POP

Core_Cali[[1]]
Homolog_info[[70]]
names(Homolog_info)

ALL_POPS_core<-list(California=Core_Cali[[1]],Haiti=Core_Haiti[[1]],Massachusetts=Core_Masa[[1]],Mexico=Core_Mex[[1]],Thailand=Core_Thai[[1]],Pandemic=Core_Pan[[1]])
ALL_POPS_core
HOMO_core<- list()
for ( i in 1:length(ALL_POPS_core)){
  Homolog_info[[70]]
  tmp_homo<-Homolog_info[grep(unlist(lapply(strsplit(as.character(ALL_POPS_core[[i]][1]),split="_"), function (x) x[2])), names(Homolog_info))]
  tmp_Gene<-as.data.frame(ALL_POPS_core[[i]])
  names(tmp_Gene)<-c("Gene")
  HOMO_core[[i]]<-merge(tmp_Gene,tmp_homo[[1]],by="Gene", all.x=F)
}

lapply(HOMO_core,function (x) dim(x))
HOMO_core[[1]]

#########


N16961HOMO_core<-cbind.data.frame(GENE=unique(unlist(lapply(HOMO_core, function (x) x[,2]))))

for (i in 1:length(HOMO_core)) {
  N16961HOMO_core<-merge(N16961HOMO_core,HOMO_core[[i]],by.x="GENE",by.y=colnames(HOMO_core[[i]])[2], all.x=T, all.y=T)
  N16961HOMO_core[,i]<-as.character(N16961HOMO_core[,i])
}
N16961HOMO_core[,i+1]<-as.character(N16961HOMO_core[,i+1])
head(N16961HOMO_core)
colnames(N16961HOMO_core)<-c("BLA",names(ALL_POPS_core))
#rownames(N16961HOMO)<-N16961HOMO[,1]


N16961HOMO_core[is.na(N16961HOMO_core)]<-0
N16961HOMO_core[N16961HOMO_core!=0]<-1
N16961HOMO_core<-N16961HOMO_core[,-1]

for (i in 1:dim(N16961HOMO_core)[2]){
  N16961HOMO_core[,i]<-as.numeric(N16961HOMO_core[,i])
}
N16961HOMO_core
dim(N16961HOMO_core)
pdf("Sets_plot_unique_core_genes_POP.pdf",height = 8,width = 20)
upset(N16961HOMO_core,  number.angles = 0, point.size = 3.5, line.size = 2, keep.order = T, sets = c("Pandemic","Mexico","California","Massachusetts","Haiti","Thailand"),
      mainbar.y.label = "Gene intersections \n between populations", sets.x.label = "Conserved \ngenes in Pop", order.by = "freq",
      text.scale = c(2, 4,2, 2, 2, 3))
dev.off()

upset(N16961HOMO_core)

#################################################################################################################

## aLL GENES IN POP


ALL_POPS_EVERYTHING<-list(California=unlist(ALL_Cali),Haiti=unlist(ALL_Haiti),Massachusetts=unlist(ALL_Masa),Mexico=unlist(ALL_Mex),Thailand=unlist(ALL_Thai),Pandemic=unlist(ALL_Pan))
ALL_POPS_EVERYTHING


HOMO_EVERYTHING<- list()
for ( i in 1:length(ALL_POPS_EVERYTHING)){
  Homolog_info[[70]]
  tmp_homo<-Homolog_info[grep(unlist(lapply(strsplit(as.character(ALL_POPS_EVERYTHING[[i]][1]),split="_"), function (x) x[2])), names(Homolog_info))]
  tmp_Gene<-as.data.frame(ALL_POPS_EVERYTHING[[i]])
  names(tmp_Gene)<-c("Gene")
  HOMO_EVERYTHING[[i]]<-merge(tmp_Gene,tmp_homo[[1]],by="Gene", all.x=F)
}

lapply(HOMO_EVERYTHING,function (x) dim(x))
HOMO_EVERYTHING[[1]]

#########


N16961HOMO_EVERYTHING<-cbind.data.frame(GENE=unique(unlist(lapply(HOMO_EVERYTHING, function (x) x[,2]))))

for (i in 1:length(HOMO_EVERYTHING)) {
  N16961HOMO_EVERYTHING<-merge(N16961HOMO_EVERYTHING,HOMO_EVERYTHING[[i]],by.x="GENE",by.y=colnames(HOMO_EVERYTHING[[i]])[2], all.x=T, all.y=T)
  N16961HOMO_EVERYTHING[,i]<-as.character(N16961HOMO_EVERYTHING[,i])
}
N16961HOMO_EVERYTHING[,i+1]<-as.character(N16961HOMO_EVERYTHING[,i+1])
head(N16961HOMO_EVERYTHING)
colnames(N16961HOMO_EVERYTHING)<-c("BLA",names(ALL_POPS_EVERYTHING))
#rownames(N16961HOMO)<-N16961HOMO[,1]


N16961HOMO_EVERYTHING[is.na(N16961HOMO_EVERYTHING)]<-0
N16961HOMO_EVERYTHING[N16961HOMO_EVERYTHING!=0]<-1
N16961HOMO_EVERYTHING<-N16961HOMO_EVERYTHING[,-1]

for (i in 1:dim(N16961HOMO_EVERYTHING)[2]){
  N16961HOMO_EVERYTHING[,i]<-as.numeric(N16961HOMO_EVERYTHING[,i])
}
N16961HOMO_EVERYTHING
dim(N16961HOMO_EVERYTHING)
pdf("Sets_plot_unique_everything_genes_POP.pdf",height = 8,width = 20)
upset(N16961HOMO_EVERYTHING,  number.angles = 0, point.size = 3.5, line.size = 2, keep.order = T, sets = c("Pandemic","Mexico","California","Massachusetts","Haiti","Thailand"),
      mainbar.y.label = "Gene intersections \n between populations", sets.x.label = "Conserved \ngenes in Pop", order.by = "freq",
      text.scale = c(2, 4,2, 2, 2, 3))
dev.off()

upset(N16961HOMO_EVERYTHING)



###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################
###########################################################################################################################################################################################################################





#################################################################################################################
#################################################################################################################
###################### # Notable genes in Pops



# in bash
#a=0;for i in $(ls *.fna); do echo $(echo $i | cut -d'_' -f3) ;makeblastdb -dbtype nucl -in $i -parse_seqids -out db_genomes/$(echo $i | cut -d'_' -f3)_db ; done
#a=0;for i in $(ls *.fna); do echo $(echo $i | cut -d'_' -f3) ;tblastn -db db_genomes/$(echo $i | cut -d'_' -f3)_db -outfmt 6 -evalue 1e-6 -show_gis -num_alignments 1 -max_hsps 20 -num_threads 30 -out db_genomes/blast_VSP-2_$(echo $i | cut -d'_' -f3).xml -query ~/Documents/EPFL/Vibrio_Chromosomal_rearrangements/Notable_Regions/VSP-2.faa ; done
#a=0;for i in $(ls *.fna); do echo $(echo $i | cut -d'_' -f3) ;tblastn -db db_genomes/$(echo $i | cut -d'_' -f3)_db -outfmt 6 -evalue 1e-6 -show_gis -num_alignments 1 -max_hsps 20 -num_threads 30 -out db_genomes/blast_VSP-1_$(echo $i | cut -d'_' -f3).xml -query ~/Documents/EPFL/Vibrio_Chromosomal_rearrangements/Notable_Regions/VSP-1.faa ; done
#a=0;for i in $(ls *.fna); do echo $(echo $i | cut -d'_' -f3) ;tblastn -db db_genomes/$(echo $i | cut -d'_' -f3)_db -outfmt 6 -evalue 1e-6 -show_gis -num_alignments 1 -max_hsps 20 -num_threads 30 -out db_genomes/blast_VPI-2_$(echo $i | cut -d'_' -f3).xml -query ~/Documents/EPFL/Vibrio_Chromosomal_rearrangements/Notable_Regions/VPI-2.faa ; done
#a=0;for i in $(ls *.fna); do echo $(echo $i | cut -d'_' -f3) ;tblastn -db db_genomes/$(echo $i | cut -d'_' -f3)_db -outfmt 6 -evalue 1e-6 -show_gis -num_alignments 1 -max_hsps 20 -num_threads 30 -out db_genomes/blast_VPI-1_$(echo $i | cut -d'_' -f3).xml -query ~/Documents/EPFL/Vibrio_Chromosomal_rearrangements/Notable_Regions/VPI-1.faa ; done
#a=0;for i in $(ls *.fna); do echo $(echo $i | cut -d'_' -f3) ;tblastn -db db_genomes/$(echo $i | cut -d'_' -f3)_db -outfmt 6 -evalue 1e-6 -show_gis -num_alignments 1 -max_hsps 20 -num_threads 30 -out db_genomes/blast_Serogroup_$(echo $i | cut -d'_' -f3).xml -query ~/Documents/EPFL/Vibrio_Chromosomal_rearrangements/Notable_Regions/Serogroup_specific_genes.faa; done
#a=0;for i in $(ls *.fna); do echo $(echo $i | cut -d'_'  -f3) ;tblastn -db db_genomes/$(echo $i | cut -d'_' -f3)_db -outfmt 6 -evalue 1e-6 -show_gis -num_alignments 1 -max_hsps 20 -num_threads 30 -out db_genomes/blast_Virulence_$(echo $i | cut -d'_' -f3).xml -query ~/Documents/EPFL/Vibrio_Chromosomal_rearrangements/Notable_Regions/Virulence_genes.faa ; done



library(pheatmap)
library(plyr)
library(ggplot2)
setwd('~/Documents/EPFL/Vibrio_Chromosomal_rearrangements/PopGenomics_V2/ALL_POPS/db_genomes/')

######################################################################################################################################
#######################################################################################################################################
############## ANALYSIS OF NOTABLE REGIONS

# Serogroup

filesToProcess <- dir(pattern = "*\\.xml$")  #files to pro# if event 3 merged
filesToProcess<-filesToProcess[grep('Serogroup',filesToProcess)]

table(gsub('.xml','',unlist(lapply(strsplit(filesToProcess,"_"), function (x) x[3]))) )

listOfFiles <- lapply(filesToProcess, function(x) tryCatch(read.table(x, header = F),
                                                           error= function (e) cbind.data.frame(V1="NA",V2="NA",V3="NA",
                                                                                                V4="NA",V5="NA",V6="NA",
                                                                                                V7=0,V8=0,V9="NA",
                                                                                                V10="NA",V11="NA",V12="NA")))
names(listOfFiles)<-gsub("blastProt_","",gsub(".xml","",gsub("Blast_","",filesToProcess))   )
colnam<-c("Gene","Chromosome","Ident_perc","Alig_length","mismatch","gapopens",
          "Start","End","Target_start","Target_end","evalue","bitscore")
listOfFiles <- lapply(listOfFiles, setNames, nm=colnam)
listOfFiles[[2]]
names(listOfFiles)

# filter 80% identity
listOfFiles<-lapply(listOfFiles, function(x) x[x$Ident_perc>80,] )


listOfFiles<-lapply(listOfFiles, function(x) table(x[,1]) )
names(listOfFiles)<-unlist(lapply(strsplit(names(listOfFiles),split = '_'),function(x) x[3]))


df<-merge(as.data.frame(listOfFiles[[1]]),as.data.frame(listOfFiles[[2]]),by='Var1',all = T)
for (i in 3:length(listOfFiles)){
  df<-merge(df,listOfFiles[[i]],by='Var1',all = T)
}
df<-df[!df$Var1=="NA",]

colnames(df)[2:dim(df)[2]]<- names(listOfFiles)
rownames(df)<-df[,1]

df<-df[,-1]
df[is.na(df)] <- 0
df[df==2] <- 1

df_serogroup<-df
pheatmap(df_serogroup,cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,legend = F,color = c('white','black'))

######################################################################################################################################
# Virulence

filesToProcess <- dir(pattern = "*\\.xml$")  #files to pro# if event 3 merged
filesToProcess<-filesToProcess[grep('Virulence',filesToProcess)]

table(gsub('.xml','',unlist(lapply(strsplit(filesToProcess,"_"), function (x) x[3]))) )

listOfFiles <- lapply(filesToProcess, function(x) tryCatch(read.table(x, header = F),
                                                           error= function (e) cbind.data.frame(V1="NA",V2="NA",V3="NA",
                                                                                                V4="NA",V5="NA",V6="NA",
                                                                                                V7=0,V8=0,V9="NA",
                                                                                                V10="NA",V11="NA",V12="NA")))
names(listOfFiles)<-gsub("blastProt_","",gsub(".xml","",gsub("Blast_","",filesToProcess))   )
colnam<-c("Gene","Chromosome","Ident_perc","Alig_length","mismatch","gapopens",
          "Start","End","Target_start","Target_end","evalue","bitscore")
listOfFiles <- lapply(listOfFiles, setNames, nm=colnam)
listOfFiles[[2]]
names(listOfFiles)

# filter 80% identity
listOfFiles<-lapply(listOfFiles, function(x) x[x$Ident_perc>80,] )


listOfFiles<-lapply(listOfFiles, function(x) table(x[,1]) )
names(listOfFiles)<-unlist(lapply(strsplit(names(listOfFiles),split = '_'),function(x) x[3]))


df<-merge(as.data.frame(listOfFiles[[1]]),as.data.frame(listOfFiles[[2]]),by='Var1',all = T)
for (i in 3:length(listOfFiles)){
  df<-merge(df,listOfFiles[[i]],by='Var1',all = T)
}
df<-df[!df$Var1=="NA",]

colnames(df)[2:dim(df)[2]]<- names(listOfFiles)
rownames(df)<-df[,1]

df<-df[,-1]
df[is.na(df)] <- 0
df[df==2] <- 1
df_virulence<-df
pheatmap(df_virulence,cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,legend = T)


######################################################################################################################################
# VIP1

filesToProcess <- dir(pattern = "*\\.xml$")  #files to pro# if event 3 merged
filesToProcess<-filesToProcess[grep('VPI-1',filesToProcess)]

table(gsub('.xml','',unlist(lapply(strsplit(filesToProcess,"_"), function (x) x[3]))) )

listOfFiles <- lapply(filesToProcess, function(x) tryCatch(read.table(x, header = F),
                                                           error= function (e) cbind.data.frame(V1="NA",V2="NA",V3="NA",
                                                                                                V4="NA",V5="NA",V6="NA",
                                                                                                V7=0,V8=0,V9="NA",
                                                                                                V10="NA",V11="NA",V12="NA")))
names(listOfFiles)<-gsub("blastProt_","",gsub(".xml","",gsub("Blast_","",filesToProcess))   )
colnam<-c("Gene","Chromosome","Ident_perc","Alig_length","mismatch","gapopens",
          "Start","End","Target_start","Target_end","evalue","bitscore")
listOfFiles <- lapply(listOfFiles, setNames, nm=colnam)
listOfFiles[[2]]
names(listOfFiles)

# filter 80% identity
listOfFiles<-lapply(listOfFiles, function(x) x[x$Ident_perc>80,] )


listOfFiles<-lapply(listOfFiles, function(x) table(x[,1]) )
names(listOfFiles)<-unlist(lapply(strsplit(names(listOfFiles),split = '_'),function(x) x[3]))


df<-merge(as.data.frame(listOfFiles[[1]]),as.data.frame(listOfFiles[[2]]),by='Var1',all = T)
for (i in 3:length(listOfFiles)){
  df<-merge(df,listOfFiles[[i]],by='Var1',all = T)
}
df<-df[!df$Var1=="NA",]

colnames(df)[2:dim(df)[2]]<- names(listOfFiles)
rownames(df)<-df[,1]

df<-df[,-1]
df[is.na(df)] <- 0
df[df==2] <- 1
df_VIP1<-df
pheatmap(df_VIP1,cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,legend = T)


######################################################################################################################################
# VIP2

filesToProcess <- dir(pattern = "*\\.xml$")  #files to pro# if event 3 merged
filesToProcess<-filesToProcess[grep('VPI-2',filesToProcess)]

table(gsub('.xml','',unlist(lapply(strsplit(filesToProcess,"_"), function (x) x[3]))) )

listOfFiles <- lapply(filesToProcess, function(x) tryCatch(read.table(x, header = F),
                                                           error= function (e) cbind.data.frame(V1="NA",V2="NA",V3="NA",
                                                                                                V4="NA",V5="NA",V6="NA",
                                                                                                V7=0,V8=0,V9="NA",
                                                                                                V10="NA",V11="NA",V12="NA")))
names(listOfFiles)<-gsub("blastProt_","",gsub(".xml","",gsub("Blast_","",filesToProcess))   )
colnam<-c("Gene","Chromosome","Ident_perc","Alig_length","mismatch","gapopens",
          "Start","End","Target_start","Target_end","evalue","bitscore")
listOfFiles <- lapply(listOfFiles, setNames, nm=colnam)
listOfFiles[[2]]
names(listOfFiles)

# filter 80% identity
listOfFiles<-lapply(listOfFiles, function(x) x[x$Ident_perc>80,] )


listOfFiles<-lapply(listOfFiles, function(x) table(x[,1]) )
names(listOfFiles)<-unlist(lapply(strsplit(names(listOfFiles),split = '_'),function(x) x[3]))


df<-merge(as.data.frame(listOfFiles[[1]]),as.data.frame(listOfFiles[[2]]),by='Var1',all = T)
for (i in 3:length(listOfFiles)){
  df<-merge(df,listOfFiles[[i]],by='Var1',all = T)
}
df<-df[!df$Var1=="NA",]

colnames(df)[2:dim(df)[2]]<- names(listOfFiles)
rownames(df)<-df[,1]

df<-df[,-1]
df[is.na(df)] <- 0
df[df==2] <- 1
df_VIP2<-df
pheatmap(df_VIP2,cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,legend = T)


######################################################################################################################################
# VSP1

filesToProcess <- dir(pattern = "*\\.xml$")  #files to pro# if event 3 merged
filesToProcess<-filesToProcess[grep('VSP-1',filesToProcess)]

table(gsub('.xml','',unlist(lapply(strsplit(filesToProcess,"_"), function (x) x[3]))) )

listOfFiles <- lapply(filesToProcess, function(x) tryCatch(read.table(x, header = F),
                                                           error= function (e) cbind.data.frame(V1="NA",V2="NA",V3="NA",
                                                                                                V4="NA",V5="NA",V6="NA",
                                                                                                V7=0,V8=0,V9="NA",
                                                                                                V10="NA",V11="NA",V12="NA")))
names(listOfFiles)<-gsub("blastProt_","",gsub(".xml","",gsub("Blast_","",filesToProcess))   )
colnam<-c("Gene","Chromosome","Ident_perc","Alig_length","mismatch","gapopens",
          "Start","End","Target_start","Target_end","evalue","bitscore")
listOfFiles <- lapply(listOfFiles, setNames, nm=colnam)
listOfFiles[[1]]
names(listOfFiles)

# filter 80% identity
listOfFiles<-lapply(listOfFiles, function(x) x[x$Ident_perc>80,] )


listOfFiles<-lapply(listOfFiles, function(x) table(x[,1]) )
names(listOfFiles)<-unlist(lapply(strsplit(names(listOfFiles),split = '_'),function(x) x[3]))


df<-merge(as.data.frame(listOfFiles[[1]]),as.data.frame(listOfFiles[[2]]),by='Var1',all = T)
for (i in 3:length(listOfFiles)){
  df<-merge(df,listOfFiles[[i]],by='Var1',all = T)
}
df<-df[!df$Var1=="NA",]
colnames(df)[2:dim(df)[2]]<- names(listOfFiles)
rownames(df)<-df[,1]

df<-df[,-1]
df[is.na(df)] <- 0
df[df==2] <- 1
df_VSP1<-df[-6,]
pheatmap(df_VSP1,cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,legend = T)


######################################################################################################################################
# VSP2

filesToProcess <- dir(pattern = "*\\.xml$")  #files to pro# if event 3 merged
filesToProcess<-filesToProcess[grep('VSP-2',filesToProcess)]

table(gsub('.xml','',unlist(lapply(strsplit(filesToProcess,"_"), function (x) x[3]))) )

listOfFiles <- lapply(filesToProcess, function(x) tryCatch(read.table(x, header = F),
                                                           error= function (e) cbind.data.frame(V1="NA",V2="NA",V3="NA",
                                                                                                V4="NA",V5="NA",V6="NA",
                                                                                                V7=0,V8=0,V9="NA",
                                                                                                V10="NA",V11="NA",V12="NA")))
names(listOfFiles)<-gsub("blastProt_","",gsub(".xml","",gsub("Blast_","",filesToProcess))   )
colnam<-c("Gene","Chromosome","Ident_perc","Alig_length","mismatch","gapopens",
          "Start","End","Target_start","Target_end","evalue","bitscore")
listOfFiles <- lapply(listOfFiles, setNames, nm=colnam)
listOfFiles[[2]]
names(listOfFiles)

# filter 80% identity
listOfFiles<-lapply(listOfFiles, function(x) x[x$Ident_perc>80,] )


listOfFiles<-lapply(listOfFiles, function(x) table(x[,1]) )
names(listOfFiles)<-unlist(lapply(strsplit(names(listOfFiles),split = '_'),function(x) x[3]))


df<-merge(as.data.frame(listOfFiles[[1]]),as.data.frame(listOfFiles[[2]]),by='Var1',all = T)
for (i in 3:length(listOfFiles)){
  df<-merge(df,listOfFiles[[i]],by='Var1',all = T)
}
df<-df[!df$Var1=="NA",]

colnames(df)[2:dim(df)[2]]<- names(listOfFiles)
rownames(df)<-df[,1]

df<-df[,-1]
df[is.na(df)] <- 0
df[df==2] <- 1
df_VSP2<-df
pheatmap(df_VSP2,cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,legend = T)


###################################################################################################################################3
# ALL TOGETHER

df_ALL<-rbind.data.frame(df_serogroup,df_virulence,df_VIP1,df_VIP2,df_VSP1,df_VSP2)

notable_genes<-data.frame(c(rep("Serogroup",dim(df_serogroup)[1]), rep("Virulence",dim(df_virulence)[1]) ,rep("VIP-1",dim(df_VIP1)[1]),rep("VIP-2",dim(df_VIP2)[1]),
                            rep("VSP-1",dim(df_VSP1)[1]),rep("VSP-2",dim(df_VSP2)[1])))
colnames(notable_genes)<-'Notable regions'
rownames(notable_genes)<-rownames(df_ALL)
colnames(df_ALL)

OrderPOP<-c(unlist(lapply(strsplit(names(Core_Cali),split="_"), function (x) x[2])),
unlist(lapply(strsplit(names(Core_Haiti),split="_"), function (x) x[2])),
unlist(lapply(strsplit(names(Core_Masa),split="_"), function (x) x[2])),
unlist(lapply(strsplit(names(Core_Mex),split="_"), function (x) x[2])),
unlist(lapply(strsplit(names(Core_Thai),split="_"), function (x) x[2])),
unlist(lapply(strsplit(names(Core_Pan),split="_"), function (x) x[2]))  )

population<-data.frame(c(rep("California",length(unlist(lapply(strsplit(names(Core_Cali),split="_"), function (x) x[2])))),
  rep("Haiti",length(unlist(lapply(strsplit(names(Core_Haiti),split="_"), function (x) x[2])))),
                          rep("Massachusetts",length(unlist(lapply(strsplit(names(Core_Masa),split="_"), function (x) x[2])))),
                                                  rep("Mexico",length(unlist(lapply(strsplit(names(Core_Mex),split="_"), function (x) x[2])))),
 rep("Thailand",length(unlist(lapply(strsplit(names(Core_Thai),split="_"), function (x) x[2])))),
  rep("Pandemic",length(unlist(lapply(strsplit(names(Core_Pan),split="_"), function (x) x[2]))))  ) )

colnames(population)<-'Populations'
rownames(population)<-OrderPOP

df_ALL<-df_ALL[OrderPOP]

## ORDER PHYLO
OrderPhylo<-c("372","504","1401","5032","601","985","1992","2007","2006-2004","211","688","971","210","204","354","IEC224","A1552","2012Env-131","2010EL-1786","2012EL-2176","2012Env-326","H1","2012Env-94","2012Env-90","MJ-1236","N16961","CRC711","E1162","E9120","C5","M66-2","M2140","NCTC9420","NCTC5395","E1320","MS6","E506","O395","2012EL-1759","Env-390","2012Env-9","YB3B05","YB4C07","YB8E08","W10G","VSH-2-22I","VHS-1-22I","IWSH488-TCY","2012Env-25","MK14","SO5Y","VID97-3-2","I-SSH2-TY6","TC183-1-1","I-WASTE-SSH2-TY1","D-SSH2-1-TY8","I-WHSH1-TY4","2012Env-2","TC8-1-1","2012Env-92","1130-W58","4T5","D-SSH2-1-TY2","TC151-1-2","1262-W278","SL4G","Sa5Y","SL5Y","YB2G07","YB4G05","YB4B03","YB3G04","YB4G06","YB4H02","YB2G01","YB2A06","SA3G","YB2A05","YB4F05","YB6A06","YB7A09","YB2G05","YB1G06","2012Env-32","SA7G","E7G","SA10G","SP6G","YB5A06","YB7A06","YB1A01","SL6Y","L6G","W7G","W6G","SP7G","VCR12","I-WASTE-HSH1-TY2","WKB-T9","PKN-T5","DT8","DT7")
df_ALL<-df_ALL[OrderPhylo]

n <- 6
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

pdf("~/Documents/EPFL/Vibrio_Chromosomal_rearrangements/PopGenomics_V2/Manuscript_v1/Notable_genes/NotableGenes_ALL.pdf",height = 12, width=18)
pheatmap(df_ALL,cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,legend = F,annotation_row = notable_genes, annotation_col = population,
         breaks=c(0,0.9,10),color = c("white","black"))
dev.off()
pheatmap(df_ALL,cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,legend = F,annotation_row = notable_genes, 
         breaks=c(0,0.9,10),color = c("white","black"))
