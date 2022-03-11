library(patchwork)
library(tidyverse)
library(qiime2R)
library(ggplot2)

# alpha diversity metrics
# shannon, faith and obversed_features

shannon<-read_qza("shannon_vector.qza")

shannon<-shannon$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

head(shannon)

p<-ggplot(data=shannon, aes(x=SampleID, y=shannon_entropy)) +
  geom_bar(stat="identity")
p

for (i in 1:length(shannon$SampleID) ){
  shannon$sample_type[i] <-  strsplit(shannon$SampleID[i], split='_')[[1]][1]
}

shannon$sample_type [ which(shannon$sample_type=='BIO') ] <- 'Biofilm'
shannon$sample_type [ which(shannon$sample_type=='CON') ] <- 'Earth control Tissue sample'
shannon$sample_type [ which(shannon$sample_type=='ISS') ] <- 'ISS Tissue sample'
shannon$sample_type [ which(shannon$sample_type=='SK') ] <- 'Skin swab'
shannon$sample_type [ which(shannon$sample_type=='TIS') ] <- 'Tissue'


p1 <- ggplot(shannon, aes(x=sample_type, y=shannon_entropy)) + 
  geom_boxplot(alpha=0.5) + geom_point(aes(color=sample_type),size=0.6)+  theme_q2r()+   theme(legend.title=element_blank())+
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ylab('Shannon diversity')+xlab('')
p1
###############################
observed_features<-read_qza("observed_features_vector.qza")

observed_features<-observed_features$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

head(observed_features)

p<-ggplot(data=observed_features, aes(x=SampleID, y=observed_features)) +
  geom_bar(stat="identity")
p

for (i in 1:length(observed_features$SampleID) ){
  
  observed_features$sample_type[i] <-  strsplit(observed_features$SampleID[i], split='_')[[1]][1]
}

observed_features$sample_type [ which(observed_features$sample_type=='BIO') ] <- 'Biofilm'
observed_features$sample_type [ which(observed_features$sample_type=='CON') ] <- 'Earth control Tissue sample'
observed_features$sample_type [ which(observed_features$sample_type=='ISS') ] <- 'ISS Tissue sample'
observed_features$sample_type [ which(observed_features$sample_type=='SK') ] <- 'Skin swab'
observed_features$sample_type [ which(observed_features$sample_type=='TIS') ] <- 'Tissue'


p2 <- ggplot(observed_features, aes(x=sample_type, y=observed_features)) + 
  geom_boxplot(alpha=0.5) + geom_point(aes(color=sample_type),size=0.6)+  theme_q2r()+   theme(legend.title=element_blank())+
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ylab('Observed features (ASV)')+xlab('')
p2
################################################################################
faith_pd<-read_qza("faith_pd_vector.qza")

faith_pd<-faith_pd$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

head(faith_pd)

p<-ggplot(data=faith_pd, aes(x=SampleID, y=faith_pd)) +
  geom_bar(stat="identity")
p

for (i in 1:length(faith_pd$SampleID) ){
  
  faith_pd$sample_type[i] <-  strsplit(faith_pd$SampleID[i], split='_')[[1]][1]
}

faith_pd$sample_type [ which(faith_pd$sample_type=='BIO') ] <- 'Biofilm'
faith_pd$sample_type [ which(faith_pd$sample_type=='CON') ] <- 'Earth control Tissue sample'
faith_pd$sample_type [ which(faith_pd$sample_type=='ISS') ] <- 'ISS Tissue sample'
faith_pd$sample_type [ which(faith_pd$sample_type=='SK') ] <- 'Skin swab'
faith_pd$sample_type [ which(faith_pd$sample_type=='TIS') ] <- 'Tissue'


p3 <- ggplot(faith_pd, aes(x=sample_type, y=faith_pd)) + 
  geom_boxplot(alpha=0.5) + geom_point(aes(color=sample_type),size=0.6)+  theme_q2r()+   theme(legend.title=element_blank())+
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ylab('Faith\'s phylogenetic diversity')+xlab('')
p3
################################################################################
p1+p2 +p3.      # samples with low sequenced depth are not shown
ggsave('alpha_diversity.pdf', height=3.2, width=11.3)




# Now we plot bray_curtis beta diversity
library(gplots)
beta <-read_qza("bray_curtis_distance_matrix.qza")

annota<-c()
for (i in 1:length( rownames(as.matrix(beta$data) )) ){
  
  annota[i] <-  strsplit( rownames(as.matrix(beta$data) )[i], split='_')[[1]][1]
}

annota [ which(annota=='BIO') ] <- 'Biofilm'
annota [ which(annota=='CON') ] <- 'Earth control Tissue sample'
annota [ which(annota=='ISS') ] <- 'ISS Tissue sample'
annota [ which(annota=='SK') ] <- 'Skin swab'
annota [ which(annota=='TIS') ] <- 'Tissue'

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 5
cols = gg_color_hue(n)

ha = HeatmapAnnotation(Sample_Type = annota,
                       annotation_legend_param = list(Sample_Type = list(title = "Sample Type")), 
                       col = list(Sample_Type = c("Biofilm" = cols[1], "Earth control Tissue sample" = cols[2], "ISS Tissue sample" = cols[3],"Skin swab" = cols[4],"Tissue" = cols[5])))

#heatmap.2(as.matrix( beta$data ) )
library(grid)
library(ComplexHeatmap)
pdf('beta_diversity.pdf', height=4, width=7.2)
Heatmap(as.matrix( beta$data ), top_annotation = ha, heatmap_legend_param = list(
  #at = c(-2, 0, 2),
  #labels = c("low", "zero", "high"),
  title = "Beta diversity",
  legend_height = unit(4, "cm"),
  title_position = "lefttop"
),
row_names_gp = grid::gpar(fontsize = 8),
column_names_gp = grid::gpar(fontsize = 8))

dev.off()

# Supp Fig. 2. 
# Barplot
####################################################
taxonomy <-read.csv("level2.txt")
head(taxonomy)
dim(taxonomy)
metadata<-read.table("merged-metadata.tsv", head=TRUE)
metadata

rownames(taxonomy) <- taxonomy$index
annota<-c()
for (i in 1:length( taxonomy$index) ){
  
  annota[i] <-  strsplit(as.character(taxonomy$index)[i] , split='_')[[1]][1]
}

SVs <- t(taxonomy[,2:22])

SVs<-apply(SVs, 2, function(x) x/sum(x, na.rm=T)*100) #convert to percent

colSums(SVs, na.rm=TRUE)

# remove rows
SVs<-SVs[ which(rownames(SVs)!='Unnamed..2'),]
SVs<-SVs[ which(rownames(SVs)!='barcode.sequence'),]

#SVs$annota<- annota
  
SVsToPlot<-  
  data.frame(MeanAbundance=rowMeans(SVs)) %>% #find the average abundance of a SV
  rownames_to_column("Feature.ID") %>%
  arrange(desc(MeanAbundance)) %>%
  top_n(19, MeanAbundance) %>%
  pull(Feature.ID) #extract only the names from the table

SVs %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  gather(-Feature.ID, key="SampleID", value="Abundance") %>%
  mutate(Feature.ID=if_else(Feature.ID %in% SVsToPlot,  Feature.ID, "Remainder")) %>% #flag features to be collapsed
  group_by(SampleID, Feature.ID) %>%
  summarize(Abundance=sum(Abundance)) %>%
  left_join(metadata) %>%
  mutate(NormAbundance=log10(Abundance+0.01)) %>% # do a log10 transformation after adding a 0.01% pseudocount. Could also add 1 read before transformation to percent
  #left_join(taxonomy) %>%
  mutate(Feature=paste(Feature.ID)) %>%
  mutate(Feature=gsub("[kpcofgs]__", "", Feature)) %>% # trim out leading text from taxonomy string
  ggplot(aes(x=SampleID, y=Feature, fill=NormAbundance)) +
  geom_tile() +
  facet_grid(~ss, scales="free_x", space="free_x") +
  theme_q2r() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_viridis_c(name="log10(% Abundance)") 

ggsave("level2.pdf", height=4, width=11, device="pdf") # save a PDF 3 inches by 4 inches


####################################################
taxonomy <-read.csv("level5.txt")
head(taxonomy)
dim(taxonomy)
metadata<-read.table("merged-metadata.tsv", head=TRUE)
metadata

rownames(taxonomy) <- taxonomy$index
SVs <- t(taxonomy[,2:113])

SVs<-apply(SVs, 2, function(x) x/sum(x, na.rm=T)*100) #convert to percent

colSums(SVs, na.rm=TRUE)

# remove rows
SVs<-SVs[ which(rownames(SVs)!='Unnamed..2'),]
SVs<-SVs[ which(rownames(SVs)!='barcode.sequence'),]


SVsToPlot<-  
  data.frame(MeanAbundance=rowMeans(SVs)) %>% #find the average abundance of a SV
  rownames_to_column("Feature.ID") %>%
  arrange(desc(MeanAbundance)) %>%
  top_n(50, MeanAbundance) %>%
  pull(Feature.ID) #extract only the names from the table

SVs %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  gather(-Feature.ID, key="SampleID", value="Abundance") %>%
  mutate(Feature.ID=if_else(Feature.ID %in% SVsToPlot,  Feature.ID, "Remainder")) %>% #flag features to be collapsed
  group_by(SampleID, Feature.ID) %>%
  summarize(Abundance=sum(Abundance)) %>%
  left_join(metadata) %>%
  mutate(NormAbundance=log10(Abundance+0.01)) %>% # do a log10 transformation after adding a 0.01% pseudocount. Could also add 1 read before transformation to percent
  #left_join(taxonomy) %>%
  mutate(Feature=paste(Feature.ID)) %>%
  mutate(Feature=gsub("[kpcofgs]__", "", Feature)) %>% # trim out leading text from taxonomy string
  ggplot(aes(x=SampleID, y=Feature, fill=NormAbundance)) +
  geom_tile() +
  facet_grid(~ss, scales="free_x", space="free_x") +
  theme_q2r() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_viridis_c(name="log10(% Abundance)")

ggsave("level5.pdf", height=7, width=12, device="pdf") # save a PDF 3 inches by 4 inches


####################################################
taxonomy <-read.csv("level7.txt")
head(taxonomy)
dim(taxonomy)
metadata<-read.table("merged-metadata.tsv", head=TRUE)
metadata

rownames(taxonomy) <- taxonomy$index
SVs <- t(taxonomy[,2:305])

SVs<-apply(SVs, 2, function(x) x/sum(x, na.rm=T)*100) #convert to percent

colSums(SVs, na.rm=TRUE)

# remove rows
SVs<-SVs[ which(rownames(SVs)!='Unnamed..2'),]
SVs<-SVs[ which(rownames(SVs)!='barcode.sequence'),]


SVsToPlot<-  
  data.frame(MeanAbundance=rowMeans(SVs)) %>% #find the average abundance of a SV
  rownames_to_column("Feature.ID") %>%
  arrange(desc(MeanAbundance)) %>%
  top_n(500, MeanAbundance) %>%
  pull(Feature.ID) #extract only the names from the table

SVs %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  gather(-Feature.ID, key="SampleID", value="Abundance") %>%
  mutate(Feature.ID=if_else(Feature.ID %in% SVsToPlot,  Feature.ID, "Remainder")) %>% #flag features to be collapsed
  group_by(SampleID, Feature.ID) %>%
  summarize(Abundance=sum(Abundance)) %>%
  left_join(metadata) %>%
  mutate(NormAbundance=log10(Abundance+0.01)) %>% # do a log10 transformation after adding a 0.01% pseudocount. Could also add 1 read before transformation to percent
  #left_join(taxonomy) %>%
  mutate(Feature=paste(Feature.ID)) %>%
  mutate(Feature=gsub("[kpcofgs]__", "", Feature)) %>% # trim out leading text from taxonomy string
  ggplot(aes(x=SampleID, y=Feature, fill=NormAbundance)) +
  geom_tile() +
  facet_grid(~ss, scales="free_x", space="free_x") +
  theme_q2r() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_viridis_c(name="log10(% Abundance)")

ggsave("level7_all.pdf", height=25, width=15, device="pdf") # save a PDF 3 inches by 4 inches

           
