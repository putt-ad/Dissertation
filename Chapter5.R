# 16S and geochemistry analysis of CPT data - Chapter 5
# 16S dataa sequenced by Yupeng available here: https://drive.google.com/drive/u/1/folders/166G-i813tvktFtEQQG8L-6FNUbd31_6X

cpt.packages <- c("ggplot2",
                  "GGally",
                  "RCurl",
                  "ggmap",
                  "mapview",
                  "sf",
                  "vegan",
                  "tidyverse",
                  "RColorBrewer",
                  "corrplot",
                  "ggcorrplot",
                  "ggsn",
                  "Hmisc",
                  "scales",
                  "plyr",
                  "dplyr",
                  "tidyr",
                  "ggpubr",
                  "ape",
                  "readr",
                  "grid",
                  "reshape2",
                  "viridis",
                  "phyloseq",
                  "ggjoy",
                  "ggtree")

new.cpt.packages <- cpt.packages[!(cpt.packages %in% installed.packages()[,"Package"])]
if(length(new.cpt.packages)) install.packages(new.cpt.packages)

# load all these
lapply(cpt.packages, require, character.only = TRUE)

# install phyloseq installer
source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",
       local = TRUE)

# install most recent version of phyloseq from github 
install_phyloseq(branch = "github")

# install phyloseq tool Biocmanager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")
install.packages('metacoder')


# make packages available in your local library
library(BiocManager)
library(microViz)
library(phyloseq)
library(metacoder)
library(microbial)
library(DESeq2)
library(DESeq2)
install.packages("Surrogate")
library(Surrogate)
install.packages("ggtern")
library(ggtern)
BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)
BiocManager::install("MatrixGenerics")

library(MatrixGenerics)

# analysis of 16S
# read in ASV matrix
CPT_otu <- read.csv("https://www.dropbox.com/s/w92c56ho1i1025t/210530_ASV.csv?dl=1", stringsAsFactors = FALSE)
CPT_otu <- data.frame(CPT_otu[,-1], row.names=CPT_otu[,1]) # set ASVID column to rownames
head(CPT_otu)
#read in tree
cpt_tree <- read_tree("https://www.dropbox.com/s/kfsq5uw7zg19a59/210530.ASV.tree.nwk?dl=1")

# read in taxa classifier
CPT_taxa <- read.csv("https://www.dropbox.com/s/f8zs2s213ecwwfa/210530.ASV.Classification.txt?dl=1", stringsAsFactors = FALSE)
glimpse(CPT_taxa)
CPT_taxa <- data.frame(CPT_taxa[,-1], row.names=CPT_taxa[,1]) # set ASVID column to rownames
CPT_taxa <- as.matrix(CPT_taxa)
head(CPT_taxa)
# sediment sample metadata
cpt_metadata <- read.csv("https://www.dropbox.com/s/mppx4lka77rw7dy/chem_transect_metadata.csv?dl=1", stringsAsFactors = FALSE, row.names=1)

#groundwater metadata
sso_metadata <- read.csv("https://www.dropbox.com/s/4dup3pkqmw5ydtg/SSO_nitrate_model_all_chem.csv?dl=1", stringsAsFactors = FALSE, row.names = 1)


glimpse(cpt_metadata)
head(cpt_metadata)

#plot soil texture
#####
ggtern(data = cpt_metadata,
       aes(x=percent_sand, y=percent_clay, z=percent_silt, color = as.factor(HydraulicConductivity))) + geom_point() +
  labs(
    yarrow = "Clay (%)",
    zarrow = "Silt (%)",
    xarrow = "Sand (%)"
  )+
  theme_showarrows()+        # Add arrows to axis titles
  theme_hidetitles()+
  theme_clockwise() +
  labs(color = "Permeability Rank")
#ggsave("~/soil_triangle_SBT.jpeg", bg = "white", width = 70, height = 45, units = "cm")
ggsave("~/soil_triangle_permeability.jpeg", bg = "white", width = 70, height = 45, units = "cm")
#####
# making phyloseq objects
#####
head(CPT_TAX)
CPT_OTU = otu_table(CPT_otu, taxa_are_rows = TRUE)
CPT_OTU

tax.clean <-CPT_taxa#taxa have been cleaned to remove NA

tax.clean 
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    Genus <- paste("Genus_", tax.clean[i,6], sep = "")
    tax.clean[i, 7:8] <- Genus
  }
}

tax.clean
CPT_TAX = phyloseq::tax_table(tax.clean)
CPT_TAX <- CPT_taxa
cpt_metadata[is.na(cpt_metadata)] = 0
CPT_sampledata = sample_data(as.data.frame(cpt_metadata))
head(CPT_sampledata)
cptphyseq <- phyloseq(CPT_OTU, CPT_TAX, CPT_sampledata)
cpt_tree = rtree(ntaxa(cptphyseq), rooted=TRUE, tip.label=taxa_names(cptphyseq))
merged_CPT<- phyloseq(CPT_OTU, CPT_TAX, CPT_sampledata, cpt_tree)
merged_CPT = merge_samples(merged_CPT, "Sed_short_ID")
merged_tree = rtree(ntaxa(merged_CPT ), rooted=TRUE, tip.label=taxa_names(merged_CPT))
merged_CPT<- phyloseq(CPT_OTU, CPT_TAX, CPT_sampledata, merged_tree)
cptphyseq
shallow_cpt <- subset_samples(cptphyseq, water_pH != 4)
sample_data(shallow_cpt)
head(sample_data(cptphyseq))

#plot(cpt_tree)
physeqCPT <- phyloseq(CPT_OTU, CPT_TAX, CPT_sampledata,cpt_tree)
physeqCPT
cpt_patesci <- subset_taxa(physeqCPT, Phylum=="Patescibacteria")
merged_cpt_patesci = merge_samples(cpt_patesci, "Sed_short_ID")
melt_cpt_patesci <- as.data.frame(phyloseq::psmelt(cpt_patesci))
write.csv(melt_cpt_patesci, "~/melt_cpt_patesci_otu.csv")
head(melt_cpt_patesci)
melt_merged_cpt_patesci 
library(tidyr)
patesci_summary <- melt_cpt_patesci  %>% # Summary by group using dplyr
  group_by(melt_merged_cpt_patesci$Sample) %>% 
  summarize(sum = sum(melt_merged_cpt_patesci$Abundance))
cpt_patesci_otu <- otu_table(cpt_patesci)
write.csv(cpt_patesci_otu, "~/cpt_patesci_otu.csv")
sum(otu_table(cpt_patesci))/1035837
cpt_proteobac <- subset_taxa(physeqCPT, Phylum=="Proteobacteria")
proteo_genus<- phyloseq::tax_glom(cpt_proteobac, "Genus", NArm = TRUE)
cpt_proteo_otu <- otu_table(cpt_proteobac)
melt_cpt_proteobac <- as.data.frame(phyloseq::psmelt(cpt_proteobac))
write.csv(melt_cpt_proteobac, "~/cpt_proteobacteria_otu.csv")
cpt_firm <- subset_taxa(physeqCPT, Phylum=="Firmicutes")
sum(otu_table(cpt_firm))
melt_cpt_firm <- as.data.frame(phyloseq::psmelt(cpt_firm))
write.csv(melt_cpt_firm, "~/melt_cpt_firmicutes_otu.csv")
cpt_arch <- subset_taxa(physeqCPT, Domain=="Archaea")
View(otu_table(cpt_arch))
melt_cpt_arch <- as.data.frame(phyloseq::psmelt(cpt_arch))
write.csv(melt_cpt_arch, "~/melt_cpt_archaea_otu.csv")
cpt_WPS <- subset_taxa(physeqCPT, Phylum=="WPS-2")
View(otu_table(cpt_WPS))
melt_cpt_WPS <- as.data.frame(phyloseq::psmelt(cpt_WPS))
write.csv(melt_cpt_WPS, "~/melt_cpt_WPS-2_otu.csv")
cpt_acido <- subset_taxa(physeqCPT, Phylum=="Acidobacteria")
melt_cpt_acido <- as.data.frame(phyloseq::psmelt(cpt_acido))
write.csv(melt_cpt_acido, "~/melt_cpt_acidobacteria_otu.csv")
cpt_nitrospirae <- subset_taxa(physeqCPT, Phylum=="Nitrospirae")
melt_cpt_nitrospirae <- as.data.frame(phyloseq::psmelt(cpt_nitrospirae))
write.csv(melt_cpt_nitrospirae, "~/melt_cpt_nitrospirae_otu.csv")


View(phyloseq::tax_table(cpt_arch))
sum(otu_table(cpt_arch))

phyloseq::tax_table(physeqCPT)
sample_data(physeqCPT)$depth_range <-as.factor(case_when((sample_data(physeqCPT)$mid_depth_cm > 250) & (sample_data(physeqCPT)$mid_depth_cm <= 450) ~ "201-450",
                                                         (sample_data(physeqCPT)$mid_depth_cm > 450) & (sample_data(physeqCPT)$mid_depth_cm <= 650) ~ "451-650",
                                                         (sample_data(physeqCPT)$mid_depth_cm > 650) & (sample_data(physeqCPT)$mid_depth_cm <= 850) ~ "601-850",
                                                         (sample_data(physeqCPT)$mid_depth_cm > 850) & (sample_data(physeqCPT)$mid_depth_cm <= 1000) ~ "851-1000"))

mergedCPT = merge_samples(physeqCPT, "Sed_short_ID")
sample_data(mergedCPT)$depth_range <-as.factor(case_when((sample_data(mergedCPT)$mid_depth_cm > 250) & (sample_data(mergedCPT)$mid_depth_cm <= 450) ~ "201-450",
                                                         (sample_data(mergedCPT)$mid_depth_cm > 450) & (sample_data(mergedCPT)$mid_depth_cm <= 650) ~ "451-650",
                                                         (sample_data(mergedCPT)$mid_depth_cm > 650) & (sample_data(mergedCPT)$mid_depth_cm <= 850) ~ "601-850",
                                                         (sample_data(mergedCPT)$mid_depth_cm > 850) & (sample_data(mergedCPT)$mid_depth_cm <= 1000) ~ "851-1000"))
cpt_otu.condensed <- as.data.frame(otu_table(mergedCPT))
glimpse(cpt_otu.condensed)
sample_data(mergedCPT)$Sed_short_ID <- row.names(sample_data(mergedCPT))
physeqCPT <- mergedCPT
phyloseq_coverage(physeqCPT)

# rarefy by even depth. Removed 9254 OTUs
set.seed(0451)
physeqCPT_rarefied<-rarefy_even_depth(physeqCPT,sample.size = min(sample_sums(physeqCPT)))
physeqCPT_rarefied

#end phyloseq
#####

#compare differences in alpha statistics
#####
install.packages('sjPlot')
library(sjPlot)
#shannon changes by day
#calculate sample alpha stats
#full 16S data set
mergedCPT_sample_data <- as.data.frame(sample_data(mergedCPT))
head(mergedCPT_sample_data)
colnames(mergedCPT_sample_data)
row.names(mergedCPT_sample_data)
rich_mergedCPT = estimate_richness(mergedCPT, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))
row.names(rich_mergedCPT) <-row.names(mergedCPT_sample_data)
merged_stat_mergedCPT <-as.data.frame(merge(rich_mergedCPT, mergedCPT_sample_data, by=0, all=TRUE))
write.csv(merged_stat_mergedCPT, "~/cpt_alpha_stats.csv")
merged_stat_mergedCPT<- read.csv( "~/cpt_alpha_stats.csv")
#patescibacteria data set
merged_cpt_patesci_sample_data <- as.data.frame(sample_data(merged_cpt_patesci))
head(merged_cpt_patesci_sample_data)
colnames(merged_cpt_patesci_sample_data)
row.names(merged_cpt_patesci_sample_data)
rich_merged_cpt_patesci = estimate_richness(merged_cpt_patesci, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))
row.names(rich_merged_cpt_patesci) <-row.names(merged_cpt_patesci_sample_data)
merged_stat_merged_cpt_patesci <-as.data.frame(merge(rich_merged_cpt_patesci, merged_cpt_patesci_sample_data, by=0, all=TRUE))
write.csv(merged_stat_merged_cpt_patesci, "~/cpt_patesci_alpha_stats.csv")


head(merged_stat_mergedCPT)
merged_stat_1$Acetate_range <-case_when((merged_stat_1$Acetate <= 1)              ~ 1,
                                        (merged_stat_1$Acetate  > 1) & (merged_stat_1$Acetate <= 10) ~ 10,
                                        (merged_stat_1$Acetate  > 10) & (merged_stat_1$Acetate <= 20) ~ 20,
                                        (merged_stat_1$Acetate  > 20) & (merged_stat_1$Acetate <= 30) ~ 30,
                                        (merged_stat_1$Acetate  > 30) & (merged_stat_1$Acetate<= 50) ~ 40)
#Simpson and Shannon stats for well types and days 
summary(aov(data=merged_stat_1, Shannon~depth_range*water_pH+Nitrate+Sulfate+fine+Coarse+))
summary(aov(data=merged_stat_1, ACE~well_type*DayCategory+nitrate_level+sulfate_level))
summary(aov(data=merged_stat_1, Shannon~well_type*Days+nitrate_level+sulfate_level))
p1_shannon <- aov(data=merged_stat_1, Shannon~well_type*Days+nitrate_level+sulfate_level)
summary(p1_shannon)
summary(aov(data=merged_stat_1, Fisher~well_type*DayCategory+nitrate_level+sulfate_level))

#0.2 only
p2017_2_sample_data <- as.data.frame(sample_data(p2017_2_filtered))
p2017_2_sample_data$N_cat <- as.factor(p2017_2_sample_data$nitrate_level);p2017_2_sample_data$S_cat <- as.factor(p2017_2_sample_data$sulfate_level)
colnames(p2017_2_sample_data)
row.names(p2017_2_sample_data)
rich_2 = estimate_richness(p2017_2_filtered, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))
plot_richness(p2017_2_filtered, measures = c( "Chao1", "ACE", "Shannon","InvSimpson"))
row.names(rich_2) <-row.names(p2017_2_sample_data)
merged_stat_2 <-as.data.frame(merge(rich_2, p2017_2_sample_data, by=0, all=TRUE))
View(merged_stat_2)
write.csv(merged_stat_2, "~/Dropbox/Andrews MacBook (Hazen MacBook Copy)/Documents/Research/EVO-UMB/merged_2_alpha_stat.csv")

merged_stat$Acetate_level <-case_when((merged_stat$Acetate <= 1)              ~ 1,
                                      (merged_stat$Acetate  > 1) & (merged_stat$Acetate <= 10) ~ 10,
                                      (merged_stat$Acetate  > 10) & (merged_stat$Acetate <= 20) ~ 20,
                                      (merged_stat$Acetate  > 20) & (merged_stat$Acetate <= 30) ~ 30,
                                      (merged_stat$Acetate  > 30) & (merged_stat$Acetate <= 50) ~ 40)
merged_2_stat <- merged_stat; head(merged_2_stat)
#Simpson and Shannon stats for well types and days 
anova(data=merged_stat_2, Shannon~Days*well_type_factor)
summary(aov(data=merged_stat_2, Shannon~well_type*DayCategory+nitrate_level+sulfate_level))
summary(aov(data=merged_stat_2, Chao1~well_type*DayCategory+nitrate_level+sulfate_level))
summary(aov(data=merged_stat_2, ACE~well_type*DayCategory+nitrate_level+sulfate_level))
summary(aov(data=merged_stat_2, Shannon~well_type*DayCategory+nitrate_level+sulfate_level))
p2_shannon <- aov(data=merged_stat_2, Shannon~well_type*DayCategory+nitrate_level+sulfate_level)
View(tidy(TukeyHSD(p2_shannon, which ='DayCategory')))

summary(aov(data=merged_stat_2, Fisher~well_type*DayCategory+nitrate_level+sulfate_level))
p2_simpson<-aov(data=merged_stat_2, InvSimpson~well_type*DayCategory+nitrate_level+sulfate_level)
View(tidy(TukeyHSD(p2_simpson, which ='well_type:DayCategory')))


#COMBINED
library(phyloseq)
p2017_combined_filtered <- p2017_combined_filtered %>%
  ps_mutate(Well = if_else(well_type == "Monitoring", true = 1, false = 0),
            Filter = if_else(FilterSize == "0.1um", true = 1, false = 2))

p2017_C_sample_data <- as.data.frame(sample_data(p2017_combined_filtered))
colnames(p2017_C_sample_data)
row.names(p2017_C_sample_data)
rich_c = estimate_richness(p2017_combined_filtered, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))
row.names(rich_c) <-row.names(p2017_C_sample_data)
merged_C_stat <-as.data.frame(merge(rich_c, p2017_C_sample_data, by=0, all=TRUE))
View(merged_C_stat)
write.csv(merged_C_stat, "~/Dropbox/Andrews MacBook (Hazen MacBook Copy)/Documents/Research/EVO-UMB/small_bacteria_alpha_stat.csv")

#merged_combined
p2017_Cm_sample_data <- as.data.frame(sample_data(p2017_combined_merged))
colnames(p2017_Cm_sample_data)
row.names(p2017_Cm_sample_data)
rich_cm = estimate_richness(p2017_combined_merged, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))
row.names(rich_cm) <-row.names(p2017_Cm_sample_data)
merged_Cm_stat <-as.data.frame(merge(rich_cm, p2017_Cm_sample_data, by=0, all=TRUE))
head(merged_Cm_stat)
write.csv(merged_Cm_stat, "~/Dropbox/Andrews MacBook (Hazen MacBook Copy)/Documents/Research/EVO-UMB/small_bacteria_mergedsamples_alpha_stat.csv")

merged_stat$Acetate_level <-case_when((merged_stat$Acetate <= 1)              ~ 1,
                                      (merged_stat$Acetate  > 1) & (merged_stat$Acetate <= 10) ~ 10,
                                      (merged_stat$Acetate  > 10) & (merged_stat$Acetate <= 20) ~ 20,
                                      (merged_stat$Acetate  > 20) & (merged_stat$Acetate <= 30) ~ 30,
                                      (merged_stat$Acetate  > 30) & (merged_stat$Acetate <= 50) ~ 40)
merged_C_stat <- merged_stat
summary(aov(data=merged_C_stat, Shannon~well_type*DayCategory+FilterSize+nitrate_level+sulfate_level))
summary(aov(data=merged_C_stat, ACE~well_type*DayCategory+FilterSize+nitrate_level+sulfate_level))
summary(aov(data=merged_C_stat, Chao1~well_type*DayCategory+FilterSize+nitrate_level+sulfate_level))
summary(aov(data=merged_C_stat, Shannon~well_type*DayCategory+FilterSize+nitrate_level+sulfate_level))
pc_shannon <-aov(data=merged_C_stat, Shannon~well_type*DayCategory+FilterSize+nitrate_level+sulfate_level)
View(tidy(TukeyHSD(pc_shannon, which ='well_type:DayCategory')))
summary(aov(data=merged_C_stat, Fisher~well_type*DayCategory+FilterSize+nitrate_level+sulfate_level))


#merged
merged_p2017_stat <- merged_p2017 %>%
  microViz::ps_mutate(Well = if_else(well_type == "Monitoring", true = 1, false = 0))

p2017_merge_sample_data <- as.data.frame(sample_data(merged_p2017))
View(p2017_merge_sample_data)
colnames(p2017_merge_sample_data)
row.names(p2017_merge_sample_data)
rich_m = estimate_richness(merged_p2017, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))
rich_m
row.names(rich_m) <-row.names(p2017_merge_sample_data)
p2017_merged_stat <-as.data.frame(merge(rich_m, p2017_merge_sample_data, by=0, all=TRUE))
p2017_merged_stat <- p2017_merged_stat %>%
  mutate(Well = if_else(well_type == "Monitoring", true = 1, false = 0))
head(p2017_merged_stat)
View(p2017_merged_stat)
write.csv(p2017_merged_stat, "~/Dropbox/Andrews MacBook (Hazen MacBook Copy)/Documents/Research/EVO-UMB/p2017_merged_alpha_stat.csv")

p2017_merged_stat$Acetate_level <-case_when((p2017_merged_stat$Acetate <= 1)              ~ 1,
                                            (p2017_merged_stat$Acetate  > 1) & (p2017_merged_stat$Acetate <= 10) ~ 10,
                                            (p2017_merged_stat$Acetate  > 10) & (p2017_merged_stat$Acetate <= 20) ~ 20,
                                            (p2017_merged_stat$Acetate  > 20) & (p2017_merged_stat$Acetate <= 30) ~ 30,
                                            (p2017_merged_stat$Acetate  > 30) & (p2017_merged_stat$Acetate <= 50) ~ 40)
head(p2017_merged_stat)
#Simpson and Shannon stats for well types and days 
summary(aov(data=p2017_merged_stat, Chao1~well_type*DayCategory+nitrate_level+sulfate_level))
summary(aov(data=p2017_merged_stat, ACE~well_type*DayCategory+nitrate_level+sulfate_level))
summary(aov(data=p2017_merged_stat, Shannon~well_type*DayCategory+nitrate_level+sulfate_level))
summary(aov(data=p2017_merged_stat, Fisher~well_type*DayCategory+nitrate_level+sulfate_level))
summary(aov(data=p2017_merged_stat, InvSimpson~well_type*DayCategory+nitrate_level+sulfate_level))

#simpson and shannon stats for well type, days 0.1 filter
head(merged_C_stat)
merged_p2017_stat_0.1 <- subset(merged_C_stat, FilterSize == "0.1um")
summary(aov(data=merged_p2017_stat_0.1, Shannon~well_type_factor*Days+nitrate_level+sulfate_level+Acetate_level +phase_cat))
summary(aov(data=merged_p2017_stat_0.1, InvSimpson~well_type_factor*Days+nitrate_level+sulfate_level+Acetate_level))

#simpson and shannon stats for well type, days 0.2 filter
merged_p2017_stat_0.2 <- subset(merged_C_stat, FilterSize == "0.2um")
summary(aov(data=merged_p2017_stat_0.2, Shannon~well_type_factor*Days+nitrate_level+sulfate_level+Acetate_level))
summary(aov(data=merged_p2017_stat_0.2, InvSimpson~well_type_factor*Days+nitrate_level+sulfate_level+Acetate_level))



#merged 0.1 and 0.2
merged_p2017_sample_data <- as.data.frame(sample_data(merged_phylo2017_rarefied))
colnames(merged_p2017_sample_data)
row.names(merged_p2017_sample_data)
rich = estimate_richness(merged_phylo2017_rarefied, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))
row.names(rich) <-row.names(merged_p2017_sample_data)
merged_stat <-as.data.frame(merge(rich, merged_p2017_sample_data, by=0, all=TRUE))
View(merged_stat)
merged_stat$Acetate_range <-case_when((merged_stat$Acetate_mgL <= 1)              ~ 1,
                                      (merged_stat$Acetate_mgL  > 1) & (merged_stat$Acetate_mgL <= 10) ~ 10,
                                      (merged_stat$Acetate_mgL  > 10) & (merged_stat$Acetate_mgL <= 20) ~ 20,
                                      (merged_stat$Acetate_mgL  > 20) & (merged_stat$Acetate_mgL <= 30) ~ 30,
                                      (merged_stat$Acetate_mgL  > 30) & (merged_stat$Acetate_mgL<= 50) ~ 40)

#summary plots for alpha stats
p2017_merged_shannon <- plot_richness(p2017_1_filtered, x="Days", color="well_type_factor", measures = "Shannon")
#p2017_merged_shannon + 
merged_C_stat$Acetate_level
hist(merged_C_stat$Acetate_level)
ggplot(data=merged_C_stat, aes(x=Acetate_level, y = Shannon, color = FilterSize))+
  geom_point() +
  geom_line()+
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  labs(y = "Shannon Index", x = "Acetate", color = "Well Type", size=" ")+
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white", color=NA),
        legend.direction="horizontal",
        legend.position = c(0.46, 0.92),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        plot.title = element_text(size=18))
ggsave("~/Dropbox/Andrews MacBook (Hazen MacBook Copy)/Documents/Research/EVO-UMB/p2017_lm_day_well_shannon_plot.pdf", height = 6, width = 6)

merged_C_stat 
library(dplyr)
#A. shannon by filter size plot
shannon_summary <- merged_C_stat %>% # Summary by group using dplyr
  group_by(FilterSize, Days) %>% 
  summarise(mean = mean(Shannon),
            std.dev = sd(Shannon),
            stdmax = mean(Shannon)+sd(Shannon),
            stdmin = mean(Shannon)-sd(Shannon),
            range = max(Shannon)-min(Shannon))
shannon_summary[is.na(shannon_summary)] = 0
head(shannon_summary)
shannon_plot <-ggplot( shannon_summary, aes(x=Days) ) + 
  geom_line( aes(y=mean,group=FilterSize,color=FilterSize) ) + 
  geom_point( aes(y=mean,color = FilterSize, shape = FilterSize),size=4, alpha = 0.8) +
  geom_linerange(aes(ymin = stdmin , ymax = stdmax, color = FilterSize), linetype = 3, alpha = 0.3) +
  geom_line( aes(x=Days,y=stdmax,group=FilterSize,color=FilterSize), linetype = 2, alpha = 0.5) + 
  geom_line( aes(x=Days,y=stdmin,group=FilterSize,color=FilterSize), linetype = 2, alpha=0.5) + 
  geom_ribbon(aes(x=Days,ymin=stdmin,ymax=stdmax, fill = FilterSize), alpha=0.1) +
  theme_test() +
  #geom_hline(yintercept = 0, linetype = 1, color = 'black', size = 0.5) +
  #geom_vline(xintercept = 0, linetype = 1, color = "gray50", size = 5.5, alpha = 0.4) +
  #geom_text(aes(x=1.8, label="Injection of EVO\n", y=23), colour="gray30", angle=90, text=element_text(size=8)) +
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = c(0.46, 0.92),
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        panel.border = element_rect(colour = "black", size=0.8))+
  scale_x_continuous( breaks = unique(shannon_summary$`Days`)) +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  #xlim(-10,250) +
  labs(x = "Days Post-Injection", y = "Shannon", fill = "",color = "", shape = "", linetype = "", title = "A. Shannon by filter size")
shannon_plot

#B. filter size inverse simpson
InvSimpson_summary <- merged_C_stat %>% # Summary by group using dplyr
  group_by(FilterSize, Days) %>% 
  summarise(mean = mean(InvSimpson),
            std.dev = sd(InvSimpson),
            stdmax = mean(InvSimpson)+sd(InvSimpson),
            stdmin = mean(InvSimpson)-sd(InvSimpson),
            range = max(InvSimpson)-min(InvSimpson))
InvSimpson_summary[is.na(InvSimpson_summary)] = 0
head(InvSimpson_summary)
InvSimpson_plot <-ggplot( InvSimpson_summary, aes(x=Days) ) + 
  geom_line( aes(y=mean,group=FilterSize,color=FilterSize) ) + 
  geom_point( aes(y=mean,color = FilterSize, shape = FilterSize),size=4, alpha = 0.8) +
  geom_linerange(aes(ymin = stdmin , ymax = stdmax, color = FilterSize), linetype = 3, alpha = 0.3) +
  geom_line( aes(x=Days,y=stdmax,group=FilterSize,color=FilterSize), linetype = 2, alpha = 0.5) + 
  geom_line( aes(x=Days,y=stdmin,group=FilterSize,color=FilterSize), linetype = 2, alpha=0.5) + 
  geom_ribbon(aes(x=Days,ymin=stdmin,ymax=stdmax, fill = FilterSize), alpha=0.1) +
  theme_test() +
  #geom_hline(yintercept = 0, linetype = 1, color = 'black', size = 0.5) +
  #geom_vline(xintercept = 0, linetype = 1, color = "gray50", size = 5.5, alpha = 0.4) +
  #geom_text(aes(x=1.8, label="Injection of EVO\n", y=23), colour="gray30", angle=90, text=element_text(size=8)) +
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = c(0.46, 0.92),
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        panel.border = element_rect(colour = "black", size=0.8))+
  scale_x_continuous( breaks = unique(InvSimpson_summary$`Days`)) +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  #xlim(-10,250) +
  labs(x = "Days Post-Injection", y = "Inverse Simpson", fill = "",color = "", shape = "", linetype = "", title = "A. Inverse Simpson by filter size")
InvSimpson_plot

#C. shannon by filter size plot
shannon_well_summary <- merged_C_stat %>% # Summary by group using dplyr
  group_by(well_type, Days) %>% 
  summarise(mean = mean(Shannon),
            std.dev = sd(Shannon),
            stdmax = mean(Shannon)+sd(Shannon),
            stdmin = mean(Shannon)-sd(Shannon),
            range = max(Shannon)-min(Shannon))
shannon_well_summary[is.na(shannon_well_summary)] = 0
head(shannon_well_summary)
shannon_well_plot <-ggplot( shannon_well_summary, aes(x=Days) ) + 
  geom_line( aes(y=mean,group=well_type,color=well_type) ) + 
  geom_point( aes(y=mean,color = well_type, shape = well_type),size=4, alpha = 0.8) +
  geom_linerange(aes(ymin = stdmin , ymax = stdmax, color = well_type), linetype = 3, alpha = 0.3) +
  geom_line( aes(x=Days,y=stdmax,group=well_type,color=well_type), linetype = 2, alpha = 0.5) + 
  geom_line( aes(x=Days,y=stdmin,group=well_type,color=well_type), linetype = 2, alpha=0.5) + 
  geom_ribbon(aes(x=Days,ymin=stdmin,ymax=stdmax, fill = well_type), alpha=0.1) +
  theme_test() +
  #geom_hline(yintercept = 0, linetype = 1, color = 'black', size = 0.5) +
  #geom_vline(xintercept = 0, linetype = 1, color = "gray50", size = 5.5, alpha = 0.4) +
  #geom_text(aes(x=1.8, label="Injection of EVO\n", y=23), colour="gray30", angle=90, text=element_text(size=8)) +
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = c(0.46, 0.92),
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        panel.border = element_rect(colour = "black", size=0.8))+
  scale_x_continuous( breaks = unique(shannon_well_summary$`Days`)) +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  #xlim(-10,250) +
  labs(x = "Days Post-Injection", y = "Shannon", fill = "",color = "", shape = "", linetype = "", title = "C. Shannon by well type")
shannon_well_plot

#D. well inverse simpson
InvSimpson_well_summary <- merged_C_stat %>% # Summary by group using dplyr
  group_by(well_type, Days) %>% 
  summarise(mean = mean(InvSimpson),
            std.dev = sd(InvSimpson),
            stdmax = mean(InvSimpson)+sd(InvSimpson),
            stdmin = mean(InvSimpson)-sd(InvSimpson),
            range = max(InvSimpson)-min(InvSimpson))
InvSimpson_well_summary[is.na(InvSimpson_well_summary)] = 0
head(InvSimpson_well_summary)
InvSimpson_well_plot <-ggplot( InvSimpson_well_summary, aes(x=Days) ) + 
  geom_line( aes(y=mean,group=well_type,color=well_type) ) + 
  geom_point( aes(y=mean,color = well_type, shape = well_type),size=4, alpha = 0.8) +
  geom_linerange(aes(ymin = stdmin , ymax = stdmax, color = well_type), linetype = 3, alpha = 0.3) +
  geom_line( aes(x=Days,y=stdmax,group=well_type,color=well_type), linetype = 2, alpha = 0.5) + 
  geom_line( aes(x=Days,y=stdmin,group=well_type,color=well_type), linetype = 2, alpha=0.5) + 
  geom_ribbon(aes(x=Days,ymin=stdmin,ymax=stdmax, fill = well_type), alpha=0.1) +
  theme_test() +
  #geom_hline(yintercept = 0, linetype = 1, color = 'black', size = 0.5) +
  #geom_vline(xintercept = 0, linetype = 1, color = "gray50", size = 5.5, alpha = 0.4) +
  #geom_text(aes(x=1.8, label="Injection of EVO\n", y=23), colour="gray30", angle=90, text=element_text(size=8)) +
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = c(0.46, 0.92),
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        panel.border = element_rect(colour = "black", size=0.8))+
  scale_x_continuous( breaks = unique(shannon_summary$`Days`)) +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  #xlim(-10,250) +
  labs(x = "Days Post-Injection", y = "Inverse Simpson", fill = "",color = "", shape = "", linetype = "", title = "D. Inverse Simpson by well type")
InvSimpson_well_plot 

#test differences in shannon means with linear model
tab_model(lm(value ~DayCategory*well_type_factor, p2017_merged_shannon$data))
ggsave("~/Dropbox/Andrews MacBook (Hazen MacBook Copy)/Documents/Research/EVO-UMB/p2017_lm_day_well_shannon.pdf", height = 6, width = 6)
#chao1 changes by day
p2017_merged_chao1 <- plot_richness(merged_phylo2017_rarefied, x="DayCategory", color="well_type_factor", measures = "Chao1")
p2017_merged_chao1
#test differences in simpson means with linear model
tab_model(lm(value ~DayCategory*well_type_factor, p2017_merged_chao1$data))
ggsave("~/Dropbox/Andrews MacBook (Hazen MacBook Copy)/Documents/Research/EVO-UMB/p2017_lm_day_well_chao1.pdf", height = 6, width = 6)
#invSimp. changes by day
p2017_merged_InvSimpson <- plot_richness(merged_phylo2017_rarefied, x="DayCategory", color="well_type_factor", measures = "InvSimpson")
p2017_merged_InvSimpson+  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  labs(y = "Inverse Simpson Index", x = "Days", color = "Well Type", size=" ")+
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white", color=NA),
        legend.direction="horizontal",
        legend.position = c(0.46, 0.92),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        plot.title = element_text(size=18))
ggsave("~/Dropbox/Andrews MacBook (Hazen MacBook Copy)/Documents/Research/EVO-UMB/p2017_lm_day_well_invsimpson_plot.pdf", height = 6, width = 6)

#test differences in simpson means with linear model
tab_model(lm(value ~DayCategory*well_type_factor, p2017_merged_InvSimpson$data))
ggsave("~/Dropbox/Andrews MacBook (Hazen MacBook Copy)/Documents/Research/EVO-UMB/p2017_lm_day_well_invsimpson.pdf", height = 6, width = 6)

#InvSimp. by nitrate
p2017_merged_NO3_shannon <- plot_richness(merged_phylo2017_rarefied, x="nitrate_range", color="well_type_factor", measures = "Shannon")
aov(value ~nitrate_range*well_type_factor, data=p2017_merged_NO3_shannon$data)
tab_model(lm(value ~nitrate_range*well_type_factor, data=p2017_merged_NO3_shannon$data))
ggsave("~/Dropbox/Andrews MacBook (Hazen MacBook Copy)/Documents/Research/EVO-UMB/p2017_lm_nitrate_well_shannon.pdf", height = 6, width = 6)

p2017_merged_SO4_shannon <- plot_richness(merged_phylo2017_rarefied, x="sulfate_range", color="well_type_factor", measures = "Shannon")
aov(value ~sulfate_range*well_type_factor, data=p2017_merged_SO4_shannon$data)
tab_model(lm(value ~sulfate_range*well_type_factor, p2017_merged_SO4_shannon$data))
ggsave("~/Dropbox/Andrews MacBook (Hazen MacBook Copy)/Documents/Research/EVO-UMB/p2017_lm_nitrate_well_shannon.pdf", height = 6, width = 6)

p2017_merged_U238_shannon <- plot_richness(merged_phylo2017_rarefied, x="uranium_reduction_cat", color="well_type_factor", measures = "Shannon")
tab_model(lm(value ~uranium_reduction_cat*well_type_factor, p2017_merged_U238_shannon$data))


p2017_merged_NO3_shannon+  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  labs(y = "Shannon Index", x = "Nitrate", color = "Well Type", size=" ")+
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white", color=NA),
        legend.direction="horizontal",
        legend.position = c(0.46, 0.92),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        plot.title = element_text(size=18))
ggsave("~/Dropbox/Andrews MacBook (Hazen MacBook Copy)/Documents/Research/EVO-UMB/p2017_lm_day_well_invsimpson_plot.pdf", height = 6, width = 6)

#test differences in simpson means with linear model
tab_model(lm(value ~nitrate_range*well_type_factor, p2017_merged_NO3_shannon$data))
ggsave("~/Dropbox/Andrews MacBook (Hazen MacBook Copy)/Documents/Research/EVO-UMB/p2017_lm_day_well_invsimpson.pdf", height = 6, width = 6)
#end of alpha stats
#####


#make abudnance plot
#####
sum(otu_table(physeqCPT))
sum(otu_table(physeqCPT)) * 0.0005
cpt_phy <-phyloseq::prune_taxa(taxa_sums(physeqCPT) > 500, physeqCPT)
table(phyloseq::tax_table(cpt_phy)[, "Phylum"])
cpt_phy <- phyloseq::tax_glom(cpt_phy, "Phylum")
melt_cpt_phy <- phyloseq::psmelt(cpt_phy)
head(melt_cpt_phy)

color_list_cpt <- list(   "pink", "plum3", "salmon", "orange3", "goldenrod3","lightgoldenrod2" , "sienna4","tan2" , "maroon"  ,"plum3",
                          "mediumpurple", "peachpuff3", "navajowhite4", "gray20",       
                          "lightblue3",  "steelblue", "lightcyan"  ,  "lightskyblue4",  "royalblue4" ,   
                          "turquoise4" ,  "slategray", "red",  "green")

cpt_phy.plot = ggplot(data = melt_cpt_phy, mapping = aes(x=Sed_short_ID, y=Abundance, fill=Phylum)) +
  #geom_bar(stat = "identity", colour = "black") + 
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, size = 17, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 18), legend.title = element_text(size = 18), 
        legend.text = element_text(colour = "black"), 
        axis.text.y = element_text(colour = "black")) + 
  scale_fill_manual(values = color_list_cpt) +
  scale_y_continuous(expand = c(0,0)) + 
  #facet_grid(FilterSize.x ~factor(Days.x, levels=c('-6','1','8','15','22', '50','78', '106', '134'))) +
  #facet_grid(~zone) +
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = "bottom",
        legend.text = element_text(size=12, color = "black"),
        axis.title = element_text(size=17, color = "black"),
        axis.text = element_text(size = 17, color = "black"),
        panel.border = element_rect(colour = "black", fill = "transparent", size=0.8))+
  labs(x = "", y = "ASV Count", fill = "", title = "Abundant Phyla in CPT samples")
cpt_phy.plot
ggsave("~/phylum_plot.tiff", bg = "white", width = 40, height = 30, units = "cm")



cpt_patesci <- subset_taxa(physeqCPT, Phylum=="Patescibacteria")
sum(otu_table(cpt_patesci))*0.005
cpt_patesci <-phyloseq::prune_taxa(taxa_sums(cpt_patesci) > 77, cpt_patesci)
#p2019_patesci <-phyloseq::prune_taxa(taxa_sums(p2019_patesci) > 367, p2019_patesci)
table(phyloseq::tax_table(cpt_patesci)[, "Order"])
cpt_patesci <- phyloseq::tax_glom(cpt_patesci, "Order")
melt_cpt_patesci <- phyloseq::psmelt(cpt_patesci)
head(melt_cpt_patesci)


cpt_patesci.plot = ggplot(data = melt_cpt_patesci, mapping = aes(x=Sed_short_ID, y=Abundance, fill=Order)) +
  #geom_bar(stat = "identity", colour = "black") + 
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, size = 17, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 18), legend.title = element_text(size = 18), 
        legend.text = element_text(colour = "black"), 
        axis.text.y = element_text(colour = "black")) + 
  scale_fill_manual(values = color_list_cpt) +
  scale_y_continuous(expand = c(0,0)) + 
  #facet_grid(FilterSize.x ~factor(Days.x, levels=c('-6','1','8','15','22', '50','78', '106', '134'))) +
  #facet_grid(FilterSize.x~well.x) +
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = "bottom",
        legend.text = element_text(size=11.7, color = "black"),
        axis.title = element_text(size=17, color = "black"),
        axis.text = element_text(size = 17, color = "black"),
        panel.border = element_rect(colour = "black", fill = "transparent", size=0.8))+
  labs(x = "", y = "ASV Count", fill = "", title = "Abundant Patescibacteria in CPT samples")
cpt_patesci.plot
ggsave("~/patesci_>0.5_plot.tiff", bg = "white", width = 40, height = 30, units = "cm")

#archaea community composition
sum(otu_table(cpt_arch))
sum(otu_table(cpt_arch)) * 0.0005
arch_ord <-phyloseq::prune_taxa(taxa_sums(cpt_arch) > 50, cpt_arch)
table(phyloseq::tax_table(arch_ord)[, "Order"])
arch_ord <- phyloseq::tax_glom(arch_ord, "Order")
melt_arch_ord <- phyloseq::psmelt(arch_ord)
#head(melt_cpt_phy)

color_list_arch <- list(   "pink", "plum3", "salmon", "orange3", "brown","goldenrod3","lightgoldenrod2" , "sienna4","tan2" , "maroon"  ,"plum3",
                           "mediumpurple", "peachpuff3", "navajowhite4", "gray20",       
                           "lightblue3",  "steelblue", "lightcyan"  ,  "lightskyblue4",  "royalblue4" ,   
                           "turquoise4" ,  "seagreen", "darkgreen","slategray", "gray30")

cpt_arch.plot = ggplot(data = melt_arch_ord, mapping = aes(x=Sed_short_ID, y=Abundance, fill=Order)) +
  #geom_bar(stat = "identity", colour = "black") + 
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, size = 15, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 18), legend.title = element_text(size = 18), 
        legend.text = element_text(colour = "black"), 
        axis.text.y = element_text(colour = "black")) + 
  scale_fill_manual(values = color_list_arch) +
  scale_y_continuous(expand = c(0,0)) + 
  #facet_grid(FilterSize.x ~factor(Days.x, levels=c('-6','1','8','15','22', '50','78', '106', '134'))) +
  #facet_grid(~zone) +
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = "bottom",
        legend.text = element_text(size=12, color = "black"),
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        panel.border = element_rect(colour = "black", fill = "transparent", size=0.8))+
  labs(x = "", y = "ASV Count", fill = "", title = "Abundant archaeal orders in CPT samples")
cpt_arch.plot
ggsave("~/achaea-plot.tiff", bg = "white", width = 40, height = 30, units = "cm")


#end of abudnance plots
#####

#taxonomic correlations with chemistry
#####

#rename col labels
colnames(sample_data(physeqCPT))
phyloseq::tax_table(physeqCPT)
phyloseq::otu_table(physeqCPT)
head(sample_data(shallow_cpt))
#shallow CPT correlations
shallow_cpt
shallow_class <- shallow_cpt
sample_data(shallow_class)$subsection <- dplyr::recode(sample_data(shallow_class)$section, "A" = 1, "B" = 2, "C"= 3)
sample_data(shallow_class)$well <- dplyr::recode(sample_data(shallow_class)$local_well, "EFPW02" = 1, "EFPW03" = 2)
colnames(sample_data(shallow_class))

shallow_class <- shallow_class %>% subset_taxa( Family!= "mitochondria" | is.na(Family) & Order!="Chloroplast Order" | is.na(Order)  & 
                                                  Class!="Chloroplast Order" | is.na(Class) & Order!="Chloroplast" | is.na(Order)& 
                                                  Class!="Chloroplast" | is.na(Class)& Genus!="Chloroplast Order" | is.na(Genus)  &
                                                  Phylum!="Chloroplast" | is.na(Phylum) & Class!="Unknown" | is.na(Class)& Class!="uncultured bacterium" | is.na(Class) &
                                                  Order!="Unknown" | is.na(Order)& Order!="uncultured bacterium" | is.na(Order) & Family!="Unknown" | is.na(Family)& 
                                                  Family!="uncultured bacterium" | is.na(Family) & Domain!="Kingdom_Eukaryota" | is.na(Domain) & 
                                                  Domain!="Eukaryota" | is.na(Domain) & Domain!="Unassigned" | is.na(Domain) &  Domain!="Kingdom_Unassigned" | is.na(Domain) &
                                                  Phylum!="Kingdom_Eukaryota" | is.na(Phylum) & Phylum!="Eukaryota" | is.na(Phylum) & Phylum!="Unassigned" | is.na(Phylum) &  
                                                  Phylum!="Kingdom_Unassigned" | is.na(Phylum) & Domain!="Eukaryota" | is.na(Domain) & Phylum!="6e1beec642e5732e97592afb9b0df55e" | is.na(Phylum) &
                                                  Phylum!="968463df5219e1bd2651952876319142"| is.na(Phylum))
shallow_class <- shallow_class %>% tax_fix(unknowns = c("uncultured bacterium", "uncultured prokaryote", "uncultured", "metagenome", "Order_Chloroplast", "Chloroplast Order", "Chloroplast",
                                                        "uncultured Chloroflexi bacterium", "uncultured Microgenomates group bacterium", "uncultured Acidimicrobidae bacterium", 
                                                        "uncultured Gemmatimonadetes bacterium", "uncultured forest soil bacterium", "uncultured Anaerolineaceae bacterium", 
                                                        "uncultured Desulfuromonadales bacterium", "uncultured archaeon", "uncultured Rhodospirillaceae bacterium",
                                                        "uncultured Rubrobacteria bacterium", "uncultured Acidobacteria bacterium", "uncultured Verrucomicrobia bacterium",
                                                        "uncultured actinobacterium", "uncultured Acidobacteriales bacterium", "Family_uncultured",
                                                        "uncultured deep-sea bacterium", "uncultured euryarchaeote", "uncultured microorganism", "uncultured delta proteobacterium", 
                                                        "uncultured soil bacterium", "uncultured Acetobacteraceae bacterium", "Order_uncultured", "uncultured Termite group 1 bacterium",
                                                        "uncultured proteobacterium", "uncultured beta proteobacterium", "uncultured Nitrosomonadaceae bacterium", 
                                                        "uncultured Syntrophobacterales bacterium", "uncultured Myxococcales bacterium", "Family_Unknown Family", "wastewater metagenome", 
                                                        "uncultured Actinomycetales bacterium", "groundwater metagenome", "uncultured alpha proteobacterium", "Family_uncultured bacterium ", 
                                                        "uncultured Bacteroidetes bacterium", "uncultured gamma proteobacterium", "uncultured Parcubacteria group bacterium", 
                                                        "uncultured planctomycete", "uncultured Acidobacterium sp.", "uncultured organism", "uncultured Holophaga sp.", "uncultured Firmicutes bacterium",
                                                        "uncultured Caldilineaceae bacterium", "uncultured Chlorobi bacterium", "uncultured Sphingobacteriales bacterium","KD4-96 CLass",
                                                        "uncultured Candidatus Dependentiae bacterium", "uncultured Conexibacteraceae bacterium", "uncultured Planctomycetales bacterium",
                                                        "Unknown Family Family", "uncultured bacterium  Family", "uncultured Sorangiineae bacterium", "unidentified","Unknown Family Family",
                                                        "uncultured bacterium  Family", "Kingdom_Eukaryota", "Kingdom_Unassigned", "Eukaryota Domain" ,"Eukaryota","Unassigned Domain","5911a83bc7642bdbb948309a94b1bf54" ,
                                                        "6e1beec642e5732e97592afb9b0df55e" ,  "968463df5219e1bd2651952876319142", "1d04ce2127c8fa8fa70da80ceb9adf63"))


#filter genera by abundance
#cpt_psq <- tax_filter(cpt_psq, min_prevalence = 1 / 100, min_sample_abundance = 1 / 100)
shallow_class <- microbiome::transform(shallow_class, 'clr')
shallow_class <- phyloseq::tax_glom(shallow_class, "Class", NArm = TRUE)
shallow_class <- tax_agg(shallow_class, "Class")
#phyloseq::tax_table(cpt_phy)
#otu_table(cpt_phy)
# select n taxa from the most abundant taxa
set.seed(0451)
taxa_shallow_class <- sample(tax_top(shallow_class, n = 50), size = 50)
taxa_shallow_class
taxa_shallow_class <-taxa_shallow_class[-15]
taxa_shallow_class
# NOTE: detection threshold set to 50 as HITchip example data seems to have background noise
ud <- 10
# chemistry correlations
jpeg("~/CPT_class_shallow_chem_cor.tiff",width=2200, height=2800, res=400, bg="white")
plot.new()
cor_heatmap(
  data = shallow_class , taxa = taxa_shallow_class, 
  vars = c("Sulfate", "Nitrate","Th","HREE","LREE","water_pH", "TOC","U" ),
  cor = "spearman",
  colors = heat_palette("Tropic", rev = FALSE, sym = TRUE),
  tax_anno = taxAnnotation(Prevalence = anno_tax_prev()),
  show_heatmap_legend = FALSE,
  var_anno = varAnnotation(
    ppm = anno_var_hist(size = grid::unit(20, "mm"))
    #Log10 = anno_var_box(function(x) log10(x + 1))
  )
)

getOption("device")
dev.set(which = dev.next())
dev.off()

# correlations and multiple annotations
jpeg("~/CPT_class_shallow_cor.tiff",width=2200, height=2800, res=400, bg="white")
plot.new()
cor_heatmap(
  data = shallow_class , taxa = taxa_shallow_class, 
  vars = c("silt","clay", "sand","large_sand","Pressure"), 
  cor = "spearman",
  colors = heat_palette("Tropic", rev = FALSE, sym = TRUE),
  tax_anno = taxAnnotation(Prevalence = anno_tax_prev()),
  show_heatmap_legend = FALSE,
  var_anno = varAnnotation(
    Percent = anno_var_hist(size = grid::unit(20, "mm"))
    #Log10 = anno_var_box(function(x) log10(x + 1))
  )
)

getOption("device")
dev.set(which = dev.next())
dev.off()

#rename col labels
colnames(sample_data(physeqCPT))
phyloseq::tax_table(physeqCPT)
phyloseq::otu_table(physeqCPT)
head(sample_data(shallow_cpt))
#shallow CPT correlations
shallow_cpt
shallow_ord <- shallow_cpt
sample_data(shallow_ord)$subsection <- dplyr::recode(sample_data(shallow_ord)$section, "A" = 1, "B" = 2, "C"= 3)
sample_data(shallow_ord)$well <- dplyr::recode(sample_data(shallow_ord)$local_well, "EFPW02" = 1, "EFPW03" = 2)
colnames(sample_data(shallow_ord))

shallow_ord <- shallow_ord %>% subset_taxa( Family!= "mitochondria" | is.na(Family) & Order!="Chloroplast Order" | is.na(Order)  & 
                                              Class!="Chloroplast Order" | is.na(Class) & Order!="Chloroplast" | is.na(Order)& 
                                              Class!="Chloroplast" | is.na(Class)& Genus!="Chloroplast Order" | is.na(Genus)  &
                                              Phylum!="Chloroplast" | is.na(Phylum) & Class!="Unknown" | is.na(Class)& Class!="uncultured bacterium" | is.na(Class) &
                                              Order!="Unknown" | is.na(Order)& Order!="uncultured bacterium" | is.na(Order) & Family!="Unknown" | is.na(Family)& 
                                              Family!="uncultured bacterium" | is.na(Family) & Domain!="Kingdom_Eukaryota" | is.na(Domain) & 
                                              Domain!="Eukaryota" | is.na(Domain) & Domain!="Unassigned" | is.na(Domain) &  Domain!="Kingdom_Unassigned" | is.na(Domain) &
                                              Phylum!="Kingdom_Eukaryota" | is.na(Phylum) & Phylum!="Eukaryota" | is.na(Phylum) & Phylum!="Unassigned" | is.na(Phylum) &  
                                              Phylum!="Kingdom_Unassigned" | is.na(Phylum) & Domain!="Eukaryota" | is.na(Domain) & Phylum!="6e1beec642e5732e97592afb9b0df55e" | is.na(Phylum) &
                                              Phylum!="968463df5219e1bd2651952876319142"| is.na(Phylum))
shallow_ord <- shallow_ord %>% tax_fix(unknowns = c("uncultured bacterium", "uncultured prokaryote", "uncultured", "metagenome", "Order_Chloroplast", "Chloroplast Order", "Chloroplast",
                                                    "uncultured Chloroflexi bacterium", "uncultured Microgenomates group bacterium", "uncultured Acidimicrobidae bacterium", 
                                                    "uncultured Gemmatimonadetes bacterium", "uncultured forest soil bacterium", "uncultured Anaerolineaceae bacterium", 
                                                    "uncultured Desulfuromonadales bacterium", "uncultured archaeon", "uncultured Rhodospirillaceae bacterium",
                                                    "uncultured Rubrobacteria bacterium", "uncultured Acidobacteria bacterium", "uncultured Verrucomicrobia bacterium",
                                                    "uncultured actinobacterium", "uncultured Acidobacteriales bacterium", "Family_uncultured",
                                                    "uncultured deep-sea bacterium", "uncultured euryarchaeote", "uncultured microorganism", "uncultured delta proteobacterium", 
                                                    "uncultured soil bacterium", "uncultured Acetobacteraceae bacterium", "Order_uncultured", "uncultured Termite group 1 bacterium",
                                                    "uncultured proteobacterium", "uncultured beta proteobacterium", "uncultured Nitrosomonadaceae bacterium", 
                                                    "uncultured Syntrophobacterales bacterium", "uncultured Myxococcales bacterium", "Family_Unknown Family", "wastewater metagenome", 
                                                    "uncultured Actinomycetales bacterium", "groundwater metagenome", "uncultured alpha proteobacterium", "Family_uncultured bacterium ", 
                                                    "uncultured Bacteroidetes bacterium", "uncultured gamma proteobacterium", "uncultured Parcubacteria group bacterium", 
                                                    "uncultured planctomycete", "uncultured Acidobacterium sp.", "uncultured organism", "uncultured Holophaga sp.", "uncultured Firmicutes bacterium",
                                                    "uncultured Caldilineaceae bacterium", "uncultured Chlorobi bacterium", "uncultured Sphingobacteriales bacterium","KD4-96 CLass",
                                                    "uncultured Candidatus Dependentiae bacterium", "uncultured Conexibacteraceae bacterium", "uncultured Planctomycetales bacterium",
                                                    "Unknown Family Family", "uncultured bacterium  Family", "uncultured Sorangiineae bacterium", "unidentified","Unknown Family Family",
                                                    "uncultured bacterium  Family", "Kingdom_Eukaryota", "Kingdom_Unassigned", "Eukaryota Domain" ,"Eukaryota","Unassigned Domain","5911a83bc7642bdbb948309a94b1bf54" ,
                                                    "6e1beec642e5732e97592afb9b0df55e" ,  "968463df5219e1bd2651952876319142", "1d04ce2127c8fa8fa70da80ceb9adf63"))


#filter genera by abundance
#cpt_psq <- tax_filter(cpt_psq, min_prevalence = 1 / 100, min_sample_abundance = 1 / 100)
shallow_ord <- microbiome::transform(shallow_ord, 'clr')
shallow_ord <- phyloseq::tax_glom(shallow_ord, "Order", NArm = TRUE)
shallow_ord <- tax_agg(shallow_ord, "Order")
#phyloseq::tax_table(cpt_phy)
#otu_table(cpt_phy)
# select n taxa from the most abundant taxa
set.seed(0451)
taxa_shallow_ord <- sample(tax_top(shallow_ord, n = 50), size = 50)
taxa_shallow_ord
taxa_shallow_ord <-taxa_shallow_ord[-15]
taxa_shallow_ord
# NOTE: detection threshold set to 50 as HITchip example data seems to have background noise
ud <- 10
# chemistry correlations
jpeg("~/CPT_class_shallow_chem_cor.tiff",width=2200, height=2800, res=400, bg="white")
plot.new()
cor_heatmap(
  data = shallow_ord , taxa = taxa_shallow_ord, 
  vars = c("Sulfate", "Nitrate","Th","HREE","LREE","water_pH", "TOC","U" ),
  cor = "spearman",
  colors = heat_palette("Tropic", rev = FALSE, sym = TRUE),
  tax_anno = taxAnnotation(Prevalence = anno_tax_prev()),
  show_heatmap_legend = FALSE,
  var_anno = varAnnotation(
    ppm = anno_var_hist(size = grid::unit(20, "mm"))
    #Log10 = anno_var_box(function(x) log10(x + 1))
  )
)

getOption("device")
dev.set(which = dev.next())
dev.off()

# correlations and multiple annotations
jpeg("~/CPT_class_shallow_cor.tiff",width=2200, height=2800, res=400, bg="white")
plot.new()
cor_heatmap(
  data = shallow_class , taxa = taxa_shallow_class, 
  vars = c("silt","clay", "sand","large_sand","Pressure","upper_pressure_difference", "lower_pressure_difference" ), 
  cor = "spearman",
  colors = heat_palette("Tropic", rev = FALSE, sym = TRUE),
  tax_anno = taxAnnotation(Prevalence = anno_tax_prev()),
  show_heatmap_legend = FALSE,
  var_anno = varAnnotation(
    Percent = anno_var_hist(size = grid::unit(20, "mm"))
    #Log10 = anno_var_box(function(x) log10(x + 1))
  )
)

getOption("device")
dev.set(which = dev.next())
dev.off()

#full community
#phylum-level
#class spearman plot
cpt_class <- cptphyseq
sample_data(cpt_class)$subsection <- dplyr::recode(sample_data(cpt_class)$section, "A" = 1, "B" = 2, "C"  = 3)
sample_data(cpt_class)$boundary <- dplyr::recode(sample_data(cpt_class)$section, "A" = 1, "B" = 0, "C"  = 1)
colnames(sample_data(cpt_class))

cpt_class <- cpt_class %>% subset_taxa( Family!= "mitochondria" | is.na(Family) & Order!="Chloroplast Order" | is.na(Order)  & 
                                          Class!="Chloroplast Order" | is.na(Class) & Order!="Chloroplast" | is.na(Order)& 
                                          Class!="Chloroplast" | is.na(Class)& Genus!="Chloroplast Order" | is.na(Genus)  &
                                          Phylum!="Chloroplast" | is.na(Phylum) & Class!="Unknown" | is.na(Class)& Class!="uncultured bacterium" | is.na(Class) &
                                          Order!="Unknown" | is.na(Order)& Order!="uncultured bacterium" | is.na(Order) & Family!="Unknown" | is.na(Family)& 
                                          Family!="uncultured bacterium" | is.na(Family) & Class!="Patescibacteria Phylum"| is.na(Class) )
cpt_class <- cpt_class%>% tax_fix(unknowns = c("uncultured bacterium", "uncultured prokaryote", "uncultured", "metagenome", "Order_Chloroplast", "Chloroplast Order", "Chloroplast",
                                               "uncultured Chloroflexi bacterium", "uncultured Microgenomates group bacterium", "uncultured Acidimicrobidae bacterium", 
                                               "uncultured Gemmatimonadetes bacterium", "uncultured forest soil bacterium", "uncultured Anaerolineaceae bacterium", 
                                               "uncultured Desulfuromonadales bacterium", "uncultured archaeon", "uncultured Rhodospirillaceae bacterium",
                                               "uncultured Rubrobacteria bacterium", "uncultured Acidobacteria bacterium", "uncultured Verrucomicrobia bacterium",
                                               "uncultured actinobacterium", "uncultured Acidobacteriales bacterium", "Family_uncultured",
                                               "uncultured deep-sea bacterium", "uncultured euryarchaeote", "uncultured microorganism", "uncultured delta proteobacterium", 
                                               "uncultured soil bacterium", "uncultured Acetobacteraceae bacterium", "Order_uncultured", "uncultured Termite group 1 bacterium",
                                               "uncultured proteobacterium", "uncultured beta proteobacterium", "uncultured Nitrosomonadaceae bacterium", 
                                               "uncultured Syntrophobacterales bacterium", "uncultured Myxococcales bacterium", "Family_Unknown Family", "wastewater metagenome", 
                                               "uncultured Actinomycetales bacterium", "groundwater metagenome", "uncultured alpha proteobacterium", "Family_uncultured bacterium ", 
                                               "uncultured Bacteroidetes bacterium", "uncultured gamma proteobacterium", "uncultured Parcubacteria group bacterium", 
                                               "uncultured planctomycete", "uncultured Acidobacterium sp.", "uncultured organism", "uncultured Holophaga sp.", "uncultured Firmicutes bacterium",
                                               "uncultured Caldilineaceae bacterium", "uncultured Chlorobi bacterium", "uncultured Sphingobacteriales bacterium","KD4-96 CLass",
                                               "uncultured Candidatus Dependentiae bacterium", "uncultured Conexibacteraceae bacterium", "uncultured Planctomycetales bacterium",
                                               "Unknown Family Family", "uncultured bacterium  Family", "uncultured Sorangiineae bacterium", "unidentified","Unknown Family Family",
                                               "uncultured bacterium  Family","uncultured Syntrophaceae bacterium", "uncultured Nitrospirae bacterium", "Unknown Family", 
                                               "uncultured Planctomycetaceae bacterium", "uncultured bacterium ", "Family XI","uncultured Rubrobacterales bacterium", "marine metagenome", 
                                               "uncultured Chlamydiales bacterium", "uncultured Acidobacteriaceae bacterium", "uncultured Rubrobacteraceae bacterium","Patescibacteria Phylum"))


#filter genera by abundance
#cpt_fam <- tax_filter(cpt_fam, min_prevalence = 1 / 100, min_sample_abundance = 1 / 100)
cpt_class <- microbiome::transform(cpt_class, 'clr')
#phyloseq::tax_table(cpt_fam)
cpt_class <- phyloseq::tax_glom(cpt_class, "Class", NArm = TRUE)
cpt_class <- tax_agg(cpt_class, "Class")
otu_table(cpt_class)
# randomly select 20 taxa from the 50 most abundant taxa
set.seed(0451)
taxa_class1 <- sample(tax_top(cpt_class, n = 50), size = 50)
taxa_class1
sample_data(cptphyseq)
# NOTE: detection threshold set to 50 as HITchip example data seems to have background noise
ud <- 10
# correlations and multiple annotations
jpeg("~/CPT_class_chem_cor.tiff",width=2800, height=2800, res=400, bg="white")
plot.new()
cor_heatmap(
  data = cpt_class , taxa = taxa_class1, 
  vars = c("Sulfate", "Nitrate","Th","HREE","LREE","water_pH", "TOC","U" ), 
  cor = "spearman",
  colors = heat_palette("Tropic", rev = FALSE, sym = TRUE),
  show_heatmap_legend = FALSE,
  tax_anno = taxAnnotation(Prevalence = anno_tax_prev()),
  var_anno = varAnnotation(
    ppm = anno_var_hist(size = grid::unit(20, "mm"))
    #Log10 = anno_var_box(function(x) log10(x + 1))
  )
)

getOption("device")
dev.set(which = dev.next())
dev.off()
colnames(sample_data(cpt_class))

jpeg("~/CPT_class_sediment_cor.tiff",width=2200, height=2800, res=400, bg="white")
plot.new()
cor_heatmap(
  data = cpt_class , taxa = taxa_class1, 
  vars = c("silt","clay", "sand","large_sand","Pressure" ), 
  cor = "spearman",
  show_heatmap_legend = FALSE,
  colors = heat_palette("Tropic", rev = FALSE, sym = TRUE),
  tax_anno = taxAnnotation(Prevalence = anno_tax_prev()),
  var_anno = varAnnotation(
    percent = anno_var_hist(size = grid::unit(20, "mm"))
    #Log10 = anno_var_box(function(x) log10(x + 1))
  )
)

getOption("device")
dev.set(which = dev.next())
dev.off()


#class spearman plot
cpt_class <- cptphyseq
sample_data(cpt_class)$subsection <- dplyr::recode(sample_data(cpt_class)$section, "A" = 1, "B" = 2, "C"  = 3)
sample_data(cpt_class)$boundary <- dplyr::recode(sample_data(cpt_class)$section, "A" = 1, "B" = 0, "C"  = 1)
colnames(sample_data(cpt_class))

cpt_class <- cpt_class %>% subset_taxa( Family!= "mitochondria" | is.na(Family) & Order!="Chloroplast Order" | is.na(Order)  & 
                                          Class!="Chloroplast Order" | is.na(Class) & Order!="Chloroplast" | is.na(Order)& 
                                          Class!="Chloroplast" | is.na(Class)& Genus!="Chloroplast Order" | is.na(Genus)  &
                                          Phylum!="Chloroplast" | is.na(Phylum) & Class!="Unknown" | is.na(Class)& Class!="uncultured bacterium" | is.na(Class) &
                                          Order!="Unknown" | is.na(Order)& Order!="uncultured bacterium" | is.na(Order) & Family!="Unknown" | is.na(Family)& 
                                          Family!="uncultured bacterium" | is.na(Family) & Class!="Patescibacteria Phylum"| is.na(Class) )
cpt_class <- cpt_class%>% tax_fix(unknowns = c("uncultured bacterium", "uncultured prokaryote", "uncultured", "metagenome", "Order_Chloroplast", "Chloroplast Order", "Chloroplast",
                                               "uncultured Chloroflexi bacterium", "uncultured Microgenomates group bacterium", "uncultured Acidimicrobidae bacterium", 
                                               "uncultured Gemmatimonadetes bacterium", "uncultured forest soil bacterium", "uncultured Anaerolineaceae bacterium", 
                                               "uncultured Desulfuromonadales bacterium", "uncultured archaeon", "uncultured Rhodospirillaceae bacterium",
                                               "uncultured Rubrobacteria bacterium", "uncultured Acidobacteria bacterium", "uncultured Verrucomicrobia bacterium",
                                               "uncultured actinobacterium", "uncultured Acidobacteriales bacterium", "Family_uncultured",
                                               "uncultured deep-sea bacterium", "uncultured euryarchaeote", "uncultured microorganism", "uncultured delta proteobacterium", 
                                               "uncultured soil bacterium", "uncultured Acetobacteraceae bacterium", "Order_uncultured", "uncultured Termite group 1 bacterium",
                                               "uncultured proteobacterium", "uncultured beta proteobacterium", "uncultured Nitrosomonadaceae bacterium", 
                                               "uncultured Syntrophobacterales bacterium", "uncultured Myxococcales bacterium", "Family_Unknown Family", "wastewater metagenome", 
                                               "uncultured Actinomycetales bacterium", "groundwater metagenome", "uncultured alpha proteobacterium", "Family_uncultured bacterium ", 
                                               "uncultured Bacteroidetes bacterium", "uncultured gamma proteobacterium", "uncultured Parcubacteria group bacterium", 
                                               "uncultured planctomycete", "uncultured Acidobacterium sp.", "uncultured organism", "uncultured Holophaga sp.", "uncultured Firmicutes bacterium",
                                               "uncultured Caldilineaceae bacterium", "uncultured Chlorobi bacterium", "uncultured Sphingobacteriales bacterium","KD4-96 CLass",
                                               "uncultured Candidatus Dependentiae bacterium", "uncultured Conexibacteraceae bacterium", "uncultured Planctomycetales bacterium",
                                               "Unknown Family Family", "uncultured bacterium  Family", "uncultured Sorangiineae bacterium", "unidentified","Unknown Family Family",
                                               "uncultured bacterium  Family","uncultured Syntrophaceae bacterium", "uncultured Nitrospirae bacterium", "Unknown Family", 
                                               "uncultured Planctomycetaceae bacterium", "uncultured bacterium ", "Family XI","uncultured Rubrobacterales bacterium", "marine metagenome", 
                                               "uncultured Chlamydiales bacterium", "uncultured Acidobacteriaceae bacterium", "uncultured Rubrobacteraceae bacterium","Patescibacteria Phylum"))


#filter genera by abundance
#cpt_fam <- tax_filter(cpt_fam, min_prevalence = 1 / 100, min_sample_abundance = 1 / 100)
cpt_class <- microbiome::transform(cpt_class, 'clr')
#phyloseq::tax_table(cpt_fam)
cpt_class <- phyloseq::tax_glom(cpt_class, "Class", NArm = TRUE)
cpt_class <- tax_agg(cpt_class, "Class")
otu_table(cpt_class)
# randomly select 20 taxa from the 50 most abundant taxa
set.seed(0451)
taxa_class1 <- sample(tax_top(cpt_class, n = 50), size = 50)
taxa_class1
taxa_class1 <-taxa_class1[-36]
sample_data(cptphyseq)
# NOTE: detection threshold set to 50 as HITchip example data seems to have background noise
ud <- 10
# correlations and multiple annotations
jpeg("~/CPT_class_chem_cor.tiff",width=2800, height=2800, res=400, bg="white")
plot.new()
cor_heatmap(
  data = cpt_class , taxa = taxa_class1, 
  vars = c("Sulfate", "Nitrate","Th","HREE","LREE","water_pH", "TOC","U" ), 
  cor = "spearman",
  colors = heat_palette("Tropic", rev = FALSE, sym = TRUE),
  show_heatmap_legend = FALSE,
  tax_anno = taxAnnotation(Prevalence = anno_tax_prev()),
  var_anno = varAnnotation(
    ppm = anno_var_hist(size = grid::unit(20, "mm"))
    #Log10 = anno_var_box(function(x) log10(x + 1))
  )
)

getOption("device")
dev.set(which = dev.next())
dev.off()
colnames(sample_data(cpt_class))

jpeg("~/CPT_class_sediment_cor.tiff",width=2200, height=2800, res=400, bg="white")
plot.new()
cor_heatmap(
  data = cpt_class , taxa = taxa_class1, 
  vars = c("silt","clay", "sand","large_sand","Pressure" ), 
  cor = "spearman",
  show_heatmap_legend = FALSE,
  colors = heat_palette("Tropic", rev = FALSE, sym = TRUE),
  tax_anno = taxAnnotation(Prevalence = anno_tax_prev()),
  var_anno = varAnnotation(
    percent = anno_var_hist(size = grid::unit(20, "mm"))
    #Log10 = anno_var_box(function(x) log10(x + 1))
  )
)

getOption("device")
dev.set(which = dev.next())
dev.off()

#Order-level correlation
cpt_psq <- cptphyseq
sample_data(cpt_psq)$subsection <- dplyr::recode(sample_data(cpt_psq)$section, "A" = 1, "B" = 2, "C"  = 3)
sample_data(cpt_psq)$boundary <- dplyr::recode(sample_data(cpt_psq)$section, "A" = 1, "B" = 0, "C"  = 1)
colnames(sample_data(cpt_psq))

cpt_psq <- cpt_psq %>% subset_taxa( Family!= "mitochondria" | is.na(Family) & Order!="Chloroplast Order" | is.na(Order)  & 
                                      Class!="Chloroplast Order" | is.na(Class) & Order!="Chloroplast" | is.na(Order)& 
                                      Class!="Chloroplast" | is.na(Class)& Genus!="Chloroplast Order" | is.na(Genus)  &
                                      Phylum!="Chloroplast" | is.na(Phylum) & Class!="Unknown" | is.na(Class)& Class!="uncultured bacterium" | is.na(Class) &
                                      Order!="Unknown" | is.na(Order)& Order!="uncultured bacterium" | is.na(Order) & Family!="Unknown" | is.na(Family)& Family!="uncultured bacterium" | is.na(Family) )
cpt_psq <- cpt_psq%>% tax_fix(unknowns = c("uncultured bacterium", "uncultured prokaryote", "uncultured", "metagenome", "Order_Chloroplast", "Chloroplast Order", "Chloroplast",
                                           "uncultured Chloroflexi bacterium", "uncultured Microgenomates group bacterium", "uncultured Acidimicrobidae bacterium", 
                                           "uncultured Gemmatimonadetes bacterium", "uncultured forest soil bacterium", "uncultured Anaerolineaceae bacterium", 
                                           "uncultured Desulfuromonadales bacterium", "uncultured archaeon", "uncultured Rhodospirillaceae bacterium",
                                           "uncultured Rubrobacteria bacterium", "uncultured Acidobacteria bacterium", "uncultured Verrucomicrobia bacterium",
                                           "uncultured actinobacterium", "uncultured Acidobacteriales bacterium", "Family_uncultured",
                                           "uncultured deep-sea bacterium", "uncultured euryarchaeote", "uncultured microorganism", "uncultured delta proteobacterium", 
                                           "uncultured soil bacterium", "uncultured Acetobacteraceae bacterium", "Order_uncultured", "uncultured Termite group 1 bacterium",
                                           "uncultured proteobacterium", "uncultured beta proteobacterium", "uncultured Nitrosomonadaceae bacterium", 
                                           "uncultured Syntrophobacterales bacterium", "uncultured Myxococcales bacterium", "Family_Unknown Family", "wastewater metagenome", 
                                           "uncultured Actinomycetales bacterium", "groundwater metagenome", "uncultured alpha proteobacterium", "Family_uncultured bacterium ", 
                                           "uncultured Bacteroidetes bacterium", "uncultured gamma proteobacterium", "uncultured Parcubacteria group bacterium", 
                                           "uncultured planctomycete", "uncultured Acidobacterium sp.", "uncultured organism", "uncultured Holophaga sp.", "uncultured Firmicutes bacterium",
                                           "uncultured Caldilineaceae bacterium", "uncultured Chlorobi bacterium", "uncultured Sphingobacteriales bacterium","KD4-96 CLass",
                                           "uncultured Candidatus Dependentiae bacterium", "uncultured Conexibacteraceae bacterium", "uncultured Planctomycetales bacterium",
                                           "Unknown Family Family", "uncultured bacterium  Family", "uncultured Sorangiineae bacterium", "unidentified","Unknown Family Family",
                                           "uncultured bacterium  Family","Phylum_FCPU426" ))


#filter genera by abundance
#cpt_psq <- tax_filter(cpt_psq, min_prevalence = 1 / 100, min_sample_abundance = 1 / 100)
cpt_psq <- microbiome::transform(cpt_psq, 'clr')
phyloseq::tax_table(cpt_psq)
cpt_psq <- phyloseq::tax_glom(cpt_psq, "Order", NArm = TRUE)
cpt_psq <- tax_agg(cpt_psq, "Order")
otu_table(cpt_psq)
# randomly select 20 taxa from the 50 most abundant taxa
set.seed(0451)
taxa_cpt <- sample(tax_top(cpt_psq, n = 50), size = 50)
taxa_cpt
# NOTE: detection threshold set to 50 as HITchip example data seems to have background noise
ud <- 10
# correlations with chemistry
jpeg("~/CPT_taxa_chem_order_cor.tiff",width=2200, height=2800, res=400, bg="white")
plot.new()
cor_heatmap(
  data = cpt_psq , taxa = taxa_cpt, 
  vars = c("Sulfate", "Nitrate","U","Th","HREE","LREE", "water_pH","TOC" ), 
  cor = "spearman",
  colors = heat_palette("Tropic", rev = FALSE, sym = TRUE),
  tax_anno = taxAnnotation(Prevalence = anno_tax_prev()),
  var_anno = varAnnotation(
    ppm = anno_var_hist(size = grid::unit(20, "mm"))
    #Log10 = anno_var_box(function(x) log10(x + 1))
  )
)

getOption("device")
dev.set(which = dev.next())
dev.off()

# correlations with sediment
jpeg("~/CPT_taxa_sed_order_cor.jpeg",width=2200, height=2800, res=400, bg="white")
plot.new()
cor_heatmap(
  data = cpt_psq , taxa = taxa_cpt, 
  vars = c("silt","clay", "sand","large_sand","Pressure" ), 
  cor = "spearman",
  colors = heat_palette("Tropic", rev = FALSE, sym = TRUE),
  tax_anno = taxAnnotation(Prevalence = anno_tax_prev()),
  var_anno = varAnnotation(
    percent = anno_var_hist(size = grid::unit(20, "mm"))
    #Log10 = anno_var_box(function(x) log10(x + 1))
  )
)

getOption("device")
dev.set(which = dev.next())
dev.off()

# correlations with sediment
jpeg("~/CPT_taxa_order_cor.tiff",width=2200, height=2800, res=400, bg="white")
plot.new()
cor_heatmap(
  data = cpt_psq , taxa = taxa_cpt, 
  vars = c("mid_depth_cm" ,"Longitude"), 
  cor = "spearman",
  colors = heat_palette("Tropic", rev = FALSE, sym = TRUE),
  tax_anno = taxAnnotation(Prevalence = anno_tax_prev())
  #var_anno = varAnnotation(
  #  meters = anno_var_hist(size = grid::unit(20, "mm"))
  #Log10 = anno_var_box(function(x) log10(x + 1))
)
)

getOption("device")
dev.set(which = dev.next())
dev.off()

#family spearman plot
cpt_fam <- cptphyseq
sample_data(cpt_fam)$subsection <- dplyr::recode(sample_data(cpt_fam)$section, "A" = 1, "B" = 2, "C"  = 3)
sample_data(cpt_fam)$boundary <- dplyr::recode(sample_data(cpt_fam)$section, "A" = 1, "B" = 0, "C"  = 1)
colnames(sample_data(cpt_fam))

cpt_fam <- cpt_fam %>% subset_taxa( Family!= "mitochondria" | is.na(Family) & Order!="Chloroplast Order" | is.na(Order)  & 
                                      Class!="Chloroplast Order" | is.na(Class) & Order!="Chloroplast" | is.na(Order)& 
                                      Class!="Chloroplast" | is.na(Class)& Genus!="Chloroplast Order" | is.na(Genus)  &
                                      Phylum!="Chloroplast" | is.na(Phylum) & Class!="Unknown" | is.na(Class)& Class!="uncultured bacterium" | is.na(Class) &
                                      Order!="Unknown" | is.na(Order)& Order!="uncultured bacterium" | is.na(Order) & Family!="Unknown" | is.na(Family)& Family!="uncultured bacterium" | is.na(Family) )
cpt_fam <- cpt_fam%>% tax_fix(unknowns = c("uncultured bacterium", "uncultured prokaryote", "uncultured", "metagenome", "Order_Chloroplast", "Chloroplast Order", "Chloroplast",
                                           "uncultured Chloroflexi bacterium", "uncultured Microgenomates group bacterium", "uncultured Acidimicrobidae bacterium", 
                                           "uncultured Gemmatimonadetes bacterium", "uncultured forest soil bacterium", "uncultured Anaerolineaceae bacterium", 
                                           "uncultured Desulfuromonadales bacterium", "uncultured archaeon", "uncultured Rhodospirillaceae bacterium",
                                           "uncultured Rubrobacteria bacterium", "uncultured Acidobacteria bacterium", "uncultured Verrucomicrobia bacterium",
                                           "uncultured actinobacterium", "uncultured Acidobacteriales bacterium", "Family_uncultured",
                                           "uncultured deep-sea bacterium", "uncultured euryarchaeote", "uncultured microorganism", "uncultured delta proteobacterium", 
                                           "uncultured soil bacterium", "uncultured Acetobacteraceae bacterium", "Order_uncultured", "uncultured Termite group 1 bacterium",
                                           "uncultured proteobacterium", "uncultured beta proteobacterium", "uncultured Nitrosomonadaceae bacterium", 
                                           "uncultured Syntrophobacterales bacterium", "uncultured Myxococcales bacterium", "Family_Unknown Family", "wastewater metagenome", 
                                           "uncultured Actinomycetales bacterium", "groundwater metagenome", "uncultured alpha proteobacterium", "Family_uncultured bacterium ", 
                                           "uncultured Bacteroidetes bacterium", "uncultured gamma proteobacterium", "uncultured Parcubacteria group bacterium", 
                                           "uncultured planctomycete", "uncultured Acidobacterium sp.", "uncultured organism", "uncultured Holophaga sp.", "uncultured Firmicutes bacterium",
                                           "uncultured Caldilineaceae bacterium", "uncultured Chlorobi bacterium", "uncultured Sphingobacteriales bacterium","KD4-96 CLass",
                                           "uncultured Candidatus Dependentiae bacterium", "uncultured Conexibacteraceae bacterium", "uncultured Planctomycetales bacterium",
                                           "Unknown Family Family", "uncultured bacterium  Family", "uncultured Sorangiineae bacterium", "unidentified","Unknown Family Family",
                                           "uncultured bacterium  Family","uncultured Syntrophaceae bacterium", "uncultured Nitrospirae bacterium", "Unknown Family", 
                                           "uncultured Planctomycetaceae bacterium", "uncultured bacterium ", "Family XI"))


#filter genera by abundance
#cpt_fam <- tax_filter(cpt_fam, min_prevalence = 1 / 100, min_sample_abundance = 1 / 100)
cpt_fam <- microbiome::transform(cpt_fam, 'clr')
#phyloseq::tax_table(cpt_fam)
cpt_fam <- phyloseq::tax_glom(cpt_fam, "Family", NArm = TRUE)
cpt_fam <- tax_agg(cpt_fam, "Family")
otu_table(cpt_fam)
# randomly select 20 taxa from the 50 most abundant taxa
set.seed(0451)
taxa_fam <- sample(tax_top(cpt_fam, n = 50), size = 50)
taxa_fam
# NOTE: detection threshold set to 50 as HITchip example data seems to have background noise
ud <- 10
# correlations and multiple annotations
jpeg("~/CPT_family_chem_cor.jpeg",width=2200, height=2800, res=400, bg="white")
plot.new()
cor_heatmap(
  data = cpt_fam , taxa = taxa_fam, 
  vars = c("Sulfate", "Nitrate","U","Th","HREE","LREE", "water_pH","TOC" ), 
  cor = "spearman",
  colors = heat_palette("Tropic", rev = FALSE, sym = TRUE),
  tax_anno = taxAnnotation(Prevalence = anno_tax_prev()),
  var_anno = varAnnotation(
    level = anno_var_hist(size = grid::unit(20, "mm"))
    #Log10 = anno_var_box(function(x) log10(x + 1))
  )
)

getOption("device")
dev.set(which = dev.next())
dev.off()

jpeg("~/CPT_family_sediment_cor.jpeg",width=2200, height=2800, res=400, bg="white")
plot.new()
cor_heatmap(
  data = cpt_fam , taxa = taxa_fam, 
  vars = c("Sulfate", "Nitrate","U_238","HREE","LREE", "water_pH","mid_depth_cm" ,"fine", "coarse", "TOC","saturation", "Longitude" ), 
  cor = "spearman",
  colors = heat_palette("Tropic", rev = FALSE, sym = TRUE),
  tax_anno = taxAnnotation(Prevalence = anno_tax_prev()),
  var_anno = varAnnotation(
    level = anno_var_hist(size = grid::unit(20, "mm"))
    #Log10 = anno_var_box(function(x) log10(x + 1))
  )
)

getOption("device")
dev.set(which = dev.next())
dev.off()

#genus spearman plot
cpt_genus <- cptphyseq
sample_data(cpt_genus)$subsection <- dplyr::recode(sample_data(cpt_genus)$section, "A" = 1, "B" = 2, "C"  = 3)
sample_data(cpt_genus)$boundary <- dplyr::recode(sample_data(cpt_genus)$section, "A" = 1, "B" = 0, "C"  = 1)
colnames(sample_data(cpt_genus))

cpt_genus <- cpt_genus %>% subset_taxa( Family!= "mitochondria" | is.na(Family) & Order!="Chloroplast Order" | is.na(Order)  & 
                                          Class!="Chloroplast Order" | is.na(Class) & Order!="Chloroplast" | is.na(Order)& 
                                          Class!="Chloroplast" | is.na(Class)& Genus!="Chloroplast Order" | is.na(Genus)  &
                                          Phylum!="Chloroplast" | is.na(Phylum) & Class!="Unknown" | is.na(Class)& Class!="uncultured bacterium" | is.na(Class) &
                                          Order!="Unknown" | is.na(Order)& Order!="uncultured bacterium" | is.na(Order) & Family!="Unknown" | is.na(Family)& Family!="uncultured bacterium" | is.na(Family) )
cpt_genus <- cpt_genus%>% tax_fix(unknowns = c("uncultured bacterium", "uncultured prokaryote", "uncultured", "metagenome", "Order_Chloroplast", "Chloroplast Order", "Chloroplast",
                                               "uncultured Chloroflexi bacterium", "uncultured Microgenomates group bacterium", "uncultured Acidimicrobidae bacterium", 
                                               "uncultured Gemmatimonadetes bacterium", "uncultured forest soil bacterium", "uncultured Anaerolineaceae bacterium", 
                                               "uncultured Desulfuromonadales bacterium", "uncultured archaeon", "uncultured Rhodospirillaceae bacterium",
                                               "uncultured Rubrobacteria bacterium", "uncultured Acidobacteria bacterium", "uncultured Verrucomicrobia bacterium",
                                               "uncultured actinobacterium", "uncultured Acidobacteriales bacterium", "Family_uncultured",
                                               "uncultured deep-sea bacterium", "uncultured euryarchaeote", "uncultured microorganism", "uncultured delta proteobacterium", 
                                               "uncultured soil bacterium", "uncultured Acetobacteraceae bacterium", "Order_uncultured", "uncultured Termite group 1 bacterium",
                                               "uncultured proteobacterium", "uncultured beta proteobacterium", "uncultured Nitrosomonadaceae bacterium", 
                                               "uncultured Syntrophobacterales bacterium", "uncultured Myxococcales bacterium", "Family_Unknown Family", "wastewater metagenome", 
                                               "uncultured Actinomycetales bacterium", "groundwater metagenome", "uncultured alpha proteobacterium", "Family_uncultured bacterium ", 
                                               "uncultured Bacteroidetes bacterium", "uncultured gamma proteobacterium", "uncultured Parcubacteria group bacterium", 
                                               "uncultured planctomycete", "uncultured Acidobacterium sp.", "uncultured organism", "uncultured Holophaga sp.", "uncultured Firmicutes bacterium",
                                               "uncultured Caldilineaceae bacterium", "uncultured Chlorobi bacterium", "uncultured Sphingobacteriales bacterium","KD4-96 CLass",
                                               "uncultured Candidatus Dependentiae bacterium", "uncultured Conexibacteraceae bacterium", "uncultured Planctomycetales bacterium",
                                               "Unknown Family Family", "uncultured bacterium  Family", "uncultured Sorangiineae bacterium", "unidentified","Unknown Family Family",
                                               "uncultured bacterium  Family","uncultured Syntrophaceae bacterium", "uncultured Nitrospirae bacterium", "Unknown Family", 
                                               "uncultured Planctomycetaceae bacterium", "uncultured bacterium ", "Family XI","uncultured Rubrobacterales bacterium", "marine metagenome", 
                                               "uncultured Chlamydiales bacterium", "uncultured Acidobacteriaceae bacterium", "uncultured Rubrobacteraceae bacterium"))


#filter genera by abundance
#cpt_fam <- tax_filter(cpt_fam, min_prevalence = 1 / 100, min_sample_abundance = 1 / 100)
cpt_genus <- microbiome::transform(cpt_genus, 'clr')
#phyloseq::tax_table(cpt_fam)
cpt_genus <- phyloseq::tax_glom(cpt_genus, "Genus", NArm = TRUE)
cpt_genus <- tax_agg(cpt_genus, "Genus")
otu_table(cpt_genus)
# randomly select 20 taxa from the 50 most abundant taxa
set.seed(0451)
taxa_genus <- sample(tax_top(cpt_genus, n = 50), size = 50)
taxa_genus
taxa_genus <- taxa_genus[-17]
# NOTE: detection threshold set to 50 as HITchip example data seems to have background noise
ud <- 10
# correlations and multiple annotations
jpeg("~/CPT_genus_chem_cor.tiff",width=2800, height=2800, res=400, bg="white")
plot.new()
cor_heatmap(
  data = cpt_genus , taxa = taxa_genus, 
  vars = c("Sulfate", "Nitrate","Th","HREE","LREE","water_pH", "TOC","U"  ), 
  cor = "spearman",
  colors = heat_palette("Tropic", rev = FALSE, sym = TRUE),
  tax_anno = taxAnnotation(Prevalence = anno_tax_prev()),
  show_heatmap_legend = FALSE,
  var_anno = varAnnotation(
    ppm = anno_var_hist(size = grid::unit(20, "mm"))
    #Log10 = anno_var_box(function(x) log10(x + 1))
  )
)

getOption("device")
dev.set(which = dev.next())
dev.off()

jpeg("~/CPT_genus_sediment_cor.tiff",width=2200, height=2800, res=400, bg="white")
plot.new()
cor_heatmap(
  data = cpt_genus , taxa = taxa_genus, 
  vars = c("silt","clay", "sand","large_sand","Pressure","Overlying.Pressure", "Deeper.Pressure"), 
  cor = "spearman",
  colors = heat_palette("Tropic", rev = FALSE, sym = TRUE),
  tax_anno = taxAnnotation(Prevalence = anno_tax_prev()),
  var_anno = varAnnotation(
    percent = anno_var_hist(size = grid::unit(20, "mm"))
    #Log10 = anno_var_box(function(x) log10(x + 1))
  )
)

getOption("device")
dev.set(which = dev.next())
dev.off()

#end taxa correlations with geochem
#####

#ordination

#ensure mitochondria, eukaryotes, chloroplasts are removed
#####
mergedCPT
mergedCPT <- mergedCPT %>% subset_taxa( Family!= "mitochondria" | is.na(Family) & Order!="Chloroplast Order" | is.na(Order)  & 
                                          Class!="Chloroplast Order" | is.na(Class) & Order!="Chloroplast" | is.na(Order)& 
                                          Class!="Chloroplast" | is.na(Class)& Genus!="Chloroplast Order" | is.na(Genus)  &
                                          Phylum!="Chloroplast" | is.na(Phylum) & Class!="Unknown" | is.na(Class)& Class!="uncultured bacterium" | is.na(Class) & Domain!="Unknown" | is.na(Domain)&
                                          Order!="Unknown" | is.na(Order)& Order!="uncultured bacterium" | is.na(Order) & Family!="Unknown" | is.na(Family)& Family!="uncultured bacterium" | is.na(Family)
                                        & Domain!="Eukaryote" | is.na(Domain) & Domain!="Eukaryota" | is.na(Domain)& Phylum!="Eukaryote" | is.na(Phylum)& Phylum!="Eukaryota" | is.na(Phylum))
mergedCPT <- mergedCPT%>% tax_fix(unknowns = c("uncultured bacterium", "uncultured prokaryote", "uncultured", "metagenome", "Order_Chloroplast", "Chloroplast Order", "Chloroplast",
                                               "unidentified","Unknown Family Family", "uncultured bacterium  Family", "Unknown Family", 
                                               "uncultured Planctomycetaceae bacterium", "Family XI","uncultured Rubrobacterales bacterium", "marine metagenome"))



#ordination with merged_p2017
#####
#CCA
merged_p2017_CCA  <- ordinate(merged_p2017~Days*well_type_factor, "CCA") #ordinate by CCA
merged_p2017_CCA
plot_scree(merged_p2017_CCA , "Scree Plot for Correspondence Analysis") #scree plot

# unifrac distance
set.seed(0451)
mergedCPT_unifrac = UniFrac(tax_glom(mergedCPT))
merged_CPT_otu <- as.matrix(otu_table(tax_glom(mergedCPT)))
View(merged_CPT_otu)
merged_p2017_genus_unifrac 

#bray curtis
merged_CPT_bray <- vegdist(wisconsin(sqrt(merged_CPT_otu)), method = "bray")
merged_CPT_bray

#NMDS of UNIFRAC
set.seed(0451)
cpt.mds <- metaMDS(mergedCPT_unifrac, zerodist=ignore, try=999 ,trymax = 1000)
p2017_merge_genus.mds; plot(p2017_merge_genus.mds)
ASV_tree <- hclust(merged_p2017_genus_unifrac, method = "average");plot(ASV_tree)
grp <- cutree(ASV_tree, k = 3)
#nmds.scores
nmds.scores <- as.data.frame(scores(p2017_merge_genus.mds, display = "sites"))
nmds.scores$grp <- grp
nrow(merged_p2017_numeric); nrow(merged_p2017_genus_unifrac)
merged_p2017_numeric
#end of ordination
#####

#merged_p2017 NMDS
#####
#ordinate with metadaata
env_ord <- envfit(p2017_merge_genus.mds, merged_p2017_numeric, na.rm = TRUE, permu = 1000)
env_ord
env_vector <- as.data.frame(scores(env_ord, display = "vectors"))
env_vector
env_vector <- env_vector[2:15,];env_vector<- env_vector[-3,]#keep only significant ordinations
env_vector
env_vector <- cbind(env_vector, Species = rownames(env_vector)); env_vector

#extract mds scores
p2017_merged_scrs <- cbind(merged_p2017_sample_data, data.frame(MDS1 = p2017_merge_genus.mds$points[,1], MDS2 = p2017_merge_genus.mds$points[,2])) 
p2017_merged_scrs$grp <- grp
env_vector
p2017_merged_scrs <- as.data.frame(p2017_merged_scrs)
head(p2017_merged_scrs)
#extract centeroid
p2017_merged_cent <- cbind(merged_p2017_sample_data, data.frame(MDS1 = p2017_merge_genus.mds$points[,1], MDS2 = p2017_merge_genus.mds$points[,2])) %>% 
  aggregate(cbind(MDS1, MDS2) ~ Days +  well_type + sulfate_level + nitrate_level,data = ., FUN = mean) 
#plot mds (not betadispersion and NOT CAP)
p2017_merged_scrs
p2017_merged_scrs$Days <- factor(p2017_merged_scrs$Days, levels = c("-6" , "1" ,  "8" ,  "15" , "22" , "50" , "78" , "106" ,"134"))

ggplot() +
  #facet_wrap(~diet) + 
  geom_hline(yintercept = 0, color = "gray50", size = 0.5, alpha=0.6, linetype="dashed")+
  geom_vline(xintercept = 0, color ="gray50", size=0.5, alpha=0.6, linetype="dashed")+
  #stat_ellipse(data=p2017_merged_scrs, aes(x = MDS1, y = MDS2, fill = as.factor(phase_cat), color = as.factor(phase_cat)),  linetype = 2, geom = "polygon", alpha = 0.6) +# add ellipse
  #geom_point(data = p2017_merged_cent, size = 5) +            # centroids
  geom_point(data = p2017_merged_scrs, aes(x = MDS1, y = MDS2,   shape=well_type),size=5.5, color = "black")+
  geom_point(data = p2017_merged_scrs, aes(x = MDS1, y = MDS2,  color = as.factor(Days), shape=well_type),size=4) +                                              # sample scores
  labs(shape = "Well Type", color = "Days") + 
  coord_fixed() + #must have a fixed centroid
  geom_segment(data = env_vector, # vector arrows
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "gray70") +
  geom_text(data = env_vector, aes(x = NMDS1, y = NMDS2, label = Species), color = "black", # vector labels
            size = 3) +# sample scores
  scale_colour_brewer(palette = "YlGnBu") + 
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = "bottom",
        axis.ticks = element_line(color = "black"),
        #legend.position = c(0.46, 0.92),
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        panel.border = element_rect(colour = "black", size=0.8, fill = "transparent"))                                           # same axis scaling
ggsave('~/Dropbox/Andrews MacBook (Hazen MacBook Copy)/Documents/Research/EVO-UMB/2017_unifrac_nmds_ordination.png', bg = "transparent", width = 40, height = 30, units = "cm")
#end NMDS
#####

#numeric dataframe for correlations and heatmap
#####
cpt_metadata
colnames(cpt_metadata)
cpt_metadata.condensed = cpt_metadata[seq(1, nrow(cpt_metadata), 2), ]#select every second row (since these are duplicate samples)
rownames(cpt_metadata.condensed) <- cpt_metadata.condensed[,1] #set row names as sediment short ID
cpt_metadata.condensed[,1] <-NULL #replace pre-existing rownames
View(cpt_metadata.condensed)
colnames(cpt_metadata.condensed)
cpt_numeric <- dplyr::select(cpt_metadata.condensed, c(8,17:22, 33,34,37,38,39,40,41,42,46,47,56,61, 67:79, 85, 87,92))
head(cpt_numeric)
# remove NA
cpt_numeric[is.na(cpt_numeric)] = 0
colnames(cpt_numeric) <- c( "depth (cm)" ,"coarse sand",   "medium/fine sand",         "silt",       "clay"     ,  ">0.05 mm"     ,  "<0.05 mm"    ,     "Nitrate"   ,   "Sulfate"    , 
                            "Na"           ,"Mg"        ,   "Al"  ,         "P"     ,       "K"   ,         "Ca"  ,         "Mn" ,          "Fe"  ,         "Rb" ,         
                            "Cd"           ,"La"         ,  "Pr"   ,        "Ce"     ,      "Nd"   ,        "Sm"   ,        "Eu"  ,         "Gd"   ,        "Dy"  ,        
                            "Ho"           ,"Er"          , "Tm"    ,       "Yb"      ,     "Lu"    ,       "Th"    ,       "U_238",        "AODC"   )
colnames(cpt_numeric)
head(cpt_numeric)
#center/scale the data
scale_numeric <- scale(cpt_numeric , scale = TRUE) 
#log-transform
log_numeric <- log10(1+cpt_numeric) 
log_numeric[is.na(log_numeric)] = 0
#spearman correlation
cc = cor(scale_numeric, method = "spearman")
cc
corrplot(cc)

#make formatted correlation plot
jpeg("Dropbox/Andrews MacBook (Hazen MacBook Copy)/Documents/Research/CPT/cpt_sed_corrplot.jpeg",width=2400, height=2200, res=400, bg="white")
plot.new()
corrplot(cc, method="circle",#addCoef.col = 'black',
         tl.col = "black", col = COL2('PuOr'),type="lower",
         #order = "hclust", hclust.method = "average", addrect = 5, tl.cex = 0.7,
         diag=TRUE,tl.cex = 0.7)
dev.off()
install.packages("dendsort")
library(dendsort)
pheatmap(log_numeric)
col.pal <- brewer.pal(9,"YlOrBr")
cluster_row <- hclust(dist(log_numeric))
dd <-dendsort::dendsort(cluster_row)

jpeg("~/cpt_sed_heatmap.jpeg",width=2800, height=2000, res=400, bg="white")
plot.new()
pheatmap(
  mat               = scale_numeric,
  color             = col.pal,
  border_color      = 'gray80',
  show_colnames     = TRUE,
  cutree_rows = 4,
  treeheight_row = 20,
  treeheight_col = 20,
  angle_col=90,
  show_rownames     = TRUE,
  drop_levels       = TRUE,
  fontsize          = 10
)

getOption("device")
dev.set(which = dev.next())
dev.off()
#end plot correlation and heatmap
#####



#NMDS distance matrix of geochemistry
#####
sso_sed_df <- subset(sso_metadata, medium == "sediment")
colnames(sso_metadata)
sso_numeric_df <- dplyr::select(sso_metadata, c(7:72))#select only numeric data
sso_sed_numeric <- dplyr::select(sso_sed_df, c(10,12,13,14, 18,19,22:72))
sso_numeric_df
sso_sed_numeric
sso_numeric_df[is.na(sso_numeric_df)] = 0
sso_sed_numeric[is.na(sso_sed_numeric)] = 0

merged_p2017_genus_bray <- vegdist(wisconsin(sqrt(merged_p2017_genus)), method = "bray")
merged_p2017_genus_bray

#NMDS of geochem
set.seed(0451)
sso_metadata.mds <- metaMDS(sso_numeric_df, zerodist="ignore", distance = "euclidean", try=999 ,trymax = 1000)
sso_metadata.mds; plot(sso_metadata.mds)
sso_nmds_merged_scrs <- cbind(sso_metadata, data.frame(MDS1 = sso_metadata.mds$points[,1], MDS2 = sso_metadata.mds$points[,2])) 
sso_nmds_merged_scrs <- as.data.frame(sso_nmds_merged_scrs)

goodness(sso_metadata.mds)
sso.euc <- vegdist(wisconsin(sqrt(sso_numeric_df)), method = "euclidean")
i = lower.tri(as.matrix(sso.euc))

hist(as.matrix(sso.euc)[i], xlab='ecological distance', main='Euclidean Dissimilarity (water and sediment)', col = "gray")
stressplot(sso_metadata.mds)
cpt_NMDS <-ggplot() +
  #facet_wrap(~diet) + 
  geom_hline(yintercept = 0, color = "gray50", size = 0.5, alpha=0.6, linetype="dashed")+
  geom_vline(xintercept = 0, color ="gray50", size=0.5, alpha=0.6, linetype="dashed")+
  #stat_ellipse(data=sso_nmds_merged_scrs, aes(x = MDS1, y = MDS2, fill = local_well, color = local_well),  linetype = 2, geom = "polygon", alpha = 0.45) +# add ellipse
  #geom_point(data = p2017_merged_cent, size = 5) +            # centroids
  geom_point(data = sso_nmds_merged_scrs, aes(x = MDS1, y = MDS2, shape = medium),size=5.5, color = "gray20")+
  geom_point(data = sso_nmds_merged_scrs, aes(x = MDS1, y = MDS2,  color = local_well, shape = medium),size=4) +                                              # sample scores
  labs( color = "Depth (cm)", fill= "Depth (cm)") + 
  coord_fixed() + #must have a fixed centroid
  #geom_segment(data = env_vector, # vector arrows
  #             aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
  #             arrow = arrow(length = unit(0.25, "cm")), colour = "gray70") +
  #geom_text(data = sso_nmds_merged_scrs, aes(x = MDS1, y = MDS2, label = sample), color = "black", # vector labels
  #          size = 3) +# sample scores
  geom_text(aes(x=-0.24, label="Stress: Type 1 weak ties\n 999 Permutations\nStress: 0.03", y=-0.042), colour="gray30", angle=0, text=element_text(size=7)) +
  scale_fill_manual(values = color_list_cpt) +
  scale_color_manual(values = color_list_cpt) +
  #scale_colour_brewer(palette = "YlGnBu") + 
  #scale_fill_brewer(palette= "YlGnBu")+
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="vertical",
        legend.position = "right",
        axis.ticks = element_line(color = "black"),
        #legend.position = c(0.46, 0.92),
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        panel.border = element_rect(colour = "black", size=0.8, fill = "transparent"))     # same axis scaling
cpt_NMDS
ggsave('~/CPT_water_sediment_geochem_NMDS.jpeg', bg = "transparent", width = 40, height = 30, units = "cm")
#end NMDS

#sed NMDS
#subset data
sso_sed_df <- subset(sso_metadata, medium == "sediment")
colnames(sso_sed_df)
sso_sed_numeric <- dplyr::select(sso_sed_df, c(10,12,13,14, 18,19,22:74))
sso_sed_numeric
sso_sed_numeric[is.na(sso_sed_numeric)] = 0
#make NMDS
sso_sed.mds <- metaMDS(sso_sed_numeric, zerodist="ignore",distance = "euclidean", try=999 ,trymax = 1000)
sso_sed.mds; plot(sso_sed.mds)
goodness(sso_sed.mds)
stressplot(sso_sed.mds)
sso.sed.euc <- vegdist(wisconsin(sqrt(sso_sed_numeric)), method = "euclidean")
i = lower.tri(as.matrix(sso.sed.euc))
hist(as.matrix(sso.sed.euc)[i], xlab='ecological distance', main='euclideanDissimilarity (sediment)', col = "gray")
cpt_sed_tree <- hclust(sso.sed.euc, method = "average");plot(cpt_sed_tree)
grp <- cutree(cpt_sed_tree, k = 3)
dend <- as.dendrogram(cpt_sed_tree)
sed_dend<- as.dendrogram(cpt_sed_tree)
install.packages("dendextend")
library(dendextend)
sed_dend  <- color_branches(sed_dend, k = 3)
sed_dend <- color_labels(sed_dend, k = 3)
plot(sed_dend)
fill_col <- c("#7fcdbb", "orange", "goldenrod")
jpeg("~/sed_dend_average_linkage_method.jpeg",width=2700, height=2700, res=400, bg="white")

sed_dend %>%
  set("labels_col", value = c("#7fcdbb", "orange", "goldenrod"), k=3) %>%
  set("branches_k_color", value = c("#7fcdbb", "orange", "goldenrod"), k = 3) %>%
  plot(horiz=TRUE, axes=TRUE) #dendrogram

getOption("device")
dev.set(which = dev.next())
dev.off()
grp <- cutree(ASV_tree, k = 3)
#extract mds scores
sso_sed_nmds_merged_scrs <- cbind(sso_sed_df, data.frame(MDS1 = sso_sed.mds$points[,1], MDS2 = sso_sed.mds$points[,2])) 
sso_sed_nmds_merged_scrs$grp <- grp
sso_sed_nmds_merged_scrs <- as.data.frame(sso_sed_nmds_merged_scrs)
sso_sed_nmds_merged_scrs$depth_range <-as.factor(case_when((sso_sed_nmds_merged_scrs$mid_depth_cm > 250) & (sso_sed_nmds_merged_scrs$mid_depth_cm <= 450) ~ "201-450",
                                                           (sso_sed_nmds_merged_scrs$mid_depth_cm > 450) & (sso_sed_nmds_merged_scrs$mid_depth_cm <= 650) ~ "451-650",
                                                           (sso_sed_nmds_merged_scrs$mid_depth_cm > 650) & (sso_sed_nmds_merged_scrs$mid_depth_cm <= 850) ~ "601-850",
                                                           (sso_sed_nmds_merged_scrs$mid_depth_cm > 850) & (sso_sed_nmds_merged_scrs$mid_depth_cm <= 1000) ~ "851-1000"))

View(sso_sed_nmds_merged_scrs)
cpt_sed_NMDS <-ggplot() +
  #facet_wrap(~diet) + 
  geom_hline(yintercept = 0, color = "gray50", size = 0.5, alpha=0.6, linetype="dashed")+
  geom_vline(xintercept = 0, color ="gray50", size=0.5, alpha=0.6, linetype="dashed")+
  #stat_ellipse(data=sso_sed_nmds_merged_scrs, aes(x = MDS1, y = MDS2, color = grp),  linetype = 2, geom = "polygon", alpha = 0.45) +# add ellipse
  #geom_point(data = p2017_merged_cent, size = 5) +            # centroids
  geom_point(data = sso_sed_nmds_merged_scrs, aes(x = MDS1, y = MDS2),size=5.5, color = "gray80")+
  geom_point(data = sso_sed_nmds_merged_scrs, aes(x = MDS1, y = MDS2,  color = depth_range),size=4) +                                              # sample scores
  labs( color = "Depth (cm)", fill= "Depth (cm)") + 
  coord_fixed() + #must have a fixed centroid
  #geom_segment(data = env_vector, # vector arrows
  #             aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
  #             arrow = arrow(length = unit(0.25, "cm")), colour = "gray70") +
  geom_text(data = sso_sed_nmds_merged_scrs, aes(x = MDS1, y = MDS2, label = sample), color = "gray20", # vector labels
            size = 3) +# sample scores
  geom_text(aes(x=-0.02, label="Stress: Type 1 weak ties\n 999 Permutations\nStress: 0.13", y=0.047), colour="gray30", angle=0, text=element_text(size=7)) +
  #scale_fill_manual(values = fill_col) +
  #scale_color_manual(values = fill_col) +
  scale_colour_brewer(palette = "YlGnBu") + 
  #scale_fill_brewer(palette= "YlGnBu")+
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="vertical",
        legend.position = "right",
        axis.ticks = element_line(color = "black"),
        #legend.position = c(0.46, 0.92),
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        panel.border = element_rect(colour = "black", size=0.8, fill = "transparent"))     # same axis scaling
cpt_sed_NMDS
ggsave('~/CPT_sediment_geochem_NMDS.jpeg', bg = "transparent", width = 40, height = 30, units = "cm")
#end NMDS

#water NMDS
#subset data
sso_water_df <- subset(sso_metadata, medium == "water")
colnames(sso_water_df)
sso_water_numeric <- dplyr::select(sso_water_df, c(7:72))
sso_water_numeric
sso_water_numeric[is.na(sso_water_numeric)] = 0
#make NMDS
sso_water.mds <- metaMDS(sso_water_numeric, zerodist="ignore",distance=="euclidean", try=999 ,trymax = 1000)
sso_water.mds; plot(sso_water.mds)
goodness(sso_water.mds)
stressplot(sso_water.mds)
sso.water.bray <- vegdist(wisconsin(sqrt(sso_water_numeric)), method = "bray")
i = lower.tri(as.matrix(sso.water.bray))
hist(as.matrix(sso.water.bray)[i], xlab='ecological distance', main='Bray Curtis Dissimilarity (water)', col = "gray")
#extract mds scores
sso_water_merged_scrs <- cbind(sso_water_df, data.frame(MDS1 = sso_water.mds$points[,1], MDS2 = sso_water.mds$points[,2])) 

sso_water_merged_scrs <- as.data.frame(sso_water_merged_scrs)
head(sso_water_merged_scrs)


sso_metadata_tree <- hclust(sso_numeric_df, method = "average");plot(sso_metadata_tree)
grp <- cutree(sso_metadata_tree, k = 3)
#nmds.scores
nmds.scores <- as.data.frame(scores(sso_metadata.mds, display = "sites"))

#end of ordination
#####

#plot NMDS
#####
#extract mds scores
sso_nmds_merged_scrs <- cbind(sso_metadata, data.frame(MDS1 = sso_metadata.mds$points[,1], MDS2 = sso_metadata.mds$points[,2])) 

sso_nmds_merged_scrs <- as.data.frame(sso_nmds_merged_scrs)
head(sso_nmds_merged_scrs)
#extract centeroid
p2017_merged_cent <- cbind(merged_p2017_sample_data, data.frame(MDS1 = p2017_merge_genus.mds$points[,1], MDS2 = p2017_merge_genus.mds$points[,2])) %>% 
  aggregate(cbind(MDS1, MDS2) ~ Days +  well_type + sulfate_level + nitrate_level,data = ., FUN = mean) 
#plot mds (not betadispersion and NOT CAP)
sso_nmds_merged_scrs
cpt_NMDS <-ggplot() +
  #facet_wrap(~diet) + 
  geom_hline(yintercept = 0, color = "gray50", size = 0.5, alpha=0.6, linetype="dashed")+
  geom_vline(xintercept = 0, color ="gray50", size=0.5, alpha=0.6, linetype="dashed")+
  #stat_ellipse(data=sso_nmds_merged_scrs, aes(x = MDS1, y = MDS2, fill = local_well, color = local_well),  linetype = 2, geom = "polygon", alpha = 0.45) +# add ellipse
  #geom_point(data = p2017_merged_cent, size = 5) +            # centroids
  geom_point(data = sso_nmds_merged_scrs, aes(x = MDS1, y = MDS2, shape = medium),size=5.5, color = "gray20")+
  geom_point(data = sso_nmds_merged_scrs, aes(x = MDS1, y = MDS2,  color = local_well, shape = medium),size=4) +                                              # sample scores
  labs( color = "Depth (cm)", fill= "Depth (cm)") + 
  coord_fixed() + #must have a fixed centroid
  #geom_segment(data = env_vector, # vector arrows
  #             aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
  #             arrow = arrow(length = unit(0.25, "cm")), colour = "gray70") +
  #geom_text(data = sso_nmds_merged_scrs, aes(x = MDS1, y = MDS2, label = sample), color = "black", # vector labels
  #          size = 3) +# sample scores
  geom_text(aes(x=-0.24, label="Stress: Type 1 weak ties\n 999 Permutations\nStress: 0.03", y=-0.042), colour="gray30", angle=0, text=element_text(size=7)) +
  scale_fill_manual(values = color_list_cpt) +
  scale_color_manual(values = color_list_cpt) +
  #scale_colour_brewer(palette = "YlGnBu") + 
  #scale_fill_brewer(palette= "YlGnBu")+
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="vertical",
        legend.position = "right",
        axis.ticks = element_line(color = "black"),
        #legend.position = c(0.46, 0.92),
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        panel.border = element_rect(colour = "black", size=0.8, fill = "transparent"))     # same axis scaling
cpt_NMDS
ggsave('~/CPT_sediment_geochem_NMDS.jpeg', bg = "transparent", width = 40, height = 30, units = "cm")
#end NMDS
#####

#geochem summary
cpt_numeric$depth_range <-as.factor(case_when((cpt_numeric$mid_depth_cm > 250) & (cpt_numeric$mid_depth_cm <= 450) ~ "201-450",
                                              (cpt_numeric$mid_depth_cm > 450) & (cpt_numeric$mid_depth_cm <= 650) ~ "451-650",
                                              (cpt_numeric$mid_depth_cm > 650) & (cpt_numeric$mid_depth_cm <= 750) ~ "601-850",
                                              (cpt_numeric$mid_depth_cm > 850) & (cpt_numeric$mid_depth_cm <= 1000) ~ "851-1000"))

colnames(cpt_numeric)
cpt_chem_mean <- cpt_numeric %>% # Summary by group using dplyr
  group_by(depth_range)%>%
  summarise_if(is.numeric, mean, na.rm = TRUE)
cpt_chem_sd <- cpt_numeric %>% # Summary by group using dplyr
  group_by(depth_range)%>%
  summarise_if(is.numeric, sd, na.rm = TRUE)
View(cpt_chem_mean)
View(cpt_chem_sd)
# alpha stat comparison
#####
p2017_Cm_sample_data <- as.data.frame(sample_data(physeqCPT))
colnames(p2017_Cm_sample_data)
row.names(p2017_Cm_sample_data)
rich_cm = estimate_richness(p2017_combined_merged, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))
row.names(rich_cm) <-row.names(p2017_Cm_sample_data)
merged_Cm_stat <-as.data.frame(merge(rich_cm, p2017_Cm_sample_data, by=0, all=TRUE))
head(merged_Cm_stat)

#anova of stats

#end alpha stat
#####

#Patesci abundance and geochem (not plotted)
head(melt_cpt_patesci)
colnames(melt_cpt_patesci)

melt_cpt_patesci_numeric <- dplyr::select(melt_cpt_patesci, c(3,12,15,17:21,26:28,30:61, 64:80, 82:91))#select numeric data
melt_cpt_patesci_numeric[is.na(melt_cpt_patesci_numeric)] <- 0
sapply(melt_cpt_patesci_numeric,class)
cpt_patesci_cor = cor(melt_cpt_patesci_numeric, method = "spearman")
corrplot(cpt_patesci_cor)
jpeg("~/patesci_chem_spearman_corrplot.jpeg",width=2700, height=2700, res=400, bg="white")
plot.new()
corrplot(cpt_patesci_cor, method="circle",#addCoef.col = 'gray40',
         tl.col = "black", col = COL2('PuOr'),type="lower",
         order = "hclust", hclust.method = "average", 
         addrect = 5, tl.cex = 0.7,diag=TRUE)
getOption("device")
dev.set(which = dev.next())

