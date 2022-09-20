# Code for Chapters 2 & 3 
# processing EVO data in phyloseq
# a.putt. 16 August 2021
# @putt-ad does not provide any warranty 
#####
# No guarantee is offered that this code will work on your version of R
# code provided here may be modified and updated to github at any time.
# the packages contained herewithin are a combination of proprietary and 
# open access code. Pipelines for these data may not apply to all data and 
# all situations. see: https://github.com/putt-ad/MemoryResponse for details.
#####

# Convert Galaxy output into phyloseq object
#otu_table - Works on any numeric matrix. You must also specify if the species are rows or columns
#sample_data - Works on any data.frame. The rownames must match the sample names in the otu_table if you plan to combine them as a phyloseq-object
#tax_table - Works on any character matrix. The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
#phyloseq - Takes as argument an otu_table and any unordered list of valid phyloseq components: sample_data, tax_table, phylo, or XStringSet. The tip labels of a phylo-object (tree) must match the OTU names of the otu_table, and similarly, the sequence names of an XStringSet object must match the OTU names of the otu_table.
#merge_phyloseq - Can take any number of phyloseq objects and/or phyloseq components, and attempts to combine them into one larger phyloseq object. This is most-useful for adding separately-imported components to an already-created phyloseq object.

# install phyloseq installer
source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",
       local = TRUE)

# install most recent version of phyloseq from github 
install_phyloseq(branch = "github")

# install phyloseq tool Biocmanager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")

# make packages available in your local library
# install phyloseq installer
source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",
       local = TRUE)
# install most recent version of phyloseq from github 
install_phyloseq(branch = "github")
library(phyloseq)
library(BiocManager)
list.of.packages <- c("plyr","dplyr","magrittr","tidyr","knitr","corrplot","ggcorrplot","vegan",
                      "tidyverse","tidytree","RColorBrewer","metacoder","ggpubr","ggplot2",
                      "ape","readr","foreach","doParallel","scales","grid","reshape2",
                      "multcompView", "viridis","VennDiagram","UpSetR","RColorBrewer",
                      "phyloseq","pheatmap", "base","plotly","vegan","ggpubr","ggtree", 
                      "ggfun", "ggjoy", "ggnewscale","RJSONIO","microeco","igraph")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages) # install
lapply(list.of.packages, require, character.only = TRUE) #  load
# install bioconductor packages
Bioc.packages <- c("microbial","microbiome","MicrobiotaProcess","DESeq2","dada2",
                   "DECIPHER","edgeR", "ANCOMBC","lattice")
new.Bioc.packages <- Bioc.packages[!(Bioc.packages %in% installed.packages()[,"Package"])]
if(length(new.Bioc.packages)) BiocManager::install(new.Bioc.packages) #install
Yes
lapply(Bioc.packages, require, character.only = TRUE) # load and print report of loaded packages
library(lattice)
library(microbiome)
library(MicrobiotaProcess)
library(broom)
library(tidyr)
library(car)
library(data.table)
library(microViz)
library(metagMisc)
BiocManager::install("dada2", version = "3.15")
install.packages("remotes")
remotes::install_github("vmikk/metagMisc",force=TRUE)

install.packages("iNEXT")
library(iNEXT)
#load Tax4Fun
install.packages(system.file("extdata", "biom_0.3.12.tar.gz", package="microeco"), repos = NULL, type = "source")
install.packages(system.file("extdata", "qiimer_0.9.4.tar.gz", package="microeco"), repos = NULL, type = "source")
install.packages(system.file("extdata", "Tax4Fun_0.3.1.tar.gz", package="microeco"), repos = NULL, type = "source")
install.packages(system.file("extdata", "Tax4Fun_0.3.1.tar.gz", package="microeco"), repos = NULL, type = "source")
library(Tax4Fun)
# Step 1. Load Data
load("merged_p2017.Rdata")
load("p2017_1_filtered.Rdata")
load("p2017_2_filtered.Rdata")
load("p2017_combined_filtered.Rdata")
load( "p2017_combined_merged.Rdata")
load("p2017_1.Rdata")
# sample data
meta_complete1 <- as.data.frame(read.csv("https://www.dropbox.com/s/rc7azdf4hsjd9qj/Unified_EVO_2017_metadatasheet.csv?dl=1", row.names=1, sep = ",", stringsAsFactors = FALSE))
meta_complete1
geochem_2017 <- as.data.frame(sample_data(phylo2017))
head(geochem_2017)

# increase memory size
memory.limit(size=60000)
memory.size()
otu_2017
plot(tree_2017)
plot(tree_2017_1)
plot(tree_2017_2)
colnames(meta_complete)

ggplot(data=meta_2017,aes(x=DO_mg.L,color=nitrate_reduction)) + 
  geom_density()

#Step 2. add / modify columns for analyses
#####
# make column for well_type
head(meta_complete)
#meta_complete$well_type <- meta_complete$WellID
#meta_complete["well_type"][meta_complete["well_type"] == c("FW215", "FW215")] <- "Control"
# set all non-control as 'monitoring'
#meta_complete["well_type"][meta_complete["well_type"] != "Control"] <- "Monitoring"
#distinct(geochem,well_type)
# make column for DO Content <0.1ppm and lower is 02 deplete, under 10ppm is slightly aerobic
#meta_complete$Oxygen_content <- meta_complete$DO_mg.L 
#meta_complete["Oxygen_content"][meta_complete["Oxygen_content"] < "0.1"] <- "O2 Deplete"
# set all non-control as 'monitoring'
#meta_complete["Oxygen_content"][meta_complete["Oxygen_content"] < "10"] <- "Slightly Aerobic"
# make filter size a text variable
#meta_complete["FilterSize"][meta_complete["FilterSize"] == "0.1"] <- "0.1um"
#meta_complete["FilterSize"][meta_complete["FilterSize"] == "0.2"] <- "0.2um"
#distinct(meta_complete,FilterSize)
# make column for categorical Nitrate
#meta_complete$nitrate_reduction <- meta_complete$Nitrate.mg.L
#meta_complete["nitrate_reduction"][meta_complete["nitrate_reduction"] >= "1"] <- "Not Reduced"
# set all non-control as 'monitoring'
#meta_complete["nitrate_reduction"][meta_complete["nitrate_reduction"] < "1"] <- "Reduced"
#distinct(meta_complete,nitrate_reduction)
#head(meta_complete)
head(meta_complete)
row.names(meta_complete) <- meta_complete$X
colnames(meta_complete) <- c('WellID','Days','FilterSize','EVOYear',
                             "pH",'DissolvedOxygen_mgL','SpecificConductivity_uScm','WaterTable_mAMSL',
                             'Magnesium_mgL', 'Aluminum_mgL', 'Potassium_mgL', "Calcium_mgL",
                             'Iron_mgL', 'Manganese_mgL', 'Uranium_mgL', 'Nitrate_mgL',
                             'Sulfate_mgL', 'Acetate_mgL', 'Ammonium_mgL','well_type',
                             'Oxygen_content','nitrate_reduction','DayCategory','sample_label',"sulfate_reduction",
                             "phase", "X.1", "X.2", "Distance")
#View(meta_complete)
# make column for categorical sulfate
meta_complete$sulfate_reduction <- meta_complete$SO4_mgL
meta_complete %>% mutate_at(vars("sulfate_reduction"), ~replace_na(.,"Reduced"))
meta_complete["sulfate_reduction"][meta_complete["sulfate_reduction"] >= "1"] <- "Not Reduced"
# set all non-control as 'monitoring'
meta_complete["sulfate_reduction"][meta_complete["sulfate_reduction"] < "1"] <- "Reduced"
meta_complete <-meta_complete %>% mutate_at(vars("sulfate_reduction"), ~replace_na(.,"Reduced"))
distinct(meta_complete,sulfate_reduction)

# update dropbox linked file
write.csv(meta_complete, "~/Unified_EVO_metadata_sheet.csv")
write.csv(meta_complete, "~/EVO_Unified_inprogress_paste-error-fix.csv")
#review the data
glimpse(meta_complete)
View(meta_complete)

# sample label
meta_complete$sample_label <-paste(meta_complete$well_type,meta_complete$DaysAfter,sep="|")
meta_complete$WellDay <- paste(meta_complete$WellID, meta_complete$Days, sep="")
head(meta_complete)


#write.csv(meta_complete, "Unified_EVO_metadata_sheet.csv")
#meta_complete <- read.csv("~/unified_evo_metadata_sheet.csv")

#geochem_2017 <- as.data.frame(read.csv("https://www.dropbox.com/s/aikypb9tlv7jt53/unified_2017_evo_metadata_sheet.csv?dl=1", sep = ",", row.names = 1, stringsAsFactors = FALSE ))
geochem_2017 <-  subset(meta_complete, EVOYear == "2017")
# subset data for phyloseq objects
geochem_2017_1 <- subset(meta_complete, FilterSize == "0.1um")
geochem_2017_1
geochem_2017_2 <- subset(geochem_2017 , FilterSize == "0.2um")
geochem_2017_2

#colnames(geochem_2017)
#change column names 
#colnames(geochem_2017) <- c('WellID','well_type',"group_sample","well_sample",
#                            'Oxygen_content','nitrate_reduction','Days','pH',
#                            'Spc.Cond_uScm','DissolvedOxygen_mgL','WaterLevel','Acetate_uM', 
#                            'Nitrate_mgL', 'Sulfate_mgL','Uranium_mgL', 'Ammonium_mgL', 
#                            "Magnesium_mgL",'Aluminum_mgL', 'Potassium_mgL', 'Calcium_mgL',
#                            'Iron_mgL','Manganese_mgL',"CellCount")
colnames(geochem_2017)
#subset categorical and numeric variables

geo2017_categorical <- select(geochem_2017, c(1,2,4, 20:26))
geochem_2017
geo2017_numeric <- geochem_2017[,5:19] #numeric
geo2017_scale_numeric <- as.data.frame(geo2017_numeric)
#geo2017_scale_numeric <- as.data.frame(geo2017_numeric[3:16]) #numeric to scale
colnames(geo2017_scale_numeric)
geo2017_summary_numeric <- geo2017_numeric
geo2017_summary_numeric$well_type<- geo2017_categorical$well_type
geo2017_summary_numeric["well_type"][geo2017_summary_numeric["well_type"] == "Monitoring"] <- 1
geo2017_summary_numeric["well_type"][geo2017_summary_numeric["well_type"] == "Control"] <- 0
geo2017_summary_numeric$Distance <- geochem_2017$WellID
geo2017_summary_numeric %>% mutate(Distance = recode(Distance, MLSB3 = 2.4, FW216 = 2.5, FW215 = -2.6, GP01 = 4.9, GP03 = 11.2))
meta_complete$Distance <- meta_complete$WellID
meta_complete %>% mutate(Distance = recode(Distance, MLSB3 = 2.4, FW216 = 2.5, FW215 = -2.6, GP01 = 4.9, GP03 = 11.2))

#site map
#####
list.map.packages <- c("ggplot2","ggmap","tidyverse","RColorBrewer","ggpubr",
                       "RCurl","corrplot","ggcorrplot", "Hmisc", "scales","qmap","plotly")

map.packages <- list.map.packages[!(list.map.packages %in% installed.packages()[,"Package"])]
if(length(map.packages)) install.packages(map.packages)

# load all these
lapply(list.map.packages, require, character.only = TRUE)

# Andrews google API
register_google(key="AIzaSyDbDh36jeyrUmPj2N9Vq0hpjkAqYLcaKAA")

#load the data
map_df <- read.csv("mapping_chemistry.csv", stringsAsFactors = FALSE)
sapply(map_df,class)
# set specific locations and elements
Area2 <- c(lon = -84.27471, lat = 35.975671)
Area2_df <- as.data.frame(t(Area2))
Area2 <-array(unlist(Area2))
plot.new()
theme_set(theme_bw(16))
EVOMap <- qmap(Area2, zoom = 22, maptype="satellite",color = "color", legend = "topright")
evowellplot <-EVOMap +
  geom_point(aes(x=Longitude, y=Latitude, colour = well_type, shape = well_type),size=7, alpha=0.8,
             data = map_df)+
  geom_text(data = map_df, aes(x=Longitude, y=Latitude, label = WellID),color = "White", position = position_dodge(width = 2),
            vjust = 1.5, hjust=0.5, size = 6, angle=355)+
  scale_color_brewer(palette= "Oranges")+
  #scale_color_distiller(palette = "YlGnBu")+
  labs(x = "Longitude", y = "Latitude", shape = "",color = "", size=" ")+
  #facet_wrap(~DaysPostinjection) +
  theme(panel.background = element_rect(fill = "gray95"), 
        legend.direction="vertical",
        legend.position = c(0.2, 0.2),
        legend.text = element_text(size=16),
        legend.title =  element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(size=18),
        strip.text.y = element_text(size=10,color="gray60"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size=10,color="gray60"))
evowellplot
ggsave('~/2017_well_map.png', bg = "transparent", width = 40, height = 30, units = "cm")

Y12 <- c(lon = -84.2735, lat = 35.9767)
Y12_df <- as.data.frame(t(Y12))
Y12 <-array(unlist(Y12))
Y12Map <- qmap(Y12, zoom = 18, maptype="satellite",color = "color", legend = "topright")
Y12plot <-Y12Map +
  geom_point(aes(x=lon, y=lat),color = "#edae49",size=6, alpha=1,
             data = Area2_df)+
  geom_text(data = Area2_df, aes(x=lon, y=lat),color = "black", label = "Area 2",position = position_dodge(width = 2),
            vjust = -1.1, hjust=0.5, size = 4.6, angle=0)+
  geom_text(data = Area2_df, aes(x=lon, y=lat),color = "white", label = "Area 2",position = position_dodge(width = 2),
            vjust = -1.1, hjust=0.5, size = 4.5, angle=0)+
  geom_text(aes(x=-84.2732, y=35.9779),color = "white", label = "S3 ponds",size = 7, angle=0)+
  #scale_color_brewer(palette= "Oranges")+
  #scale_color_distiller(palette = "YlGnBu")+
  labs(x = "Longitude", y = "Latitude", shape = "",color = "Area 2", size=" ")+
  #facet_wrap(~DaysPostinjection) +
  theme(panel.background = element_rect(fill = "gray95"), 
        legend.direction="vertical",
        legend.position = c(0.3, 0.26),
        legend.text = element_text(size=13),
        legend.title = element_text(size=10),
        plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(size=18),
        strip.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.4),
        axis.text = element_blank())
Y12plot
install.packages("patchwork")
library(patchwork)
library(cowplot)
evowellplot+ inset_element(Y12plot, left = 0.5, bottom = 0.5, right = 0.9, top = 1)
map.with.inset <-
  ggdraw() +
  draw_plot(evowellplot) +
  draw_plot(Y12plot, x = 0.115, y = .7, width = .3, height = .3)
map.with.inset
ggsave('~/2017_well_map.png', bg = "transparent", width = 40, height = 30, units = "cm")

#end of map
#####

#Figure 2.2 Correlations
######
# Geochem Correlations
#run with scaled data
geo2017_aov$Distance <- geo2017_summary_numeric$Distance
geo2017_aov$phase
colnames(geo2017_aov)
geo2017_cor <- select(geo2017_aov, c(2,8:16))
colnames(geo2017_cor)
colnames(geo2017_cor) <- c('Specific Conductivity','Iron',
                           'Manganese','Uranium', 'Nitrate', 'Sulfate', 'Acetate', 
                           'Ammonium', 'pH',"Days")
geo17_cor = cor(geo2017_cor, method = "spearman")
corrplot(geo17_cor)
jpeg("~/2017_chem_kendall_corrplot.jpg",width=2400, height=2200, res=400, bg="white")
plot.new()
corrplot(geo17_cor, method="circle",#addCoef.col = 'gray40',
         tl.col = "black", col = COL2('PuOr'),type="lower",
         order = "hclust", hclust.method = "average", addrect = 5, tl.cex = 0.7,diag=TRUE)
dev.off()
######

#Figure 2.1 Geochem Heatmap
#####
pheatmap(geo2017_scale_center)
colnames(geo2017_aov[1:14])
pheatmap(geo2017_aov[1:14])
pheatmap(geo2017_cor[4:10])
pheatmap(geo2017_log_numeric[8:12])
head(geo2017_log_numeric)
col.pal <- brewer.pal(9,"YlOrBr")
cluster_row <- hclust(dist(geo2017_scale_numeric))
dd <-dendsort::dendsort(cluster_row)
category_df = data.frame("Well" = geo2017_categorical$well_type)
#category_df$Filter <- geochem_2017$FilterSize
category_df$Phase <- geo2017_categorical$phase
unique(category_df$phase)
category_df$Days <- geo2017_numeric$Days


rownames(category_df) = rownames(geo2017_scale_center) # name matching
head(category_df)
category_df$Days <- factor(category_df$Days, levels = c("-6" , "1" ,  "8" ,  "15" , "22" , "50" , "78" , "106" ,"134"))
category_df$Phase <- factor(category_df$Phase, levels = c("Control", "Phase0",  "Phase1", "Phase2", "Phase3" , "Phase4",  "Phase5",  "Phase6"   ))

anno_colors = list(
  Days = c("-6" = "#ffffd9", "1" = "#edf8b1",  "8"  = "#c7e9b4",  "15" = "#7fcdbb", "22" = "#41b6c4" , "50"  = "#1d91c0", "78" = "#225ea8", "106" = "#253494","134" = "#081d58"),
  #Well = c(Control = "tomato3", Monitoring = "pink2"),
  #Well = c(Control = "darkred", Monitoring = "darksalmon"),
  Well = c(Control = "black", Monitoring = "blueviolet"),
  #Filter = c("0.1um" = "black", "0.2um" = "blueviolet"),
  Phase = c("Control" = "white","Phase0"=  "#ffffe5","Phase1"= "#f7fcb9" ,"Phase2"= "#d9f0a3","Phase3"= "#41ab5d" ,"Phase4"= "#238b45" ,"Phase5"= "#006d2c" ,"Phase6"= "#00441b" )
  #Filter = c("0.1um" = "palegreen", "0.2um" = "darkgreen")
)
#make pheatmap with day annotation
jpeg("~/2017_chem_ScaleCenter_heatmap.jpeg",width=2200, height=2800, res=400, bg="white")
plot.new()
pheatmap(
  mat               = geo2017_scale_center,
  color             = col.pal,
  border_color      = 'gray80',
  show_colnames     = TRUE,
  cutree_rows = 7,
  cutree_cols = 4,
  treeheight_row = 20,
  treeheight_col = 20,
  annotation_row = category_df,
  annotation_colors = anno_colors,
  angle_col=315,
  show_rownames     = FALSE,
  drop_levels       = TRUE,
  fontsize          = 12
)

getOption("device")
dev.set(which = dev.next())
dev.off()
######

#Data set for tables 2.1, 2.2 Geochem summary
#####
# summarize geochemistry means
geo2017_summary_type_output <- geo2017_summary_numeric %>%
  group_by(well_type) %>%
  summarize(mean_nitrate = mean(Nitrate_mgL),
            mean_sulfate = mean(Sulfate_mgL),
            mean_uranium = mean(Uranium_mgL),
            mean_iron = mean(Iron_mgL),
            mean_acetate = mean(Acetate_mgL))
geo2017_summary_typeDays_output
write.csv(geo2017_summary_typeDays_output, "~/2017_summary_typeDays_output.csv")
geo2017_summary_type_output
write.csv(geo2017_summary_type_output, "~/2017_summary_type_output.csv")
geo2017_table2.1.2 <- geo2017_summary_numeric %>%
  group_by(WellID, Days)
######

#dataset and testing data suitability for ANOVA - Supplementary QQ plot 
######
library(dplyr)
# transform
geo2017_scale_numeric[geo2017_scale_numeric == 0] <- NA
geo2017_log_numeric <- as.data.frame(log(geo2017_scale_numeric) ) 
geo2017_log_numeric[is.na(geo2017_log_numeric)] = 0
geo2017_log_numeric
geo2017_log_numeric$pH <- geochem_2017$pH
geo2017_log_numeric$Days <- as.character(geochem_2017$Days)
geo2017_log_numeric$well_type<- geochem_2017$well_type
geo2017_log_numeric["well_type"][geo2017_log_numeric["well_type"] == "Monitoring"] <- 1
geo2017_log_numeric["well_type"][geo2017_log_numeric["well_type"] == "Control"] <- 0
geo2017_log_numeric$Distance <- geochem_2017$WellID
geo2017_log_numeric %>% mutate(Distance = dplyr::recode(Distance, MLSB3 = 2.4, FW216 = 2.5, FW215 = -2.6, GP01 = 4.9, GP03 = 11.2))

geo2017_log_numeric <- as.data.frame(geo2017_log_numeric)
geo2017_aov <- geo2017_log_numeric
geo2017_aov$well_type<- geochem_2017$well_type
geo2017_aov["well_type"][geo2017_aov["well_type"] == "Monitoring"] <- 1
geo2017_aov["well_type"][geo2017_aov["well_type"] == "Control"] <- 0

sapply(geo2017_aov,class)
head(geo2017_aov)
#scale and center
colnames(geo2017_scale_numeric)
geo2017_scale_numeric$pH <- geo2017_scale_numeric$pH
geo2017_scale_numeric$Days <- geo2017_scale_numeric$Days
geo2017_scale_center <- select(geo2017_scale_numeric, c(1,2,3,6:10,15,16))
geo2017_scale_center <- scale(geo2017_scale_center, scale = TRUE) 
colnames(geo2017_scale_center)
colnames(geo2017_scale_center) <- c('pH','Specific Conductivity','Acetate', 
                                    'Nitrate', 'Sulfate','Uranium', 'Ammonium', 
                                    'Iron','Manganese')

#test difference in group variance
leveneTest(Nitrate_mgL~well_type,data = geo2017_aov)
#nitrate test of normality
NO3_aov <- aov(Nitrate_mgL ~ Days,data = geo2017_aov)
par(mfrow = c(1, 2));hist(NO3_aov$residuals);qqPlot(NO3_aov$residuals,id = FALSE) # id = FALSE to remove point identification
shapiro.test(NO3_aov$residuals) #p~0.05. NEAR NORMAL DISTRIBUTION
geo2017_aov %>% group_by(well_type) %>% shapiro_test(Nitrate_mgL) 
ggqqplot(geo2017_aov, "Nitrate_mgL", facet.by = "well_type")
#sulfate test of normality
SO4_aov <- aov(Sulfate_mgL ~ Days,data = geo2017_aov)
par(mfrow = c(1, 2));hist(SO4_aov$residuals);qqPlot(SO4_aov$residuals,id = FALSE) # id = FALSE to remove point identification
shapiro.test(SO4_aov$residuals) #p<0.05. NOT NORMAL DISTRIBUTION
geo2017_aov %>% group_by(well_type) %>% shapiro_test(Sulfate_mgL) 
ggqqplot(geo2017_aov, "Sulfate_mgL", facet.by = "well_type")
#uranium test of normality
U_aov <- aov(Uranium_mgL ~ Days,data = geo2017_aov)
par(mfrow = c(1, 2));hist(U_aov$residuals);qqPlot(U_aov$residuals,id = FALSE) # id = FALSE to remove point identification
shapiro.test(U_aov$residuals) #p<0.05. NOT NORMAL DISTRIBUTION
geo2017_aov %>% group_by(well_type) %>% shapiro_test(Uranium_mgL) 
ggqqplot(geo2017_aov, "Uranium_mgL", facet.by = "well_type")
#pHtest of normality
pH_aov <- aov(pH ~ Days,data = geo2017_aov)
par(mfrow = c(1, 2));hist(pH_aov$residuals);qqPlot(pH_aov$residuals,id = FALSE) # id = FALSE to remove point identification
shapiro.test(pH_aov$residuals) #p>>0.05. IS NORMAL DISTRIBUTION
#Iron test of normality
Fe_aov <- aov(Iron_mgL ~ Days,data = geo2017_aov)
par(mfrow = c(1, 2));hist(Fe_aov$residuals);qqPlot(Fe_aov$residuals,id = FALSE) # id = FALSE to remove point identification
shapiro.test(Fe_aov$residuals) #p<0.05. NOT NORMAL DISTRIBUTION
#Ammonium test of normality
NH4_aov <- aov(Ammonium_mgL ~ Days,data = geo2017_aov)
par(mfrow = c(1, 2));hist(NH4_aov$residuals);qqPlot(NH4_aov$residuals,id = FALSE) # id = FALSE to remove point identification
shapiro.test(NH4_aov$residuals) #p~0.05. NEAR NORMAL DISTRIBUTION
#Acetate test of normality
Ac_aov <- aov(Acetate_uM ~ Days,data = geo2017_aov)
par(mfrow = c(1, 2));hist(Ac_aov$residuals);qqPlot(Ac_aov$residuals,id = FALSE) # id = FALSE to remove point identification
shapiro.test(Ac_aov$residuals) #p<0.05. NOT NORMAL DISTRIBUTION
#Manganese test of normality
Mn_aov <- aov(Manganese_mgL ~ Days,data = geo2017_aov)
par(mfrow = c(1, 2));hist(Mn_aov$residuals);qqPlot(Mn_aov$residuals,id = FALSE) # id = FALSE to remove point identification
shapiro.test(Mn_aov$residuals) #p<0.05. NOT NORMAL DISTRIBUTION
#Manganese test of normality
Mg_aov <- aov(Magnesium_mgL ~ Days,data = geo2017_aov)
par(mfrow = c(1, 2));hist(Mg_aov$residuals);qqPlot(Mg_aov$residuals,id = FALSE) # id = FALSE to remove point identification
shapiro.test(Mg_aov$residuals) #p<0.05. NOT NORMAL DISTRIBUTION
#conductivity test of normality
spc_aov <- aov(Spc.Cond_uScm ~ Days,data = geo2017_aov)
par(mfrow = c(1, 2));hist(spc_aov$residuals);qqPlot(spc_aov$residuals,id = FALSE) # id = FALSE to remove point identification
shapiro.test(spc_aov$residuals) #p~0.05. NEAR NORMAL DISTRIBUTION
#Al test of normality
Al_aov <- aov(Aluminum_mgL ~ Days,data = geo2017_aov)
par(mfrow = c(1, 2));hist(Al_aov$residuals);qqPlot(Al_aov$residuals,id = FALSE) # id = FALSE to remove point identification
shapiro.test(Al_aov$residuals) #p<0.05. NOT NORMAL DISTRIBUTION
#DTW test of normality
DTW_aov <- aov(WaterLevel ~ Days,data = geo2017_aov)
par(mfrow = c(1, 2));hist(DTW_aov$residuals);qqPlot(DTW_aov$residuals,id = FALSE) # id = FALSE to remove point identification
shapiro.test(DTW_aov$residuals) #p<0.05. NOT NORMAL DISTRIBUTION
library(rstatix)
geo2017_aov %>% group_by(well_type) %>% shapiro_test(Nitrate_mgL)

#plot histograms and QQ plots of residuals
png('~/2017_chem_hist_by_day_plot.png',
    width = 850,
    height = 1100,
    units = "px", 
    pointsize = 12,
    bg = "white",
    res = 100)
par(mfrow = c(4, 3))
hist(NO3_aov$residuals);hist(SO4_aov$residuals);hist(U_aov$residuals);hist(Fe_aov$residuals);hist(NH4_aov$residuals);hist(Ac_aov$residuals);hist(Mn_aov$residuals);hist(Mg_aov$residuals);hist(Al_aov$residuals);hist(DTW_aov$residuals);hist(pH_aov$residuals);hist(spc_aov$residuals)
dev.off()
png('~/2017_chem_QQ-plot_by_day_plot.png',
    width = 850,
    height = 1100,
    units = "px", 
    pointsize = 12,
    bg = "white",
    res = 100)
par(mfrow = c(4, 3))
qqPlot(NO3_aov$residuals,id = FALSE);qqPlot(SO4_aov$residuals,id = FALSE);qqPlot(U_aov$residuals,id = FALSE);qqPlot(Fe_aov$residuals,id = FALSE);qqPlot(NH4_aov$residuals,id = FALSE);qqPlot(Ac_aov$residuals,id = FALSE);qqPlot(Mn_aov$residuals,id = FALSE);qqPlot(Mg_aov$residuals,id = FALSE);qqPlot(Al_aov$residuals,id = FALSE);qqPlot(DTW_aov$residuals,id = FALSE);qqPlot(pH_aov$residuals,id = FALSE);qqPlot(spc_aov$residuals,id = FALSE)
dev.off()
#####

#ANOVA Table 2.3 and Tukey 2.4
######
#add categories to log-transformed data
#convert categories to factors

head(geo2017_log_numeric)
geo2017_log_numeric$Days<- geo2017_categorical$Days
geo2017_log_numeric$well_type<- geo2017_categorical$well_type
geo2017_log_numeric$Distance <- geo2017_categorical$WellID
geo2017_log_numeric$Well <- geo2017_categorical$WellID
geo2017_log_numeric$sample_label <- geo2017_categorical$sample_label
head(geo2017_log_numeric)
geo2017_log_numeric <- geo2017_log_numeric %>% mutate(well_type = if_else(well_type == "Monitoring", true = 1, false = 0),
                                                      Distance = recode(Distance, MLSB3 = 2.4, FW216 = 2.5, FW215 = -2.6, GP01 = 4.9, GP03 = 11.2),
                                                      Well = recode(Well, MLSB3 = 1, FW216 = 2, FW215 = 0, GP01 = 3, GP03 = 4))
geo2017_log_numeric <- geo2017_log_numeric %>% mutate(well_type = factor(well_type, labels = c( "0", "1")),
                                                      Well = factor(Well, labels = c( "0", "1", "2", "3", "4")),
                                                      Distance = factor(Distance, c("-2.6", "2.4", "2.5", "4.9", "11.2")),
                                                      Days = factor(Days, c("-6", "1", "8", "15", "22", "50", "78", "106", "134")),
                                                      sample_label = factor(sample_label, label=c( "Control|-6"    , "Control|1"    ,  "Control|8"    ,  "Control|15"   ,  "Control|22"   ,  "Control|50"   ,  "Control|78"   ,  "Control|106"  ,  "Control|134"   ,
                                                                                                   "Monitoring|-6" , "Monitoring|1" ,  "Monitoring|8" ,  "Monitoring|15"  ,"Monitoring|22" , "Monitoring|50" , "Monitoring|78"  ,"Monitoring|106", "Monitoring|134"))
)
head(geo2017_log_numeric)
sulfate_aov_model <- aov(Sulfate  ~ Days   +well_type + Distance , data=geo2017_log_numeric)
summary(sulfate_aov_model)
nitrate_aov_model <- aov(Nitrate  ~  Days + well_type + Distance, data=geo2017_log_numeric)
summary(nitrate_aov_model)
iron_aov_model <- aov(Iron  ~ Days +well_type + Distance, data=geo2017_log_numeric)
summary(iron_aov_model)
uranium_aov_model <- aov(Uranium  ~ Days +well_type + Distance, data=geo2017_log_numeric)
summary(uranium_aov_model)
nh4_aov_model <- aov(Ammonium  ~ Days + well_type + Distance, data=geo2017_log_numeric)
summary(nh4_aov_model)
pH_aov_model <- aov(pH  ~ Days * well_type * Distance, data=geo2017_log_numeric)
spc_aov_model <- aov(SpecificConductivity_uScm  ~ Days + well_type +Distance, data=geo2017_log_numeric)
summary(spc_aov_model)
#view summary of three-way ANOVA
library(broom)
summary(sulfate_aov_model)
tidy_sulfate_type_tukey <- tidy(TukeyHSD(sulfate_aov_model, which ='well_type'))
tidy_sulfate_type_tukey
tidy_sulfate_days_tukey <- tidy(TukeyHSD(sulfate_aov_model, which ='Days'))
View(tidy_sulfate_days_tukey )
tidy_sulfate_dist_tukey <- tidy(TukeyHSD(sulfate_aov_model, which ='Distance'))
View(tidy_sulfate_dist_tukey)
summary(Nitrate_aov_model)
tidy_nitrate_well_tukey <- tidy(TukeyHSD(nitrate_aov_model, which ='well_type'))
tidy_nitrate_well_tukey
tidy_nitrate_day_tukey <- tidy(TukeyHSD(nitrate_aov_model, which ='Days'))
View(tidy_nitrate_day_tukey)
tidy_nitrate_dist_tukey <- tidy(TukeyHSD(nitrate_aov_model, which ='Distance'))
View(tidy_nitrate_dist_tukey)
summary(tidy_nitrate_day_tukey)
summary(iron_aov_model)
tidy_iron_well_tukey <- tidy(TukeyHSD(iron_aov_model, which ='well_type'))
tidy_iron_well_tukey
tidy_iron_day_tukey <- tidy(TukeyHSD(iron_aov_model, which ='Days'))
View(tidy_iron_day_tukey)
tidy_iron_dist_tukey <- tidy(TukeyHSD(iron_aov_model, which ='Distance'))
View(tidy_iron_dist_tukey)
summary(uranium_aov_model)
tidy_uranium_well_tukey <- tidy(TukeyHSD(uranium_aov_model, which ='well_type'))
tidy_uranium_well_tukey
tidy_uranium_dist_tukey <- tidy(TukeyHSD(uranium_aov_model, which ='Distance'))
View(tidy_uranium_dist_tukey)
summary(nh4_aov_model)
tidy_nh4_well_tukey <- tidy(TukeyHSD(nh4_aov_model, which ='well_type'))
tidy_nh4_well_tukey
tidy_nh4_dist_tukey <- tidy(TukeyHSD(nh4_aov_model, which ='Distance'))
View(tidy_nh4_dist_tukey)
summary(pH_aov_model)
tidy_pH_well_tukey <- tidy(TukeyHSD(pH_aov_model, which ='well_type'))
tidy_pH_well_tukey
summary(spc_aov_model)
tidy_spc_type_tukey <- tidy(TukeyHSD(spc_aov_model, which ='well_type'))
tidy_spc_type_tukey
tidy_spc_dist_tukey <- tidy(TukeyHSD(spc_aov_model, which ='Distance'))
View(tidy_spc_dist_tukey)

#subset only monitoring
geo2017_log_monitoring <-  subset(geo2017_log_numeric, FilterSize == "0.2um")

geo2017_log_numeric$Days<- geo2017_categorical$Days
geo2017_log_numeric$well_type<- geo2017_categorical$well_type
geo2017_log_numeric$Distance <- geo2017_categorical$WellID
geo2017_log_numeric$Well <- geo2017_categorical$WellID
geo2017_log_numeric$sample_label <- geo2017_categorical$sample_label
head(geo2017_log_numeric)
geo2017_log_monitoring <-  subset(geo2017_log_numeric, well_type == "1")
geo2017_log_monitoring <- geo2017_log_monitoring %>% mutate(well_type = if_else(well_type == "Monitoring", true = 1, false = 0),
                                                            Distance = recode(Distance, MLSB3 = 2.4, FW216 = 2.5, FW215 = -2.6, GP01 = 4.9, GP03 = 11.2),
                                                            Well = recode(Well, MLSB3 = 1, FW216 = 2, FW215 = 0, GP01 = 3, GP03 = 4))
geo2017_log_monitoring <- geo2017_log_monitoring %>% mutate(#well_type = factor(well_type, labels = c( "0", "1")),
  #Well = factor(Well, labels = c( "0", "1", "2", "3", "4")),
  #Distance = factor(Distance, c("-2.6", "2.4", "2.5", "4.9", "11.2")),
  Days = factor(Days, c("-6", "1", "8", "15", "22", "50", "78", "106", "134")),
  sample_label = factor(sample_label, label=c("Monitoring|-6" , "Monitoring|1" ,  "Monitoring|8" ,  "Monitoring|15"  ,"Monitoring|22" , "Monitoring|50" , "Monitoring|78"  ,"Monitoring|106", "Monitoring|134"))
)


unique(geo2017_log_monitoring$Days)
#sulfate anova - monitoring wells
sulfate_monitoring_model <- aov(Sulfate  ~ Days  +  Distance, data=geo2017_log_monitoring)
summary(sulfate_monitoring_model)
tidy_sulfate_type_monitor <- tidy(TukeyHSD(sulfate_monitoring_model, which ='Days'))
View(tidy_sulfate_type_monitor)
tidy_sulfate_dist_tukey <- tidy(TukeyHSD(sulfate_aov_model, which ='Distance'))
View(tidy_sulfate_dist_tukey)
#nitrate_monitoring
Nitrate_monitoring_model <- aov(Nitrate  ~  Days +  Distance, data=geo2017_log_monitoring)
summary(Nitrate_monitoring_model)
tidy_nitrate_model_tukey <- tidy(TukeyHSD(Nitrate_monitoring_model, which ='Days'))
View(tidy_nitrate_model_tukey)
#iron_modeling
summary(iron_monitor_model)
tidy_iron_day_tukey <- tidy(TukeyHSD(iron_aov_model, which ='Days'))
View(tidy_iron_day_tukey)
tidy_iron_model_tukey <- tidy(TukeyHSD(iron_monitor_model, which ='Days'))
View(tidy_iron_model_tukey)
#uranium_modeling
summary(iron_monitor_model)
tidy_iron_day_tukey <- tidy(TukeyHSD(iron_aov_model, which ='Days'))
View(tidy_iron_day_tukey)
uranium_monitoring_model <- aov(Uranium  ~ Days  +  Distance, data=geo2017_log_monitoring)
summary(uranium_monitoring_model)
tidy_iron_model_tukey <- tidy(TukeyHSD(uranium_monitoring_model, which ='Days'))
View(tidy_iron_model_tukey)
#ammonium_modeling
nh4_aov_model <- aov(Ammonium  ~ Days + well_type + Distance, data=geo2017_log_numeric)
summary(nh4_aov_model)
nh4_monitoring_model <- aov(Ammonium  ~ Days  +  Distance, data=geo2017_log_monitoring)
summary(nh4_monitoring_model)
tidy_nh4_model_tukey <- tidy(TukeyHSD(nh4_monitoring_model, which ='Days'))
View(tidy_nh4_model_tukey)
#spc_modeling
spc_monitor_model <- aov(SpecificConductivity_uScm  ~ Days  +Distance, data=geo2017_log_monitoring)
summary(spc_monitor_model)
nh4_monitoring_model <- aov(Ammonium  ~ Days  +  Distance, data=geo2017_log_monitoring)
summary(nh4_monitoring_model)
tidy_nh4_model_tukey <- tidy(TukeyHSD(nh4_monitoring_model, which ='Days'))
View(tidy_nh4_model_tukey)


iron_monitor_model <- aov(Iron  ~ Days + Distance, data=geo2017_log_monitoring)
uranium_aov_model <- aov(Uranium_mgL  ~ Days *well_type * Distance, data=geo2017_log_monitoring)
nh4_aov_model <- aov(Ammonium_mgL  ~ Days * well_type * Distance, data=geo2017_log_monitoring)
pH_aov_model <- aov(pH  ~ Days * well_type * Distance, data=geo2017_log_monitoring)
spc_aov_model <- aov(SpecificConductivity_uScm  ~ Days * well_type *Distance, data=geo2017_log_monitoring)
#anova
summary(sulfate_aov_model)
tidy_sulfate_days_tukey <- tidy(TukeyHSD(sulfate_aov_model, which ='Days'))
View(tidy_sulfate_days_tukey)
tidy_sulfate_dist_tukey <- tidy(TukeyHSD(sulfate_aov_model, which ='Distance'))
View(tidy_sulfate_dist_tukey)
summary(Nitrate_aov_model)
tidy_nitrate_well_tukey <- tidy(TukeyHSD(Nitrate_aov_model, which ='well_type'))
tidy_nitrate_well_tukey
tidy_nitrate_day_tukey <- tidy(TukeyHSD(Nitrate_aov_model, which ='Days'))
View(tidy_nitrate_day_tukey)
tidy_nitrate_dist_tukey <- tidy(TukeyHSD(Nitrate_aov_model, which ='Distance'))
View(tidy_nitrate_dist_tukey)
summary(tidy_nitrate_day_tukey)
summary(iron_aov_model)


tidy_iron_well_tukey <- tidy(TukeyHSD(iron_aov_model, which ='well_type'))
tidy_iron_well_tukey
tidy_iron_day_tukey <- tidy(TukeyHSD(iron_aov_model, which ='Days'))
View(tidy_iron_day_tukey)
tidy_iron_dist_tukey <- tidy(TukeyHSD(iron_aov_model, which ='Distance'))
View(tidy_iron_dist_tukey)
summary(uranium_aov_model)
tidy_uranium_well_tukey <- tidy(TukeyHSD(uranium_aov_model, which ='well_type'))
tidy_uranium_well_tukey
tidy_uranium_days_tukey <- tidy(TukeyHSD(uranium_aov_model, which ='Days'))
View(tidy_uranium_days_tukey)
summary(nh4_aov_model)
tidy_nh4_well_tukey <- tidy(TukeyHSD(nh4_aov_model, which ='well_type'))
tidy_nh4_well_tukey
tidy_nh4_days_tukey <- tidy(TukeyHSD(nh4_aov_model, which ='Days'))
View(tidy_nh4_days_tukey)
summary(pH_aov_model)
tidy_pH_well_tukey <- tidy(TukeyHSD(pH_aov_model, which ='well_type'))
tidy_pH_well_tukey
summary(spc_aov_model)
tidy_spc_days_tukey <- tidy(TukeyHSD(spc_aov_model, which ='Days'))
View(tidy_spc_days_tukey)
View(geo2017_numeric)
#####

#Figure 2.3 Geochem Plots
#####
#nitrate plot
NO3_summary <- geochem_2017 %>% # Summary by group using dplyr
  group_by(well_type, Days) %>% 
  summarize(mean = mean(Nitrate_mgL),
            std.dev = sd(Nitrate_mgL),
            stdmax = mean(Nitrate_mgL)+sd(Nitrate_mgL),
            stdmin = mean(Nitrate_mgL)-sd(Nitrate_mgL),
            range = max(Nitrate_mgL)-min(Nitrate_mgL))
NO3_summary[is.na(NO3_summary)] = 0
View(NO3_summary)
nitrate_plot <-ggplot( NO3_summary, aes(x=Days) ) + 
  geom_line( aes(y=mean,group=well_type,color=well_type) ) + 
  geom_point( aes(y=mean,color = well_type, shape = well_type),size=4, alpha = 0.8) +
  geom_linerange(aes(ymin = stdmin , ymax = stdmax, color = well_type), linetype = 3, alpha = 0.3) +
  geom_line( aes(x=Days,y=stdmax,group=well_type,color=well_type), linetype = 2, alpha = 0.5) + 
  geom_line( aes(x=Days,y=stdmin,group=well_type,color=well_type), linetype = 2, alpha=0.5) + 
  geom_ribbon(aes(x=Days,ymin=stdmin,ymax=stdmax, fill = well_type), alpha=0.1) +
  theme_test() +
  geom_hline(yintercept = 0, linetype = 1, color = 'black', size = 0.5) +
  geom_vline(xintercept = 0, linetype = 1, color = "gray50", size = 5.5, alpha = 0.4) +
  geom_text(aes(x=1.8, label="Injection of EVO\n", y=23), colour="gray30", angle=90, text=element_text(size=8)) +
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = c(0.46, 0.92),
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        panel.border = element_rect(colour = "black", size=0.8))+
  scale_x_continuous( breaks = unique(NO3_summary$`Days`)) +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  #xlim(-10,250) +
  labs(x = "Days Post-Injection", y = "Nitrate (mg/L)", fill = "",color = "", shape = "", linetype = "", title = "A. Nitrate")
nitrate_plot

#sulfate plot
SO4_summary <- geochem_2017 %>% # Summary by group using dplyr
  group_by(well_type, Days) %>% 
  summarize(mean = mean(Sulfate_mgL),
            std.dev = sd(Sulfate_mgL),
            stdmax = mean(Sulfate_mgL)+sd(Sulfate_mgL),
            stdmin = mean(Sulfate_mgL)-sd(Sulfate_mgL),
            range = max(Sulfate_mgL)-min(Sulfate_mgL))
SO4_summary[is.na(SO4_summary)] = 0
sulfate_plot <- ggplot( SO4_summary, aes(x=Days) ) + 
  geom_ribbon(aes(x=Days,ymin=stdmin,ymax=stdmax, fill = well_type), alpha=0.1) +
  geom_line( aes(x=Days,y=stdmax,group=well_type,color=well_type), linetype = 2, alpha = 0.5) + 
  geom_line( aes(x=Days,y=stdmin,group=well_type,color=well_type), linetype = 2, alpha=0.5) + 
  geom_linerange(aes(ymin = stdmin , ymax = stdmax, color = well_type), linetype = 3, alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = 1, color = 'black', size = 0.5) +
  geom_vline(xintercept = 0, linetype = 1, color = "gray50", size = 5.5, alpha = 0.2) +
  geom_text(aes(x=1.8, label="Injection of EVO\n", y=37), colour="gray30", angle=90, text=element_text(size=8)) +
  geom_line( aes(y=mean,group=well_type,color=well_type) ) + 
  geom_point( aes(y=mean,color = well_type, shape = well_type),size=4, alpha = 0.8) +
  theme_test() +
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = c(0.8, 0.96),
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        panel.border = element_rect(colour = "black", size=0.8))+
  scale_x_continuous( breaks = unique(SO4_summary$`Days`)) +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  #xlim(-10,250) +
  labs(x = "Days Post-Injection", y = "Sulfate (mg/L)", fill = "",color = "", shape = "", linetype = "", title = "B. Sulfate")
sulfate_plot
colnames(sample_data(phylo2017_rarefied))

#uranium plot
uranium_summary <- geochem_2017 %>% # Summary by group using dplyr
  group_by(well_type, Days) %>% 
  summarize(mean = mean(Uranium_mgL),
            std.dev = sd(Uranium_mgL),
            stdmax = mean(Uranium_mgL)+sd(Uranium_mgL),
            stdmin = mean(Uranium_mgL)-sd(Uranium_mgL),
            range = max(Uranium_mgL)-min(Uranium_mgL))
uranium_summary[is.na(uranium_summary)] = 0
uranium_plot <- ggplot( uranium_summary, aes(x=Days) ) + 
  geom_ribbon(aes(x=Days,ymin=stdmin,ymax=stdmax, fill = well_type), alpha=0.1) +
  geom_line( aes(x=Days,y=stdmax,group=well_type,color=well_type), linetype = 2, alpha = 0.5) + 
  geom_line( aes(x=Days,y=stdmin,group=well_type,color=well_type), linetype = 2, alpha=0.5) + 
  geom_linerange(aes(ymin = stdmin , ymax = stdmax, color = well_type), linetype = 3, alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = 1, color = 'black', size = 0.5) +
  geom_vline(xintercept = 0, linetype = 1, color = "gray50", size = 5.5, alpha = 0.2) +
  geom_text(aes(x=1.8, label="Injection of EVO\n", y=0.65), colour="gray30", angle=90, text=element_text(size=8)) +
  geom_line( aes(y=mean,group=well_type,color=well_type) ) + 
  geom_point( aes(y=mean,color = well_type, shape = well_type),size=4, alpha = 0.8) +
  theme_test() +
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = c(0.8, 0.96),
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        panel.border = element_rect(colour = "black", size=0.8))+
  scale_x_continuous( breaks = unique(uranium_summary$`Days`)) +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  #xlim(-10,250) +
  labs(x = "Days Post-Injection", y = "Uranium (mg/L)", fill = "",color = "", shape = "", linetype = "", title = "C. Uranium")
uranium_plot

#iron plot
ac_summary <- geochem_2017 %>% # Summary by group using dplyr
  group_by(well_type, Days) %>% 
  summarize(mean = mean(Acetate_mgL),
            std.dev = sd(Acetate_mgL),
            stdmax = mean(Acetate_mgL)+sd(Acetate_mgL),
            stdmin = mean(Acetate_mgL)-sd(Acetate_mgL),
            range = max(Acetate_mgL)-min(Acetate_mgL))
ac_summary[is.na(ac_summary)] = 0
ac_plot <- ggplot( ac_summary, aes(x=Days) ) + 
  geom_ribbon(aes(x=Days,ymin=stdmin,ymax=stdmax, fill = well_type), alpha=0.1) +
  geom_line( aes(x=Days,y=stdmax,group=well_type,color=well_type), linetype = 2, alpha = 0.5) + 
  geom_line( aes(x=Days,y=stdmin,group=well_type,color=well_type), linetype = 2, alpha=0.5) + 
  geom_linerange(aes(ymin = stdmin , ymax = stdmax, color = well_type), linetype = 3, alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = 1, color = 'black', size = 0.5) +
  geom_vline(xintercept = 0, linetype = 1, color = "gray50", size = 5.5, alpha = 0.2) +
  geom_text(aes(x=1.8, label="Injection of EVO\n", y=29.9), colour="gray30", angle=90, text=element_text(size=8)) +
  geom_line( aes(y=mean,group=well_type,color=well_type) ) + 
  geom_point( aes(y=mean,color = well_type, shape = well_type),size=4, alpha = 0.8) +
  theme_test() +
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = c(0.8, 0.96),
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        panel.border = element_rect(colour = "black", size=0.8))+
  scale_x_continuous( breaks = unique(ac_summary$`Days`)) +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  #xlim(-10,250) +
  labs(x = "Days Post-Injection", y = "Acetate (mg/L)", fill = "",color = "", shape = "", linetype = "", title = "D. Acetate")
ac_plot

#ammonium plot
nh4_summary <- geochem_2017 %>% # Summary by group using dplyr
  group_by(well_type, Days) %>% 
  summarize(mean = mean(Ammonium_mgL),
            std.dev = sd(Ammonium_mgL),
            stdmax = mean(Ammonium_mgL)+sd(Ammonium_mgL),
            stdmin = mean(Ammonium_mgL)-sd(Ammonium_mgL),
            range = max(Ammonium_mgL)-min(Ammonium_mgL))
nh4_summary[is.na(nh4_summary)] = 0
nh4_plot <- ggplot( nh4_summary, aes(x=Days) ) + 
  geom_ribbon(aes(x=Days,ymin=stdmin,ymax=stdmax, fill = well_type), alpha=0.1) +
  geom_line( aes(x=Days,y=stdmax,group=well_type,color=well_type), linetype = 2, alpha = 0.5) + 
  geom_line( aes(x=Days,y=stdmin,group=well_type,color=well_type), linetype = 2, alpha=0.5) + 
  geom_linerange(aes(ymin = stdmin , ymax = stdmax, color = well_type), linetype = 3, alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = 1, color = 'black', size = 0.5) +
  geom_vline(xintercept = 0, linetype = 1, color = "gray50", size = 5.5, alpha = 0.2) +
  geom_text(aes(x=1.8, label="Injection of EVO\n", y=0.3), colour="gray30", angle=90, text=element_text(size=8)) +
  geom_line( aes(y=mean,group=well_type,color=well_type) ) + 
  geom_point( aes(y=mean,color = well_type, shape = well_type),size=4, alpha = 0.8) +
  theme_test() +
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = c(0.8, 0.96),
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        panel.border = element_rect(colour = "black", size=0.8))+
  scale_x_continuous( breaks = unique(nh4_summary$`Days`)) +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  #xlim(-10,250) +
  labs(x = "Days Post-Injection", y = "Ammonium (mg/L)", fill = "",color = "", shape = "", linetype = "", title = "E. Ammonium")
nh4_plot

#spc plot
spc_summary <- geochem_2017 %>% # Summary by group using dplyr
  group_by(well_type, Days) %>% 
  summarize(mean = mean(SpecificConductivity_uScm),
            std.dev = sd(SpecificConductivity_uScm),
            stdmax = mean(SpecificConductivity_uScm)+sd(SpecificConductivity_uScm),
            stdmin = mean(SpecificConductivity_uScm)-sd(SpecificConductivity_uScm),
            range = max(SpecificConductivity_uScm)-min(SpecificConductivity_uScm))
spc_summary[is.na(spc_summary)] = 0
spc_plot <- ggplot( spc_summary, aes(x=Days) ) + 
  geom_ribbon(aes(x=Days,ymin=stdmin,ymax=stdmax, fill = well_type), alpha=0.1) +
  geom_line( aes(x=Days,y=stdmax,group=well_type,color=well_type), linetype = 2, alpha = 0.5) + 
  geom_line( aes(x=Days,y=stdmin,group=well_type,color=well_type), linetype = 2, alpha=0.5) + 
  geom_linerange(aes(ymin = stdmin , ymax = stdmax, color = well_type), linetype = 3, alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = 1, color = 'black', size = 0.5) +
  geom_vline(xintercept = 0, linetype = 1, color = "gray50", size = 5.5, alpha = 0.2) +
  geom_text(aes(x=1.8, label="Injection of EVO\n", y=485), colour="gray30", angle=90, text=element_text(size=8)) +
  geom_line( aes(y=mean,group=well_type,color=well_type) ) + 
  geom_point( aes(y=mean,color = well_type, shape = well_type),size=4, alpha = 0.8) +
  theme_test() +
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = c(0.8, 0.96),
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        panel.border = element_rect(colour = "black", size=0.8))+
  scale_x_continuous( breaks = unique(spc_summary$`Days`)) +
  scale_y_continuous(limits = c(0, 800), breaks = c(0, 200, 400, 600, 800))+
  #scale_y_continuous(breaks=seq(600,800,by=200)) +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  #xlim(-10,250) +
  labs(x = "Days Post-Injection", y = "Specific Conductivity (µS/cm)", fill = "",color = "", shape = "", linetype = "", title = "F. Specific Conductivity")
spc_plot

#make geochem plots
ggarrange(nitrate_plot, sulfate_plot, uranium_plot, ac_plot, nh4_plot, spc_plot ,widths = c(3,3,3))
ggsave('~/2017_chem_line_plot.png', bg = "transparent", width = 40, height = 30, units = "cm")
#####


# Step 4. Make Phyloseq object and check effects of rarefying on taxa
#####
# coerce taxonomy data into a phyloseq format
taxamat <- dplyr::select(taxa_2017, -ends_with(c(".Confidence"))) # remove confidence scores
names(taxamat)[names(taxamat) == "Kingdom"] <- "Domain" # phyloseq uses 'Domain' instead of 'Kingdom'
taxamat <-as.matrix(taxamat)
#Remove NA genus
taxamat_NA<- as.matrix(taxamat[!is.na(taxamat_df$Genus),])
#taxamat_NA[is.na(taxamat_NA)] = 0
#remove unknown phyla
taxamat <- as.matrix(taxamat[!is.na(taxamat$Phylum),])
# sample data
otu_2017
taxamat
geochem_2017 <- geochem_2017_2
phylo2017 <- phyloseq(otu_table(otu_2017, taxa_are_rows = TRUE), phyloseq::tax_table(taxamat), 
                      sample_data(meta_complete), tree_2017)
phylo2017
#rarefy based on coverage
set.seed(0451)
phylo2017_cov_raref <- phyloseq_coverage_raref(physeq = phylo2017, iter = 1, coverage = 0.98)
phyloseq_coverage(phylo2017_cov_raref)
phylo2017_cov_raref
# rarefy by even depth. Removed 2866 OTUs
set.seed(0451)
phylo2017_rarefied<-rarefy_even_depth(phylo2017,sample.size = min(sample_sums(phylo2017)))
phylo2017_rarefied
phylo2017_rarefied_tree <- phyloseq::phy_tree(phylo2017_rarefied)
#line plot of alpha diversity
ggrare <- function(physeq_object, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  
  x <- methods::as(phyloseq::otu_table(physeq_object), "matrix")
  if (phyloseq::taxa_are_rows(physeq_object)) { x <- t(x) }
  
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  
  # Get sample data
  if (!is.null(phyloseq::sample_data(physeq_object, FALSE))) {
    sdf <- methods::as(phyloseq::sample_data(physeq_object), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  # Add, any custom-supplied plot-mapped variables
  if ( length(color) > 1 ) {
    data$color <- color
    names(data)[names(data) == "color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  
  if ( length(label) > 1 ) {
    labels$label <- label
    names(labels)[names(labels) == "label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot2::ggplot(data = data,
                       ggplot2::aes_string(x = "Size",
                                           y = ".S",
                                           group = "Sample",
                                           color = color))
  
  p <- p + ggplot2::labs(x = "Sequence Sample Size", y = "Species Richness")
  
  if (!is.null(label)) {
    p <- p + ggplot2::geom_text(data = labels,
                                ggplot2::aes_string(x = "x",
                                                    y = "y",
                                                    label = label,
                                                    color = color),
                                size = 4, hjust = 0)
  }
  
  p <- p + ggplot2::geom_line()
  if (se) { ## add standard error if available
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = ".S - .se",
                                               ymax = ".S + .se",
                                               color = NULL,
                                               fill = color),
                           alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}
p2017 <- ggrare(phylo2017, step = 1000, color = "WellID", se = FALSE)
(p2017 <- p2017 + facet_wrap(~FilterSize) +theme(axis.text.x = element_text(angle = 90))+ labs(title = "Non-rarefied"))

cov2017 <- ggrare(phylo2017_cov_raref, step = 1000, color = "WellID",se = FALSE)
(cov2017 <- cov2017 + facet_wrap(~FilterSize) +theme(axis.text.x = element_text(angle = 90))+ labs(title = "Coverage-rarefied"))

rare2017 <- ggrare(phylo2017_rarefied, step = 1000, color = "WellID",se = FALSE)
(rare2017 <- rare2017 + facet_wrap(~FilterSize) +theme(axis.text.x = element_text(angle = 90))+ labs(title = "Even-depth-rarefied"))

ggarrange(p2017, cov2017, rare2017, widths = c(3,3,3))
ggsave('~/2017_rarefied_taxa_lineplot.png', bg = "transparent", width = 40, height = 30, units = "cm")
dev.off()
plot.new()

#step 5. venn diagram
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)
library(MicrobiotaProcess)
library(VennDiagram)
vennlist_filt <- get_vennlist(obj=phylo2017, factorNames="FilterSize")
View(vennlist_filt)
glimpse(vennlist_filt)
venn_filt_diagram <- venn.diagram(vennlist_filt,
                                  height=5,
                                  width=5, 
                                  filename=NULL, 
                                  fill=c("salmon", "steelblue"),
                                  cat.col=c("salmon", "steelblue"),
                                  alpha = 0.4, 
                                  fontfamily = "serif",
                                  fontface = "bold",
                                  main = "Even-depth rarefied",
                                  cex = 1.2,
                                  cat.cex = 1.3,
                                  cat.default.pos = "outer",
                                  cat.dist=0.1,
                                  margin = 0.1, 
                                  lwd = 3,
                                  lty ='dotted',
                                  imagetype = "svg")

venn_filt_diagram
vennlist_filt <- get_vennlist(obj=phylo2017, factorNames="FilterSize")
glimpse(vennlist_filt)
venn_cov <- venn.diagram(vennlist_filt_cov,
                         height=5,
                         width=5, 
                         filename=NULL, 
                         fill=c("salmon", "steelblue"),
                         cat.col=c("salmon", "steelblue"),
                         alpha = 0.4, 
                         fontfamily = "serif",
                         fontface = "bold",
                         main = "Coverage-rarefied ASVs",
                         cex = 1.2,
                         cat.cex = 1.3,
                         cat.default.pos = "outer",
                         cat.dist=0.1,
                         margin = 0.1, 
                         lwd = 3,
                         lty ='dotted',
                         imagetype = "svg")
vennlist_filt <- get_vennlist(obj=phylo2017, factorNames="FilterSize")
vennlist_filt
#extract data from the venn diagram output
taxamat
venn_filt_1<- as.data.frame(vennlist_filt$`0.1um`)
row.names(venn_filt_1) <- venn_filt_1[,1]
venn_filt_1 <- as.data.frame(venn_filt_1)
venn_filt_2<- as.data.frame(vennlist_filt$`0.2um`)
row.names(venn_filt_2) <- venn_filt_2[,1]
venn_filt_2 <- as.data.frame(venn_filt_2)
View(venn_filt_1)
View(venn_filt_2)
dplyr::inner_join(venn_filt_1 , taxamat)
venn_filt_both <- as.data.frame(generics::intersect(venn_filt_1 , venn_filt_2))
row.names(venn_filt_1) <- venn_filt_1[,1]
taxamat_df <- as.data.frame(taxamat)
taxamat_df
filtered_0.1_taxa <- as.data.frame(merge(x=venn_filt_1, y=taxamat_df, by = 0, all.x = TRUE, all.y = FALSE))
filtered_0.1_taxa 
filtered_0.1_taxa <- filtered_0.1_taxa[,-1]; filt_0.1_taxa <- filtered_0.1_taxa[,-1]
row.names(filt_0.1_taxa) <- filtered_0.1_taxa[,1]
filt_0.1_taxa
write.csv(filt_0.1_taxa, '~/2017_EVO_0.1_venn_taxa.csv')
filt_0.1_taxa  %>% count(Phylum, sort = TRUE)
filt_0.1_taxa  %>% count(Phylum, Class, sort = TRUE)

venn_filt_2<- as.data.frame(vennlist_filt$`0.2um`)
filtered_0.2_taxa <- merge(venn_filt_2, taxamat_df, by = 0, all.x = TRUE, all.y = FALSE)
filtered_0.2_taxa <- filtered_0.2_taxa[,-1]; filt_0.2_taxa <- filtered_0.2_taxa[,-1]
row.names(filt_0.2_taxa) <- filtered_0.2_taxa[,1]
View(filt_0.2_taxa)
write.csv(filt_0.2_taxa, '~/2017_EVO_0.2_venn_taxa.csv')
filt_0.2_taxa  %>% count(Phylum, sort = TRUE)
filt_0.2_taxa  %>% count(Phylum, Class, sort = TRUE)
filtered_both_taxa <- merge(venn_filt_1, venn_filt_2, by = 0, all.x = FALSE, all.y = FALSE)
filtered_both_taxa <- merge(filtered_both_taxa, taxamat, by.x=2, by.y= 0, all.x = FALSE, all.y = FALSE)
write.csv(filtered_both_taxa, '~/2017_EVO_both_venn_taxa.csv')
filtered_both_taxa  %>% count(Phylum, sort = TRUE)
filtered_both_taxa  %>% count(Phylum, Class, sort = TRUE)
venn_p2017 <- venn.diagram(vennlist_filt,
                           height=5,
                           width=5, 
                           filename=NULL, 
                           fill=c("salmon", "steelblue"),
                           cat.col=c("salmon", "steelblue"),
                           alpha = 0.4, 
                           fontfamily = "serif",
                           fontface = "bold",
                           main = "Coverage-rarefied ASVs",
                           cex = 1.2,
                           cat.cex = 1.3,
                           cat.default.pos = "outer",
                           cat.dist=0.1,
                           margin = 0.1, 
                           lwd = 3,
                           lty ='dotted',
                           imagetype = "svg")
#subset the phylo2017 
p2017_taxtable <- as.data.frame(phyloseq::tax_table(phylo2017))
p2017_otu <- as.data.frame(phyloseq::otu_table(phylo2017))
p2017_sample_data <- sample_data(phylo2017)
head(sample_data(p2017_1))

p2017_1 <- subset_samples(phylo2017, FilterSize=="0.1um")
p2017_1 = subset_samples(p2017_1, Date != "8/16/18")#remove post-injection samples
p2017_1 <- prune_taxa(taxa_sums(p2017_1) > 0, p2017_1)
save(p2017_1, file = "~/p2017_1.Rdata")
p2017_1_filtered = prune_taxa(taxa_sums(p2017_1_filtered) >= 1, p2017_1_filtered) 
p2017_1_top_filtered = prune_taxa(taxa_sums(p2017_1_filtered) > 20, p2017_1_filtered) 
taxtable
#0.2
p2017_2_filtered <- subset_samples(phylo2017, FilterSize=="0.2um")
p2017_2_filtered = prune_taxa(taxa_sums(p2017_2_filtered) >= 1, p2017_2_filtered) 
p2017_2_top_filtered = prune_taxa(taxa_sums(p2017_2_filtered) > 20, p2017_2_filtered) 
p2017_2_filtered
#taxonomy table for phyloseq
p2017_1_taxtable <- phyloseq::tax_table(p2017_1_filtered)
p2017_1_taxtable_df <-as.data.frame(p2017_1_taxtable)
p2017_1_taxtable_rows <- as.data.frame(row.names(p2017_1_taxtable_df))
rownames(p2017_1_taxtable_rows) <- p2017_1_taxtable_rows[,1]
p2017_1_taxtable_rows 
p2017_1_taxtable_df
p2017_1_otutable <- phyloseq::otu_table(p2017_1_filtered)
p2017_1_otutable_df <-as.data.frame(p2017_1_otutable)
p2017_1_otutable
p2017_1_samples <- phyloseq::sample_data(p2017_1_filtered)
#subset OTU and taxonomy for 0.2
p2017_2_taxtable <- phyloseq::tax_table(p2017_2_filtered)
p2017_2_taxtable_df <-as.data.frame(p2017_2_taxtable)
p2017_2_taxtable_df
p2017_2_taxtable_rows <- as.data.frame(row.names(p2017_2_taxtable_df))
rownames(p2017_2_taxtable_rows) <- p2017_2_taxtable_rows[,1]
p2017_2_taxtable_rows 
p2017_2_otutable <- phyloseq::otu_table(p2017_2_filtered)
p2017_2_otutable_df <-as.data.frame(p2017_2_otutable)
p2017_2_filtered
p2017_2_samples <- phyloseq::sample_data(p2017_2_filtered)

#combined taxa matrix
#merge all ASV IDs into a single dataframe
total_tax_id <- p2017_2_taxtable_rows
total_tax_id<- merge(total_tax_id, p2017_1_taxtable_rows, by=0, all=TRUE)
row.names(total_tax_id) <- total_tax_id$Row.names
head(total_tax_id)
combined_taxID <- subset(total_tax_id, total_tax_id$`row.names(p2017_2_taxtable_df)` == total_tax_id$`row.names(p2017_1_taxtable_df)`)
row.names(combined_taxID) <- combined_taxID$`row.names(p2017_2_taxtable_df)`
combined_taxID <- as.data.frame(combined_taxID[2]);head(combined_taxID); nrow(combined_taxID)
View(combined_taxID)
head(p2017_taxtable)
combined_taxamat1<- merge(combined_taxID, p2017_taxtable, by=0, all = TRUE)#add taxonomy information
head(combined_taxamat1);nrow(combined_taxamat1)
combined_taxamat1 <-combined_taxamat1[!is.na(combined_taxamat1$`row.names(p2017_2_taxtable_df)`),]
head(combined_taxamat1);nrow(combined_taxamat1)
row.names(combined_taxamat1) <- combined_taxamat1$Row.names;combined_taxamat1 <- as.matrix(combined_taxamat1[3:9]);head(combined_taxamat1); nrow(combined_taxamat1)

#combined otu
combined_otu<- as.data.frame(merge(combined_taxID, p2017_otu, by=0, all = TRUE))#add taxonomy information
combined_otu <-combined_otu[!is.na(combined_otu$`row.names(p2017_2_taxtable_df)`),]
head(combined_otu);ncol(combined_otu);nrow(combined_otu)
row.names(combined_otu) <- combined_otu$`row.names(p2017_2_taxtable_df)`;combined_otu <- as.matrix(combined_otu[3:77]);head(combined_otu); nrow(combined_otu)

#combined phyloseq
p2017_combined_filtered<- phyloseq(otu_table(combined_otu, taxa_are_rows = TRUE), phyloseq::tax_table(combined_taxamat1), 
                                   sample_data(p2017_sample_data))
#combined sequence tree
combined_tree = rtree(ntaxa(p2017_combined_filtered), rooted=TRUE, tip.label=taxa_names(p2017_combined_filtered))
plot(combined_tree)
p2017_combined_filtered<- phyloseq(otu_table(combined_otu, taxa_are_rows = TRUE), phyloseq::tax_table(combined_taxamat1), 
                                   sample_data(p2017_sample_data), combined_tree)
p2017_combined_merged  <- merge_samples(p2017_combined_filtered, "WellDay")
head(sample_data(p2017_combined_merged))
save(p2017_combined_filtered, file = "~/p2017_combined_filtered.Rdata")
save(p2017_combined_merged, file = "~/p2017_combined_merged.Rdata")



#taxa merge (known taxa)
combined_intersecting_taxa <- as.data.frame(generics::intersect(p2017_1_taxtable_df, p2017_2_taxtable_df))
combined_intersecting_taxa_matrix <- as.matrix(combined_intersecting_taxa)

#0.1um taxa matrix
UMB_taxID1 <- subset(total_tax_id,  is.na(total_tax_id$`row.names(p2017_2_taxtable_df)`))#exclude missing 0.2
UMB_taxID <- UMB_taxID1[3]
head(UMB_taxID);nrow(UMB_taxID)#review output with new rownames
UMB_taxamat<- as.data.frame(merge(UMB_taxID, p2017_1_taxtable_df, by=0, all = TRUE))#add taxonomy information
UMB_taxamat <-UMB_taxamat[!is.na(UMB_taxamat$`row.names(p2017_1_taxtable_df)`),];head(UMB_taxamat);ncol(UMB_taxamat)
row.names(UMB_taxamat) <- UMB_taxamat$`row.names(p2017_1_taxtable_df)`;UMB_taxamat <- as.matrix(UMB_taxamat[3:9]);head(UMB_taxamat); nrow(UMB_taxamat)
#0.1um otu dataframe
UMB_otu<- as.data.frame(merge(UMB_taxID, p2017_1_otutable_df, by=0, all = TRUE))#add taxonomy information
head(UMB_otu);ncol(UMB_otu)
UMB_otu <-UMB_otu[!is.na(UMB_otu$`row.names(p2017_1_taxtable_df)`),]
row.names(UMB_otu) <- UMB_otu$`row.names(p2017_1_taxtable_df)`;UMB_otu <- as.data.frame(UMB_otu[3:33]);head(UMB_otu); nrow(UMB_otu)

#0.1 phyloseq
p2017_1_filtered<- phyloseq(otu_table(UMB_otu, taxa_are_rows = TRUE), phyloseq::tax_table(UMB_taxamat), 
                            sample_data(p2017_1_samples))
library(ape)
p1_tree = rtree(ntaxa(p2017_1_filtered), rooted=TRUE, tip.label=taxa_names(p2017_1_filtered))
plot(p1_tree)
p2017_1_filtered<- phyloseq(otu_table(UMB_otu, taxa_are_rows = TRUE), phyloseq::tax_table(UMB_taxamat), 
                            sample_data(p2017_1_samples), p1_tree)
save(p2017_1_filtered, file = "~/p2017_1_filtered.Rdata")

#0.2um taxa matrix
large_taxID <- subset(total_tax_id, is.na(total_tax_id$`row.names(p2017_1_taxtable_df)`) )
row.names(large_taxID) <- large_taxID$`row.names(p2017_2_taxtable_df)`;large_taxID <- large_taxID[-1] ;head(large_taxID); nrow(large_taxID);ncol(large_taxID)#review output with new rownmanes
large_taxamat<- as.data.frame(merge(large_taxID, p2017_2_taxtable_df, by=0, all = TRUE))#add taxonomy information
large_taxamat <-large_taxamat[!is.na(large_taxamat$`row.names(p2017_2_taxtable_df)`),]
head(large_taxamat);ncol(large_taxamat)
row.names(large_taxamat) <- large_taxamat$`row.names(p2017_2_taxtable_df)`;large_taxamat <- as.matrix(large_taxamat[4:10]);head(large_taxamat); nrow(large_taxamat)
#0.2um otu dataframe
large_otu<- as.data.frame(merge(large_taxID, p2017_2_otutable_df, by=0, all = TRUE))#add taxonomy information
head(large_otu);ncol(large_otu)
large_otu <-large_otu[!is.na(large_otu$`row.names(p2017_2_taxtable_df)`),]
head(large_otu);ncol(large_otu)
row.names(large_otu) <- large_otu$`row.names(p2017_2_taxtable_df)`;large_otu <- as.data.frame(large_otu[4:47]);head(large_otu); nrow(large_otu)
#0.2 phyloseq
p2017_2_filtered<- phyloseq(otu_table(large_otu, taxa_are_rows = TRUE), phyloseq::tax_table(large_taxamat), 
                            sample_data(p2017_2_samples))
p2_tree = rtree(ntaxa(p2017_2_filtered), rooted=TRUE, tip.label=taxa_names(p2017_2_filtered))
plot(p2_tree)
p2017_2_filtered<- phyloseq(otu_table(large_otu, taxa_are_rows = TRUE), phyloseq::tax_table(large_taxamat), 
                            sample_data(p2017_2_samples),p2_tree)
p2017_2_filtered
save(p2017_2_filtered, file = "~/p2017_2_filtered.Rdata")

#merge0.1 and 0.2 to well day
sample_data(phylo2017)
merged_p2017 <- merge_samples(phylo2017, "WellDay")
#check sample data before saving the file
head(sample_data(merged_p2017)) #proceeed to 
save(merged_p2017, file = "~/merged_p2017.Rdata")

#0.1 taxa
umb_1_taxa <- as.matrix(dplyr::setdiff(p2017_1_taxtable_df,combined_intersecting_taxa))
umb_1_taxa_df <- as.data.frame(umb_1_taxa)
#0.2 taxa
large_taxa <- as.matrix(dplyr::setdiff(p2017_2_taxtable_df,combined_intersecting_taxa))
large_taxa_df <- as.data.frame(large_taxa)

#otu merge
combined_otu<-as.data.frame(merge(x=combined_intersecting_taxa, y=otu_2017,by=0,all.x=TRUE))
row.names(combined_otu) <- combined_otu$Row.names
combined_otu <- combined_otu[,9:91]; combined_otu
#0.1 OTU
umb_otu<-as.data.frame(merge(x=umb_1_taxa_df, y=p2017_1_otutable_df,by=0,all.x=TRUE))
row.names(umb_otu) <- umb_otu$Row.names
umb_otu <- umb_otu[,9:39] ;umb_otu
#0.2 OTU
large_otu<-as.data.frame(merge(x=large_taxa_df, y=p2017_2_otutable_df,by=0,all.x=TRUE))
row.names(large_otu) <- large_otu$Row.names
large_otu <- large_otu[,9:52] ;large_otu

#phyloseq objects for filtered OTUs
#0.2um
p2017_2_known_filtered <- phyloseq(otu_table(large_otu, taxa_are_rows = TRUE), phyloseq::tax_table(large_taxa), 
                                   sample_data(p2017_2_samples))

#0.1 um
p2017_1_known_filtered <- phyloseq(otu_table(umb_otu, taxa_are_rows = TRUE), phyloseq::tax_table(umb_1_taxa), 
                                   sample_data(p2017_1_samples))
#combined taxa
p2017_combined_known_filtered <- phyloseq(otu_table(combined_otu, taxa_are_rows = TRUE), phyloseq::tax_table(combined_intersecting_taxa_matrix), 
                                          sample_data(p2017_sample_data))


write.csv(filt_0.2_taxa, '~/2017_EVO_0.2_venn_taxa.csv')


#make rarefaction plots
ggarrange(venn_p2017, venn_cov, venn_rare, widths = c(3,3,3))
ggsave('~/2017_rarefied_taxa_venn.png', bg = "transparent", width = 40, height = 30, units = "cm")
dev.off()
plot.new()


# 2017 phyloseq object with all NA taxa removed
phylo2017_NA <- phyloseq::phyloseq(otu_table(otu_2017, taxa_are_rows = TRUE), phyloseq::tax_table(taxamat_NA), 
                                   sample_data(meta_complete), 
                                   tree_2017)
# rarefy the NA_phyloseq object. Removed 2198 OTUs
set.seed(0451)
phylo2017_NA_rarefied<-rarefy_even_depth(phylo2017_NA,sample.size = min(sample_sums(phylo2017_NA)))

#phyloseq merged by WellDay
sample_data(phylo2017)
merged_phylo2017 <- merge_samples(phylo2017, "WellDay")

# 2017 0.2 phyloseq object 
phylo2017_2 <- phyloseq::phyloseq(otu_table(otu_2017, taxa_are_rows = TRUE), phyloseq::tax_table(taxamat), 
                                  sample_data(geochem_2017_2), 
                                  tree_2017)
# rarefy the 0.2 um_phyloseq object. Removed 3116 OTUs
set.seed(0451)
phylo2017_2_rarefied<-rarefy_even_depth(phylo2017_2,sample.size = min(sample_sums(phylo2017_2)))

# 2017 0.1 phyloseq object 
phylo2017_1 <- phyloseq::phyloseq(otu_table(otu_2017, taxa_are_rows = TRUE), phyloseq::tax_table(taxamat), 
                                  sample_data(geochem_2017_1), 
                                  tree_2017)
phylo2017_1 <- phyloseq::phyloseq(otu_table(otu_2017, taxa_are_rows = TRUE), phyloseq::tax_table(taxamat), 
                                  sample_data(geochem_2017_1), 
                                  tree_2017)
phylo2017
phylo2017_1
# rarefy the 01 um phyloseq object. Removed 9498 OTUs
set.seed(0451)
phylo2017_1_rarefied<-rarefy_even_depth(phylo2017_1,sample.size = min(sample_sums(phylo2017_1)))
#Phyloseq Objects produced

filtered_0.1_taxa  #filtered o.1 um taxa from non-reduced phyloseq object _nrow:3950 (unique: 1023)
filtered_0.2_taxa  #filtered 0.2 um taxa from non-reduced phyloseq object_nrow:10967 (unique: 8040)
filtered_both_taxa #taxa overlap in both 0.2 and 0.1 from non-reduced phyloseq object_nrow:2927
phylo2017
phylo2017_rarefied
phylo2017_NA_rarefied
merged_phylo2017_rarefied
p2017_1_filtered
tax_1_filtered <- as.data.frame(phyloseq::tax_table(p2017_1_filtered))
tax_1_filtered  %>% count(Phylum, Class, Order, Family, sort = TRUE)
tax_comb_filtered <- as.data.frame(phyloseq::tax_table(p2017_combined_filtered))
tax_comb_filtered  %>% count(Phylum, Class,Order, sort = TRUE)
tax_2_filtered <- as.data.frame(phyloseq::tax_table(p2017_2_filtered))
tax_2_filtered  %>% count(Phylum, Class, sort = TRUE)
p2017_2_filtered
p2017_combined_filtered
p2017_1_top_filtered = prune_taxa(taxa_sums(p2017_1_filtered) > 20, p2017_1_filtered) 
p2017_2_top_filtered = prune_taxa(taxa_sums(p2017_2_filtered) > 20, p2017_2_filtered) 
p2017_combined_top_filtered = prune_taxa(taxa_sums(p2017_combined_filtered) > 20, p2017_combined_filtered) 
p2017_1_top_filtered
p2017_2_top_filtered
p2017_combined_top_filtered
ggarrange(venn_p2017, venn_cov, venn_rare, widths = c(3,3,3))
save(merged_phylo2017, file = "~/merged_phylo2017.Rdata")
save(phylo2017, file = "~/phylo2017.Rdata")
save(phylo2017_1, file = "~/phylo2017_1.Rdata")
save(phylo2017_2, file = "~/phylo2017_2.Rdata")

merged_p2017_sample_data <- as.data.frame((sample_data(merged_phylo2017_rarefied)))

#restructure sample data in merged phyloseq
#####
library(dplyr)
head(sample_data(merged_p2017))
head(sample_data(p2017_combined_merged))
sample_data(p2017_combined_merged)$DayCategory <- as.factor(sample_data(p2017_combined_merged)$Days)
sample_data(merged_p2017)$sample_label <- row.names(sample_data(merged_p2017))
sample_data(merged_p2017)$WellDay <- row.names(sample_data(merged_p2017))
sample_data(merged_p2017)$WellID <- gsub("\\Day.*", "", sample_data(merged_p2017)$WellDay)  
sample_data(merged_p2017)$DayCategory <- gsub(".*Day", "Day", sample_data(merged_p2017)$WellDay)  
sample_data(merged_p2017)$sample_label <-ifelse(sample_data(merged_p2017)$sample_label=="FW215","Control","Monitoring")
sample_data(p2017_combined_merged)$well_type_factor <-as.factor(ifelse(sample_data(p2017_combined_merged)$well_type=="Control",0,1))
sample_data(p2017_combined_merged)$Distance <-sample_data(p2017_combined_merged)$WellID
sample_data(p2017_combined_merged)$Distance <- dplyr::recode(sample_data(p2017_combined_merged)$Distance, "MLSB3" = 2.4, "FW216" = 2.5, "FW215" = -2.6, "GP01" = 4.9, "GP03" = 11.2)
sample_data(merged_p2017)$type_day <- paste(sample_data(merged_p2017)$well_type,sample_data(merged_p2017)$DayCategory)
sample_data(merged_p2017)$type_day <- sub("Control Day.*", "Control", sample_data(merged_p2017)$type_day)  
sample_data(p2017_combined_merged)$nitrate_reduction <-ifelse(sample_data(p2017_combined_merged)$Nitrate<=5,"Reduced","Not Reduced")
sample_data(p2017_combined_merged)$nitrate_reduction_factor <-as.factor(ifelse(sample_data(p2017_combined_merged)$nitrate_reduction=="Reduced",0,1))
sample_data(p2017_combined_merged)$sulfate_reduction <-ifelse(sample_data(p2017_combined_merged)$Sulfate<=50,"Reduced","Not Reduced")
sample_data(p2017_combined_merged)$sulfate_reduction_factor <-as.factor(ifelse(sample_data(p2017_combined_merged)$sulfate_reduction=="Reduced",0,1))
sample_data(p2017_combined_merged)$uranium_reduction <-ifelse(sample_data(p2017_combined_merged)$Uranium<=1,"Reduced","Not Reduced")
sample_data(p2017_combined_merged)$uranium_reduction_factor <-as.factor(ifelse(sample_data(p2017_combined_merged)$uranium_reduction=="Reduced",0,1))
sample_data(p2017_combined_merged)$phase <-dplyr::recode(sample_data(p2017_combined_merged)$sample_label, "FW215Day-6" = "Control",  "FW215Day1" = "Control",  "FW215Day106"="Control","FW215Day134"="Control","FW215Day15"="Control","FW215Day22"="Control","FW215Day50"="Control","FW215Day78"="Control", "FW215Day8"="Control",
                                                         "FW216Day-6"="Phase0",  "FW216Day1" = "Phase0", "FW216Day8" = "Phase2", "FW216Day15" = "Phase2","FW216Day22" = "Phase2", "FW216Day50" = "Phase3", "FW216Day78" = "Phase4","FW216Day106"= "Phase4","FW216Day134"= "Phase4",
                                                         "GP01Day-6" = "Phase0","GP01Day8"="Phase1", "GP01Day15"= "Phase2", "GP01Day22"="Phase2", "GP01Day50"="Phase3","GP01Day78"="Phase4","GP01Day106"="Phase4", "GP01Day134"="Phase4",
                                                         "GP03Day-6"="Phase0", "GP03Day1"="Phase0","GP03Day8"="Phase1", "GP03Day15" = "Phase2", "GP03Day22"="Phase2", "GP03Day50" = "Phase3", "GP03Day78"="Phase4","GP03Day106"="Phase4","GP03Day134"="Phase4",
                                                         "MLSB3Day-6"="Phase0", "MLSB3Day1"="Phase0", "MLSB3Day8"="Phase1", "MLSB3Day15"="Phase2", "MLSB3Day22"="Phase2", "MLSB3Day50"="Phase3", "MLSB3Day78" = "Phase4", "MLSB3Day106" = "Phase4", "MLSB3Day134"="Phase4")
sample_data(p2017_combined_merged)$phase_C134 <-dplyr::recode(sample_data(p2017_combined_merged)$sample_label, "FW215Day-6" = "Control",  "FW215Day1" = "Control",  "FW215Day106"="Control","FW215Day134"="Control","FW215Day15"="Control","FW215Day22"="Control","FW215Day50"="Control","FW215Day78"="Control", "FW215Day8"="Control",
                                                              "FW216Day-6"="Phase1",  "FW216Day1" = "Phase1", "FW216Day8" = "Phase3", "FW216Day15" = "Phase3","FW216Day22" = "Phase3", "FW216Day50" = "Phase3", "FW216Day78" = "Phase4","FW216Day106"= "Phase4","FW216Day134"= "Phase4",
                                                              "GP01Day-6" = "Phase1","GP01Day8"="Phase3", "GP01Day15"= "Phase3", "GP01Day22"="Phase3", "GP01Day50"="Phase3","GP01Day78"="Phase4","GP01Day106"="Phase4", "GP01Day134"="Phase4",
                                                              "GP03Day-6"="Phase1", "GP03Day1"="Phase1","GP03Day8"="Phase1", "GP03Day15" = "Phase3", "GP03Day22"="Phase3", "GP03Day50" = "Phase3", "GP03Day78"="Phase4","GP03Day106"="Phase4","GP03Day134"="Phase4",
                                                              "MLSB3Day-6"="Phase1", "MLSB3Day1"="Phase1", "MLSB3Day8"="Phase3", "MLSB3Day15"="Phase3", "MLSB3Day22"="Phase3", "MLSB3Day50"="Phase3", "MLSB3Day78" = "Phase4", "MLSB3Day106" = "Phase4", "MLSB3Day134"="Phase4")
sample_data(p2017_combined_merged)$phase_cat <- dplyr::recode(sample_data(p2017_combined_merged)$phase_C134, "Control" = 0, "Phase1" = 1, "Phase3" = 2, "Phase4" = 3)
#0.1
p2017_combined_0.1_filtered <-subset_samples(p2017_combined_filtered, FilterSize=="0.1um")
head(sample_data(p2017_combined_0.1_filtered))
p2017_combined_0.1_filtered <- prune_taxa(taxa_sums(p2017_combined_0.1_filtered) > 0, p2017_combined_0.1_filtered)
sum(otu_table(p2017_combined_0.1_filtered))
sample_data(p2017_combined_0.1_filtered)$acetate_presence <-as.factor(case_when((sample_data(p2017_combined_0.1_filtered)$Acetate <= 0.1)              ~ "≤0.1",
                                                                                (sample_data(p2017_combined_0.1_filtered)$Acetate   > 0.1) ~ ">0.1"))
save(p2017_combined_0.1_filtered, file = "~/p2017_combined_0.1_filtered.Rdata")

head(sample_data(p2017_1_filtered))
sample_data(merged_p2017)$phase_C134 <-dplyr::recode(sample_data(merged_p2017)$sample_label, "FW215Day-6" = "Control",  "FW215Day1" = "Control",  "FW215Day106"="Control","FW215Day134"="Control","FW215Day15"="Control","FW215Day22"="Control","FW215Day50"="Control","FW215Day78"="Control", "FW215Day8"="Control",
                                                     "FW216Day-6"="Phase1",  "FW216Day1" = "Phase1", "FW216Day8" = "Phase3", "FW216Day15" = "Phase3","FW216Day22" = "Phase3", "FW216Day50" = "Phase3", "FW216Day78" = "Phase4","FW216Day106"= "Phase4","FW216Day134"= "Phase4",
                                                     "GP01Day-6" = "Phase1","GP01Day8"="Phase3", "GP01Day15"= "Phase3", "GP01Day22"="Phase3", "GP01Day50"="Phase3","GP01Day78"="Phase4","GP01Day106"="Phase4", "GP01Day134"="Phase4",
                                                     "GP03Day-6"="Phase1", "GP03Day1"="Phase1","GP03Day8"="Phase1", "GP03Day15" = "Phase3", "GP03Day22"="Phase3", "GP03Day50" = "Phase3", "GP03Day78"="Phase4","GP03Day106"="Phase4","GP03Day134"="Phase4",
                                                     "MLSB3Day-6"="Phase1", "MLSB3Day1"="Phase1", "MLSB3Day8"="Phase3", "MLSB3Day15"="Phase3", "MLSB3Day22"="Phase3", "MLSB3Day50"="Phase3", "MLSB3Day78" = "Phase4", "MLSB3Day106" = "Phase4", "MLSB3Day134"="Phase4")
sample_data(merged_p2017)$phase_cat <- dplyr::recode(sample_data(merged_p2017)$phase_C134, "Control" = 0, "Phase1" = 1, "Phase3" = 2, "Phase4" = 3)

sample_data(p2017_1_filtered)$phase_C134 <-dplyr::recode(sample_data(p2017_1_filtered)$WellDay, "FW215Day-6" = "Control",  "FW215Day1" = "Control",  "FW215Day106"="Control","FW215Day134"="Control","FW215Day15"="Control","FW215Day22"="Control","FW215Day50"="Control","FW215Day78"="Control", "FW215Day8"="Control",
                                                         "FW216Day-6"="Phase1",  "FW216Day1" = "Phase1", "FW216Day8" = "Phase3", "FW216Day15" = "Phase3","FW216Day22" = "Phase3", "FW216Day50" = "Phase3", "FW216Day78" = "Phase4","FW216Day106"= "Phase4","FW216Day134"= "Phase4",
                                                         "GP01Day-6" = "Phase1","GP01Day8"="Phase3", "GP01Day15"= "Phase3", "GP01Day22"="Phase3", "GP01Day50"="Phase3","GP01Day78"="Phase4","GP01Day106"="Phase4", "GP01Day134"="Phase4",
                                                         "GP03Day-6"="Phase1", "GP03Day1"="Phase1","GP03Day8"="Phase1", "GP03Day15" = "Phase3", "GP03Day22"="Phase3", "GP03Day50" = "Phase3", "GP03Day78"="Phase4","GP03Day106"="Phase4","GP03Day134"="Phase4",
                                                         "MLSB3Day-6"="Phase1", "MLSB3Day1"="Phase1", "MLSB3Day8"="Phase3", "MLSB3Day15"="Phase3", "MLSB3Day22"="Phase3", "MLSB3Day50"="Phase3", "MLSB3Day78" = "Phase4", "MLSB3Day106" = "Phase4", "MLSB3Day134"="Phase4")
sample_data(p2017_1_filtered)$phase_cat <- dplyr::recode(sample_data(p2017_1_filtered)$phase_C134, "Control" = 0, "Phase1" = 1, "Phase3" = 2, "Phase4" = 3)


head(sample_data(p2017_combined_merged))
sample_data(p2017_1_filtered)$WellID <- gsub("\\Day.*", "", sample_data(p2017_1_filtered)$WellDay) 
sample_data(p2017_2_filtered)$WellID <- gsub("\\Day.*", "", sample_data(p2017_2_filtered)$WellDay) 
sample_data(p2017_combined_filtered)$WellID <- gsub("\\Day.*", "", sample_data(p2017_combined_filtered)$WellDay) 
sample_data(p2017_1_filtered)$well_type <-ifelse(sample_data(p2017_1_filtered)$WellID=="FW215","Control","Monitoring")
sample_data(p2017_1_filtered)$well_type_factor <-as.factor(ifelse(sample_data(p2017_1_filtered)$well_type=="Control",0,1))
sample_data(p2017_2_filtered)$well_type <-ifelse(sample_data(p2017_2_filtered)$WellID=="FW215","Control","Monitoring")
sample_data(p2017_2_filtered)$well_type_factor <-as.factor(ifelse(sample_data(p2017_2_filtered)$well_type=="Control",0,1))
sample_data(p2017_combined_filtered)$well_type <-ifelse(sample_data(p2017_combined_filtered)$WellID=="FW215","Control","Monitoring")
sample_data(p2017_combined_filtered)$well_type_factor <-as.factor(ifelse(sample_data(p2017_combined_filtered)$well_type=="Control",0,1))
head(sample_data(p2017_combined_filtered))
sample_data(p2017_combined_filtered)$phase_C134 <-dplyr::recode(sample_data(p2017_combined_filtered)$WellDay, "FW215Day-6" = "Control",  "FW215Day1" = "Control",  "FW215Day106"="Control","FW215Day134"="Control","FW215Day15"="Control","FW215Day22"="Control","FW215Day50"="Control","FW215Day78"="Control", "FW215Day8"="Control",
                                                                "FW216Day-6"="Phase1",  "FW216Day1" = "Phase1", "FW216Day8" = "Phase3", "FW216Day15" = "Phase3","FW216Day22" = "Phase3", "FW216Day50" = "Phase3", "FW216Day78" = "Phase4","FW216Day106"= "Phase4","FW216Day134"= "Phase4",
                                                                "GP01Day-6" = "Phase1","GP01Day8"="Phase3", "GP01Day15"= "Phase3", "GP01Day22"="Phase3", "GP01Day50"="Phase3","GP01Day78"="Phase4","GP01Day106"="Phase4", "GP01Day134"="Phase4",
                                                                "GP03Day-6"="Phase1", "GP03Day1"="Phase1","GP03Day8"="Phase1", "GP03Day15" = "Phase3", "GP03Day22"="Phase3", "GP03Day50" = "Phase3", "GP03Day78"="Phase4","GP03Day106"="Phase4","GP03Day134"="Phase4",
                                                                "MLSB3Day-6"="Phase1", "MLSB3Day1"="Phase1", "MLSB3Day8"="Phase3", "MLSB3Day15"="Phase3", "MLSB3Day22"="Phase3", "MLSB3Day50"="Phase3", "MLSB3Day78" = "Phase4", "MLSB3Day106" = "Phase4", "MLSB3Day134"="Phase4")


sample_data(p2017_combined_filtered)$phase_cat <- dplyr::recode(sample_data(p2017_combined_filtered)$phase_C134, "Control" = 0, "Phase1" = 1, "Phase3" = 2, "Phase4" = 3)

#assign injection phases

#give levels to geochem
sample_data(phylo2017)$nitrate_range <-case_when((sample_data(phylo2017)$NO3_mgL <= 1)              ~ 1,
                                                 (sample_data(phylo2017)$NO3_mgL  > 1) & (sample_data(phylo2017)$NO3_mgL <= 10) ~ 10,
                                                 (sample_data(phylo2017)$NO3_mgL  > 10) & (sample_data(phylo2017)$NO3_mgL <= 20) ~ 20,
                                                 (sample_data(phylo2017)$NO3_mgL  > 20) & (sample_data(phylo2017)$NO3_mgL <= 30) ~ 30)
sample_data(phylo2017)$uranium_range <-as.factor(case_when((sample_data(phylo2017)$U238_mgL <= 1)              ~ 1,
                                                           (sample_data(phylo2017)$U238_mgL   > 1) & (sample_data(phylo2017)$U238_mgL  <= 3) ~ 2))
sample_data(merged_p2017)$acetate_presence <-as.factor(case_when((sample_data(merged_p2017)$Acetate <= 20)              ~ "low",
                                                                 (sample_data(merged_p2017)$Acetate   > 20) ~ "high"))

sample_data(merged_p2017)$ammonium_presence <-as.factor(case_when((sample_data(merged_p2017)$Ammonium <= 0.4)              ~ "≤0.4",
                                                                  (sample_data(merged_p2017)$Ammonium   > 0.4) ~ ">0.4"))
sample_data(p2017_1)$ammonium_presence <-as.factor(case_when((sample_data(p2017_1)$NH4_mgL <= 0.5)              ~ "≤0.5",
                                                             (sample_data(p2017_1)$NH4_mgL   > 0.5) ~ ">0.5"))


View(sample_data(phylo2017))
head(sample_data(p2017_1))
sample_data(p2017_1)$acetate_presence <-as.factor(case_when((sample_data(p2017_1)$Acetate_uM <= 1)              ~ "low",
                                                            (sample_data(p2017_1)$Acetate_uM   > 1) ~ "high"))
sample_data(phylo2017)$sulfate_range <-case_when((sample_data(phylo2017)$SO4_mgL <= 20)              ~ 20,
                                                 (sample_data(phylo2017)$SO4_mgL  > 20) & (sample_data(phylo2017)$SO4_mgL <= 40) ~ 40,
                                                 (sample_data(phylo2017)$SO4_mgL  > 40) & (sample_data(phylo2017)$SO4_mgL <= 60) ~ 60,                                                                 (sample_data(phylo2017)$SO4_mgL  > 40) & (sample_data(phylo2017)$SO4_mgL <= 60) ~ 60,
                                                 (sample_data(phylo2017)$SO4_mgL  > 60) & (sample_data(phylo2017)$SO4_mgL <= 70) ~ 70)
#add type_day and levels to 0.1 phyloseq objects
sample_data(phylo2017_2)$type_day <- paste(sample_data(phylo2017_2)$well_type,sample_data(phylo2017_2)$DayCategory)
#give levels to geochem
sample_data(phylo2017_2)$nitrate_range <-case_when((sample_data(phylo2017_2)$NO3_mgL <= 1)              ~ 1,
                                                   (sample_data(phylo2017_2)$NO3_mgL  > 1) & (sample_data(phylo2017_2)$NO3_mgL <= 10) ~ 10,
                                                   (sample_data(phylo2017_2)$NO3_mgL  > 10) & (sample_data(phylo2017_2)$NO3_mgL <= 20) ~ 20,
                                                   (sample_data(phylo2017_2)$NO3_mgL  > 20) & (sample_data(phylo2017_2)$NO3_mgL <= 30) ~ 30)
sample_data(phylo2017_2)$uranium_range <-as.factor(case_when((sample_data(phylo2017_2)$U238_mgL <= 1)              ~ 1,
                                                             (sample_data(phylo2017_2)$U238_mgL   > 1) & (sample_data(phylo2017_2)$U238_mgL  <= 3) ~ 2))
sample_data(phylo2017_2)$sulfate_range <-case_when((sample_data(phylo2017_2)$SO4_mgL <= 20)              ~ 20,
                                                   (sample_data(phylo2017_2)$SO4_mgL  > 20) & (sample_data(phylo2017_2)$SO4_mgL <= 40) ~ 40,
                                                   (sample_data(phylo2017_2)$SO4_mgL  > 40) & (sample_data(phylo2017_2)$SO4_mgL <= 60) ~ 60,
                                                   (sample_data(phylo2017_2)$SO4_mgL  > 60) & (sample_data(phylo2017_2)$SO4_mgL <= 70) ~ 70)
#add type_day and levels to 0.2 phyloseq objects
sample_data(phylo2017_1)$type_day <- paste(sample_data(phylo2017_1)$well_type,sample_data(phylo2017_1)$DayCategory)
#give levels to geochem
sample_data(phylo2017_1)$nitrate_range <-case_when((sample_data(phylo2017_1)$NO3_mgL <= 1)              ~ 1,
                                                   (sample_data(phylo2017_1)$NO3_mgL  > 1) & (sample_data(phylo2017_1)$NO3_mgL <= 10) ~ 10,
                                                   (sample_data(phylo2017_1)$NO3_mgL  > 10) & (sample_data(phylo2017_1)$NO3_mgL <= 20) ~ 20,
                                                   (sample_data(phylo2017_1)$NO3_mgL  > 20) & (sample_data(phylo2017_1)$NO3_mgL <= 30) ~ 30)
sample_data(phylo2017_1)$uranium_range <-as.factor(case_when((sample_data(phylo2017_1)$U238_mgL <= 1)              ~ 1,
                                                             (sample_data(phylo2017_1)$U238_mgL   > 1) & (sample_data(phylo2017_1)$U238_mgL  <= 3) ~ 2))
sample_data(phylo2017_1)$sulfate_range <-case_when((sample_data(phylo2017_1)$SO4_mgL <= 20)              ~ 20,
                                                   (sample_data(phylo2017_1)$SO4_mgL  > 20) & (sample_data(phylo2017_1)$SO4_mgL <= 40) ~ 40,
                                                   (sample_data(phylo2017_1)$SO4_mgL  > 40) & (sample_data(phylo2017_1)$SO4_mgL <= 60) ~ 60,
                                                   (sample_data(phylo2017_1)$SO4_mgL  > 60) & (sample_data(phylo2017_1)$SO4_mgL <= 70) ~ 70)

#mean(meta_complete$NO3_mgL)
merged_p2017_sample_data <- sample_data(merged_phylo2017)
head(merged_p2017_sample_data)
#####

#abundance plots
#####
#0.1

prune_p1 <-phyloseq::prune_taxa(taxa_sums(p2017_1_filtered) > 7, p2017_1_filtered)
table(phyloseq::tax_table(prune_p1)[, "Phylum"])
ps1_phylum <- phyloseq::tax_glom(prune_p1, "Phylum")
phyloseq::taxa_names(ps1_phylum) <- phyloseq::tax_table(ps1_phylum)[, "Phylum"]
melt_p1_phylo2017 <- phyloseq::psmelt(prune_p1)
head(melt_p1_phylo2017)
ps1_samdat.df <- as.data.frame(sample_data(p2017_1_filtered))
head(ps1_samdat.df)
ps.1_phylum_merge <- merge(x=melt_p1_phylo2017, y=ps1_samdat.df,by.x= "Sample",
                           by.y =  0, all = TRUE)
ps.1_phylum_merge$WellID <- gsub("\\Day.*", "", ps.1_phylum_merge$WellDay.x)  
ps.1_phylum_merge

color_list_1 <- list(   "steelblue", "pink", "plum3", "orange3", "lightblue", "lightgoldenrod2" , "black","lightskyblue4" , "green2",
                        "maroon1", "mediumpurple", "lightcyan", "navajowhite4", "mediumpurple3", "blue",  "salmon",        
                        "tan2" ,"lightblue3", "peachpuff3" , "plum3",  "royalblue4" ,   "goldenrod3",
                        "turquoise4" ,  "seagreen3"  ,"sienna4","slategray", "red",  "green", 
                        "purple", "black", "salmon")
asv.count.plot.1 = ggplot(data = ps.1_phylum_merge, mapping = aes(x=WellID.x, y=Abundance, fill=Phylum)) +
  #geom_bar(stat = "identity", colour = "black") + 
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 10), legend.title = element_text(size = 10), 
        legend.text = element_text(colour = "black"), 
        axis.text.y = element_text(colour = "black")) + 
  scale_fill_manual(values = color_list_1) +
  scale_y_continuous(expand = c(0,0)) + 
  facet_grid(~factor(Days.x, levels=c('1','8','15', '22', '78', '106', '134'))) +
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = "bottom",
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        panel.border = element_rect(colour = "black", fill = "transparent", size=0.8))+
  labs(x = "Samples", y = "ASV Count", fill = "", title = "0.1um 2017 EVO")
asv.count.plot.1
ggsave("~/2017_facet_0.1count_phylum_plot.tiff", bg = "transparent", width = 30, height = 20, units = "cm")

#0.1 Abundance

prune_p1 <-phyloseq::prune_taxa(taxa_sums(p2017_1_filtered) > 7, p2017_1_filtered)
sample_data(prune_p1)
ps1 <- phyloseq::merge_samples(prune_p1, "WellDay")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
ps2_otu <- (phyloseq::otu_table(ps2)) ;ps2_sample <- (sample_data(ps2));ps2_tax <- (phyloseq::tax_table(ps2))
ps2_otu[is.na(ps2_otu)] = 0
ps2 <- phyloseq(otu_table(ps2_otu, taxa_are_rows = TRUE), phyloseq::tax_table(ps2_tax), 
                sample_data(ps2_sample))
table(phyloseq::tax_table(ps1)[, "Phylum"])
ps1_phylum <- phyloseq::tax_glom(ps1, "Phylum")
phyloseq::taxa_names(ps1_phylum) <- phyloseq::tax_table(ps1_phylum)[, "Phylum"]
melt_p1_phylo2017 <- phyloseq::psmelt(ps1)
head(melt_p1_phylo2017)
ps_samdat.df <- as.data.frame(sample_data(p2017_1_filtered))
head(ps_samdat.df)
ps.1_abun_merge <- merge(x=melt_p1_phylo2017, y=ps_samdat.df,
                         by.x= "Sample",
                         by.y =  "WellDay", all = TRUE)
ps.1_abun_merge
ps.1_abun_merge$WellID <- gsub("\\Day.*", "", ps.1_abun_merge$Sample)  
ps.1_abun_merge$Days.x <- gsub(".....\\Day*", "", ps.1_abun_merge$Sample)  
ps.1_abun_merge <- as.data.frame(ps.1_abun_merge)
head(ps.1_abun_merge)
asv.abund.plot.1 = ggplot(data = melt_p1_phylo2017, mapping = aes(x=WellID, y=Abundance, fill=Phylum)) +
  #geom_bar(stat = "identity", colour = "black") + 
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 10), legend.title = element_text(size = 10), 
        legend.text = element_text(colour = "black"), 
        axis.text.y = element_text(colour = "black")) + 
  scale_fill_manual(values = color_list_2) +
  scale_y_continuous(expand = c(0,0)) + 
  facet_grid(~factor(Days, levels=c('1','8','15', '22', '78', '106', '134'))) +
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = "bottom",
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        panel.border = element_rect(colour = "black", fill = "transparent", size=0.8))+
  labs(x = "Samples", y = "ASV Abundance (%)", fill = "", title = "0.1um 2017 EVO")
asv.abund.plot.1
ggsave("~/2017_facet_0.1count_phylum_plot.tiff", bg = "transparent", width = 30, height = 20, units = "cm")

#0.2
prune_p2 <-phyloseq::prune_taxa(taxa_sums(p2017_2_filtered) > 100, p2017_2_filtered)
table(phyloseq::tax_table(prune_p2)[, "Phylum"])
ps_phylum <- phyloseq::tax_glom(prune_p2, "Phylum")
phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
melt_p2_phylo2017 <- phyloseq::psmelt(prune_p2)
head(melt_p2_phylo2017)
ps2_samdat.df <- as.data.frame(sample_data(p2017_2_filtered))
head(ps2_samdat.df)
ps.2_phylum_merge <- merge(x=melt_p2_phylo2017, y=ps2_samdat.df,
                           by.x= "Sample",
                           by.y =  0, all = TRUE)
ps.2_phylum_merge
ps.2_phylum_merge$WellID <- gsub("\\Day.*", "", ps.2_phylum_merge$WellDay.x)  
ps.2_phylum_merge


color_list_2 <- list(   "steelblue", "pink", "plum3", "lightblue", "lightgoldenrod2" , "black","lightskyblue4" , "green2",
                        "maroon1", "mediumpurple", "yellow","lightcyan", "plum2","navajowhite4", "mediumpurple3", "blue", "purple3" ,"tan2", "salmon",        
                        "orchid3" ,"lightblue3", "peachpuff3" , "plum3", "red3", "royalblue4" ,  "paleturquoise2" ,  "goldenrod3",
                        "seagreen3"  ,"turquoise4" , "sienna4","slategray", "red",  "green", 
                        "purple", "black", "salmon")

asv.count.plot.2 = ggplot(data = ps.2_phylum_merge, mapping = aes(x=WellID, y=Abundance, fill=Phylum)) +
  #geom_bar(stat = "identity", colour = "black") + 
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 10), legend.title = element_text(size = 10), 
        legend.text = element_text(colour = "black"), 
        axis.text.y = element_text(colour = "black")) + 
  scale_fill_manual(values = color_list_2) +
  scale_y_continuous(expand = c(0,0)) + 
  facet_grid(~factor(Days.x, levels=c('-6','1','8','15','22', '50', '78', '106', '134'))) +
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = "bottom",
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        panel.border = element_rect(colour = "black", fill = "transparent", size=0.8))+
  labs(x = "Samples", y = "ASV Count", fill = "", title = "0.2um 2017 EVO")
asv.count.plot.2
ggsave("~/2017_facet_0.2count_phylum_plot.tiff", bg = "transparent", width = 30, height = 20, units = "cm")

#0.2 Abundance
prune_p2 <-phyloseq::prune_taxa(taxa_sums(p2017_2_filtered) > 200, p2017_2_filtered)
ps1 <- merge_samples(prune_p2, "WellDay")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
table(phyloseq::tax_table(ps2)[, "Phylum"])
ps_phylum <- phyloseq::tax_glom(ps2, "Phylum")
phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
melt_p2_phylo2017 <- phyloseq::psmelt(ps2)
ps_samdat.df <- as.data.frame(sample_data(p2017_2_filtered))
ps.2_abun_merge <- merge(x=melt_p2_phylo2017, y=ps_samdat.df,
                         by = 'row.names', all = TRUE)
ps.2_abun_merge
ps.2_abun_merge$WellID <- gsub("\\Day.*", "", ps.2_abun_merge$Sample)  
ps.2_abun_merge
dev.off()
asv.abund.plot.2 = ggplot(data = ps.2_abun_merge, mapping = aes(x=WellID, y=Abundance, fill=Phylum)) +
  #geom_bar(stat = "identity", colour = "black") + 
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 10), legend.title = element_text(size = 10), 
        legend.text = element_text(colour = "black"), 
        axis.text.y = element_text(colour = "black")) + 
  scale_fill_manual(values = color_list) +
  scale_y_continuous(expand = c(0,0)) + 
  facet_grid(~factor(Days.x, levels=c('-6','1','8','15','50', '22', '78', '106', '134'))) +
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = "bottom",
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        panel.border = element_rect(colour = "black", fill = "transparent", size=0.8))+
  labs(x = "Samples", y = "ASV Count", fill = "", title = "0.2um 2017 EVO")
asv.abund.plot.2
ggsave("~/2017_facet_0.2count_phylum_plot.tiff", bg = "transparent", width = 30, height = 20, units = "cm")

#small bacteria
prune_pc <-phyloseq::prune_taxa(taxa_sums(p2017_combined_filtered) > 1000, p2017_combined_filtered)
table(phyloseq::tax_table(prune_pc)[, "Phylum"])
psc_phylum <- phyloseq::tax_glom(prune_pc, "Phylum")
phyloseq::taxa_names(psc_phylum) <- phyloseq::tax_table(psc_phylum)[, "Phylum"]
melt_pcombined_phylo2017 <- phyloseq::psmelt(prune_pc)
head(melt_pcombined_phylo2017)
psc_samdat.df <- as.data.frame(sample_data(p2017_combined_filtered))
head(psc_samdat.df)
ps.c_phylum_merge <- as.data.frame(merge(x=melt_pcombined_phylo2017, y=psc_samdat.df, by.x= "Sample",by.y =  0, all = TRUE))
ps.c_phylum_merge
ps.c_phylum_merge$WellID <- gsub("\\Day.*", "", ps.c_phylum_merge$WellDay.x)  
head(ps.c_phylum_merge)
View(ps.c_phylum_merge)

color <- c(colors())
color_list_complete <- color[seq(379, length(color), 11)]
color_list_complete
asv.count.plot.count = ggplot(data = ps.c_phylum_merge, mapping = aes(x=WellID, y=Abundance, fill=Phylum)) +
  #geom_bar(stat = "identity", colour = "black") + 
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 10), legend.title = element_text(size = 10), 
        legend.text = element_text(colour = "black"), 
        axis.text.y = element_text(colour = "black")) + 
  scale_fill_manual(values = color_list) +
  scale_y_continuous(expand = c(0,0)) + 
  #facet_grid(FilterSize.x ~factor(Days.x, levels=c('-6','1','8','15','22', '50','78', '106', '134'))) +
  facet_grid(FilterSize.x ~Days.x) +
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = "bottom",
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        panel.border = element_rect(colour = "black", fill = "transparent", size=0.8))+
  labs(x = "", y = "ASV Count", fill = "", title = "Small Bacteria (0.1um and 0.2um overlap)")
asv.count.plot.count
ggsave("~/2017_facet_small_bacteria_phylum_plot.tiff", bg = "transparent", width = 30, height = 20, units = "cm")

#merged
prune_pc <-phyloseq::prune_taxa(taxa_sums(merged_p2017) > 1000, merged_p2017)
table(phyloseq::tax_table(prune_pc)[, "Phylum"])
psc_phylum <- phyloseq::tax_glom(prune_pc, "Phylum")
phyloseq::taxa_names(psc_phylum) <- phyloseq::tax_table(psc_phylum)[, "Phylum"]
melt_pcombined_phylo2017 <- phyloseq::psmelt(prune_pc)
head(melt_pcombined_phylo2017)
psc_samdat.df <- as.data.frame(sample_data(merged_p2017))
head(psc_samdat.df)
ps.c_phylum_merge <- as.data.frame(merge(x=melt_pcombined_phylo2017, y=psc_samdat.df, by.x= "Sample",by.y =  0, all = TRUE))
ps.c_phylum_merge
ps.c_phylum_merge$WellID <- gsub("\\Day.*", "", ps.c_phylum_merge$WellDay.x)  
head(ps.c_phylum_merge)
View(ps.c_phylum_merge)

color <- c(colors())
color_list <- color[seq(379, length(color), 11)]
asv.count.plot.count = ggplot(data = ps.c_phylum_merge, mapping = aes(x=WellID, y=Abundance, fill=Phylum)) +
  #geom_bar(stat = "identity", colour = "black") + 
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 10), legend.title = element_text(size = 10), 
        legend.text = element_text(colour = "black"), 
        axis.text.y = element_text(colour = "black")) + 
  scale_fill_manual(values = color_list) +
  scale_y_continuous(expand = c(0,0)) + 
  #facet_grid(FilterSize.x ~factor(Days.x, levels=c('-6','1','8','15','22', '50','78', '106', '134'))) +
  facet_grid(~Days.x) +
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = "bottom",
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        panel.border = element_rect(colour = "black", fill = "transparent", size=0.8))+
  labs(x = "", y = "ASV Count", fill = "", title = "Complete Community")
asv.count.plot.count
ggsave("~/2017_facet_complete_community_phylum_plot.tiff", bg = "transparent", width = 30, height = 20, units = "cm")


#end of phylum plots
#####

#Figure 2.5 - t-test boxplots of changes in abundance with geochem
#####
#assign new categories in phylo2017_rarefied
head(sample_data(phylo2017))
sample_data(phylo2017)$nitrate_level <-case_when((sample_data(phylo2017)$Nitrate_mgL <= 1)              ~ 1,
                                                 (sample_data(phylo2017)$Nitrate_mgL  > 1) & (sample_data(phylo2017)$Nitrate_mgL <= 10) ~ 10,
                                                 (sample_data(phylo2017)$Nitrate_mgL  > 10) & (sample_data(phylo2017)$Nitrate_mgL <= 20) ~ 20,
                                                 (sample_data(phylo2017)$Nitrate_mgL  > 20) & (sample_data(phylo2017)$Nitrate_mgL <= 30) ~ 30)
sample_data(phylo2017)$uranium_level <-as.factor(case_when((sample_data(phylo2017)$Uranium_mgL <= 1)              ~ 1,
                                                           (sample_data(phylo2017)$Uranium_mgL   > 1) & (sample_data(phylo2017)$Uranium_mgL  <= 3) ~ 2))
sample_data(phylo2017)$sulfate_level <-case_when((sample_data(phylo2017)$Sulfate_mgL <= 20)              ~ 20,
                                                 (sample_data(phylo2017)$Sulfate_mgL  > 20) & (sample_data(phylo2017)$Sulfate_mgL <= 40) ~ 40,
                                                 (sample_data(phylo2017)$Sulfate_mgL  > 40) & (sample_data(phylo2017)$Sulfate_mgL <= 60) ~ 60,
                                                 (sample_data(phylo2017)$Sulfate_mgL  > 60) & (sample_data(phylo2017)$Sulfate_mgL <= 70) ~ 70)

#select only most abundant phyla
prune_phylo2017 <-phyloseq::prune_taxa(taxa_sums(phylo2017_rarefied) > 2000, phylo2017_rarefied)
table(phyloseq::tax_table(prune_phylo2017)[, "Phylum"])
ps_phylum <- phyloseq::tax_glom(prune_phylo2017, "Phylum")
phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
melt_pr_phylo2017_rarefied <- phyloseq::psmelt(prune_phylo2017)
ps_phylum
ps_phylum.df <- as.data.frame(t(otu_table(ps_phylum)))
ps_phylum.df
colnames(ps_phylum.df)
ps_samdat.df <- as.data.frame(sample_data(phylo2017_rarefied))
ps_phylum_merge <- merge(x=ps_phylum.df, y=ps_samdat.df,
                         by = 'row.names', all = TRUE)
ps_phylum_merge

phyloseq::psmelt(ps_phylum)
compare_means(Patescibacteria~nitrate_level, data=ps_phylum_merge)
nitrate_comparisons = list(c("1","10"), c("1","20"), c("1","30"), c("10","20"), c("10","30"), c("20", "30"))
phyloseq::psmelt(ps_phylum) %>%
  ggplot(data = ., aes(x = as.character(nitrate_level), y = Abundance)) +
  #geom_boxplot(outlier.shape  = NA) +
  theme_bw()+
  geom_jitter(aes(color=OTU), height = 0, width = .2, size = 5, alpha =0.8) +
  labs(x = "Nitrate Levels (mg/L)", y = "Number of ASVs\n", color = "Phyla") +
  facet_grid(cols = vars(OTU), rows = vars(FilterSize),scales = "free") +
  #ggtitle("Changes in taxonomic abundance with nitrate content") +
  theme(axis.title = element_text(size = 45),
        legend.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 20),
        axis.text = element_text(size = 18))+ stat_compare_means(comparisons = nitrate_comparisons , method = "t.test", aes(label = ..p.signif..), ref.group = "0.5", hide.ns = TRUE, size = 12,step.increase = 0.09,label.y.npc = "bottom")+
  stat_summary(
    geom = "point",
    fun.y = "mean",
    col = "black",
    size = 7,
    shape = 24,
    fill = "black"
  )
ggsave("~/2017_nitrate_boxplot_facet.tiff", 
       bg = "transparent", width = 110, height = 70, units = "cm")

compare_means(Patescibacteria~sulfate_level, data=ps_phylum_merge)
sulfate_comparisons = list(c("70","20"),c("70","40"), c("70","60"), c("60","20"), c("60","40"))
phyloseq::psmelt(ps_phylum) %>%
  ggplot(data = ., aes(x = as.character(sulfate_level), y = Abundance)) +
  #geom_boxplot(outlier.shape  = NA) +
  theme_bw()+
  geom_jitter(aes(color=OTU), height = 0, width = .2, size = 5, alpha =0.8) +
  labs(x = "Sulfate Levels (mg/L)", y = "Number of ASVs\n", color = "Phyla") +
  facet_grid(cols = vars(OTU), rows = vars(FilterSize),scales = "free") +
  #ggtitle("Changes in taxonomic abundance with nitrate content") +
  theme(axis.title = element_text(size = 45),
        legend.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 20),
        axis.text = element_text(size = 18))+ stat_compare_means(comparisons = sulfate_comparisons, method = "t.test", aes(label = ..p.signif..), ref.group = "0.5" , hide.ns = TRUE ,size = 12,step.increase = 0.09,label.y.npc = "bottom")+
  stat_summary(
    geom = "point",
    fun.y = "mean",
    col = "black",
    size = 7,
    shape = 24,
    fill = "black"
  )

ggsave("~/2017_sulfate_boxplot_facet.tiff", 
       bg = "transparent", width = 110, height = 70, units = "cm")

#uranium

compare_means(Patescibacteria~uranium_level, data=ps_phylum_merge)
uranium_comparisons = list(c("1","2"))
phyloseq::psmelt(ps_phylum) %>%
  ggplot(data = ., aes(x = as.character(uranium_level), y = Abundance)) +
  #geom_boxplot(outlier.shape  = NA) +
  theme_bw()+
  geom_jitter(aes(color=OTU), height = 0, width = .2, size = 5, alpha =0.8) +
  labs(x = "Uranium Levels (mg/L)", y = "Number of ASVs\n", color = "Phyla") +
  facet_grid(cols = vars(OTU), rows = vars(FilterSize),scales = "free") +
  #ggtitle("Changes in taxonomic abundance with nitrate content") +
  theme(axis.title = element_text(size = 45),
        legend.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 20),
        axis.text = element_text(size = 18))+ stat_compare_means(comparisons = uranium_comparisons , method = "t.test", aes(label = ..p.signif..), ref.group = "0.5", hide.ns = TRUE, 
                                                                 size = 12,step.increase = 0.09,label.y.npc = "bottom")+
  stat_summary(
    geom = "point",
    fun.y = "mean",
    col = "black",
    size = 7,
    shape = 24,
    fill = "black"
  )
ggsave("~/2017_uranium_boxplot_facet.png", 
       bg = "transparent", width = 110, height = 70, units = "cm")

#end of t-test boxplots
#####

#t-test boxplots for bacteria in the 0.1 only 
#####
#subset most abundant phyla
prune_p1 <-phyloseq::prune_taxa(taxa_sums(p2017_1_filtered) > 7, p2017_1_filtered)
sum(otu_table(prune_p1))
table(phyloseq::tax_table(prune_p1)[, "Phylum"])
p1_phylum <- phyloseq::tax_glom(prune_p1, "Phylum")
phyloseq::taxa_names(p1_phylum) <- phyloseq::tax_table(p1_phylum)[, "Phylum"]
melt_p1_phylo2017 <- phyloseq::psmelt(prune_p1)
#export datatable as csv
View(melt_p1_phylo2017)
write.csv(melt_p1_phylo2017, "~/melt_p1_phylo2017.csv")

#continue preparing data for plots
melt.1_phylo2017_rarefied <- phyloseq::psmelt(prune_p1)
melt.1_phylo2017_rarefied
ps.1_order
ps.1_order.df <- as.data.frame(t(otu_table(ps.1_order)))
ps.1_order.df
colnames(ps.1_order.df)
ps.1_samdat.df <- as.data.frame(sample_data(p2017_1_filtered))
ps.1_samdat.df
ps.1_order_merge <- merge(x=ps.1_order.df, y=ps.1_samdat.df,
                          by = 'row.names', all = TRUE)
ps.1_order_merge

#0.1 phylum
table(phyloseq::tax_table(prune_p1)[, "Phylum"])
ps.1_Phylum <- phyloseq::tax_glom(prune_p1, "Phylum")
phyloseq::taxa_names(ps.1_Phylum) <- phyloseq::tax_table(ps.1_Phylum)[, "Phylum"]
melt.1_phylo2017_rarefied <- phyloseq::psmelt(prune_p1)
ps.1_Phylum
ps.1_Phylum.df <- as.data.frame(t(otu_table(ps.1_Phylum)))
ps.1_Phylum.df
colnames(ps.1_Phylum.df)
ps.1_samdat.df <- as.data.frame(sample_data(p2017_1_filtered))
ps.1_Phylum_merge <- merge(x=ps.1_Phylum.df, y=ps.1_samdat.df,
                           by = 'row.names', all = TRUE)
ps.1_Phylum_melt_merge <- merge(x=melt.1_phylo2017_rarefied, y=ps.1_samdat.df,
                                by = 'row.names', all = TRUE)
ps.1_Phylum_melt_merge
psmelt(p2017_1_filtered)
phyloseq::psmelt(p2017_1_filtered)
compare_means(Patescibacteria~nitrate_level, data=ps_phylum_merge)
nitrate_comparisons = list(c("1","10"), c("1","20"), c("1","30"), c("10","20"), c("10","30"), c("20", "30"))
phyloseq::psmelt(ps.1_Phylum) %>%
  ggplot(data = ., aes(x = as.character(nitrate_level), y = Abundance)) +
  #geom_boxplot(outlier.shape  = NA) +
  theme_bw()+
  geom_jitter(aes(color=Phylum), height = 0, width = .2, size = 5, alpha =0.8) +
  labs(x = "Nitrate Levels (mg/L)", y = "Number of ASVs\n", color = "Phyla") +
  facet_grid(cols = vars(Phylum),scales = "free") +
  #ggtitle("Changes in taxonomic abundance with nitrate content") +
  theme(axis.title = element_text(size = 45),
        legend.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 20),
        axis.text = element_text(size = 18))+ stat_compare_means(comparisons = nitrate_comparisons , method = "t.test", aes(label = ..p.signif..), ref.group = "0.5", hide.ns = TRUE, size = 12,step.increase = 0.09,label.y.npc = "bottom")+
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 7,
    shape = 24,
    fill = "black"
  )
ggsave("~/2017_nitrate_boxplot_0.1ONLY_facet.tiff", 
       bg = "transparent", width = 110, height = 70, units = "cm")

compare_means(Patescibacteria~sulfate_level, data=ps_phylum_merge)
sulfate_comparisons = list(c("70","20"),c("70","40"), c("70","60"), c("60","20"), c("60","40"))
phyloseq::psmelt(ps.1_Phylum) %>%
  ggplot(data = ., aes(x = as.character(sulfate_level), y = Abundance)) +
  #geom_boxplot(outlier.shape  = NA) +
  theme_bw()+
  geom_jitter(aes(color=Phylum), height = 0, width = .2, size = 5, alpha =0.8) +
  labs(x = "Sulfate Levels (mg/L)", y = "Number of ASVs\n", color = "Phyla") +
  facet_grid(cols = vars(Phylum),scales = "free") +
  #ggtitle("Changes in taxonomic abundance with nitrate content") +
  theme(axis.title = element_text(size = 45),
        legend.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 20),
        axis.text = element_text(size = 18))+ stat_compare_means(comparisons = sulfate_comparisons, method = "t.test", aes(label = ..p.signif..), ref.group = "0.5" , hide.ns = TRUE ,size = 12,step.increase = 0.09,label.y.npc = "bottom")+
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 7,
    shape = 24,
    fill = "black"
  )

ggsave("~/2017_sulfate_boxplot_0.1ONLYfacet.tiff", 
       bg = "transparent", width = 110, height = 70, units = "cm")

#uranium

compare_means(Patescibacteria~uranium_level, data=ps_phylum_merge)
uranium_comparisons = list(c("1","2"))
phyloseq::psmelt(ps.1_Phylum) %>%
  ggplot(data = ., aes(x = as.character(uranium_level), y = Abundance)) +
  #geom_boxplot(outlier.shape  = NA) +
  theme_bw()+
  geom_jitter(aes(color=Phylum), height = 0, width = .2, size = 5, alpha =0.8) +
  labs(x = "Uranium Levels (mg/L)", y = "Number of ASVs\n", color = "Phyla") +
  facet_grid(cols = vars(Phylum),scales = "free") +
  #ggtitle("Changes in taxonomic abundance with nitrate content") +
  theme(axis.title = element_text(size = 45),
        legend.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 20),
        axis.text = element_text(size = 18))+ stat_compare_means(comparisons = uranium_comparisons , method = "t.test", aes(label = ..p.signif..), ref.group = "0.5", hide.ns = TRUE, 
                                                                 size = 12,step.increase = 0.09,label.y.npc = "bottom")+
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 7,
    shape = 24,
    fill = "black"
  )
ggsave("~/2017_uranium_boxplot_0.1ONLYfacet.png", 
       bg = "transparent", width = 110, height = 70, units = "cm")

#end 0.1 taxa
#####

#t-test boxplots for bacteria in the 0.2 only 
#####
#subset most abundant phyla in 0.2
prune_p2 <-phyloseq::prune_taxa(taxa_sums(p2017_2_filtered) > 300, p2017_2_filtered)
table(phyloseq::tax_table(prune_p2)[, "Phylum"])
ps.2_Phylum <- phyloseq::tax_glom(prune_p2, "Phylum")
phyloseq::taxa_names(ps.2_Phylum) <- phyloseq::tax_table(ps.2_Phylum)[, "Phylum"]
melt.2_phylo2017_rarefied <- phyloseq::psmelt(prune_p2)
ps.2_Phylum
ps.2_Phylum.df <- as.data.frame(t(otu_table(ps.2_Phylum)))
ps.2_Phylum.df
colnames(ps.2_Phylum.df)
ps.2_samdat.df <- as.data.frame(sample_data(p2017_2_filtered))
ps.2_Phylum_merge <- merge(x=ps.2_Phylum.df, y=ps.2_samdat.df,
                           by = 'row.names', all = TRUE)
ps.2_Phylum_melt_merge <- merge(x=melt.2_phylo2017_rarefied, y=ps.2_samdat.df,
                                by = 'row.names', all = TRUE)
ps.2_Phylum_melt_merge

phyloseq::psmelt(p2017_2_filtered)
compare_means(Patescibacteria~nitrate_level, data=ps_phylum_merge)
nitrate_comparisons = list(c("1","10"), c("1","20"), c("1","30"), c("10","20"), c("10","30"), c("20", "30"))
phyloseq::psmelt(ps.2_Phylum) %>%
  ggplot(data = ., aes(x = as.character(nitrate_level), y = Abundance)) +
  #geom_boxplot(outlier.shape  = NA) +
  theme_bw()+
  geom_jitter(aes(color=Phylum), height = 0, width = .2, size = 5, alpha =0.8) +
  labs(x = "Nitrate Levels (mg/L)", y = "Number of ASVs\n", color = "Phyla") +
  facet_grid(cols = vars(Phylum),scales = "free") +
  #ggtitle("Changes in taxonomic abundance with nitrate content") +
  theme(axis.title = element_text(size = 45),
        legend.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 20),
        axis.text = element_text(size = 18))+ stat_compare_means(comparisons = nitrate_comparisons , method = "t.test", aes(label = ..p.signif..), ref.group = "0.5", hide.ns = TRUE, size = 12,step.increase = 0.09,label.y.npc = "bottom")+
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 7,
    shape = 24,
    fill = "black"
  )
ggsave("~/2017_nitrate_boxplot_0.2ONLY_facet.tiff", 
       bg = "transparent", width = 110, height = 70, units = "cm")

#Days
day_comparisons = list(c("-6","1"), c("-6","8"), c("-6","15"), c("-6","22"), c("-6","50"), c("-6","78"), c("-6","106"), c("-6","134"),
                       c("1","8"), c("1","15"), c("1","22"), c("1","50"), c("1","78"), c("1", "106"), c("1","134"), 
                       c("8","15"), c("8","22"), c("8","50"), c("8","78"), c("8", "106"), c("8","134"),
                       c("15","22"), c("15","50"), c("15","78"), c("15", "106"), c("15","134"), 
                       c("22","50"), c("22","78"), c("22", "106"), c("22","134"),
                       c("50","78"), c("50", "106"), c("50","134"),
                       c("78", "106"), c("78","134"),
                       c("106","134") )
psmelt(psc_phylum)
phyloseq::psmelt(ps.2_Phylum) %>%
  ggplot(data = ., aes(x = as.character(Days), y = Abundance)) +
  #geom_boxplot(outlier.shape  = NA) +
  theme_bw()+
  geom_jitter(aes(color=Phylum), height = 0, width = .2, size = 5, alpha =0.8) +
  labs(x = "Days", y = "Number of ASVs\n", color = "Phyla") +
  facet_grid(cols = vars(Phylum), scales = "free") +
  #ggtitle("Changes in taxonomic abundance with nitrate content") +
  theme(axis.title = element_text(size = 45),
        legend.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 20),
        axis.text = element_text(size = 18))+ 
  stat_compare_means(comparisons = day_comparisons , method = "t.test", aes(label = ..p.signif..), ref.group = "0.5", hide.ns = TRUE, size = 9)+
  
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 7,
    shape = 24,
    fill = "black"
  )
ggsave("~/2017_day_0.2ONLY_boxplot_combinedfacet.tiff", 
       bg = "transparent", width = 110, height = 70, units = "cm")


compare_means(Patescibacteria~sulfate_level, data=ps_phylum_merge)
sulfate_comparisons = list(c("70","20"),c("70","40"), c("70","60"), c("60","20"), c("60","40"))
phyloseq::psmelt(ps.2_Phylum) %>%
  ggplot(data = ., aes(x = as.character(sulfate_level), y = Abundance)) +
  #geom_boxplot(outlier.shape  = NA) +
  theme_bw()+
  geom_jitter(aes(color=Phylum), height = 0, width = .2, size = 5, alpha =0.8) +
  labs(x = "Sulfate Levels (mg/L)", y = "Number of ASVs\n", color = "Phyla") +
  facet_grid(cols = vars(Phylum),scales = "free") +
  #ggtitle("Changes in taxonomic abundance with nitrate content") +
  theme(axis.title = element_text(size = 45),
        legend.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 20),
        axis.text = element_text(size = 18))+ stat_compare_means(comparisons = sulfate_comparisons, method = "t.test", aes(label = ..p.signif..), ref.group = "0.5" , hide.ns = TRUE ,size = 12,step.increase = 0.09,label.y.npc = "bottom")+
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 7,
    shape = 24,
    fill = "black"
  )

ggsave("~/2017_sulfate_boxplot_0.2ONLYfacet.tiff", 
       bg = "transparent", width = 110, height = 70, units = "cm")

#uranium

compare_means(Patescibacteria~uranium_level, data=ps_phylum_merge)
uranium_comparisons = list(c("1","2"))
phyloseq::psmelt(ps.2_Phylum) %>%
  ggplot(data = ., aes(x = as.character(uranium_level), y = Abundance)) +
  #geom_boxplot(outlier.shape  = NA) +
  theme_bw()+
  geom_jitter(aes(color=Phylum), height = 0, width = .2, size = 5, alpha =0.8) +
  labs(x = "Uranium Levels (mg/L)", y = "Number of ASVs\n", color = "Phyla") +
  facet_grid(cols = vars(Phylum),scales = "free") +
  #ggtitle("Changes in taxonomic abundance with nitrate content") +
  theme(axis.title = element_text(size = 45),
        legend.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 20),
        axis.text = element_text(size = 18))+ stat_compare_means(comparisons = uranium_comparisons , method = "t.test", aes(label = ..p.signif..), ref.group = "0.5", hide.ns = TRUE, 
                                                                 size = 12,step.increase = 0.09,label.y.npc = "bottom")+
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 7,
    shape = 24,
    fill = "black"
  )
ggsave("~/2017_uranium_boxplot_0.2ONLYfacet.png", 
       bg = "transparent", width = 110, height = 70, units = "cm")

#end 0.2 taxa
#####

#t-test boxplots for combined 0.2 and 0.1 bacteria
#####

prune_pc <-phyloseq::prune_taxa(taxa_sums(p2017_combined_filtered) > 900, p2017_combined_filtered)
table(phyloseq::tax_table(prune_pc)[, "Phylum"])
psc_phylum <- phyloseq::tax_glom(prune_pc, "Phylum")
phyloseq::taxa_names(psc_phylum) <- phyloseq::tax_table(psc_phylum)[, "Phylum"]
melt_pcombined_phylo2017 <- phyloseq::psmelt(prune_pc)

View(melt_pcombined_phylo2017)
geochem_2017_1 <- subset(meta_complete, FilterSize == "0.1um")
write.csv(melt_pcombined_phylo2017, "~/melt_pcombined_phylo2017.csv")
psmerge_samdat.df <- as.data.frame(sample_data(merged_p2017))

head(melt_pcombined_phylo2017)
psc_samdat.df <- as.data.frame(sample_data(p2017_combined_filtered))
head(psc_samdat.df)
ps.c_phylum_merge <- as.data.frame(merge(x=melt_pcombined_phylo2017, y=psc_samdat.df, by.x= "Sample",by.y =  0, all = TRUE))
ps.c_phylum_merge
ps.c_phylum_merge$WellID <- gsub("\\Day.*", "", ps.c_phylum_merge$WellDay.x)  
head(ps.c_phylum_merge)
View(ps.c_phylum_merge)

color <- c(colors())
color_list <- color[seq(379, length(color), 11)]

#select only most abundant phyla
p2017_combined_0.1_filtered <-subset_samples(p2017_combined_filtered, FilterSize=="0.1um")
p2017_combined_filtered
p2017_combined_0.1_filtered
p2017_combined_0.1_filtered <- prune_taxa(taxa_sums(p2017_combined_0.1_filtered) > 0, p2017_combined_0.1_filtered)
sum(otu_table(p2017_combined_0.1_filtered))

prune_p1combined <-phyloseq::prune_taxa(taxa_sums(p2017_combined_0.1_filtered) > 350, p2017_combined_0.1_filtered)
table(phyloseq::tax_table(prune_p1combined)[, "Phylum"])
#table(phyloseq::tax_table(p2017_combined_top_filtered)[, "Phylum"])
psc1_phylum <- phyloseq::tax_glom(prune_p1combined, "Phylum")
phyloseq::taxa_names(psc1_phylum) <- phyloseq::tax_table(psc1_phylum)[, "Phylum"]
melt_pc1_phylo2017 <- phyloseq::psmelt(prune_p1combined)
melt_pc1_phylo2017
write.csv(melt_pc1_phylo2017, "~/melt_0.1_combined_asvs.csv")


psc_phylum
psc_phylum.df <- as.data.frame(t(otu_table(psc_phylum)))
psc_phylum.df
colnames(psc_phylum.df)
psc_samdat.df <- as.data.frame(sample_data(p2017_combined_filtered))
head(melt_pc_phylo2017_rarefied)
psc_phylum_merge <- merge(x=melt_pc_phylo2017_rarefied, y=psc_samdat.df, 
                          by.x= "Sample",by.y =  0, all = TRUE)

ps_phylum_merge
psmelt(psc_phylum)
phyloseq::psmelt(psc_phylum)
compare_means(Patescibacteria~nitrate_level, data=ps_phylum_merge)
nitrate_comparisons = list(c("1","10"), c("1","20"), c("1","30"), c("10","20"), c("10","30"), c("20", "30"))
phyloseq::psmelt(psc_phylum) %>%
  ggplot(data = ., aes(x = as.character(nitrate_level), y = Abundance)) +
  #geom_boxplot(outlier.shape  = NA) +
  theme_bw()+
  geom_jitter(aes(color=OTU), height = 0, width = .2, size = 5, alpha =0.8) +
  labs(x = "Nitrate Levels (mg/L)", y = "Number of ASVs\n", color = "Phyla") +
  facet_grid(cols = vars(OTU), rows = vars(FilterSize),scales = "free") +
  #ggtitle("Changes in taxonomic abundance with nitrate content") +
  theme(axis.title = element_text(size = 45),
        legend.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 20),
        axis.text = element_text(size = 18)) + stat_compare_means(comparisons = nitrate_comparisons , method = "t.test", aes(label = ..p.signif..), ref.group = "0.5", hide.ns = TRUE, size = 9)+
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 7,
    shape = 24,
    fill = "black"
  )
ggsave("~/2017_nitrate_Both_boxplot_facet.tiff", 
       bg = "transparent", width = 110, height = 70, units = "cm")

#Days
day_comparisons = list(c("-6","1"), c("-6","8"), c("-6","15"), c("-6","22"), c("-6","50"), c("-6","78"), c("-6","106"), c("-6","134"),
                       c("1","8"), c("1","15"), c("1","22"), c("1","50"), c("1","78"), c("1", "106"), c("1","134"), 
                       c("8","15"), c("8","22"), c("8","50"), c("8","78"), c("8", "106"), c("8","134"),
                       c("15","22"), c("15","50"), c("15","78"), c("15", "106"), c("15","134"), 
                       c("22","50"), c("22","78"), c("22", "106"), c("22","134"),
                       c("50","78"), c("50", "106"), c("50","134"),
                       c("78", "106"), c("78","134"),
                       c("106","134") )
psmelt(psc_phylum)
phyloseq::psmelt(psc_phylum) %>%
  ggplot(data = ., aes(x = as.character(Days), y = Abundance)) +
  #geom_boxplot(outlier.shape  = NA) +
  theme_bw()+
  geom_jitter(aes(color=Phylum), height = 0, width = .2, size = 5, alpha =0.8) +
  labs(x = "Days", y = "Number of ASVs\n", color = "Phyla") +
  facet_grid(cols = vars(Phylum), rows = vars(FilterSize),scales = "free") +
  #ggtitle("Changes in taxonomic abundance with nitrate content") +
  theme(axis.title = element_text(size = 45),
        legend.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 20),
        axis.text = element_text(size = 18))+ 
  #stat_compare_means(comparisons = day_comparisons , method = "t.test", aes(label = ..p.signif..), ref.group = "0.5", hide.ns = TRUE, size = 12,step.increase = 0.09,label.y.npc = "bottom")+
  
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 7,
    shape = 24,
    fill = "black"
  )
ggsave("~/2017_day_Both_boxplot_combinedfacet.tiff", 
       bg = "transparent", width = 110, height = 70, units = "cm")


compare_means(Patescibacteria~sulfate_level, data=ps_phylum_merge)
sulfate_comparisons = list(c("70","20"),c("70","40"), c("70","60"), c("60","20"), c("60","40"))
phyloseq::psmelt(psc_phylum) %>%
  ggplot(data = ., aes(x = as.character(sulfate_level), y = Abundance)) +
  #geom_boxplot(outlier.shape  = NA) +
  theme_bw()+
  geom_jitter(aes(color=OTU), height = 0, width = .2, size = 5, alpha =0.8) +
  labs(x = "Sulfate Levels (mg/L)", y = "Number of ASVs\n", color = "Phyla") +
  facet_grid(cols = vars(OTU), rows = vars(FilterSize),scales = "free") +
  #ggtitle("Changes in taxonomic abundance with nitrate content") +
  theme(axis.title = element_text(size = 45),
        legend.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 20),
        axis.text = element_text(size = 18))+ stat_compare_means(comparisons = sulfate_comparisons, method = "t.test", aes(label = ..p.signif..), ref.group = "0.5" , hide.ns = TRUE ,size = 12,step.increase = 0.09,label.y.npc = "bottom")+
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 7,
    shape = 24,
    fill = "black"
  )

ggsave("~/2017_sulfate_Both_boxplot_facet.tiff", 
       bg = "transparent", width = 110, height = 70, units = "cm")

#uranium

compare_means(Patescibacteria~uranium_level, data=ps_phylum_merge)
uranium_comparisons = list(c("1","2"))
phyloseq::psmelt(ps_phylum) %>%
  ggplot(data = ., aes(x = as.character(uranium_level), y = Abundance)) +
  #geom_boxplot(outlier.shape  = NA) +
  theme_bw()+
  geom_jitter(aes(color=OTU), height = 0, width = .2, size = 5, alpha =0.8) +
  labs(x = "Uranium Levels (mg/L)", y = "Number of ASVs\n", color = "Phyla") +
  facet_grid(cols = vars(OTU), rows = vars(FilterSize),scales = "free") +
  #ggtitle("Changes in taxonomic abundance with nitrate content") +
  theme(axis.title = element_text(size = 45),
        legend.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 20),
        axis.text = element_text(size = 18))+ stat_compare_means(comparisons = uranium_comparisons , method = "t.test", aes(label = ..p.signif..), ref.group = "0.5", hide.ns = TRUE, 
                                                                 size = 12,step.increase = 0.09,label.y.npc = "bottom")+
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 7,
    shape = 24,
    fill = "black"
  )
ggsave("~/2017_uranium_Both_boxplot_facet.png", 
       bg = "transparent", width = 110, height = 70, units = "cm")

#end t-test plots for the both-filter size-fraction
#####

#t-test boxplots for ccomplete community
#####

prune_pmerge <-phyloseq::prune_taxa(taxa_sums(merged_p2017) > 2000, merged_p2017)
table(phyloseq::tax_table(prune_pmerge)[, "Phylum"])
psmerge_phylum <- phyloseq::tax_glom(prune_pmerge, "Phylum")
phyloseq::taxa_names(psmerge_phylum) <- phyloseq::tax_table(psmerge_phylum)[, "Phylum"]
melt_pmerge2017 <- phyloseq::psmelt(prune_pmerge)

View(melt_pmerge2017)
write.csv(melt_pmerge2017, "~/melt_pmerge2017.csv")
psmerge_samdat.df <- as.data.frame(sample_data(merged_p2017))
head(psmerge_samdat.df)
ps.merge_phylum_merge <- as.data.frame(merge(x=melt_pmerge2017, y=psmerge_samdat.df, by.x= "Sample",by.y =  0, all = TRUE))
ps.merge_phylum_merge
ps.merge_phylum_merge$WellID <- gsub("\\Day.*", "", ps.merge_phylum_merge$WellDay.x)  
head(ps.merge_phylum_merge)
View(ps.merge_phylum_merge)

color <- c(colors())
color_list <- color[seq(379, length(color), 11)]

ps_phylum_merge
psmelt(psmerge_phylum)
phyloseq::psmelt(psmerge_phylum)
compare_means(Patescibacteria~nitrate_level.y, data=ps.merge_phylum_merge)
nitrate_comparisons = list(c("1","10"), c("1","20"), c("1","30"), c("10","20"), c("10","30"), c("20", "30"))
phyloseq::psmelt(psmerge_phylum) %>%
  ggplot(data = ., aes(x = as.character(nitrate_level), y = Abundance)) +
  #geom_boxplot(outlier.shape  = NA) +
  theme_bw()+
  geom_jitter(aes(color=OTU), height = 0, width = .2, size = 5, alpha =0.8) +
  labs(x = "Nitrate Levels (mg/L)", y = "Number of ASVs\n", color = "Phyla") +
  facet_grid(cols = vars(OTU),scales = "free") +
  #ggtitle("Changes in taxonomic abundance with nitrate content") +
  theme(axis.title = element_text(size = 35),
        legend.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20),
        axis.text = element_text(size = 18)) + stat_compare_means(comparisons = nitrate_comparisons , method = "t.test", aes(label = ..p.signif..), ref.group = "0.5", hide.ns = TRUE, size = 9)+
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 7,
    shape = 24,
    fill = "black"
  )
ggsave("~/2017_nitrate_Complete_microbial_community_facet.tiff", 
       bg = "transparent", width = 150, height = 80, units = "cm", limitsize=FALSE)

#Days
day_comparisons = list(c("-6","1"), c("-6","8"), c("-6","15"), c("-6","22"), c("-6","50"), c("-6","78"), c("-6","106"), c("-6","134"),
                       c("1","8"), c("1","15"), c("1","22"), c("1","50"), c("1","78"), c("1", "106"), c("1","134"), 
                       c("8","15"), c("8","22"), c("8","50"), c("8","78"), c("8", "106"), c("8","134"),
                       c("15","22"), c("15","50"), c("15","78"), c("15", "106"), c("15","134"), 
                       c("22","50"), c("22","78"), c("22", "106"), c("22","134"),
                       c("50","78"), c("50", "106"), c("50","134"),
                       c("78", "106"), c("78","134"),
                       c("106","134") )
psmelt(psc_phylum)
phyloseq::psmelt(psmerge_phylum) %>%
  ggplot(data = ., aes(x = as.character(Days), y = Abundance)) +
  #geom_boxplot(outlier.shape  = NA) +
  theme_bw()+
  geom_jitter(aes(color=Phylum), height = 0, width = .2, size = 5, alpha =0.8) +
  geom_smooth()+
  labs(x = "Days", y = "Number of ASVs\n", color = "Phyla") +
  facet_grid(cols = vars(OTU),scales = "free") +
  #ggtitle("Changes in taxonomic abundance with nitrate content") +
  theme(axis.title = element_text(size = 45),
        legend.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 20),
        axis.text = element_text(size = 18))+ 
  #stat_compare_means(comparisons = day_comparisons , method = "t.test", aes(label = ..p.signif..), ref.group = "0.5", hide.ns = TRUE, size = 12,step.increase = 0.09,label.y.npc = "bottom")+
  
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 7,
    shape = 24,
    fill = "black"
  )
ggsave("~/2017_day_Complete_microbial_community.tiff", 
       bg = "transparent", width = 150, height = 80, units = "cm", limitsize=FALSE)


compare_means(Patescibacteria~sulfate_level, data=ps_phylum_merge)
sulfate_comparisons = list(c("70","20"),c("70","40"), c("70","60"), c("60","20"), c("60","40"))
phyloseq::psmelt(psmerge_phylum) %>%
  ggplot(data = ., aes(x = as.character(sulfate_level), y = Abundance)) +
  #geom_boxplot(outlier.shape  = NA) +
  theme_bw()+
  geom_jitter(aes(color=OTU), height = 0, width = .2, size = 5, alpha =0.8) +
  labs(x = "Sulfate Levels (mg/L)", y = "Number of ASVs\n", color = "Phyla") +
  facet_grid(cols = vars(OTU),scales = "free") +
  #ggtitle("Changes in taxonomic abundance with nitrate content") +
  theme(axis.title = element_text(size = 35),
        legend.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20),
        axis.text = element_text(size = 18))+ stat_compare_means(comparisons = sulfate_comparisons, method = "t.test", aes(label = ..p.signif..), ref.group = "0.5" , hide.ns = TRUE ,size = 12,step.increase = 0.09,label.y.npc = "bottom")+
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 7,
    shape = 24,
    fill = "black"
  )

ggsave("~/2017_sulfate_Complete_microbial_community.tiff", 
       bg = "transparent", width = 110, height = 70, units = "cm")

#uranium

compare_means(Patescibacteria~uranium_level, data=ps_phylum_merge)
uranium_comparisons = list(c("1","2"))
phyloseq::psmelt(psmerge_phylum) %>%
  ggplot(data = ., aes(x = as.character(uranium_level), y = Abundance)) +
  #geom_boxplot(outlier.shape  = NA) +
  theme_bw()+
  geom_jitter(aes(color=OTU), height = 0, width = .2, size = 5, alpha =0.8) +
  labs(x = "Uranium Levels (mg/L)", y = "Number of ASVs\n", color = "Phyla") +
  facet_grid(cols = vars(OTU), scales = "free") +
  #ggtitle("Changes in taxonomic abundance with nitrate content") +
  theme(axis.title = element_text(size = 45),
        legend.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 20),
        axis.text = element_text(size = 18))+ stat_compare_means(comparisons = uranium_comparisons , method = "t.test", aes(label = ..p.signif..), ref.group = "0.5", hide.ns = TRUE, 
                                                                 size = 12,step.increase = 0.09,label.y.npc = "bottom")+
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 7,
    shape = 24,
    fill = "black"
  )
ggsave("~/2017_uranium_Complete_microbial_community.tiff", 
       bg = "transparent", width = 110, height = 70, units = "cm")

#end t-test plots for the complete microbial community
#####


# Figure 2.6 - taxa and geochem correlations
#####
# taxa correlation plots
# create a couple of numerical variables to use
colnames(sample_data(phylo2017))
colnames(sample_data(p2017_1_filtered))
phylo2017_NA_rarefied
#rename if not named as desired
colnames(sample_data(phylo2017)) <- c('WellID','Days','FilterSize','EVOYear',
                                      "pH",'DissolvedOxygen','SpecificConductivity_uScm','WaterTable_mAMSL',
                                      'Magnesium', 'Aluminum', 'Potassium', "Calcium",
                                      'Iron', 'Manganese', 'Uranium', 'Nitrate',
                                      'Sulfate', 'Acetate', 'Ammonium','well_type',
                                      'Oxygen_content','nitrate_reduction','DayCategory','sample_label',
                                      "sulfate_reduction", "phase", "WellDay")

#complete microbial community
psq_comb_2017 <- p2017_combined_filtered %>%
  ps_mutate(Well = if_else(well_type == "Monitoring", true = 1, false = 0),
            Filter = if_else(FilterSize == "0.1um", true = 1, false = 2))

psq_comb_2017 <- tax_filter(psq_comb_2017, min_prevalence = 1 / 10, min_sample_abundance = 1 / 10)
psq_comb_2017 <- phyloseq::tax_glom(psq_comb_2017, "Genus", NArm = TRUE)
psq_comb_2017 <- tax_agg(psq_comb_2017, "Genus")

# randomly select 20 taxa from the 50 most abundant taxa
set.seed(0451)
taxa <- sample(tax_top(psq_comb_2017, n = 20), size = 20)

# NOTE: detection threshold set to 50 as HITchip example data seems to have background noise
ud <- 20

# correlations and multiple annotations
jpeg("~/2017_taxa_genus_notreduced_cor_heatmap.tiff",width=2200, height=2800, res=400, bg="white")
plot.new()
cor_heatmap(
  data = psq_comb_2017, taxa = taxa[1:20], 
  vars = c("Sulfate", "Nitrate","Uranium","Iron","Acetate", "Ammonium"), 
  cor = "spearman",
  colors = heat_palette("Tropic", rev = FALSE, sym = TRUE),
  tax_anno = taxAnnotation(
    Abundance = anno_tax_box(undetected = ud)
  ),
  var_anno = varAnnotation(
    ppm = anno_var_hist(size = grid::unit(20, "mm"))
    #Log10 = anno_var_box(function(x) log10(x + 1))
  )
)

getOption("device")
dev.set(which = dev.next())
dev.off()

#0.1um community
#rename if not named as desired
colnames(sample_data(phylo2017_1))
colnames(sample_data(phylo2017_1)) <- c('WellID','Days','FilterSize','EVOYear',"pH",
                                        'DissolvedOxygen','SpecificConductivity_uScm','WaterTable_mAMSL','Magnesium', 'Aluminum',
                                        'Potassium', "Calcium",'Iron', 'Manganese', 'Uranium',
                                        'Nitrate', 'Sulfate', 'Acetate', 'Ammonium','well_type',
                                        'Oxygen_content','nitrate_reduction','DayCategory',"Days",'sample_label',"sulfate_reduction", 
                                        "phase", "type_day","nitrate_range", "uranium_range", "sulfate_range")
#0.1 OVERLAP microbial community
psq_1_2017 <- p2017_combined_0.1_filtered
#filter genera by abundance
psq_1_2017 <- tax_filter(psq_1_2017, min_prevalence = 1 / 10, min_sample_abundance = 1 / 10)
psq_1_2017 <- phyloseq::tax_glom(psq_1_2017, "Genus", NArm = TRUE)
psq_1_2017 <- tax_agg(psq_1_2017, "Genus")
# randomly select 20 taxa from the 50 most abundant taxa
set.seed(0451)
taxa1 <- sample(tax_top(psq_1_2017, n = 20), size = 20)
taxa1
# NOTE: detection threshold set to 50 as HITchip example data seems to have background noise
ud <- 10
# correlations and multiple annotations
jpeg("~/2017_taxa_0.1_combined_genus_cor_heatmap.jpeg",width=2200, height=2800, res=400, bg="white")
plot.new()
cor_heatmap(
  data = psq_1_2017, taxa = taxa1, 
  vars = c("Sulfate", "Nitrate","Uranium","Iron","Acetate", "Ammonium"), 
  cor = "spearman",
  colors = heat_palette("Tropic", rev = FALSE, sym = TRUE),
  tax_anno = taxAnnotation(Abundance = anno_tax_box(undetected = ud)),
  var_anno = varAnnotation(
    ppm = anno_var_hist(size = grid::unit(20, "mm"))
    #Log10 = anno_var_box(function(x) log10(x + 1))
  )
)

getOption("device")
dev.set(which = dev.next())
dev.off()

#order-level plot
#0.1 OVERLAP microbial community
psq_1ord_2017 <- p2017_combined_0.1_filtered
#filter genera by abundance
psq_1ord_2017 <- tax_filter(psq_1ord_2017, min_prevalence = 1 / 10, min_sample_abundance = 1 / 10)
psq_1ord_2017 <- phyloseq::tax_glom(psq_1ord_2017, "Order", NArm = TRUE)
psq_1ord_2017 <- tax_agg(psq_1ord_2017, "Order")
# randomly select 20 taxa from the 50 most abundant taxa
set.seed(0451)
taxa1ord <- sample(tax_top(psq_1ord_2017, n = 20), size = 20)
taxa1ord
# NOTE: detection threshold set to 50 as HITchip example data seems to have background noise
ud <- 10

# correlations and multiple annotations
jpeg("~/2017_taxa_0.1_combined_order_cor_heatmap.jpeg",width=2200, height=2800, res=400, bg="white")
plot.new()
cor_heatmap(
  data = psq_1ord_2017, taxa = taxa1ord, 
  vars = c("Sulfate", "Nitrate","Uranium","Iron","Acetate", "Ammonium"), 
  cor = "spearman",
  colors = heat_palette("Tropic", rev = FALSE, sym = TRUE),
  tax_anno = taxAnnotation(Abundance = anno_tax_box(undetected = ud)),
  var_anno = varAnnotation(
    ppm = anno_var_hist(size = grid::unit(20, "mm"))
    #Log10 = anno_var_box(function(x) log10(x + 1))
  )
)

getOption("device")
dev.set(which = dev.next())
dev.off()

#0.1 EXCLUSIVE microbial community
psq_1E_2017 <- p2017_1_filtered
#filter genera by abundance
psq_1E_2017 <- tax_filter(psq_1E_2017, min_prevalence = 1 / 100, min_sample_abundance = 1 / 100)
psq_1E_2017 <- phyloseq::tax_glom(psq_1E_2017, "Genus", NArm = TRUE)
psq_1E_2017 <- tax_agg(psq_1E_2017, "Genus")
psq_1E_2017
# randomly select 20 taxa from the 50 most abundant taxa
set.seed(0451)
taxa1E <- sample(tax_top(psq_1E_2017, n = 20), size = 20)
taxa1E
# NOTE: detection threshold set to 50 as HITchip example data seems to have background noise
ud <- 10
# correlations and multiple annotations
jpeg("~/2017_taxa_0.1_exclusive_genus_cor_heatmap.jpeg",width=2200, height=2800, res=400, bg="white")
plot.new()
cor_heatmap(
  data = psq_1E_2017, taxa = taxa1E, 
  vars = c("Sulfate", "Nitrate","Uranium","Iron","Acetate", "Ammonium"), 
  cor = "spearman",
  colors = heat_palette("Tropic", rev = FALSE, sym = TRUE),
  tax_anno = taxAnnotation(Abundance = anno_tax_box(undetected = ud)),
  var_anno = varAnnotation(
    ppm = anno_var_hist(size = grid::unit(20, "mm"))
    #Log10 = anno_var_box(function(x) log10(x + 1))
  )
)

getOption("device")
dev.set(which = dev.next())
dev.off()

#order-level plot
#0.1 EXCLUSIVE microbial community
psq_1Eord_2017 <- p2017_1_filtered
#filter genera by abundance
psq_1Eord_2017 <- tax_filter(psq_1Eord_2017, min_prevalence = 1 / 100, min_sample_abundance = 1 / 100)
psq_1Eord_2017 <- phyloseq::tax_glom(psq_1Eord_2017, "Order", NArm = TRUE)
psq_1Eord_2017 <- tax_agg(psq_1Eord_2017, "Order")
psq_1Eord_2017
# randomly select 20 taxa from the 50 most abundant taxa
set.seed(0451)
taxa1Eord <- sample(tax_top(psq_1Eord_2017, n = 20), size = 20)
taxa1Eord
# NOTE: detection threshold set to 50 as HITchip example data seems to have background noise
ud <- 10

# correlations and multiple annotations
jpeg("~/2017_taxa_0.1_exclusive_order_cor_heatmap.jpeg",width=2200, height=2800, res=400, bg="white")
plot.new()
cor_heatmap(
  data = psq_1Eord_2017, taxa = taxa1Eord, 
  vars = c("Sulfate", "Nitrate","Uranium","Iron","Acetate", "Ammonium"), 
  cor = "spearman",
  colors = heat_palette("Tropic", rev = FALSE, sym = TRUE),
  tax_anno = taxAnnotation(Abundance = anno_tax_box(undetected = ud)),
  var_anno = varAnnotation(
    ppm = anno_var_hist(size = grid::unit(20, "mm"))
    #Log10 = anno_var_box(function(x) log10(x + 1))
  )
)

getOption("device")
dev.set(which = dev.next())
dev.off()

#0.2um community
#rename col labels
colnames(sample_data(phylo2017_2_rarefied))
colnames(sample_data(phylo2017_2_rarefied)) <- c('WellID','Days','FilterSize','EVOYear',"pH",
                                                 'DissolvedOxygen','SpecificConductivity_uScm','WaterTable_mAMSL','Magnesium', 'Aluminum',
                                                 'Potassium', "Calcium",'Iron', 'Manganese', 'Uranium',
                                                 'Nitrate', 'Sulfate', 'Acetate', 'Ammonium','well_type',
                                                 'Oxygen_content','nitrate_reduction','DayCategory','sample_label',"sulfate_reduction", 
                                                 "phase", "type_day","nitrate_range", "uranium_range", "sulfate_range")

#complete microbial community
psq_2_2017 <- p2017_2_filtered %>%
  ps_mutate(Well = if_else(well_type == "Monitoring", true = 1, false = 0))
phyloseq::tax_table(p2017_2_filtered)
phyloseq::otu_table(p2017_2_filtered)
#filter abudannte genera
psq_2_2017 <- tax_filter(psq_2_2017, min_prevalence = 1 / 10, min_sample_abundance = 1 / 10)
psq_2_2017 <- phyloseq::tax_glom(psq_2_2017, "Genus", NArm = TRUE)
psq_2_2017 <- tax_agg(psq_2_2017, "Genus")

# randomly select 20 taxa from the 50 most abundant taxa
set.seed(0451)
taxa2 <- sample(tax_top(psq_2_2017, n = 20), size = 20)

# NOTE: detection threshold set to 50 as HITchip example data seems to have background noise
ud <- 20

# correlations and multiple annotations
jpeg("~/2017_taxa_0.2_genus_notreduced_cor_heatmap.tiff",width=2200, height=2800, res=400, bg="white")
plot.new()
cor_heatmap(
  data = psq_2_2017, taxa = taxa2[1:20], 
  vars = c("Sulfate", "Nitrate","Uranium","Iron","Acetate", "Ammonium"), 
  cor = "spearman",
  colors = heat_palette("Tropic", rev = FALSE, sym = TRUE),
  tax_anno = taxAnnotation(
    Abundance = anno_tax_box(undetected = ud)
  ),
  var_anno = varAnnotation(
    ppm = anno_var_hist(size = grid::unit(20, "mm"))
    #Log10_ppm = anno_var_box(function(x) log10(x + 1))
  )
)

getOption("device")
dev.set(which = dev.next())
dev.off()

#reset column names back
colnames(sample_data(phylo2017_NA_rarefied)) <- c('WellID','Days','FilterSize','EVOYear',
                                                  "pH",'DissolvedOxygen_mgL','SpecificConductivity_uScm','WaterTable_mAMSL',
                                                  'Magnesium_mgL', 'Aluminum_mgL', 'Potassium_mgL', "Calcium_mgL",
                                                  'Iron_mgL', 'Manganese_mgL', 'Uranium_mgL', 'Nitrate_mgL',
                                                  'Sulfate_mgL', 'Acetate_mgL', 'Ammonium_mgL','well_type',
                                                  'Oxygen_content','nitrate_reduction','DayCategory','sample_label',
                                                  "sulfate_reduction", "phase", "WellDay")
colnames(sample_data(phylo2017_1_rarefied)) <- c('WellID','Days','FilterSize','EVOYear',
                                                 "pH",'DissolvedOxygen_mgL','SpecificConductivity_uScm','WaterTable_mAMSL',
                                                 'Magnesium_mgL', 'Aluminum_mgL', 'Potassium_mgL', "Calcium_mgL",
                                                 'Iron_mgL', 'Manganese_mgL', 'Uranium_mgL', 'Nitrate_mgL',
                                                 'Sulfate_mgL', 'Acetate_mgL', 'Ammonium_mgL','well_type',
                                                 'Oxygen_content','nitrate_reduction','DayCategory','sample_label',
                                                 "sulfate_reduction", "phase", "type_day","nitrate_range", "uranium_range", "sulfate_range")
colnames(sample_data(phylo2017_2_rarefied)) <- c('WellID','Days','FilterSize','EVOYear',
                                                 "pH",'DissolvedOxygen_mgL','SpecificConductivity_uScm','WaterTable_mAMSL',
                                                 'Magnesium_mgL', 'Aluminum_mgL', 'Potassium_mgL', "Calcium_mgL",
                                                 'Iron_mgL', 'Manganese_mgL', 'Uranium_mgL', 'Nitrate_mgL',
                                                 'Sulfate_mgL', 'Acetate_mgL', 'Ammonium_mgL','well_type',
                                                 'Oxygen_content','nitrate_reduction','DayCategory','sample_label',
                                                 "sulfate_reduction", "phase", "type_day","nitrate_range", "uranium_range", "sulfate_range")
#####end figure 2.6

#taxa correlation plots
#merged phyla heatmap plots
#####

#taxonomic correlations
#####
#correlation analysis 
#extract genus-level ids
library(MicrobiotaProcess)
library(phyloseq)
#extract genera 
genustab <- MicrobiotaProcess::get_taxadf(merged_p2017, taxlevel=6)
genustab <- as.data.frame(t(phyloseq::otu_table(genustab)), check.names=FALSE)
genustab <- as.data.frame(apply(genustab, 2, function(x)x/sum(x)), check.names=FALSE)
View(genustab)
genustab  <- na.omit(t(genustab)) 
genustab <-t(genustab);View(genustab)

genustab_cor = cor(as.matrix(genustab), method = "spearman")
genustab_cor <- as.data.frame(genustab_cor)
View(genustab_cor)



install.packages("lares")
library(lares)

top_cross_corr<- corr_cross(genustab_cor, # name of dataset
                            plot= TRUE,
                            ignore = NA,
                            max_pvalue = 0.05, # max correlation value
                            top = 50 # display top 20 couples of variables (by correlation coefficient)
)
top_cross_corr
top_cross_corr.95<- corr_cross(genustab_cor, # name of dataset
                               plot= TRUE,
                               ignore = NA,
                               max=0.95, # max correlation value
                               top = 100 # display top 20 couples of variables (by correlation coefficient)
)
View(top_cross_corr.95$data)
cross_corr.99<- corr_cross(genustab_cor, # name of dataset
                           plot= TRUE,
                           ignore = NA,
                           max=0.99, # max correlation value
                           top = 1000 # display top 20 couples of variables (by correlation coefficient)
)
View(cross_corr.99$data)

#test correlations on families
#extract genera 
famtab <- MicrobiotaProcess::get_taxadf(merged_p2017, taxlevel=5)
famtab <- as.data.frame(t(phyloseq::otu_table(famtab)), check.names=FALSE)
famtab <- as.data.frame(apply(famtab, 2, function(x)x/sum(x)), check.names=FALSE)
View(famtab)
famtab  <- na.omit(t(famtab)) 
famtab <-t(famtab);View(famtab)

famtab_cor = cor(as.matrix(famtab), method = "spearman")
famtab_cor <- as.data.frame(famtab_cor)
View(famtab_cor)
corr_cross(famtab_cor, # name of dataset
           plot= TRUE,
           ignore = NA,
           max=0.99, # max correlation value
           top = 50 # display top 20 couples of variables (by correlation coefficient)
)

cross_corr_fam.99<- corr_cross(famtab_cor, # name of dataset
                               plot= TRUE,
                               ignore = NA,
                               max=0.99, # max correlation value
                               top = 6000 # display top 20 couples of variables (by correlation coefficient)
)
View(cross_corr_fam.99$data)

#correlations
genustab_cor = cor(as.matrix(genustab), method = "spearman")
View(genustab_cor)
genustab_cor[genustab_cor < 0.6 | genustab_cor ==1] <- NA
View(genustab_cor)
genustab_cor_filter <-genustab_cor[rowSums(is.na(genustab_cor)) != ncol(genustab_cor), ]
genustab_cor %>% filter_all(all_vars(.> 0.6))
View(genustab_cor_filter)
#filter correlations
genustab_cor_filter <- filter(genustab_cor <= -0.5 | genustab_cor >= 0.8)
corrplot(genustab_cor)
jpeg("~/2017_chem_kendall_corrplot.jpg",width=2400, height=2200, res=400, bg="white")
plot.new()
corrplot(geo17_cor, method="circle",#addCoef.col = 'gray40',
         tl.col = "black", col = COL2('PuOr'),type="lower",
         order = "hclust", hclust.method = "average", addrect = 5, tl.cex = 0.7,diag=TRUE)
dev.off()


#spearman correlation analysis
cortest <- WGCNA::corAndPvalue(genustab, method="spearman", alternative="two.sided")
cortest$cor[upper.tri(cortest$cor, diag = TRUE)] <- NA
cortest$p[upper.tri(cortest$p, diag = TRUE)] <- NA
cortab1 <- na.omit(melt(t(cortest$cor))) %>% rename(from=Var1,to=Var2,cor=value)
corptab1 <- na.omit(melt(t(cortest$p))) %>% rename(pvalue=value)
cortab1$fdr <- p.adjust(corptab1$pvalue, method="fdr")

cortab1 <- cortab1 %>% mutate(correlation=case_when(cor>0 ~ "positive",cor < 0 ~ "negative",TRUE ~ "No"))
cortab2 <- cortab1 %>% filter(fdr <= 0.05) %>% filter(cor <= -0.5 | cor >= 0.8)
#end of taxonomic correlations
#####

# phylum-specific heatmap
#####
phylo2017_rarefied
prune_phylo2017 <-phyloseq::prune_taxa(taxa_sums(p2017_1_filtered) > 500, p2017_1_filtered)
table(phyloseq::tax_table(prune_phylo2017)[, "Phylum"])
ps_phylum <- phyloseq::tax_glom(prune_phylo2017, "Phylum")
phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
melt_pr_phylo2017_rarefied <- phyloseq::psmelt(prune_phylo2017)
ps_phylum
ps_phylum.df <- as.data.frame(t(otu_table(ps_phylum)))
ps_phylum.df

colnames(ps_phylum.df)
ps_samdat.df <- as.data.frame(sample_data(phylo2017_rarefied))
ps_phylum_merge <- merge(x=ps_phylum.df, y=ps_samdat.df,
                         by = 'row.names', all = TRUE)
#design elements for heatmap
col.pal <- brewer.pal(9,"YlOrBr")
category_df = data.frame("Well" = ps_phylum_merge$well_type)
category_df$Days <- ps_phylum_merge$Days
category_df$Filter <- ps_phylum_merge$FilterSize
category_df$Phase <- ps_phylum_merge$phase
category_df
unique(category_df$Filter)

rownames(category_df) <- rownames(ps_samdat.df) # name matching
head(category_df)
category_df$Days <- factor(category_df$Days, levels = c("-6","1" ,  "8" ,  "15" , "22" , "50", "78" , "106" ,"134"))
category_df$Phase <- factor(category_df$Phase, levels = c("Control", "Phase0",  "Phase1", "Phase2", "Phase3" , "Phase4",  "Phase5",  "Phase6"   ))

anno_colors = list(
  Days = c("-6" = "#ffffd9", "1" = "#edf8b1",  "8"  = "#c7e9b4",  "15" = "#7fcdbb", "22" = "#41b6c4" , "50"  = "#1d91c0", "78" = "#225ea8", "106" = "#253494","134" = "#081d58"),
  Well = c(Control = "tomato3", Monitoring = "pink2"),
  #Well = c(Control = "darkred", Monitoring = "darksalmon"),
  #Well = c(Control = "black", Monitoring = "blueviolet"),
  Filter = c("0.1um" = "black", "0.2um" = "blueviolet"),
  Phase = c("Control" = "white","Phase0"=  "#ffffe5","Phase1"= "#f7fcb9" ,"Phase2"= "#d9f0a3","Phase3"= "#41ab5d" ,"Phase4"= "#238b45" ,"Phase5"= "#006d2c" ,"Phase6"= "#00441b" )
  #Filter = c("0.1um" = "palegreen", "0.2um" = "darkgreen")
)

library(microbiomeutilities)
heat.phylo2017 <- plot_taxa_heatmap(phylo2017_rarefied,
                                    subset.top = 20,
                                    VariableA = c("Days","well_type", "FilterSize"),
                                    heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                                    transformation = "log10"
)

ps_phylum.df[ps_phylum.df == 0] <- NA
ps_phylum.df_log <- as.data.frame(log(ps_phylum.df) ) 
ps_phylum.df_log[is.na(ps_phylum.df_log)] = 0
rownames(category_df) <- rownames(ps_phylum.df_log) # name matching
head(ps_phylum.df_log)
head(ps_phylum.df)
geo2017_log_numeric
#make pheatmap with day annotation
jpeg("~/2017_heatmap_phylum_top500.jpeg",width=2200, height=2800, res=400, bg="white")
plot.new()
ps_phylum.df_log
p2017_tax_heatmap <-pheatmap(
  mat               = ps_phylum.df_log,
  color             = col.pal,
  border_color      = 'gray80',
  show_colnames     = TRUE,
  cutree_rows = 7,
  cutree_cols = 4,
  treeheight_row = 20,
  treeheight_col = 20,
  annotation_row = category_df,
  annotation_colors = anno_colors,
  angle_col=315,
  show_rownames     = FALSE,
  drop_levels       = TRUE,
  fontsize          = 10
)
getOption("device")
dev.set(which = dev.next())
dev.off()

#merged heatmap
phylo2017_rarefied
prune_merged_phylo2017 <-phyloseq::prune_taxa(taxa_sums(merged_phylo2017_rarefied) > 500, merged_phylo2017_rarefied)
table(phyloseq::tax_table(prune_merged_phylo2017)[, "Phylum"])
ps_merged_phylum <- phyloseq::tax_glom(prune_merged_phylo2017, "Phylum")
phyloseq::taxa_names(ps_merged_phylum) <- phyloseq::tax_table(ps_merged_phylum)[, "Phylum"]
melt_merged_phylo2017_rarefied <- phyloseq::psmelt(prune_merged_phylo2017)
#ps_phylum
ps_merged_phylum.df <- as.data.frame(otu_table(ps_merged_phylum))
ps_merged_phylum.df

colnames(ps_merged_phylum.df)
head(sample_data(merged_phylo2017_rarefied))
sample_data()
ps_merged_samdat.df <- as.data.frame(sample_data(merged_phylo2017_rarefied))
ps_merged_samdat.df
merged_ps_phylum_merge <- merge(x=ps_merged_phylum.df, y=ps_merged_samdat.df,
                                by = 'row.names', all = TRUE)
merged_ps_phylum_merge
ps_merged_phylum.df[ps_merged_phylum.df== 0] <- NA
ps_merged_phylum.df_log <- as.data.frame(log(ps_merged_phylum.df) ) 
ps_merged_phylum.df_log[is.na(ps_merged_phylum.df_log)] = 0

category_merge_df = data.frame("Days" = ps_merged_samdat.df$DayCategory)
category_merge_df$Well <- as.character(ps_merged_samdat.df$well_type)
#category_merge_df$Phase <- ps_phylum_merge$phase
category_merge_df
unique(category_merge_df$Days)

rownames(category_merge_df) <- rownames(ps_merged_phylum.df_log) # name matching
head(category_df)
category_merge_df$Days <- factor(category_merge_df$Days, levels = c("-6","1" ,  "8" ,  "15" , "22" , "50", "78" , "106" ,"134", "NA"))

anno_merge_colors = list(
  Days = c("-6" = "#ffffd9", "1" = "#edf8b1",  "8"  = "#c7e9b4",  "15" = "#7fcdbb", "22" = "#41b6c4" , "50"  = "#1d91c0", "78" = "#225ea8", "106" = "#253494","134" = "#081d58", "NA" = "White"),
  #Well = c(Control = "tomato3", Monitoring = "pink2"),
  #Well = c(Control = "darkred", Monitoring = "darksalmon"),
  Well = c("Control" = "black", "Monitoring" = "blueviolet")
  #Filter = c("0.1um" = "black", "0.2um" = "blueviolet"),
  #Phase = c("Control" = "white","Phase0"=  "#ffffe5","Phase1"= "#f7fcb9" ,"Phase2"= "#d9f0a3","Phase3"= "#41ab5d" ,"Phase4"= "#238b45" ,"Phase5"= "#006d2c" ,"Phase6"= "#00441b" )
  #Filter = c("0.1um" = "palegreen", "0.2um" = "darkgreen")
)
category_merge_df

rownames(category_df) <- rownames(ps_merged_phylum.df_log) # name matching
head(ps_merged_phylum.df_log)
head(ps_merged_phylum.df)
geo2017_log_numeric
#make pheatmap with day annotation
jpeg("~/2017_heatmap_phylum_combined_top500.jpeg",width=2200, height=2800, res=400, bg="white")
plot.new()
pheatmap(
  mat               = ps_merged_phylum.df_log,
  color             = col.pal,
  border_color      = 'gray80',
  show_colnames     = TRUE,
  cutree_rows = 7,
  cutree_cols = 4,
  treeheight_row = 20,
  treeheight_col = 20,
  annotation_row = category_merge_df,
  annotation_colors = anno_merge_colors,
  angle_col=315,
  show_rownames     = TRUE,
  drop_levels       = TRUE,
  fontsize          = 10
)
getOption("device")
dev.set(which = dev.next())
dev.off()
#if dev.off() is not working run this to reset the device.next()
#while (!is.null(dev.list()))  dev.off()

#combined taxa heatmap
p2017_combined_filtered
prune_merged_p2017_combined <-phyloseq::prune_taxa(taxa_sums(p2017_combined_filtered) > 200, p2017_combined_filtered)
table(phyloseq::tax_table(prune_merged_p2017_combined)[, "Phylum"])
ps_merged_phylum <- phyloseq::tax_glom(prune_merged_p2017_combined, "Phylum")
phyloseq::taxa_names(ps_merged_phylum) <- phyloseq::tax_table(ps_merged_phylum)[, "Phylum"]
melt_merged_p2017_combined <- phyloseq::psmelt(prune_merged_p2017_combined)
#ps_phylum
ps_merged_phylum.df <- as.data.frame(otu_table(ps_merged_phylum))
ps_merged_phylum.df <- t(ps_merged_phylum.df)
head(ps_merged_phylum.df)

colnames(ps_merged_phylum.df)
head(sample_data(p2017_combined_filtered))
sample_data()
ps_combined_samdat.df <- as.data.frame(sample_data(p2017_combined_filtered))
head(ps_combined_samdat.df)
combined_ps_phylum_merge <- merge(x=ps_merged_phylum.df, y=ps_combined_samdat.df,
                                  by = 'row.names', all = TRUE)
combined_ps_phylum_merge
ps_merged_phylum.df[ps_merged_phylum.df== 0] <- NA
ps_merged_phylum.df_log <- as.data.frame(log(ps_merged_phylum.df) ) 
ps_merged_phylum.df_log[is.na(ps_merged_phylum.df_log)] = 0
head(ps_merged_phylum.df_log)
category_merge_df = data.frame("Days" = ps_combined_samdat.df$DayCategory)
category_merge_df$Well <- as.character(ps_combined_samdat.df$well_type)
category_merge_df$Filter <- as.character(ps_combined_samdat.df$FilterSize)
#category_merge_df$Phase <- ps_phylum_merge$phase
category_merge_df
unique(category_merge_df$Days)

rownames(category_merge_df) <- rownames(ps_merged_phylum.df_log) # name matching
head(category_merge_df)
category_merge_df$Days <- factor(category_merge_df$Days, levels = c("-6","1" ,  "8" ,  "15" , "22" , "50", "78" , "106" ,"134", "NA"))

anno_merge_colors = list(
  Days = c("-6" = "#ffffd9", "1" = "#edf8b1",  "8"  = "#c7e9b4",  "15" = "#7fcdbb", "22" = "#41b6c4" , "50"  = "#1d91c0", "78" = "#225ea8", "106" = "#253494","134" = "#081d58", "NA" = "White"),
  #Well = c(Control = "tomato3", Monitoring = "pink2"),
  Well = c(Control = "darkred", Monitoring = "darksalmon"),
  #Well = c("Control" = "black", "Monitoring" = "blueviolet")
  Filter = c("0.1um" = "black", "0.2um" = "blueviolet")
  #Phase = c("Control" = "white","Phase0"=  "#ffffe5","Phase1"= "#f7fcb9" ,"Phase2"= "#d9f0a3","Phase3"= "#41ab5d" ,"Phase4"= "#238b45" ,"Phase5"= "#006d2c" ,"Phase6"= "#00441b" )
  #Filter = c("0.1um" = "palegreen", "0.2um" = "darkgreen")
)
category_merge_df

rownames(category_df) <- rownames(ps_merged_phylum.df_log) # name matching
head(ps_merged_phylum.df_log)
head(ps_merged_phylum.df)
geo2017_log_numeric
#make pheatmap with day annotation
jpeg("~/2017_heatmap_phylum_0.1+0.2_combined_top200.jpeg",width=2200, height=2800, res=400, bg="white")
plot.new()
pheatmap(
  mat               = ps_merged_phylum.df_log,
  color             = col.pal,
  border_color      = 'gray80',
  show_colnames     = TRUE,
  cutree_rows = 5,
  cutree_cols = 5,
  treeheight_row = 20,
  treeheight_col = 20,
  annotation_row = category_merge_df,
  annotation_colors = anno_merge_colors,
  angle_col=315,
  show_rownames     = FALSE,
  drop_levels       = TRUE,
  fontsize          = 10
)
getOption("device")
dev.set(which = dev.next())
dev.off()


#end of heatmaps
#####

#compare differences in alpha statistics
#####
install.packages('sjPlot')
library(sjPlot)
#shannon changes by day
#calculate sample alpha stats
#0.1 only
p2017_1_sample_data <- as.data.frame(sample_data(p2017_1_filtered))
head(p2017_1_sample_data)
p2017_1_sample_data$N_cat <- as.factor(p2017_1_sample_data$nitrate_level);p2017_1_sample_data$S_cat <- as.factor(p2017_1_sample_data$sulfate_level)
colnames(p2017_1_sample_data)
row.names(p2017_1_sample_data)
rich_1 = estimate_richness(p2017_1_filtered, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))
row.names(rich_1) <-row.names(p2017_1_sample_data)
merged_stat_1 <-as.data.frame(merge(rich_1, p2017_1_sample_data, by=0, all=TRUE))
head(merged_stat_1)
merged_stat_1$Acetate_range <-case_when((merged_stat_1$Acetate <= 1)              ~ 1,
                                        (merged_stat_1$Acetate  > 1) & (merged_stat_1$Acetate <= 10) ~ 10,
                                        (merged_stat_1$Acetate  > 10) & (merged_stat_1$Acetate <= 20) ~ 20,
                                        (merged_stat_1$Acetate  > 20) & (merged_stat_1$Acetate <= 30) ~ 30,
                                        (merged_stat_1$Acetate  > 30) & (merged_stat_1$Acetate<= 50) ~ 40)
#Simpson and Shannon stats for well types and days 
summary(aov(data=merged_stat_1, Shannon~well_type*DayCategory+nitrate_level+sulfate_level))
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
write.csv(merged_stat_2, "~/merged_2_alpha_stat.csv")

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
write.csv(merged_C_stat, "~/small_bacteria_alpha_stat.csv")

#merged_combined
p2017_Cm_sample_data <- as.data.frame(sample_data(p2017_combined_merged))
colnames(p2017_Cm_sample_data)
row.names(p2017_Cm_sample_data)
rich_cm = estimate_richness(p2017_combined_merged, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))
row.names(rich_cm) <-row.names(p2017_Cm_sample_data)
merged_Cm_stat <-as.data.frame(merge(rich_cm, p2017_Cm_sample_data, by=0, all=TRUE))
head(merged_Cm_stat)
write.csv(merged_Cm_stat, "~/small_bacteria_mergedsamples_alpha_stat.csv")

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
write.csv(p2017_merged_stat, "~/p2017_merged_alpha_stat.csv")

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
ggsave("~/p2017_lm_day_well_shannon_plot.pdf", height = 6, width = 6)

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
ggsave("~/p2017_lm_day_well_shannon.pdf", height = 6, width = 6)
#chao1 changes by day
p2017_merged_chao1 <- plot_richness(merged_phylo2017_rarefied, x="DayCategory", color="well_type_factor", measures = "Chao1")
p2017_merged_chao1
#test differences in simpson means with linear model
tab_model(lm(value ~DayCategory*well_type_factor, p2017_merged_chao1$data))
ggsave("~/p2017_lm_day_well_chao1.pdf", height = 6, width = 6)
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
ggsave("~/p2017_lm_day_well_invsimpson_plot.pdf", height = 6, width = 6)

#test differences in simpson means with linear model
tab_model(lm(value ~DayCategory*well_type_factor, p2017_merged_InvSimpson$data))
ggsave("~/p2017_lm_day_well_invsimpson.pdf", height = 6, width = 6)

#InvSimp. by nitrate
p2017_merged_NO3_shannon <- plot_richness(merged_phylo2017_rarefied, x="nitrate_range", color="well_type_factor", measures = "Shannon")
aov(value ~nitrate_range*well_type_factor, data=p2017_merged_NO3_shannon$data)
tab_model(lm(value ~nitrate_range*well_type_factor, data=p2017_merged_NO3_shannon$data))
ggsave("~/p2017_lm_nitrate_well_shannon.pdf", height = 6, width = 6)

p2017_merged_SO4_shannon <- plot_richness(merged_phylo2017_rarefied, x="sulfate_range", color="well_type_factor", measures = "Shannon")
aov(value ~sulfate_range*well_type_factor, data=p2017_merged_SO4_shannon$data)
tab_model(lm(value ~sulfate_range*well_type_factor, p2017_merged_SO4_shannon$data))
ggsave("~/p2017_lm_nitrate_well_shannon.pdf", height = 6, width = 6)

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
ggsave("~/p2017_lm_day_well_invsimpson_plot.pdf", height = 6, width = 6)

#test differences in simpson means with linear model
tab_model(lm(value ~nitrate_range*well_type_factor, p2017_merged_NO3_shannon$data))
ggsave("~/p2017_lm_day_well_invsimpson.pdf", height = 6, width = 6)
#end of alpha stats
#####

#differential abundance of the full community using data transformed in a distance matrix
#####
library("nlme")
library("reshape2")

#CCA plot and analysis
require(vegan) 
data(varespec) 
data(varechem)

#merged_p2017 meta data
#####
#extract merged sample data
merged_p2017_metadata <- as(sample_data(merged_p2017), "data.frame") ## convert sample_data to data.frame
merged_p2017_sample_data <-sample_data(merged_p2017)
merged_p2017_sample_df <- as.data.frame(merged_p2017_sample_data)
#numeric data only
merged_p2017_numeric <- as.matrix(merged_p2017_sample_data);merged_p2017_numeric  <- merged_p2017_numeric[,5:19]
merged_p2017_numeric
#end of making metadata frame
#####
head(merged_p2017_sample_data)

#ordination with merged_p2017
#####
#CCA
merged_p2017_CCA  <- ordinate(merged_p2017~Days*well_type_factor, "CCA") #ordinate by CCA
merged_p2017_CCA
plot_scree(merged_p2017_CCA , "Scree Plot for Correspondence Analysis") #scree plot

# unifrac distance
set.seed(0451)
merged_p2017_genus_unifrac = UniFrac(tax_glom(merged_p2017, taxrank = "Genus"))
merged_p2017_genus <- as.matrix(otu_table(tax_glom(merged_p2017, taxrank = "Genus")))
View(merged_p2017_genus)
merged_p2017_genus_unifrac 

#bray curtis
merged_p2017_genus_bray <- vegdist(wisconsin(sqrt(merged_p2017_genus)), method = "bray")
merged_p2017_genus_bray

#NMDS of UNIFRAC
set.seed(0451)
p2017_merge_genus.mds <- metaMDS(merged_p2017_genus_unifrac, zerodist=ignore, try=999 ,trymax = 1000)
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
ggsave('~/2017_unifrac_nmds_ordination.png', bg = "transparent", width = 40, height = 30, units = "cm")
#end NMDS
#####

#merged_p2017 Canonicanl Analysis (CAP)
#####
merged_p2017_metadata
#multivariate anova to determine what percent of total variation is explained  by covariates
adonis2(merged_p2017_genus_unifrac ~ Days + well_type_factor + nitrate_level + sulfate_level + acetate_presence+ammonium_presence, data = merged_p2017_metadata, perm = 999)
p2017_merged_CCA  <- ordinate(merged_p2017, "CCA") #ordinate by CCA
p2017_merged_CAP <-p2017_merged_CCA$CA$u
p2017_merged_CCA_metadata; plot_ordination(merged_p2017, ordinate(merged_p2017, "CCA"), "samples", color = "Days", shape = "well_type")
#Legendre, P. & Anderson, M.J. (1999). Distance-based redundancy analysis: testing multispecies responses in multifactorial ecological experiments. Ecological Monographs 69: 1-24.
#Anderson, M.J. & Willis, T.J. (2003). Canonical analysis of principal coordinates: a useful method of constrained ordination for ecology. Ecology 84: 511-525.
cap_genus <- capscale(merged_p2017_genus_unifrac ~  +Nitrate + Sulfate + Uranium + Acetate + Ammonium, data = merged_p2017_metadata)
cap_genus
plot(cap_genus)
anova(cap_genus)
cap_genus_cols <-cap_genus$CCA$wa
cap_genus_cols[,1:2]
#export vector as CCA
cap_genus_vector <-as.data.frame(cap_genus$CCA$biplot)
cap_genus_vector <- (cap_genus_vector/2)
cap_genus_vector$sites <-row.names(cap_genus_vector);cap_genus_vector
#merge cap_genus data with metadata
cap_genus_metadata <- cbind(merged_p2017_metadata, data.frame(CAP1 = cap_genus_cols[,1], CAP2=cap_genus_cols[,2])) 
cap_genus_metadata
head(cap_genus_metadata)
head(p2017_merged_CCA_metadata)
plot(cap_genus)

#CAP PLOT
#set color levels in data
cap_genus_metadata$Days <- factor(cap_genus_metadata$Days, levels = c("-6" , "1" ,  "8" ,  "15" , "22" , "50" , "78" , "106" ,"134"))
cap_genus_metadata$Phase <- factor(cap_genus_metadata$phase_C12, levels = c("Control", "Phase0",  "Phase1", "Phase2", "Phase3" ))
cap_genus_metadata$type_day <- factor(cap_genus_metadata$type_day, levels = c("Control", "Monitoring Day-6",  "Monitoring Day1", "Monitoring Day8", "Monitoring Day15", "Monitoring Day22", "Monitoring Day50", "Monitoring Day78", "Monitoring Day106", "Monitoring Day134"))
cap_genus_metadata$phase_cat <- factor(cap_genus_metadata$phase_cat, levels = c("0", "1",  "2", "3"))

#plot CAP analysis
ggplot() +
  #facet_wrap(~diet) + 
  scale_colour_brewer(palette = "YlGnBu") + 
  scale_fill_brewer(palette = "YlGnBu") + 
  geom_hline(yintercept = 0, color = "gray50", size = 0.5, alpha=0.6, linetype="dashed")+
  geom_vline(xintercept = 0, color ="gray50", size=0.5, alpha=0.6, linetype="dashed")+
  stat_ellipse(data=cap_genus_metadata, aes(x = CAP1, y = CAP2, fill = phase_cat),  linetype = 2, geom = "polygon", alpha = 0.6) +# add ellipse
  geom_segment(data = cap_genus_vector, # vector arrows
               aes(x = 0, xend = CAP1, y = 0, yend = CAP2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "gray35") + #scores
  #geom_point(data = p2017_merged_cent, size = 5) +            # centroids
  geom_point(data = cap_genus_metadata, aes(x = CAP1, y = CAP2,   shape=well_type),size=5.5, color = "black")+
  geom_point(data = cap_genus_metadata, aes(x = CAP1, y = CAP2,  color = as.factor(type_day), shape=well_type),size=4) +   
  geom_label(data = cap_genus_vector, aes(x = CAP1, y = CAP2, label = sites), # vector labels
             size = 3.5, check_overlap = TRUE) +# sample scores
  geom_text(aes(x=0.25, label="Significance of constraints:\n Sum Sq=2.18\n F=3.45\n p=0.001", y=-0.6), colour="gray30", angle=0, text=element_text(size=7)) +
  labs(shape = "Well Type", color = "Well and Day", fill="Injection Phase", y="CAP2 [6.21%]", x="CAP1 [19.8%]") + 
  coord_fixed() + #must have a fixed centroid
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="vertical",
        legend.position = "right",
        axis.ticks = element_line(color = "gray50"),
        #legend.position = c(0.46, 0.92),
        axis.title = element_text(size=25, color = "black"),
        plot.title = element_text(size=25, color = "black"),
        axis.text = element_text(size = 18, color = "black"),
        panel.border = element_rect(colour = "gray50", size=0.8, fill = "transparent"))                                        # same axis scaling
ggsave('~/2017_unifrac_CAP_ordination.jpeg', bg = "white", width = 30, height = 25, units = "cm")
#end CAP
#####

#merged_p2017 anosim and adonis
#####
#anosim testing p.133 ADP-003 to test if composition differs between groups
#groups DO NOT differ in composition with well_type
anosim(merged_p2017_genus_unifrac, merged_p2017_sample_data$well_type)
#groups DO!! differ in composition with NITRATE. R: 0.2202; p: 0.001
anosim(merged_p2017_genus_unifrac, merged_p2017_sample_data$nitrate_level)
#groups DO!! differ in composition with SULFATE. R: 0.3341; p: 0.001
anosim(merged_p2017_genus_unifrac, merged_p2017_sample_data$sulfate_level)
#groups DO!! differ in composition with DAYS. R: 0.5724; p: 0.001
anosim(merged_p2017_genus_unifrac, merged_p2017_sample_data$Days)
#groups DO!! differ in composition with URANIUM. R: 0.2528; p: 0.003
anosim(merged_p2017_genus_unifrac, merged_p2017_sample_data$uranium_level)
#groups DO!! differ in composition by phase. R2: 0.6621; p:0.001
anosim(merged_p2017_genus_unifrac, merged_p2017_sample_data$phase_cat)
#ammonium 
anosim(merged_p2017_genus_unifrac, merged_p2017_sample_data$ammonium_presence)
(sample_data(merged_p2017)$Ammonium)
#acetate
anosim(merged_p2017_genus_unifrac, merged_p2017_sample_data$acetate_presence)
#adonis test (p.133-134 ADP-003)
head(sample_data(merged_p2017))
(microbes.adonis <- adonis2(merged_p2017_genus_unifrac ~ well_type_factor *Days + nitrate_level +sulfate_level+ acetate_presence+ammonium_presence, as(sample_data(merged_p2017), "data.frame"), 
                            permutations = 999, by=NULL))
(sulfate.adonis <- adonis2(merged_p2017_genus_unifrac ~ well_type_factor *Days, as(sample_data(merged_p2017), "data.frame"), 
                           permutations = 999, by=NULL))
(phase.adonis <- adonis2(merged_p2017_genus_unifrac ~ phase_cat, as(sample_data(merged_p2017), "data.frame"), 
                         permutations = 999, by=NULL))
(microbes.adonis <- adonis2(p2017_rarefied_genus_unifrac ~  nitrate_range * sulfate_range  *type_day , as(sample_data(merged_phylo2017_rarefied), "data.frame"), 
                            permutations = 10000, by="terms"))

plot(microbes.adonis)
#end adonis and anosim
#####

#betadispersion test (p.135 ADP-003)
sample_data(merged_p2017)
#####
anova(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$phase))
anova(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$well_type))
plot(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$well_type))
anova(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$nitrate_level))
plot(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$nitrate_level))
anova(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$sulfate_level))
plot(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$sulfate_level))
anova(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$uranium_level))
plot(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$uranium_level))
anova(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$sample_label))
plot(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$type_day))
betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$type_day)
anova(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$type_day))
anova(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$phase_cat))
betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$phase_cat)
anova(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$phase_C12))
betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$phase_C12)
anova(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$phase_C134))
betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$phase_C134)
plot(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$ammonium_presence))
betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$ammonium_presence)
anova(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$ammonium_presence))
plot(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$acetate_presence))
betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$acetate_presence)
anova(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$acetate_presence))
#extract centroids and sites from phase, type_day and well_type beta dispersion tests
#type_day
merged_betadisp_type_day <-plot(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$type_day))
anova(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$type_day))
merged_betadisp_type_day_sites <- as.data.frame(merged_betadisp_type_day$sites)
merged_betadisp_type_day_sites$label <- row.names(merged_betadisp_type_day_sites)
merged_betadisp_type_day_centroid <- as.data.frame(merged_betadisp_type_day$centroids)
merged_betadisp_type_day_centroid$label <- row.names(merged_betadisp_type_day_centroid)
merged_betadisp_type_day_centroid
merged_betadisp_type_day_sites <- cbind(merged_p2017_sample_data, data.frame(PCoA1 = merged_betadisp_type_day$sites[,1], PCoA2 = merged_betadisp_type_day$sites[,2]))
merged_betadisp_type_day_sites $grp <-grp
merged_betadisp_type_day_sites$phase_cat <- factor(merged_betadisp_type_day_sites$phase_cat, levels = c("0", "1",  "2", "3"))

#phase
merged_betadisp_phase_cat <-plot(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$phase_cat))
merged_betadisp_phase_cat_sites <- as.data.frame(merged_betadisp_phase_cat$sites)
merged_betadisp_phase_cat_sites$label <- row.names(merged_betadisp_phase_cat_sites)
merged_betadisp_phase_cat_centroid <- as.data.frame(merged_betadisp_phase_cat$centroids)
merged_betadisp_phase_cat_centroid$label <- row.names(merged_betadisp_phase_cat_centroid)
head(merged_betadisp_phase_cat_centroid)
#make dataframe for plot
merged_betadisp_phase_cat_sites <- cbind(merged_p2017_sample_data, data.frame(PCoA1 = merged_betadisp_phase_cat$sites[,1], PCoA2 = merged_betadisp_phase_cat$sites[,2]))
merged_betadisp_phase_cat_sites $grp <-grp
merged_betadisp_phase_cat_sites$phase_cat <- factor(merged_betadisp_phase_cat_sites$phase_cat, levels = c("0", "1",  "2", "3"))
head(merged_betadisp_phase_cat_sites)
#well_type
anova(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$well_type))
betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$well_type)
merged_betadisp_well_type <-plot(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$well_type))
merged_betadisp_well_type_sites <- as.data.frame(merged_betadisp_well_type$sites)
merged_betadisp_well_type_sites$label <- row.names(merged_betadisp_well_type_sites)
merged_betadisp_well_type_centroid <- as.data.frame(merged_betadisp_well_type$centroids)
merged_betadisp_well_type_centroid$label <- row.names(merged_betadisp_well_type_centroid)

merged_betadisp_type_day_sites <- cbind(merged_p2017_sample_data, data.frame(PCoA1 = merged_betadisp_type_day$sites[,1], PCoA2 = merged_betadisp_type_day$sites[,2]))
merged_betadisp_type_day_sites $grp <-grp
merged_betadisp_type_day_sites$phase_cat <- factor(merged_betadisp_type_day_sites$phase_cat, levels = c("0", "1",  "2", "3"))

#permutation tests of dispersion within the monitoring wells # indicates a significant community shift in monitoring wells related to days
merged_p2017_genus.monitoring <- tax_glom(subset_samples(merged_p2017, well_type == "Monitoring"), taxrank = "Genus")  
merged_p2017_genus_unifrac.monitoring = UniFrac(merged_p2017_genus.monitoring)
(merged_p2017_genus.monitoring.adonis <- adonis2(merged_p2017_genus_unifrac.monitoring ~DayCategory *nitrate_level * phase_cat*sulfate_level * uranium_level , as(sample_data(merged_p2017_genus.monitoring), "data.frame"), permutations = 999, by="term"))
#permutation tests of dispersion within the CONTROL well # indicates a significant community shift in control wells related to days
merged_p2017_genus.control <- tax_glom(subset_samples(merged_p2017, well_type == "Control"), taxrank = "Genus")  
sample_data(merged_p2017_genus.control)
merged_p2017_genus_unifrac.control = UniFrac(merged_p2017_genus.control)
(merged_p2017_genus.control.adonis <- adonis2(merged_p2017_genus_unifrac.control ~Days *nitrate_level + phase_cat+sulfate_level + uranium_level , as(sample_data(merged_p2017_genus.control), "data.frame"), permutations = 999, by="term"))
merged_p2017_genus.control.adonis
head(merged_betadisp_type_day_sites)

#plot dispersion eigen values
ggplot(merged_betadisp_phase_cat_sites, aes(x = PCoA1, y = PCoA2,  color = phase_cat, shape = well_type)) +
  #facet_wrap(~diet) + 
  scale_colour_brewer(palette = "YlGnBu") + 
  scale_fill_brewer(palette = "YlGnBu") + 
  geom_hline(yintercept = 0, color = "gray50", size = 0.5, alpha=0.6, linetype="dashed")+
  geom_vline(xintercept = 0, color ="gray50", size=0.5, alpha=0.6, linetype="dashed")+
  stat_ellipse(data=merged_betadisp_phase_cat_sites, aes(x=PCoA1,y=PCoA2, fill = phase_cat, color = phase_cat),  linetype = 2, geom = "polygon", alpha = 0.6) +# add ellipse
  geom_point(size=5.5, color = "black")+
  geom_text(aes(x=0.3, label="Variance in group dispersion:\n F=14.78\n p≤0.0001", y=-0.32), colour="gray30", angle=0, text=element_text(size=10)) +
  geom_point(size=4) +   
  #geom_label(data = merged_betadisp_phase_cat_centroid, aes(x=PCoA1,y=PCoA2, label = label),size = 3, check_overlap = TRUE) +# sample scores
  #geom_point(data = merged_betadisp_type_day_centroid, aes(x = PCoA1, y = PCoA2),color = "blue",size = 6) +
  labs(shape = "", color = "", fill="Phase", title = "Unweighted UniFrac of the merged microbial community ") + 
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = "bottom",
        axis.ticks = element_line(color = "gray50"),
        #legend.position = c(0.46, 0.92),
        axis.title = element_text(size=25, color = "black"),
        plot.title = element_text(size=25, color = "black"),
        axis.text = element_text(size = 18, color = "black"),
        panel.border = element_rect(colour = "gray50", size=0.8, fill = "transparent"))                                       # same axis scaling
ggsave('~/2017_betadispersion_phase.jpeg', bg = "white", width = 20, height = 17, units = "cm")

#testing unweighted distance and betadispersion in the overlapping community
ggplot(merged_betadisp_sites, aes(x = PCoA1, y = PCoA2,  color = as.factor(grp), shape = well_type)) +
  #facet_wrap(~diet) + 
  #scale_colour_brewer(palette = "YlGnBu") + 
  geom_hline(yintercept = 0, color = "gray50", size = 0.5, alpha=0.6, linetype="dashed")+
  geom_vline(xintercept = 0, color ="gray50", size=0.5, alpha=0.6, linetype="dashed")+
  stat_ellipse(aes(x=PCoA1,y=PCoA2, fill = as.factor(grp), color = as.factor(grp)),  linetype = 2, geom = "polygon", alpha = 0.6) +# add ellipse
  geom_point(size=5.5, color = "black")+
  geom_text(aes(x=0.35, label="Variance across groups:\n F=2.47\n p≤0.05", y=-0.32), colour="gray30", angle=0, text=element_text(size=8)) +
  geom_point(size=4) +                                              # sample scores
  #geom_point(data = merged_betadisp_type_day_centroid, aes(x = PCoA1, y = PCoA2),color = "blue",size = 6) +
  labs(shape = "Well Type", color = "", fill="Phase") + 
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = "bottom",
        axis.ticks = element_line(color = "gray50"),
        #legend.position = c(0.46, 0.92),
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        panel.border = element_rect(colour = "gray50", size=0.8, fill = "transparent"))                                        # same axis scaling
ggsave("~/2017_betadispersion_type_day.png", bg = "transparent", width = 60, height = 35, units = "cm")

plot(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$phase_cat))
disp_p12 <- plot(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$phase_C12))
disp_p134 <- plot(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$phase_C134))
pdf(disp_p12,"~/2017_p2017_disp_p12.pdf")
pdf(disp_p134,"~/2017_p2017_disp_p134.pdf")


#end beta diversity testing - merged
#####

#p2017_1meta data
#####
#extract merged sample data
p2017_1 = subset_samples(p2017_1, Date != "8/16/18")#remove post-injection samples
p2017_1 <- prune_taxa(taxa_sums(p2017_1) > 0, p2017_1)
p2017_1_metadata <- as(sample_data(p2017_1), "data.frame") ## convert sample_data to data.frame
p2017_1_sample_data <-sample_data(p2017_1)
p2017_1_sample_df <- as.data.frame(p2017_1_sample_data)
#numeric data only
p2017_1_numeric <- as.matrix(p2017_1_sample_data)
p2017_1_numeric
p2017_1_numeric  <- p2017_1_numeric[,8:36]
p2017_1_numeric
#end of making metadata frame
#####

#ordination with p2017_1
#####
#CCA
p2017_1_CCA  <- ordinate(p2017_1~Days*well_type, "CCA") #ordinate by CCA
p2017_1_CCA
plot_scree(p2017_1_CCA , "Scree Plot for Correspondence Analysis") #scree plot

# unifrac distance
set.seed(0451)
p2017_1_genus_unifrac = UniFrac(tax_glom(p2017_1, taxrank = "Genus"))
p2017_1_genus <- as.matrix(otu_table(tax_glom(p2017_1, taxrank = "Genus")))
p2017_1_genus
p2017_1_genus_unifrac 

#bray curtis
p2017_1_genus_bray <- vegdist(wisconsin(sqrt(p2017_1_genus)), method = "bray")
p2017_1_genus_bray

#NMDS of UNIFRAC
set.seed(0451)
p2017_1_genus.mds <- metaMDS(p2017_1_genus_unifrac, zerodist=ignore, try=999 ,trymax = 1000)
p2017_1_genus.mds; plot(p2017_1_genus.mds)
ASV_tree <- hclust(p2017_1_genus_unifrac, method = "average");plot(ASV_tree)
grp <- cutree(ASV_tree, k = 3)
#nmds.scores
nmds.scores <- as.data.frame(scores(p2017_1_genus.mds, display = "sites"))
nmds.scores$grp <- grp
nrow(p2017_1_numeric); nrow(p2017_1_genus_unifrac)
p2017_1_numeric
#end of ordination
#####

#p2017_1 NMDS
#####
#ordinate with metadaata
env_ord <- envfit(p2017_1_genus.mds, p2017_1_numeric, na.rm = TRUE, permu = 1000)
env_ord
env_vector <- as.data.frame(scores(env_ord, display = "vectors"))
env_vector
env_vector <- env_vector[2:15,];env_vector<- env_vector[-3,]#keep only significant ordinations
env_vector
env_vector <- cbind(env_vector, Species = rownames(env_vector)); env_vector

#extract mds scores
p2017_1_scrs <- cbind(p2017_1_sample_data, data.frame(MDS1 = p2017_1_genus.mds$points[,1], MDS2 = p2017_1_genus.mds$points[,2])) 
p2017_1_scrs$grp <- grp
env_vector
p2017_1_scrs <- as.data.frame(p2017_1_scrs)
head(p2017_1_scrs)
#extract centeroid
p2017_merged_cent <- cbind(merged_p2017_sample_data, data.frame(MDS1 = p2017_merge_genus.mds$points[,1], MDS2 = p2017_merge_genus.mds$points[,2])) %>% 
  aggregate(cbind(MDS1, MDS2) ~ Days +  well_type + sulfate_level + nitrate_level,data = ., FUN = mean) 
#plot mds (not betadispersion and NOT CAP)
p2017_1_scrs
p2017_1_scrs$Days <- factor(p2017_1_scrs$Days, levels = c("-6" , "1" ,  "8" ,  "15" , "22" , "50" , "78" , "106" ,"134"))

ggplot() +
  #facet_wrap(~diet) + 
  geom_hline(yintercept = 0, color = "gray50", size = 0.5, alpha=0.6, linetype="dashed")+
  geom_vline(xintercept = 0, color ="gray50", size=0.5, alpha=0.6, linetype="dashed")+
  #stat_ellipse(data=p2017_merged_scrs, aes(x = MDS1, y = MDS2, fill = as.factor(phase_cat), color = as.factor(phase_cat)),  linetype = 2, geom = "polygon", alpha = 0.6) +# add ellipse
  #geom_point(data = p2017_merged_cent, size = 5) +            # centroids
  geom_point(data = p2017_1_scrs, aes(x = MDS1, y = MDS2,   shape=well_type),size=5.5, color = "black")+
  geom_point(data = p2017_1_scrs, aes(x = MDS1, y = MDS2,  color = as.factor(Days), shape=well_type),size=4) +                                              # sample scores
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
ggsave('~/2017_0.1_unifrac_nmds_ordination.png', bg = "transparent", width = 40, height = 30, units = "cm")
#end NMDS
#####

#merged_p2017 Canonicanl Analysis (CAP)
#####
View(p2017_1_metadata)
#multivariate anova to determine what percent of total variation is explained  by covariates
View(sample_data(p2017_1))
adonis2(p2017_1_genus_unifrac ~ Days+ well_type + nitrate_range + sulfate_range + acetate_presence + NH4_mgL, data = p2017_1_metadata, perm = 999)
p2017_1_CCA  <- ordinate(p2017_1, "CCA") #ordinate by CCA
p2017_1_CAP <-p2017_1_CCA$CA$u
plot_ordination(p2017_1, ordinate(p2017_1, "CCA"), "samples", color = "Days", shape = "well_type")
#Legendre, P. & Anderson, M.J. (1999). Distance-based redundancy analysis: testing multispecies responses in multifactorial ecological experiments. Ecological Monographs 69: 1-24.
#Anderson, M.J. & Willis, T.J. (2003). Canonical analysis of principal coordinates: a useful method of constrained ordination for ecology. Ecology 84: 511-525.
cap_genus <- capscale(p2017_1_genus_unifrac ~ Days * well_type +NO3_mgL + SO4_mgL +Acetate_uM + NH4_mgL, data = p2017_1_metadata)
cap_genus
plot(cap_genus)
anova(cap_genus)
cap_genus_cols <-cap_genus$CCA$wa
cap_genus_cols[,1:2]
#export vector as CCA
cap_genus_vector <-as.data.frame(cap_genus$CCA$biplot)
cap_genus_vector <- (cap_genus_vector/2)
cap_genus_vector$sites <-row.names(cap_genus_vector) ;cap_genus_vector 
row.names(cap_genus_vector) <-c("DaysDay106", "DaysDay134","DaysDay15", "DaysDay22", "DaysDay78", "DaysDay8",  "Monitoring Wells", "Nitrate", "Sulfate", "Acetate", "Ammonium", "Uranium" )
cap_genus_vector$sites <- row.names(cap_genus_vector) 
cap_genus_vector <- cap_genus_vector[8:12,]
cap_genus_vector 
#merge cap_genus data with metadata
cap_genus_metadata <- cbind(p2017_1_metadata, data.frame(CAP1 = cap_genus_cols[,1], CAP2=cap_genus_cols[,2])) 
cap_genus_metadata
head(cap_genus_metadata)
head(p2017_1_CCA_metadata)
plot(cap_genus)

#CAP PLOT
#set color levels in data
cap_genus_metadata$DaysPostinjection <- factor(cap_genus_metadata$DaysPostinjection, levels = c( "1" ,  "8" ,  "15" , "22" ,  "78" , "106" ,"134"))
cap_genus_metadata$Phase <- factor(cap_genus_metadata$phase_C12, levels = c("Control", "Phase0",  "Phase1", "Phase2", "Phase3" ))
cap_genus_metadata$type_day <- factor(cap_genus_metadata$type_day, levels = c("Control", "Monitoring Day-6",  "Monitoring Day1", "Monitoring Day8", "Monitoring Day15", "Monitoring Day22", "Monitoring Day50", "Monitoring Day78", "Monitoring Day106", "Monitoring Day134"))
cap_genus_metadata$phase_cat <- factor(cap_genus_metadata$phase_cat, levels = c("0", "1",  "2", "3"))

#plot CAP analysis
ggplot() +
  #facet_wrap(~diet) + 
  scale_colour_brewer(palette = "YlGnBu") + 
  scale_fill_brewer(palette = "YlGnBu") + 
  geom_hline(yintercept = 0, color = "gray50", size = 0.5, alpha=0.6, linetype="dashed")+
  geom_vline(xintercept = 0, color ="gray50", size=0.5, alpha=0.6, linetype="dashed")+
  #stat_ellipse(data=cap_genus_metadata, aes(x = CAP1, y = CAP2, fill = phase_cat),  linetype = 2, geom = "polygon", alpha = 0.6) +# add ellipse
  geom_segment(data = cap_genus_vector, # vector arrows
               aes(x = 0, xend = CAP1, y = 0, yend = CAP2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "gray35") + #scores
  #geom_point(data = p2017_merged_cent, size = 5) +            # centroids
  geom_point(data = cap_genus_metadata, aes(x = CAP1, y = CAP2,   shape=well_type),size=5.5, color = "black")+
  geom_point(data = cap_genus_metadata, aes(x = CAP1, y = CAP2,  color = DaysPostinjection, shape=well_type),size=4) +   
  geom_label(data = cap_genus_vector, aes(x = CAP1, y = CAP2, label = sites), # vector labels
             size = 3.3, check_overlap = TRUE) +# sample scores
  geom_text(aes(x=0.3, label="Significance of constraints:\n Sum Sq=3.73\n F=2.25\n p=0.001", y=0.47), colour="gray30", angle=0, text=element_text(size=7)) +
  labs(shape = "Well Type", color = "Day", fill="Injection Phase", y="CAP2 [10%]", x="CAP1 [23.4%]", title = "A. Canonical Analysis") + 
  coord_fixed() + #must have a fixed centroid
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="vertical",
        legend.position = "right",
        axis.ticks = element_line(color = "gray50"),
        #legend.position = c(0.46, 0.92),
        axis.title = element_text(size=25, color = "black"),
        plot.title = element_text(size=25, color = "black"),
        axis.text = element_text(size = 18, color = "black"),
        panel.border = element_rect(colour = "gray50", size=0.8, fill = "transparent"))                                       # same axis scaling
ggsave('~/2017_1_unifrac_CAP_ordination.jpeg', bg = "white", width = 20, height = 20, units = "cm")
#end CAP
#####

#p2017_1 anosim and adonis
#####
#anosim testing p.133 ADP-003 to test if composition differs between groups
#groups DO NOT differ in composition with well_type
anosim(p2017_1_genus_unifrac, (p2017_1_sample_data$well_type))
#groups DO!! differ in composition with NITRATE. R: 0.2202; p: 0.001
anosim(p2017_1_genus_unifrac, p2017_1_sample_data$nitrate_range)
#groups DO!! differ in composition with SULFATE. R: 0.3341; p: 0.001
anosim(p2017_1_genus_unifrac, p2017_1_sample_data$sulfate_range)
#groups DO!! differ in composition with DAYS. R: 0.5724; p: 0.001
anosim(p2017_1_genus_unifrac, p2017_1_sample_data$Days)
#groups DO!! differ in composition with URANIUM. R: 0.2528; p: 0.003
anosim(p2017_1_genus_unifrac, p2017_1_sample_data$uranium_range)
#groups DO!! differ in composition by acetate. R2: 0.104; p:0.03
anosim(p2017_1_genus_unifrac, p2017_1_sample_data$acetate_presence)
#ammonium
anosim(p2017_1_genus_unifrac, p2017_1_sample_data$ammonium_presence)

#adonis test (p.133-134 ADP-003)
head(sample_data(merged_p2017))
(microbes.adonis <- adonis2(p2017_1_genus_unifrac ~ well_type +Days +NH4_mgL+ nitrate_range +acetate_presence, as(sample_data(p2017_1), "data.frame"), 
                            permutations = 999, by=NULL))
microbes.adonis$R
(sulfate.adonis <- adonis2(merged_p2017_genus_unifrac ~ well_type_factor *Days, as(sample_data(merged_p2017), "data.frame"), 
                           permutations = 999, by=NULL))
(phase.adonis <- adonis2(merged_p2017_genus_unifrac ~ phase_cat, as(sample_data(merged_p2017), "data.frame"), 
                         permutations = 999, by=NULL))
(microbes.adonis <- adonis2(p2017_rarefied_genus_unifrac ~  nitrate_range * sulfate_range  *type_day , as(sample_data(merged_phylo2017_rarefied), "data.frame"), 
                            permutations = 10000, by="terms"))

plot(microbes.adonis)
#end adonis and anosim
#####

#betadispersion test 0.1 um (p.140 ADP-003)
#####
#test significant groups: sulfate
betadisper(p2017_1_genus_unifrac, sample_data(p2017_1)$sulfate_range)
anova(betadisper(p2017_1_genus_unifrac, sample_data(p2017_1)$sulfate_range))
plot(betadisper(p2017_1_genus_unifrac, sample_data(p2017_1)$sulfate_range))
#re-assign sulfate levels
sample_data(p2017_1)$sulfate_level <-as.factor(case_when((sample_data(p2017_1)$sulfate_range <= 20)              ~ "low",
                                                         (sample_data(p2017_1)$sulfate_range   > 20) ~ "high"))
betadisper(p2017_1_genus_unifrac, sample_data(p2017_1)$sulfate_level)
plot(betadisper(p2017_1_genus_unifrac, sample_data(p2017_1)$sulfate_level))
anova(betadisper(p2017_1_genus_unifrac, sample_data(p2017_1)$sulfate_level))
p2017_1_betadisp_sulfate <-plot(betadisper(p2017_1_genus_unifrac, sample_data(p2017_1)$sulfate_level))
#Days
betadisper(p2017_1_genus_unifrac, sample_data(p2017_1)$Days)
anova(betadisper(p2017_1_genus_unifrac, sample_data(p2017_1)$Days))
plot(betadisper(p2017_1_genus_unifrac, sample_data(p2017_1)$Days))
#re-assign days
head(sample_data(p2017_1))
sample_data(p2017_1)$day_level <-as.factor(case_when((sample_data(p2017_1)$DaysPostinjection <= 22)              ~ "Days 1 - 22",
                                                     (sample_data(p2017_1)$DaysPostinjection   > 22) ~ "Days 78 - 134"))
betadisper(p2017_1_genus_unifrac, sample_data(p2017_1)$day_level)
anova(betadisper(p2017_1_genus_unifrac, sample_data(p2017_1)$day_level))
p2017_1_betadisp_daylevel <-plot(betadisper(p2017_1_genus_unifrac, sample_data(p2017_1)$day_level))
#Acetate
betadisper(p2017_1_genus_unifrac, sample_data(p2017_1)$acetate_presence)
anova(betadisper(p2017_1_genus_unifrac, sample_data(p2017_1)$acetate_presence))
plot(betadisper(p2017_1_genus_unifrac, sample_data(p2017_1)$acetate_presence))
#Ammonium
sample_data(p2017_1)$NH4_mgL
betadisper(p2017_1_genus_unifrac, sample_data(p2017_1)$ammonium_presence)
anova(betadisper(p2017_1_genus_unifrac, sample_data(p2017_1)$ammonium_presence))
plot(betadisper(p2017_1_genus_unifrac, sample_data(p2017_1)$ammonium_presence))

#extract centroids and sites from DAY dispersion test (for plotting)
p2017_1_betadisp_daylevel_sites <- as.data.frame(p2017_1_betadisp_daylevel$sites) #see lines above for betadispersion testing
p2017_1_betadisp_daylevel_sites$label <- row.names(p2017_1_betadisp_daylevel_sites) #add label column
p2017_1_betadisp_daylevel_centroid <- as.data.frame(p2017_1_betadisp_daylevel$centroids)
p2017_1_betadisp_daylevel_centroid$label <- row.names(p2017_1_betadisp_daylevel_centroid)#add label column

#extract centroids and sites from sulfate dispersion test (for plotting)
p2017_1_betadisp_sulfate_sites <- as.data.frame(p2017_1_betadisp_sulfate$sites) #see lines above for betadispersion testing
p2017_1_betadisp_sulfate_sites$label <- row.names(p2017_1_betadisp_sulfate_sites) #add label column
p2017_1_betadisp_sulfate_centroid <- as.data.frame(p2017_1_betadisp_sulfate$centroids)
p2017_1_betadisp_sulfate_centroid$label <- row.names(p2017_1_betadisp_sulfate_centroid)#add label column

#make data frames for plotting - daylevel
p2017_1_betadisp_daylevel_sites <- cbind(p2017_1_sample_data, data.frame(PCoA1 = p2017_1_betadisp_daylevel$sites[,1], PCoA2 = p2017_1_betadisp_daylevel$sites[,2]))
p2017_1_betadisp_daylevel_sites $grp <-grp
p2017_1_betadisp_daylevel_sites$DaysPostinjection <- factor(p2017_1_betadisp_daylevel_sites$DaysPostinjection, levels = c("1", "8",  "15", "22", "78", "106", "134"))
p2017_1_betadisp_daylevel_sites$DaysPostinjection <- factor(p2017_1_betadisp_daylevel_sites$DaysPostinjection, levels = c( "1" ,  "8" ,  "15" , "22" ,  "78" , "106" ,"134"))

#make dataframe for plot - sulfate
p2017_1_betadisp_sulfate_sites <- cbind(p2017_1_sample_data, data.frame(PCoA1 = p2017_1_betadisp_sulfate$sites[,1], PCoA2 = p2017_1_betadisp_sulfate$sites[,2]))
p2017_1_betadisp_sulfate_sites $grp <-grp
p2017_1_betadisp_sulfate_sites$DaysPostinjection <- factor(p2017_1_betadisp_sulfate_sites$DaysPostinjection, levels = c("1", "8",  "15", "22", "78", "106", "134"))


head(p2017_1_betadisp_sulfate_sites)
head(p2017_1_betadisp_daylevel_centroid)

#plot dispersion eigen values # changes by day
ggplot(p2017_1_betadisp_daylevel_sites, aes(x = PCoA1, y = PCoA2,  color = day_level, shape = well_type)) +
  #facet_wrap(~diet) + 
  scale_colour_brewer(palette = "YlGnBu") + 
  scale_fill_brewer(palette = "YlGnBu") + 
  geom_hline(yintercept = 0, color = "gray50", size = 0.5, alpha=0.6, linetype="dashed")+
  geom_vline(xintercept = 0, color ="gray50", size=0.5, alpha=0.6, linetype="dashed")+
  stat_ellipse(data=p2017_1_betadisp_daylevel_sites, aes(x=PCoA1,y=PCoA2, fill = day_level), inherit.aes = FALSE, geom = "polygon", alpha = 0.6) +# add ellipse
  #stat_ellipse(data=p2017_1_betadisp_sulfate_sites, aes(x=PCoA1,y=PCoA2, color = sulfate_level), size = 1.2, inherit.aes = FALSE, linetype = 2,  alpha = 1) +# add ellipse
  geom_point(size=5.5, color = "black")+
  geom_text(aes(x=-0.4, label="Variance in group dispersion:\n Sum Sq=0.046\n F=6.25\n p<0.05", y=0.32), colour="gray30", angle=0, text=element_text(size=10)) +
  #geom_text(aes(x=0.28, label="Sulfate >20mg/L", y=-0.35), colour="#2c7fb8", angle=0, text=element_text(size=10)) +
  #geom_text(aes(x=-0.5, label="Sulfate <20 mg/L", y=-0.35), colour="#253494", angle=0, text=element_text(size=10)) +
  
  geom_point(size=4) +   
  geom_label(data = p2017_1_betadisp_daylevel_centroid, aes(x=PCoA1,y=PCoA2, label = label),inherit.aes = FALSE,size = 3, check_overlap = TRUE) +# sample scores
  #geom_point(data = merged_betadisp_type_day_centroid, aes(x = PCoA1, y = PCoA2),color = "blue",size = 6) +
  labs(shape = "Well Type", color = "", fill="", title = "B. Dispersion by day") + 
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="vertical",
        legend.position = "right",
        axis.ticks = element_line(color = "gray50"),
        #legend.position = c(0.46, 0.92),
        axis.title = element_text(size=25, color = "black"),
        plot.title = element_text(size=25, color = "black"),
        axis.text = element_text(size = 18, color = "black"),
        panel.border = element_rect(colour = "gray50", size=0.8, fill = "transparent"))                                       # same axis scaling
ggsave('~/2017_1_betadispersion_days.jpeg', bg = "white", width = 20, height = 17, units = "cm")

#changes by sulfate
ggplot(p2017_1_betadisp_sulfate_sites, aes(x = PCoA1, y = PCoA2,  color = DaysPostinjection, shape = well_type, size=SO4_mgL)) +
  #facet_wrap(~diet) + 
  scale_colour_brewer(palette = "YlGnBu") + 
  scale_fill_brewer(palette = "YlGnBu") + 
  geom_hline(yintercept = 0, color = "gray50", size = 0.5, alpha=0.6, linetype="dashed")+
  geom_vline(xintercept = 0, color ="gray50", size=0.5, alpha=0.6, linetype="dashed")+
  stat_ellipse(data=p2017_1_betadisp_sulfate_sites, aes(x=PCoA1,y=PCoA2, fill = sulfate_level, color = sulfate_level), inherit.aes = FALSE, linetype = 2, geom = "polygon", alpha = 0.6) +# add ellipse
  geom_point(color = "black")+
  geom_text(aes(x=-0.4, label="Variance in group dispersion:\n Sum Sq=0.046\n F=6.25\n p<0.05", y=0.32), inherit.aes = FALSE, colour="gray30", angle=0, text=element_text(size=10)) +
  geom_point() +   
  geom_label(data = p2017_1_betadisp_sulfate_centroid, aes(x=PCoA1,y=PCoA2, label = label),inherit.aes = FALSE,size = 3, check_overlap = TRUE) +# sample scores
  #geom_point(data = merged_betadisp_type_day_centroid, aes(x = PCoA1, y = PCoA2),color = "blue",size = 6) +
  labs(shape = "", color = "Days", fill="", title = "C. Group Dispersion by Sulfate ") + 
  scale_size(range = c(2.2, 6)) +
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = "bottom",
        axis.ticks = element_line(color = "gray50"),
        #legend.position = c(0.46, 0.92),
        axis.title = element_text(size=25, color = "black"),
        plot.title = element_text(size=25, color = "black"),
        axis.text = element_text(size = 18, color = "black"),
        panel.border = element_rect(colour = "gray50", size=0.8, fill = "transparent"))                                       # same axis scaling
ggsave('~/2017_betadispersion_phase.jpeg', bg = "white", width = 20, height = 17, units = "cm")


#testing unweighted distance and betadispersion in the overlapping community
ggplot(merged_betadisp_sites, aes(x = PCoA1, y = PCoA2,  color = as.factor(grp), shape = well_type)) +
  #facet_wrap(~diet) + 
  #scale_colour_brewer(palette = "YlGnBu") + 
  geom_hline(yintercept = 0, color = "gray50", size = 0.5, alpha=0.6, linetype="dashed")+
  geom_vline(xintercept = 0, color ="gray50", size=0.5, alpha=0.6, linetype="dashed")+
  stat_ellipse(aes(x=PCoA1,y=PCoA2, fill = as.factor(grp), color = as.factor(grp)),  linetype = 2, geom = "polygon", alpha = 0.6) +# add ellipse
  geom_point(size=5.5, color = "black")+
  geom_text(aes(x=0.35, label="Variance across groups:\n F=2.47\n p≤0.05", y=-0.32), colour="gray30", angle=0, text=element_text(size=8)) +
  geom_point(size=4) +                                              # sample scores
  #geom_point(data = merged_betadisp_type_day_centroid, aes(x = PCoA1, y = PCoA2),color = "blue",size = 6) +
  labs(shape = "Well Type", color = "", fill="Phase") + 
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = "bottom",
        axis.ticks = element_line(color = "gray50"),
        #legend.position = c(0.46, 0.92),
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        panel.border = element_rect(colour = "gray50", size=0.8, fill = "transparent"))                                        # same axis scaling
ggsave("~/2017_betadispersion_type_day.png", bg = "transparent", width = 60, height = 35, units = "cm")

plot(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$phase_cat))
disp_p12 <- plot(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$phase_C12))
disp_p134 <- plot(betadisper(merged_p2017_genus_unifrac, sample_data(merged_p2017)$phase_C134))
pdf(disp_p12,"~/2017_p2017_disp_p12.pdf")
pdf(disp_p134,"~/2017_p2017_disp_p134.pdf")


#end beta diversity testing - merged
#####


#Mantel test - statistical correlations between 
install.packages("geosphere")
library(geosphere)
#OTU distance matrix
merged_p2017_genus_unifrac
#Environmental distance matrix
Nitrate_merged <- merged_p2017_sample_df$Nitrate
N_euclidean <- dist(Nitrate_merged, method="euclidean")
Days_merged <- merged_p2017_sample_df$Days
N_euclidean <- dist(Days_merged, method="euclidean")
p2017_euclidean <- dist(merged_p2017_numeric, method="euclidean")
#perform mantel test
N_merge_mantel <- vegan::mantel(merged_p2017_genus_unifrac, N_euclidean, method="spearman", permutations = 999, na.rm = TRUE);N_merge_mantel
N_merge_bray_mantel <- vegan::mantel(merged_p2017_genus_bray, N_euclidean, method="spearman", permutations = 999, na.rm = TRUE);N_merge_bray_mantel
p2017_bray_mantel <- vegan::mantel(merged_p2017_genus_bray, p2017_euclidean, method="spearman", permutations = 999, na.rm = TRUE);p2017_bray_mantel
p2017_unifrac_mantel <- vegan::mantel(merged_p2017_genus_unifrac, p2017_euclidean, method="spearman", permutations = 999, na.rm = TRUE);p2017_bray_mantel




#combined 0.1 and 0.2 beta diversity
#####
# microbial community differentiation based on unifrac distance
merged_p2017_sample_data <-sample_data(merged_phylo2017_rarefied)
set.seed(0451)
p2017_rarefied_genus_unifrac = UniFrac(tax_glom(merged_phylo2017_rarefied, taxrank = "Genus"))
p2017_rarefied_genus_unifrac
set.seed(0451)
p2017_merge_genus.mds <- metaMDS(p2017_rarefied_genus_unifrac, zerodist=ignore, try=500 ,trymax = 1000)
p2017_merge_genus.mds
plot(p2017_merge_genus.mds)
ASV_tree <- hclust(p2017_rarefied_genus_unifrac, method = "average")
plot(ASV_tree)
grp <- cutree(ASV_tree, k = 4)

#extract mds scores
p2017_merged_scrs <- cbind(merged_p2017_sample_data, data.frame(MDS1 = p2017_merge_genus.mds$points[,1], MDS2 = p2017_merge_genus.mds$points[,2])) 
p2017_merged_scrs$grp <- grp
head(p2017_merged_scrs)
#extract centeroid
p2017_merged_cent <- cbind(merged_p2017_sample_data, data.frame(MDS1 = p2017_merge_genus.mds$points[,1], MDS2 = p2017_merge_genus.mds$points[,2])) %>% 
  aggregate(cbind(MDS1, MDS2) ~ DaysPostinjection + DayCategory + well_type_factor + uranium_range +uranium_range + sulfate_range + nitrate_range,data = ., FUN = mean) 
#anosim 
#groups DO NOT differ in composition with well_type
anosim(p2017_rarefied_genus_unifrac, merged_p2017_sample_data$well_type)
#groups DO!! differ in composition with NITRATE. R: 0.2202; p: 0.001
anosim(p2017_rarefied_genus_unifrac, merged_p2017_sample_data$nitrate_range)
#groups DO!! differ in composition with SULFATE. R: 0.3341; p: 0.001
anosim(p2017_rarefied_genus_unifrac, merged_p2017_sample_data$sulfate_range)
#groups DO!! differ in composition with DAYS. R: 0.5724; p: 0.001
anosim(p2017_rarefied_genus_unifrac, merged_p2017_sample_data$DayCategory)
#groups DO!! differ in composition with URANIUM. R: 0.2528; p: 0.003
anosim(p2017_rarefied_genus_unifrac, merged_p2017_sample_data$uranium_range)

cbn <- combn(x=unique(metadata$FilterSize), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  ps.subs <- phyloseq::subset_samples(merged_p2017, phase_cat %in% cbn[,i])
  metadata_sub <- data.frame(phyloseq::sample_data(ps.subs))
  permanova_pairwise <- anosim(phyloseq::distance(ps.subs, method="unifrac", binary = TRUE), 
                               metadata_sub$FilterSize)
  p <- c(p, permanova_pairwise$signif[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(t(cbn), p=p, p.adj=p.adj)

p.table
#      1     2     p p.adj
#1 0.1um 0.2um 0.001 0.001

# ANOSIM by Day
metadata <- data.frame(sample_data(phylo2017))
anosim(dist_17_jacc, metadata$Days)

cbn <- combn(x=unique(metadata$Days), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  ps.subs <- subset_samples(phylo2017, Days %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  permanova_pairwise <- anosim(phyloseq::distance(ps.subs, method="jaccard", binary = TRUE), 
                               metadata_sub$Days)
  p <- c(p, permanova_pairwise$signif[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(t(cbn), p=p, p.adj=p.adj)
# ANOSIM results 
p.table <- as.data.frame(p.table)
sig_by_day_ANOSIM <- subset(p.table, p.adj < 0.01)
sig_by_day_ANOSIM

#Repeat distance and ANOSIM with 0.2 um data
ASV_table <- as.data.frame(t(otu_2017_2))
row.names(ASV_table)
colnames(ASV_table)


#plot mds
ggplot(p2017_merged_scrs, aes(x = MDS1, y = MDS2,  color = DayCategory, shape = well_type_factor)) +
  #facet_wrap(~diet) + 
  scale_colour_brewer(palette = "YlGnBu") + 
  geom_hline(yintercept = 0, color = "gray50", size = 0.5, alpha=0.6, linetype="dashed")+
  geom_vline(xintercept = 0, color ="gray50", size=0.5, alpha=0.6, linetype="dashed")+
  #geom_point(data = p2017_merged_cent, size = 5) +            # centroids
  geom_point(size=4, color = "black")+
  geom_point(size=3) +                                              # sample scores
  labs(shape = "Well Type", color = "Days") + 
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = "bottom",
        axis.ticks = element_line(color = "gray50"),
        #legend.position = c(0.46, 0.92),
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        panel.border = element_rect(colour = "gray50", size=0.8, fill = "transparent"))                                           # same axis scaling
ggsave('~/2017_EVO_unifrac_nmds.png', bg = "transparent", width = 40, height = 30, units = "cm")


(microbes.adonis <- adonis2(p2017_rarefied_genus_unifrac ~ well_type_factor *DayCategory * nitrate_range +sulfate_range + type_day, as(sample_data(merged_phylo2017_rarefied), "data.frame"), 
                            permutations = 10000, by=NULL))
microbes.adonis
(microbes.adonis <- adonis2(p2017_rarefied_genus_unifrac ~  nitrate_range * sulfate_range  *type_day , as(sample_data(merged_phylo2017_rarefied), "data.frame"), 
                            permutations = 10000, by="terms"))

plot(microbes.adonis)
anova(betadisper(p2017_rarefied_genus_unifrac, sample_data(merged_phylo2017_rarefied)$type_day))
p2017_betadisper_rarefied_genus_unifrac<- betadisper(p2017_rarefied_genus_unifrac, sample_data(merged_phylo2017_rarefied)$type_day)
p2017_betadisper_rarefied_genus_unifrac
plot(betadisper(p2017_rarefied_genus_unifrac, sample_data(merged_phylo2017_rarefied)$type_day))
ggsave("~/2017_p2017_diff_abundance_betadispersion_type_day.png", bg = "transparent", width = 60, height = 35, units = "cm")

#permutation tests of dispersion within the monitoring wells # indicates a significant community shift in monitoring wells related to days
merged_p2017_genus.monitoring <- tax_glom(subset_samples(merged_phylo2017_rarefied, well_type == "Monitoring"), taxrank = "Genus")  
merged_p2017_genus_unifrac.monitoring = UniFrac(merged_p2017_genus.monitoring)
(merged_p2017_genus.monitoring.adonis <- adonis2(merged_p2017_genus_unifrac.monitoring ~nitrate_range * DayCategory*sulfate_range * uranium_range , as(sample_data(merged_p2017_genus.monitoring), "data.frame"), permutations = 10000, by="term"))
#permutation tests of dispersion within the CONTROL well # indicates a significant community shift in control wells related to days
merged_p2017_genus.control <- tax_glom(subset_samples(merged_phylo2017_rarefied, well_type == "Control"), taxrank = "Genus")  
merged_p2017_genus_unifrac.control = UniFrac(merged_p2017_genus.control)
(merged_p2017_genus.control.adonis <- adonis2(merged_p2017_genus_unifrac.control ~  Days * nitrate_range , as(sample_data(merged_p2017_genus.control), "data.frame"), permutations = 10000, strata = nitrate_reduction_cat))

#end of combined 0.1 and 0.2
#####

#0.1 beta diversity
#####
# microbial community differentiation based on unifrac distance
p2017_1 <- prune_taxa(taxa_sums(p2017_1) > 0, p2017_1)
p2017_1_sample_data <-sample_data(p2017_1)
p2017_1_metadata <- p2017_1_sample_data 
set.seed(0451)
p2017_1_genus_unifrac = UniFrac(tax_glom(p2017_1, taxrank = "Genus"))
p2017_1_genus_unifrac
p2017_1_order_unifrac = UniFrac(tax_glom(p2017_1, taxrank = "Order"))
p2017_1_order_unifrac
set.seed(0451)
p2017_1_genus.mds <- metaMDS(p2017_1_genus_unifrac, zerodist=ignore, try=500 ,trymax = 1000)
p2017_1_genus.mds
plot(p2017_1_genus.mds)
p2017_1_order.mds <- metaMDS(p2017_1_order_unifrac, zerodist=ignore, try=500 ,trymax = 1000)
p2017_1_order.mds
plot(p2017_1_order.mds)
ASV_tree <- hclust(p2017_1_genus_unifrac, method = "average")
plot(ASV_tree)
grp <- cutree(ASV_tree, k = 4)

#extract mds scores
p2017_1_merged_scrs <- cbind(p2017_1_sample_data, data.frame(MDS1 = p2017_1_genus.mds$points[,1], MDS2 = p2017_1_genus.mds$points[,2])) 
p2017_1_merged_scrs$grp <- grp
head(p2017_1_merged_scrs)
p2017_1_sample_data[is.na(p2017_1_sample_data)] <- 0
View(p2017_1_sample_data)
#extract centeroid
p2017_merged_cent <- cbind(merged_p2017_sample_data, data.frame(MDS1 = p2017_merge_genus.mds$points[,1], MDS2 = p2017_merge_genus.mds$points[,2])) %>% 
  aggregate(cbind(MDS1, MDS2) ~ DaysPostinjection + DayCategory + well_type_factor + uranium_range +uranium_range + sulfate_range + nitrate_range,data = ., FUN = mean) 
#groups DO NOT differ in composition with well_type
anosim(p2017_1_genus_unifrac, p2017_1_sample_data$well_type)
#groups DO!! differ in composition with NITRATE. R: 0.2757; p: 0.004
anosim(p2017_1_genus_unifrac, p2017_1_sample_data$nitrate_range)
#groups DO!! differ in composition with SULFATE. R: 0.2875; p: 0.002
anosim(p2017_1_genus_unifrac, p2017_1_sample_data$sulfate_range)
#groups DO!! differ in composition with DAYS. R: 0.597 p: 0.001
anosim(p2017_1_genus_unifrac, p2017_1_sample_data$Days)
#groups DO!! differ in composition with URANIUM. R: 0.2528; p: 0.003
anosim(p2017_rarefied_genus_unifrac, merged_p2017_sample_data$uranium_range)

#plot mds
ggplot(p2017_1_merged_scrs, aes(x = MDS1, y = MDS2,  color = as.factor(sulfate_range), shape = well_type)) +
  #facet_wrap(~diet) + 
  scale_colour_brewer(palette = "YlGnBu") + 
  geom_hline(yintercept = 0, color = "gray50", size = 0.5, alpha=0.6, linetype="dashed")+
  geom_vline(xintercept = 0, color ="gray50", size=0.5, alpha=0.6, linetype="dashed")+
  #geom_point(data = p2017_merged_cent, size = 5) +            # centroids
  geom_point(size=4, color = "black")+
  geom_point(size=3) +                                              # sample scores
  labs(shape = "Well Type", color = "Days") + 
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="horizontal",
        legend.position = "bottom",
        axis.ticks = element_line(color = "gray50"),
        #legend.position = c(0.46, 0.92),
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        panel.border = element_rect(colour = "gray50", size=0.8, fill = "transparent"))                                           # same axis scaling
ggsave('~/2017_EVO_unifrac_nmds.png', bg = "transparent", width = 40, height = 30, units = "cm")


(microbes.adonis <- adonis2(p2017_rarefied_genus_unifrac ~ well_type_factor *DayCategory * nitrate_range +sulfate_range + type_day, as(sample_data(merged_phylo2017_rarefied), "data.frame"), 
                            permutations = 10000, by=NULL))
microbes.adonis
(microbes.adonis <- adonis2(p2017_rarefied_genus_unifrac ~  nitrate_range * sulfate_range  *type_day , as(sample_data(merged_phylo2017_rarefied), "data.frame"), 
                            permutations = 10000, by="terms"))

plot(microbes.adonis)
anova(betadisper(p2017_rarefied_genus_unifrac, sample_data(merged_phylo2017_rarefied)$type_day))
p2017_betadisper_rarefied_genus_unifrac<- betadisper(p2017_rarefied_genus_unifrac, sample_data(merged_phylo2017_rarefied)$type_day)
p2017_betadisper_rarefied_genus_unifrac
plot(betadisper(p2017_rarefied_genus_unifrac, sample_data(merged_phylo2017_rarefied)$type_day))
ggsave("~/2017_p2017_diff_abundance_betadispersion_type_day.png", bg = "transparent", width = 60, height = 35, units = "cm")

#permutation tests of dispersion within the monitoring wells # indicates a significant community shift in monitoring wells related to days
merged_p2017_genus.monitoring <- tax_glom(subset_samples(merged_phylo2017_rarefied, well_type == "Monitoring"), taxrank = "Genus")  
merged_p2017_genus_unifrac.monitoring = UniFrac(merged_p2017_genus.monitoring)
(merged_p2017_genus.monitoring.adonis <- adonis2(merged_p2017_genus_unifrac.monitoring ~nitrate_range * DayCategory*sulfate_range * uranium_range , as(sample_data(merged_p2017_genus.monitoring), "data.frame"), permutations = 10000, by="term"))
#permutation tests of dispersion within the CONTROL well # indicates a significant community shift in control wells related to days
merged_p2017_genus.control <- tax_glom(subset_samples(merged_phylo2017_rarefied, well_type == "Control"), taxrank = "Genus")  
merged_p2017_genus_unifrac.control = UniFrac(merged_p2017_genus.control)
(merged_p2017_genus.control.adonis <- adonis2(merged_p2017_genus_unifrac.control ~  Days * nitrate_range , as(sample_data(merged_p2017_genus.control), "data.frame"), permutations = 10000, strata = nitrate_reduction_cat))

#end of 0.2 only
#####

# 0.1 Canonicanl Analysis (CAP)
#####
#multivariate anova to determine what percent of total variation is explained  by covariates
adonis2(p2017_1_genus_unifrac ~ Days + well_type + nitrate_level + sulfate_level, data = p2017_1_sample_data, perm = 9999)
p2017_merged_CCA  <- ordinate(merged_p2017, "CCA") #ordinate by CCA
p2017_merged_CAP <-p2017_merged_CCA$CA$u
p2017_merged_CCA_metadata; plot_ordination(merged_p2017, ordinate(merged_p2017, "CCA"), "samples", color = "Days", shape = "well_type")
#Legendre, P. & Anderson, M.J. (1999). Distance-based redundancy analysis: testing multispecies responses in multifactorial ecological experiments. Ecological Monographs 69: 1-24.
#Anderson, M.J. & Willis, T.J. (2003). Canonical analysis of principal coordinates: a useful method of constrained ordination for ecology. Ecology 84: 511-525.
cap_genus <- capscale(merged_p2017_genus_unifrac ~ well_type +Nitrate + Sulfate + Uranium + Acetate + Ammonium, data = merged_p2017_metadata)
cap_genus
plot(cap_genus)
anova(cap_genus)
cap_genus_cols <-cap_genus$CCA$wa
cap_genus_cols[,1:2]
#export vector as CCA
cap_genus_vector <-as.data.frame(cap_genus$CCA$biplot)
cap_genus_vector <- (cap_genus_vector/2)
cap_genus_vector$sites <-row.names(cap_genus_vector);cap_genus_vector
#merge cap_genus data with metadata
cap_genus_metadata <- cbind(merged_p2017_metadata, data.frame(CAP1 = cap_genus_cols[,1], CAP2=cap_genus_cols[,2])) 
cap_genus_metadata
head(cap_genus_metadata)
head(p2017_merged_CCA_metadata)
plot(cap_genus)

#CAP PLOT
#set color levels in data
cap_genus_metadata$Days <- factor(cap_genus_metadata$Days, levels = c("-6" , "1" ,  "8" ,  "15" , "22" , "50" , "78" , "106" ,"134"))
cap_genus_metadata$Phase <- factor(cap_genus_metadata$phase_C12, levels = c("Control", "Phase0",  "Phase1", "Phase2", "Phase3" ))
cap_genus_metadata$type_day <- factor(cap_genus_metadata$type_day, levels = c("Control", "Monitoring Day-6",  "Monitoring Day1", "Monitoring Day8", "Monitoring Day15", "Monitoring Day22", "Monitoring Day50", "Monitoring Day78", "Monitoring Day106", "Monitoring Day134"))
cap_genus_metadata$phase_cat <- factor(cap_genus_metadata$phase_cat, levels = c("0", "1",  "2", "3"))

#plot CAP analysis
ggplot() +
  #facet_wrap(~diet) + 
  scale_colour_brewer(palette = "YlGnBu") + 
  scale_fill_brewer(palette = "YlGnBu") + 
  geom_hline(yintercept = 0, color = "gray50", size = 0.5, alpha=0.6, linetype="dashed")+
  geom_vline(xintercept = 0, color ="gray50", size=0.5, alpha=0.6, linetype="dashed")+
  #stat_ellipse(data=cap_genus_metadata, aes(x = CAP1, y = CAP2, fill = phase_cat),  linetype = 2, geom = "polygon", alpha = 0.6) +# add ellipse
  geom_segment(data = cap_genus_vector, # vector arrows
               aes(x = 0, xend = CAP1, y = 0, yend = CAP2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "gray35") + #scores
  #geom_point(data = p2017_merged_cent, size = 5) +            # centroids
  geom_point(data = cap_genus_metadata, aes(x = CAP1, y = CAP2,   shape=well_type),size=5.5, color = "black")+
  geom_point(data = cap_genus_metadata, aes(x = CAP1, y = CAP2,  color = as.factor(type_day), shape=well_type),size=4) +   
  geom_label(data = cap_genus_vector, aes(x = CAP1, y = CAP2, label = sites), # vector labels
             size = 3, check_overlap = TRUE) +# sample scores
  geom_text(aes(x=0.25, label="Significance of constraints:\n F=3.43\n p≤0.001", y=-0.4), colour="gray30", angle=0, text=element_text(size=7)) +
  labs(shape = "Well Type", color = "Well and Day", fill="Injection Phase", y="CAP2 [6.21%]", x="CAP1 [19.46%]") + 
  coord_fixed() + #must have a fixed centroid
  theme(panel.background = element_rect(fill = "white"), 
        legend.background = element_rect(fill="white"),
        legend.direction="vertical",
        legend.position = "right",
        axis.ticks = element_line(color = "gray50"),
        #legend.position = c(0.46, 0.92),
        axis.title = element_text(size=15, color = "black"),
        axis.text = element_text(size = 12.5, color = "black"),
        panel.border = element_rect(colour = "gray50", size=0.8, fill = "transparent"))                                           # same axis scaling
ggsave('~/2017_unifrac_CAP_ordination.png', bg = "transparent", width = 40, height = 30, units = "cm")
#end CAP
#####






#
#13. differential analysis using ggdifclade 
# differential analysis using linear discriminate analysis - 0.2 µm community
######
# mergedby well-type
set.seed(0451)# the seed should be set for reproducibly results.
normTSS <- function(physeq)
{
  aux <- data.frame(sample_data(physeq))
  aux$"NF.TSS" <- sample_sums(physeq)#/exp(mean(log(sample_sums(physeq))))
  sample_data(physeq) <- aux
  physeq
}
remotes::install_github("yiluheihei/microbiomeMarker")
library(microbiomeMarker)
merged_p2017
p2017_merged_TSS <- normTSS(merged_p2017)
p2017_2_TSS <- normalize(merged_p2017, method = "TSS")
head(sample_data(p2017_2_TSS))
sample_data(p2017_2_TSS)
deres_merged_WT <- diff_analysis(obj = p2017_merged_TSS, classgroup = "well_type",
                                 mlfun = "lda", filtermod = "pvalue", firstcomfun = "wilcox.test",
                                 firstalpha = 0.05, strictmod = TRUE, secondcomfun = "wilcox.test",
                                 subclmin = 3, subclwilc = TRUE, secondalpha = 0.01, lda=3)
deres_rf_merged_well <- diff_analysis(obj = p2017_merged_TSS, classgroup = "WellID",
                                      mlfun = "rf", filtermod = "pvalue", firstcomfun = "wilcox.test",
                                      firstalpha = 0.05, strictmod = TRUE, secondcomfun = "wilcox.test",
                                      subclmin = 3, subclwilc = TRUE, secondalpha = 0.01, lda=3)
# print out results of the diff_analysis                     
deres_2_WT

#boxplot
diffbox_merged_deres <- ggdiffbox(obj=deres_merged_WT, box_notch=FALSE, 
                                  colorlist=c("darkgoldenrod4", "slateblue"), l_xlabtext="relative abundance") 
diffbox_merged_deres
ggsave("~/2017_welltype_color_es_boxplot.png", 
       bg = "transparent", width = 60, height = 35, units = "cm")
#plot the linear discriminate analysis effect size
WT_merged_effect <- ggeffectsize(obj=deres_merged_WT, 
                                 lineheight=0.1, linewidth=0.2) + 
  scale_color_manual(values=c("darkgoldenrod4", "slateblue")) 
WT_merged_effect
ggsave("~/2017_welltype_color_es_boxplot.png", 
       bg = "transparent", width = 60, height = 35, units = "cm")

#0.2 by nitrate
set.seed(0451)
#p2017_2_filtered_transform <- transform(p2017_2_filtered, 'log10p') with sulfate reduction
colnames(sample_data(p2017_2_TSS))
deres_2_SR <- diff_analysis(obj = p2017_2_filtered_transform, classgroup = "sulfate_reduction",
                            mlfun = "lda", filtermod = "pvalue", firstcomfun = "wilcox.test",
                            firstalpha = 0.05, strictmod = TRUE, secondcomfun = "wilcox.test",
                            subclmin = 3, subclwilc = TRUE, secondalpha = 0.01, lda=3.5)#modified from 3 to 3.5
deres_2_SR# print out results of the diff_analysis   
#boxplot
diffbox_2SR_deres <- ggdiffbox(obj=deres_2_SR, box_notch=FALSE,colorlist=c("darkgoldenrod4", "slateblue"), l_xlabtext="relative abundance") 
diffbox_2SR_deres
#plot the linear discriminate analysis effect size
WT_2SR_effect <- ggeffectsize(obj=deres_2_SR, lineheight=0.1, linewidth=0.2) + scale_color_manual(values=c("darkgoldenrod4", "slateblue")) 
WT_2SR_effect
ggsave("~/2017_filtersize_color_es_boxplot.png", 
       bg = "transparent", width = 60, height = 35, units = "cm")

#0.2 by nitrate
head(sample_data(merged_p2017))

deres_2_NR <- diff_analysis(obj = merged_p2017, classgroup = "phase_C134",
                            mlfun = "lda", filtermod = "pvalue", firstcomfun = "kruskal.test",
                            firstalpha = 0.05, strictmod = TRUE, secondcomfun = "kruskal.test",
                            subclmin = 3, subclwilc = TRUE, secondalpha = 0.01, lda=3)#modified from 3 to 3.5
deres_2_NR# print out results of the diff_analysis   
#boxplot
diffbox_2NR_deres <- ggdiffbox(obj=deres_2_NR, box_notch=FALSE,colorlist=c("darkgoldenrod4", "slateblue", "orange", "brown"), l_xlabtext="relative abundance") ;diffbox_2NR_deres
#plot the linear discriminate analysis effect size
NR_2_effect <- ggeffectsize(obj=deres_2_NR, lineheight=0.1, linewidth=0.2) + scale_color_manual(values=c("darkgoldenrod4", "slateblue"));NR_2_effect
ggsave("~/2017_filtersize_color_es_boxplot.png", 
       bg = "transparent", width = 60, height = 35, units = "cm")

#0.2 by uranium
head(sample_data(p2017_2_filtered_transform))
deres_2u <- diff_analysis(obj = p2017_2_filtered_transform, classgroup = "uranium_level",
                          mlfun = "lda", filtermod = "pvalue", firstcomfun = "kruskal.test",
                          firstalpha = 0.05, strictmod = TRUE, secondcomfun = "wilcox.test",
                          subclmin = 3, subclwilc = TRUE, secondalpha = 0.01, lda=3.5)#modified from 3 to 3.5
deres_2u# print out results of the diff_analysis   
#boxplot
diffbox_2u_deres <- ggdiffbox(obj=deres_2u, box_notch=FALSE,colorlist=c("darkgoldenrod4", "slateblue"), l_xlabtext="relative abundance") ;diffbox_2u_deres
#plot the linear discriminate analysis effect size
WT_2u_effect <- ggeffectsize(obj=deres_2u, lineheight=0.1, linewidth=0.2) + scale_color_manual(values=c("darkgoldenrod4", "slateblue"));WT_2u_effect
ggsave("~/2017_filtersize_color_es_boxplot.png", 
       bg = "transparent", width = 60, height = 35, units = "cm")

#end of differential abundance for 0.2 by geochem categories 
#####

# differential analysis using linear discriminate analysis - merged community
######
# merged by well-type
set.seed(0451)# the seed should be set for reproducibly results.
normTSS <- function(physeq)
{
  aux <- data.frame(sample_data(physeq))
  aux$"NF.TSS" <- sample_sums(physeq)#/exp(mean(log(sample_sums(physeq))))
  sample_data(physeq) <- aux
  physeq
}
remotes::install_github("yiluheihei/microbiomeMarker")
library(microbiomeMarker)
p2017_merged_CLR <- merged_p2017
p2017_merged_CLR %>% tax_transform("clr", rank = NA)
p2017_merged_CLR
p2017_merged_TSS = transform_sample_counts(merged_p2017, function(x) x/sum(x))
head(sample_data(p2017_2_TSS))
sample_data(p2017_2_TSS)
deres_merged_WT_TSS <- diff_analysis(obj = p2017_merged_TSS, classgroup = "well_type",
                                     mlfun = "lda", filtermod = "pvalue", firstcomfun = "wilcox.test",
                                     firstalpha = 0.05, strictmod = TRUE, secondcomfun = "wilcox.test",
                                     subclmin = 3, subclwilc = TRUE, secondalpha = 0.01, lda=3)
deres_1_filtered_WT_CLR <- diff_analysis(obj = p2017_1_filtered_CLR, classgroup = "well_type",
                                         mlfun = "lda", filtermod = "pvalue", firstcomfun = "wilcox.test",
                                         firstalpha = 0.05, strictmod = FALSE, secondcomfun = "wilcox.test",
                                         subclmin = 3, subclwilc = TRUE, secondalpha = 0.05, lda=3)
deres_1_filtered_WT_rare <- diff_analysis(obj = p2017_combined_0.1_filtered, classgroup = "well_type",
                                          mlfun = "lda", filtermod = "pvalue", firstcomfun = "wilcox.test",
                                          firstalpha = 0.05, strictmod = TRUE, secondcomfun = "wilcox.test",
                                          subclmin = 3, subclwilc = TRUE, secondalpha = 0.05, lda=3)


# print out results of the diff_analysis                     
deres_2_WT
deres_merged_WT_TSS

#boxplot
diffbox_merged_deres <- ggdiffbox(obj=deres_merged_WT_TSS, box_notch=FALSE, 
                                  colorlist=c("orange", "slateblue"), l_xlabtext="Relative Abundance") 
diffbox_merged_deres
ggsave("~/2017_LDA+ABUN_kruskal_welltype_merged.tiff", 
       bg = "white", width = 35, height = 35, units = "cm")
#plot the linear discriminate analysis effect size
WT_merged_effect <- ggeffectsize(obj=deres_merged_WT_TSS, 
                                 lineheight=0.1, linewidth=0.2) + 
  scale_color_manual(values=c("orange", "slateblue")) 
WT_merged_effect

ggsave("~/2017_LDA_kruskal_welltype_merged.tiff", 
       bg = "white", width = 35, height = 65, units = "cm")

#sulfate
set.seed(0451)
#p2017_2_filtered_transform <- transform(p2017_2_filtered, 'log10p') with sulfate reduction
colnames(sample_data(p2017_merged_TSS))
deres_merged_SO4 <- diff_analysis(obj = p2017_merged_TSS, classgroup = "sulfate_reduction",
                                  mlfun = "lda", filtermod = "pvalue", firstcomfun = "wilcox.test",
                                  firstalpha = 0.05, strictmod = TRUE, secondcomfun = "wilcox.test",
                                  subclmin = 3, subclwilc = TRUE, secondalpha = 0.01, lda=3)#modified from 3 to 3.5
View(deres_merged_SO4)# print out results of the diff_analysis   
so4_complete_result <- deres_merged_SO4@result
deres_merged_SO4

#boxplot
diffbox_merged_SO4_deres <- ggdiffbox(obj=deres_merged_SO4, box_notch=FALSE,colorlist=c("orange", "slateblue"), l_xlabtext="Relative Abundance") 
diffbox_merged_SO4_deres
ggsave("~/2017_LDA+ABUN_wilcox_sulfate_merged.tiff", 
       bg = "transparent", width = 35, height = 75, units = "cm")
#plot the linear discriminate analysis effect size
WT_merged_SO4_effect <- ggeffectsize(obj=deres_merged_SO4, lineheight=0.1, linewidth=0.2) + scale_color_manual(values=c("orange", "slateblue")) 
WT_merged_SO4_effect
ggsave("~/2017_LDA_wilcox_sulfate_merged.tiff", 
       bg = "transparent", width = 35, height = 85, units = "cm")

#nitrate
head(sample_data(p2017_2_filtered_transform))
deres_merged_NR <- diff_analysis(obj = p2017_merged_TSS, classgroup = "nitrate_reduction",
                                 mlfun = "lda", filtermod = "pvalue", firstcomfun = "wilcox.test",
                                 firstalpha = 0.05, strictmod = TRUE, secondcomfun = "wilcox.test",
                                 subclmin = 3, subclwilc = TRUE, secondalpha = 0.01, lda=3)#modified from 3 to 3.5
deres_merged_NR# print out results of the diff_analysis   
#boxplot
diffbox_merged_NR_deres <- ggdiffbox(obj=deres_merged_NR, box_notch=FALSE,colorlist=c("orange", "slateblue"), l_xlabtext="Relative Abundance") ;diffbox_merged_NR_deres
ggsave("~/2017_LDA+ABUN_wilcox_nitrate_merged.tiff", 
       bg = "transparent", width = 35, height = 75, units = "cm")
#plot the linear discriminate analysis effect size
merged_NR_effect <- ggeffectsize(obj=deres_merged_NR, lineheight=0.1, linewidth=0.2) + scale_color_manual(values=c("orange", "slateblue"));merged_NR_effect
ggsave("~/2017_LDA_wilcox_nitrate_merged.tiff", 
       bg = "white", width = 35, height = 85, units = "cm")

# acetate
deres_merged_acetate <- diff_analysis(obj = p2017_merged_TSS, classgroup = "acetate_presence",
                                      mlfun = "lda", filtermod = "pvalue", firstcomfun = "wilcox.test",
                                      firstalpha = 0.05, strictmod = TRUE, secondcomfun = "wilcox.test",
                                      subclmin = 3, subclwilc = TRUE, secondalpha = 0.01, lda=3)#modified from 3 to 3.5
deres_merged_acetate
#boxplot
diffbox_merged_acetate <- ggdiffbox(obj=deres_merged_acetate, box_notch=FALSE,colorlist=c("orange", "slateblue"), l_xlabtext="relative abundance") ;diffbox_merged_acetate
ggsave("~/2017_LDA+ABUN_wilcox_acetate_merged.tiff", 
       bg = "transparent", width = 35, height = 75, units = "cm")
#plot the linear discriminate analysis effect size
merged_acetate_effect <- ggeffectsize(obj=deres_merged_acetate, lineheight=0.1, linewidth=0.2) + scale_color_manual(values=c("orange", "slateblue"));merged_acetate_effect
ggsave("~/2017_LDA_wilcox_acetate_merged.tiff", 
       bg = "white", width = 35, height = 95, units = "cm")

#phase
head(sample_data(merged_p2017))

deres_merged_phase <- diff_analysis(obj = p2017_merged_TSS, classgroup = "phase_C134",
                                    mlfun = "lda", filtermod = "pvalue", firstcomfun = "kruskal.test",
                                    firstalpha = 0.05, strictmod = TRUE, secondcomfun = "kruskal.test",
                                    subclmin = 3, subclwilc = TRUE, secondalpha = 0.01, lda=3)#modified from 3 to 3.5
deres_merged_phase# print out results of the diff_analysis   
#boxplot
diffbox_merged_phase_deres <- ggdiffbox(obj=deres_merged_phase, box_notch=FALSE,colorlist=c("darkgoldenrod4", "slateblue", "orange", "brown"), l_xlabtext="relative abundance") ;diffbox_merged_phase_deres
ggsave("~/2017_LDA+ABUN_kruskal_phase_merged.tiff", 
       bg = "transparent", width = 35, height = 45, units = "cm")
#plot the linear discriminate analysis effect size
merged_phase_effect <- ggeffectsize(obj=deres_merged_phase, lineheight=0.1, linewidth=0.2) + scale_color_manual(values=c("orange", "brown"));merged_phase_effect
ggsave("~/2017_LDA_kruskal_phase_merged.tiff", 
       bg = "transparent", width = 35, height = 45, units = "cm")

#end of differential abundance for complete community by geochem categories 
#####

# differential analysis using linear discriminate analysis - 0.1 OVERLAP
######
# merged by well-type
set.seed(0451)# the seed should be set for reproducibly results.
normTSS <- function(physeq)
{
  aux <- data.frame(sample_data(physeq))
  aux$"NF.TSS" <- sample_sums(physeq)#/exp(mean(log(sample_sums(physeq))))
  sample_data(physeq) <- aux
  physeq
}
remotes::install_github("yiluheihei/microbiomeMarker")
library(microbiomeMarker)
p2017_1_filtered_CLR <- p2017_combined_0.1_filtered
p2017_1_filtered_CLR %>% tax_transform("clr", rank = NA)
p2017_1_filtered_CLR
sample_data(p2017_combined_0.1_filtered)
sample_data(p2017_combined_0.1_filtered)$day_level <-as.factor(case_when((sample_data(p2017_combined_0.1_filtered)$Days <= 22)              ~ "Days 1 - 22",
                                                                         (sample_data(p2017_combined_0.1_filtered)$Days   > 22) ~ "Days 78 - 134"))

p2017_1_filtered_TSS = transform_sample_counts(p2017_combined_0.1_filtered, function(x) x/sum(x))
head(sample_data(p2017_1_filtered_TSS))

deres_1_filtered_WT_TSS <- diff_analysis(obj = p2017_1_filtered_TSS, classgroup = "well_type",
                                         mlfun = "lda", filtermod = "pvalue", firstcomfun = "wilcox.test",
                                         firstalpha = 0.05, strictmod = TRUE, secondcomfun = "wilcox.test",
                                         subclmin = 3, subclwilc = TRUE, secondalpha = 0.01, lda=3)
deres_1_filtered_WT_CLR <- diff_analysis(obj = p2017_1_filtered_CLR, classgroup = "well_type",
                                         mlfun = "lda", filtermod = "pvalue", firstcomfun = "wilcox.test",
                                         firstalpha = 0.05, strictmod = TRUE, secondcomfun = "wilcox.test",
                                         subclmin = 3, subclwilc = TRUE, secondalpha = 0.05, lda=2)
deres_1_filtered_WT_rare <- diff_analysis(obj = p2017_combined_0.1_filtered, classgroup = "well_type",
                                          mlfun = "lda", filtermod = "pvalue", firstcomfun = "wilcox.test",
                                          firstalpha = 0.05, strictmod = TRUE, secondcomfun = "wilcox.test",
                                          subclmin = 3, subclwilc = TRUE, secondalpha = 0.05, lda=2)

# print out results of the diff_analysis                     
deres_1_filtered_WT_TSS
deres_1_filtered_WT_CLR
deres_1_filtered_WT_rare

#boxplot
diffbox_1_filtered_deres <- ggdiffbox(obj=deres_1_filtered_WT_TSS, box_notch=FALSE, 
                                      colorlist=c("darkgoldenrod4", "slateblue"), l_xlabtext="Relative Abundance") 
diffbox_1_filtered_deres
ggsave("~/2017_LDA+ABUN_wilcox_welltype_0.1_filtered.tiff", 
       bg = "white", width = 35, height = 35, units = "cm")
#plot the linear discriminate analysis effect size
WT_1_filtered_effect <- ggeffectsize(obj=deres_1_filtered_WT, 
                                     lineheight=0.1, linewidth=0.2) + 
  scale_color_manual(values=c("darkgoldenrod4", "slateblue")) 
WT_1_filtered_effect

ggsave("~/2017_LDA_wilcox_welltype_0.1_filtered.tiff", 
       bg = "white", width = 35, height = 35, units = "cm")

#sulfate
set.seed(0451)
#p2017_2_filtered_transform <- transform(p2017_2_filtered, 'log10p') with sulfate reduction
deres_0.1_filtered_SO4 <- diff_analysis(obj = p2017_1_filtered_TSS, classgroup = "sulfate_reduction",
                                        mlfun = "lda", filtermod = "pvalue", firstcomfun = "wilcox.test",
                                        firstalpha = 0.05, strictmod = TRUE, secondcomfun = "wilcox.test",
                                        subclmin = 3, subclwilc = TRUE, secondalpha = 0.01, lda=3)#modified from 3 to 3.5
deres_1_filtered_SO4_TSS <- diff_analysis(obj = p2017_1_filtered_TSS, classgroup = "sulfate_reduction",
                                          mlfun = "lda", filtermod = "pvalue", firstcomfun = "wilcox.test",
                                          firstalpha = 0.05, strictmod = TRUE, secondcomfun = "wilcox.test",
                                          subclmin = 3, subclwilc = TRUE, secondalpha = 0.01, lda=2)
deres_1_filtered_SO4_CLR <- diff_analysis(obj = p2017_1_filtered_CLR, classgroup = "sulfate_reduction",
                                          mlfun = "lda", filtermod = "pvalue", firstcomfun = "wilcox.test",
                                          firstalpha = 0.05, strictmod = TRUE, secondcomfun = "wilcox.test",
                                          subclmin = 3, subclwilc = TRUE, secondalpha = 0.01, lda=2)
deres_1_filtered_SO4_rare <- diff_analysis(obj = p2017_combined_0.1_filtered, classgroup = "sulfate_reduction",
                                           mlfun = "lda", filtermod = "pvalue", firstcomfun = "wilcox.test",
                                           firstalpha = 0.05, strictmod = TRUE, secondcomfun = "wilcox.test",
                                           subclmin = 3, subclwilc = TRUE, secondalpha = 0.01, lda=2)
deres_0.1_filtered_SO4
deres_1_filtered_SO4_TSS
deres_1_filtered_SO4_CLR
deres_1_filtered_SO4_rare# print out results of the diff_analysis   
#boxplot
diffbox_0.1_filtered_SO4_TSS <- ggdiffbox(obj=deres_0.1_filtered_SO4, box_notch=FALSE,colorlist=c("orange", "slateblue"), l_xlabtext="Relative Abundance") 
diffbox_0.1_filtered_SO4_TSS
diffbox_0.1_filtered_SO4_CLR <- ggdiffbox(obj=deres_1_filtered_SO4_CLR, box_notch=FALSE,colorlist=c("darkgoldenrod4", "slateblue"), l_xlabtext="Relative Abundance") 
diffbox_0.1_filtered_SO4_CLR
diffbox_0.1_filtered_SO4_rare <- ggdiffbox(obj=deres_1_filtered_SO4_rare, box_notch=FALSE,colorlist=c("darkgoldenrod4", "slateblue"), l_xlabtext="Relative Abundance") 
diffbox_0.1_filtered_SO4_rare

ggsave("~/2017_LDA+ABUN_wilcox_sulfate_0.1_filtered.tiff", 
       bg = "white", width = 35, height = 35, units = "cm")
#plot the linear discriminate analysis effect size
filtered_0.1_SO4_effect <- ggeffectsize(obj=deres_0.1_filtered_SO4, lineheight=0.1, linewidth=0.2) + scale_color_manual(values=c("orange", "slateblue")) 
filtered_0.1_SO4_effect
ggsave("~/2017_LDA_wilcox_sulfate_0.1_filtered.tiff", 
       bg = "transparent", width = 35, height = 35, units = "cm")

#nitrate
head(sample_data(p2017_2_filtered_transform))
deres_0.1_filtered_NR <- diff_analysis(obj = p2017_1_filtered_TSS, classgroup = "nitrate_reduction",
                                       mlfun = "lda", filtermod = "pvalue", firstcomfun = "wilcox.test",
                                       firstalpha = 0.05, strictmod = TRUE, secondcomfun = "wilcox.test",
                                       subclmin = 3, subclwilc = TRUE, secondalpha = 0.01, lda=3)#modified from 3 to 3.5
deres_1_filtered_NR_TSS <- diff_analysis(obj = p2017_1_filtered_TSS, classgroup = "nitrate_reduction",
                                         mlfun = "lda", filtermod = "pvalue", firstcomfun = "wilcox.test",
                                         firstalpha = 0.05, strictmod = TRUE, secondcomfun = "wilcox.test",
                                         subclmin = 3, subclwilc = TRUE, secondalpha = 0.01, lda=2)
deres_1_filtered_NR_CLR <- diff_analysis(obj = p2017_1_filtered_CLR, classgroup = "nitrate_reduction",
                                         mlfun = "lda", filtermod = "pvalue", firstcomfun = "wilcox.test",
                                         firstalpha = 0.05, strictmod = TRUE, secondcomfun = "wilcox.test",
                                         subclmin = 3, subclwilc = TRUE, secondalpha = 0.01, lda=2)
deres_1_filtered_NR_rare <- diff_analysis(obj = p2017_combined_0.1_filtered, classgroup = "nitrate_reduction",
                                          mlfun = "lda", filtermod = "pvalue", firstcomfun = "wilcox.test",
                                          firstalpha = 0.05, strictmod = TRUE, secondcomfun = "wilcox.test",
                                          subclmin = 3, subclwilc = TRUE, secondalpha = 0.01, lda=2)

deres_1_filtered_NR_TSS
deres_1_filtered_NR_CLR
deres_1_filtered_NR_rare
deres_0.1_filtered_NR# print out results of the diff_analysis   
#boxplot
diffbox_0.1_filtered_NR_TSS <- ggdiffbox(obj=deres_0.1_filtered_NR, box_notch=FALSE,colorlist=c("orange", "slateblue"), l_xlabtext="Relative Abundance") ;diffbox_0.1_filtered_NR_TSS
diffbox_0.1_filtered_NR_CLR <- ggdiffbox(obj=deres_1_filtered_NR_CLR, box_notch=FALSE,colorlist=c("darkgoldenrod4", "slateblue"), l_xlabtext="Relative Abundance") ;diffbox_0.1_filtered_NR_CLR
diffbox_0.1_filtered_NR_rare <- ggdiffbox(obj=deres_1_filtered_NR_rare, box_notch=FALSE,colorlist=c("darkgoldenrod4", "slateblue"), l_xlabtext="Relative Abundance") ;diffbox_0.1_filtered_NR_rare
ggsave("~/2017_LDA+ABUN_wilcox_nitrate_0.1_filtered.tiff", 
       bg = "transparent", width = 35, height = 35, units = "cm")
#plot the linear discriminate analysis effect size
filtered_0.1_NR_effect <- ggeffectsize(obj=deres_0.1_filtered_NR, lineheight=0.1, linewidth=0.2) + scale_color_manual(values=c("orange", "slateblue"));filtered_0.1_NR_effect
ggsave("~/2017_LDA_wilcox_nitrate_0.1_filtered.tiff", 
       bg = "white", width = 35, height = 35, units = "cm")

# acetate
deres_0.1_filtered_acetate <- diff_analysis(obj = p2017_1_filtered_TSS, classgroup = "acetate_presence",
                                            mlfun = "lda", filtermod = "pvalue", firstcomfun = "wilcox.test",
                                            firstalpha = 0.05, strictmod = TRUE, secondcomfun = "wilcox.test",
                                            subclmin = 3, subclwilc = TRUE, secondalpha = 0.01, lda=2)#modified from 3 to 3.5

#boxplot
diffbox_0.1_filtered_acetate <- ggdiffbox(obj=deres_0.1_filtered_acetate, box_notch=FALSE,colorlist=c("slateblue" ,"orange"), l_xlabtext="relative abundance") ;diffbox_0.1_filtered_acetate
ggsave("~/2017_LDA+ABUN_wilcox_acetate_0.1_filtered.tiff", 
       bg = "transparent", width = 35, height = 85, units = "cm")
#plot the linear discriminate analysis effect size
filtered_0.1_acetate_effect <- ggeffectsize(obj=deres_0.1_filtered_acetate, lineheight=0.1, linewidth=0.2) + scale_color_manual(values=c("orange","slateblue" ));filtered_0.1_acetate_effect
ggsave("~/2017_LDA_wilcox_acetate_0.1_filtered.tiff", 
       bg = "white", width = 35, height = 95, units = "cm")

#days 
deres_1_day <- diff_analysis(obj = p2017_1_filtered_TSS, classgroup = "day_level",
                             mlfun = "lda", filtermod = "pvalue", firstcomfun = "wilcox.test",
                             firstalpha = 0.05, strictmod = TRUE, secondcomfun = "wilcox.test",
                             subclmin = 3, subclwilc = TRUE, secondalpha = 0.01, lda=2)#modified from 3 to 3.5
deres_1_day# print out results of the diff_analysis   
#boxplot
diffbox_1_day_deres <- ggdiffbox(obj=deres_1_day, box_notch=FALSE,colorlist=c("orange", "slateblue","darkgoldenrod4",  "brown"), l_xlabtext="relative abundance") ;diffbox_1_day_deres
ggsave("~/2017_LDA+ABUN_wilcox_day_0.1_filtered.tiff", 
       bg = "transparent", width = 35, height = 45, units = "cm")
#plot the linear discriminate analysis effect size
day_effect <- ggeffectsize(obj=deres_1_day, lineheight=0.1, linewidth=0.2) + scale_color_manual(values=c("orange", "slateblue"));day_effect
ggsave("~/2017_LDA_wilcox_day_0.1_filtered.tiff", 
       bg = "transparent", width = 35, height = 75, units = "cm")

#end of differential abundance for 0.1 OVERLAP bacteria by geochem categories 
#####

# differential analysis using linear discriminate analysis - 0.1 EXCLUSIVE
######
# merged by well-type
set.seed(0451)# the seed should be set for reproducibly results.
#Add acetate category
sample_data(p2017_1_filtered)$acetate_presence <-as.factor(case_when((sample_data(p2017_1_filtered)$Acetate <= 0.1)              ~ "≤0.1",
                                                                     (sample_data(p2017_1_filtered)$Acetate   > 0.1) ~ ">0.1"))
sample_data(p2017_1_filtered)$day_level <-as.factor(case_when((sample_data(p2017_1_filtered)$Days <= 22)              ~ "Days 1 - 22",
                                                              (sample_data(p2017_1_filtered)$Days   > 22) ~ "Days 78 - 134"))


library(microViz)
p2017_1_E_CLR = tax_transform(p2017_1_filtered, trans = "clr")
p2017_1_E_TSS = transform_sample_counts(p2017_1_filtered, function(x) x/sum(x))
head(sample_data(p2017_1_E_TSS))
deres_1_E_WT <- diff_analysis(obj = p2017_1_E_TSS, classgroup = "well_type",
                              mlfun = "lda", filtermod = "pvalue", firstcomfun = "wilcox.test",
                              firstalpha = 0.05, strictmod = TRUE, secondcomfun = "wilcox.test",
                              subclmin = 3, subclwilc = TRUE, secondalpha = 0.01, lda=3)


# print out results of the diff_analysis                     
deres_1_E_WT

#boxplot
diffbox_1_E_deres <- ggdiffbox(obj=deres_1_E_WT, box_notch=FALSE, 
                               colorlist=c("darkgoldenrod4", "slateblue"), l_xlabtext="Relative Abundance") 
diffbox_1_E_deres
ggsave("~/2017_LDA+ABUN_kruskal_welltype_0.1_E.tiff", 
       bg = "white", width = 35, height = 35, units = "cm")
#plot the linear discriminate analysis effect size
WT_1_E_effect <- ggeffectsize(obj=deres_1_E_WT, lineheight=0.1, linewidth=0.2) + 
  scale_color_manual(values=c("orange", "slateblue")) 
WT_1_E_effect

ggsave("~/2017_LDA_kruskal_welltype_0.1_E.tiff", 
       bg = "white", width = 35, height = 35, units = "cm")

#sulfate
set.seed(0451)
#p2017_2_E_transform <- transform(p2017_2_E, 'log10p') with sulfate reduction
deres_0.1_E_SO4 <- diff_analysis(obj = p2017_1_E_TSS, classgroup = "sulfate_reduction",
                                 mlfun = "lda", filtermod = "pvalue", firstcomfun = "wilcox.test",
                                 firstalpha = 0.05, strictmod = TRUE, secondcomfun = "wilcox.test",
                                 subclmin = 2, subclwilc = TRUE, secondalpha = 0.01, lda=2)#modified from 3 to 3.5
deres_0.1_E_SO4# print out results of the diff_analysis   
#boxplot
diffbox_0.1_E_SO4_deres <- ggdiffbox(obj=deres_0.1_E_SO4, box_notch=FALSE,colorlist=c("darkgoldenrod4", "slateblue"), l_xlabtext="Relative Abundance") 
diffbox_0.1_E_SO4_deres
ggsave("~/2017_LDA+ABUN_kruskal_sulfate_0.1_E.tiff", 
       bg = "white", width = 35, height = 35, units = "cm")
#plot the linear discriminate analysis effect size
E_0.1_SO4_effect <- ggeffectsize(obj=deres_0.1_E_SO4, lineheight=0.1, linewidth=0.2) + scale_color_manual(values=c("orange", "slateblue")) 
E_0.1_SO4_effect
ggsave("~/2017_LDA_kruskal_sulfate_0.1_E.tiff", 
       bg = "transparent", width = 35, height = 35, units = "cm")

#nitrate
head(sample_data(p2017_2_E_transform))
deres_0.1_E_NR <- diff_analysis(obj = p2017_1_E_TSS, classgroup = "nitrate_reduction",
                                mlfun = "lda", filtermod = "pvalue", firstcomfun = "kruskal.test",
                                firstalpha = 0.05, strictmod = TRUE, secondcomfun = "kruskal.test",
                                subclmin = 2, subclwilc = TRUE, secondalpha = 0.05, lda=2)#modified from 3 to 3.5
deres_0.1_E_NR# print out results of the diff_analysis   
#boxplot
diffbox_0.1_E_NR_deres <- ggdiffbox(obj=deres_0.1_E_NR, box_notch=FALSE,colorlist=c("darkgoldenrod4", "slateblue"), l_xlabtext="Relative Abundance") ;diffbox_0.1_E_NR_deres
ggsave("~/2017_LDA+ABUN_kruskal_nitrate_0.1_E.tiff", 
       bg = "transparent", width = 35, height = 35, units = "cm")
#plot the linear discriminate analysis effect size
E_0.1_NR_effect <- ggeffectsize(obj=deres_0.1_E_NR, lineheight=0.1, linewidth=0.2) + scale_color_manual(values=c("darkgoldenrod4", "slateblue"));E_0.1_NR_effect
ggsave("~/2017_LDA_kruskal_nitrate_0.1_E.tiff", 
       bg = "white", width = 35, height = 35, units = "cm")

# acetate
deres_0.1_E_acetate <- diff_analysis(obj = p2017_1_E_TSS, classgroup = "acetate_presence",
                                     mlfun = "lda", filtermod = "pvalue", firstcomfun = "kruskal.test",
                                     firstalpha = 0.05, strictmod = TRUE, secondcomfun = "kruskal.test",
                                     subclmin = 2, subclwilc = TRUE, secondalpha = 0.05, lda=2)#modified from 3 to 3.5
deres_0.1_E_acetate

#boxplot
diffbox_0.1_E_acetate <- ggdiffbox(obj=deres_0.1_E_acetate, box_notch=FALSE,colorlist=c("darkgoldenrod4", "slateblue"), l_xlabtext="relative abundance") ;diffbox_0.1_E_acetate
ggsave("~/2017_LDA+ABUN_kruskal_acetate_0.1_E.tiff", 
       bg = "transparent", width = 35, height = 35, units = "cm")
#plot the linear discriminate analysis effect size
E_0.1_acetate_effect <- ggeffectsize(obj=deres_0.1_E_acetate, lineheight=0.1, linewidth=0.2) + scale_color_manual(values=c("darkgoldenrod4", "slateblue"));E_0.1_acetate_effect
ggsave("~/2017_LDA_kruskal_acetate_0.1_E.tiff", 
       bg = "white", width = 35, height = 95, units = "cm")
# day_level
deres_0.1_E_day_level <- diff_analysis(obj = p2017_1_E_TSS, classgroup = "day_level",
                                       mlfun = "lda", filtermod = "pvalue", firstcomfun = "wilcox.test",
                                       firstalpha = 0.05, strictmod = FALSE, secondcomfun = "wilcox.test",
                                       subclmin = 2, subclwilc = FALSE, secondalpha = 0.05, lda=1)#modified from 3 to 3.5
deres_0.1_E_day_level

#boxplot
diffbox_0.1_E_day_level <- ggdiffbox(obj=deres_0.1_E_day_level, box_notch=FALSE,colorlist=c("darkgoldenrod4", "slateblue"), l_xlabtext="relative abundance") ;diffbox_0.1_E_day_level
ggsave("~/2017_LDA+ABUN_kruskal_day_level_0.1_E.tiff", 
       bg = "transparent", width = 35, height = 35, units = "cm")
#plot the linear discriminate analysis effect size
E_0.1_day_level_effect <- ggeffectsize(obj=deres_0.1_E_day_level, lineheight=0.1, linewidth=0.2) + scale_color_manual(values=c("darkgoldenrod4", "slateblue"));E_0.1_day_level_effect
ggsave("~/2017_LDA_kruskal_day_level_0.1_E.tiff", 
       bg = "white", width = 35, height = 95, units = "cm")

#end of differential abundance for 0.1 E bacteria by geochem categories 
#####

# 14. functional data from 16S using Tax2Fun
#example data
#####
data(dataset)
dataset
library(microeco)
# load the example data; 16S rRNA gene amplicon sequencing dataset
# metadata table (this is categorical data, but can include all data)
data(sample_info_16S)
view(sample_info_16S)
class(sample_info_16S)
# feature table
data(otu_table_16S)
otu_table_16S[1:5, 1:5]
otu_2017[1:5, 1:5]
class(otu_2017)
# taxonomic assignment table
data(taxonomy_table_16S)
taxonomy_table_16S[1:5, 1:5]
view(taxonomy_table_16S)
taxa_2017[1:5,1:5]
# make the taxonomic information unified, very important
taxa_2017 %<>% tidy_taxonomy
# use phylogenetic tree to calculate phylogeny-based alpha and beta metrics
data(phylo_tree_16S)
class(phylo_tree_16S)
class(tree_2017)
# load the environmental data table if it is not in sample table
data(env_data_16S)
class(env_data_16S)
View(env_data_16S)
meta_2017
class(meta_2017)
# Let's create a microtable object with more information
head(sample_data(merged_p2017))
micro2017 <- microtable$new(sample_table = meta_2017, otu_table = otu_2017, tax_table = taxa_2017, phylo_tree = tree_2017)
#end test data
#####

#funtion of phase 1 complete community
#####
#phase 1 function
m2017_phase1 <- subset_samples(merged_p2017, phase_cat == "1") #subset phyloseq
m2017_phase1 <- prune_taxa(taxa_sums(m2017_phase1) > 0, m2017_phase1)
otu_phase1 <-as.data.frame(t(otu_table(m2017_phase1)))
tax_phase1 <- as.data.frame(phyloseq::tax_table(m2017_phase1)); colnames(tax_phase1) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
meta_phase1 <- as.data.frame(phyloseq::sample_data(m2017_phase1))
tree_phase1 = rtree(ntaxa(m2017_phase1), rooted=TRUE, tip.label=taxa_names(m2017_phase1))
merged_phase1 <- microtable$new(sample_table = meta_phase1, otu_table = otu_phase1, tax_table = tax_phase1, phylo_tree = tree_phase1)

dataset
# create object of trans_func
t2 <- trans_func$new(dataset) #example
t2017 <- trans_func$new(micro2017)
t2017_phase1 <- trans_func$new(merged_phase1)
t2017
# mapping the taxonomy to the database
# this can recognize prokaryotes or fungi automatically if the names of taxonomic levels are standard.
# for fungi example, see https://chiliubio.github.io/microeco_tutorial/other-dataset.html#fungi-data
# default database for prokaryotes is FAPROTAX database
t2$cal_spe_func(prok_database = "FAPROTAX")
t2$res_spe_func[1:5, 1:2]

t2017$cal_spe_func(prok_database = "FAPROTAX")
t2017$res_spe_func[1:5, 1:2]

t2017_phase1$cal_spe_func(prok_database = "FAPROTAX")
t2017_phase1$res_spe_func[1:5, 1:2]
# calculate the percentages for communities
# here do not consider the abundance
t2$cal_spe_func_perc(abundance_weighted = FALSE)
t2$res_spe_func_perc[1:5, 1:2]

t2017$cal_spe_func_perc(abundance_weighted = FALSE)
t2017$res_spe_func_perc[1:5, 1:2]

t2017_phase1$cal_spe_func_perc(abundance_weighted = FALSE)
t2017_phase1$res_spe_func_perc[1:5, 1:2]
t2017_phase1$cal_spe_func_perc(abundance_weighted = TRUE)
t2017_phase1$res_spe_func_perc[1:5, 1:2]

# construct a network for the example
network <- trans_network$new(dataset = micro2017, cal_cor = "base", taxa_level = "OTU", filter_thres = 0.0001, cor_method = "spearman")
network$cal_network(p_thres = 0.01, COR_cut = 0.7)
network$cal_module()

network_merge2017_phase1 <- trans_network$new(dataset = merged_phase1, cal_cor = "base", taxa_level = "OTU", filter_thres = 0.0001, cor_method = "spearman")
network_merge2017_phase1$cal_network(p_thres = 0.01, COR_cut = 0.7)
network_merge2017_phase1$cal_module()
# convert module info to microtable object
meco_phase1_module <- network_merge2017_phase1$trans_comm(use_col = "module")
meco_phase1_module_func <- trans_func$new(meco_phase1_module)
meco_phase1_module_func$cal_spe_func(prok_database = "FAPROTAX")
meco_phase1_module_func$cal_spe_func_perc(abundance_weighted = FALSE)
meco_phase1_module_func$plot_spe_func_perc()
meco_phase1_module_func$cal_spe_func_perc(abundance_weighted = TRUE)
meco_phase1_module_func$plot_spe_func_perc()
jpeg("~/phase1_complete_function.tiff",width=3200, height=2000, res=400, bg="white")
plot.new()
meco_phase1_module_func$plot_spe_func_perc()
getOption("device")
dev.set(which = dev.next())
dev.off()

#end of phase 1 complete community function
#####

#Phase 2 complete community function
#####
m2017_phase2 <- subset_samples(merged_p2017, phase_cat == "2") #subset phyloseq
m2017_phase2 <- prune_taxa(taxa_sums(m2017_phase2) > 0, m2017_phase2)
otu_phase2 <-as.data.frame(t(otu_table(m2017_phase2)))
tax_phase2 <- as.data.frame(phyloseq::tax_table(m2017_phase2)); colnames(tax_phase2) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
meta_phase2 <- as.data.frame(phyloseq::sample_data(m2017_phase2))
tree_phase2 = rtree(ntaxa(m2017_phase2), rooted=TRUE, tip.label=taxa_names(m2017_phase2))
merged_phase2 <- microtable$new(sample_table = meta_phase2, otu_table = otu_phase2, tax_table = tax_phase2, phylo_tree = tree_phase2)

# create object of trans_func
t2017_phase2 <- trans_func$new(merged_phase2)
# mapping the taxonomy to the database
# this can recognize prokaryotes or fungi automatically if the names of taxonomic levels are standard.
# for fungi example, see https://chiliubio.github.io/microeco_tutorial/other-dataset.html#fungi-data
# default database for prokaryotes is FAPROTAX database

t2017_phase2$cal_spe_func(prok_database = "FAPROTAX")
t2017_phase2$res_spe_func[1:5, 1:2]
# calculate the percentages for communities
# here do not consider the abundance

t2017_phase2$cal_spe_func_perc(abundance_weighted = FALSE)
t2017_phase2$res_spe_func_perc[1:5, 1:2]
t2017_phase2$cal_spe_func_perc(abundance_weighted = TRUE)
t2017_phase2$res_spe_func_perc[1:5, 1:2]

# construct a network
network_merge2017_phase2 <- trans_network$new(dataset = merged_phase2, cal_cor = "base", taxa_level = "OTU", filter_thres = 0.0001, cor_method = "spearman")
network_merge2017_phase2$cal_network(p_thres = 0.01, COR_cut = 0.7)
network_merge2017_phase2$cal_module()
# convert module info to microtable object
meco_phase2_module <- network_merge2017_phase2$trans_comm(use_col = "module")
meco_phase2_module_func <- trans_func$new(meco_phase2_module)
meco_phase2_module_func$cal_spe_func(prok_database = "FAPROTAX")
meco_phase2_module_func$cal_spe_func_perc(abundance_weighted = FALSE)
meco_phase2_module_func$plot_spe_func_perc()
meco_phase2_module_func$cal_spe_func_perc(abundance_weighted = TRUE)
meco_phase2_module_func$plot_spe_func_perc()
jpeg("~/phase2_complete_function.tiff",width=3200, height=2000, res=400, bg="white")
plot.new()
meco_phase2_module_func$plot_spe_func_perc()
getOption("device")
dev.set(which = dev.next())
dev.off()

#end complete community phase2 function
#####

#phase3 complete community function
#####
m2017_phase3 <- subset_samples(merged_p2017, phase_cat == "3") #subset phyloseq
m2017_phase3 <- prune_taxa(taxa_sums(m2017_phase3) > 0, m2017_phase3)
otu_phase3 <-as.data.frame(t(otu_table(m2017_phase3)))
tax_phase3 <- as.data.frame(phyloseq::tax_table(m2017_phase3)); colnames(tax_phase3) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
meta_phase3 <- as.data.frame(phyloseq::sample_data(m2017_phase3))
tree_phase3 = rtree(ntaxa(m2017_phase3), rooted=TRUE, tip.label=taxa_names(m2017_phase3))
merged_phase3 <- microtable$new(sample_table = meta_phase3, otu_table = otu_phase3, tax_table = tax_phase3, phylo_tree = tree_phase3)

# create object of trans_func
t2017_phase3 <- trans_func$new(merged_phase3)
# mapping the taxonomy to the database
# this can recognize prokaryotes or fungi automatically if the names of taxonomic levels are standard.
# for fungi example, see https://chiliubio.github.io/microeco_tutorial/other-dataset.html#fungi-data
# default database for prokaryotes is FAPROTAX database

t2017_phase3$cal_spe_func(prok_database = "FAPROTAX")
t2017_phase3$res_spe_func[1:5, 1:2]
# calculate the percentages for communities
# here do not consider the abundance

t2017_phase3$cal_spe_func_perc(abundance_weighted = FALSE)
t2017_phase3$res_spe_func_perc[1:5, 1:2]
t2017_phase3$cal_spe_func_perc(abundance_weighted = TRUE)
t2017_phase3$res_spe_func_perc[1:5, 1:2]

# construct a network 
network_merge2017_phase3 <- trans_network$new(dataset = merged_phase3, cal_cor = "base", taxa_level = "OTU", filter_thres = 0.0001, cor_method = "spearman")
network_merge2017_phase3$cal_network(p_thres = 0.01, COR_cut = 0.7)
network_merge2017_phase3$cal_module()
# convert module info to microtable object
meco_phase3_module <- network_merge2017_phase3$trans_comm(use_col = "module")
meco_phase3_module_func <- trans_func$new(meco_phase3_module)
meco_phase3_module_func$cal_spe_func(prok_database = "FAPROTAX")
meco_phase3_module_func$cal_spe_func_perc(abundance_weighted = FALSE)
meco_phase3_module_func$plot_spe_func_perc()
meco_phase3_module_func$cal_spe_func_perc(abundance_weighted = TRUE)
meco_phase3_module_func$plot_spe_func_perc()
jpeg("~/phase3_complete_function.tiff",width=3200, height=2000, res=400, bg="white")
plot.new()
meco_phase3_module_func$plot_spe_func_perc()
getOption("device")
dev.set(which = dev.next())
dev.off()
#end phase 3 function
#####

#control phase function
#####
m2017_phase0 <- subset_samples(merged_p2017, phase_cat == "0") #subset phyloseq
m2017_phase0 <- prune_taxa(taxa_sums(m2017_phase0) > 0, m2017_phase0)
otu_phase0 <-as.data.frame(t(otu_table(m2017_phase0)))
tax_phase0 <- as.data.frame(phyloseq::tax_table(m2017_phase0)); colnames(tax_phase0) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
meta_phase0 <- as.data.frame(phyloseq::sample_data(m2017_phase0))
tree_phase0 = rtree(ntaxa(m2017_phase0), rooted=TRUE, tip.label=taxa_names(m2017_phase0))
merged_phase0 <- microtable$new(sample_table = meta_phase0, otu_table = otu_phase0, tax_table = tax_phase0, phylo_tree = tree_phase0)

# create object of trans_func
t2017_phase0 <- trans_func$new(merged_phase0)
# mapping the taxonomy to the database
# this can recognize prokaryotes or fungi automatically if the names of taxonomic levels are standard.
# for fungi example, see https://chiliubio.github.io/microeco_tutorial/other-dataset.html#fungi-data
# default database for prokaryotes is FAPROTAX database

t2017_phase0$cal_spe_func(prok_database = "FAPROTAX")
t2017_phase0$res_spe_func[1:5, 1:2]
# calculate the percentages for communities
# here do not consider the abundance

#t2017_phase3$cal_spe_func_perc(abundance_weighted = FALSE)
#t2017_phase3$res_spe_func_perc[1:5, 1:2]
t2017_phase0$cal_spe_func_perc(abundance_weighted = TRUE)
t2017_phase0$res_spe_func_perc[1:5, 1:2]

# construct a network 
network_merge2017_phase0 <- trans_network$new(dataset = merged_phase0, cal_cor = "base", taxa_level = "OTU", filter_thres = 0.0001, cor_method = "spearman")
network_merge2017_phase0$cal_network(p_thres = 0.01, COR_cut = 0.7)
network_merge2017_phase0$cal_module()
# convert module info to microtable object
meco_phase0_module <- network_merge2017_phase0$trans_comm(use_col = "module")
meco_phase0_module_func <- trans_func$new(meco_phase0_module)
meco_phase0_module_func$cal_spe_func(prok_database = "FAPROTAX")
meco_phase0_module_func$cal_spe_func_perc(abundance_weighted = FALSE)
meco_phase0_module_func$plot_spe_func_perc()
meco_phase0_module_func$cal_spe_func_perc(abundance_weighted = TRUE)
meco_phase0_module_func$plot_spe_func_perc()
jpeg("~/phase0_complete_function.tiff",width=3200, height=2000, res=400, bg="white")
plot.new()
meco_phase0_module_func$plot_spe_func_perc()
getOption("device")
dev.set(which = dev.next())
dev.off()

#end phase 0 function
#####

#Phase 2 - 0.1 OVERLAP
#####
head(sample_data(p2017_combined_0.1_filtered))
overlap_1_phase2 <- subset_samples(p2017_combined_0.1_filtered, phase_cat == "3") #subset phyloseq
overlap_1_phase2
overlap_1_phase2 <- prune_taxa(taxa_sums(overlap_1_phase2) > 0, overlap_1_phase2)
otu_1_phase2 <-as.data.frame(t(otu_table(overlap_1_phase2)))
tax_1_phase2 <- as.data.frame(phyloseq::tax_table(overlap_1_phase2)); colnames(tax_1_phase2) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
meta_1_phase2 <- as.data.frame(phyloseq::sample_data(overlap_1_phase2))
tree_1_phase2 = rtree(ntaxa(overlap_1_phase2), rooted=TRUE, tip.label=taxa_names(overlap_1_phase2))
merged_1_phase2 <- microtable$new(sample_table = meta_1_phase2, otu_table = otu_1_phase2, tax_table = tax_1_phase2, phylo_tree = tree_1_phase2)


otu_1 <-as.data.frame(t(otu_table(p2017_combined_0.1_filtered)))
tax_1 <- as.data.frame(phyloseq::tax_table(p2017_combined_0.1_filtered)); colnames(tax_1) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
meta_1 <- as.data.frame(phyloseq::sample_data(p2017_combined_0.1_filtered))
tree_1 = rtree(ntaxa(p2017_combined_0.1_filtered), rooted=TRUE, tip.label=taxa_names(p2017_combined_0.1_filtered))
merged_1 <- microtable$new(sample_table = meta_1, otu_table = otu_1, tax_table = tax_1, phylo_tree = tree_1)

# create object of trans_func
t2017_1 <- trans_func$new(merged_1)
# mapping the taxonomy to the database
# this can recognize prokaryotes or fungi automatically if the names of taxonomic levels are standard.
# for fungi example, see https://chiliubio.github.io/microeco_tutorial/other-dataset.html#fungi-data
# default database for prokaryotes is FAPROTAX database

t2017_1$cal_spe_func(prok_database = "FAPROTAX")
t2017_1$res_spe_func[1:5, 1:2]
# calculate the percentages for communities
# here do not consider the abundance

t2017_1$cal_spe_func_perc(abundance_weighted = FALSE)
t2017_1$res_spe_func_perc[1:5, 1:2]
t2017_1$cal_spe_func_perc(abundance_weighted = TRUE)
t2017_1$res_spe_func_perc[1:5, 1:2]

# construct a network
network_merge2017_1 <- trans_network$new(dataset = merged_1, cal_cor = "base", taxa_level = "OTU", filter_thres = 0.00001, cor_method = "spearman")
network_merge2017_1$cal_network(p_thres = 0.05, COR_cut = 0.5)
network_merge2017_1$cal_module()
# convert module info to microtable object
meco_1 <- network_merge2017_1$trans_comm(use_col = "module")
meco_1_phase2_module_func <- trans_func$new(meco_1_phase2_module)
meco_1_phase2_module_func$cal_spe_func(prok_database = "FAPROTAX")
meco_phase2_module_func$cal_spe_func_perc(abundance_weighted = FALSE)
meco_phase2_module_func$plot_spe_func_perc()
meco_phase2_module_func$cal_spe_func_perc(abundance_weighted = TRUE)
meco_phase2_module_func$plot_spe_func_perc()
jpeg("~/phase2_complete_function.tiff",width=3200, height=2000, res=400, bg="white")
plot.new()
meco_phase2_module_func$plot_spe_func_perc()
getOption("device")
dev.set(which = dev.next())
dev.off()

#end complete community phase2 function
#####

foo()
# M represents module, ordered by the nodes number from high to low
# If you want to change the group list, reset the list t2$func_group_list
t2$func_group_list
# use show_prok_func to see the detailed information of prokaryotic traits
t2$show_prok_func("fermentation")

# then we try to correlate the res_spe_func_perc of communities to environmental variables
head(meta_numeric)
head(meta_phase1)
phase1_numeric <- meta_phase1[,9:19]
phase1_numeric$specific_conductivity <- meta_phase1$SpecificConductivity_uScm
View(phase1_numeric)
t3 <- trans_env$new(dataset = merged_phase1, add_data = phase1_numeric)
t3   
t3$cal_cor(add_abund_table = t2017_phase1$res_spe_func_perc, cor_method = "spearman")
t3$plot_cor(pheatmap=TRUE)
faprotax_heatmap

t.merged_phase1 <- trans_func$new(dataset)
# https://chiliubio.github.io/microeco_tutorial/intro.html#tax4fun for the installation description
# and provide the file path of SILVA123
setwd("C:/Users/User/Documents")

getwd()
t.merged_phase1$cal_tax4fun(folderReferenceData = "~/Downloads/SILVA123")
t.merged_phase1 <- trans_func$new(merged_phase1)
class(t.merged_phase1)
t.merged_phase1$cal_tax4fun(folderReferenceData = "~/Downloads/SILVA123")
# return two files: t1$tax4fun_KO: KO file; t1$tax4fun_path: pathway file.
#view the output
t1$tax4fun_KO$Tax4FunProfile[1:5, 1:2]
# must transpose to taxa row, sample column
pathway_file <- t1$tax4fun_path$Tax4FunProfile %>% t %>% as.data.frame
# filter rownames, only keep ko+number
rownames(pathway_file) %<>% gsub("(^.*);\\s.*", "\\1", .)
# load the pathway hierarchical metadata
data(Tax4Fun2_KEGG)
View(Tax4Fun2_KEGG)
# further create a microtable object, familiar?
func1 <- microtable$new(otu_table = pathway_file, tax_table = Tax4Fun2_KEGG$ptw_desc, sample_table = t1$sample_table)
print(func1)
#trim& calculate abundance
func1$tidy_dataset()
# calculate abundance automatically at three functionaal levels: Level.1, Level.2, Level.3
func1$cal_abund()
print(func1)
#make bar plot
func2 <- trans_abund$new(func1, taxrank = "Level.1", groupmean = "Group")
func2$plot_bar(legend_text_italic = FALSE)
#test abundances to find most abundant pathways
func2 <- trans_diff$new(dataset = func1, method = "lefse", group = "Saline", alpha = 0.05, lefse_subgroup = NULL)
func2$plot_diff_bar(threshold = 3, width = 0.8)


#0.2 differential abundance by categorical phases
#####
set.seed(0451) # the seed should be set for reproducibly results.
p2017_2_filtered_transform <- transform(p2017_2_filtered, 'log10p') #(x+1) log transform
p2017_2_monitoring <- subset_samples(p2017_2_filtered_transform, well_type== "Monitoring"))
p2017_2_phase1.phase2 <- subset_samples(p2017_2_monitoring, phase== c("Phase1","Phase2")))
+)
)+
  
  
  colnames(sample_data(p2017_2_phase1.Control))
deres_p1C_WT <- diff_analysis(obj = p2017_2_phase1.phase2, classgroup = "well_type",
                              mlfun = "lda", filtermod = "pvalue", firstcomfun = "kruskal.test",
                              firstalpha = 0.05, strictmod = TRUE, secondcomfun = "wilcox.test",
                              subclmin = 3, subclwilc = TRUE, secondalpha = 0.01, lda=3.5)#modified from 3 to 3.5
deres_p1C_WT# print out results of the diff_analysis  
#boxplot
diffbox_2_deres <- ggdiffbox(obj=deres_2_WT, box_notch=FALSE, 
                             colorlist=c("darkgoldenrod4", "slateblue"), l_xlabtext="relative abundance") 
diffbox_2_deres
#plot the linear discriminate analysis effect size
WT_2_effect <- ggeffectsize(obj=deres_2_WT, 
                            lineheight=0.1, linewidth=0.2) + 
  scale_color_manual(values=c("darkgoldenrod4", "slateblue")) 
WT_2_effect
ggsave("~/2017_filtersize_color_es_boxplot.png", 
       bg = "transparent", width = 60, height = 35, units = "cm")

