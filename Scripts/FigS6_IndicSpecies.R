# Teresita M Porter, November 10, 2021

library(stringr) # str_split
library(reshape2) # dcast
library(vegan) # rrarefy
library(ggplot2) # ggplot
library(data.table) # setDT
library(indicspecies)

###########################
# Repeat at species level
##########################

# read in MetaWorks results
a <- read.table("diatom_invert_combined.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE)
head(a)

# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
a.1 <- data.frame(a, do.call(rbind, str_split(a$SampleName,"_")), stringsAsFactors = FALSE)
names(a.1)[33:37] <- c("Name","Marker","Site","Replicate","Illumina sample")

# Create new column for dcast
a.1$SiteMarkerReplicate <- paste(a.1$Site, a.1$Marker, a.1$Replicate, sep="_")

# create taxon field for species and add a zotu if needed
# rbcL diatom 200bp frag, 90% accuracy sBP > 0.90
# COI 200bp frag, 95% accuracy, sBP > 0.70
a.1$taxon <- ifelse(a.1$Marker=="Diatom", 
                    ifelse(a.1$sBP >= 0.90, paste(a.1$Phylum, a.1$Species, sep="_"), 
                           paste(a.1$Phylum, a.1$Genus, a.1$GlobalESV, sep="_")), 
                    "")
a.1$taxon <- ifelse(a.1$Marker=="COI", 
                    ifelse(a.1$sBP >= 0.70, paste(a.1$Phylum, a.1$Species, sep="_"), 
                           paste(a.1$Phylum, a.1$Genus, a.1$GlobalESV, sep="_")), 
                    a.1$taxon)
                    
# pivot to make sample x taxon matrix
A.species <- reshape2::dcast(a.1, SiteMarkerReplicate ~ taxon, value.var = "ESVsize", fun.aggregate = sum)

# move sample to rownames then delete
rownames(A.species) <- A.species$SiteMarkerReplicate
A.species$SiteMarkerReplicate <- NULL

#remove columns with only zeros
esv.notnull<-A.species[,colSums(A.species) !=0]

#remove rows with only zeros & edit rownames
esv.notnull2<-esv.notnull[rowSums(esv.notnull) !=0,]

# check distribution of library sizes - ok
# sums <- rowSums(esv.notnull2)
# hist(sums)
# sort(sums)

#calculate 15th percentile for rrarefy function
esv.percentile<-quantile(rowSums(esv.notnull2), prob=0.15)
esv.percentile
# 15% 
# 45937.3

# set random seed
set.seed(1234)

# Rarefy the dataset down to the 15th percentile
rare.mat <- rrarefy(esv.notnull2, sample=esv.percentile)

# process diatoms and macroinverts separately
diat <- rare.mat[grepl("Diatom", rownames(rare.mat)), grepl("rbcL", colnames(rare.mat))]
mi <- rare.mat[grepl("COI", rownames(rare.mat)), !grepl("rbcL", colnames(rare.mat))]

# Add site status 
S <- read.csv("Sites.csv", header=TRUE)
# site      lat       lon Condition
# 1  Clair15 43.46290 -80.58468      Good
# 2  Laurel7 43.47073 -80.55627      Fair
# 3 Beaver18 43.49200 -80.60978      Good
# 4  Clair12 43.46541 -80.57132      Fair


############################################
# Find diatom indicators
############################################

# set up vector to represent clusters of sites
status <- as.numeric(factor(ifelse(grepl("COWB18", rownames(diat)), "GOOD",
                 ifelse(grepl("COWC12", rownames(diat)), "FAIR",
                       ifelse(grepl("COWC15", rownames(diat)), "GOOD", "FAIR")))))
# 2==Good, 1==Fair

#Run indicator species command for status, sites in rows
diat_r_g <- multipatt(diat, status, func = "r.g", duleg = TRUE, control = how(nperm=9999))

#view results
summary(diat_r_g)

# Grab the data frame containing the p-values from each analysis
diat_sign2 <- data.frame(diat_r_g$sign)

# Grab group indicator OTUs with pvalues <= 0.05
diat_good2 <- diat_sign2[which(diat_sign2$p.value<=0.05 & 
                                diat_sign2$index==2),]
diat_fair2 <- diat_sign2[which(diat_sign2$p.value<=0.05 & 
                                diat_sign2$index==1),]

# for the significant OTUs grab A and B components
b <- data.frame(diat_r_g[[7]])
a <- data.frame(diat_r_g[[6]])

sig <- rownames(diat_good2)
a.sig <- a[rownames(a) %in% sig,]
b.sig <- b[rownames(b) %in% sig,]
diat_good2$A <- a.sig[,1]
diat_good2$B <- b.sig[,1]

sig <- rownames(diat_fair2)
a.sig <- a[rownames(a) %in% sig,]
b.sig <- b[rownames(b) %in% sig,]
diat_fair2$A <- a.sig[,2]
diat_fair2$B <- b.sig[,2]

# Add Marker to the data frame
diat_good2$Marker <- "Diatom"
diat_fair2$Marker <- "Diatom"

# Add Condition to the data frame
diat_good2$Status <- "GOOD"
diat_fair2$Status <- "FAIR"

# Move rownames to column then reset row names
diat_good2$species <- rownames(diat_good2)
diat_fair2$species <- rownames(diat_fair2)
rownames(diat_good2) <- NULL
rownames(diat_fair2) <- NULL

# Combine the data frames
indic2 <- rbind(
  diat_good2, 
  diat_fair2)
# add a placeholder for empty A and B cols
indic2$A <- ""
indic2$B <- ""
indic2$test <- "r.g"

#Run indicator species command for status, sites in rows
diat_IndVal_g <- multipatt(diat, status, func = "IndVal.g", duleg = TRUE, control = how(nperm=9999))

#view results
summary(diat_IndVal_g)

# Grab the data frame containing the p-values from each analysis
diat_sign <- data.frame(diat_IndVal_g$sign)

# Grab group indicator OTUs with pvalues <= 0.05
diat_good <- diat_sign[which(diat_sign$p.value<0.05 & 
                              diat_sign$index==2),]
diat_fair <- diat_sign[which(diat_sign$p.value<0.05 & 
                               diat_sign$index==1),]

# for the significant OTUs grab A and B components
b <- data.frame(diat_IndVal_g[[7]])
a <- data.frame(diat_IndVal_g[[6]])

sig <- rownames(diat_good)
a.sig <- a[rownames(a) %in% sig,]
b.sig <- b[rownames(b) %in% sig,]
diat_good$A <- a.sig[,2]
diat_good$B <- b.sig[,2]

sig <- rownames(diat_fair)
a.sig <- a[rownames(a) %in% sig,]
b.sig <- b[rownames(b) %in% sig,]
diat_fair$A <- a.sig[,1]
diat_fair$B <- b.sig[,1]

# Add Marker to the data frame
diat_good$Marker <- "Diatom"
diat_fair$Marker <- "Diatom"

# Add Condition to the data frame
diat_good$Status <- "GOOD"
diat_fair$Status <- "FAIR"

# Move rownames to column then reset row names
diat_good$species <- rownames(diat_good)
diat_fair$species <- rownames(diat_fair)
rownames(diat_good) <- NULL
rownames(diat_fair) <- NULL

# Combine the data frames
indic <- rbind(
  diat_good, 
  diat_fair)
indic$test <- "IndVal.g"

# put it together
diat.trad.indic <- rbind(indic, indic2)

############################################
# Repeat for macroinverts
############################################

# set up vector to represent clusters of sites
status <- as.numeric(factor(ifelse(grepl("COWB18", rownames(mi)), "GOOD",
                                   ifelse(grepl("COWC12", rownames(mi)), "FAIR",
                                          ifelse(grepl("COWC15", rownames(mi)), "GOOD", "FAIR")))))
# 2==Good, 1==Fair

#Run indicator species command for status, sites in rows
mi_r_g <- multipatt(mi, status, func = "r.g", duleg = TRUE, control = how(nperm=9999))

#view results
summary(mi_r_g)

# Grab the data frame containing the p-values from each analysis
mi_sign2 <- data.frame(mi_r_g$sign)

# Grab group indicator OTUs with pvalues <= 0.05
mi_good2 <- mi_sign2[which(mi_sign2$p.value<=0.05 & 
                                 mi_sign2$index==2),]
mi_fair2 <- mi_sign2[which(mi_sign2$p.value<=0.05 & 
                                 mi_sign2$index==1),]

# for the significant OTUs grab A and B components
b <- data.frame(mi_r_g[[7]])
a <- data.frame(mi_r_g[[6]])

sig <- rownames(mi_good2)
a.sig <- a[rownames(a) %in% sig,]
b.sig <- b[rownames(b) %in% sig,]
mi_good2$A <- a.sig[,2]
mi_good2$B <- b.sig[,2]

sig <- rownames(mi_fair2)
a.sig <- a[rownames(a) %in% sig,]
b.sig <- b[rownames(b) %in% sig,]
mi_fair2$A <- a.sig[,1]
mi_fair2$B <- b.sig[,1]

# Add Marker to the data frame
mi_good2$Marker <- "Macroinvertebrate"
mi_fair2$Marker <- "Macroinvertebrate"

# Add Condition to the data frame
mi_good2$Status <- "GOOD"
mi_fair2$Status <- "FAIR"

# Move rownames to column then reset row names
mi_good2$species <- rownames(mi_good2)
mi_fair2$species <- rownames(mi_fair2)
rownames(mi_good2) <- NULL
rownames(mi_fair2) <- NULL

# Combine the data frames
indic2 <- rbind(
  mi_good2, 
  mi_fair2)
# add a placeholder for empty A and B cols
indic2$A <- ""
indic2$B <- ""
indic2$test <- "r.g"

#Run indicator species command for status, sites in rows
mi_IndVal_g <- multipatt(mi, status, func = "IndVal.g", duleg = TRUE, control = how(nperm=9999))

#view results
summary(mi_IndVal_g)

# Grab the data frame containing the p-values from each analysis
mi_sign <- data.frame(mi_IndVal_g$sign)

# Grab group indicator OTUs with pvalues <= 0.05
mi_good <- mi_sign[which(mi_sign$p.value<0.05 & 
                               mi_sign$index==2),]
mi_fair <- mi_sign[which(mi_sign$p.value<0.05 & 
                               mi_sign$index==1),]

# for the significant OTUs grab A and B components
b <- data.frame(mi_IndVal_g[[7]])
a <- data.frame(mi_IndVal_g[[6]])

sig <- rownames(mi_good)
a.sig <- a[rownames(a) %in% sig,]
b.sig <- b[rownames(b) %in% sig,]
mi_good$A <- a.sig[,2]
mi_good$B <- b.sig[,2]

sig <- rownames(mi_fair)
a.sig <- a[rownames(a) %in% sig,]
b.sig <- b[rownames(b) %in% sig,]
mi_fair$A <- a.sig[,1]
mi_fair$B <- b.sig[,1]

# Add Marker to the data frame
mi_good$Marker <- "Macroinvertebrate"
mi_fair$Marker <- "Macroinvertebrate"

# Add Condition to the data frame
mi_good$Status <- "GOOD"
mi_fair$Status <- "FAIR"

# Move rownames to column then reset row names
mi_good$species <- rownames(mi_good)
mi_fair$species <- rownames(mi_fair)
rownames(mi_good) <- NULL
rownames(mi_fair) <- NULL

# Combine the data frames
indic <- rbind(
  mi_good, 
  mi_fair)
indic$test <- "IndVal.g"
names(indic)[1:2] <- c("s.1", "s.2")

# put it together
mi.trad.indic <- rbind(indic, indic2)




# plot diatom indicators to visualize IndVal/A/B and r stats
mi <- mi.trad.indic[,-c(1:3,5,8)]
mi$taxon <- paste(mi$Status, mi$species, sep=" ")
mi$Status <- NULL
mi$species <- NULL
mi$taxon <- gsub("_", " ", mi$taxon)

# melt for ggplot
mi.gg <- gather(mi, key, value, -c(test, taxon))
mi.gg$key <- ifelse(mi.gg$key=="stat" & mi.gg$test=="IndVal.g", "IndVal.g", mi.gg$key)
mi.gg$key <- ifelse(mi.gg$test=="r.g", "r.g", mi.gg$key)
mi.gg$test <- NULL
mi.gg2 <- mi.gg[!mi.gg$value=="",]

# create factors
mi.gg2$taxon <- factor(mi.gg2$taxon, levels=rev(levels(as.factor(mi.gg2$taxon))))
mi.gg2$key <- factor(mi.gg2$key, levels=c("A", "B", "IndVal.g", "r.g"))
mi.gg2$value <- as.numeric(mi.gg2$value)
mi.gg2$stat <- mi.gg2$value

my_breaks = c(0, 0.20, 0.4, 0.6, 0.8, 1.0)

# plot it
p1 <- ggplot(mi.gg2, aes(x=key, y=taxon)) +
  geom_tile(aes(fill=stat), stat="identity") +
  scale_fill_viridis_c(begin=0, end=1,
                       limits = c(0, 1),
                       breaks = my_breaks, labels = my_breaks
                    ) +
  coord_equal() +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'grey'),
        axis.line = element_line(colour = "black"),
        text = element_text(size=9),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")
p1



# plot diatom indicators to visualize IndVal/A/B and r stats
diat <- diat.trad.indic[,-c(1:3,5,8)]
diat$taxon <- paste(diat$Status, diat$species, sep=" ")
diat$Status <- NULL
diat$species <- NULL
diat$taxon <- gsub("_", " ", diat$taxon)

# melt for ggplot
diat.gg <- gather(diat, key, value, -c(test, taxon))
diat.gg$key <- ifelse(diat.gg$key=="stat" & diat.gg$test=="IndVal.g", "IndVal.g", diat.gg$key)
diat.gg$key <- ifelse(diat.gg$test=="r.g", "r.g", diat.gg$key)
diat.gg$test <- NULL
diat.gg2 <- diat.gg[!diat.gg$value=="",]

# create factors
diat.gg2$taxon <- factor(diat.gg2$taxon, levels=rev(levels(as.factor(diat.gg2$taxon))))
diat.gg2$key <- factor(diat.gg2$key, levels=c("A", "B", "IndVal.g", "r.g"))
diat.gg2$value <- as.numeric(diat.gg2$value)
diat.gg2$stat <- diat.gg2$value

# plot it
p2 <- ggplot(diat.gg2, aes(x=key, y=taxon)) +
  geom_tile(aes( fill=stat), stat="identity") +
  scale_fill_viridis_c(begin=0, end=1,
                       limits = c(0, 1),
                       breaks = my_breaks, labels = my_breaks
  ) +
  coord_equal() +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'grey'),
        axis.line = element_line(colour = "black"),
        text = element_text(size=9),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())
p2


g <- plot_grid(p1, p2, labels = c('A', 'B'), ncol=2, label_size = 12)

ggsave("FigS6_Indicators.jpg", g, width = 8, height = 6)




