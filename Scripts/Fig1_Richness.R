# Teresita M. Porter, Nov. 11/21
# Compare richness (ESV) and effective number of ESVs

library(stringr) # str_split
library(reshape2) # dcast
library(vegan) # rrarefy
library(ggplot2) # ggplot
library(data.table) # setDT
library(cowplot) # get_legend
library(ggpubr) # normality
library(vegetarian) # calc species equivalents

#####################################################################
# Look at richness
#####################################################################

a <- read.table("diatom_invert_combined.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE)

# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
a.1 <- data.frame(a, do.call(rbind, str_split(a$SampleName,"_")), stringsAsFactors = FALSE)
names(a.1)[33:37] <- c("Name","Marker","Site","Replicate","Illumina sample")

# Create new column for dcast
a.1$SiteMarkerReplicate <- paste(a.1$Site, a.1$Marker, a.1$Replicate, sep="_")

# calculate stats for each primer for Supplementary Table
BR5_esvs <- length(unique(a.1$GlobalESV[grepl("BR5_", a.1$GlobalESV)]))
# 575 ESVs
BR5_reads <- sum(a.1$ESVsize[grepl("BR5_", a.1$GlobalESV)])
# 131,576 reads

F230R_esvs <- length(unique(a.1$GlobalESV[grepl("F230", a.1$GlobalESV)]))
# 923 ESVs
F230R_reads <- sum(a.1$ESVsize[grepl("F230", a.1$GlobalESV)])
# 338612

mljg_esvs <- length(unique(a.1$GlobalESV[grepl("ml-jg", a.1$GlobalESV)]))
# 955 ESVs
mljg_reads <- sum(a.1$ESVsize[grepl("ml-jg", a.1$GlobalESV)])
# 259419

rbcL_esvs <- length(unique(a.1$GlobalESV[grepl("rbcL", a.1$GlobalESV)]))
# 1573 ESVs
rbcL_reads <- sum(a.1$ESVsize[grepl("rbcL", a.1$GlobalESV)])
# 574866

# filter to only retain Arthropoda and Bacillariophyta
a.2 <- a.1[a.1$Phylum=="Arthropoda" | a.1$Phylum=="Bacillariophyta",]

# pivot to make esv matrix (pool across versions, keep only substrate + sites separate)
A.esv <- reshape2::dcast(a.2, SiteMarkerReplicate ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# move sample to rownames then delete
rownames(A.esv) <- A.esv$SiteMarkerReplicate
A.esv$SiteMarkerReplicate <- NULL

#remove columns with only zeros
esv.notnull<-A.esv[,colSums(A.esv) !=0]

#remove rows with only zeros & edit rownames
esv.notnull2<-esv.notnull[rowSums(esv.notnull) !=0,]

# check distribution of library sizes - ok
#sums <- rowSums(esv.notnull2)
#hist(sums)
#sort(sums)

#calculate 15th percentile for rrarefy function
esv.percentile<-quantile(rowSums(esv.notnull2), prob=0.15)
esv.percentile
# 15% 
# 43419.05

# set random seed
set.seed(1234)

# Rarefy the dataset down to the 15th percentile
rare.mat.abund <- rrarefy(esv.notnull2, sample=esv.percentile)

# Convert to df
df <- data.frame(rare.mat.abund, stringsAsFactors = FALSE)  

# Move rownames to first column
setDT(df, keep.rownames = TRUE)[]

# Get separate method and siterep cols
setDT(df)[, paste0("S", 1:3) := tstrsplit(rn, "_")]
colnames(df)[colnames(df)=="S1"] <- "Site"
colnames(df)[colnames(df)=="S2"] <- "Taxon"
colnames(df)[colnames(df)=="S3"] <- "Replicate"

# Create a site status column
df$Condition <- df$Site
df$Condition <- gsub("COWB18","Good", df$Condition)
df$Condition <- gsub("COWC12","Fair", df$Condition)
df$Condition <- gsub("COWC15","Good", df$Condition)
df$Condition <- gsub("COWL07","Fair", df$Condition)

# now create abund and pa df's
# omit non-numeric, add back later
df2 <- data.frame(df[,c(2:3740)])
df2.abund <- df2
df2.abund <- cbind(df$rn, df2.abund, df$Site, df$Taxon, df$Replicate, df$Condition)
names(df2.abund)[1] <- "rn"
names(df2.abund)[3741:3744] <- c("Site", "Taxon", "Replicate", "Condition")

# create factors
df2.abund$rn <- factor(df2.abund$rn, levels = unique(df2.abund$rn))
df2.abund$Taxon <- factor(df2.abund$Taxon,
                       levels = c("COI","Diatom"),
                       labels = c("Macroinvertebrate","Diatom"))
df2.abund$Replicate <- factor(df2.abund$Replicate, 
                           levels = c("1", "2","3"),
                           labels = c("1", "2","3"))
df2.abund$Condition <- factor(df2.abund$Condition, levels=c("Fair", "Good"))

# use VEGETARIAN package to calculate richness and effective number of species according to Jost
# be sure to only include numeric columns

# figure out for each site_taxon set
df2.abund$Site_Taxon_Rep <- paste(df2.abund$Site, df2.abund$Taxon, df2.abund$Replicate, sep="_")
samples <- unique(df2.abund$Site_Taxon_Rep)

# Calculate richness and std error, q=0, sensitive to rare 'species'
richness.list <- list()
for (i in 1:length(samples)) {
  sample_name <- samples[i]
  richness <- d(df2.abund[df2.abund$Site_Taxon_Rep==sample_name, c(2:3740)], q=0, boot=TRUE, boot.arg = list(s.sizes = NULL, num.iter = 100))
  richness.list[[i]] <- richness
  names(richness.list)[[i]] <- sample_name
}
# turn list into df, then add to existing df
df3 <-do.call(rbind.data.frame, richness.list)

# Calculate effective number of species, q=1, exponent(x) of Shannon, less sensitive to rare 'species'
effrichness.list <- list()
for (i in 1:length(samples)) {
  sample_name <- samples[i]
  effrichness <- d(df2.abund[df2.abund$Site_Taxon_Rep==sample_name, c(2:3740)], q=1, boot=TRUE, boot.arg = list(s.sizes = NULL, num.iter = 100))
  effrichness.list[[i]] <- effrichness
  names(effrichness.list)[[i]] <- sample_name
}
# turn list into df, then add to existing df
df4 <-do.call(rbind.data.frame, effrichness.list)

# put richness and effective number species together
df5 <- data.frame(cbind(df3$D.Value, df3$StdErr, df4$D.Value, df4$StdErr))
names(df5) <- c("Richness", "RichnessStdErr", "EffectiveNumberSpecies", "EffectiveNumberSpeciesStdErr")
rownames(df5) <- rownames(df3)

# Add cols for condition and taxon
df5$Condition <- c(rep("Good",6), rep("Fair", 6), rep("Good",6), rep("Fair", 6))
df5$Taxon <- rep(c(rep("Macroinvert",3), rep("Diatom", 3)), 4)

# create factors
df5$Taxon <- factor(df5$Taxon,
                    levels = c("Macroinvert","Diatom"),
                    labels = c("Macroinvertebrate","Diatom"))
df5$Condition <- factor(df5$Condition, levels=c("Fair", "Good"))

# summarize total richness
df2.abund.fair <- df2.abund[df2.abund$Condition=="Fair",]
d(df2.abund.fair[,c(2:3740)], q=0)
# [1] 346.75
df2.abund.good <- df2.abund[df2.abund$Condition=="Good",]
d(df2.abund.good[,c(2:3740)], q=0)
# [1] 167.8333
# 346.75/167.833 = 2.06

# summarize effective number of species
d(df2.abund.fair[,c(2:3740)], q=1)
# [1] 44.9365
d(df2.abund.good[,c(2:3740)], q=1)
# [1] 21.28261
# 44.9365/21.28261 = 2.111419

# plot it
p1.tmp <- ggplot(df5, aes(x=Condition, y=Richness, color=Taxon)) +
  geom_boxplot(position=position_dodge(width = 0.8), show.legend=FALSE) +
  geom_dotplot(aes(fill=Taxon), binaxis='y', stackdir='center', position=position_dodge(0.8)) +
  ggtitle("A)") +
  labs(x="Condition", y="Richness") +
  theme(legend.title=element_blank()) +
  scale_fill_manual(values = c("#4d82ff", "#33bd2b")) +
  scale_color_manual(values = c("#4d82ff", "#33bd2b")) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(hjust = 0.5,vjust = 0),
        legend.title = element_blank(),
        legend.position = "bottom")
p1.tmp

l <- get_legend(p1.tmp)

p1 <- ggplot(df5, aes(x=Condition, y=Richness, color=Taxon)) +
  geom_boxplot(position=position_dodge(width = 0.8), show.legend=FALSE) +
  geom_dotplot(aes(fill=Taxon), binaxis='y', stackdir='center', position=position_dodge(0.8), alpha=0.5) +
  ggtitle("A)") +
  labs(x="Condition", y="Richness") +
  theme(legend.title=element_blank()) +
  scale_fill_manual(values = c("#4d82ff", "#33bd2b")) +
  scale_color_manual(values = c("#4d82ff", "#33bd2b")) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(hjust = 0.5,vjust = 0),
        legend.title = element_blank(),
        legend.position = "none")
p1

# plot it
p2 <- ggplot(df5, aes(x=Condition, y=EffectiveNumberSpecies, color=Taxon)) +
  geom_boxplot(position=position_dodge(width = 0.8), show.legend=FALSE) +
  geom_dotplot(aes(fill=Taxon), binaxis='y', stackdir='center', position=position_dodge(0.8), alpha=0.5) +
  ggtitle("B)") +
  labs(x="Condition", y="Effective Number of ESVs") +
  theme(legend.title=element_blank()) +
  scale_fill_manual(values = c("#4d82ff", "#33bd2b")) +
  scale_color_manual(values = c("#4d82ff", "#33bd2b")) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(hjust = 0.5,vjust = 0),
        legend.title = element_blank(),
        legend.position = "none")
p2

g <- plot_grid(p1, p2, l, ncol=2, rel_heights = c(1, 0.2))

ggsave("Fig1_Richness.pdf", g, height = 4, width = 8)

# test for normality pkg "ggpubr"
ggdensity(df5$Richness, 
          main = "Density plot",
          xlab = "Sample ESV richness")
# not normal

ggqqplot(df5$Richness)
# normal

#Shapiro-Wilk test of normality
shapiro.test(df5$Richness)
# data:  df2$sums
# W = 0.95026, p-value = 0.2745
# fail to reject null hypothesis of normality (normal)

# test for normality pkg "ggpubr"
ggdensity(df5$EffectiveNumberSpecies, 
          main = "Density plot",
          xlab = "Sample Effective Number of Species")
# right-skewed normal 

ggqqplot(df5$EffectiveNumberSpecies)
# normal

#Shapiro-Wilk test of normality
shapiro.test(df5$EffectiveNumberSpecies)
# data:  df2$sums
# W = 0.95159, p-value = 0.2932
# fail to reject null hypothesis of normality (normal)

# t.test to compare group means (2 conditions)
# uses Holm method for adjustment
data.frame(compare_means(Richness ~ Condition, df5, method = "t.test", paired = FALSE,
              group.by = NULL, ref.group = NULL))
# .y.         group1  group2  p         p.adj     p.format p.signif method
# 1 Richness   Fair   Good 1.121895e-05 1.1e-05  1.1e-05     **** T-test

data.frame(compare_means(EffectiveNumberSpecies ~ Condition, df5, method = "t.test", paired = FALSE,
                         group.by = NULL, ref.group = NULL))
# .y.                       group1  group2  p         p.adj p.format p.signif
# 1 EffectiveNumberSpecies   Fair   Good 0.02624493 0.026    0.026        *
#   method
# 1 T-test

# Create condition_taxon field
df5$Condition_Taxon <- paste(df5$Condition, df5$Taxon, sep="_")

# t.test to compare group means (2 conditions)
# uses Holm method for adjustment
data.frame(compare_means(Richness ~ Condition_Taxon, df5, method = "t.test", paired = FALSE,
                         group.by = NULL, ref.group = NULL))
# .y.                 group1                 group2            p  p.adj p.format p.signif method
# 1 Richness Good_Macroinvertebrate            Good_Diatom 0.7057203419 0.7100  0.70572       ns T-test
# 2 Richness Good_Macroinvertebrate Fair_Macroinvertebrate 0.0007795487 0.0039  0.00078      *** T-test
# 3 Richness Good_Macroinvertebrate            Fair_Diatom 0.0189419367 0.0570  0.01894        * T-test
# 4 Richness            Good_Diatom Fair_Macroinvertebrate 0.0001664121 0.0010  0.00017      *** T-test
# 5 Richness            Good_Diatom            Fair_Diatom 0.0059805888 0.0240  0.00598       ** T-test
# 6 Richness Fair_Macroinvertebrate            Fair_Diatom 0.0902747861 0.1800  0.09027       ns T-test

data.frame(compare_means(EffectiveNumberSpecies ~ Condition_Taxon, df5, method = "t.test", paired = FALSE,
                         group.by = NULL, ref.group = NULL))
# .y.                     group1                 group2                 p           p.adj p.format p.signif method
# 1 EffectiveNumberSpecies Good_Macroinvertebrate            Good_Diatom 0.27588734  0.92    0.276       ns T-test
# 2 EffectiveNumberSpecies Good_Macroinvertebrate Fair_Macroinvertebrate 0.22923031  0.92    0.229       ns T-test
# 3 EffectiveNumberSpecies Good_Macroinvertebrate            Fair_Diatom 0.28931908  0.92    0.289       ns T-test
# 4 EffectiveNumberSpecies            Good_Diatom Fair_Macroinvertebrate 0.06461844  0.39    0.065       ns T-test
# 5 EffectiveNumberSpecies            Good_Diatom            Fair_Diatom 0.06904210  0.39    0.069       ns T-test
# 6 EffectiveNumberSpecies Fair_Macroinvertebrate            Fair_Diatom 0.71181113  0.92    0.712       ns T-test
