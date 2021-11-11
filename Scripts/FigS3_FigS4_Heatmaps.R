# Teresita M. Porter, Apr. 2/21

library(stringr) # str_split
library(reshape2) # dcast
library(vegan) # rrarefy
library(ggplot2) # ggplot
library(data.table) # setDT
library(cowplot) # get_legend
library(ggpubr) # normality
library(dplyr) # summarize_all
library(tidyr) # gather
library(scales) # scales

#####################################################################
# Look at richness
#####################################################################

a <- read.table("diatom_invert_combined.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE)
head(a)

# replace the . in ml.jg with _
a$GlobalESV <- gsub("ml-jg", "mljg", a$GlobalESV)

# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
a.1 <- data.frame(a, do.call(rbind, str_split(a$SampleName,"_")), stringsAsFactors = FALSE)
names(a.1)[33:37] <- c("Name","Marker","Site","Replicate","Illumina sample")

# Create new column for dcast
a.1$SiteMarkerReplicate <- paste(a.1$Site, a.1$Marker, a.1$Replicate, sep="_")

# pivot to make esv matrix (pool across verions, keep only substrate + sites separate)
A.esv <- reshape2::dcast(a.1, SiteMarkerReplicate ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# move sample to rownames then delete
rownames(A.esv) <- A.esv$SiteMarkerReplicate
A.esv$SiteMarkerReplicate <- NULL

#remove columns with only zeros
esv.notnull<-A.esv[,colSums(A.esv) !=0]

#remove rows with only zeros & edit rownames
esv.notnull2<-esv.notnull[rowSums(esv.notnull) !=0,]

# check distribution of library sizes - ok
# sums <- rowSums(esv.notnull2)
# hist(sums)
# sort(sums)

#calculate 15th percentile for rrarefy function
esv.percentile <- quantile(rowSums(esv.notnull2), prob=0.15)
esv.percentile
# 15% 
# 45937.3

# set random seed
set.seed(1234)

# Rarefy the dataset down to the 15th percentile
rare.mat <- rrarefy(esv.notnull2, sample=esv.percentile)

# Convert to df
df<-data.frame(rare.mat, stringsAsFactors = FALSE)  

# Get total ESVs per sample
df$sums<-rowSums(df)

# Move rownames to first column
df2<-data.frame(df, stringsAsFactors = FALSE)
setDT(df2, keep.rownames = TRUE)[]

# Get separate method and siterep cols
setDT(df2)[, paste0("S", 1:3) := tstrsplit(rn, "_")]
colnames(df2)[colnames(df2)=="S1"] <- "Site"
colnames(df2)[colnames(df2)=="S2"] <- "Marker"
colnames(df2)[colnames(df2)=="S3"] <- "Replicate"

# Create a site status column
df2$Condition <- df2$Site
df2$Condition <- gsub("COWB18","Beaver18\nGood", df2$Condition)
df2$Condition <- gsub("COWC12","Clair12\nFair", df2$Condition)
df2$Condition <- gsub("COWC15","Clair15\nGood", df2$Condition)
df2$Condition <- gsub("COWL07","Laurel7\nFair", df2$Condition)

df2$Condition <- factor(df2$Condition)

# put rownames back
df2 <- data.frame(df2)
rownames(df2) <- df2$rn
df2$rn <- NULL

# append taxonomic info
df3 <- data.frame(t(df2))
setDT(df3, keep.rownames = TRUE)[]
names(df3)[1] <- "GlobalESV"
# remove extra non GlobalESV rows
df3 <- df3[-c(4027:4031),]

# get taxonomy for each GlobalESV
################ prefix with phylum, then add order ###############
tax <- unique(a[,c(1,15,21)])
tax$taxon <- paste(tax$Phylum, tax$Order, sep=";")
tax <- tax[,-c(2:3)]

# put it together
df4 <- merge(df3, tax, by="GlobalESV", all.x = TRUE)

# aim for 95% accuracy with COI, 200bp frag, no cutoff needed for order rank
# https://github.com/terrimporter/CO1Classifier
# aim for 95% accuracy with rbcL, 200bp frag, no cutoff needed for order rank
# https://github.com/terrimporter/rbcLdiatomClassifier

# create df for order vs site status
df5 <- gather(df4[,-1], key, value, -taxon)

# get condition
cond <- cbind.data.frame(rownames(df2), factor(df2$Replicate), factor(df2$Marker), df2$Condition)
names(cond) <- c("Site", "Replicate", "Taxon", "Condition")

# put it together 
df6 <- merge(df5, cond, by.x="key", by.y="Site", all.x=TRUE)
df7 <- df6[,-1]
df8 <- df7 %>% group_by(taxon,Replicate,Taxon,Condition) %>% summarize(n=sum(as.numeric(value)))

# add names
names(df8)[5] <- "Reads" 

# plot it
my_breaks = c(0, 10, 100, 1000, 10000)

# get diatoms
d.df <- df8[df8$Taxon=="Diatom" & df8$Reads > 0,]
d.df.count <- d.df %>% group_by(Replicate, Condition) %>% tally()
# 1 1         "Beaver18\nGood"    16
# 2 1         "Clair12\nFair"     12
# 3 1         "Clair15\nGood"      6
# 4 1         "Laurel7\nFair"     17
# 5 2         "Beaver18\nGood"     9
# 6 2         "Clair12\nFair"     12
# 7 2         "Clair15\nGood"     10
# 8 2         "Laurel7\nFair"     13
# 9 3         "Beaver18\nGood"     6
# 10 3         "Clair12\nFair"     13
# 11 3         "Clair15\nGood"     11
# 12 3         "Laurel7\nFair"     14
d.df.count$map <- paste(d.df.count$Replicate, d.df.count$Condition, sep="\n")
d.df.count$labels <- paste(d.df.count$map, d.df.count$n, sep="\n")

df8$map <- paste(df8$Replicate, df8$Condition, sep="\n")
df9 <- merge(df8, d.df.count, by="map", all.x = TRUE)

# create factors, reverse taxon on y axis
df9$taxon <- factor(df9$taxon, levels=rev(levels(as.factor(df9$taxon))))
df9$Samples <- factor(df9$labels, levels=c("1\nClair12\nFair\n12", "2\nClair12\nFair\n12", "3\nClair12\nFair\n13",
                                              "1\nLaurel7\nFair\n17", "2\nLaurel7\nFair\n13", "3\nLaurel7\nFair\n14",
                                              "1\nBeaver18\nGood\n16", "2\nBeaver18\nGood\n9", "3\nBeaver18\nGood\n6",
                                              "1\nClair15\nGood\n6", "2\nClair15\nGood\n10", "3\nClair15\nGood\n11"))

tmp <- ggplot(df9[df9$Taxon=="Diatom" & df9$Reads > 0,], aes(x=Samples, y=taxon)) +
  geom_tile(aes(fill=Reads), stat="identity", show.legend = TRUE) +
  ggtitle("A) Diatoms") +
  scale_fill_viridis_c(begin=0,
                       end=1, 
                       breaks = my_breaks, labels = my_breaks,
                       trans=scales::pseudo_log_trans(sigma = 0.001)) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'grey'),
        axis.line = element_line(colour = "black"),
        text = element_text(size=12),
        axis.text.x = element_text(size=9),
        legend.text = element_text(size=8),
        legend.position = "bottom",
        axis.title.y = element_blank())
tmp
l <- get_legend(tmp)

d <- ggplot(df9[df9$Taxon=="Diatom" & df9$Reads > 0,], aes(x=Samples, y=taxon)) +
  geom_tile(aes(fill=Reads), stat="identity", show.legend = FALSE) +
  scale_fill_viridis_c(begin=0,
                       end=1, 
                       breaks = my_breaks, labels = my_breaks,
                       trans=scales::pseudo_log_trans(sigma = 0.001)) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'grey'),
        axis.line = element_line(colour = "black"),
        text = element_text(size=12),
        axis.text.x = element_text(size=8),
        axis.title.y = element_blank())
d

g <- plot_grid(d, l, ncol=1, rel_heights = c(7,1))
ggsave("FigS3_Diatom_Heatmap.jpg", g, width = 8, height = 10)

# get macroinverts
mi.df <- df8[df8$Taxon=="COI" & df8$Reads > 0,]
mi.df.count <- mi.df %>% group_by(Replicate, Condition) %>% tally()
# 1 1         "Beaver18\nGood"    29
# 2 1         "Clair12\nFair"     30
# 3 1         "Clair15\nGood"     17
# 4 1         "Laurel7\nFair"     42
# 5 2         "Beaver18\nGood"    17
# 6 2         "Clair12\nFair"     32
# 7 2         "Clair15\nGood"     17
# 8 2         "Laurel7\nFair"     43
# 9 3         "Beaver18\nGood"    12
# 10 3         "Clair12\nFair"     32
# 11 3         "Clair15\nGood"     26
# 12 3         "Laurel7\nFair"     51
mi.df.count$map <- paste(mi.df.count$Replicate, mi.df.count$Condition, sep="\n")
mi.df.count$labels <- paste(mi.df.count$map, mi.df.count$n, sep="\n")

df8$map <- paste(df8$Replicate, df8$Condition, sep="\n")
df9 <- merge(df8, mi.df.count, by="map", all.x = TRUE)

# create factors, reverse taxon on y axis
df9$taxon <- factor(df9$taxon, levels=rev(levels(as.factor(df9$taxon))))
df9$Samples <- factor(df9$labels, levels=c("1\nClair12\nFair\n30", "2\nClair12\nFair\n32", "3\nClair12\nFair\n32",
                                           "1\nLaurel7\nFair\n42", "2\nLaurel7\nFair\n43", "3\nLaurel7\nFair\n51",
                                           "1\nBeaver18\nGood\n29", "2\nBeaver18\nGood\n17", "3\nBeaver18\nGood\n12",
                                           "1\nClair15\nGood\n17", "2\nClair15\nGood\n17", "3\nClair15\nGood\n26"))



# plot macroinverts 
mi <- ggplot(df9[df9$Taxon=="COI" & df9$Reads > 0,], aes(x=Samples, y=taxon)) +
  geom_tile(aes(fill=Reads), stat="identity", show.legend = FALSE) +
  scale_fill_viridis_c(begin=0,
                       end=1, 
                       breaks = my_breaks, labels = my_breaks,
                       trans=scales::pseudo_log_trans(sigma = 0.001)) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'grey'),
        axis.line = element_line(colour = "black"),
        text = element_text(size=12),
        axis.text.x = element_text(size=8),
        axis.title.y = element_blank())
mi

g <- plot_grid(mi, l, ncol=1, rel_heights = c(7,1))
ggsave("FigS4_Macroinvertebrate_Heatmap.jpg", g, width = 8, height = 10)

