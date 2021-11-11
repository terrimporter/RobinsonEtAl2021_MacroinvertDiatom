# Teresita M. Porter, Nov. 11/21
# NMDS plots with fitted vars
# Also plots of pair-wise dissimilarities for the supplement

library(stringr) # str_split
library(reshape2) # dcast
library(vegan) # rrarefy
library(ggplot2) # ggplot
library(goeveg) # scree
library(plyr) # ddply
library(gridExtra) #grid.arrange
library(otuSummary) #matrixConvert

a <- read.table("diatom_invert_combined.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE)
head(a)

# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
a.1 <- data.frame(a, do.call(rbind, str_split(a$SampleName,"_")), stringsAsFactors = FALSE)
names(a.1)[33:37]<- c("Name","Marker","Site","Replicate","Illumina sample")
# map site codes to site name
a.1$site <- a.1$Site
a.1$site <- gsub("COWL07", "Laurel7", a.1$site)
a.1$site <- gsub("COWC12", "Clair12", a.1$site)
a.1$site <- gsub("COWB18", "Beaver18", a.1$site)
a.1$site <- gsub("COWC15", "Clair15", a.1$site)

# add site status
status <- read.csv("Sites.csv", header=TRUE)
a.2 <- merge(a.1, status, by="site", all.x=TRUE)

# Create new column for dcast
a.2$SiteMarkerReplicateStatus <- paste(a.2$Site, a.2$Marker, a.2$Replicate, a.2$Condition, sep="_")

# pivot to make esv matrix
A.esv <- reshape2::dcast(a.2, SiteMarkerReplicateStatus ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# move sample to rownames then delete
rownames(A.esv) <- A.esv$SiteMarkerReplicateStatus
A.esv$SiteMarkerReplicateStatus <- NULL

#remove columns with only zeros
esv.notnull<-A.esv[,colSums(A.esv) !=0]

#remove rows with only zeros & edit rownames
esv.notnull2<-esv.notnull[rowSums(esv.notnull) !=0,]

# # check out library size variation # ok
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

# Convert to presence-absence matrix
rare.mat[rare.mat>0] <-1

# analyze diatoms and macroinvertebrates separately
diat <- rare.mat[grepl("Diatom", rownames(rare.mat)),
                 grepl("rbcL", colnames(rare.mat))]
#dim(diat)
#12 1573
mi <- rare.mat[grepl("COI", rownames(rare.mat)),
                 !grepl("rbcL", colnames(rare.mat))]
#dim(mi)
# 12 2453

# Scree plots to determine number of dimensions to use for NMDS
pdf("Scree_diat.pdf")
# check dims
dimcheckMDS(diat)
dev.off()
# use k=2
pdf("Scree_macroinvert.pdf")
# check dims
dimcheckMDS(mi)
dev.off()
# use k=2

# Do 2 dimensional NMDS
nmds2_diat <- metaMDS(diat, k=2, trymax=100)
# stress = 0.04207974
nmds2_mi <- metaMDS(mi, k=2, trymax=100)
# stress = 0.05513536

# Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
pdf("stressplot_diat.pdf")
stressplot(nmds2_diat)
gof <-goodness(nmds2_diat)
gof
plot(nmds2_diat, display = "sites", type="n", main="SSU")
points(nmds2_diat, display="sites",cex=2*gof/mean(gof))
dev.off()
# linear R2 = 0.994
pdf("stressplot_macroinvert.pdf")
stressplot(nmds2_mi)
gof <-goodness(nmds2_mi)
gof
plot(nmds2_mi, display = "sites", type="n", main="SSU")
points(nmds2_mi, display="sites",cex=2*gof/mean(gof))
dev.off()
# linear R2 = 0.983


# Create grouping matrix for samples by grabbing row names from above matrix
names_diat <- data.frame(row.names(diat), stringsAsFactors = FALSE)
# Rename the column
names(names_diat)<-"sample"
# Copy column to row names
row.names(names_diat)<-names_diat$sample
# Split first column into their own fields
names_diat.1<-data.frame(names_diat, do.call(rbind, strsplit(names_diat$sample,'_')), stringsAsFactors = FALSE)
names(names_diat.1)[2:5]<-c("Site", "Marker", "Replicate","Status")
# Remove first column
names_diat.1 <- names_diat.1[,-1]
# Grab sites/species scores from NMDS output
df <- data.frame(scores(nmds2_diat, display = "sites"))
# Put it all in one df for ggplot
gg_diat <- merge(df, names_diat.1, by="row.names")


# Create grouping matrix for samples by grabbing row names from above matrix
names_mi <- data.frame(row.names(mi), stringsAsFactors = FALSE)
# Rename the column
names(names_mi)<-"sample"
# Copy column to row names
row.names(names_mi)<-names_mi$sample
# Split first column into their own fields
names_mi.1<-data.frame(names_mi, do.call(rbind, strsplit(names_mi$sample,'_')), stringsAsFactors = FALSE)
names(names_mi.1)[2:5]<-c("Site", "Marker", "Replicate","Status")
# Remove first column
names_mi.1 <- names_mi.1[,-1]
# Grab sites/species scores from NMDS output
df <- data.frame(scores(nmds2_mi, display = "sites"))
# Put it all in one df for ggplot
gg_mi <- merge(df, names_mi.1, by="row.names")

# put it together
gg <- rbind(gg_diat, gg_mi)

# create factors
gg$Site <- factor(gg$Site, 
                   levels = c("COWB18", "COWC12", "COWC15","COWL07"),
                   labels = c("Beaver18", "Clair12", "Clair15","Laurel7"))
gg$Marker <- factor(gg$Marker,
                    levels = c("Diatom", "COI"),
                    labels = c("Diatom", "Macroinvertebrate"))
gg$Replicate <- factor(gg$Replicate,
                       levels = c("1","2","3"),
                       labels = c("1","2","3"))
gg$Status <- factor(gg$Status,
                       levels = c("Fair","Good"),
                       labels = c("Fair","Good"))

# Just the diatoms
gg_diat <- gg[gg$Marker=="Diatom",]

# Just the macroinverts
gg_mi <- gg[gg$Marker=="Macroinvertebrate",]

# Create metadata from rownames 'sample'
env_diat <- gg_diat[,c(1,4:7)]
# Assess dispersion (variance) using ANOVA
# Create distance matrix based on P-A data using binary Bray Curtis (Sorensen) dissimilarity
sor_diat <- vegdist(diat, "bray", binary=TRUE)

# Calculate beta dispersion (homogeneity needed for adonis)
# Break it down by study to ensure balanced design
bd.site<-betadisper(sor_diat, as.factor(env_diat$Site))
bd.status<-betadisper(sor_diat, as.factor(env_diat$Status))

# check for heterogeneity of beta dispersions within groups BALANCED DESIGN
set.seed(1234)
anova(bd.site) # 0.0007665 ***
anova(bd.status) # 0.0006997 ***
# PERMANOVA differences could simply be due to significant differences in beta dispersion,
# but NMDS also shows location effect (this is common) 

pdf("BetaDispersion_diatom.pdf")
par(mfrow=c(2,2))
boxplot(bd.site, main="Site", las=2)
boxplot(bd.status, main="Status")
dev.off()

# Use ADONIS to check for significant interactions between site and status
adonis(sor_diat~Status*Site, data=env_diat, permutations=999)
            # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Status     1    0.8103 0.81032  3.3119 0.21859  0.001 ***
# Site       2    0.9394 0.46971  1.9197 0.25341  0.004 ** 
# Residuals  8    1.9574 0.24467         0.52801           
# Total     11    3.7071                 1.00000   

# Create metadata from rownames 'sample'
env_mi <- gg_mi[,c(1,4:7)]
# Assess dispersion (variance) using ANOVA
# Create distance matrix based on P-A data using Bray Curtis (Sorensen) dissimilarity
sor_mi <- vegdist(mi, "bray", binary=TRUE)

# Calculate beta dispersion (homogeneity needed for adonis)
# Break it down by study to ensure balanced design
bd.site<-betadisper(sor_mi, as.factor(env_mi$Site))
bd.status<-betadisper(sor_mi, as.factor(env_mi$Status))

# check for heterogeneity of beta dispersions within groups BALANCED DESIGN
set.seed(1234)
anova(bd.site) # 0.03058 *
anova(bd.status) # 0.3285 n/s
# PERMANOVA differences could simply be due to significant differences in beta dispersion,
# but NMDS also shows location effect (this is common) 

pdf("BetaDispersion_macroinvertebrate.pdf")
par(mfrow=c(2,2))
boxplot(bd.site, main="Site", las=2)
boxplot(bd.status, main="Status")
dev.off()

# Use ADONIS to check for significant interactions between site and status
adonis(sor_mi~Status*Site, data=env_mi, permutations=999)
          # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Status     1    0.7966 0.79664  3.3021 0.19006  0.001 ***
# Site       2    1.4649 0.73245  3.0360 0.34949  0.001 ***
# Residuals  8    1.9300 0.24125         0.46045           
# Total     11    4.1916                 1.00000  

#############################################
# Fit environmental variables for diatom
# Read in metadata
D <-read.csv(file='metadata_diatom.csv',head=TRUE)
names(D)[1] <- "SampleName"

fit <- envfit(nmds2_diat, D, perm = 999)
fit.vectors <- fit[[1]]
fit.pvals <- fit[[1]]$pvals
fit.df <- as.data.frame(fit.vectors[[1]])
fit.df$pvals <- fit.pvals

fit.df.sig.diat <- fit.df[fit.df$pvals <0.05,]

fit.df.sig.diat$x <- 0
fit.df.sig.diat$y <- 0

# color by site
chulls.status.diat <- ddply(gg_diat, .(Status), function(gg_diat) gg_diat[chull(gg_diat$NMDS1, gg_diat$NMDS2), ])

# NMDS plot, add sig env vars
diatom_fitted <- ggplot(data=gg_diat, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.status.diat, aes(x=NMDS1, y=NMDS2, fill=Status), alpha=0.25) +
  geom_segment(fit.df.sig.diat, 
               mapping=aes(x=x,y=y,xend=NMDS1*0.90, yend=NMDS2*0.90), 
               arrow=arrow(angle=25, length=unit(0.10, "inches")), 
               size=0.25,
               color="black") +
  geom_text(fit.df.sig.diat,
            mapping=aes(x=NMDS1*0.95, y=NMDS2*0.95, label=rownames(fit.df.sig.diat)),
            hjust="outward", vjust="outward", size=2.5, color="black") +
  geom_point(data=gg_diat, size = 2.5, aes(color = Status, shape=Site)) +
  scale_shape_manual(values=c(21:24))+
  ggtitle("Diatoms") +
  theme_bw() +
  theme(
    text = element_text(size=8),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=8),
    legend.position = "bottom")+
  guides(color="none",
         fill="none",
         shape=guide_legend(nrow=2, byrow=TRUE, keyheight=0.1)) +
  scale_color_viridis_d(begin=0.1, end=0.5)+
  scale_fill_viridis_d(begin=0.1, end=0.5)
diatom_fitted

#############################################
# Fit environmental variables for inverts
# Read in metadata
C <-read.csv(file='metadata_invert.csv',head=TRUE)
names(C)[1] <- "SampleName"

fit <- envfit(nmds2_mi, C, perm = 999)
fit.vectors <- fit[[1]]
fit.pvals <- fit[[1]]$pvals
fit.df <- as.data.frame(fit.vectors[[1]])
fit.df$pvals <- fit.pvals

fit.df.sig.mi <- fit.df[fit.df$pvals <0.05,]

fit.df.sig.mi$x <- 0
fit.df.sig.mi$y <- 0

# color by site
chulls.status.mi <- ddply(gg_mi, .(Status), function(gg_mi) gg_mi[chull(gg_mi$NMDS1, gg_mi$NMDS2), ])

set.seed(12)
# NMDS plot, add sig env vars
invert_fitted <- ggplot(gg_mi, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.status.mi, aes(x=NMDS1, y=NMDS2, fill=Status), alpha=0.25) +
    geom_segment(fit.df.sig.mi, 
               mapping=aes(x=x,y=y,xend=NMDS1*0.90, yend=NMDS2*0.90), 
               arrow=arrow(angle=25, length=unit(0.10, "inches")), 
               size=0.25,
               color="black") +
  geom_text(fit.df.sig.mi,
            mapping=aes(x=NMDS1*0.95, y=NMDS2*0.95, label=rownames(fit.df.sig.mi)),
             size=2.5, color="black", position=position_jitter(width=ifelse(rownames(fit.df.sig.mi)=='Temp',0.15,0.01),
                                                               height=ifelse(rownames(fit.df.sig.mi)=='Temp',0.15,0.01))) +
  geom_point(data=gg_mi, size = 2.5, aes(color = Status, shape=Site)) +
  scale_shape_manual(values=c(21:24)) +
  ggtitle("Macroinvertebrates") +
  theme_bw() +
  theme(
    text = element_text(size=8),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=8),
    legend.position = "bottom"
    ) +
  guides(
         fill=FALSE,
         shape=FALSE,
         color=guide_legend(nrow=2, byrow=TRUE, keyheight=0.1)) +
  scale_color_viridis_d(begin=0.1, end=0.5)+
  scale_fill_viridis_d(begin=0.1, end=0.5)

invert_fitted

g <- grid.arrange(invert_fitted, diatom_fitted, ncol=2 )

ggsave("Fig2_NMDS.jpeg", g, width = 6, height = 4)
# only plotted fitted environmental variables with p-val <0.05





# plot range of beta dissimilarities

# turn pairwise dissimilarities to df list
sor_diat.df <- matrixConvert(sor_diat, colname = c("row", "col", "diss"))
sor_diat.df$org <- "Diatoms"
sor_diat.df$cond <- ifelse(grepl("Good", sor_diat.df$row) & 
                            grepl("Good", sor_diat.df$col), "Good", "NA")
sor_diat.df$cond <- ifelse(grepl("Fair", sor_diat.df$row) & 
                             grepl("Fair", sor_diat.df$col), "Fair", sor_diat.df$cond)
sor_diat.df2 <- sor_diat.df[!sor_diat.df$cond=="NA",]

sor_mi.df <- matrixConvert(sor_mi, colname = c("row", "col", "diss"))
sor_mi.df$org <- "Macroinvertebrates"
sor_mi.df$cond <- ifelse(grepl("Good", sor_mi.df$row) & 
                             grepl("Good", sor_mi.df$col), "Good", "NA")
sor_mi.df$cond <- ifelse(grepl("Fair", sor_mi.df$row) & 
                             grepl("Fair", sor_mi.df$col), "Fair", sor_mi.df$cond)
sor_mi.df2 <- sor_mi.df[!sor_mi.df$cond=="NA",]

# combine dfs for ggpot
sor.df <- rbind(sor_diat.df2, sor_mi.df2)

p <- ggplot(sor.df, aes(x=diss, fill=cond)) +
  geom_histogram(alpha=0.4, position = 'identity', bins=10) +
  labs(x="Binary Bray Curtis Dissimilarities",
       y = "Count") +
  scale_fill_viridis_d() +
  facet_wrap(~org) +
  theme_bw() +
  theme(
    text = element_text(size=7),
    strip.background = element_blank(),
   strip.text.x = element_text(size=7, hjust=0),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    legend.title = element_blank(),
     legend.position = "bottom") 
p

ggsave("FigS5_dissimilarities.jpeg", p, width = 6, height = 4)
# only plotted fitted environmental variables with p-val <0.05










# calc ave beta diversity with each site for aic.R
COWB18.sor <- sor.df[grepl("COWB18", sor.df$row) &
                         grepl("COWB18", sor.df$col),]
mean(COWB18.sor$diss)
# 0.816694

COWC15.sor <- sor.df[grepl("COWC15", sor.df$row) &
                       grepl("COWC15", sor.df$col),]
mean(COWC15.sor$diss)
# 0.7879363


COWC12.sor <- sor.df[grepl("COWC12", sor.df$row) &
                       grepl("COWC12", sor.df$col),]
mean(COWC12.sor$diss)
# 0.6461458


