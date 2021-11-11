# Teresita M. Porter, Nov. 11/21
# Do cheddar trophic analysis (food webs)
# combine this with network analysis based on resource-consumer relationships

library(cheddar)
library(stringr) # str_split
library(ggplot2) # ggplot
library(data.table) # setDT
library(plyr) # ddply
library(ggpubr) #text_grob
library(igraph) # for more plotting options in a circular format
library(tidyr) # gather
library(cowplot) # get_legend

#Now load data by individual sites (B18,C12,C15,L07)
B18 <- LoadCommunity('waterloo/communities/Beaver18')
B18

C12 <- LoadCommunity('waterloo/communities/Clair12')
C12

C15 <- LoadCommunity('waterloo/communities/Clair15')
C15

L7 <- LoadCommunity('waterloo/communities/Laurel7')
L7

# Number of nodes
NumberOfNodes(B18) #83
NumberOfNodes(C12) #108
NumberOfNodes(C15) #45
NumberOfNodes(L7) #105

# Number of Trophic Links
NumberOfTrophicLinks(B18)
# 316
NumberOfTrophicLinks(C12)
# 402
NumberOfTrophicLinks(C15)
# 106
NumberOfTrophicLinks(L7)
# 390

# Density
LinkageDensity(B18)
#3.8
LinkageDensity(C12)
#3.7
LinkageDensity(C15)
#2.4
LinkageDensity(L7)
#3.7

# Connectance
DirectedConnectance(B18)
#0.05
DirectedConnectance(C12)
#0.03
DirectedConnectance(C15)
#0.05
DirectedConnectance(L7)
#0.04

# path lengths (small-world)
CharacteristicPathLength(B18) #2.4
CharacteristicPathLength(C12) #2.3
CharacteristicPathLength(C15) #2.2
CharacteristicPathLength(L7) #3.0

# max Trophic height
max(TrophicHeight(B18)) # 3.7
max(TrophicHeight(C12)) # 4.8
max(TrophicHeight(C15)) # 2.9
max(TrophicHeight(L7)) # 4.4

# identify top predators
sort(TrophicHeight(B18))
sort(TrophicHeight(C12))
sort(TrophicHeight(C15))
sort(TrophicHeight(L7))

# Plot cheddar food webs
pdf("Fig3_foodwebs.pdf", width = 8, height = 10)
par(mfrow=c(3,2))

p1 <- PlotWebByLevel(B18, 
                     main="Beaver18 - Good",
                     cex.main = 1.5,
                     font.main = 1,
                     xlab="",
                     level="ChainAveragedTrophicLevel", 
                     show.level.labels=TRUE,
                     ylim = c(1,5),
                     weight.by=NULL,
                     cex=1)

p2 <- PlotWebByLevel(C12, 
                     main="Clair12 - Fair",
                     cex.main = 1.5,
                     font.main = 1,
                     xlab="",
                     level="ChainAveragedTrophicLevel", 
                     show.level.labels=TRUE,
                     ylim = c(1,5),
                     weight.by=NULL,
                     cex=1)

p3 <- PlotWebByLevel(C15, 
                     main="Clair15 - Good",
                     cex.main = 1.5,
                     font.main = 1,
                     xlab="",
                     level="ChainAveragedTrophicLevel", 
                     show.level.labels=TRUE,
                     ylim = c(1,5),
                     weight.by=NULL,
                     cex=1)

p4 <- PlotWebByLevel(L7, 
                     main="Laurel7 - Fair",
                     cex.main = 1.5,
                     font.main = 1,
                     xlab="",
                     level="ChainAveragedTrophicLevel", 
                     show.level.labels=TRUE,
                     ylim = c(1,5),
                     weight.by=NULL,
                     cex=1)

plot.new()
plot.new()

# get cheddar style defaults
DefaultCategoryColours()
DefaultCategoryLabelColours()
DefaultCategorySymbols()
DefaultLinkColour()

op <- par(cex = 0.8)
legend("top", 
       legend=c("Producers (Diatoms)", "Macroinvertebrates", "Macroinvertebrates (Cannibals)"), 
       pch = c(21,22,21), 
       col = c("#33bd2b", "#4d82ff", "#26417fff"),
       pt.bg = c("#33bd2b", "#4d82ff", "#a6c0ffff"),
      y.intersp=1,
      bty = "n",
      cex=1)

dev.off()









# Circularize networks to calculate addional network properties
pdf("Fig4_circlegraphs.pdf", width = 8, height = 8)
par(mfrow=c(2,2))

# Calculate stats and create circularized network for Beaver18
# turn df into graph get hubscores/groups etc to find keystone taxa
B18.df <- TLPS(B18)

# turn TLPS(B18) into a directed graph (resource -> consumer)
B18.g <- graph.data.frame(d = B18.df, directed = TRUE)

# add category from NPS(B18) to graph
B18.g = set_vertex_attr(B18.g, "category", index=V(B18.g), NPS(B18)$category)
# set color using category
V(B18.g)$color <- ifelse(V(B18.g)$category == "producer", "#33bd2b", "#4d82ff")
V(B18.g)$frame.color <- ifelse(V(B18.g)$category == "producer", "#33bd2b", "#4d82ff")
V(B18.g)$vertex.label.color <- ifelse(degree(B18.g)>3, "black", "transparent")
V(B18.g)$label <- as.character(seq(1, length(V(B18.g))))

# find clusters
wtc <- cluster_walktrap(B18.g)
# circular layout
coords <- layout_in_circle(B18.g, order=order(membership(wtc)))

#get groups for plotting
grps<-groups(wtc)

# set up groups for OTUs part of clusters
grps <- lapply(1:length(groups(wtc)), function(x) which(V(B18.g)$name %in% grps[[x]]))

# calculate modularity # add to table above
mod <- modularity(wtc)
# 0.07987023

# don't bother to size vertices by centrality or labels to reduce visual clutter
# plot circle graph
plot(B18.g, layout=coords, 
     vertex.size=5,
     edge.arrow.size=0, 
     vertex.label.color = "transparent",
     vertex.label.dist = 0, vertex.label.family="Helvetica", vertex.label.font=1,
     vertex.label.cex = 0.75,
     margin=c(0.125, 0.125, 0.125, 0.125),
     mark.groups = grps, mark.col = "transparent", mark.border ="black")
mtext("Beaver18 - Good", side=3, cex=0.8, font=2)

# track groups
l <- list()
for (i in 1:length(grps)) {
  species <- grps[[i]]
  l[[i]] <- data.frame(Group=rep(i, length(species)), Taxon=species)
}
df.groups <- do.call(rbind.data.frame, l)

# track taxa
df.V <- data.frame(VertexNumber=V(B18.g)$label, Taxon=V(B18.g)$name)
df.groups.V <- merge(df.groups, df.V, by.x="Taxon", by.y="VertexNumber", all.x=TRUE)
names(df.groups.V) <- c("VertexNumber", "Group", "Taxon")
df.groups.V <- df.groups.V[,c("Group", "VertexNumber", "Taxon")]

# add degree
df.degree <- data.frame(degree(B18.g))
setDT(df.degree, keep.rownames = TRUE)[]
names(df.degree)[1:2] <- c("Taxon", "Degree")
df.degree <- data.frame(df.degree)
df.groups.V.degree <- merge(df.groups.V, df.degree, by="Taxon", all.x=TRUE)

# add hubscores
df.hubscore <- data.frame(hub_score(B18.g)$vector)
setDT(df.hubscore, keep.rownames = TRUE)[]
names(df.hubscore)[1:2] <- c("Taxon", "HubScore")
df.hubscore <- data.frame(df.hubscore)
df.hubscore$HubScore <- lapply(df.hubscore$HubScore, format, scientific=FALSE)
df.hubscore$HubScore <- round(as.numeric(df.hubscore$HubScore), 2)
df.groups.V.degree.hubscore <- merge(df.groups.V.degree, df.hubscore, by="Taxon", all.x=TRUE)
df.groups.V.degree.hubscore <- df.groups.V.degree.hubscore[order(df.groups.V.degree.hubscore$Group,
                                                                 df.groups.V.degree.hubscore$Taxon),]

x<-order(membership(wtc))
df.groups.V.degree.hubscore <- df.groups.V.degree.hubscore %>%
  arrange(sapply(VertexNumber, function(y) which(y == x)))

df.final <- merge(df.groups.V.degree.hubscore, B18$nodes, by.x="Taxon", by.y="node", all.x=TRUE)
df.final <- df.final[,c(6,1:5)]
names(df.final)[1] <- "Category"
df.final$Category <- gsub("producer", "Diatom", df.final$Category)
df.final$Category <- gsub("invertebrate", "Macroinvertebrate", df.final$Category)
write.csv(df.final, "Beaver18centrality.csv", quote=FALSE, row.names = FALSE)

df.final2 <- gather(df.final, key, value, -c(Category, Taxon, Group, VertexNumber))
df.final2$Category <- factor(df.final2$Category, levels=c("Macroinvertebrate","Diatom"))

# plot centrlaity scores
b18.g.tmp <- ggplot(df.final2, aes(x=Category, y=value, color=Category))+
  geom_boxplot(show.legend=FALSE) +
  geom_jitter(aes(fill=Category), size=2, alpha=1, width = 0.2, shape=16) +
  ggtitle("A)") +
  labs(x="", y="Centrlaity") +
  scale_fill_manual(values = c("#4d82ff", "#33bd2b")) +
  scale_color_manual(values = c("#4d82ff", "#33bd2b")) +
  facet_wrap(~key, scales = "free_y")+
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        legend.title = element_blank())
l <- get_legend(b18.g.tmp)

b18.g <- ggplot(df.final2, aes(x=Category, y=value, color=Category))+
  geom_boxplot(show.legend=FALSE) +
  geom_jitter(aes(fill=Category), size=1.5, alpha=0.5, width = 0.2, shape=16) +
  ggtitle("A)") +
  labs(x="", y="Centrlaity") +
  scale_fill_manual(values = c("#4d82ff", "#33bd2b")) +
  scale_color_manual(values = c("#4d82ff", "#33bd2b")) +
  facet_wrap(~key, scales = "free_y")+
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        legend.position = "none",
        legend.title = element_blank())







# Calculate stats and create circularized network for Clair12
# turn df into graph get hubscores/groups etc to find keystone taxa
C12.df <- TLPS(C12)

# turn TLPS(C12) into a directed graph (resource -> consumer)
C12.g <- graph.data.frame(d = C12.df, directed = TRUE)

# add category from NPS(C12) to graph
C12.g = set_vertex_attr(C12.g, "category", index=V(C12.g), NPS(C12)$category)
# set color using category
V(C12.g)$color <- ifelse(V(C12.g)$category == "producer", "#33bd2b", "#4d82ff")
V(C12.g)$frame.color <- ifelse(V(C12.g)$category == "producer", "#33bd2b", "#4d82ff")
V(C12.g)$vertex.label.color <- ifelse(degree(C12.g)>3, "black", "transparent")
V(C12.g)$label <- as.character(seq(1, length(V(C12.g))))

# find clusters
wtc <- cluster_walktrap(C12.g)
# circular layout
coords <- layout_in_circle(C12.g, order=order(membership(wtc)))

#get groups for plotting
grps<-groups(wtc)

# set up groups for OTUs part of clusters
grps <- lapply(1:length(groups(wtc)), function(x) which(V(C12.g)$name %in% grps[[x]]))

# calculate modularity
mod <- modularity(wtc)
# 0.06343061

# plot circle graph
plot(C12.g, layout=coords, 
     vertex.size=5,
     edge.arrow.size=0, 
     vertex.label.color = "transparent",
     vertex.label.dist = 0, vertex.label.family="Helvetica", vertex.label.font=1,
     vertex.label.cex = 0.75,
     margin=c(0.125, 0.125, 0.125, 0.125),
     mark.groups = grps, mark.col = "transparent", mark.border ="black")
mtext("Claire12 - Fair", side=3, cex=0.8, font=2)

# track groups
l <- list()
for (i in 1:length(grps)) {
  species <- grps[[i]]
  l[[i]] <- data.frame(Group=rep(i, length(species)), Taxon=species)
}
df.groups <- do.call(rbind.data.frame, l)

# track taxa
df.V <- data.frame(VertexNumber=V(C12.g)$label, Taxon=V(C12.g)$name)
df.groups.V <- merge(df.groups, df.V, by.x="Taxon", by.y="VertexNumber", all.x=TRUE)
names(df.groups.V) <- c("VertexNumber", "Group", "Taxon")
df.groups.V <- df.groups.V[,c("Group", "VertexNumber", "Taxon")]

# add degree
df.degree <- data.frame(degree(C12.g))
setDT(df.degree, keep.rownames = TRUE)[]
names(df.degree)[1:2] <- c("Taxon", "Degree")
df.degree <- data.frame(df.degree)
df.groups.V.degree <- merge(df.groups.V, df.degree, by="Taxon", all.x=TRUE)

# add hubscores
df.hubscore <- data.frame(hub_score(C12.g)$vector)
setDT(df.hubscore, keep.rownames = TRUE)[]
names(df.hubscore)[1:2] <- c("Taxon", "HubScore")
df.hubscore <- data.frame(df.hubscore)
df.hubscore$HubScore <- lapply(df.hubscore$HubScore, format, scientific=FALSE)
df.hubscore$HubScore <- round(as.numeric(df.hubscore$HubScore), 2)
df.groups.V.degree.hubscore <- merge(df.groups.V.degree, df.hubscore, by="Taxon", all.x=TRUE)
df.groups.V.degree.hubscore <- df.groups.V.degree.hubscore[order(df.groups.V.degree.hubscore$Group,
                                                                 df.groups.V.degree.hubscore$Taxon),]

x<-order(membership(wtc))
df.groups.V.degree.hubscore <- df.groups.V.degree.hubscore %>%
  arrange(sapply(VertexNumber, function(y) which(y == x)))

df.final <- merge(df.groups.V.degree.hubscore, C12$nodes, by.x="Taxon", by.y="node", all.x=TRUE)
df.final <- df.final[,c(6,1:5)]
names(df.final)[1] <- "Category"
df.final$Category <- gsub("producer", "Diatom", df.final$Category)
df.final$Category <- gsub("invertebrate", "Macroinvertebrate", df.final$Category)
write.csv(df.final, "Clair12centrality.csv", quote=FALSE, row.names = FALSE)

df.final2 <- gather(df.final, key, value, -c(Category, Taxon, Group, VertexNumber))
df.final2$Category <- factor(df.final2$Category, levels=c("Macroinvertebrate", "Diatom"))

# plot centrality scores
c12.g <- ggplot(df.final2, aes(x=Category, y=value, color=Category))+
  geom_boxplot(show.legend=FALSE) +
  geom_jitter(aes(fill=Category), size=1.5, alpha=0.5, width = 0.2, shape=16) +
  ggtitle("C)") +
  labs(x="", y="Centrality") +
  scale_fill_manual(values = c("#4d82ff", "#33bd2b")) +
  scale_color_manual(values = c("#4d82ff", "#33bd2b")) +
  facet_wrap(~key, scales = "free_y")+
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        legend.position = "none",
        legend.title = element_blank())






# Calculate stats and create circularized network for Clair15
# turn df into graph get hubscores/groups etc to find keystone taxa
C15.df <- TLPS(C15)

# turn TLPS(C15) into a directed graph (resource -> consumer)
C15.g <- graph.data.frame(d = C15.df, directed = TRUE)

# add category from NPS(C15) to graph
C15.g = set_vertex_attr(C15.g, "category", index=V(C15.g), NPS(C15)$category)
# set color using category
V(C15.g)$color <- ifelse(V(C15.g)$category == "producer", "#33bd2b", "#4d82ff")
V(C15.g)$frame.color <- ifelse(V(C15.g)$category == "producer", "#33bd2b", "#4d82ff")
V(C15.g)$vertex.label.color <- ifelse(degree(C15.g)>3, "black", "transparent")
V(C15.g)$label <- as.character(seq(1, length(V(C15.g))))

# find clusters
wtc <- cluster_walktrap(C15.g)
# circular layout
coords <- layout_in_circle(C15.g, order=order(membership(wtc)))

#get groups for plotting
grps<-groups(wtc)

# set up groups for OTUs part of clusters
grps <- lapply(1:length(groups(wtc)), function(x) which(V(C15.g)$name %in% grps[[x]]))

# calculate modularity
mod <- modularity(wtc)
# 0.1638

# plot circle graph
plot(C15.g, layout=coords, 
     vertex.size=5,
     edge.arrow.size=0, 
     vertex.label.color = "transparent",
     vertex.label.dist = 0, vertex.label.family="Helvetica", vertex.label.font=1,
     vertex.label.cex = 0.75,
     margin=c(0.125, 0.125, 0.125, 0.125),
     mark.groups = grps, mark.col = "transparent", mark.border ="black")
mtext("Claire15 - Good", side=3, cex=0.8, font=2)

# track groups
l <- list()
for (i in 1:length(grps)) {
  species <- grps[[i]]
  l[[i]] <- data.frame(Group=rep(i, length(species)), Taxon=species)
}
df.groups <- do.call(rbind.data.frame, l)

# track taxa
df.V <- data.frame(VertexNumber=V(C15.g)$label, Taxon=V(C15.g)$name)
df.groups.V <- merge(df.groups, df.V, by.x="Taxon", by.y="VertexNumber", all.x=TRUE)
names(df.groups.V) <- c("VertexNumber", "Group", "Taxon")
df.groups.V <- df.groups.V[,c("Group", "VertexNumber", "Taxon")]

# add degree
df.degree <- data.frame(degree(C15.g))
setDT(df.degree, keep.rownames = TRUE)[]
names(df.degree)[1:2] <- c("Taxon", "Degree")
df.degree <- data.frame(df.degree)
df.groups.V.degree <- merge(df.groups.V, df.degree, by="Taxon", all.x=TRUE)

# add hubscores
df.hubscore <- data.frame(hub_score(C15.g)$vector)
setDT(df.hubscore, keep.rownames = TRUE)[]
names(df.hubscore)[1:2] <- c("Taxon", "HubScore")
df.hubscore <- data.frame(df.hubscore)
df.hubscore$HubScore <- lapply(df.hubscore$HubScore, format, scientific=FALSE)
df.hubscore$HubScore <- round(as.numeric(df.hubscore$HubScore), 2)
df.groups.V.degree.hubscore <- merge(df.groups.V.degree, df.hubscore, by="Taxon", all.x=TRUE)
df.groups.V.degree.hubscore <- df.groups.V.degree.hubscore[order(df.groups.V.degree.hubscore$Group,
                                                                 df.groups.V.degree.hubscore$Taxon),]

x<-order(membership(wtc))
df.groups.V.degree.hubscore <- df.groups.V.degree.hubscore %>%
  arrange(sapply(VertexNumber, function(y) which(y == x)))

df.final <- merge(df.groups.V.degree.hubscore, C15$nodes, by.x="Taxon", by.y="node", all.x=TRUE)
df.final <- df.final[,c(6,1:5)]
names(df.final)[1] <- "Category"
df.final$Category <- gsub("producer", "Diatom", df.final$Category)
df.final$Category <- gsub("invertebrate", "Macroinvertebrate", df.final$Category)
write.csv(df.final, "Clair15centrality.csv", quote=FALSE, row.names = FALSE)

df.final2 <- gather(df.final, key, value, -c(Category, Taxon, Group, VertexNumber))
df.final2$Category <- factor(df.final2$Category, levels=c("Macroinvertebrate","Diatom"))

# plot centrality
c15.g <- ggplot(df.final2, aes(x=Category, y=value, color=Category))+
  geom_boxplot() +
  geom_jitter(aes(fill=Category), size=1.5, alpha=0.5, width = 0.2, shape=16) +
  ggtitle("B)") +
  labs(x="", y="Centrality") +
  facet_wrap(~key, scales = "free_y") +
  scale_fill_manual(values = c("#4d82ff", "#33bd2b")) +
  scale_color_manual(values = c("#4d82ff", "#33bd2b")) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        legend.position = "none",
        legend.title = element_blank())





# Calculate stats and create circularized network for Laurel7
# turn df into graph get hubscores/groups etc to find keystone taxa
L7.df <- TLPS(L7)

# turn TLPS(L7) into a directed graph (resource -> consumer)
L7.g <- graph.data.frame(d = L7.df, directed = TRUE)

# add category from NPS(L7) to graph
L7.g = set_vertex_attr(L7.g, "category", index=V(L7.g), NPS(L7)$category)
# set color using category
V(L7.g)$color <- ifelse(V(L7.g)$category == "producer", "#33bd2b", "#4d82ff")
V(L7.g)$frame.color <- ifelse(V(L7.g)$category == "producer", "#33bd2b", "#4d82ff")
V(L7.g)$vertex.label.color <- ifelse(degree(L7.g)>3, "black", "transparent")
V(L7.g)$label <- as.character(seq(1, length(V(L7.g))))

# find clusters
wtc <- cluster_walktrap(L7.g)
# circular layout
coords <- layout_in_circle(L7.g, order=order(membership(wtc)))

#get groups for plotting
grps<-groups(wtc)

# set up groups for OTUs part of clusters
grps <- lapply(1:length(groups(wtc)), function(x) which(V(L7.g)$name %in% grps[[x]]))

# calculate modularity
mod <- modularity(wtc)
# 0.06482168

# plot circle graph
plot(L7.g, layout=coords, 
     vertex.size=5,
     edge.arrow.size=0, 
     vertex.label.color = "transparent",
     vertex.label.dist = 0, vertex.label.family="Helvetica", vertex.label.font=1,
     vertex.label.cex = 0.75,
     margin=c(0.125, 0.125, 0.125, 0.125),
     mark.groups = grps, mark.col = "transparent", mark.border ="black")
mtext("Laurel7 - Fair", side=3, cex=0.8, font=2)

# track groups
l <- list()
for (i in 1:length(grps)) {
  species <- grps[[i]]
  l[[i]] <- data.frame(Group=rep(i, length(species)), Taxon=species)
}
df.groups <- do.call(rbind.data.frame, l)

# track taxa
df.V <- data.frame(VertexNumber=V(L7.g)$label, Taxon=V(L7.g)$name)
df.groups.V <- merge(df.groups, df.V, by.x="Taxon", by.y="VertexNumber", all.x=TRUE)
names(df.groups.V) <- c("VertexNumber", "Group", "Taxon")
df.groups.V <- df.groups.V[,c("Group", "VertexNumber", "Taxon")]

# add degree
df.degree <- data.frame(degree(L7.g))
setDT(df.degree, keep.rownames = TRUE)[]
names(df.degree)[1:2] <- c("Taxon", "Degree")
df.degree <- data.frame(df.degree)
df.groups.V.degree <- merge(df.groups.V, df.degree, by="Taxon", all.x=TRUE)

# add hubscores
df.hubscore <- data.frame(hub_score(L7.g)$vector)
setDT(df.hubscore, keep.rownames = TRUE)[]
names(df.hubscore)[1:2] <- c("Taxon", "HubScore")
df.hubscore <- data.frame(df.hubscore)
df.hubscore$HubScore <- lapply(df.hubscore$HubScore, format, scientific=FALSE)
df.hubscore$HubScore <- round(as.numeric(df.hubscore$HubScore), 2)
df.groups.V.degree.hubscore <- merge(df.groups.V.degree, df.hubscore, by="Taxon", all.x=TRUE)
df.groups.V.degree.hubscore <- df.groups.V.degree.hubscore[order(df.groups.V.degree.hubscore$Group,
                                                                 df.groups.V.degree.hubscore$Taxon),]

x<-order(membership(wtc))
df.groups.V.degree.hubscore <- df.groups.V.degree.hubscore %>%
  arrange(sapply(VertexNumber, function(y) which(y == x)))

df.final <- merge(df.groups.V.degree.hubscore, L7$nodes, by.x="Taxon", by.y="node", all.x=TRUE)
df.final <- df.final[,c(6,1:5)]
names(df.final)[1] <- "Category"
df.final$Category <- gsub("producer", "Diatom", df.final$Category)
df.final$Category <- gsub("invertebrate", "Macroinvertebrate", df.final$Category)
write.csv(df.final, "Laurel7centrality.csv", quote=FALSE, row.names = FALSE)

df.final2 <- gather(df.final, key, value, -c(Category, Taxon, Group, VertexNumber))
df.final2$Category <- factor(df.final2$Category, levels=c("Macroinvertebrate","Diatom"))

# plot centrality scores
l7.g <- ggplot(df.final2, aes(x=Category, y=value, color=Category))+
  geom_boxplot(show.legend=FALSE) +
  geom_jitter(aes(fill=Category), size=1.5, alpha=0.5, width = 0.2, shape=16) +
  ggtitle("D)") +
  labs(x="", y="Centrality") +
  facet_wrap(~key, scales = "free_y") +
  scale_fill_manual(values = c("#4d82ff", "#33bd2b")) +
  scale_color_manual(values = c("#4d82ff", "#33bd2b")) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        legend.position = "none",
        legend.title = element_blank())

dev.off()

# plot centrality figures in one plate
g <- plot_grid(b18.g, c15.g, c12.g, l7.g, l, nrow=3, rel_heights=c(1,1,0.2))
ggsave("Fig5_centrality.jpg", g, width = 8, height = 8)


