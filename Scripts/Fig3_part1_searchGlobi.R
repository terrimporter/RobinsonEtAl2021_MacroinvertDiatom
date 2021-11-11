# Teresita M. Porter, Nov. 11/21
# Script to query GloBI

# install.packages("rglobi")
library(rglobi)
library(stringr) # str_split, word

# read in MetaWorks results
a <- read.table("diatom_invert_combined.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE)
head(a)

# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
a.1 <- data.frame(a, do.call(rbind, str_split(a$SampleName,"_")), stringsAsFactors = FALSE)
names(a.1)[33:37] <- c("Name","Marker","Site","Replicate","Illumina sample")

# Create new column for dcast
a.1$SiteMarkerReplicate <- paste(a.1$Site, a.1$Marker, a.1$Replicate, sep="_")

# Figure out good species / all species
# create taxon field for species and add a zotu if needed
# rbcL diatom 200bp frag, 90% accuracy sBP > 0.90
# COI 200bp frag, 95% accuracy, sBP > 0.70
a.1$goodspecies <- ifelse(a.1$Marker=="Diatom", 
                    ifelse(a.1$sBP >= 0.90, a.1$Species, ""), 
                    "")
a.1$goodspecies <- ifelse(a.1$Marker=="COI", 
                    ifelse(a.1$sBP >= 0.70, a.1$Species, ""), 
                    a.1$goodspecies)
length(unique(a.1$goodspecies))
# 289 / 1,186 = 0.24, so 76% of ESVs can't be confidently assigned at species rank

# check genus rank instead
a.1$goodgenera <- ifelse(a.1$Marker=="Diatom",
                         ifelse(a.1$gBP >= 0, a.1$Genus, ""),
                         "")
a.1$goodgenera <- ifelse(a.1$Marker=="COI",
                         ifelse(a.1$gBP >= 0, a.1$Genus, ""),
                         a.1$goodgenera)
length(unique(a.1$goodgenera))
# 704 genera (COI - 95% accuracy, rbcL diatom - 90% accuracy) no cutoff needed for 200bp frag
# so we can work with ALL genera assignments with at 90-95% accuracy at least

# note Eisenia_Lumbricidae to just Eisenia (genus of segmented worms AND brown algae) 
a.1$Genus <- gsub("Eisenia_Lumbricidae", "Eisenia", a.1$Genus)

# get species or genus if not possible
a.1$taxon <- ifelse(a.1$Marker=="Diatom", 
                    ifelse(a.1$sBP >= 0.90, a.1$Species, a.1$Genus), #it''s rbcL
                    ifelse(a.1$sBP >= 0.70, a.1$Species, a.1$Genus)) # it's COI
# create taxonList of target species + genera + breakdown
a.1$taxon <- gsub("_", " ", a.1$taxon)
# needs manual editing because some are varieties, 'cf', 'sp', etc.
a.1$taxon <- gsub("Gomphonema pumilum var rigidum", "Gomphonema pumilum", a.1$taxon)
a.1$taxon <- gsub("Nitzschia dissipata var media", "Nitzschia dissipata", a.1$taxon)
a.1$taxon <- gsub("Nitzschia cf pusilla", "Nitzschia", a.1$taxon)
a.1$taxon <- gsub("Nitzschia cf bulnheimiana", "Nitzschia", a.1$taxon)
a.1$taxon <- gsub("Caloneis sp NE-S01", "Caloneis", a.1$taxon)
a.1$taxon <- gsub("Nitzschia cf microcephala", "Nitzschia", a.1$taxon)
a.1$taxon <- gsub("Stegopterna mutata/diplomutata", "Stegopterna", a.1$taxon)
a.1$taxon <- gsub("Fallacia sp", "Fallacia", a.1$taxon)
a.1$taxon <- gsub("Discostella sp", "Discostella", a.1$taxon)
a.1$taxon <- gsub("Tryblionella sp", "Tryblionella", a.1$taxon)
taxonList <- unique(a.1$taxon)
# 911 taxa
# count number of taxa that are genera (no spaces) or species (has space)
length(grep(" ", taxonList))
# 279 species
length(grep(" ", taxonList, invert=TRUE))
# 632 genera

# sanity check, try another way just in case
# finds first word (one word)
length(unique(word(taxonList, 1))) #911-279 = 632 genera
#finds second word, filters out NAs, count (two words)
length(taxonList[!is.na(word(taxonList, 2))]) #279 species
speciesList <- taxonList[!is.na(word(taxonList, 2))]

# #################
# First run.  Be patient, it takes a while
interactions <- list()
for (i in 1:length(taxonList)) {
  interactions[[i]]<- get_interactions_by_taxa(sourcetaxon = taxonList[[i]],
                                               showfield = c(
                                                 "source_taxon_name",
                                                 "source_taxon_path",
                                                 "interaction_type",
                                                 "target_taxon_name",
                                                 "target_taxon_path"
                                                            )

                                              )
  Sys.sleep(1)
}

# save to file so don't have to re-run globi query again later
df <- do.call(rbind.data.frame, interactions)
write.csv(df, "new_globi.csv", row.names = FALSE)
################

# use new file from above for a new search, or the original file to reproduce manuscript results
df <- read.csv("original_globi.csv", header=TRUE)

# keep species results (create mapping)
#speciesList <- taxonList[str_count(taxonList, '\\w+') == 2] #281
df2 <- df[df$source_taxon_name %in% speciesList |
            df$target_taxon_name %in% speciesList,]

# keep genera results (merge across species)
df3 <- df[!(df$source_taxon_name %in% speciesList |
            df$target_taxon_name %in% speciesList),]
# collapse names to genus rank to pool results across all species in genus
df3$source_taxon_name <- word(df3$source_taxon_name,1)
df3$target_taxon_name <- word(df3$target_taxon_name,1)

df_thisOne <- rbind(unique(df2), unique(df3))

# figure out missing species to collapse to genus rank for second search
uniqueSpeciesInResults <- unique(cbind(df2$source_taxon_name, df2$target_taxon_name))
speciesFound <- speciesList[speciesList %in% uniqueSpeciesInResults]
length(unique(speciesFound))
# 105
speciesNotFound <- speciesList[!speciesList %in% uniqueSpeciesInResults]
length(unique(speciesNotFound))
# 174
generaToFind <- unique(word(speciesNotFound, 1))
length(unique(generaToFind))
# 112

# figure out number of Genera found from the first search
uniqueGeneraInResults <- unique(cbind(df3$source_taxon_name, df3$target_taxon_name))
generaList <- taxonList[!str_count(taxonList, '\\w+') == 2]
generaFound <- generaList[generaList %in% uniqueGeneraInResults]
length(unique(generaFound))
# 412
generaNotFound <- generaList[!generaList %in% uniqueGeneraInResults]
length(unique(generaNotFound))
# 219

# #################
# # Second run.  Be patient, it takes a while
 interactions2 <- list()
 for (i in 1:length(generaToFind)) {
  interactions2[[i]]<- get_interactions_by_taxa(sourcetaxon = generaToFind[[i]],
                                               showfield = c(
                                                 "source_taxon_name",
                                                 "source_taxon_path",
                                                 "interaction_type",
                                                 "target_taxon_name",
                                                 "target_taxon_path"
                                                            )

                                              )
  Sys.sleep(1)
}

# save to file so don't have to re-run globi query again later
df4 <- do.call(rbind.data.frame, interactions2)
write.csv(df4, "new_globi2.csv", row.names = FALSE)
################

# use file from above for a new search, or original file to reproduce results from the manuscript
df4 <- read.csv("original_globi2.csv", header=TRUE)

# collapse names to genus rank to pool results across all species in genus
df4$source_taxon_name <- word(df4$source_taxon_name,1)
df4$target_taxon_name <- word(df4$target_taxon_name,1)

# put together all genera found
uniqueGeneraInResults2 <- unique(cbind(df4$source_taxon_name, df4$target_taxon_name))
generaFound2 <- generaToFind[generaToFind %in% uniqueGeneraInResults2]
length(unique(generaFound2))
# 90
generaNotFound2 <- generaToFind[!generaToFind %in% uniqueGeneraInResults2]
length(unique(generaNotFound2))
# 22

# put search results 1 and 2 together
df5 <- unique(rbind(df_thisOne, df4))

# number of unique taxa used to search GloBI
taxonList2 <- c(speciesFound, generaToFind, generaFound, generaNotFound)
length(unique(taxonList2))
# 777 unique target taxa

# number of these unique taxa found in the two-step search
length(unique(c(speciesFound, generaFound, generaFound2)))
# 548 71%

# create mapping file
# here figure out which taxa are macroinvertebrates, the result by default will be diatoms
inverts <- unique(a.1$taxon[!a.1$Phylum=="Bacillariophyta"])
# also collapse to genus to make sure I don't miss any
inverts.genus <- word(inverts, 1)
inverts.all <- unique(c(inverts, inverts.genus))
# manually edit to account for unresolved species

# for inverts (from results.csv) in taxonList2 (all unique taxa used to search globi)
inverts2 <- inverts.all[inverts.all %in% taxonList2]

# edit original to reflect what was in the results csv
mapping <- data.frame(original=c(speciesFound, speciesNotFound, generaFound, generaNotFound), 
                      new=c(speciesFound, word(speciesNotFound, 1), generaFound, generaNotFound))
mapping$taxon <- ifelse(mapping$new %in% inverts2, "Macroinvertebrate", "Diatom")

# create infiles for cheddar
# first get a list of genera in each site
# for that list, keep relevant interactions
sites <- unique(a.1$Site)
all <- list()
for (i in 1:length(sites)) {
  # process sites one at a time
  s <- a.1[a.1$Site==sites[[i]],]
  keep <- data.frame(original=unique(s$taxon)) # query species & genera
  keepNew.df <- merge(keep, mapping, by="original", all.x = TRUE)
  keepNew <- unique(keepNew.df$new)
  
  # subset by taxa in site
  # ensure source and target are both found in the site
  all[[i]] <- unique(df5[df5$source_taxon_name %in% keepNew &
                         df5$target_taxon_name %in% keepNew,])
  all[[i]]$site <- sites[[i]]
  
}
# put it together
allRecords <- do.call(rbind.data.frame, all)

# format for cheddar
df.list <- list()
for (i in 1:length(sites)) {
  allRecords.site <- allRecords[allRecords$site==sites[[i]],]
  df.list[[i]] <- data.frame(resource=allRecords.site$source_taxon_name, 
                             interaction_type=allRecords.site$interaction_type, 
                             consumer=allRecords.site$target_taxon_name,
                             site=sites[[i]])
}
trophic.links <- unique(do.call(rbind.data.frame, df.list))
length(unique(c(trophic.links$resource, trophic.links$consumer)))
#266, i.e. 266/777 = 34% retained for further analysis after filtering out interactions with non-target taxa

# filter by interaction type in the resource -> consumer direction
keepInteractions <- c("eatenBy", "preyedUponBy",
                      "hostOf", "hasParasitoid", "hadEctoparasite",
                      "hasEndoparasitoid", "hasEctoparasitoid",
                      "vectorOf", "livedInsideOfBy", "hasKleptoparasite",
                      "dispersalVectorOf")
trophic.links <- trophic.links[trophic.links$interaction_type %in% keepInteractions,]
length(unique(c(trophic.links$resource, trophic.links$consumer)))
# 171 / 777 = 22% retained after filtering out off-target interactions

# see how many taxa retained now
resource.consumer2 <- unique(c(trophic.links$resource, trophic.links$consumer)) #171
# so we found interactions for 77.5% of target taxa in Globi

# L7
original.L7.taxa.df <- data.frame(original=unique(a.1$taxon[a.1$Site=="COWL07"]))
new.L7.taxa <- merge(original.L7.taxa.df, mapping, by.x="original", by.y="original", all.x = TRUE)
length(unique(new.L7.taxa$new))
# 407 original taxa
# now figure out how many interactions were retained for site
countL7.df <- trophic.links[trophic.links$site=="COWL07",]
# count number of unique taxa represented
length(unique(c(countL7.df$resource, countL7.df$consumer)))
# 105 included in network
# 25.8% unique taxa represented in network

# C12
original.C12.taxa.df <- data.frame(original=unique(a.1$taxon[a.1$Site=="COWC12"]))
# 477 original unique taxa, but need to consider species that were collapsed to genera
new.C12.taxa <- merge(original.C12.taxa.df, mapping, by.x="original", by.y="original", all.x = TRUE)
length(unique(new.C12.taxa$new))
# 407 original taxa
# now figure out how many interactions were retained for site
countC12.df <- trophic.links[trophic.links$site=="COWC12",]
# count number of unique taxa represented
length(unique(c(countC12.df$resource, countC12.df$consumer)))
# 113 included in network
# 27.8% unique taxa represented in network

# B18
original.B18.taxa.df <- data.frame(original=unique(a.1$taxon[a.1$Site=="COWB18"]))
# 315 original unique taxa, but need to consider species that were collapsed to genera
new.B18.taxa <- merge(original.B18.taxa.df, mapping, by.x="original", by.y="original", all.x = TRUE)
length(unique(new.B18.taxa$new))
# 267 original taxa
# now figure out how many interactions were retained for site
countB18.df <- trophic.links[trophic.links$site=="COWB18",]
# count number of unique taxa represented
length(unique(c(countB18.df$resource, countB18.df$consumer)))
# 86 included in network
# 32.2% unique taxa represented in network

# C15
original.C15.taxa.df <- data.frame(original=unique(a.1$taxon[a.1$Site=="COWC15"]))
# 193 original unique taxa, but need to consider species that were collapsed to genera
new.C15.taxa <- merge(original.C15.taxa.df, mapping, by.x="original", by.y="original", all.x = TRUE)
length(unique(new.C15.taxa$new))
# 178 original taxa
# now figure out how many interactions were retained for site
countC15.df <- trophic.links[trophic.links$site=="COWC15",]
# count number of unique taxa represented
length(unique(c(countC15.df$resource, countC15.df$consumer)))
# 48 included in network
# 27.0% unique taxa represented in network



# reorder columns
trophic.links <- trophic.links[,c("resource","consumer","interaction_type","site")]
trophic.links$interaction.type <- NULL
trophic.links <- unique(trophic.links)

# Create Beaver18 files for Cheddar (limited interaction types for now)
Beaver18.t.l <- trophic.links[trophic.links$site=="COWB18",]
Beaver18.t.l$site <- NULL
b18 <- data.frame(node=unique(c(Beaver18.t.l$resource, Beaver18.t.l$consumer)), 
                  category="")
b18 <- unique(b18)
b18$category <- ifelse(b18$node %in% inverts2, "invertebrate", "producer")
b18.p <- data.frame(title="Beaver18")
dir.create("waterloo")
dir.create("waterloo/communities")
dir.create("waterloo/communities/Beaver18")
write.csv(b18, "waterloo/communities/Beaver18/nodes.csv", row.names = FALSE, quote=FALSE)
write.csv(Beaver18.t.l, "waterloo/communities/Beaver18/trophic.links.csv", row.names = FALSE, quote=FALSE)
write.csv(b18.p, "waterloo/communities/Beaver18/properties.csv", row.names = FALSE, quote=FALSE)

# Create Clair12 files for Cheddar
Clair12.t.l <- trophic.links[trophic.links$site=="COWC12",]
Clair12.t.l$site <- NULL
c12 <- data.frame(node=unique(c(Clair12.t.l$resource, Clair12.t.l$consumer)), 
                  category="")
c12 <- unique(c12)
c12$category <- ifelse(c12$node %in% inverts2, "invertebrate", "producer")
c12.p <- data.frame(title="Clair12")
dir.create("waterloo/communities/Clair12")
write.csv(c12, "waterloo/communities/Clair12/nodes.csv", row.names = FALSE, quote=FALSE)
write.csv(Clair12.t.l, "waterloo/communities/Clair12/trophic.links.csv", row.names = FALSE, quote=FALSE)
write.csv(c12.p, "waterloo/communities/Clair12/properties.csv", row.names = FALSE, quote=FALSE)

# Create Clair15 files for Cheddar
Clair15.t.l <- trophic.links[trophic.links$site=="COWC15",]
Clair15.t.l$site <- NULL
c15 <- data.frame(node=unique(c(Clair15.t.l$resource, Clair15.t.l$consumer)), 
                  category="")
c15 <- unique(c15)
c15$category <- ifelse(c15$node %in% inverts2, "invertebrate", "producer")
c15.p <- data.frame(title="Clair15")
dir.create("waterloo/communities/Clair15")
write.csv(c15, "waterloo/communities/Clair15/nodes.csv", row.names = FALSE, quote=FALSE)
write.csv(Clair15.t.l, "waterloo/communities/Clair15/trophic.links.csv", row.names = FALSE, quote=FALSE)
write.csv(c15.p, "waterloo/communities/Clair15/properties.csv", row.names = FALSE, quote=FALSE)

# Create laurel7 files for Cheddar
Laurel7.t.l <- trophic.links[trophic.links$site=="COWL07",]
Laurel7.t.l$site <- NULL
l7 <- data.frame(node=unique(c(Laurel7.t.l$resource, Laurel7.t.l$consumer)), 
                  category="")
l7 <- unique(l7)
l7$category <- ifelse(l7$node %in% inverts2, "invertebrate", "producer")
l7.p <- data.frame(title="Laurel7")
dir.create("waterloo/communities/Laurel7")
write.csv(l7, "waterloo/communities/Laurel7/nodes.csv", row.names = FALSE, quote=FALSE)
write.csv(Laurel7.t.l, "waterloo/communities/Laurel7/trophic.links.csv", row.names = FALSE, quote=FALSE)
write.csv(l7.p, "waterloo/communities/Laurel7/properties.csv", row.names = FALSE, quote=FALSE)



