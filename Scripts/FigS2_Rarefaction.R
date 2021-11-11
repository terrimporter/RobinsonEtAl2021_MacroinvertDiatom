# Teresita M. Porter, November 10, 2021

library(stringr) # str_split
library(reshape2) # dcast
library(vegan) # rarecurve
library(purrr) # for map_dfr
library(ggplot2) # ggplot
library(cowplot) # get_legend
library(scales) # comma

###################################################################
# Edit rarecurve function to remove the horizontal lines
###################################################################

rarecurve2 <- function (x, step = 1, sample, xlab = "Sample Size", ylab = "Species", 
                        label = TRUE, col, lty, ...) 
{
  x <- as.matrix(x)
  if (!identical(all.equal(x, round(x)), TRUE)) 
    stop("function accepts only integers (counts)")
  if (missing(col)) 
    col <- par("col")
  if (missing(lty)) 
    lty <- par("lty")
  tot <- rowSums(x)
  S <- specnumber(x)
  if (any(S <= 0)) {
    message("empty rows removed")
    x <- x[S > 0, , drop = FALSE]
    tot <- tot[S > 0]
    S <- S[S > 0]
  }
  nr <- nrow(x)
  col <- rep(col, length.out = nr)
  lty <- rep(lty, length.out = nr)
  out <- lapply(seq_len(nr), function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) 
      n <- c(n, tot[i])
    drop(rarefy(x[i, ], n))
  })
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)
  plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab, 
       type = "n", ...)
  if (!missing(sample)) {
    abline(v = sample) 
    rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"), 
                                           y = z, xout = sample, rule = 1)$y)
    #    abline(h = rare, lwd = 0.5) #turn off horizontal lines
  }
  for (ln in seq_along(out)) {
    N <- attr(out[[ln]], "Subsample")
    lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)
  }
  if (label) { 
    ordilabel(cbind(tot, S), labels = rownames(x), ...)
  }
  invisible(out)
}

#####################################################################

# Read in cat.csv
a <- read.table("diatom_invert_combined.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE)

# total number of ESVs
length(unique(a$GlobalESV))
# 4026

# total number of reads
sum(a$ESVsize)
# 1,304,473

# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
a.1 <- data.frame(a, do.call(rbind, str_split(a$SampleName,"_")), stringsAsFactors = FALSE)
names(a.1)[33:37] <- c("Name","Marker","Site","Replicate","Illumina Sample")

# pivot to make esv matrix (pool across versions, keep only substrate + sites separate)
A.esv <- reshape2::dcast(a.1, SampleName ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# move sample to rownames then delete
rownames(A.esv) <- A.esv$SampleName
A.esv$SampleName <- NULL

#remove columns with only zeros
esv.notnull <- A.esv[,colSums(A.esv) !=0]

#remove rows with only zeros & edit rownames
esv.notnull2 <- esv.notnull[rowSums(esv.notnull) !=0,]

# check read depth per sample - ok
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

# Do rarefection with pkg 'vegan'
rarecurveout <- rarecurve2(esv.notnull2, 
                           sample=esv.percentile, 
                           step=250, 
                           label=T)

# Reformat vegan list as df (cols OTU, raw.read)
rare.df <- lapply(rarecurveout, function(x){
  b <- as.data.frame(x)
  b <- data.frame(OTU = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})

# Add sample names to vegan output (df) (rownames)
sample_names <- rownames(esv.notnull2)
names(rare.df) <- sample_names

# Map rownames to vegan output (df)
rare.df <- map_dfr(rare.df, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")

# Parse out metadata from sample
rare.df <- data.frame(rare.df, do.call(rbind, str_split(rare.df$sample,"_")), stringsAsFactors = FALSE)
names(rare.df)[4:8]<-c("Name","Marker","Site","Replicate","IlluminaSample")

# create factors
rare.df$Marker <- factor(rare.df$Marker,
                         levels = c("Diatom", "COI"),
                         labels = c("rbcL", "COI"))
rare.df$Site <- factor(rare.df$Site, 
                       levels = c("COWB18", "COWC12", "COWC15","COWL07"),
                       labels = c("Beaver18", "Clair12", "Clair15","Laurel7"))
rare.df$Replicate <- factor(rare.df$Replicate, 
                         levels = c("1", "2", "3"),
                         labels = c("Replicate 1", "Replicate 2", "Replicate 3"))

# color by site
p1 <- ggplot(data = rare.df) +
  ggtitle("Sites") +
  labs(x="Reads", y="ESVs") +
  geom_point(aes(x = raw.read, y = OTU, color = Site), size=0.1) +
  geom_vline(xintercept = esv.percentile, linetype = "dashed") +
  scale_x_continuous(label = comma) +
  scale_color_viridis_d()+
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle=90),
        text = element_text(size=10),
        plot.title = element_text(size=10),
        legend.title = element_blank(),
        legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size=5)))
p1

# color by marker
p2 <- ggplot(data = rare.df) +
  ggtitle("Markers") +
  labs(x="Reads", y="ESVs") +
  geom_point(aes(x = raw.read, y = OTU, color = Marker), size=0.1) +
  geom_vline(xintercept = esv.percentile, linetype = "dashed") +
  scale_x_continuous(label = comma) +
  scale_color_viridis_d()+
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=10),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle=90),
        plot.title = element_text(size=10),
        legend.title = element_blank(),
        legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size=5)))
p2

# color by site status
sites <- read.csv("Sites.csv", header=TRUE)
rare.df2 <- merge(rare.df, sites, by.x="Site", by.y="site", all.x=TRUE)
names(rare.df2)[11] <- "Site status"
p3 <- ggplot(data = rare.df2) +
  ggtitle("Site status") +
  labs(x="Reads", y="ESVs") +
  geom_point(aes(x = raw.read, y = OTU, color = `Site status`), size=0.1) +
  geom_vline(xintercept = esv.percentile, linetype = "dashed") +
  scale_x_continuous(label = comma) +
  scale_color_viridis_d()+
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=10),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle=90),
        plot.title = element_text(size=10),
        legend.title = element_blank(),
        legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size=5)))
p3

#combine plots
g <- plot_grid(p1, p2, p3)
g
ggsave("FigS2_Rarefaction.jpeg", g, width = 8, height = 8)
# vertical line represents 15% percentile of read depth across sample (excluding controls)
