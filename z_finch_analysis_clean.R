#Analysis 2020:
#Process so far on the cluster:
#1) Quality control using TrimGalore,
#2) Demultiplex and assemble ref with epiGBS scripts,
#3) Align using bwa-meth amd call variants using methyldackl
#4) Export snp and methylation files from cluster.

#Goals for this script: analyze zebra finches

# Packages --------------------------------------------------------------------
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#library(BiocManager)
#install_github("al2na/methylKit", build_vignettes=FALSE,
#               repos=BiocManager::repositories(),
#               dependencies=TRUE)
#BiocManager::install("genomation")
#BiocManager::install("GenomicFeatures")
#BiocManager::install("rtracklayer")
#BiocManager::install("GenomicRanges")

library(methylKit)
library(genomation)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(stringr)
library(chromPlot)
library(dplyr)

# Data --------------------------------------------------------------------
finchpredict <-
  read.csv("./processed_data/finch_pheno.csv") #finch metadata
ref_metadata <-
  read.csv("Reference_data/z_finch_name_bTaeGut2.pat.W.v2.csv")

# Methylkit ---------------------------------------------------------------
finch <- finchpredict$sample %>% as.character()
fl <- "cluster_output_finch_June2020/methyld_methylkit/"
file.list <-
  list() #create a list of files, one per sample, from methyldackl
for (i in 1:length(finch)) {
  file.list[[i]] <- paste(fl, finch[i], "_CpG.methylKit", sep = "")
}

myobj <- methRead (
  file.list,
  sample.id = as.list(finch),
  assembly = "z_finch",
  treatment = finchpredict$treatment,
  context = "CpG",
  mincov = 10
)

#some summary stats
for (i in 1:11) {
  getCoverageStats(myobj[[i]], plot = T, both.strands = F)
}
table(finchpredict$treatment)
sapply(myobj, nrow) %>% fivenum() #Number of smps called for all samples
head(myobj[[1]])
sapply(1:11, function(x)
  min(myobj[[x]]$coverage)) # mincov arg in methRead

#how many SMPs with coverage > 20?
sapply(1:11, function(x)
  myobj[[x]]$coverage[myobj[[x]]$coverage > 20] %>% length())


# Filtering coverage: depth at least 10, discards bases that have too much coverage
filtered.myobj <-
  filterByCoverage(
    myobj,
    lo.count = 10,
    lo.perc = NULL,
    hi.count = NULL,
    hi.perc = 99.9
  )

#normalize coverage
filtered.myobj <- normalizeCoverage(filtered.myobj)

# Merge samples, make sure that we have at least 4 per group.
# Decided to destrand because focus is on CpG methylation (see documentation)
meth2 <- unite(filtered.myobj,
               destrand = T,
               min.per.group = 4L)

dim(meth2)  # number of SMPs seen in at least 4 per group
percMethylation(meth2) %>% apply(1, function(x)
  sd(x, na.rm = T)) %>% mean() #mean variance = 1059
percMethylation(meth2) %>% apply(1, function(x)
  (max(x, na.rm = T) - min(x, na.rm = T))) %>%  .[. > 0] %>% length()
meth2 %>% data.frame %>% select(chr) %>% unique %>% dim

#PCA
methpca <- PCASamples(meth2, obj.return = T)
summary(methpca)
head(data.frame(methpca$x))

#Pca plot doesn't distinguish treatments
pdf(
  "output_plots/pca_finch.pdf",
  height = 4,
  width = 8,
  useDingbats = F
)
ggplot(data.frame(methpca$x)) +
  aes(x = PC1,
      y = PC2,
      color = finchpredict$treatment) +
  geom_point(size = 3,) +
  scale_color_manual(values = c("#DDAA33", "#4477AA")) +
  theme_classic(base_size = 15) +
  stat_ellipse()
dev.off()
#Are PCA axes significantly associated with covariates?


# Calculates difference of water - permethrin (control according to methylkit)
myDiff <- calculateDiffMeth(meth2,
                            #covariates= dplyr::select(finchpredict, nest),
                            overdispersion = "MN", test = "Chisq")

diffMethPerChr(
  myDiff,
  plot = TRUE,
  qvalue.cutoff = 0.01,
  meth.cutoff = 25
)

myDiff25 <- getMethylDiff(myDiff, difference = 25, qvalue = 0.01)
nrow(myDiff25)
colnames(myDiff25)[1] <- "Accession"

sig_dif <-
  merge(getData(myDiff25),
        dplyr::select(ref_metadata, Chr, Accession),
        all.x = T)

write.csv(sig_dif, "results_output/sig_dif.csv")
write.table(
  sig_dif[, 1:4],
  "results_output/sig_dif.bed",
  sep = "\t",
  col.names = F,
  row.names = F,
  quote = F
)
# Mapping genes -----------------------------------------------------------
# format gff file
# cut -f9 GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.gff  | awk -F';' '{print $1, $0}' | awk '{gsub(/ID=/, "", $1)}1'| grep "gene"| awk '{gsub(/^.*gene=/, "", $2)}1' | awk '{gsub(/;.*$/, "", $2)}1' | cut -d ' ' -f1,2  > finch_gff_forR.txt
gene_key <- read.table("Reference_data/finch_gff_forR.txt")
names(gene_key) <- c("feature.name", "gene")
gene.obj <-
  readTranscriptFeatures(
    "Reference_data/GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.bed",
    remove.unusual = F,
    up.flank = 2000,
    down.flank = 0
  ) #default estimates promoters as 1000 bp downstream

# All sequenced regions

all_regions <-
  annotateWithGeneParts(as(myDiff, "GRanges"), gene.obj)

# Just diff meth regions
regions <- annotateWithGeneParts(as(myDiff25, "GRanges"), gene.obj)
dist_tss <-
  getAssociationWithTSS(regions) # target.row is the row number in myDiff25
head(regions@dist.to.TSS)

finch_anno <-
  regions@members %>% as.data.frame %>% mutate(target.row = 1:182) %>% merge(., regions@dist.to.TSS)
finch_anno <- merge(finch_anno, gene_key)

# numbers of genes associated with DMCs
finch_anno %>% filter(prom == 1) %>% select(gene) %>% unique %>% dim()
finch_anno %>% filter(exon == 1) %>% select(gene) %>% unique %>% dim()

finch_anno %>% filter(exon == 1 |
                        prom == 1 | intron == 1) %>% select(gene) %>% unique %>% dim()
finch_anno %>% filter(exon == 1 | prom == 1 | intron == 1) %>% dim()

finch_anno_merged <-
  myDiff25 %>% data.frame %>% mutate(target.row = 1:182) %>% merge(., finch_anno)
write.csv(finch_anno_merged,
          "results_output/finch_annotation_granges_nov.csv")


#Percentage of differentially methylated regions that overlap with regions etc.
getTargetAnnotationStats(regions, percentage = T, precedence = T)

# Randomization -----------------------------------------------------------


#Randomize treatments
test2 <- NULL

for (i in 1:1000) {
  randt <- sample(finchpredict$treatment,
                  size = 11,
                  replace = F)
  myobj_reorg <- reorganize(myobj, sample.ids = finchpredict$sample,
                            treatment = randt) %>%
    filterByCoverage(
      .,
      lo.count = 10,
      lo.perc = NULL,
      hi.count = NULL,
      hi.perc = 99.9
    ) %>%
    normalizeCoverage(.) %>%
    unite(., destrand = T, min.per.group = 4L) %>%
    calculateDiffMeth(., overdispersion = "MN", test = "Chisq") %>%
    getMethylDiff(., difference = 25, qvalue = 0.01) %>%
    nrow() -> test2[i]
}

pdf("output_plots/randomization_hist_finch.pdf")
hist(test2,
     xlab = "Differentially methylated cytosines",
     col = "#EEDD88",
     main = "")
abline(v = nrow(myDiff25), lwd = 2)
dev.off()


test2[test2 >= nrow(myDiff25)] %>% length() # only 3 simulations out of 1000 (p, 0.003)
median(test2)
3 / 1000
head(finchmeth)



# Plots -------------------------------------------------------------------

finch_genome_coord <-
  filter(ref_metadata, Sequence.Role == "assembled-molecule") %>%
  mutate(Start = 1,
         End = Sequence.Length,
         Chrom = as.character(Chr)) %>%
  dplyr::select(Chrom, Start, End) %>% distinct(Chrom, .keep_all = T) %>%
  filter(Chrom != "na", Chrom != "MT", Chrom != "W")


head(myDiff)
sites <- myDiff %>% data.frame %>% mutate(Accession = chr) %>%
  merge(., dplyr::select(ref_metadata, Chr, Accession), all.x = T) %>%
  dplyr::select(Chrom = Chr, Start = start, End = end) %>% filter(Chrom != "na", Chrom !=
                                                                    "MT", Chrom != "W") %>%
  mutate(Chrom = as.character(Chrom))

sitesx2 <- rbind(sites, sites, sites, sites, sites)

diff <- myDiff25 %>% data.frame %>% mutate(Accession = chr) %>%
  merge(., dplyr::select(ref_metadata, Chr, Accession), all.x = T) %>%
  dplyr::select(Chrom = Chr, Start = start, End = end) %>% filter(Chrom != "na", Chrom !=
                                                                    "MT", Chrom != "W") %>%
  mutate(Chrom = as.character(Chrom))
head(diff)
diffx2 <- diff



write.csv(diff, "finch_diff.csv", row.names = F)
pdf("output_plots/chrom_zfinch.pdf")
chromPlot(
  gaps  = finch_genome_coord,
  cex = 1.5,
  #segment2 = sites, colSegments2 = "#332288",
  annot1 = diff,
  colAnnot1 = "#882255",
  colStat = "grey",
  chrSide = c(1, 1, 1, 1, -1, 1, 1, 1)
)
dev.off()
