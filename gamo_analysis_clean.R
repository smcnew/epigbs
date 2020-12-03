#Analysis 2020 mockingbird:
#Process so far on the cluster:
#1) Quality control using TrimGalore,
#2) Demultiplex and assemble ref with epiGBS scripts,
#3) Align and call variants using bwa-meth+methyldackl
#4) Export snp and methylation files from cluster.
#
#Goals for this script: analyze mockingbirds


# Packages --------------------------------------------------------------------
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#library(BiocManager)
#install_github("al2na/methylKit", build_vignettes=FALSE,
#               repos=BiocManager::repositories(),
#               dependencies=TRUE)
#BiocManager::install("genomation")
#BiocManager::install("GenomicFeatures")
#BiocManager::install("rtracklayer")
#BiocManager::install("GenomicRanges")
#BiocManager::install("chromPlot")

library(methylKit)
library(chromPlot)
library(genomation)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(ggplot2)
library(scales)
library(stringr)
library(dplyr)
# Data --------------------------------------------------------------------
# New output, June 2020

gamopredict <- read.csv("./processed_data/gamo_metadata.csv") #experiment metadata
genes <- read.csv("processed_data/genes_gamo_finch_gamode.csv", stringsAsFactors = F) # list of genes from knutie's experiment

# Methylkit ---------------------------------------------------------------


gamo <- gamopredict$sample %>% as.character()
fl <- "cluster_output_gamo_June2020/methyld_methylkit/"
file.list <- list() #create a list of files, one per sample, from methyldackl
for (i in 1:length(gamo)) {
  file.list[[i]] <- paste(fl, gamo[i], "_CpG.methylKit", sep="")}

myobj <- methRead (file.list,
                   sample.id = as.list(gamo),
                   assembly = "gamo",
                   treatment = gamopredict$treatment,
                   context = "CpG",
                   mincov = 10)

#some summary stats
for ( i in 1:71){
  getCoverageStats(myobj[[i]], plot=T, both.strands=F)
}

filtered.myobj <- filterByCoverage(myobj, lo.count=10, lo.perc=NULL, hi.count=NULL,
                                   hi.perc=99.9) %>% normalizeCoverage()
sapply(filtered.myobj, nrow) %>% fivenum() #Number of smps called for all samples

#how many SMPs with coverage > x ?
sapply(1:length(myobj_reorg), function(x) myobj_reorg[[x]]$coverage[myobj_reorg[[x]]$coverage > 3] %>% length())

#Some samples were sequenced really poorly. Let's remove samples with < X smps

gamopredict$smps <- sapply(filtered.myobj, nrow)

gamopredict %>% arrange(desc(smps))
smpm <- mean(gamopredict$smps) *.75 #Min number of smps we want ~ 19K
hist(gamopredict$smps)
table(gamopredict$treatment, gamopredict$year)


filtered_predict <- filter(gamopredict, smps > smpm, sample !="DC686")

myobj_reorg <- reorganize(filtered.myobj, sample.ids = filtered_predict$sample,
                          treatment = filtered_predict$treatment)

sapply(myobj_reorg, nrow) %>% max()

# Merge samples, make sure that we have at least 5 per group.
meth2 <- unite(myobj_reorg, destrand=TRUE, min.per.group = 5L)
percMethylation(meth2) %>% apply(1, function(x) sd(x, na.rm =T)) %>% mean() #mean variance = 1059

percMethylation(meth2) %>% apply(1, function(x) (max(x, na.rm = T) - min(x, na.rm =T))) %>%  .[.>0] %>% length() #sites 11419, number range = 0? 11117/11419

dim (meth2)
percMethylation(meth2) %>% head()
#PCA
methpca <- PCASamples(meth2, obj.return = T)
summary(methpca)
tail(methpca$x)

#Pca plot doesn't distinguish treatments

pdf("output_plots/gamo_pca.pdf", width = 8, height = 4, useDingbats = F )
ggplot(data.frame(methpca$x))+
  aes(x=PC1,y=PC2,color = as.factor(filtered_predict$treatment))+
  geom_point(size = 3) +
  scale_color_manual(values = c("#4477AA", "#DDAA33"))+
  theme_classic(base_size = 15)+
  stat_ellipse()
dev.off()

# Calculates difference of water - permethrin (control according to methylk)
myDiff <- calculateDiffMeth(meth2,
                            covariates= dplyr::select(filtered_predict, year),
                            overdispersion="MN", test= "Chisq")
myDiff25 <- getMethylDiff(myDiff, difference=25, qvalue=0.01)


myDiff25 %>% data.frame %>% write.csv(., "results_output/sig_dif_gamo.csv")
myDiff25 %>% data.frame %>% dplyr::select(chr, start, end, strand) %>%
  write.table(
    .,
    file = "results_output/sig_dif_gamo.bed",
    sep = "\t",
    col.names = F,
    row.names = F,
    quote = F
  )
write.csv(filtered_predict, "results_output/filtered_predict_gamo_metadata.csv",
          row.names = FALSE)
# Export all cpgs so we can compare with Knutie's gene expression results
myDiff %>% data.frame %>% dplyr::select(chr, start, end, strand) %>%
  write.table(
    .,
    file = "results_output/all_cpg_gamo.bed",
    sep = "\t",
    col.names = F,
    row.names = F,
    quote = F
  )

#

# Varying N per group  ----------------------------------------------------


# how much does the number per group matter?

unite(filtered.myobj,
      destrand = T,
      min.per.group = 6L) %>% #dim
  calculateDiffMeth(
    .,
    overdispersion = "MN",
    test = "Chisq",
    covariates = dplyr::select(filtered_predict, year)
  ) %>%
  getMethylDiff(., difference = 25, qvalue = 0.01) %>%
  nrow()

unite(filtered.myobj,
      destrand = T,
      min.per.group = 8L) %>% #dim()
  calculateDiffMeth(
    .,
    overdispersion = "MN",
    test = "Chisq",
    covariates = dplyr::select(filtered_predict, year)
  ) %>%
  getMethylDiff(., difference = 25, qvalue = 0.01) %>%
  nrow()

unite(filtered.myobj,
      destrand = T,
      min.per.group = 10L) %>%
  calculateDiffMeth(
    .,
    overdispersion = "MN",
    test = "Chisq",
    covariates = dplyr::select(filtered_predict, year)
  ) %>%
  getMethylDiff(., difference = 25, qvalue = 0.01) %>%
  nrow()



# Mapping genes -----------------------------------------------------------
# Create a gene map
# cut -f9 mmel_ns_up_onlygene.ipr.gff| awk -F';' '{print $1, $0}' | grep "Note=Similar to" | awk '{gsub(/^.*Similar/, "", $2)}1' | awk '{gsub(/:.*$/, "")}1' | sed 's/to //g' | sed 's/ID=//g' | awk -F '  ' '{print $1, $2}' > gamo_gene_listR.txt
# cat gamo_gene_listR.txt | awk '{print $1, $2}' > gamo_gene_listR_corrected.txt
gene_key <- read.table("Reference_data/gamo_gene_listR_corrected.txt")
head(gene_key)
names(gene_key) <- c("feature.name", "gene")
gene.obj <- readTranscriptFeatures("processed_data/gamo.bed",remove.unusual = F, up.flank = 2000, down.flank = 0)

# Annotation of all cpgs
all_regions <- annotateWithGeneParts(as(myDiff, "GRanges"), gene.obj)
all_regions
gamo_all_anno <- all_regions@members %>% as.data.frame %>% mutate(target.row = 1:nrow(all_regions@members)) %>% merge(., all_regions@dist.to.TSS)
gamo_all_anno <- merge(gamo_all_anno, gene_key, all.x = T)
head(gamo_all_anno)

#is there any overlap with knutie's genes?
gamo_all_anno %>% mutate(ingene = prom + exon + intron) %>%
  filter(ingene > 0) %>% select(gene) %>% unique() %>% dim()


gamo_all_anno$gene %in% genes$gene[genes$species=="gamo_ge"] %>% sum()
genes$gene[genes$species=="gamo_ge"] %in% gamo_all_anno$gene %>% sum()


# Annotation of significant cpgs
regions <- annotateWithGeneParts(as(myDiff25,"GRanges"), gene.obj)
regions@members %>% colSums()
filter(gamo_anno, exon ==1 ) %>% dim()

gamo_anno <- regions@members %>% as.data.frame %>% mutate(target.row = 1:196) %>% merge(., regions@dist.to.TSS)
gamo_anno <- merge(gamo_anno, gene_key, all.x = T )
gamo_anno[,4:6] %>% colSums()

gamo_anno %>% filter(intron ==1) %>% nrow()


gamo_anno_merged <- myDiff25 %>% data.frame %>% mutate(target.row = 1:196) %>% merge(., gamo_anno)
write.csv(gamo_anno_merged, "results_output/gamo_annotation_granges_nov.csv")

filter(gamo_anno, exon ==1)
dist_tss <- getAssociationWithTSS(regions) # target.row is the row number in myDiff25p

# Randomization -----------------------------------------------------------
#Randomize treatments
filtered_predict
test_4_per <- NULL

for (i in 1:1000){
  randt <- sample(filtered_predict$treatment, size = length(filtered_predict$treatment), replace =F )
  myobj_reorg <- reorganize(filtered.myobj, sample.ids = filtered_predict$sample,
                            treatment = randt) %>%
    filterByCoverage(., lo.count=10, lo.perc=NULL, hi.count=NULL,
                     hi.perc=99.9) %>%
    normalizeCoverage(.) %>%
    unite(., destrand=T, min.per.group = 4L) %>% #-> gamometh
    calculateDiffMeth(., overdispersion="MN", test= "Chisq",
                      covariates= dplyr::select(filtered_predict, year)) %>%
    getMethylDiff(., difference=25, qvalue=0.01) %>%
    nrow() -> test_4_per[i]
}
hist(test_4_per)

# median = 193, true = 194 out of 11,419

# Plots -------------------------------------------------------------------

# need chromosome lengths. Extract them from the ref genome fasta using the
# .fai file (cut -f1,2 gamo_ref.fa.fai > sizes.gamo.genome). If the .fai file
# doesn't exist it can be made with samtools faidx.

gamo_genome_coord <- read.table("processed_data/sizes.gamo.genome", col.names =
                                  c("Chrom", "End"), stringsAsFactors = F)
gamo_genome_coord <- mutate(gamo_genome_coord, Start = 1) %>%
  dplyr::select(.,Chrom, Start, End) %>%
  mutate(Chrom = str_replace(Chrom, "Mmel_Chr", "" ))


sites <- myDiff %>% data.frame %>% dplyr::select(chr, start, end) %>%
  dplyr::select(Chrom = chr, Start = start, End = end) %>%
  mutate(Chrom = str_replace(Chrom, "Mmel_Chr", "" ))

diff <- myDiff25 %>% data.frame %>% dplyr::select(chr, start, end) %>%
  dplyr::select(Chrom = chr, Start = start, End = end) %>%
  mutate(Chrom = str_replace(Chrom, "Mmel_Chr", "" ))
dim(diff)
dim(sites)

str(sites)
pdf("output_plots/chrom_gamo.pdf")
chromPlot(gaps = gamo_genome_coord,
          cex = 1.5,
          # annot1 = sites, colAnnot1 = "#332288",
          segment = diff, colSegments = "#882255",
          colStat = "grey",
          chrSide=c(1,1,1,1,-1,1,1,1))
dev.off()



# Overlap mockingbird finch  ----------------------------------------------
go <- read.csv("processed_data/go_finch_gamo.csv", stringsAsFactors = F)


filter(genes, species == "gamo_ge") %>% dim

tail(genes)
gamo_genes <- genes$gene[genes$species=="gamo"]
finch_genes <- genes$gene[genes$species=="finch"]
gamode_genes <- genes$gene[genes$species=="gamo_ge"] #differentially expressed genes
gamo_cpg_genes <- genes$gene[genes$species=="gamo_cpg"]
length(gamode_genes)
gamode_genes[gamode_genes %in% gamo_cpg_genes]
gamode_genes %in% gamo_genes %>% sum() # no overlap

go$go[go$species=="finch"] %in% go$go[go$species=="gamo"] # no overlap
go$go[go$species=="gamo"] %in% go$go[go$species=="finch"] # no overlap


