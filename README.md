# epigbs

## Scripts for analysis of McNew et al. 2020 
### Methylation of Galapagos mockingbird and zebra finch nestlings


### In this repo: 

- slurm scripts for alignment and variant calling 
- R scripts for testing for differential methylation and plotting

## Tools used 
**Command line**  
bwameth Version: 0.2.2

MethylDackel: A tool for processing bisulfite sequencing alignments.
Version: 0.3.0-3-g084d926-dirty (using HTSlib version 1.2.1)

Other dependencies: samtools, python, bwa

**R**
```
> sessionInfo()
R version 3.6.2 (2019-12-12)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Catalina 10.15.7

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
 [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] dplyr_1.0.1            stringr_1.4.0          scales_1.1.1           ggplot2_3.3.2          rtracklayer_1.46.0     GenomicFeatures_1.38.2
 [7] AnnotationDbi_1.48.0   Biobase_2.46.0         genomation_1.18.0      chromPlot_1.14.0       biomaRt_2.42.1         methylKit_1.12.0      
[13] GenomicRanges_1.38.0   GenomeInfoDb_1.22.1    IRanges_2.20.2         S4Vectors_0.24.4       BiocGenerics_0.32.0   
```

### Notes: 

- Mockingbird and zebra finch datasets were analyzed separately but using the same methods. 
- Demultiplexed fastqs have been deposited in the NCBI SRA (F and R for each individual). 
