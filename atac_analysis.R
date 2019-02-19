## Analysis of DelToro ATAC-Seq data from human neutrophils 
## Date: 02.14.2019
## Author: Michael Chimenti
## Organism: hg19
## Aligners/peak callers: bwa-mem / genrich 
## Design:  7 samples, no replicates 
## Reps: 0

# 1. PMN 0h
# 2. PMN 8h
# 3. PMN + H. pylori 8h
# 4. PMN + killed H. pylori 8h
# 5. PMN + H. pylori 24h 
# 6. PMN 24h
# 7. PMN + killed H. pylori 24h

##########
## Imports
##########

setwd('~/iihg/ChIPSeq/allen_lab/del_toro_feb2019/')

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)
library(tidyverse)
library(ReactomePA)
library(DOSE)
#library(rtracklayer)

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream = 3000)

###########
## Function defs
###########

get_peaks_ucsc <- function(filename) {
    peaks <- readPeakFile(filename, header = FALSE)
    newStyle <- mapSeqlevels(seqlevels(peaks), "UCSC")
    newStyle <- newStyle[!is.na(newStyle)]
    peaks <- renameSeqlevels(peaks, newStyle)
    return(peaks)
}

write_narrowPeak <- function(granges_obj, samplename){
    granges_obj %>% 
    as.tibble() %>% 
    select(-c(width,strand)) %>% 
    filter(!grepl("GL*", seqnames)) %>%
    write.table(file = paste(samplename,"ucsc.narrowPeak", sep=''), 
                                 row.names=FALSE, quote=FALSE, sep='\t', col.names=FALSE)
}

width_histogram <- function(granges_obj, binsize=20){
    width_df <- as.data.frame(width(granges_obj))
    colnames(width_df) <- "region"
    ggplot(width_df) + geom_histogram(aes(x=region), bins = 100, binwidth = binsize) + 
    coord_cartesian(xlim = c(0,1000)) + 
    ggtitle("ATAC peak width (bp)")
    #ggtitle(sprintf("Sample %s ATAC peak width histogram", deparse(substitute(granges_obj))))
  }

make_plots <- function(granges_obj, tag_matrix){
  png(sprintf("%s_widthplot.png", deparse(substitute(granges_obj))), 1200, 1200, pointsize=20, res=150)
  p1 = width_histogram(granges_obj, 15)
  show(p1)
  dev.off()
  
  png(sprintf("%s_covplot.png", deparse(substitute(granges_obj))), 1200, 1200, pointsize=20, res=150)
  p2 = covplot(granges_obj, weightCol='V5')
  show(p2)
  dev.off()
  
  png(sprintf("%s_avgTssProf.png", deparse(substitute(granges_obj))), 1200, 1200, pointsize=20, res=150)
  p3 = plotAvgProf(tag_matrix, xlim=c(-3000,3000))
  show(p3)
  dev.off()
  
}


###########
## Analysis
###########

####################### PMN 0hrs
pmn0 <- get_peaks_ucsc('1_20190121000_S59_L006_sort.narrowPeak')
make_plots(pmn0, promoter)

pmn0_anno <- annotatePeak(pmn0, tssRegion=c(-3000,3000), TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoBar(pmn0_anno)
upsetplot(pmn0_anno)

pmn0_genes <- as.data.frame(pmn0_anno)$geneId
pmn0_path <- enrichPathway(pmn0_genes)
pmn0_dose <- enrichDO(gene = pmn0_genes)

#################### PMN 8hrs
pmn8 <- get_peaks_ucsc('2_20190121000_S60_L006_sort.narrowPeak')
pmn8_tag <- getTagMatrix(pmn8, windows=promoter)
make_plots(pmn8, pmn8_tag)

pmn8_anno <- annotatePeak(pmn8, tssRegion=c(-3000,3000), TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoBar(pmn8_anno)
upsetplot(pmn8_anno)

pmn8_genes <- as.data.frame(pmn8_anno)$geneId
pmn8_path <- enrichPathway(pmn8_genes)

png("pmn8_dotplot.png", 1500, 800, pointsize=20, res=150)
dotplot(pmn8_path)
dev.off()

################## PMN + Hpylori 8hrs
pmn8_hp <- get_peaks_ucsc('3_20190121000_S61_L006_sort.narrowPeak')
pmn8_hp_tag <- getTagMatrix(pmn8_hp, windows=promoter)
make_plots(pmn8_hp, pmn8_hp_tag)

pmn8_hp_anno <- annotatePeak(pmn8_hp, tssRegion=c(-3000,3000), TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoBar(pmn8_hp_anno)
upsetplot(pmn8_hp_anno)

pmn8_hp_genes <- as.data.frame(pmn8_hp_anno)$geneId
pmn8_hp_path <- enrichPathway(pmn8_hp_genes)



#################### PMN + killedHpylori 8hrs
pmn8_hpk <- get_peaks_ucsc('4_20190121000_S62_L006_sort.narrowPeak')
pmn8_hpk_tag <- getTagMatrix(pmn8_hpk, windows = promoter)
make_plots(pmn8_hpk, pmn8_hp_tag)

pmn8_hpk_anno <- annotatePeak(pmn8_hp, tssRegion=c(-3000,3000), TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoBar(pmn8_hpk_anno)
upsetplot(pmn8_hpk_anno)



#################### PMN + Hpylori 24hrs
pmn24_hp <- get_peaks_ucsc('5_20190121000_S63_L006_sort.narrowPeak')
pmn24_hp_tag <- getTagMatrix(pmn24_hp, windows=promoter)
make_plots(pmn24_hp, pmn24_hp_tag)

pmn24_hp_anno <- annotatePeak(pmn24_hp, tssRegion=c(-3000,3000), TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoBar(pmn24_hp_anno)
upsetplot(pmn24_hp_anno)

pmn24_hp_genes <- as.data.frame(pmn24_hp_anno)$geneId
pmn24_hp_path <- enrichPathway(pmn24_hp_genes)

png("pmn24_hp_dotplot.png", 1500, 800, pointsize=20, res=150)
dotplot(pmn24_hp_path)
dev.off()

##################### PMN 24hrs
pmn24 <- get_peaks_ucsc('6_20190121000_S64_L006_sort.narrowPeak')
pmn24_tag <- getTagMatrix(pmn24, windows=promoter)
make_plots(pmn24, pmn24_tag)

pmn24_anno <- annotatePeak(pmn24, tssRegion=c(-3000,3000), TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoBar(pmn24_anno)
upsetplot(pmn24_anno)

pmn24_genes <- as.data.frame(pmn24_anno)$geneId
pmn24_path <- enrichPathway(pmn24_genes)

png("pmn24_dotplot.png", 1500, 800, pointsize=20, res=150)
dotplot(pmn24_path)
dev.off()


##################### PMN + killedHpylori 24hrs

pmn24_hpk <- get_peaks_ucsc('7_20190121000_S65_L006_sort.narrowPeak')
pmn24_hpk_tag <- getTagMatrix(pmn24_hpk, windows = promoter)
make_plots(pmn24_hpk, pmn24_hpk_tag)

pmn24_hpk_anno <- annotatePeak(pmn24_hpk,tssRegion=c(-3000,3000), TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoBar(pmn24_hpk_anno)
upsetplot(pmn24_hpk_anno)

