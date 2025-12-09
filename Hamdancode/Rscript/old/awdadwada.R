# install.packages("karyoploteR")
# 1.  Make sure BiocManager is present
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# 2.  Install karyoploteR from the current Bioconductor release
BiocManager::install("karyoploteR")



install.packages("karyoploteR")
library(ggplot2)
library(karyoploteR)
library(dplyr)
head(var_cpg)
## ----- 1. tidy your data frame -----------------------------------------------
# var_cpg must have columns: Chr (1-22), Pos (bp), var_beta
vc <- var_cpg %>% 
  filter(Chr %in% 1:22) %>% 
  mutate(chrom = paste0("chr", Chr))

## ----- 2. start an hg38 ideogram ---------------------------------------------
kp <- plotKaryotype(genome = "hg38", plot.type = 2)   # 2 = all chromosomes stacked

## optional grey background grid
kpAddBaseNumbers(kp, tick.dist = 5e7, add.units = TRUE, tick.len = 0)

## ----- 3. scatter the points --------------------------------------------------
# normal points in light grey
kpPoints(kp,
         chr   = vc$chrom,
         x     = vc$Pos,
         y     = vc$var_beta,
         pch   = 20,
         cex   = 0.3,
         col   = "#BBBBBB")

# highlight the top 1 % in red
q99 <- quantile(vc$var_beta, 0.99)
kpPoints(kp,
         chr   = vc$chrom[vc$var_beta > q99],
         x     = vc$Pos [vc$var_beta > q99],
         y     = vc$var_beta[vc$var_beta > q99],
         pch   = 20,
         cex   = 0.5,
         col   = "red")

 


# Bioconductor installation, if needed
# if (!require("BiocManager")) install.packages("BiocManager")
 BiocManager::install("ggbio")

 nn
 # ─────────────────────────────  SET-UP  ──────────────────────────────
 # Install once if needed:
 # if (!requireNamespace("BiocManager", quietly = TRUE))
 #     install.packages("BiocManager")
 # BiocManager::install(c("ggbio", "GenomicRanges"))
 
 library(ggbio)
 library(GenomicRanges)
 library(ggplot2)
 
 # ── 1.  GRanges ----------------------------------------------------------------
 gr <- GRanges(seqnames = paste0("chr", var_cpg$Chr),
               ranges   = IRanges(start = var_cpg$Pos, width = 1),
               var_beta = var_cpg$var_beta)
 seqlevelsStyle(gr) <- "UCSC"           # 'chr1', 'chr2', …
 
 # small data-frame for the scatter layer
 df <- data.frame(seqnames = seqnames(gr),
                  start    = start(gr),
                  var_beta = mcols(gr)$var_beta)
 
 q99 <- quantile(df$var_beta, 0.99)     # 99th-percentile cut-off
 
 # ── 2.  PLOT -------------------------------------------------------------------
 ggplot() +
   layout_karyogram(gr, cytobands = FALSE) +   # chromosome bars only
   geom_point(data = df,                       # ← pass the data frame here
              aes(x = start, y = var_beta),
              colour = "grey60", size = 0.3) +
   geom_hline(yintercept = q99, colour = "red", linetype = 2) +
   facet_wrap(~ seqnames, ncol = 1, scales = "free_x") +
   labs(x = "Genomic position (bp)", y = "Variance of methylation (β)") +
   theme_bw() +
   theme(strip.text.y = element_text(angle = 0))
 
 
