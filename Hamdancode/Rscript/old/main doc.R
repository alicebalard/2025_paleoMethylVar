install.packages("ggplot2")
install.packages("openxlsx")
install.packages("tidyverse")
install.packages("VennDiagram")
install.packages("ggmap")
install.packages("plotly")
install.packages("dplyr")
install.packages("ggmap")
install.packages("tmaptools")
install.packages("car")
install.packages(" hrbrthemes")
install.packages("limma")
install.packages("Rtools")
install.packages("UpSetR")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
library(BiocManager)
library(Rtools)
library(dplyr)
library(tidyverse)
library(openxlsx)
library(ggplot2)
library(VennDiagram) 
library(ggmap)
library(plotly)
library(ggmap)
library(tmaptools)
library(car)
library(tmaptools)
library(scales)
library(UpSetR)   
library(limma)

library(tibble)

save.image("C:/Users/hamda/Desktop/UGI_thingy/data/.RData")

setwd("/Users/hamda/Desktop/UGI_thingy")
Mdata <- read.xlsx("data/gkac503_supplemental_files/SUPPLEMENTARY_TABLES_resubmission.xlsx", sheet = 6)
ITlocation <- read.delim("data/filereport_read_run_PRJEB32566_tsv.txt")
MChrNum <- Mdata$chr
chrpos <- Mdata$hg19.pos 

 





 
#the interesting ones are they also found in mdata
# List which cites are found most often in IT data, Than 
#start with best you can get than shrink down 
# work out which ones have the most variation / which ones occur the most .
#This code will asign all of the files in the folder to the ITdata"number"



#files <- list.files("../UGI_thingy/data/ALL_MethylMaps_Antonio2019", full.names = TRUE)
files <- list.files("../UGI_thingy/data/newmethyl/", full.names = TRUE)


files
for (i in seq_along(files)) {
  file_path <- files[i]
  df_name <- paste0("ITdata", i)
  
  # Check file is not empty
  if (file.info(file_path)$size == 0) {
    message("Skipping empty file: ", file_path)
    next
  }
  
  # Try reading the file safely
  tryCatch({
    if (grepl("\\.gz$", file_path)) {
      assign(df_name, read.delim(gzfile(file_path, header = FALSE, skip = 1)))
    } else {
      assign(df_name, read.delim(file_path, header = FALSE, skip = 1))
    }
  }, error = function(e) {
    message("Skipping file due to error: ", file_path)
    message("Error: ", e$message)
  })
}


ITdatamerge <- NULL

for (i in seq_along(files)) {
  df_name  <- paste0("ITdata", i)
  ITdata_i <- get(df_name)
  
  # add the source BEFORE binding
  ITdata_i$source <- df_name   # or use df_name if you prefer "ITdata1", etc.
  
  ITdatamerge <- rbind(ITdatamerge, ITdata_i)
}
 
table(ITdatamerge$source)
# Unsure what this does 
for (i in seq_along(files)) {
  df_name <- paste0("ITdata", i)
  ITdata_i <- get(df_name)
  print(paste("DataFrame:", df_name))
  print(names(ITdata_i))
  print(ncol(ITdata_i))
}

# this now contains all rows from ITdata1 to ITdata63
# UNZ is for unzip pls rember 
# ../../ for navigating files showing directory like ls  
#read gz for zipped files 
 
#since now this stage has worked the next step is to make a function where the
#data set is given and its compared to maria data

#from the data 
#hvCpG annotations from Illumina450K manifest and hvCpG clusters considered in our study.
#n.clust cpgs = number of hvCpGs that are within cluster		

#save.image("C:/Users/hamda/Desktop/UGI_thingy/data/dataaaa.RData")
#Mmatches <- data.frame() #run this when you want to clearn it 

#Insert all of the files, write a function where you can pick a file and compare ----
#how many of the chromasomal positions are the same between maria and italy for all the positions 


all_ITmatches$BaseID <- sub("\\..*$", "", all_ITmatches$V1)

count_df <- all_ITmatches %>%
  group_by(BaseID) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))

 
print(head(count_df, 10))

top_counts <- count_df %>% slice_max(Count, n = 20)

ggplot(top_counts, aes(x = reorder(BaseID, -Count), y = Count)) +
  geom_bar(stat = "identity") +
  labs(title = "Top 20 ITdata Groups by Count", x = "ITdata ID", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


### again 
 
all_ITmatches$BaseID <- sub("\\..*$", "", all_ITmatches$V1)

count_df <- all_ITmatches %>%
  group_by(BaseID) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))


print(head(count_df, 10))
#So essentially, In both of the data sets I counted how many CPGS were matching on each chhromasome
#After this I put it onto a graph and this were the results. chromasome 19 is weirdly having alot matching 
#even though it is one of the smallest chromasomes.Therefore one of the next steps is to do an anova to count 
#if there is a significant difference in the number of entries for each chromasome. 

#This was run for chroamsomes 20-16 and showed that the number of values for chromasomes was significantly increasing 
#depending on the lower chr number (lower = larger)
#Even in ancient DNA CPG sites are very high on the chromasomes 19 
#find the percentage of cpg coverage for chromosome 19 specifically and compare to all other chromasomes
 
#Is it that chromasome 19 has more genes compared to the others so it has a greater conc of cpg sites 
#possible that the genes are those more influenced by methylations
#https://www.sciencedirect.com/topics/immunology-and-microbiology/chromosome-19 This shows that it may be that it has more genes 
#cpg sites are found in areas which are more gene rich 
 
#put ITdata together, Merge is for if they have a column in common, It will also put all rows together 
#Rbind is for just mashing all of it together 


#### matching as a funciton for individulas and matching specifically for ITdatamerge 
blah <- function(datar) {
  
  Mmatches <- data.frame()
  ITmatches <- data.frame() 
  
  
  for (i in 1:nrow(datar)) {
    
    pos_i <- datar$V2[i]  # genomic position from HGdata
    
    
    if (pos_i %in% chrpos) {   #checking them one by one by eachtother
      
      ITmatches <- rbind(ITmatches, datar[i, ])
      
      # Find the matching in mdata chromasome pos and HG data pos, if same adding to the thing
      match_row <- Mdata[Mdata$hg19.pos == pos_i, ]
      
      # Add the matching Mdata row(s) to Mmatches 
      Mmatches <- rbind(Mmatches, match_row)
    }
  }
  
  return(list(ITmatches = ITmatches, Mmatches = Mmatches)) 
  
}
blah(ITdatamerge)
## Add a column for the ID for each person 
## list of data frames, keep column of samples # how many times each hvcpg apear as next step 

Mmatches <- data.frame()
ITmatches <- data.frame()

for (i in 1:nrow(ITdatamerge)) {
  ITpos <- ITdatamerge$V2[i]  # or use correct column name like 'hg19.pos' or similar
  if (ITpos %in% chrpos) {
    # Add current IT row
    ITmatches <- rbind(ITmatches, ITdatamerge[i, ])
    # Add matching Mdata row(s)
    match_row <- Mdata[Mdata$hg19.pos == ITpos, ]
    Mmatches <- rbind(Mmatches, match_row)
  }
}
 
#table(ITmatches$V1)

ITmatches$V6 <- paste(ITmatches$V1, ITmatches$V2, sep = "_")
#ITmatches
ITmatches$V1 <- ITmatches$V6 
ITmatches$V1 <- NULL 
ITmatches$V2 <- NULL 
ITmatches$V5 <- NULL

ITmatches
ITmatches$V6 <- NULL
ITmatches
  
TableITmatches <- table(ITmatches$V1)
TableITmatches <- as.data.frame(TableITmatches)
TableITmatches
TableITmatches <- TableITmatches[order(-TableITmatches$Freq), ]
head(TableITmatches)

#regression of 
head(ITdatamerge)
ITdatamerge$V6 <- paste(ITdatamerge$V1, ITdatamerge$V2, sep = "_")
head(ITdatamerge)
ITdatamerge$V1 <- ITdatamerge$V6 
ITdatamerge$V1 <- NULL 
ITdatamerge$V2 <- NULL 
ITdatamerge$V5 <- NULL
ITdatamerge


head(ITdatamerge)
ITtable <- table(ITdatamerge$V1)
ITtable
head(ITdatamerge$V1, X = 300)
hist(ITtable[ITtable>6])
ITdf <- as.data.frame(ITtable)

table(ITtable)
ggplot( data=ITdf, aes(x = Var1, y = Freq)) + geom_histogram( ) +   
head(ITtable)

 
#ggplot(ITdf), aes( x = Freq) + geom_histogram(binwidth = 1)  + scale_y_log10()
 

ITdf$Freq <- ITdf$Freq + 1  

ggplot(data =ITdf, aes(x = Freq) ) +
  geom_histogram(binwidth = 1)  +
  scale_y_log10()
 
#you should add the variance row before doing the changes to the table !!! Also need to add person row for ITdatamerge
 
 
library(ggplot2)
head(ITdf)
# Remove rows with Freq == 0 or NA
ITdf_filtered <- ITdf[ITdf$Freq > 0 & !is.na(ITdf$Freq), ]

ggplot(data = ITdf_filtered, aes(x = Freq)) +
  geom_histogram(binwidth = 1) +
  scale_y_log10(labels = label_comma())  +
  theme_minimal()

sum
more3pos <- names(ITtable[ITtable >= 3]) ## 3 threshold 

dfmore3 <- ITdatamerge[ITdatamerge$V1 %in% more3pos,]
nrow(dfmore3)
dfmore3
## more than 3 positions !!! 

row.names(dfmore3) <- NULL 
head(dfmore3)
 

head(ITmatches$V1)
ITmatchtable <- table(ITmatches$V1)
ITmatchtable <- as.data.frame(ITmatchtable)
head(ITmatchtable)

ITmatchtable <- ITmatchtable[order(-ITmatchtable$Freq), ]
ITmatchtable

 head(ITmatchtable)
hist(ITmatchtable$Freq)
hist(ITmatchtable$Freq[ITmatchtable$Freq > 1], )



ggplot(data = ITmatchtable, aes(x = Freq))  +
  geom_histogram(binwidth = 1)


  





library(ggplot2)
head(ITdatamerge)
testdatamerge <- ITdatamerge

testdatamerge$V5 <- testdatamerge$V3 / testdatamerge$V4

ggplot(subset(testdatamerge, V5 > 0 ), aes(x = V5)) +
  geom_histogram(binwidth = 0.05) +
  theme_minimal()

ggplot(testdatamerge, aes(x = V5)) +
  geom_histogram(binwidth = 0.05) +
  theme_minimal()

variance_results


ggplot(variance_results, aes(x = Variance )) +
  geom_histogram(binwidth = 0.0005) +
  theme_minimal() 

ggplot(subset(variance_results, Variance > 0), aes(x = Variance )) +
  geom_histogram(binwidth = 0.0005) +
  theme_minimal()   

p1 <- ggplot(variance_results, aes(x = Variance, y = Freq)) +
  geom_point(color = "#69b3a2") +
  geom_smooth(method = "lm",            
              se = FALSE,               
              color = "black",    
              linewidth = 1)           
 
 
 #addd the unique 



site_period <- testmerge %>% distinct(V1, period4ana = period4ana) 
site_period
 #list for the cpg things 

cpg_sets <- site_period %>% group_by(period4ana) %>%  
  summarise(CpGs = list(V1)) %>% 
  deframe() 

actual_upset <- fromList(cpg_sets)

upset(actual_upset, nsets = length(cpg_sets), 
      nintersects = NA,
       
      order.by = "freq")


#make the upset for the values present in atleast 3 people 

newtbl <- table(testmerge$V1) 
newdf1 <- as.data.frame(newtbl)  
newdf1 <- newdf1[newdf1$Freq > 3, ]
more3testmerge <- testmerge[testmerge$V1 %in% newdf1$Var1, ]
more3testmerge

site_period_1 <- more3testmerge %>% distinct(V1, period4ana = period4ana) 
site_period

#list for the cpg things 
cpg_sets <- site_period_1 %>% group_by(period4ana) %>%  
  summarise(CpGs = list(V1)) %>% 
  deframe() 

actual_upset <- fromList(cpg_sets)

upset(actual_upset, nsets = length(cpg_sets), 
      nintersects = NA,
      
      order.by = "freq")




 
 


#-------------------- Upset for presentation 1 and 2 


# 1) Keep only rows where the CpG is methylated at least once
meth_rows <- testmerge %>%
  filter(!is.na(V3), V3 > 0, !is.na(period4ana))

# 2) Build per-period sets of CpGs (unique V1 values) that are methylated ≥1 time
cpg_sets <- meth_rows %>%
  distinct(period4ana, V1) %>%                 # one CpG per period
  group_by(period4ana) %>%
  summarise(CpGs = list(unique(V1)), .groups = "drop") %>%
  deframe()                                    # named list: period -> vector of CpGs

# 3) UpSet input + plot
actual_upset <- fromList(cpg_sets)

upset(
  actual_upset,
  nsets        = length(cpg_sets),
  nintersects  = NA,
  order.by     = "freq"
)



# --- 1) Keep rows where CpG is methylated at least once (V3 > 0) ---
meth_rows <- testmerge %>%
  filter(!is.na(V3), V3 > 0, !is.na(period4ana))

# --- 2) Build per-period CpG sets (unique V1) ---
cpg_sets <- meth_rows %>%
  distinct(period4ana, V1) %>%
  group_by(period4ana) %>%
  summarise(CpGs = list(unique(V1)), .groups = "drop") %>%
  deframe()

# Drop empty sets (if any)
cpg_sets <- cpg_sets[lengths(cpg_sets) > 0]

# Order sets by size (largest first) for a cleaner layout
set_sizes <- sort(sapply(cpg_sets, length), decreasing = TRUE)
set_order <- names(set_sizes)

# --- 3) UpSet input + neat plot settings ---
actual_upset <- fromList(cpg_sets)
pdf("C:/Users/hamda/Desktop/UGI_thingy/fig3.pdf", width = 10, height = 5)

upset(
  actual_upset,
  sets              = set_order,
  keep.order        = TRUE,
  nsets             = length(cpg_sets),
  nintersects       = min(30, nrow(actual_upset)),  # cap to top 30 intersections
  order.by          = "freq",
  decreasing        = c(TRUE, FALSE),
  # Aesthetics
  main.bar.color    = "grey20",
  sets.bar.color    = "grey35",
  
  matrix.color      = "grey20",
  shade.color       = "grey90",
  shade.alpha       = 0.3,
  matrix.dot.alpha  = 0.7,
  text.scale        = 1.15,
  mb.ratio          = c(0.65, 0.35),
  mainbar.y.label   = "Intersecting CpGs (≥1 methylation)",
  sets.x.label      = "CpGs per period (≥1 methylation)",
  show.numbers      = "yes",
  number.angles     = 0
)

dev.off()

 


#-------------------- Upset for presentation 1 and 2 

















# make df all the hvcpg > 3 people, 4 people than 5 and 6 
# so that if rand pick in hvcpg do they 
# 10k increase **** 
# range for 3.(sumary in tbl) repeat for 4 5 6 

IT_MT_more3 <- data.frame(ITmatchtable$Var1, ITmatchtable$Freq[ITmatchtable$Freq >= 3]) 
IT_MT_more3 <- ITmatchtable$Freq[ITmatchtable$Freq >= 3] 
IT_MT_more3 <- subset(ITmatchtable, Freq >= 3)
write.xlsx(IT_MT_more3, file = "morethan3.xlsx")

 
 
methyl_matrix <- testmerge %>%
  select(SiteID = V1, SampleID = run_accession, Beta = V5) %>%
  pivot_wider(names_from = SampleID, values_from = Beta) %>%
  column_to_rownames("SiteID") %>%
  as.matrix()
 
sample_metadata <- testmerge %>%
  select(SampleID = run_accession, period4ana) %>%
  distinct() %>%
  filter(SampleID %in% colnames(methyl_matrix)) %>%
  arrange(match(SampleID, colnames(methyl_matrix)))

population <- factor(sample_metadata$period4ana)

design <- model.matrix(~ 0 + population)
colnames(design) <- levels(population)
rownames(design) <- sample_metadata$SampleID

# Sanity check
ncol(methyl_matrix) == nrow(design)  # should be TRUE
all(colnames(methyl_matrix) == rownames(design))  # should be TRUE

# Run limma

# Fix column names in design matrix
colnames(design) <- make.names(colnames(design))


library(limma)
fit <- lmFit(methyl_matrix, design)

fitcontrast_matrix <- makeContrasts(Imperial - Iron.Republic, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
resultsy <- topTable(fit2, adjust = "BH", number = Inf)

library(ggplot2)
resultsy

# Add significance flag
resultsy$Significant <- resultsy$adj.P.Val < 0.05

# Basic volcano plot
ggplot(resultsy, aes(x = logFC, y = -log10(P.Value), color = Significant)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("grey70", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Differential Methylation",
       x = "Log2 Fold Change (Imperial vs Iron/Republic)",
       y = "-log10(P-value)")

library(limma)
library(ggplot2)

# 1. Clean and order metadata
sample_metadata <- sample_metadata %>%
  filter(SampleID %in% colnames(methyl_matrix)) %>%
  arrange(match(SampleID, colnames(methyl_matrix)))

# 2. Build design matrix (fix column names for contrasts)
population <- factor(sample_metadata$period4ana)
design <- model.matrix(~ 0 + population)
colnames(design) <- make.names(levels(population))
rownames(design) <- sample_metadata$SampleID

# 3. Fit model once
fit <- lmFit(methyl_matrix, design)

# 4. Create all pairwise combinations
groups <- colnames(design)
pairs <- combn(groups, 2, simplify = FALSE)

# 5. Loop through each pair and make volcano plot
for (pair in pairs) {
  group1 <- pair[1]
  group2 <- pair[2]
  contrast_name <- paste0(group1, "_vs_", group2)
  
  # Create contrast
  contrast_matrix <- makeContrasts(contrasts = paste0(group1, " - ", group2), levels = design)
  
  # Fit and get results
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  results <- topTable(fit2, adjust = "BH", number = Inf)
  results$Significant <- results$adj.P.Val < 0.05
  
  # Volcano plot
  p <- ggplot(results, aes(x = logFC, y = -log10(P.Value), color = Significant)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("grey70", "red")) +
    theme_minimal() +
    labs(
      title = paste("Volcano Plot:", group1, "vs", group2),
      x = "log2 Fold Change",
      y = "-log10(P-value)"
    )
  
  print(p)  # Display plot
}



library(ggplot2)
library(dplyr)

# Assuming your metadata is in `testmerge` or similar
sample_qc <- testmerge %>%
  select(run_accession, period4ana, `%Endogenous`, Coverage) %>%
  distinct()  # One row per sample

# Plot 1: % Endogenous DNA
ggplot(sample_qc, aes(x = period4ana, y = `%Endogenous`, fill = period4ana)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Endogenous DNA by Time Period", y = "% Endogenous") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 2: Coverage
ggplot(sample_qc, aes(x = period4ana, y = Coverage, fill = period4ana)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Sample Coverage by Time Period") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

sample_qc

# Assuming your methylation matrix is `methyl_matrix` with CpGs as rows, samples as columns

 


# Calculate variance difference and p-value (F-test) for same CpG sites
var_test <- testmerge2 %>%
  filter(period4ana %in% period_pair) %>%
  group_by(SiteID) %>%
  nest() %>%
  mutate(
    p_val = map_dbl(data, ~ {
      d <- .
      v1 <- d$V5[d$period4ana == period_pair[1]]
      v2 <- d$V5[d$period4ana == period_pair[2]]
      if (length(v1) > 1 && length(v2) > 1) {
        var.test(v1, v2)$p.value
      } else NA_real_
    }),
    var_diff = map_dbl(data, ~ {
      d <- .
      var(d$V5[d$period4ana == period_pair[1]], na.rm = TRUE) -
        var(d$V5[d$period4ana == period_pair[2]], na.rm = TRUE)
    })
  ) %>%
  ungroup() %>%
  filter(!is.na(p_val)) %>%
  mutate(logp = -log10(p_val))

# Volcano plot for variance
ggplot(var_test, aes(x = var_diff, y = logp)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "red") +
  geom_vline(xintercept = c(-0.01, 0.01), linetype = "dashed", colour = "grey") +
  labs(
    x = "Variance Difference",
    y = "-log10(p-value)",
    title = paste("Volcano Plot: Variance (", period_pair[1], "vs", period_pair[2], ")")
  ) +
  theme_minimal()

View(testmerge)
testmerge <- testmerge %>%
  mutate(SiteID = paste0("chr", Chr, "_", Pos))
df <- testmerge

# Get all unique time periods
time_periods <- unique(df$period4ana)

# Generate all pairwise combinations
combinations <- combn(time_periods, 2, simplify = FALSE)

# Create output folder
dir.create("volcano_plots", showWarnings = FALSE)

# Loop for β-values (mean methylation)
for (pair in combinations) {
  p1 <- pair[1]
  p2 <- pair[2]
  pair_df <- df %>%
    filter(period4ana %in% c(p1, p2)) %>%
    group_by(SiteID) %>%
    nest() %>%
    mutate(t_result = map(data, ~ {
      dat <- .
      val1 <- dat$V5[dat$period4ana == p1]
      val2 <- dat$V5[dat$period4ana == p2]
      if(length(val1) > 1 && length(val2) > 1 && (sd(val1) > 0 || sd(val2) > 0)) {
        t.test(val1, val2)
      } else NA
    })) %>%
    filter(!is.na(t_result)) %>%
    mutate(p_val = map_dbl(t_result, ~ .x$p.value),
           logFC = map_dbl(t_result, ~ mean(.x$data[[1]]$V5[.x$data[[1]]$period4ana == p1], na.rm=TRUE) -
                             mean(.x$data[[1]]$V5[.x$data[[1]]$period4ana == p2], na.rm=TRUE)),
           neg_log_p = -log10(p_val),
           Significant = p_val < 0.05 & abs(logFC) > 0.1)
  
  # Plot
  ggplot(pair_df, aes(x = logFC, y = neg_log_p, color = Significant)) +
    geom_point(size = 1) +
    scale_color_manual(values = c("grey", "red")) +
    labs(title = paste("Volcano Plot: ", p1, "vs", p2),
         x = "log2 Fold Change",
         y = "-log10(P-value)") +
    theme_minimal() +
    ggsave(filename = paste0("volcano_plots/volcano_beta_", p1, "_vs_", p2, ".png"),
           width = 8, height = 6)
}


# Loop for Variance comparison
for (pair in combinations) {
  p1 <- pair[1]
  p2 <- pair[2]
  pair_df <- df %>%
    filter(period4ana %in% c(p1, p2)) %>%
    group_by(SiteID) %>%
    nest() %>%
    mutate(t_result = map(data, ~ {
      dat <- .
      val1 <- dat$var_site[dat$period4ana == p1]
      val2 <- dat$var_site[dat$period4ana == p2]
      if(length(val1) > 1 && length(val2) > 1 && (sd(val1) > 0 || sd(val2) > 0)) {
        t.test(val1, val2)
      } else NA
    })) %>%
    filter(!is.na(t_result)) %>%
    mutate(p_val = map_dbl(t_result, ~ .x$p.value),
           logFC = map_dbl(t_result, ~ mean(.x$data[[1]]$var_site[.x$data[[1]]$period4ana == p1], na.rm=TRUE) -
                             mean(.x$data[[1]]$var_site[.x$data[[1]]$period4ana == p2], na.rm=TRUE)),
           neg_log_p = -log10(p_val),
           Significant = p_val < 0.05 & abs(logFC) > 0.005)
  
  # Plot
  ggplot(pair_df, aes(x = logFC, y = neg_log_p, color = Significant)) +
    geom_point(size = 1) +
    scale_color_manual(values = c("grey", "red")) +
    labs(title = paste("Volcano Plot (Variance): ", p1, "vs", p2),
         x = "log2 Fold Change (variance)",
         y = "-log10(P-value)") +
    theme_minimal() +
    ggsave(filename = paste0("volcano_plots/volcano_var_", p1, "_vs_", p2, ".png"),
           width = 8, height = 6)
}

 
 
periods <- unique(testmerge$period4ana)
period_combinations <- combn(periods, 2, simplify = FALSE)

# Loop through all combinations
for (pair in period_combinations) {
  period1 <- pair[1]
  period2 <- pair[2]
  
  cat("Processing:", period1, "vs", period2, "\n")
  
  # Subset data
  subset_df <- testmerge %>%
    filter(period4ana %in% c(period1, period2)) %>%
    select(V1, period4ana, V5)
  
  # Group and test
  results <- subset_df %>%
    group_by(V1) %>%
    nest() %>%
    mutate(stats = map(data, ~ {
      d <- .
      v1 <- d$V5[d$period4ana == period1]
      v2 <- d$V5[d$period4ana == period2]
      if (length(v1) > 1 && length(v2) > 1 && sd(v1) > 0 && sd(v2) > 0) {
        ttest <- t.test(v1, v2)
        list(p_val = ttest$p.value,
             log2FC = log2(mean(v1, na.rm = TRUE) + 1e-6) - log2(mean(v2, na.rm = TRUE) + 1e-6))
      } else {
        list(p_val = NA, log2FC = NA)
      }
    })) %>%
    mutate(
      p_val = map_dbl(stats, "p_val"),
      log2FC = map_dbl(stats, "log2FC"),
      neg_log10_p = -log10(p_val),
      Significant = ifelse(!is.na(p_val) & p_val < 0.05 & abs(log2FC) > 0.2, TRUE, FALSE)
    )
  
  # Plot
  p <- ggplot(results, aes(x = log2FC, y = neg_log10_p, color = Significant)) +
    geom_point(alpha = 0.7, size = 1.1) +
    scale_color_manual(values = c("grey", "red")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "red") +
    geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", colour = "black") +
    labs(
      title = paste("Volcano Plot:", period1, "vs", period2),
      x = "log2 Fold Change (β-value)",
      y = "-log10(p-value)"
    ) +
    theme_minimal()
  
  print(p)
  
}



sum(table(ITdatamerge$V4))



(2581003/3222209)*100









# How many individuals have at least 1 coverage in each CPG for summary 
# CGPS on y axis 
# x is number of individuals 
# 1, 2, 3, 4 
# instead of bar plot do lines for different coverage values. 
# how many cpgs have atleast 
#plot of individuals against atleast one methylation at an amount of cpgs


coveragedf <- ITdatamerge



#---------------------------- The first plot 
hist(coveragedf$V4[coveragedf$V4 >= 4], ylab = "Number of CpGs", xlab = "Coverage", 
     main = "CpG count by Coverage")
#----------------------------  





# Y is number of cpgs, X is number of individuals

plotset <- ITdatamerge
plotset$V7 <- plotset$V3 >= 1
View(plotset)
# how many individuals have at least 1 methylation at CPG cites 
# you need the individual, also how many are true for this person 


newtable <- c()


cpg_set <- as.data.frame(table(plotset$V6))
cpg_set_methy <- as.data.frame(table(plotset$V6[ plotset$V7 == TRUE  ] ))
cpg_set$TrueMethyl <- cpg_set_methy$Freq
cpg_set

# match cpg_set rows to cpg_set_methy by CpG ID
idx <- match(cpg_set$Var1 , cpg_set_methy$Var1 )

cpg_set$TrueMethyl <- cpg_set_methy$Freq[idx]
cpg_set$TrueMethyl[is.na(cpg_set$TrueMethyl)] <- 0L
cpg_set





#--------------------------------------------- The second plot 

hist(cpg_set$TrueMethyl[cpg_set$TrueMethyl >= 10])


plot(cpg_set$Freq, cpg_set$TrueMethyl, main = , xlab = "" )

#--------------------------------------------




#i need it to loop through the plotset table, and count how many times each CPG occurs 
#count how many times it has a methylation 

 

# also make one for coverage also 
# essentially combining the plots 


# Run the it data where the 
# one plot for matching hvcpg and another for all cpgs 
# plot for the array CPGs 
# If there are enough cpgs in ppls u can run pca
 
# Dont just look coverage also the proportion of methylation across individuals. 
# The amount of methylation data was v low 
# lack of variability at hvcpgs (cus they are 0) 

# bar plot 0 1, or how many people have a methylation of 1 
# how many cpgs have at least one methylated 

# how many people have atleast 1 read 
# once for every cpg make a histogram 
# make a loop  to check for how many people have atleast 1 in each cpg ! 

ITdatamerge
# The number of cpg which have a coverage about a certain number in a number of individuals 
library(ggplot2)
trail <- ITdatamerge
trailwitsource <- as.data.frame(table(trail$source))
trail$V5 <- NULL
 

ggplot(trailwitsource, aes(x = reorder(Var1, -Freq), y = Freq)) +
  geom_col() +
  labs(x = NULL, y = "Count", title = "Counts by Var1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ITdatamerge$V5
View(ITdatamerge)
# how many ind have atleast coverage 1 at 4 cpgs 
# prep for pres !!! 

#x num of ind 
#y num of cpg 

summary(lm(data = cpg_set, TrueMethyl ~ Freq))

# Depth of data, how many cpgs are common in more nubmer of people ** # Do the same with CpGs with at least one methylation 

names(trail) <- c("chr", "pos", "mthyl", "cov", "ID")

head(trail) 



trail <- trail %>%
  mutate(SiteID = paste0(chr, "_", pos))


 
write.csv(trail, file = "trail.csv")


library(paletteer)

 


#----------------------------------------------------------------- playground 


library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(scales)
 
# Build cumulative "at least N people" distribution for a given min coverage
cum_dist <- function(dat, min_cov) {
  site_counts <- dat %>%
    filter(!is.na(cov), cov >= min_cov) %>%           # cumulative: ≥ min_cov
    distinct(chr, pos, ID) %>%                        # one row per site/person
    count(chr, pos, name = "n_people")                # people per site
  
  dist <- site_counts %>%
    count(n_people, name = "n_sites") %>%             # sites with exactly N people
    arrange(n_people) %>%
    complete(n_people = full_seq(n_people, 1), fill = list(n_sites = 0)) %>%
    mutate(
      n_sites_at_least = rev(cumsum(rev(n_sites))),   # sites with ≥ N people
      n_people = as.integer(n_people),
      min_cov  = min_cov,
      label    = paste0("≥", min_cov)
    )
  
  dist %>% select(n_people, n_sites_at_least, min_cov, label)
}

# Compute for thresholds 4..10
dist_all_cum <- map_df(4:10, ~cum_dist(trail, .x))

# make sure 'label' has ordered factor levels
dist_all_cum <- dist_all_cum %>%
  mutate(label = factor(label, levels = paste0("≥", 4:10)))

#dist_all_cum$min_cov <- as.factor(dist_all_cum$min_cov)
# Colours per threshold (you can swap these to your palette of choice)
Hcolors <- c("≥4"="grey90", "≥5"="#1b9e77","≥6"="#d95f02","≥7"="#7570b3",
             "≥8"="#e7298a","≥9"="#66a61e","≥10"="#e6ab02")

pdf("C:/Users/hamda/Desktop/UGI_thingy/fig1.pdf", width = 10, height = 5)
ggplot() +
  # overlays for ≥5…≥10
  geom_col(
    data = dist_all_cum,
    aes(x = n_people, y = n_sites_at_least, fill = label),
    position = "identity", width = 0.9, alpha = 0.65, colour = "black"
  ) +
  scale_fill_manual(values = Hcolors, name = "Min coverage") +
  scale_y_continuous(
    trans = "log10",
    labels = label_number(scale_cut = cut_short_scale()),
    breaks = 10^(0:6),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_x_continuous(breaks = pretty_breaks(n = 12),
                     expand = expansion(mult = c(0.01, 0.01))) +
  labs(
    title = "Sites shared by at least N people at different coverage thresholds",
    x = "N people (threshold)",
    y = "CpG sites (log scale)"
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title.position = "plot",
        legend.position = "right")
dev.off()
# how much data tdo we have # can we compare variability!!! 
# drops of quite fats ! # treshold - enough people enough sites with enough coveraage 
# exact same with is it methylated, do we have methylation signal 



# Only count CpG sites that are methylated, then compute "≥ N people" per site
cum_dist_methylated <- function(dat, min_cov) {
  site_counts <- dat %>%
    filter(!is.na(cov), cov >= min_cov,
           !is.na(mthyl), mthyl > 0) %>%         # keep only methylated CpGs
    distinct(chr, pos, ID) %>%                   # one row per site/person
    count(chr, pos, name = "n_people")           # how many people have this site methylated
  
  dist <- site_counts %>%
    count(n_people, name = "n_sites") %>%        # sites with exactly N people methylated
    arrange(n_people) %>%
    complete(n_people = full_seq(n_people, 1), fill = list(n_sites = 0)) %>%
    mutate(
      n_sites_at_least = rev(cumsum(rev(n_sites))),  # sites with ≥ N people methylated
      n_people = as.integer(n_people),
      min_cov  = min_cov,
      label    = paste0("≥", min_cov)
    )
  
  dist %>% select(n_people, n_sites_at_least, min_cov, label)
}

# Compute for thresholds 4..10 using your 'trail' table
dist_all_cum <- map_df(4:10, ~cum_dist_methylated(trail, .x)) %>%
  mutate(label = factor(label, levels = paste0("≥", 4:10)))

# Palette (same as before)
Hcolors <- c("≥4"="grey90", "≥5"="#1b9e77","≥6"="#d95f02","≥7"="#7570b3",
             "≥8"="#e7298a","≥9"="#66a61e","≥10"="#e6ab02")


pdf("C:/Users/hamda/Desktop/UGI_thingy/fig2.pdf", width = 10, height = 5)
ggplot() +
  geom_col(
    data = dist_all_cum,
    aes(x = n_people, y = n_sites_at_least, fill = label),
    position = "identity", width = 0.9, alpha = 0.65, colour = "black"
  ) +
  scale_fill_manual(values = Hcolors, name = "Min coverage") +
  scale_y_continuous(
    trans = "log10",
    labels = label_number(scale_cut = cut_short_scale()),
    breaks = 10^(0:6),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_x_continuous(breaks = pretty_breaks(n = 12),
                     expand = expansion(mult = c(0.01, 0.01))) +
  labs(
    title = "Methylated CpG sites shared by at least N people (by coverage threshold)",
    x = "N people",
    y = "Methylated CpG sites (log scale)"
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title.position = "plot",
        legend.position = "right")


dev.off()


# --- setup ---------------------------------------------------------------
library(dplyr)
library(tidyr)
library(stringr)

# Ensure numeric and build a robust methylation proportion column
# Uses V5 if present; otherwise falls back to V3/V4 when available.
testmerge <- testmerge %>%
  mutate(
    V3   = suppressWarnings(as.numeric(V3)),
    V4   = suppressWarnings(as.numeric(V4)),
    V5   = suppressWarnings(as.numeric(V5)),
    meth = dplyr::coalesce(V5, ifelse(!is.na(V3) & !is.na(V4) & V4 > 0, V3 / V4, NA_real_))
  )

# --- 1) Keep CpGs with frequency > 10 -----------------------------------
cpg_counts <- testmerge %>% count(V1, name = "n_obs")
cpg_keep   <- cpg_counts %>% filter(n_obs > 10)

tm_sub <- testmerge %>%
  semi_join(cpg_keep, by = "V1") %>%
  filter(!is.na(meth))                      # keep rows with a valid methylation proportion

# (Optional: add any extra QC filters here, e.g. minimum coverage V4 >= 4)

# --- 2) Overall variance per CpG ----------------------------------------
cpg_var_overall <- tm_sub %>%
  group_by(V1) %>%
  summarise(
    n_obs       = n(),
    mean_meth   = mean(meth),
    var_meth    = if (n() >= 2) var(meth) else NA_real_,
    .groups = "drop"
  ) %>%
  arrange(desc(var_meth))

# --- 3) Per-period variance per CpG (within-period) ---------------------
cpg_var_by_period_long <- tm_sub %>%
  group_by(V1, period4ana) %>%
  summarise(
    n_obs_period = n(),
    mean_meth    = mean(meth),
    var_meth     = if (n() >= 2) var(meth) else NA_real_,
    .groups = "drop"
  )

# Wide view: one row per CpG, columns are var per period
cpg_var_by_period_wide <- cpg_var_by_period_long %>%
  select(V1, period4ana, var_meth) %>%
  pivot_wider(
    names_from  = period4ana,
    values_from = var_meth,
    names_prefix = "var_"
  )

# --- 4) (Optional) Between-period variance of period means ---------------
# This captures how much the *mean methylation* differs across periods for each CpG.
cpg_between_period_var <- cpg_var_by_period_long %>%
  group_by(V1) %>%
  summarise(
    var_between_period_means =
      if (n() >= 2) var(mean_meth, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )

# --- 5) Final merged table ----------------------------------------------
cpg_variance_summary <- cpg_var_overall %>%
  left_join(cpg_var_by_period_wide, by = "V1") %>%
  left_join(cpg_between_period_var, by = "V1")

# Objects you now have:
# - cpg_var_overall: V1, n_obs, mean_meth, var_meth (overall)
# - cpg_var_by_period_long: long per-period stats per CpG
# - cpg_var_by_period_wide: wide table of per-period variances
# - cpg_variance_summary: overall + per-period + between-period variance in one DF

# Example: peek at the top high-variance CpGs
head(cpg_variance_summary, 20)
View(cpg_variance_summary)







# --- packages -------------------------------------------------------------
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

# --- prep: methylation proportion and CpG filtering ----------------------
testmerge <- testmerge %>%
  mutate(
    V3   = suppressWarnings(as.numeric(V3)),       # methylated reads
    V4   = suppressWarnings(as.numeric(V4)),       # total coverage
    V5   = suppressWarnings(as.numeric(V5)),       # precomputed proportion (if present)
    meth = dplyr::coalesce(V5, ifelse(!is.na(V3) & !is.na(V4) & V4 > 0, V3 / V4, NA_real_))
  )

# Keep CpGs seen > 10 times
cpg_keep <- testmerge %>%
  count(V1, name = "n_obs") %>%
  filter(n_obs > 10)

# Subset to valid rows (optionally enforce coverage threshold, e.g. V4 >= 4)
tm_sub <- testmerge %>%
  semi_join(cpg_keep, by = "V1") %>%
  filter(!is.na(meth), !is.na(period4ana), !is.na(V4), V4 > 0)

# --- overall + per-period variance (as before) ---------------------------
cpg_var_overall <- tm_sub %>%
  group_by(V1) %>%
  summarise(
    n_obs     = n(),
    mean_meth = mean(meth),
    var_meth  = if (n() >= 2) var(meth) else NA_real_,
    .groups = "drop"
  )

cpg_var_by_period_long <- tm_sub %>%
  group_by(V1, period4ana) %>%
  summarise(
    n_obs_period = n(),
    mean_meth_p  = mean(meth),
    var_meth_p   = if (n() >= 2) var(meth) else NA_real_,
    .groups = "drop"
  )

cpg_var_by_period_wide <- cpg_var_by_period_long %>%
  select(V1, period4ana, var_meth_p) %>%
  pivot_wider(names_from = period4ana, values_from = var_meth_p, names_prefix = "var_")

cpg_between_period_var <- cpg_var_by_period_long %>%
  group_by(V1) %>%
  summarise(
    var_between_period_means = if (n() >= 2) var(mean_meth_p, na.rm = TRUE) else NA_real_,
    n_periods = n_distinct(period4ana),
    .groups = "drop"
  )

cpg_variance_summary <- cpg_var_overall %>%
  left_join(cpg_var_by_period_wide, by = "V1") %>%
  left_join(cpg_between_period_var, by = "V1")

# --- per-CpG ANOVA tests --------------------------------------------------
run_tests_for_cpg <- function(df) {
  # df is rows for a single CpG (columns: meth, V3, V4, period4ana, etc.)
  out <- list(
    aov_p = NA_real_, aov_df1 = NA_integer_, aov_df2 = NA_integer_, aov_eta2 = NA_real_,
    glm_p = NA_real_, glm_dev_null = NA_real_, glm_dev_full = NA_real_, glm_r2_mcfadden = NA_real_
  )
  
  n_periods <- dplyr::n_distinct(df$period4ana)
  # ---------- weighted ANOVA on proportion ----------
  if (n_periods >= 2 && length(unique(df$meth)) > 1) {
    mod_aov <- try(aov(meth ~ period4ana, data = df, weights = V4), silent = TRUE)
    if (!inherits(mod_aov, "try-error")) {
      an <- summary(mod_aov)[[1]]
      if ("period4ana" %in% rownames(an)) {
        out$aov_p   <- an["period4ana", "Pr(>F)"]
        out$aov_df1 <- an["period4ana", "Df"]
        out$aov_df2 <- an["Residuals",  "Df"]
        ss_eff <- an["period4ana", "Sum Sq"]
        ss_err <- an["Residuals",  "Sum Sq"]
        out$aov_eta2 <- if (is.finite(ss_eff) && is.finite(ss_err) && (ss_eff + ss_err) > 0) {
          ss_eff / (ss_eff + ss_err)     # partial eta^2
        } else NA_real_
      }
    }
  }
  
  # ---------- binomial GLM (analysis of deviance) ----------
  dfg <- df %>% filter(!is.na(V3), !is.na(V4), V4 > 0, V3 >= 0, V3 <= V4)
  if (nrow(dfg) >= 2 && dplyr::n_distinct(dfg$period4ana) >= 2) {
    m0 <- try(glm(cbind(V3, V4 - V3) ~ 1,           family = binomial, data = dfg), silent = TRUE)
    m1 <- try(glm(cbind(V3, V4 - V3) ~ period4ana,  family = binomial, data = dfg), silent = TRUE)
    if (!inherits(m0, "try-error") && !inherits(m1, "try-error")) {
      lr <- anova(m0, m1, test = "Chisq")
      out$glm_p          <- lr$`Pr(>Chi)`[2]
      out$glm_dev_null   <- m0$deviance
      out$glm_dev_full   <- m1$deviance
      out$glm_r2_mcfadden <- 1 - (m1$deviance / m0$deviance)
    }
  }
  
  tibble::tibble(
    aov_p = out$aov_p, aov_df1 = out$aov_df1, aov_df2 = out$aov_df2, aov_eta2 = out$aov_eta2,
    glm_p = out$glm_p, glm_dev_null = out$glm_dev_null, glm_dev_full = out$glm_dev_full,
    glm_r2_mcfadden = out$glm_r2_mcfadden
  )
}

cpg_anova <- tm_sub %>%
  group_by(V1) %>%
  group_modify(~ run_tests_for_cpg(.x)) %>%
  ungroup() %>%
  mutate(
    aov_p_adj = p.adjust(aov_p, method = "BH"),
    glm_p_adj = p.adjust(glm_p, method = "BH")
  )

# --- final merged table: variances + ANOVA/GLM tests ----------------------
cpg_results <- cpg_variance_summary %>%
  left_join(cpg_anova, by = "V1") %>%
  arrange(dplyr::coalesce(glm_p, 1), dplyr::coalesce(aov_p, 1))

# Example: inspect top CpGs by GLM p-value
head(cpg_results, 20)


install.packages("leaflet")
library(leaflet)
leaflet() %>% addTiles() %>% setView(lng = -0.1276, lat = 51.5074, zoom = 12)
















# packages
install.packages(c("leaflet","dplyr","scales"))
library(dplyr)
library(scales)
library(leaflet)

# 1) keep one row per sample/coordinate (avoids stacked duplicates)
dat <- testmerge %>%
  filter(!is.na(Latitude), !is.na(Longitude)) %>%
  distinct(sample_alias, Latitude, Longitude, .keep_all = TRUE) %>%
  mutate(
    meth = as.numeric(meth),
    # circle radius scaled from small to big
    radius = rescale(meth, to = c(4, 18), from = range(meth, na.rm = TRUE)),
    colour_val = Latitude     # <-- change to Longitude if you prefer
  )

# 2) colour palette (by latitude; change title if you switch to Longitude)
pal <- colorNumeric(palette = "viridis", domain = dat$colour_val, na.color = "#BDBDBD")

# 3) map
leaflet(dat) %>%
  addProviderTiles("CartoDB.Positron") %>%
   

addCircleMarkers(
  lng = ~Longitude, lat = ~Latitude,
  radius = ~radius,
  fillColor = ~pal(colour_val),
  fillOpacity = 0.75,   # ← was 0.85; lower = more see-through
  stroke = TRUE, color = "#333333",
  opacity = 0.5,        # ← border transparency (optional)
  weight = 1,
  label = ~paste0(sample_alias, " — ", Site),
  popup = ~paste0("<b>", sample_alias, "</b><br/>", "Site: ", Site, "<br/>",
                  "Period: ", period4ana, "<br/>",
                  "Meth: ", scales::percent(meth, accuracy = 0.1), "<br/>",
                  "Lat/Lon: ", round(Latitude, 5), ", ", round(Longitude, 5))
)




#0.558, 0.641, 0.631



ITdatamerge$V5 <- ITdatamerge$V3 / ITdatamerge$V4

library(dplyr)

 

# (Optional) normalise source names like "ITdata1T" -> "ITdata1"
ITdatamerge2 <- ITdatamerge %>%
  mutate(source = sub("^([A-Za-z]+\\d+).*", "\\1", source))

mean_by_source <- ITdatamerge2 %>%
  group_by(source) %>%
  summarise(
    n_rows       = n(),
    n_with_V4    = sum(!is.na(V4)),
    mean_V4      = mean(V4, na.rm = TRUE)  # mean methylation per person
  ) %>%
  arrange(as.integer(gsub("\\D", "", source)))

mean_by_source
# If you want it rounded:
# mean_by_source %>% mutate(mean_V3 = round(mean_V3, 3))

# (Optional) save to CSV for PowerPoint/Excel
# write.csv(mean_by_source, "mean_methylation_by_source.csv", row.names = FALSE)



library(dplyr)

# (Optional) normalise source names like "ITdata1T" -> "ITdata1"
ITdatamerge2 <- ITdatamerge %>%
  mutate(source = sub("^([A-Za-z]+\\d+).*", "\\1", source))

by_source <- ITdatamerge2 %>%
  group_by(source) %>%
  summarise(
    n_rows        = n(),
    n_with_V5     = sum(!is.na(V5)),
    mean_V5       = mean(V5, na.rm = TRUE),  # unweighted mean of site-level methylation proportions
    sd_V5         = sd(V5, na.rm = TRUE),
    # coverage-weighted methylation proportion across all sites for this person
    weighted_mean = if (sum(V4, na.rm = TRUE) > 0) sum(V3, na.rm = TRUE) / sum(V4, na.rm = TRUE) else NA_real_
  ) %>%
  arrange(as.integer(gsub("\\D", "", source)))

by_source
# Optional: rounding / export
# by_source %>% mutate(across(c(mean_V5, weighted_mean), ~round(.x, 3)))
# write.csv(by_source, "mean_methylation_by_source_V5.csv", row.names = FALSE)

values <- c(0.558, 0.641, 0.631)

library(dplyr)
library(ggplot2)
library(scales)

# Ensure numeric order for individuals
plot_df <- by_source %>%
  mutate(source_num = as.integer(gsub("\\D", "", source))) %>%
  arrange(source_num)

# Reference "modern human methylation rate" lines
ref <- data.frame(value = c(0.558, 0.641, 0.631))
ref_breaks <- sort(unique(ref$value))
ref_labels <- paste0(formatC(100 * ref_breaks, format = "f", digits = 1), "%")

# Y-axis upper limit to include both your means and the reference lines
ymax <- max(c(plot_df$mean_V5, ref$value), na.rm = TRUE) * 1.05

# Reasonable x breaks/labels (every 5th individual)
x_breaks <- seq(min(plot_df$source_num), max(plot_df$source_num), by = 5)
x_labels <- paste0("ITdata", x_breaks)

pdf("C:/Users/hamda/Desktop/UGI_thingy/fig5.pdf", width = 10, height = 5)

ggplot(plot_df, aes(x = source_num, y = mean_V5)) +
  geom_point(size = 1.8, alpha = 0.9) +
  geom_hline(
    data = ref,
    aes(yintercept = value, colour = factor(value)),
    linetype = "dashed", linewidth = 0.7
  ) +
  scale_colour_discrete(
    name = "Modern human methylation rate",
    breaks = ref_breaks, labels = ref_labels
  ) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), limits = c(0, ymax)) +
  labs(
    x = "Individual (source)",
    y = "Mean methylation (V5)",
    title = "Mean site-level methylation per individual with modern human reference rates"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

# Optional: save
# ggsave("mean_methylation_scatter_with_refs.png", width = 10, height = 5, dpi = 300)

dev.off()

