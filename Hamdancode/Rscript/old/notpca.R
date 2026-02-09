
# Prepare data: remove rows with too many NAs, impute, then transpose
install.packages("FactoMineR")
install.packages("factoextra")
library(FactoMineR)
library(factoextra)

# Impute simple row-wise mean
methyl_imputed <- apply(methyl_matrix, 1, function(row) {
  row[is.na(row)] <- mean(row, na.rm = TRUE)
  return(row)
})
methyl_imputed <- t(methyl_imputed)

# Transpose: samples as rows
pca_input <- t(methyl_imputed)

# Run PCA
res.pca <- PCA(pca_input, graph = FALSE)

# Plot PCA with group annotation
fviz_pca_ind(res.pca, 
             habillage = population,
             addEllipses = TRUE,
             title = "PCA of Methylation Profiles")


View(head(site_var))





library(pheatmap)
var_matrix
# Step 1: Get variance across time periods for each CpG
row_var <- apply(var_matrix, 1, var)

# Step 2: Select top 100 CpGs with most variance
top_cpgs <- names(sort(row_var, decreasing = TRUE)[1:50])

# Step 3: Subset matrix
heatmap_data <- var_matrix[top_cpgs, ]

# Step 4: Z-score scaling (optional)
heatmap_scaled <- t(scale(t(heatmap_data)))

# Step 5: Plot heatmap
pheatmap(heatmap_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Top 100 CpGs with Highest Variance Across Time Periods")

# If using 450K array IDs (adjust for your actual probe IDs)
# If you have genomic coordinates instead, use annotatr

# For example using annotatr with hg19 regions
BiocManager::install("annotatr")
BiocManager::install("GRanges")
library(annotatr)
library(GRanges)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
library(GenomicRanges)




# Build GRanges object
gr <- stringr::str_split_fixed(rownames(var_matrix), "_", 2)
gr <- GRanges(seqnames = paste0("chr", gr[,1]),
              ranges = IRanges(start = as.integer(gr[,2]), width = 1))

# Annotate using CpG features
annotations <- build_annotations(genome = 'hg19', annotations = 'hg19_cpgs')
annotated <- annotate_regions(regions = gr, annotations = annotations, ignore.strand = TRUE)

View(as.data.frame(annotated))







site_var <- site_var %>%
  mutate(SiteID = paste0(Chr, "_", Pos))


var_matrix <- site_var %>%
  select(SiteID, period4ana, var_site) %>%
  pivot_wider(names_from = period4ana, values_from = var_site) %>%
  column_to_rownames("SiteID") %>%
  as.matrix()


var_matrix[is.na(var_matrix)] <- rowMeans(var_matrix, na.rm = TRUE)[row(var_matrix)[is.na(var_matrix)]]


pca_input <- t(var_matrix)


res.pca <- PCA(pca_input, graph = FALSE)


group_metadata <- data.frame(
  period = rownames(pca_input)
)


time_palette <- c(
  "Copper Age" = "#E69F00",
  "Imperial" = "#56B4E9",
  "Iron/Republic" = "#009E73",
  "LateAntiquity" = "#F0E442",
  "Medieval/EarlyModern" = "#0072B2",
  "Mesolithic" = "#D55E00",
  "Neolithic" = "#CC79A7",
  "(not included in analyses)" = "grey60"
)


fviz_pca_ind(res.pca,
             geom.ind = "point",
             pointshape = 21,
             pointsize = 6,
             fill.ind = group_metadata$period,
             col.ind = "black",
             palette = time_palette,
             repel = TRUE,
             addEllipses = FALSE,
             title = "PCA of CpG Methylation Variance by Time Period") +
  theme_minimal(base_size = 14) +
  guides(fill = guide_legend(title = "Time Period"))





library(dplyr)
library(tidyr)
library(FactoMineR)
library(factoextra)
library(ggplot2)

# Step 1: Build CpG site ID
site_var <- site_var %>%
  mutate(SiteID = paste0(Chr, "_", Pos))

# Step 2: Pivot wider
site_var_wide <- site_var %>%
  select(SiteID, period4ana, var_site) %>%
  pivot_wider(names_from = period4ana, values_from = var_site)

# Step 3: Build matrix
var_matrix <- site_var_wide %>%
  column_to_rownames("SiteID") %>%
  as.matrix()

# Step 4: Impute missing values with row mean
var_matrix[is.na(var_matrix)] <- rowMeans(var_matrix, na.rm = TRUE)[row(var_matrix)[is.na(var_matrix)]]

# Step 5: Transpose — PCA on time periods
pca_input <- t(var_matrix)

# Step 6: Run PCA
res.pca <- PCA(pca_input, graph = FALSE)

# Step 7: Build metadata for aesthetic labels
groups_df <- data.frame(
  period = rownames(pca_input)
)

# Step 8: Styled PCA Plot (matching earlier style)
fviz_pca_ind(res.pca,
             geom.ind = "point",
             pointshape = 21,
             pointsize = 5,
             fill.ind = "tomato",
             col.ind = "black",
             repel = TRUE,
             title = "PCA of CpG Site Variance by Time Period") +
  theme_minimal(base_size = 14)

group_colors <- c(
  "Copper Age" = "#E69F00",
  "Imperial" = "#56B4E9",
  "Iron/Republic" = "#009E73",
  "LateAntiquity" = "#F0E442",
  "Medieval/EarlyModern" = "#0072B2",
  "Mesolithic" = "#D55E00",
  "Neolithic" = "#CC79A7",
  "(not included in analyses)" = "grey70"
)

fviz_pca_ind(res.pca,
             geom.ind = "point",
             pointshape = 21,
             pointsize = 6,
             fill.ind = groups_df$period,
             col.ind = "black",
             palette = group_colors,
             addEllipses = TRUE,
             label = "none",
             repel = TRUE,
             title = "PCA of CpG Site Variance by Time Period") +
  theme_minimal(base_size = 14) +
  guides(fill = guide_legend(title = "Time Period"))




# Get top contributing variables to PC1 and PC2
var_contrib <- res.pca$var$contrib

# Top 10 for Dim1
top_pc1 <- sort(var_contrib[, "Dim.1"], decreasing = TRUE)[1:10]
print("Top CpGs driving PC1:")
print(top_pc1)

# Top 10 for Dim2
top_pc2 <- sort(var_contrib[, "Dim.2"], decreasing = TRUE)[1:10]
print("Top CpGs driving PC2:")
print(top_pc2)

#fviz_pca_var(res.pca,
#col.var = "contrib",
# gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
# repel = TRUE,
# title = "CpG Contributions to PCA Dimensions")


fviz_pca_var(res.pca,
             select.var = list(contrib = 50),  # show only top 50 CpGs
             col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,
             title = "Top 50 CpG Contributions to PCA")




Step 1: Filter to top 1000 most variable sites
site_sd <- apply(methyl_matrix, 1, sd, na.rm = TRUE)
top_sites <- names(sort(site_sd, decreasing = TRUE))[1:1000]
methyl_topvar <- methyl_matrix[top_sites, ]

# Step 2: Transpose for PCA (samples as rows)
pca_input <- t(methyl_topvar)

# Step 3: Run PCA

library(FactoMineR)
library(factoextra)
res.pca <- PCA(pca_input, graph = FALSE)
res.pca <- FactoMineR::PCA(pca_input, graph = FALSE)
# Step 4: Prepare sample metadata
sample_metadata <- testmerge %>%
  select(run_accession, period4ana) %>%
  distinct() %>%
  filter(run_accession %in% rownames(pca_input))

# Step 5: Plot PCA
fviz_pca_ind(res.pca,
             geom.ind = "point",
             pointshape = 21,
             pointsize = 5,
             fill.ind = sample_metadata$period4ana,
             col.ind = "black",
             palette = "jco",
             repel = TRUE,
             title = "PCA on Top 1000 Most Variable CpGs") +
  theme_minimal(base_size = 14)


# Pick outlier samples from PCA manually if needed
outliers <- c("ERR3546000", "ERR3559881")  # Replace with actual run IDs

# Extract methylation values
outlier_profiles <- methyl_matrix[, colnames(methyl_matrix) %in% outliers]
other_profiles <- methyl_matrix[, !colnames(methyl_matrix) %in% outliers]

# Plot density distributions
outlier_df <- data.frame(value = as.vector(outlier_profiles), type = "Outlier")
other_df   <- data.frame(value = as.vector(other_profiles), type = "Other")
combined_df <- rbind(outlier_df, other_df)

ggplot(combined_df, aes(x = value, fill = type)) +
  geom_density(alpha = 0.4) +
  labs(title = "Methylation Value Distributions: Outliers vs Others",
       x = "Methylation β-value", y = "Density") +
  theme_minimal()

View(head(testmerge2))

# Required libraries
library(tidyverse)

# Define periods
period_pair <- c("Imperial", "Iron/Republic")

# Perform t-test with filtering
meth_test <- testmerge2 %>%
  filter(period4ana %in% period_pair) %>%
  group_by(SiteID) %>%
  nest() %>%
  mutate(
    stats = map(data, ~ {
      d <- .
      v1 <- d$V5[d$period4ana == period_pair[1]]
      v2 <- d$V5[d$period4ana == period_pair[2]]
      if (length(unique(v1)) > 1 && length(unique(v2)) > 1) {
        list(p_val = t.test(v1, v2)$p.value,
             mean_diff = mean(v1, na.rm = TRUE) - mean(v2, na.rm = TRUE))
      } else {
        list(p_val = NA, mean_diff = NA)
      }
    }),
    p_val = map_dbl(stats, "p_val"),
    mean_diff = map_dbl(stats, "mean_diff")
  ) %>%
  ungroup() %>%
  filter(!is.na(p_val)) %>%
  mutate(logp = -log10(p_val))

# Volcano plot
ggplot(meth_test, aes(x = mean_diff, y = logp)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "red") +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", colour = "grey") +
  labs(
    x = "Mean β-Value Difference",
    y = "-log10(p-value)",
    title = paste("Volcano Plot: Methylation Level (", period_pair[1], "vs", period_pair[2], ")")
  ) +
  theme_minimal()


