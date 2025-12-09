#This page is for the ANOVA calculations for the positions 
#In order to complete this 


ITdatamerge
testmerge


library(dplyr)
library(ggplot2)
library(forcats)

 

#–– 1A. Compute the variance of V5 for every Chr × Period
var_chr_period <- testmerge %>% 
  group_by(Period = period4ana, Chr = V1) %>% 
  summarise(var_beta   = var(V5, na.rm = TRUE),  # variance of methylation
            n_datapts  = n(),                    # rows that went into it
            .groups    = "drop")

#–– 1B. Quick diagnostic glance
table(var_chr_period_new$Chr)
View(testmerge)
head(var_chr_period)
var_chr_period_new <- subset(var_chr_period, n_datapts >= 5)  
var_chr_period_new <- var_chr_period_new[order(-var_chr_period_new$var_beta), ]
head(var_chr_period_new)

###
 
ggplot(var_chr_period_new, aes(x = factor(Chr),
                           y = var_beta,
                           fill = Period)) +
  geom_col(position = "dodge") +
  labs(x = "Chromosome", y = "Variance of methylation (V5)") +
  theme_minimal()



ggplot(var_chr_period_new, aes(x = factor(Chr), y = var_beta)) +
  geom_col() +
  facet_wrap(~Period, scales = "free_x") +
  labs(x = "Chromosome", y = "Variance of methylation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size = 6))













#### testing
library(dplyr)
library(ggplot2)
library(tidyr)
var_chr_period_new <- var_chr_period_new %>%
  separate(Chr, into = c("Chromosome", "Position"), sep = "_") %>%
  mutate(
    Chromosome = as.numeric(Chromosome),
    Position = as.numeric(Position),
    Chromosome = factor(Chromosome, levels = 1:22)
  )

ggplot(var_chr_period_new, aes(x = Position, y = var_beta)) +
  geom_point(size = 0.4, alpha = 0.6) +
  facet_wrap(~Chromosome, scales = "free_x") +
  labs(x = "Position on Chromosome", y = "Variance of methylation") +
  theme_minimal()
 












####test2 
head(var_chr_period_new)

top_cut <- 0.99  # top 1 %
df <- df %>% 
  mutate(top = var_beta  > quantile(var_beta , top_cut))

ggplot(var_chr_period_new, aes(pos/1e6, var_beta )) +
  geom_point(aes(colour = top), alpha = 0.6, size = 1) +
  scale_colour_manual(values = c("grey60", "red"), guide = "none") +
  facet_wrap(~ Chromosome , scales = "free_x", ncol = 4) +
  scale_x_continuous(labels = scales::label_number(suffix = " Mb")) +
  labs(x = "Genomic position (Mb)", y = "Variance of methylation (β)") +
  theme_bw(base_size = 12) +
  theme(panel.spacing = unit(0.8, "lines"))

library(dplyr)

var_chr_period_new <- var_chr_period_new %>% 
  rename(
    pos  = Position,   # lowercase name the code expects
    var  = var_beta    # only needed if you want the shorter name
  ) %>% 
  mutate(top = var > quantile(var, 0.99))

library(ggplot2); library(scales)

ggplot(var_chr_period_new, aes(pos/1e6, var)) +
  geom_point(aes(colour = top), alpha = 0.6, size = 1) +
  scale_colour_manual(values = c("grey60", "red"), guide = "none") +
  facet_wrap(~ Chromosome, scales = "free_x", ncol = 4) +
  scale_x_continuous(labels = label_number(suffix = " Mb")) +
  labs(x = "Genomic position (Mb)", y = "Variance of methylation (β)") +
  theme_bw(base_size = 12) +
  theme(panel.spacing = unit(0.8, "lines"))

library(dplyr); library(ggplot2); library(scales)

top_cut <- 0.99
q99     <- quantile(var_chr_period_new$var_beta, top_cut)

ggplot(df, aes(Position/1e6, var_beta)) +                 # <- real names
  geom_point(aes(colour = var_beta > q99), alpha = 0.6, size = 1) +
  scale_colour_manual(values = c("grey60", "red"), guide = "none") +
  facet_wrap(~ Chromosome, scales = "free_x", ncol = 4) +
  scale_x_continuous(labels = label_number(suffix = " Mb")) +
  labs(x = "Genomic position (Mb)", y = "Variance of methylation (β)") +
  theme_bw(base_size = 12) +
  theme(panel.spacing = unit(0.8, "lines"))

head(var_chr_period_new)
ggplot(var_chr_period_new,
       aes(Position/1e6, var_beta, colour = Period)) +
  geom_point(alpha = 0.7, size = 0.8) +
  facet_wrap(~ Chromosome, scales = "free_x", ncol = 4) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Genomic position (Mb)",
       y = "Variance of methylation (β)") +
  theme_bw()


ggplot(var_chr_period_new,
       aes(Position/1e6, var_beta)) +
  geom_count(alpha = 0.5) +                # size ∝ number of rows
  facet_wrap(~ Chromosome, scales = "free_x", ncol = 4) +
  scale_size(range = c(0.5, 3), guide = "none") +
  labs(x = "Genomic position (Mb)",
       y = "Variance of methylation (β)") +
  theme_bw()


var_chr_period_uni <- var_chr_period_new %>% 
  group_by(Period, Chromosome, Position) %>% 
  summarise(var_beta = mean(var_beta),
            n_datapts = sum(n_datapts),
            .groups   = "drop")


  
var_cpg <- var_chr_period_new %>% 
  group_by(Chromosome, Position) %>% 
  summarise(var_beta = mean(var_beta),   # or max(), median(), etc.
            .groups   = "drop")

ggplot(var_cpg, aes(Position/1e6, var_beta)) +
  geom_point(alpha = 0.6, size = 0.8) +
  facet_wrap(~ Chromosome, scales = "free_x", ncol = 4) +
  labs(x = "Genomic position (Mb)",
       y = "Variance of methylation (β)") +
  theme_bw()

df_plot <- var_chr_period_uni      # choose the one you want
top_cut <- 0.99
q99     <- quantile(df_plot$var_beta, top_cut)

ggplot(df_plot,
       aes(Position/1e6, var_beta,
           colour = var_beta > q99)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_colour_manual(values = c("grey60", "red"), guide = "none") +
  facet_wrap(~ Chromosome, scales = "free_x", ncol = 4) +
  scale_x_continuous(labels = scales::label_number(suffix = " Mb")) +
  labs(x = "Genomic position (Mb)",
       y = "Variance of methylation (β)") +
  theme_bw(base_size = 12) +
  theme(panel.spacing = unit(0.8, "lines"))


###

imperial


ggplot(var_chr_period_new,
       aes(fct_reorder(factor(Chr), var_beta),
           var_beta)) +
  geom_point(size = 2, colour = "#69b3a2") +
  facet_wrap(~ Period, nrow = 1) +
  labs(x = "Chromosome", y = "Variance of methylation") +
  theme_minimal()

head(var_chr_period_new)
##


ggplot(var_chr_period_new,
       aes(x = factor(Chr),          # each chromosome gets its own slot
           y = var_beta)) +
  geom_point(size = 2, colour = "#2ca25f") +
  facet_wrap(~ Period, nrow = 1) +   # keeps your period facets
  labs(x = "Chromosome",
       y = "Variance of methylation") +
  theme_minimal()


##


period_rank <- var_chr_period_new %>%
  group_by(Period) %>%
  summarise(med = median(var_beta)) %>%
  arrange(desc(med)) %>%
  pull(Period)

var_chr_period_new <- var_chr_period_new %>%
  mutate(Period = factor(Period, levels = period_rank))
period_rank

##


ggplot(var_chr_period_new,
       aes(factor(Chr), var_beta)) +
  geom_point(alpha = 0.6, colour = "#2ca25f") +
  stat_summary(fun = median,
               geom = "segment",
               aes(xend = after_stat(x),
                   yend = after_stat(y)),
               linewidth = 2,
               colour = "black") +
  facet_wrap(~ Period, nrow = 1)
head(testmerge)


###################################################################################
##############################################################################
# 1. Split the “chr_position” column (V1) into Chr and Pos -------------------
#    convert = TRUE turns Chr and Pos into numeric; drop it if you want Chr
#    to stay as a character (e.g. to keep “X”, “Y”).
testmerge2 <- testmerge %>% 
  separate(V1, into = c("Chr", "Pos"), sep = "_", convert = TRUE)

##############################################################################
# 2. OPTIONAL: remove duplicate calls of the same CpG in the same individual -
#    Keeps just one row per sample × CpG × period.
testmerge2 <- testmerge2 %>%
  distinct(sample_alias, Chr, Pos, period4ana, .keep_all = TRUE)

##############################################################################
# 3.  VARIANCE ACROSS CpG SITES *WITHIN* EACH CHR × PERIOD -------------------
#    (i.e. treat CpGs as the replicates)
var_chr_period <- testmerge2 %>% 
  group_by(period4ana, Chr) %>% 
  summarise(var_beta = var(V5, na.rm = TRUE),
            n_sites  = n(),              # how many CpGs went into it
            .groups  = "drop")

##############################################################################
# 4.  QUICK SCATTER PLOT (facets = periods) ----------------------------------
ggplot(var_chr_period,
       aes(x = factor(Chr), y = var_beta)) +
  geom_point(size = 2, colour = "#2ca25f") +
  facet_wrap(~ period4ana, nrow = 1) +
  labs(x = "Chromosome",
       y = "Variance of methylation (V5)") +
  theme_minimal()

##############################################################################
# 5.  (ALTERNATIVE) VARIANCE ACROSS INDIVIDUALS *AT EACH SITE* ---------------
#     then summarise those site-level variances per Chr × Period
site_var <- testmerge2 %>% 
  group_by(Chr, Pos, period4ana) %>% 
  filter(n() > 1) %>%                    # <-- keep sites seen in ≥2 individuals
  summarise(var_site = var(V5, na.rm = TRUE), .groups = "drop")


var_chr_period_site <- site_var %>% 
  group_by(period4ana, Chr) %>% 
  summarise(mean_site_var = mean(var_site, na.rm = TRUE),
            med_site_var  = median(var_site, na.rm = TRUE),
            n_sites       = n(),
            .groups       = "drop")

# Plot median of site-level variances
ggplot(var_chr_period_site,
       aes(x = factor(Chr), y = med_site_var)) +
  geom_col(fill = "#69b3a2") +
  facet_wrap(~ period4ana, nrow = 1) +
  labs(x = "Chromosome",
       y = "Median per-site variance of methylation") +
  theme_minimal()
