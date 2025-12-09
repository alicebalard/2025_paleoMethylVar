
# Package intsallation
#-------------------------------------------------------------------------------    

install.packages("magrittr")
install.packages("cli")
install.packages("glue")   # to fix the permission denied one, after restart
install.packages("openxlsx")
install.packages("remotes")
remotes::install_github("al2na/methylKit")


# libraries 
#-------------------------------------------------------------------------------     

library(ggplot2)
library(methylKit)
library(openxlsx)
  
#-------------------------------------------------------------------------------     

# Loading Bisulfite sequencing (methylation reads)
#-------------------------------------------------------------------------------     

Mdata <- read.xlsx("C:/Users/hamda/Desktop/UGI_thingy/data/gkac503_supplemental_files/SUPPLEMENTARY_TABLES_resubmission.xlsx", sheet = 6)

file1 <- "C:/Users/hamda/Desktop/UGI_thingy/3rdyear/genome/quality/methylation_hg38/Zvej16_bisulfite.trim_bismark_bt2.deduplicated.bismark.cov.gz"

Zvej16 <- methRead(
  location   = file1,
  sample.id  = "Zvej16",
  assembly   = "hg38",
  pipeline   = "bismarkCoverage",
  mincov     = 0,
  treatment  = 0      
)

View(Zvej16)


file2 <- "C:/Users/hamda/Desktop/UGI_thingy/3rdyear/genome/quality/methylation_hg38/SP75_bisulfite.trim_bismark_bt2.deduplicated.bismark.cov.gz"

SP75 <- methRead(
  location   = file2,
  sample.id  = "SP75",
  assembly   = "hg38",
  pipeline   = "bismarkCoverage",
  mincov     = 0,
  treatment  = 0      
)

SP75$methy <- SP75$numCs / SP75$coverage
Zvej16$methy <- Zvej16$numCs / Zvej16$coverage

file1 <- "C:/Users/hamda/Desktop/UGI_thingy/3rdyear/genome/quality/methylation_hg38/Zvej16_bisulfite.trim_bismark_bt2.deduplicated.bismark.cov.gz"

write.csv(Zvej16, file = "Zvej16file.csv")
write.csv(SP75, file = "SP75file.csv")


#-------------------------------------------------------------------------------     

# Exploratory analysis
#-------------------------------------------------------------------------------     

# density of methylation (for positions)

ggplot(Zvej16, aes(x=methy)) + 
  geom_density()

ggplot(SP75, aes(x=methy)) + 
  geom_density()


# histogram of methylation (for positions)
 
ggplot(Zvej16, aes(x=coverage)) + 
  geom_histogram()

ggplot(SP75, aes(x=methy)) + 
  geom_histogram()


#-------------------------------------------------------------------------------     

 


