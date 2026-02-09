install.packages("openxlsx")
library(openxlsx)

ITdatamerge

Allfiles <- lapply(files, read.delim)
Allfiles <- lapply(Allfiles, shorterningfunc)

Allfiles[[130]]
names(Allfiles)

excelLOC1 <- read.xlsx("../UGI_thingy/data/aay6826_tables_s1_to_s4.xlsx", sheet = 1)
excelLOC2 <- read.xlsx("../UGI_thingy/data/aay6826_tables_s1_to_s4.xlsx", sheet = 2)
excelLOC3 <- read.xlsx("../UGI_thingy/data/aay6826_tables_s1_to_s4.xlsx", sheet = 3)
excelLOC4 <- read.xlsx("../UGI_thingy/data/aay6826_tables_s1_to_s4.xlsx", sheet = 4)
excelLOC5 <- read.xlsx("../UGI_thingy/data/aay6826_tables_s1_to_s4.xlsx", sheet = 5)
filetemp <- read.delim("data/filereport_read_run_PRJEB32566_tsv.txt")
filetemp
head(excelLOC1)
nrow(excelLOC1)
 
# 1) List all your files
files <- list.files(
  "../UGI_thingy/data/ALL_MethylMaps_Antonio2019",
  full.names = TRUE
)

# 2) Prepare a list to hold each data.frame
sample_list <- list()

for (file_path in files) {
  # 3) Skip empty files
  if (file.info(file_path)$size == 0) {
    message("Skipping empty file: ", file_path)
    next
  }
  
  # 4) Derive the sample ID == first 10 chars of the filename
  fn <- basename(file_path)
  samp <- substr(fn, 1, 10)
  
  # 5) Read the data (skipping header line)
  dat <- tryCatch({
    if (grepl("\\.gz$", file_path)) {
      read.delim(gzfile(file_path), header = FALSE, skip = 1)
    } else {
      read.delim(file_path,   header = FALSE, skip = 1)
    }
  }, error = function(e) {
    message("Error reading ", file_path, ":", e$message)
    return(NULL)
  })
  if (is.null(dat)) next
  
  # 6) Add metadata columns
  dat$sample <- samp          # first 10 chars
   
  # 7) Store in our list
  sample_list[[samp]] <- dat
}

# 8) Merge them all by row
dfmore3 <- do.call(rbind, sample_list)

# 9) Inspect
head(dfmore3)
nrow(dfmore3)

ITdatamerge$V5 <- NULL
ITdatamerge$Vchr <- NULL
ITdatamerge

dfmore3$V1 <- paste(ITdatamerge$V1, ITdatamerge$V2, sep = "_")  
dfmore3$V2 <- NULL

head(dfmore3)


 
# filetemp has run_accession and all the other run‐metadata columns
# ITdatamerge has sample plus your methylation columns

# Rename the key in dfmore3 so it matches filetemp
names(dfmore3)[names(dfmore3)=="sample"] <- "run_accession"


merged_all <- merge(
  dfmore3,
  filetemp,
  by    = "run_accession",
  all.x = TRUE,
  sort  = FALSE
)
 

# now merge
merged_all <- merge(
  merged_all, 
  excelLOC2,
  by.x  = "sample_alias",   # in ITdatamerge
  by.y  = "Sample",   # in excelLOC
  all.x = TRUE,       # keep all rows of ITdatamerge
  sort  = FALSE
)


# inspect
head(merged_all)
  
 
row.names(dfmore3) <- NULL 
dfmore3
View(merged_all)
 
#noice

testmerge <- merged_all
write.csv(testmerge, file = "testmerge2.csv")

 
testmerge$Cultural.Affiliation <- NULL 
testmerge$Skeletal.Element <- NULL 
hist(ITmatchtable$Freq)

ggplot(testmerge, aes(x = Period.Label.for.Analyses)) +
  geom_bar(fill = "#69b3a2") +       # geom_bar defaults to counting rows per level
  labs(x = "Period", y = "Count") +
  theme_minimal()

head(testmerge)

testmerge <- merge(
  testmerge, 
  excelLOC1,
  by.x  = "Site",   # in ITdatamerge
  #by.y  = "Table.S1..Archaeological.site.information",   # in excelLOC
  all.x = TRUE,       # keep all rows of ITdatamerge
  sort  = FALSE
)

# inspect



head(testmerge)
v




View(head(testmerge))
testmerge$V5 <- testmerge$V3 / testmerge$V4
testmerge$mtDNA <- NULL 
testmerge$Radiocarbon.Dated <- NULL 
testmerge$SNP.Coverage <- NULL 
testmerge$sample_accession <- NULL 
testmerge$UDG.treatment <- NULL 
testmerge$experiment_accession <- NULL 
testmerge$study_accession <- NULL
testmerge$Contamination_MT <- NULL
testmerge$Contamination_Xchr <- NULL
testmerge$`Y-chr` <- NULL 
testmerge$`Damage.rate.on.last.base.(5p.G>A)` <- NULL 
testmerge$Archaeological.ID.and.Notes <- NULL 

names(testmerge)[names(testmerge) == "Period.according.to.archaeological.context"] <- "Period"
names(testmerge)[names(testmerge) == "Period.Label.for.Analyses"] <- "period4ana"
names(testmerge)[names(testmerge) 
== "Date.(Direct.radiocarbon.date.on.the.individual.calibrated.95%.confidence.interval.or.date.range.based.on.the.archaeological.context,.including.AMS.dating.of.stratigraphic.unit)"] <- "Date"

write.csv(testmerge, file = "testmergeclean.csv")
View(head(testmerge)) 

coverage_sum <- testmerge %>% group_by(sample_alias, period4ana) %>% summarise(Coverage = first(Coverage), .groups = "drop") %>% group_by(period4ana) %>% summarise(mean_cov = mean(Coverage))
coverage_sum
ggplot(coverage_sum, aes(period4ana, mean_cov)) + geom_boxplot()


testiinngggg <- testmerge %>% group_by(sample_alias, period4ana)
testiinngggg


View( testmerge )
 




