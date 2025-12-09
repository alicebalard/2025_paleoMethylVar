#turns out you can just use table to see how many of each there are  
MDperchr <- table(Mdata$chr)

# Convert it to a data frame
MDperchr_df <- as.data.frame(MDperchr)
colnames(MDperchr_df) <- c("chr", "count")

# Plot with ggplot2
ggplot(data = MDperchr_df, aes(x = chr, y = count)) +
  geom_bar(stat = "identity") +
  labs(title = "Count per Chromosome in Mdata",
       x = "Chromosome",
       y = "Count") +
  theme_minimal()

ggplot(data = MDperchr_df, aes(x = reorder(chr, -count), y = count)) +
  geom_bar(stat = "identity") +
  labs(title = "Count per Chromosome in Mdata",
       x = "Chromosome",
       y = "Count") +
  theme_minimal()


# cpg sites are conserved - conservation of cpg sites especially on chromasome 19 
# to further investigate there were 4125 cpgs included in the top 5% variance and involved in 60% of the data sets 
# potentially chromasome 19 cpg count can give information on past histories in methylation 

#first combine the ITaly data and then test the counts of chromasomal number 


MDperchr <- table(Mdata$chr)

# Convert it to a data frame
MDperchr_df <- as.data.frame(MDperchr)
colnames(MDperchr_df) <- c("chr", "count")

# Plot with ggplot2
ggplot(data = MDperchr_df, aes(x = chr, y = count)) +
  geom_bar(stat = "identity") +
  labs(title = "Count per Chromosome in Mdata",
       x = "Chromosome",
       y = "Count") +
  theme_minimal()

ggplot(data = MDperchr_df, aes(x = reorder(chr, -count), y = count)) +
  geom_bar(stat = "identity") +
  labs(title = "Count per Chromosome in Mdata",
       x = "Chromosome",
       y = "Count") +
  theme_minimal()

MDperchr <- table(Mdata$chr)

# Convert it to a data frame
MDperchr_df <- as.data.frame(MDperchr)
colnames(MDperchr_df) <- c("chr", "count")

# Plot with ggplot2
ggplot(data = MDperchr_df, aes(x = chr, y = count)) +
  geom_bar(stat = "identity") +
  labs(title = "Count per Chromosome in Mdata",
       x = "Chromosome",
       y = "Count") +
  theme_minimal()

ggplot(data = MDperchr_df, aes(x = reorder(chr, -count), y = count)) +
  geom_bar(stat = "identity") +
  labs(title = "Count per Chromosome in Mdata",
       x = "Chromosome",
       y = "Count") +
  theme_minimal()


plot(density(variance_results$Variance), main = "", xlab = "Value")

plot_ly(data = variance_results,
        x = ~Freq, y = ~Variance,
        type = 'scatter', mode = 'markers',
        text = ~paste("Position:", Var1),
        marker = list(size = 6, color = 'blue'))




 


plot(variance_results$Freq, variance_results$Variance,
     xlab = "Frequency", ylab = "Methylation Variance",
     main = "Methylation Variance vs Frequency",
     pch = 20, col = "darkblue")



hist(variance_results$Variance,
     breaks = 50,
     main = "Distribution of Methylation Variance",
     xlab = "Variance",
     col = "skyblue")


plot(density(variance_results$Variance, na.rm = TRUE),
     main = "Density of Methylation Variance",
     xlab = "Variance", col = "purple")


variance_results$FreqGroup <- cut(variance_results$Freq,
                                  breaks = c(0, 10, 20, 50, 100, Inf),
                                  labels = c("0-10", "11-20", "21-50", "51-100", "100+"))

boxplot(Variance ~ FreqGroup, data = variance_results,
        main = "Variance by Frequency Group",
        xlab = "Frequency Range", ylab = "Variance",
        col = "orange") 


top_vars <- variance_results[order(-variance_results$Variance), ][1:20, ]
barplot(top_vars$Variance, names.arg = top_vars$Var1,
        las = 2, col = "tomato", main = "Top 20 Variable Positions",
        ylab = "Variance", cex.names = 0.6)
 

top_vars

 
testo 

top_vars <- variance_results[order(-variance_results$Variance), ][1:20, ]
barplot(top_vars$Variance, names.arg = top_vars$Var1,
        las = 2, col = "tomato", main = "Top 20 variance",
        ylab = "Variance", cex.names = 0.6)


top_vars
var_res_1_head

# with linear trend
p2 <- ggplot(var_res_1_head, aes(x=Variance , y=Freq)) +
  geom_point() +
  geom_smooth(method=lm , color="red", se=FALSE) 
p2      
 


#CPG count for each chromasome in marias data   
MDperchr <- table(Mdata$chr)

# Convert it to a data frame
MDperchr_df <- as.data.frame(MDperchr)
colnames(MDperchr_df) <- c("chr", "count")

# Plot with ggplot2
ggplot(data = MDperchr_df, aes(x = chr, y = count)) +
  geom_bar(stat = "identity") +
  labs(title = "Count per Chromosome in Mdata",
       x = "Chromosome",
       y = "Count") +
  theme_minimal()

ggplot(data = MDperchr_df, aes(x = reorder(chr, -count), y = count)) +
  geom_bar(stat = "identity") +
  labs(title = "Count per Chromosome in Mdata",
       x = "Chromosome",
       y = "Count") +
  theme_minimal()
# cpg sites are conserved - conservation of cpg sites especially on chromasome 19 
# to further investigate there were 4125 cpgs included in the top 5% variance and involved in 60% of the data sets 
# potentially chromasome 19 cpg count can give information on past histories in methylation 

#first combine the ITaly data and then test the counts of chromasomal number 

top_counts

top_counts <- count_df %>% slice_max(Count, n = 20)

ggplot(top_counts, aes(x = reorder(BaseID, -Count), y = Count)) +
  geom_bar(stat = "identity") +
  labs(title = "Top 20 ITdata Groups by Count", x = "ITdata ID", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




var_chr_period_new
 
Varrank <- var_chr_period_new[order(-var_chr_period_new$var_beta), ]
Varrank
dprank <- var_chr_period_new[order(-var_chr_period_new$n_datapts), ]
dprank

library(ggplot2)
ggplot(data = var_chr_period_new) + 
  geom_point(
    mapping = aes(x = n_datapts, y = var_beta)
  )


nrow(var_chr_period)
head(var_chr_period)
var_chr_period
 
  


