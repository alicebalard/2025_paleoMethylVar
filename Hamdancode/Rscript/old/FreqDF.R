#Everything is aprox for meth ratio 

# take the first value from ITloccount 
# cycle through ITdatamerge and find all of the same cgp location 


# calculate the variance for each chromasome
# I need it to do it for all of the Var1, going down the table
# find all instances of the Var1, turn methylation into a value between 0  - 1 


#SLIGHT mistake, the positions could be counting those on all chromasomes
#tmrw change this so it makes more sense 

 


# List which cites are found most often in IT data, Than 
# work out which ones have the most variation / which ones occur the most .

ITdatamerge
  

ITdatamerge$V2
ITlocCount <- table(ITdatamerge$V2)
ITlocCount

ITlocCount <- as.data.frame(ITlocCount)
ITlocCount

ITlocCount <- ITlocCount[order(-ITlocCount$Freq), ]
 

##REDO to include the chromasomal number aswell 
 


ITGreatCount <- ITlocCount2[freq_df$frequency > 10, ]

calc_position_variances <- function(ITGreatCount, ITdatamerge) {
  # Initialize the results data frame
  results <- data.frame(Var1 = numeric(), Variance = numeric(), Freq = numeric(), stringsAsFactors = FALSE)
  
  for (i in 1:nrow(ITGreatCount)) {
    position <- ITGreatCount$Var1[i]  # Position from ITGreatCount
    
    # Match positions in ITdatamerge using v2
    subset_data <- ITdatamerge[ITdatamerge$V2 == position, ]
    
    if (nrow(subset_data) > 0) {
      # Calculate methylation ratios safely
      ratios <- with(subset_data, ifelse(V4 != 0, V3 / V4, NA))
      var_val <- var(ratios, na.rm = TRUE)
    } else {
      var_val <- NA
    }
    
    # Store result
    results <- rbind(results, data.frame(
      Var1 = position,
      Variance = var_val,
      Freq = ITGreatCount$Freq[i]
    ))
  }
  
  # Sort by descending frequency
  results <- results[order(-results$Freq), ]
  return(results)
}

# Example usage:
variance_results <- calc_position_variances(ITGreatCount, ITdatamerge)
head(variance_results)
variance_results


all_ITmatches

sum(ITGreatCount$Var1 %in% ITdatamerge$V2)

nrow(ITGreatCount)


#Check that its not clumping them together regardless of chromasome number 
#you want it to be both chromasome and position consistently 

#REDO incling for ITlocCount to include chromasome
ITlocCount2 <- data.frame(ITdatamerge$V1, ITdatamerge$V2)
ITlocCount2

ITlocCount2 <- as.data.frame(ITlocCount2)

ITlocCount2

freq_df <- ITlocCount2 %>%
  group_by(ITdatamerge.V1, ITdatamerge.V2) %>%
  summarise(frequency = n(), .groups = 'drop')

 
#make it 3 instead of 10 
 
 
FreqGreatCount <- freq_df[freq_df$frequency > 4, ]

calc_position_variances <- function(FreqGreatCount, ITdatamerge) {
  # Initialize the results data frame
  results <- data.frame(Var1 = numeric(), Variance = numeric(), Freq = numeric(), CHR1 = numeric(),  stringsAsFactors = FALSE)
  
  for (i in 1:nrow(FreqGreatCount)) {
    position <- FreqGreatCount$ITdatamerge.V2[i]  # Position from ITGreatCount
    
    # Match positions in ITdatamerge using v2
    subset_data <- ITdatamerge[ITdatamerge$V2 == position, ]
    
    if (nrow(subset_data) > 0) {
      # Calculate methylation ratios safely
      ratios <- with(subset_data, ifelse(V4 != 0, V3 / V4, NA))
      var_val <- var(ratios, na.rm = TRUE)
    } else {
      var_val <- NA
    }
    
    # Store result
    results <- rbind(results, data.frame(
      CHR1 = FreqGreatCount$ITdatamerge.V1[i],
      Var1 = position,
      Variance = var_val,
      Freq = FreqGreatCount$frequency[i]
      
    ))
  }
  
  # Sort by descending frequency
  results <- results[order(-results$Freq), ]
  return(results)
}

# Example usage:
variance_results <- calc_position_variances(FreqGreatCount, ITdatamerge)
head(variance_results)
 
#find how many sites in marias 

#leveneTest
testo <- subset(ITdatamerge, V2 %in% variance_results$Var1)
head(testo)
ITdatamerge
variance_results
for( i in 1:nrow(testo)){ 
  testo$V6 <- testo$V3 / testo$V4  
  
  
  }
testo
leveneTest(testo$V6 ~ as.factor(testo$V2))
as.character(testo$V2)

#Significant result as the P value was really small 

#By chromasome 

testo$V1 <- as.factor(testo$V1)
leveneTest(V6 ~ V1, data = testo)

#there was significance between both the variance of 


#This was variance within now calculate the variance between all of the 
#positions on the chromosomes. 
variance_results
hist( x = variance_results$Freq, breaks = 50 )
hist( x = variance_results$Variance, breaks = 50 )



all_ITmatches
all_ITmatches$chr_pos <-  paste(all_ITmatches$V1, all_ITmatches$V2, sep = "_")  
all_ITmatches
nrow(all_ITmatches)



nrow(variance_results_1)
#you want to make a row with chrnum, actual position, a new row with the freq and further another 
# GODMC 
 
#scatter of the variance vs frequency 
#Variance results one will rank them based on the highest variance 
variance_results_1 <- variance_results[order(-variance_results$Variance ), ]
variance_results_1

variance_results_1$num_pos <- paste(variance_results_1$CHR, variance_results_1$Var1, sep = "_")  
 

var_res_1_head <- head(variance_results_1, 4000)
var_res_1_head

head(testmerge)
 


