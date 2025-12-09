
results_list <- list()  # store all checker results
results_list
for (i in 1:134) {
  df_name <- paste0("ITdata", i)
  df <- get(df_name)  # get the actual data frame object
  
  results <- checker(df)  # run the checker function
  
  # Save results in a list for later use
  results_list[[df_name]] <- results
}

results

results_list[["ITdata2"]]$ITmatches
results_list[["ITdata2"]]$Mmatches

all_ITmatches <- do.call(rbind, lapply(results_list, function(x) x$ITmatches))
all_Mmatches  <- do.call(rbind, lapply(results_list, function(x) x$Mmatches)) 


all_ITmatches
all_Mmatches 
nrow(all_ITmatches)
nrow(all_Mmatches)
