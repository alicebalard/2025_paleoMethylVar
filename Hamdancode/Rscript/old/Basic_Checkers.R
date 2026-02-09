#All the different positions of the chromasomes 
#count numb of chromasome 1 before adding to count   

run = 0 
check = 0
same1 <- c()
x <- 1


count <- 0
for (i in 1:nrow(Mdata)) {
  if (Mdata$chr[i] == "chr1") {
    count <- count + 1 
  }  
}

count

#this one was done with an and/or to get values and stuff
#you can change the values depending on what chromasomes you want to count 
count <- 0

for (i in 1:nrow(HGdata)) {
  if (isTRUE(Mdata$chr[i] == "chr1" & HGdata$V1[i] == "1")) {
    count <- count + 1 
  }  
}
print(count)

#Checker 2

HGmatches <- data.frame()

HGmatches

chrpos

Mmatches <- data.frame()

check = 0 
for (i in 1:nrow(ITdata63)) {
  if (ITdata63$V2[i] %in% chrpos) {
    check <- check + 1
    #print(HGdata$V2[i])
    Mmatches <- rbind(Mmatches, ITdata63[i, ])
    
  }
}

Mmatches
nrow(Mmatches)

Mmatches <- data.frame()
ITmatches <- data.frame()


checker <- function(datar) {
  
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

checker(ITdatamerge) 
nrow(ITmatches)
