# Standardise function
# From Mottl et al (2021) RRatepol -- https://rdrr.io/github/HOPE-UIB-BIO/R-Ratepol-package/src/R/fc_standardise_community_data.R
fc_standardise_community_data <- 
  function (sequence_pollen, sequence_sample_id, N_individuals, Debug = FALSE){
    
    Samples <-  1:nrow(sequence_pollen) 
    Samples <-  Samples[!is.na(sequence_sample_id)]
    
    for(i in Samples){  # for each row(sample)
      
      select_row <-  sequence_pollen[i, ] # selected row
      
      n1 <-  1:ncol(sequence_pollen) #number for each species name in community data
      ab1 <-  as.vector(select_row) #frequencies of species in each sample
      
      vec1 <-  NULL    #a vector for the species pool
      
      # create a vector with species numbers replicated X times, where X is number of individuals
      for(j in 1:length(n1))  {
        
        #1: number of species 
        v1 <-  rep(n1[j], ab1[j]) #repeat species names 
        vec1 <-  c(vec1, v1) #a vector that repeat the species for the occurrences
      }
      
      # sample species X time (150 times)
      rsample <-  sample(vec1, size = N_individuals, replace = FALSE)
      
      # replace all values in community data by 0
      sequence_pollen[i, ]<-  rep(0,length(select_row))
      
      # replace individuals by new randomised values
      sequence_pollen[i,as.numeric(names(table(rsample)))] <- as.numeric(table(rsample))
      
    }
    
    return (sequence_pollen) 
    
  }
