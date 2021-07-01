###################################################
## comparing the hmm aa df to the gene name method 
###################################################

require(dplyr)


gene_name_aa <- read.csv("~/Dropbox/phd/PMEN3/test_pen/updated_whole_cdc_aa/aa_df.csv",
                         stringsAsFactors = FALSE)
hmm_aa <- read.csv("~/Dropbox/phd/PMEN3/test_pen/aa_df.csv", 
                   stringsAsFactors = FALSE)


## which ones missing from the hmm_aa df?

length(which(gene_name_aa$id %in% hmm_aa$id))
length(which(!(gene_name_aa$id %in% hmm_aa$id)))
hmm_aa[which(!(hmm_aa$id %in% gene_name_aa$id)),"id"]
gene_name_aa[which(!(gene_name_aa$id %in% hmm_aa$id)),"id"]

### So there are 3 missing from the hmm aa df and one in the hmm aa 
## For 6187_7#15 no pbp2b hit
##     6569_6#20 no pbp1a hit
##     6678_3#11 no pbp2b hit
##     6736_4#14 no 1a hit      
## All these were too small to extract properly at the moment. 
### df that isn't in the gene name one

#############################################################################
## Now lets look into whether the AA predictions are the same for the #######
## different isolates #######################################################
#############################################################################


df_compo <- function(df_one, df_two){
  ## Use the set diff function to compare joint rows 
  ## Assume df_one is the smaller hmm one 
  
  odd_ids <- NULL
  
  for(k in 1:nrow(df_one)){
    current_id <- df_one[k, "id"]
    current_one <- df_one[k,]
    current_two <- df_two[df_two$id == current_id,]
    if(!nrow(current_two)){
      next
    }
    
    if(nrow(setdiff(current_one, current_two)))
      odd_ids <- append(odd_ids, current_id)
    
  }
  
  
  return(odd_ids)
  
}


oddities <- df_compo(hmm_aa, gene_name_aa)

hmm_odd <- hmm_aa[hmm_aa$id == oddities,] %>% mutate(dataset = "hmm")
gene_name_odd <- gene_name_aa[gene_name_aa$id == oddities,] %>% mutate(dataset = "gene_name")
compo <- bind_rows(hmm_odd, gene_name_odd) %>% t() %>% as.data.frame() %>% rename(HMM = V1) %>%
  rename(gene_name = V2) %>% mutate(diff = ifelse(HMM == gene_name, "same","differ"))
dplyr::count(compo, diff)
compo[compo$diff == "differ",]




