###############################################################################
## Extract only the gffs and fasta for the 1220 in the pmen_full data set #####
###############################################################################

whole_aa_data_pmen <- read.csv(file = "~/Dropbox/phd/cdc_mic_predictions/all_isolates_with_three_pbp_profiles_db_updated.csv",
                               stringsAsFactors = TRUE)


pmen_gff <- readLines("~/Dropbox/phd/insertion_site_analysis/data/pmen_run/pmen_gff_list.txt")
pmen_fasta <- readLines("~/Dropbox/phd/insertion_site_analysis/data/pmen_run/pmen_fastas_list")

pmen_bases_gff <- basename(pmen_gff)
pmen_bases_fasta <- basename(pmen_fasta)

pmen_gff_isos <- sub("\\..*$","",pmen_bases_gff)
pmen_fasta_isos <- sub("\\..*$","",pmen_bases_fasta)

missing_fasta <- (which(!(pmen_fasta_isos %in% pmen_gff_isos)))
pmen_fasta_isos <- pmen_fasta_isos[-missing_fasta]
pmen_fasta <- pmen_fasta[-missing_fasta]

## Split on the dash, then split on the dot or a combination of dot and grep

ordered_gff <- sort(pmen_gff_isos)
ordered_fasta <- sort(pmen_fasta_isos)

orderer_function <- function(ordered_1, unordered_1){
  place_list <- NULL
  for(k in 1:length(ordered_1)){
    place_list <- append(place_list,which(unordered_1 == ordered_1[k]))
  }
  
  return(place_list)
}




gff_place_list <- orderer_function(ordered_gff, pmen_gff_isos)
fasta_place_list <- orderer_function(ordered_fasta, pmen_fasta_isos)

identical(pmen_gff_isos[gff_place_list], pmen_fasta_isos[fasta_place_list])

writeLines(pmen_gff[gff_place_list], con = "~/Dropbox/phd/insertion_site_analysis/data/pmen_run/pmen_gff_list.txt")
writeLines(pmen_fasta[fasta_place_list], con = "~/Dropbox/phd/insertion_site_analysis/data/pmen_run/pmen_fasta_list_narrow.txt")

missing_places <- c(grep("6187_7#15$",ordered_fasta), grep("6678_3#11$",ordered_fasta),
                    grep("6187_5#19$",ordered_fasta), grep("11679_1#89$",ordered_fasta))

writeLines(pmen_gff[gff_place_list][missing_places], con = "~/Dropbox/phd/insertion_site_analysis/data/pmen_run/pmen_gff_list_4.txt")
writeLines(pmen_fasta[fasta_place_list][missing_places], con = "~/Dropbox/phd/insertion_site_analysis/data/pmen_run/pmen_fasta_list_narrow_4.txt")

test_6187_5_19_seq <- whole_aa_data_pmen[whole_aa_data_pmen$id == "6187_5#19",2:253][1,]
test_6187_5_19_seq <- apply(test_6187_5_19_seq,1,as.character) 
test_6187_5_19_seq <- ifelse(test_6187_5_19_seq == "TRUE","T",test_6187_5_19_seq)
test_6187_5_19_seq <- paste(test_6187_5_19_seq, collapse = "")

test_aa_seq <- as.AAbin(test_6187_5_19_seq)
write.fasta(test_aa_seq, names = "6187_5#19_pbp1a_seq",file.out = "~/Dropbox/phd/cdc_mic_predictions/pmen_pbp_extraction/testing_on_4/6187_5#19_old_pbp1a.prot")

missing_fasta <- (which(!(pmen_fasta_isos %in% pmen_gff_isos))) 
pmen_fasta_isos <- pmen_fasta_isos[-missing_fasta]

identical(sort(pmen_gff_isos), sort(pmen_fasta_isos))



fasta_update <- pmen_fasta[-missing_fasta]


pmen_gff <- sort(pmen_gff)



writeLines(pmen_fasta[-missing_fasta], con = "~/Dropbox/phd/insertion_site_analysis/data/pmen_run/pmen_fasta_list_narrow.txt")

identical(pmen_gff_isos, pmen_fasta_isos[-missing_fasta])


## read in aa df 

true_changer <- function(cdc_dataset){
  for(k in 2:(ncol(cdc_dataset) - 1)){
    current_col <- cdc_dataset[,k]
    true_vals <- which(current_col == TRUE)
    if(length(true_vals) > 0){
      levels(current_col) <- c(levels(current_col), "T")
      current_col[true_vals] <- "T"
      current_col <- as.factor(current_col)
      
      current_col <- droplevels(current_col)
      
    }
    
    false_vals <- which(current_col == FALSE)
    if(length(false_vals) != 0){
      levels(current_col) <- c(levels(current_col), "F")
      current_col[false_vals] <- "F"
      current_col <- as.factor(current_col)
      
      current_col <- droplevels(current_col)
      
    }
    cdc_dataset[,k] <- current_col
  }
  return(cdc_dataset)
}


aa_df <- read.csv("~/Dropbox/phd/cdc_mic_predictions/pmen_pbp_extraction/pmen_aa_df.csv",
                  stringsAsFactors = TRUE)
aa_df_2 <- read.csv("~/Dropbox/phd/cdc_mic_predictions/pmen_pbp_extraction/aa_df.csv",
                  stringsAsFactors = TRUE)


aa_df_2 <- aa_df_2[complete.cases(aa_df_2),]
aa_df_2 <- true_changer(aa_df_2)

whole_aa_data_pmen <- true_changer(whole_aa_data_pmen)

## test this against the old whole_aa_data_pmen 

tester_func <- function(new_aa_df, old_aa_df){
  
  old_a1_names <- colnames(old_aa_df)[grep("a1",colnames(old_aa_df))]
  old_b2_names <- colnames(old_aa_df)[grep("b2",colnames(old_aa_df))]
  old_x2_names <- colnames(old_aa_df)[grep("x2",colnames(old_aa_df))]
  
  match_df <- as.data.frame(matrix(ncol = 4, nrow = nrow(new_aa_df)))
  colnames(match_df) <- c("id","a1","b2","x2")
  
  for(k in 1:nrow(new_aa_df)){
    current_id <- as.character(new_aa_df$id[k])
    old_row <- old_aa_df[old_aa_df$id == current_id,]
    if(nrow(old_row) > 0){
      a1_match <- identical(sapply(old_row[,old_a1_names], as.character), sapply(new_aa_df[k,old_a1_names], as.character))
      b2_match <- identical(sapply(old_row[,old_b2_names], as.character), sapply(new_aa_df[k,old_b2_names], as.character))
      x2_match <- identical(sapply(old_row[,old_x2_names], as.character), sapply(new_aa_df[k,old_x2_names], as.character))
      
      print(paste("pbp 1a match:", a1_match))
      print(paste("pbp 2b match:", b2_match))
      print(paste("pbp 2x match:", x2_match))
      
      match_df[k,1] <- current_id
      match_df[k,2] <- a1_match
      match_df[k,3] <- b2_match
      match_df[k,4] <- x2_match
    }else{
      match_df[k,1] <- current_id
    }
    
    
  }

  return(match_df)
}

tester_output <- tester_func(aa_df_2, whole_aa_data_pmen)

cdc_in_mod <- read.csv("~/Dropbox/phd/PMEN3/test_pen/pbp_tpd_extraction/data/cdc_seqs_df.csv")

ncol(cdc_data)






