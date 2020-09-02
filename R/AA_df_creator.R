###############################################################################
## Rscript to create output csv of AAs for future RF and MN model fitting #####
###############################################################################
require(seqinr)
require(stringr)

aa_adder <- function(df, list_of_files, gene){
  missing_mic_data <- NULL
  
  if(gene == "1a"){
    col_start = 2
    col_end = 253
    tpd_length = 252
    col_vec = seq(2,253)
  }else if(gene == "2b"){
    col_start = 254
    col_end = 530
    tpd_length = 277
    col_vec = seq(col_start, col_end)
  }else if (gene == "2x"){
    col_start = 531
    col_end = 829
    tpd_length = 299
    col_vec = seq(col_start, col_end)
  }
  
  for (k in 1:length(list_of_files)){
    current_file <- list_of_files[k]
    basename_o_file <- basename(current_file)
    isolate_name <- str_split_fixed(basename_o_file, "_pb", 2)[,1]
    # if (k == 582){
    #    browser()
    # }
    db_index <- which(df[,1] == isolate_name)
    
    if(length(db_index) == 1){
      
      current_fasta <- read.fasta(file = current_file,
                                  seqtype = "AA")
      
      if(length(current_fasta[[1]]) == tpd_length){
        for(j in 1:length(current_fasta[[1]])){
          df[db_index,col_vec[j]] <- current_fasta[[1]][[j]] 
          
        }
      }
    }else{
      print(isolate_name)
      if (!(isolate_name %in% missing_mic_data)){
        missing_mic_data <- append(isolate_name, missing_mic_data)
      }
      
    }
    
    
    print(k)
    # Sys.sleep(0.1)
  }
  
  return(df)
}

###############################################################################
## Load up the input data strings #############################################
###############################################################################
input_args <- commandArgs(trailingOnly = TRUE)
pbp1a_files <- readLines(input_args[1])
pbp2b_files <- readLines(input_args[2])
pbp2x_files <- readLines(input_args[3])

out_csv_name <- input_args[4]

## Get the total number of sequences to run through (some might be missing some genes)

pbp1a_basename <- basename(pbp1a_files)
pbp2b_basename <- basename(pbp2b_files)
pbp2x_basename <- basename(pbp2x_files)

pbp1a_isolates <- stringr::str_split_fixed(string = pbp1a_basename, "_pbp",2)[,1]
pbp2b_isolates <- stringr::str_split_fixed(string = pbp2b_basename, "_pbp",2)[,1]
pbp2x_isolates <- stringr::str_split_fixed(string = pbp2x_basename, "_pbp",2)[,1]

tot_isolates <- union(pbp1a_isolates, pbp2b_isolates)
tot_isolates <- union(tot_isolates, pbp2x_isolates)

###############################################################################
## Set up df to collate AA tpds ###############################################
###############################################################################

one_a <- paste("a1_", seq(1,252), sep = "_")
two_b <- paste("b2_", seq(1, 277), sep = "_")
two_x <- paste("x2_", seq(61,359), sep = "_")

amino_acid_db <- data.frame(data = matrix(data = NA, ncol = 829, nrow = length(tot_isolates)))
colnames(amino_acid_db) <- c("id",one_a, two_b, two_x)
amino_acid_db$id <- tot_isolates


###############################################################################
## Ok, now lets run that through for all the isolates #########################
###############################################################################

amino_acid_db <- aa_adder(amino_acid_db, pbp1a_files, "1a")
amino_acid_db <- aa_adder(amino_acid_db, pbp2b_files, "2b")
amino_acid_db <- aa_adder(amino_acid_db, pbp2x_files, "2x")


write.csv(amino_acid_db, file = out_csv_name)

