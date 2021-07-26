#input_args <- c("./full_seqs_aa_dir/pbp1a_prot_list","./full_seqs_aa_dir/pbp2b_prot_list","./full_seqs_aa_dir/pbp2x_prot_list","wasssup.csv")
require(stringr, quietly = TRUE)
require(seqinr, quietly = TRUE)

input_args <- commandArgs(trailingOnly = TRUE)
pbp1a_files <- readLines(input_args[1])
pbp2b_files <- readLines(input_args[2])
pbp2x_files <- readLines(input_args[3])
out_csv_name <- input_args[4]
pbp1a_basename <- basename(pbp1a_files)
pbp2b_basename <- basename(pbp2b_files)
pbp2x_basename <- basename(pbp2x_files)
pbp1a_isolates <- stringr::str_split_fixed(string = pbp1a_basename, "_pbp",2)[,1]
pbp2b_isolates <- stringr::str_split_fixed(string = pbp2b_basename, "_pbp",2)[,1]
pbp2x_isolates <- stringr::str_split_fixed(string = pbp2x_basename, "_pbp",2)[,1]

tot_isolates <- union(pbp1a_isolates, pbp2b_isolates)
tot_isolates <- union(tot_isolates, pbp2x_isolates)

# a1_lengs <- NULL
# for(isolate in pbp1a_files){
#   print(isolate)
#   current_fasta <- read.fasta(file = isolate,
#                               seqtype = "AA")
#   a1_lengs <- append(a1_lengs, length(current_fasta[[1]]))
#   
# }
# hist.data = hist(a1_lengs, breaks = seq(651,728,1))
# hist.data$counts[hist.data$counts != 0] = log10(hist.data$counts[hist.data$counts != 0])
# plot(hist.data)
# 
# b2_lengs <- NULL
# for(isolate in pbp2b_files){
#   print(isolate)
#   current_fasta <- read.fasta(file = isolate,
#                               seqtype = "AA")
#   b2_lengs <- append(b2_lengs, length(current_fasta[[1]]))
#   
# }
# range(b2_lengs)
# hist.data = hist(b2_lengs, breaks = seq(610,681,1))
# hist.data$counts[hist.data$counts != 0] = log10(hist.data$counts[hist.data$counts != 0])
# plot(hist.data)
# 
# 
# x2_lengs <- NULL
# for(isolate in pbp2x_files){
#   print(isolate)
#   current_fasta <- read.fasta(file = isolate,
#                               seqtype = "AA")
#   x2_lengs <- append(x2_lengs, length(current_fasta[[1]]))
#   
# }
# range(x2_lengs)
# hist.data = hist(x2_lengs, breaks = seq(682,750,1))
# hist.data$counts[hist.data$counts != 0] = log10(hist.data$counts[hist.data$counts != 0])
# plot(hist.data)
# 

one_a <- paste("a1_", seq(1,719), sep = "_")
two_b <- paste("b2_", seq(1, 680), sep = "_")
two_x <- paste("x2_", seq(1,750), sep = "_")



amino_acid_db <- data.frame(data = matrix(data = NA, ncol = 2150, nrow = length(tot_isolates)))
colnames(amino_acid_db) <- c("id",one_a, two_b, two_x)
amino_acid_db$id <- tot_isolates


## Lets just take 719 for a1, 680 for 2b and 750 for 2x
aa_adder <- function(df, list_of_files, gene){
  
  missing_mic_data <- NULL
  #browser()
  if(gene == "1a"){
    col_start = 2
    col_end = 720
    tpd_length = 719
    col_vec = seq(2,720)
  }else if(gene == "2b"){
    col_start = 721
    col_end = 1400
    tpd_length = 680
    col_vec = seq(col_start, col_end)
  }else if (gene == "2x"){
    col_start = 1401
    col_end = 2150
    tpd_length = 750
    col_vec = seq(col_start, col_end)
  }
  
  for (k in 1:length(list_of_files)){
    current_file <- list_of_files[k]
    basename_o_file <- basename(current_file)
    isolate_name <- stringr::str_split_fixed(string = basename_o_file, "_pbp",2)[,1]
    # if (k == 582){
    #    browser()
    # }
    
    
    
    db_index <- which(df[,1] == isolate_name)
    
    if(length(db_index) == 1){
      
      current_fasta <- read.fasta(file = current_file,
                                  seqtype = "AA")
      
      if(length(current_fasta[[1]]) >= tpd_length){
        for(j in 1:tpd_length){
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

amino_acid_db <- aa_adder(amino_acid_db, pbp1a_files, "1a")
amino_acid_db <- aa_adder(amino_acid_db, pbp2b_files, "2b")
amino_acid_db <- aa_adder(amino_acid_db, pbp2x_files, "2x")

missing_rows <- amino_acid_db[rowSums(is.na(amino_acid_db)) > 0,]
amino_acid_db <- amino_acid_db[complete.cases(amino_acid_db),]

if(nrow(missing_rows) > 0){
  writeLines(missing_rows$id, con = "missing_pbps_aa_df_isolates.txt")  
}

write.csv(amino_acid_db, file = out_csv_name, row.names = FALSE)




