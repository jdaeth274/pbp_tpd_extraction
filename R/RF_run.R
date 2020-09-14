###############################################################################
## Run RF on cdc aa data ######################################################
###############################################################################


require(randomForest)
require(seqinr)
require(protr)
require(BiocManager)
if(!("Biostrings" %in% rownames(installed.packages())))
  BiocManager::install("Biostrings")
require(Biostrings)
require(ape)
require(dplyr)
require(tibble)
require(ggplot2)
require(tidyr)
require(tictoc)
require(reshape2)
require(ranger)

###############################################################################
## First we'll code out the functions required ################################
###############################################################################

catergory_assessment <- function(mic_vector, S_upper, R_lower, transform_mic = FALSE){
  if(transform_mic)
    back_once_again <- 2^(mic_vector)
  else
    back_once_again <- mic_vector
  out_vec <- rep(0, length(back_once_again))
  for(k in 1:length(back_once_again)){
    current_val <- as.numeric(back_once_again[k])
    if (current_val <= S_upper | (all.equal(current_val, S_upper) == TRUE)){
      out_vec[k] <- "S"
    }else if(current_val >= R_lower | (all.equal(current_val, R_lower) == TRUE)){
      out_vec[k] <- "R"
    }else{
      out_vec[k] <- "I"
    }
    
  }
  return(out_vec)
  
}

CA_assessment <- function(actual_cat, predicted_cat){
  agree_vec <- 0
  disagree_vec <- 0
  major_erros <- 0
  very_major_errors <- 0
  nd_changer <- 0
  for(k in 1:length(actual_cat)){
    current_actual <- actual_cat[k]
    current_pred <- predicted_cat[k]
    if(current_actual == current_pred){
      agree_vec <- agree_vec + 1
      
    }else{
      if(current_actual == "ND"){
        nd_changer <- nd_changer + 1
      }else{
        disagree_vec <- disagree_vec + 1
        if(current_actual == "S" & grepl("R",current_pred)){
          major_erros <- major_erros + 1
        }else if (grepl("R",current_actual) & current_pred == "S"){
          very_major_errors <- very_major_errors + 1
        }
      }
    }
  }
  CA_val <- agree_vec / (agree_vec + disagree_vec) * 100
  majjor_erro <- major_erros / (agree_vec + disagree_vec) * 100
  very_major_errors <- very_major_errors / (agree_vec + disagree_vec) * 100
  
  return(list(CA = CA_val, maj = majjor_erro, vmj = very_major_errors,
              ND_changes = nd_changer))
}

data("BLOSUM62")
blossy <- as.data.frame(BLOSUM62)


blossy_changer <- function(validation, training_purp, blossy){
  # for(k in 1:(ncol(training_purp) - 1)){
  #   training_purp[,k] <- as.character(training_purp[,k])
  # }
  # for(k in 1:(ncol(validation) - 1)){
  #   validation[,k] <- as.character(validation[,k])
  # }
  # 
  
  changed_row <- NULL
  
  
  for (k in 1:(ncol(validation) - 1)){
    num_changed <- 0
    validation[,k] <- droplevels(validation[,k])
    current_col <- training_purp[,k]
    train_vals <- plyr::count(current_col)
    valid_vals <- plyr::count(validation[,k])
    
    for(l in 1:nrow(valid_vals)){
      current_aa <- as.character(valid_vals[l,1])
      if(!(current_aa %in% train_vals[,1])){
        num_changed <- num_changed + 1
        current_aa <- as.character(current_aa)
        training_rows <- which(row.names(blossy) %in% train_vals[,1])
        valid_col <- which(colnames(blossy) == current_aa)
        
        mat_vals <- blossy[training_rows, valid_col]
        single_val <- which.max(mat_vals)
        
        new_aa <- row.names(blossy)[training_rows[single_val]]
        
        levels(validation[,k]) <- c(levels(validation[,k]), new_aa)
        
        
        validation[validation[,k] == current_aa,k] <- new_aa
        
        validation[,k] <- droplevels(validation[,k])
        
        
        
      }
    }
    
    changed_row <- append(changed_row, num_changed)
    
  }
  
  
  common <- intersect(names(training_purp), names(validation)) 
  for (p in common) { 
    # if (p == "a1__16"){
    #   browser()
    # }
    # 
    if (class(training_purp[[p]]) == "factor") { 
      training_levels <- levels(training_purp[[p]]) 
      valid_levels <- levels(validation[[p]])
      for(k in 1:length(training_levels)){
        if(!(training_levels[k] %in% valid_levels))
          levels(validation[[p]]) <- c(levels(validation[[p]]), training_levels[k])
      }
      
      validation[[p]] <- factor(validation[[p]], levels = training_levels)
      
    } 
  }
  
  
  # for(k in 1:(ncol(training_purp) - 1)){
  #   training_purp[,k] <- as.factor(training_purp[,k])
  # }
  # for(k in 1:(ncol(validation) - 1)){
  #   validation[,k] <- as.factor(validation[,k])
  #   validation[,k] <- droplevels(validation[,k])
  #   levels(validation[,k]) <- levels(training_purp[,k])
  # }
  # 
  
  
  
  return(list(data = validation, changes = changed_row))
}

blossy_changer_old <- function(validation, training_purp, blossy){
  # for(k in 1:(ncol(training_purp) - 1)){
  #   training_purp[,k] <- as.character(training_purp[,k])
  # }
  # for(k in 1:(ncol(validation) - 1)){
  #   validation[,k] <- as.character(validation[,k])
  # }
  # 
  
  changed_row <- NULL
  
  
  for (k in 1:(ncol(validation) - 1)){
    num_changed <- 0
    current_col <- training_purp[,k]
    train_vals <- plyr::count(current_col)
    valid_vals <- plyr::count(validation[,k])
    
    for(l in 1:nrow(valid_vals)){
      current_aa <- as.character(valid_vals[l,1])
      if(!(current_aa %in% train_vals[,1])){
        num_changed <- num_changed + 1
        current_aa <- as.character(current_aa)
        training_rows <- which(row.names(blossy) %in% train_vals[,1])
        valid_col <- which(colnames(blossy) == current_aa)
        
        mat_vals <- blossy[training_rows, valid_col]
        single_val <- which.max(mat_vals)
        
        new_aa <- row.names(blossy)[training_rows[single_val]]
        
        levels(validation[,k]) <- c(levels(validation[,k]), new_aa)
        
        
        validation[validation[,k] == current_aa,k] <- new_aa
        
        validation[,k] <- droplevels(validation[,k])
        
        
        
      }
    }
    changed_row <- append(changed_row, num_changed)
    
  }
  
  
  common <- intersect(names(training_purp), names(validation)) 
  for (p in common) { 
    # if (p == "a1__16"){
    #   browser()
    # }
    # 
    if (class(training_purp[[p]]) == "factor") { 
      training_levels <- levels(training_purp[[p]]) 
      valid_levels <- levels(validation[[p]])
      for(k in 1:length(training_levels)){
        if(!(training_levels[k] %in% valid_levels))
          levels(validation[[p]]) <- c(levels(validation[[p]]), training_levels[k])
      }
    } 
  }
  
  
  for(k in 1:(ncol(training_purp) - 1)){
    training_purp[,k] <- as.factor(training_purp[,k])
  }
  for(k in 1:(ncol(validation) - 1)){
    validation[,k] <- as.factor(validation[,k])
  }
  
  
  
  return(list(data = validation, changes = changed_row))
}


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

levels_checker <- function(valid_set, train_set){
  counter <- 0
  indexers <- NULL
  for(k in 1:(ncol(valid_set) -1 )){
    val_levels <- levels(valid_set[,k])
    train_levels <- levels(train_set[,k])
    if (all(train_levels %in% val_levels)){
      counter <- counter + 1
    }else{
      print(k)
      indexers <- append(indexers, k)
    }
    
  }
  return(indexers)
}

catergory_assessment_multi <- function(mic_vector, S_upper, R_lower, R2_lower, R3_lower, R4_lower, R5_lower){
  back_once_again <- 2^(mic_vector)
  out_vec <- rep(0, length(back_once_again))
  for(k in 1:length(back_once_again)){
    current_val <- as.numeric(back_once_again[k])
    if (current_val <= S_upper | (all.equal(current_val, S_upper) == TRUE)){
      out_vec[k] <- "S"
    }else if(current_val >= R_lower & current_val < R2_lower | (all.equal(current_val, R_lower) == TRUE)){
      out_vec[k] <- "R"
    }else if (current_val >= R2_lower & current_val < R3_lower | (all.equal(current_val, R2_lower) == TRUE)){
      out_vec[k] <- "R2"
    }else if (current_val >= R3_lower & current_val < R4_lower | (all.equal(current_val, R3_lower) == TRUE)){
      out_vec[k] <- "R3"
    }else if (current_val >= R4_lower & current_val < R5_lower | (all.equal(current_val, R4_lower) == TRUE)){
      out_vec[k] <- "R4"
    }else if(current_val >= R5_lower | (all.equal(current_val, R5_lower) == TRUE)){
      out_vec[k] <- "R5"
    }
    else{
      out_vec[k] <- "I"
    }
    
  }
  return(out_vec)
  
}

catergory_assessment_multi_smol <- function(mic_vector, S_upper, R_lower, R2_lower = 1, R3_lower = 2,  transform_mic = FALSE){
  
  if(transform_mic){
    back_once_again <- 2^(mic_vector)
  }else{
    back_once_again <- mic_vector
  }
  out_vec <- rep(0, length(back_once_again))
  for(k in 1:length(back_once_again)){
    current_val <- as.numeric(back_once_again[k])
    if (current_val <= S_upper | (all.equal(current_val, S_upper) == TRUE)){
      out_vec[k] <- "S"
    }else if(current_val >= R_lower & current_val < R2_lower | (all.equal(current_val, R_lower) == TRUE)){
      out_vec[k] <- "R"
    }else if (current_val >= R2_lower & current_val < R3_lower | (all.equal(current_val, R2_lower) == TRUE)){
      out_vec[k] <- "R2"
    }else if (current_val >= R3_lower | (all.equal(current_val, R3_lower) == TRUE)){
      out_vec[k] <- "R3"
    }
    else{
      out_vec[k] <- "I"
    }
    
  }
  return(out_vec)
  
}

difference_checker <- function(old_data, new_data, valid_data){
  
  change_rows <- old_data$changes
  out_list <- list()
  mismatched_old_new <- NULL
  levels_list <- NULL
  for(k in 1:ncol(old_data[[1]])){
    initial_dat <- as.data.frame(plyr::count(valid_data[,k]))
    colnames(initial_dat) <- c("id","initial")
    old_dat <- as.data.frame(plyr::count(old_data[[1]][,k]))
    colnames(old_dat) <- c("id","old")
    new_dat <- as.data.frame(plyr::count(new_data[[1]][,k]))
    colnames(new_dat) <- c("id","new")
    
    old_initial <- dplyr::left_join(initial_dat, old_dat, by = "id")
    old_initial_new <- dplyr::left_join(old_initial, new_dat, by = "id")
    out_list[[k]] <- old_initial_new
    
    if(!(identical(old_dat$old, new_dat$new))){
      mismatched_old_new <- append(k, mismatched_old_new)
    }
    if(!(identical(levels(old_data[[1]][,k]), levels(new_data[[1]][,k])))){
      levels_list <- append(k, levels_list)
    }
    
  }
  return(list(out_list, mismatched_old_new, levels_list))
}



rf_model_fitter <- function(input_data, test_data = NULL, test_known = FALSE, input_model = NULL, ranger_use = FALSE,
                            categories = 3){
  
  start_time <- Sys.time()
  
  if(categories == 3){
    cat_func <- catergory_assessment 
  }else if(categories == 5){
    cat_func <- catergory_assessment_multi_smol
  }
  
  tic("ensuring factors")
  
  for(i in 1:(ncol(input_data) - 1)){
    input_data[,i] <- as.factor(input_data[,i])
  }
  if(!is.null(test_data)){
    for(i in 1:(ncol(test_data) - 1)){
      test_data[,i] <- as.factor(test_data[,i])
    } 
  }
  
  toc()
  set.seed(5252)
  tic("prepping training and test data")
  if(is.null(test_data)){
    pmen_train_mic_rows <- sample(nrow(input_data), 0.7*nrow(input_data), replace = FALSE)
    pmen_train_mic <- input_data[pmen_train_mic_rows,]
    pmen_valid_mic <- input_data[-pmen_train_mic_rows,]
  }else{
    pmen_train_mic <- input_data
    pmen_valid_mic <- test_data
  }
  toc()
  tic("Running rf fit")
  if(is.null(input_model)){
    if(ranger_use){
      test_model_run <- ranger::ranger(mic ~., data = pmen_train_mic, num.trees = 500)
    }else{
      system.time(test_model_run <- randomForest(mic ~ ., data = pmen_train_mic,
                                                 ntree = 500, mtry = 828))
    }
  }else{
    test_model_run <- input_model
  }
  toc()
  tic("comparing to training set")
  
  pred_train <- predict(test_model_run, pmen_train_mic)
  if(ranger_use){
    actual_cats <- cat_func(pmen_train_mic$mic, S_upper = 0.06, R_lower = 0.12, transform_mic = TRUE)
    pred_cats <- cat_func(pred_train$predictions, S_upper = 0.06, R_lower = 0.12, transform_mic = TRUE)
  }else{
    actual_cats <- cat_func(pmen_train_mic$mic, S_upper = 0.06, R_lower = 0.12, transform_mic = TRUE)
    pred_cats <- cat_func(pred_train, S_upper = 0.06, R_lower = 0.12, transform_mic = TRUE)
  }
  CA_value_trained <- CA_assessment(actual_cats, pred_cats)
  toc()
  tic("Predicting on test data")
  CA_value_valid <- NULL
  if(is.null(test_data)){
    pmen_valid_mic <- blossy_changer(pmen_valid_mic, pmen_train_mic, blossy = BLOSUM62)
    pred_valid <- predict(test_model_run, pmen_valid_mic[[1]])
    if(ranger_use){
      actual_cats <- cat_func(pmen_valid_mic[[1]]$mic, S_upper = 0.06, R_lower = 0.12, transform_mic = TRUE)
      pred_cats <- cat_func(pred_valid$predictions, S_upper = 0.06, R_lower = 0.12, transform_mic = TRUE)
    }else{
      actual_cats <- cat_func(pmen_valid_mic[[1]]$mic, S_upper = 0.06, R_lower = 0.12, transform_mic = TRUE)
      pred_cats <- cat_func(pred_valid, S_upper = 0.06, R_lower = 0.12, transform_mic = TRUE)
    }
    CA_value_valid <- CA_assessment(actual_cats, pred_cats)
    
  }else{
    pmen_valid_mic <- blossy_changer(pmen_valid_mic, pmen_train_mic, blossy = BLOSUM62)
    pred_valid <- predict(test_model_run, pmen_valid_mic[[1]])
    if(ranger_use){
      pred_cats <- cat_func(pred_valid$predictions, S_upper = 0.06, R_lower = 0.12, transform_mic = TRUE)
    }else{
      pred_cats <- cat_func(pred_valid, S_upper = 0.06, R_lower = 0.12, transform_mic = TRUE)
    }
    if(test_known){
      actual_cats <- cat_func(pmen_valid_mic[[1]]$mic, S_upper = 0.06, R_lower = 0.12, transform_mic = TRUE)
      CA_value_valid <- CA_assessment(actual_cats, pred_cats)
    }
    
  }
  toc()
  end_time <- Sys.time()
  
  trained <- bind_rows(CA_value_trained)
  if(is.null(test_data) | test_known){
    valid_d <- bind_rows(CA_value_valid)
    trained$data <- "Train"
    valid_d$data <- "Valid"
    trained <- bind_rows(trained, valid_d)
  }else{
    trained$data <- "Train"
  }
  
  cat(paste("Total time for run:",(end_time - start_time)))
  return(list(model = test_model_run, tests = trained, valid_preds = pred_cats))
  
}

writing_out_res <- function(pmen_res, preds, pmen_data, new_colname){
  
  pmen_3_epi_data <- pmen_res[[1]]
  pmen_9_epi_data <- pmen_res[[2]]
  
  if(new_colname %in% colnames(pmen_3_epi_data)){
    print("Warning: Deleting old colname that is the same as input colname PMEN3")
    pmen_3_epi_data <- pmen_3_epi_data[,-which(colnames(pmen_3_epi_data) ==new_colname)]
  }
  
  if(new_colname %in% colnames(pmen_9_epi_data)){
    print("Warning: Deleting old colname that is the same as input colname PMEN9")
    pmen_9_epi_data <- pmen_9_epi_data[,-which(colnames(pmen_9_epi_data) ==new_colname)]
  }
  
  
  pmen_preds <- cbind.data.frame(pmen_data$id[which(is.na(pmen_data$mic))], preds)
  colnames(pmen_preds) <- c("id","place_holder")
  
  pmen_3_epi_data$place_holder <- pmen_3_epi_data$pen_SIR_pmen
  pmen_3_epi_data <- pmen_3_epi_data %>%
    left_join(pmen_preds, by = "id") %>%
    mutate(place_holder = coalesce(place_holder.y, place_holder.x)) %>%
    select(-place_holder.y, -place_holder.x)
  
  pmen_9_epi_data$place_holder <- pmen_9_epi_data$pen_SIR_pmen__autocolour
  pmen_9_epi_data <- pmen_9_epi_data %>%
    left_join(pmen_preds, by = "id") %>%
    mutate(place_holder = coalesce(place_holder.y, place_holder.x)) %>%
    select(-place_holder.y, -place_holder.x) 
  
  colnames(pmen_3_epi_data)[which(colnames(pmen_3_epi_data) == "place_holder")] <- new_colname
  colnames(pmen_9_epi_data)[which(colnames(pmen_9_epi_data) == "place_holder")] <- new_colname
  
  return(list(pmen3 = pmen_3_epi_data, pmen9 = pmen_9_epi_data))
  
}

compo_df_create <- function(epi_data, col_names, pdf_out = "~/Dropbox/phd/PMEN3/pen_res/levels_for_models.pdf"){
  
  out_df <- NULL
  for(k in 1:length(col_names)){
    current_colname <- col_names[k]
    current_count <- as.data.frame(plyr::count(epi_data[,current_colname]))
    colnames(current_count) <- c("id", current_colname)
    if(k == 1){
      out_df <- current_count
    }else{
      out_df <- dplyr::left_join(out_df, current_count, by = "id")
    }
    
    
  }
  
  ## plot changes 
  
  plot_df <- melt(out_df, id.vars = "id")
  
  pdf(file = pdf_out, paper = "A4r", width = 10, height = 8)
  
  graph_out <- ggplot(data = plot_df, aes(x = variable, y = value, colour = id)) + geom_point() +
    theme(axis.text = element_text(angle = 90))
  print(graph_out)
  
  dev.off()
  
  return(list(df = out_df, ggplot_graph = graph_out))
  
}

testing_results_df <- function(list_of_df, list_of_trained, list_of_tests, model_type){
  
  out_df <- NULL
  for(k in 1:length(list_of_df)){
    current_df <- list_of_df[[k]]
    current_df$trained <- list_of_trained[k]
    current_df$test <- list_of_tests[k]
    current_df$model <- model_type
    out_df <- bind_rows(out_df, current_df)
  }
  
  return(out_df)
  
}


results_testing_plotter <- function(results_df, column_to_use,y_lab, pdf_out_file){
  
  pdf(file = pdf_out_file, paper = "A4r", width = 10, height = 8)
  col_name_of_col <- colnames(results_df)[column_to_use]
  min_val <- min(results_df[,column_to_use])
  max_val <- ifelse(column_to_use ==1 , 100, 10)
  unique_datsets <- unique(results_df$trained)
  for(k in unique_datsets){
    current_training_set <- results_df[results_df$trained == k,]
    current_training_set <- current_training_set %>% filter(!(data == "Train" & test != k))
    current_training_set$grouping_val <- paste(current_training_set$test, current_training_set$data, sep = "-")
    if(nrow(current_training_set) == 9){
      non_k_val <- current_training_set$test[which(current_training_set$test != k)][1]
      k_trains <- which(current_training_set$grouping_val == paste(k,"-Train",sep = ""))
      k_valids <- which(current_training_set$grouping_val == paste(k,"-Valid",sep = ""))
      non_k_trains <- which(current_training_set$grouping_val == paste(non_k_val,"-Valid",sep = "")) 
      current_training_set$ordering_vec <- 1
      current_training_set$ordering_vec[k_valids] <- 2
      current_training_set$ordering_vec[non_k_trains] <- 3
      
    }else{
      k_trains <- which(current_training_set$grouping_val == paste(k,"-Train",sep = ""))
      k_valids <- which(current_training_set$grouping_val == paste(k,"-Valid",sep = ""))
      current_training_set$ordering_vec <- 1
      current_training_set$ordering_vec[k_valids] <- 2
      
    }
    
    print(ggplot(data = current_training_set,
                 aes(x = model, y = current_training_set[,column_to_use],
                     group = ordering_vec, color = grouping_val, fill = grouping_val, order = ordering_vec)) +
            geom_bar(stat = "identity", position = "dodge") +
            xlab("Model type") + ylab(y_lab) + ggtitle(paste("CA on model trained on:",k)) +
            coord_cartesian(ylim = c(min_val, max_val)))
    
    
  }
  
  dev.off()
  
}

###############################################################################
## Right, now we'll run on the pmen known data to test fitting to the data ####
###############################################################################

input_args <- commandArgs(trailingOnly = TRUE)
test_set <- input_args[1]
training_set <- input_args[2]
number_of_categories <- input_args[3]
out_csv_name <- input_args[4]

## load up the AA csv from the bash script here 
whole_aa_data_sample <- read.csv(file = test_set,
                               stringsAsFactors = TRUE)

## Ensure T & F treated as characters not boolean
whole_aa_data_sample <- true_changer(whole_aa_data_sample)

## Remove any isolates with missing values 
whole_aa_data <- whole_aa_data_sample[complete.cases(whole_aa_data),]

cdc_data <- read.csv(training_set)

cdc_data <- true_changer(cdc_data)

###############################################################################
## Lets run this on the CDC data ##############################################
###############################################################################

## remove isolate id 
cdc_data_no_isolate <- cdc_data[, -1]
## get mic vals in log format 
cdc_data_no_isolate$mic <- log(cdc_data_no_isolate$mic, base = 2)

## Test on CDC data only 
cdc_known_mic <- rf_model_fitter(cdc_data_no_isolate, ranger_use = TRUE,
                                 categories = number_of_categories)

## Predict Unknown data now
cdc_preds <- rf_model_fitter(cdc_data_no_isolate,test_data = whole_aa_data,
                             ranger_use = TRUE, categories = number_of_categories)

### Predicted categories for penicillin below

predicted_categories <- cdc_preds$valid_preds

out_df <- cbind.data.frame(whole_aa_data$id, predicted_categories)
colnames(out_df) <- c("id","penicillin_cat")

write.csv(out_df,
          file = out_csv_name)
