##updating 121125 we will set an option to re-run and only change the threshold
##updating 120325 we will set the group
##updating 051925 this is for the peak calling of the treatment
##updating 011525 refer from Pablo script to call the ct ACR


# load libraries
library(dplyr)
library("edgeR")
library(Matrix)
library(gplots)
library(RColorBrewer)
library(irlba)
library(proxy)
library(png)
library(tidyverse)
library(ComplexHeatmap)
library(here)
library(rsample)
library(purrr)
library(data.table)
library(future)
library(furrr)
library(parallel)
library(argparse)

# load arguments
args <- commandArgs(T)

message('define the arguments')

# load arguments
args <- commandArgs(T)

input_data <- as.character(args[1])
##peak tn5 file

meta <- as.character(args[2])
##meta file

peak_file <- as.character(args[3])
##peak file

input_config_fl <- as.character(args[4])

input_output_dir <- as.character(args[5])

source(input_config_fl)

##updating 051925
##here we use the meta_slot as the treat group
#meta_slot <- target_cluster



print(input_data)
print(meta)
print(peak_file)




calculating_specificity <- function(x, threads=30){
    
    # add pseudo-count
    x <- x+1
    
    # convert to probability distribution
    p <- t(apply(x, 1, function(z){
        z/sum(z)
    }))
    hp <- apply(p, 1, function(z){
        z <- z[z > 0]
        -1*sum(z*(log2(z)))
        
    })
    sp <- apply(p, 2, function(z){
        hp - log2(z)
    })
    
    return(sp)
}

convert_to_sparse_matrix <- function(three_col_tribble, meta_slot_var) {
    
    CPM_matrix_prep <- three_col_tribble  %>% 
        dplyr::select(!!sym(meta_slot_var), geneID, grouped_CPM)

    three_col_prep <- CPM_matrix_prep  %>% 
        dplyr::rename("V2" = "geneID")  %>% 
        dplyr::rename("V1" = !!sym(meta_slot_var))


    # make sure bins/cells are factors
    three_col_prep$V1 <- factor(three_col_prep$V1)
    three_col_prep$V2 <- factor(three_col_prep$V2)


    # convert to sparseMatrix format
    sparse_count_matrix <- Matrix::sparseMatrix(i=as.numeric(three_col_prep$V1),
        j=as.numeric(three_col_prep$V2),
        x=as.numeric(three_col_prep$grouped_CPM),
        dimnames=list(levels(three_col_prep$V1),levels(three_col_prep$V2)))


    return(sparse_count_matrix)
    
}

# Function to generate null distribution
generate_null_distribution <- function(data, col_name) {
  # Determine the number of rows per class (equal for all classes)
  n_rows_per_class <- min(250)

  ## Replicate splitting is old and being abandoned in place of bootstrapping
  # Split the data into replicate1 and replicate2 groups based on the rep_values column
  #data_rep1 <- data %>% filter(data[[rep_values]] == "rep1")
  #data_rep2 <- data %>% filter(data[[rep_values]] == "rep2")

  # Group data by the specified column and select n_rows_per_class for each class 
  sampled_data <- data %>%
    group_by_at(col_name) %>%
    sample_n(n_rows_per_class,replace = TRUE) %>%
    ungroup()

  # Shuffle the column and make sure less than 20% of cells retain their original cluster ID
  #shuffled_column <- shuffle_and_assign_v3(sampled_data, col_name)
  #sampled_data[[col_name]] <- shuffled_column
   sampled_data[[col_name]] <- sampled_data[[col_name]][sample(nrow(sampled_data))]

  return(sampled_data)
}


generate_null_dist_values <- function(meta_data, meta_slot_var, raw_cpm_counts_all_genes,thread_num){
#    meta_slot_var <- c("final_annotation")
    options(dplyr.summarise.inform = FALSE)
    merged_meta_cpm_information <- left_join(meta_data, raw_cpm_counts_all_genes, by = c("cellID"), relationship = "many-to-many")  %>%
        group_by(!!sym(meta_slot_var), geneID)  %>%
        summarise(counts = sum(accessability, na.rm = TRUE))

    ### Alt CPM Calc
    merged_meta_cpm_information_copied <- merged_meta_cpm_information
    catch <- merged_meta_cpm_information_copied  %>%
        group_by(!!sym(meta_slot_var)) %>%
        group_map(~(edgeR::cpm(.x$counts, log = FALSE, group = .f)), .keep = TRUE)  %>%
        unlist()

    caught_values <- as_tibble(catch)
    see <- ungroup(merged_meta_cpm_information_copied)
    merged_meta_cpm_information_copied<- bind_cols(merged_meta_cpm_information_copied,caught_values)  %>% 
        rename(grouped_CPM = value)
    
    sparse_null_dist <- convert_to_sparse_matrix(merged_meta_cpm_information_copied, meta_slot_var)
    transposed_ACRs_by_ct <- as.matrix(t(sparse_null_dist))

    message("Generating Null Distribution ...")

    calculate_specificity <- calculating_specificity(transposed_ACRs_by_ct,threads=thread_num)
       
    return(calculate_specificity)
}

generate_null_dist_values_optimized <- function(meta_data, meta_slot_var, raw_cpm_counts_all_genes) {

  #message(paste0("Working on bootstrap ", counter))
    print(paste("Processing permutation:", Sys.getpid()))

  # Convert dataframes to data.tables
  setDT(meta_data)
  setDT(raw_cpm_counts_all_genes)
  
  # Replace left_join with merge
  merged_meta_cpm_information <- merge(meta_data, raw_cpm_counts_all_genes, by = "cellID", all.x = TRUE, allow.cartesian=TRUE)
  
  # Group by and summarise
  merged_meta_cpm_information[, counts := sum(accessability, na.rm = TRUE), by = c(meta_slot_var, "geneID")]

  # Calculate CPM
  merged_meta_cpm_information[, grouped_CPM := edgeR::cpm(counts, log = FALSE, group = get(meta_slot_var)), by = meta_slot_var]
  
  sparse_null_dist <- convert_to_sparse_matrix(merged_meta_cpm_information, meta_slot_var)
  transposed_ACRs_by_ct <- as.matrix(t(sparse_null_dist))

  message("Generating Null Distribution ...")

  calculate_specificity <- calculating_specificity(transposed_ACRs_by_ct,threads=thread_num)
       
  return(calculate_specificity)
}


generate_pvalues_bootstraps_fast <- function(meta_data, raw_cpm_counts_all_genes, meta_slot_var) {
  
  message(paste0("Working on bootstrap..."))
  #counter <<- counter +1 

  setDT(meta_data)
  setDT(raw_cpm_counts_all_genes)
  
  merged_meta_cpm_information <- merge(meta_data, raw_cpm_counts_all_genes, by = "cellID", all.x = TRUE, allow.cartesian=TRUE)
  
  merged_meta_cpm_information[, counts := sum(accessability, na.rm = TRUE), by = c(meta_slot_var, "geneID")]
  
  #message("generating the CPM values")
  
  merged_meta_cpm_information[, grouped_CPM := edgeR::cpm(counts, log = FALSE, group = get(meta_slot_var)), by = meta_slot_var]
  
  sparse_null_dist <- convert_to_sparse_matrix(merged_meta_cpm_information, meta_slot_var)
  transposed_ACRs_by_ct <- as.matrix(t(sparse_null_dist))
  
  #message("Generating P values ...")
  
  calculate_specificity <- calculating_specificity(transposed_ACRs_by_ct,threads = thread_num)
  return(calculate_specificity)
}

#parser <- ArgumentParser(description = "Command Line Args....")

is_valid_matrix <- function(mat) {
  # Check if it's a matrix
  if(!is.matrix(mat)) return(FALSE)
  
  # Check if it's empty
  if(any(dim(mat) == 0)) return(FALSE)
  
  # Check for negative values
  if(any(mat < 0, na.rm = TRUE)) return(FALSE)
  
  # If passed all checks
  return(TRUE)
}

calc_pvals <- function(qp, mean_val, sd) {
  obs <- qp[is.finite(qp)]
  #ave <- mean(all_null_values_array, na.rm = TRUE)
  #sd <- sd(all_null_values_array, na.rm = TRUE)
  pvals <- pnorm(obs, mean = mean_val, sd = sd, lower.tail = TRUE)
  return(pvals)
}

calculate_pvalues_per_ACR <- function(ACR_tribble_list) {
  results <- ACR_tribble_list %>%
    ungroup() %>%
    rowwise() %>%
    mutate(
      list_len = lengths(null_dist) + 1,
      median_val = mean(unlist(distribution), na.rm = TRUE),
      conf_interval_lower = quantile(unlist(distribution), probs = c(0.025), na.rm = TRUE),           
      conf_interval_upper = quantile(unlist(distribution), probs = c(0.975), na.rm = TRUE),           
      med_val_perm = sum(median_val > unlist(null_dist)),
      perm_pval = med_val_perm/list_len, # Permutation p-value for the median
      perm_pval_lower = sum(unlist(conf_interval_lower[1][1]) > unlist(null_dist))/list_len, # Permutation p-value for the lower CI bound
      perm_pval_upper = sum(unlist(conf_interval_upper[1][1]) > unlist(null_dist))/list_len, # Permutation p-value for the upper CI bound
      null_dist_mean = mean(unlist(null_dist), na.rm = TRUE),
      null_dist_sd = sd(unlist(null_dist), na.rm = TRUE),
      pnorm_pval = map(null_dist, ~calc_pvals(median_val, null_dist_mean, null_dist_sd)),
      z_score = (median_val - null_dist_mean) / null_dist_sd,
      zscore_pval = pnorm(z_score, lower.tail = TRUE)
    ) %>% # calculate Z score 
    mutate(pnorm_pval = map_dbl(pnorm_pval, ~ .x[[1]])) %>%
    ungroup()
  return(results)
} 

quantify_cell_type_specific_acrs <- function(signifigant_matrix, pval_slot, pval_filter, prefix) {
  
  # Save histogram of p-values
  hist_save <- paste0(opt_ct_dir,'/',prefix, ".pvalue_dist.png")
  png(hist_save)
  hist(as.numeric(as.matrix(signifigant_matrix[,pval_slot])), breaks = 100)
  dev.off()
  
  
  signifigant_matrix <- signifigant_matrix %>% 
    rename(pval = !!sym(pval_slot))
  
  print(head(signifigant_matrix))
  ct_specific_acrs_filtered_cell_types <- signifigant_matrix %>% 
    dplyr::filter(z_score < -1) %>%
    dplyr::filter(pval < pval_filter)
  print(head(ct_specific_acrs_filtered_cell_types))
  
  
  counts_ct <- ct_specific_acrs_filtered_cell_types  %>% 
    group_by(ACR_values) %>% 
    summarise(count_n_cell_types = n())  %>% 
    arrange(desc(count_n_cell_types))  %>% 
    mutate(class_acr = case_when(count_n_cell_types == 1 ~ "cts_acr",
                                 count_n_cell_types > 1 & count_n_cell_types <=3 ~ "ctr_acr", 
                                 count_n_cell_types > 3 ~ "broadly_accessible_acr"))
  
  
  
  # Check 1
  if (any(duplicated(counts_ct$ACR_values))) {
    stop("Error: Some ACRs are assigned to multiple categories.")
  }
  
  
  ## Different filtering to assign and filter different classes of ACRs
  cell_type_specific_restricted_ACRs <- counts_ct  %>% 
    dplyr::filter(class_acr == "cts_acr" | class_acr == "ctr_acr" )
  
  cell_type_specific_ACRs <- counts_ct  %>% 
    dplyr::filter(class_acr == "cts_acr")
  
  cell_type_restricted_ACRs <- counts_ct  %>% 
    dplyr::filter(class_acr == "ctr_acr" )
  
  
  ct_specific_acrs_filtered_cell_types.filtered <- ct_specific_acrs_filtered_cell_types  %>% 
    dplyr::filter(ACR_values %in% cell_type_specific_restricted_ACRs$ACR_values)
  
  
  write_delim(ct_specific_acrs_filtered_cell_types.filtered, file = paste0(opt_ct_dir,'/',prefix, ".pvalues.csv", collapse = "."), delim=",")
  
  
  ##########################################
  ###Isoalte write cell-type specific ACRs
  ##########################################
  ct_specific_acrs_filtered_cell_types_joined_names <- ct_specific_acrs_filtered_cell_types  %>% 
    ungroup() %>% 
    dplyr::filter(ACR_values %in% cell_type_specific_ACRs$ACR_values)  %>% 
    mutate(sc_acr_name = str_c(ACR_values, cell_type, sep = ";"))
  
  
  # Check 2
  if (nrow(unique(ct_specific_acrs_filtered_cell_types_joined_names["ACR_values"])) != nrow(cell_type_specific_ACRs)) {
    stop("Error: Mismatch in number of ACRs assigned to cell-type specific category.")
  }
  
  
  combined <- left_join(ct_specific_acrs_filtered_cell_types_joined_names, bed_file_read, by = c("ACR_values" = "acr_number"))
  
  # Write Cell Type Specific Peaks ALl
  cell_type_specific_acrs <- combined  %>% 
    dplyr::ungroup() %>% 
    #mutate_at(vars(!!sym(pval_slot)), character)  %>% 
    dplyr::select(chrom, start, stop, sc_acr_name, pval) 
  
  write_tsv(cell_type_specific_acrs, file = paste0(opt_ct_dir,'/',prefix, ".all_cts.ACRs.bed", collapse = "."), col_names = FALSE)
  
  # Write Cell Type Specific Peaks by Cell Type 
  combined  %>% 
    group_by(cell_type)  %>% 
    dplyr::select(chrom, start, stop, sc_acr_name, pval)  %>% 
    group_walk(~ write_delim(.x, paste0(opt_ct_dir,'/',prefix,".", .y$cell_type, ".cts.ACRs.bed", collapse = "."), delim = "\t", col_names = FALSE))
  
  
  ##########################################
  #Select cell-type restricted ACRs and Write
  ##########################################
  ctr_acrs_filtered_cell_types_joined_names <- ct_specific_acrs_filtered_cell_types  %>% 
    dplyr::filter(ACR_values %in% cell_type_restricted_ACRs$ACR_values) %>%
    group_by(ACR_values) %>% 
    summarize(combined_cell_type = paste(cell_type, collapse = ','),
              pval = paste(pval, collapse = ","))  %>% 
    ungroup()  %>% 
    mutate(sc_acr_name = str_c(ACR_values, combined_cell_type, sep = ";"))
  
  combined_ctr_acrs <- left_join(ctr_acrs_filtered_cell_types_joined_names, bed_file_read, by = c("ACR_values" = "acr_number"))
  
  # Write Cell Type restricted Peaks ALl
  cell_type_restricted_acrs <- combined_ctr_acrs  %>% 
    dplyr::select(chrom, start, stop, sc_acr_name, pval)
  
  write_tsv(cell_type_restricted_acrs, file = paste0(opt_ct_dir,'/',prefix, ".all_ctr.ACRs.bed", collapse = "."), col_names = FALSE)
  
  
  ##########################################   
  #Select Broadly Accessible ACRs and write 
  ##########################################
  '%ni%' <- Negate("%in%")
  broadly_accessible_acrs <- bed_file_read  %>% 
    dplyr::filter(acr_number %ni% cell_type_specific_restricted_ACRs$ACR_values)  %>% 
    dplyr::mutate(sc_acr_name = str_c(acr_number, "broadly_accessible", sep = ";"))  %>% 
    dplyr::select(-acr_number) %>%
    dplyr::select(-accessability) %>%
    mutate(pval = NA) %>% 
    dplyr::select(chrom, start, stop, sc_acr_name, pval)
  
  #write_tsv(broadly_accessible_acrs, file = paste0(prefix, "broadly_accessible.ACRs.bed", collapse = "."), col_names = FALSE)
  
  #Combine All ACRs and Write
  
  cell_type_specific_acrs <- cell_type_specific_acrs %>% 
    mutate(pval = as.character(pval))
  
  cell_type_restricted_acrs <- cell_type_restricted_acrs %>% 
    mutate(pval = as.character(pval))
  
  all_acrs_combined <- bind_rows(broadly_accessible_acrs, cell_type_specific_acrs, cell_type_restricted_acrs)
  write_tsv(all_acrs_combined, file = paste0(opt_ct_dir,'/',prefix, ".all_ACRs.classified.bed", collapse = "."), col_names = FALSE)
  
  return(all_acrs_combined)
}




## Read Inputs 
input <- input_data
message("Reading Peak File")
bed_file_read <- read_delim(peak_file, col_names = c("chrom", "start", "stop", "acr_number", "accessability"))

message("Reading Meta Data...")
meta_data <- read.delim(meta,row.names = 1)



##updating 120325
##the pattern is the
##heat:leaf1,leaf2;control:leaf3,leaf4
if (group_pattern != 'na'){
  
  group_map <- unlist(
    lapply(strsplit(group_pattern, ";")[[1]], function(x) {
      tmp <- strsplit(x, ":")[[1]]
      libs <- strsplit(tmp[2], ",")[[1]]
      setNames(rep(tmp[1], length(libs)), libs)
    })
  )
  
  meta_data$group <- group_map[meta_data[[target_treat_colnm]]]
  meta_slot <- 'group'
  message(paste0('the group pattern is ',group_pattern))
  
}else{

  ##target_treat_colnm
  meta_slot <- target_treat_colnm
  message('there is no group pattern')
}

##updating 120825
write.table(meta_data,paste0(input_output_dir,'/temp_add_group_meta.txt'),quote = F,sep = '\t')



##updating 121125
##do not run the group that contain 0 for either one of group
meta_data$count <- 1
aggregated_dt <- aggregate(count ~ cell_identity + group,data = meta_data,sum)
all_groups <- unique(aggregated_dt$group)
cell_ids <- unique(aggregated_dt$cell_identity)
missing_cell_types <- cell_ids[
  sapply(cell_ids, function(cid) {
    groups_present <- aggregated_dt$group[aggregated_dt$cell_identity == cid]
    length(setdiff(all_groups, groups_present)) > 0
  })
]

#message(paste0('the missing cell types is ',missing_cell_types))






##updating 051925
##we will split the meta based on treat group 
##use the for loop
all_celltype_all_list <- unique(meta_data[[target_cluster]])

##updating 121125 remove the missing cell types
all_celltype_list <- setdiff(all_celltype_all_list, missing_cell_types)

for (i in 1:length(all_celltype_list)) {
  
  ipt_target_celltype <- all_celltype_list[i]
  
  ##create an directory to run for each cell type
  opt_ct_dir <- paste0(input_output_dir,'/',ipt_target_celltype)
  if (!dir.exists(opt_ct_dir)) {
    dir.create(opt_ct_dir, recursive = TRUE)
    message("Directory created at: ", opt_ct_dir)
  } else {
    message("Directory already exists: ", opt_ct_dir)
  }
  
  
  ##updating 121125
  ##check if the opt_pvalues.csv exist
  if (file.exists(paste0(opt_ct_dir,'/','opt', ".all_pvalues.csv"))) {
    
    message("File found and we will re-run to test different threshold")
    
    final_prefix <- paste0(prefix,'.',threshold)
    
    calculated_Pvals_all_quantify <- read.csv(paste0(opt_ct_dir,'/','opt', ".all_pvalues.csv"))
    
    if (stat_test == "perm") {
      quantify_cell_type_specific_acrs(calculated_Pvals_all_quantify, "perm_pval", threshold, final_prefix)
    } else if (stat_test == "pnorm") {
      quantify_cell_type_specific_acrs(calculated_Pvals_all_quantify, "pnorm_pval", threshold, final_prefix)
    } else {
      stop("Invalid stat_test value")
    }
    
    message("Done! Check your outputs Dork! :D ")
    
    
    
  } else {
    
    final_prefix <- paste0(prefix,'.',threshold)
  
    meta_ct_data <- meta_data[meta_data[[target_cluster]] == ipt_target_celltype,]
    
    ## Use the column for meta_data for cell type ACR calling 
    meta_slot_var <- c(meta_slot)
    
    message("Reading Input Data...")
    raw_cpm_counts_all_genes <- read_delim(input, delim="\t", col_names = c("gene_name", "barcode", "accessability")) %>%
        dplyr::mutate(cellID = barcode)  %>%
        dplyr::mutate(geneID = gene_name)
    
    
    ## Generate the Null distribution of values... 
    null_distributions <- replicate(null_permutations, generate_null_distribution(meta_ct_data, meta_slot_var), simplify = FALSE)
    #Older slower implementation 
    #null_dist_values <- lapply(null_distributions, generate_null_dist_values_optimized, "final_annotation_n", raw_cpm_counts_all_genes)
    
    
    counter <- 1
    null_dist_values <- mclapply(null_distributions, generate_null_dist_values_optimized, 
             meta_slot_var, raw_cpm_counts_all_genes,
             mc.preschedule = FALSE, mc.set.seed = TRUE,
             mc.silent = FALSE, mc.cores = 25,
             mc.cleanup = TRUE, mc.allow.recursive = TRUE, affinity.list = NULL)
    
    
    # Apply the function to the list of matrices
    valid_indices <- sapply(null_dist_values, is_valid_matrix)
    # Filter out invalid matrices
    null_dist_values <- null_dist_values[valid_indices]
    
    
    
    
    message("Generating 1000 Bootstraps")
    y <- bootstraps(meta_ct_data, times = entropy_bootstraps, strata = !!sym(meta_slot_var))
    
    # Set up parallel processing
    message("Running Bootstraps...")
    
    ## Older Slower Version
    # counter <- 1
    # tic()
    # results <- y %>%
    #   mutate(p_values = map(splits, ~ generate_pvalues_bootstraps_fast(analysis(.x), raw_cpm_counts_all_genes, "final_annotation_n"))) %>%
    #   select(id, p_values)
    # toc()
    
    #Set to 3gb
    #options(future.globals.maxSize= 9e+9)
    ## Increase speed by threading...
    options(future.globals.maxSize= 100000 * 1024^2)
    
    ##debug
    ##debug
    #saveRDS(y,'temp_y.rds')
    #results <- y %>%
    #  mutate(safe_output = future_map(splits, ~ generate_pvalues_bootstraps_fast(analysis(.x), raw_cpm_counts_all_genes, meta_slot_var)))
    
    #saveRDS(results,'temp_results_1.rds')
    
    
    
    generate_pvalues_bootstraps_fast_safe_function <- purrr::safely(generate_pvalues_bootstraps_fast)
    counter <- 1
    plan(multicore)
    results <- y %>%
     mutate(safe_output = future_map(splits, ~ generate_pvalues_bootstraps_fast_safe_function(analysis(.x), raw_cpm_counts_all_genes, meta_slot_var))) %>%
     mutate(has_error = map_lgl(safe_output, ~ !is.null(.x$error))) %>%
     # Filter out rows with errors
     filter(!has_error) %>%
     # Extract p_values from the safe_output
     mutate(p_values = map(safe_output, "result")) %>%
     select(id, p_values)
    
    saveRDS(results,'temp_results.rds')
    
    expaneded_bootstraps <- results %>%
      mutate(p_values = map(p_values, ~ as_tibble(.x, rownames = "ACR_values"))) %>% # Convert matrices to tibbles and include row names
      unnest(cols = p_values) %>%                          # Unnest the tibbles
      pivot_longer(cols = -c(id, ACR_values),              # Keep 'id' and 'ACR_values' fixed
                   names_to = "cell_type",                 # Assign the column names to 'cell_type'
                   values_to = "value") %>%                # Assign the values to 'value'
      rename(BootstrapID = id)                             # Rename the columns to the desired names
    
    nested_data <- expaneded_bootstraps %>% 
        select(-BootstrapID) %>%
        group_by(ACR_values, cell_type) %>%
        nest() %>% 
        rename(distribution = "data")
    
    
    message("Merging Bootstraps and Null...")
    
    melted_nested_nulls <- imap_dfr(null_dist_values, function(matrix, index) {
      df_matrix <- reshape2::melt(matrix)
      df_matrix <- df_matrix %>%
        rename(
          row_name = Var1,
          column_name = Var2,
          value = value
        ) %>%
        mutate(matrix_index = index)
      return(df_matrix)
    })
    
    
    rm(null_dist_values)
    gc()
    
    null_dist_generation <- melted_nested_nulls %>% 
        rename(ACR_values = row_name) %>%
        select(ACR_values, value) %>%
        group_by(ACR_values) %>%
        nest() %>%
        rename(null_dist = data)
    
    saveRDS(null_dist_generation,'temp_null_dist_generation.rds')
    
    merged_bootstraps_nulls <- left_join(nested_data, null_dist_generation, by = c("ACR_values"))
    
    
    
    
    message("Calculating Pvalues....")
    calculated_Pvals_all <- calculate_pvalues_per_ACR(merged_bootstraps_nulls)
    ## Save the RDS object here after all the hard work is done... 
    
    message("Saving Intermediate Data")
    save_pval_null_dists <- paste0(opt_ct_dir,'/',final_prefix, ".combined_data.rds")
    saveRDS(calculated_Pvals_all, file = save_pval_null_dists)
    
    message("Filtering Pvalues...")
    calculated_Pvals_all_quantify <- calculated_Pvals_all %>% 
        dplyr::select(ACR_values, cell_type, perm_pval, perm_pval_lower, perm_pval_upper, pnorm_pval, z_score) 
    
    ##change the prefix to opt
    save_pval_null_dists <- paste0(opt_ct_dir,'/','opt', ".all_pvalues.csv")
    write_delim(calculated_Pvals_all_quantify, file = save_pval_null_dists, delim=",")
    
    
    ## Clean up happen later after script is finalized. Stops from StackOverflow
    rm(calculated_Pvals_all)
    rm(merged_bootstraps_nulls)
    gc()
    
    #setwd("/scratch/jpm73279/comparative_single_cell/dev_location/entropy_final")
    
    
    
    if (stat_test == "perm") {
      quantify_cell_type_specific_acrs(calculated_Pvals_all_quantify, "perm_pval", threshold, final_prefix)
    } else if (stat_test == "pnorm") {
      quantify_cell_type_specific_acrs(calculated_Pvals_all_quantify, "pnorm_pval", threshold, final_prefix)
    } else {
      stop("Invalid stat_test value")
    }
    
    message("Done! Check your outputs Dork! :D ")

  
  }
  
  
}
