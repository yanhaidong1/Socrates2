##this script is to build the quality control plots for the article
library(ggplot2)
library(dplyr)
library(viridis)
library(ggrepel)

##updating 022826 we will use the two sd to perform the cutoff

args <- commandArgs(T)

input_obj_fl <- as.character(args[1])

input_output_dir <- as.character(args[2])

input_prefix <- as.character(args[3])

ipt_config_fl <- as.character(args[4])

source(ipt_config_fl)
##it includes the min_effect_value

input_obj <- readRDS(input_obj_fl)

merged_dt <- input_obj$meta

merged_dt$FRIP <- merged_dt$acrs/merged_dt$total
merged_dt$TSS <- merged_dt$tss/merged_dt$total

merged_dt <- merged_dt[c('cellID','FRIP','TSS','doubletscore','LouvainClusters','library','umap1','umap2')]


##updating 112025
##plot replications
library_list <- unique(merged_dt$library)

library_colors <- setNames(
  scales::hue_pal()(length(library_list)), 
  library_list
)


for (i in 1:length(library_list)){
  
  lib_name <- library_list[i]
  
  p <- ggplot(merged_dt %>% filter(library==lib_name),
              aes(x=umap1, y=umap2, color=library)) +
    geom_point(alpha=0.7, size=1.5) +
    scale_color_manual(values = library_colors) +
    #scale_color_manual(values = setNames("skyblue", lib_name))+
    #scale_color_manual(values=c(lib_name="skyblue")) +
    labs(title=paste0("UMAP of ",lib_name)) +
    theme_minimal()+
    theme_classic()
  pdf(paste0(input_output_dir,'/','opt_',input_prefix,'.UMAP.library.',lib_name,'.pdf'),height = 8,width = 8)
  print(p)
  dev.off()
  
}



##############
##for the FRIP
df <- merged_dt
df <- df %>%
  mutate(FRIP_scaled = (FRIP - min(FRIP)) / (max(FRIP) - min(FRIP)))

p <- ggplot(df, aes(x = umap1, y = umap2, color = FRIP_scaled)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_viridis(option = "magma", direction = -1) +
  #scale_color_gradient(low = "blue", high = "red") +
  theme_classic() +
  labs(color = "FRIP (scaled)",
       x = "UMAP1",
       y = "UMAP2")

pdf(paste0(input_output_dir,'/','opt_',input_prefix,'.UMAP.FRiP.pdf'),width = 8,height = 8)
print(p)
dev.off()

##assign the df to the ipt_dt
ipt_dt <- df

test_each_cluster_strict <- function(df, sd_fold = 2) {
  
  clusters <- unique(df$LouvainClusters)
  
  #df$FRIP <- df$Peak.Tn5/df$Total.Tn5
  
  ## 1. 先计算每个 cluster 的均值
  cluster_means <- tapply(df$FRIP,
                          df$LouvainClusters,
                          mean,
                          na.rm = TRUE)
  
  ## 2. cluster 均值之间的 SD
  cluster_sd <- sd(cluster_means, na.rm = TRUE)
  
  results <- data.frame(
    LouvainClusters = clusters,
    p_value = NA_real_,
    mean_FRIP = NA_real_,
    delta = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (cl in clusters) {
    x <- df$FRIP[df$LouvainClusters == cl]
    y <- df$FRIP[df$LouvainClusters != cl]
    
    res <- wilcox.test(x, y, alternative = "two.sided")
    
    mean_x <- mean(x, na.rm = TRUE)
    mean_y <- mean(y, na.rm = TRUE)
    
    ##we will use the global - cluster
    delta <- mean_y - mean_x
    #delta <- mean_x - mean_y
    
    results[results$LouvainClusters == cl,
            c("p_value", "mean_FRIP", "delta")] <-
      c(res$p.value, mean_x, delta)
  }
  
  ## 3. 严格判定：基于 cluster mean SD
  results$significant <- results$p_value < 0.05 &
    results$delta > sd_fold * cluster_sd
  
  ## 可选：保留阈值，方便画图和解释
  results$cluster_sd <- cluster_sd
  results$threshold  <- sd_fold * cluster_sd
  
  results
}

##for the original without filtering the xylem with high doublet scores
res_table <- test_each_cluster_strict(ipt_dt, sd_fold = sd_fold_value)
#default 2

res_table$mlog10pvalue <- -log10(res_table$p_value)

##replace the Inf to the biggest value
res_table[] <- lapply(res_table, function(x) {
  if (is.numeric(x)) {
    max_val <- max(x[is.finite(x)], na.rm = TRUE)
    x[is.infinite(x)] <- max_val
  }
  x
})

##set the range 
x_range <- range(res_table$mlog10pvalue, na.rm = TRUE)
y_range <- range(res_table$delta, na.rm = TRUE)

cutoff_threshold <-unique(res_table$threshold)

p <- ggplot(res_table, aes(x = mlog10pvalue, y = delta)) +
  geom_point(aes(color = significant), size = 8, alpha = 1) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red")) +
  
  geom_hline(yintercept = cutoff_threshold,
             linetype = "dotted", color = "blue", linewidth = 0.7) +
  geom_vline(xintercept = -log10(0.05),
             linetype = "dotted", color = "darkgreen", linewidth = 0.7) +
  
  
  geom_text_repel(
    aes(label = LouvainClusters),
    size = 6,
    max.overlaps = Inf,      # 确保每个点都有标签
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey70",
    segment.size = 0.3
  ) +
  coord_cartesian(
    xlim = x_range + c(-0.5, 0.5) * diff(x_range),
    ylim = y_range + c(-0.5, 0.5) * diff(y_range)
  ) +
  
  labs(
    x = expression(-log[10](pvalue)),
    y = expression(Delta~"(Global mean - Cluster mean FRIP)"),
    color = "Significant",
    title = "Cluster-level FRIP Changes"
  ) +
  theme_classic(base_size = 14)

pdf(paste0(input_output_dir,'/','opt_',input_prefix,'.UMAP.FRiP.cluster.test.pdf'),width = 8,height = 6 )
p
dev.off()



#############
##for the TSS
df <- merged_dt
df <- df %>%
  mutate(FRIP_scaled = (TSS - min(TSS)) / (max(TSS) - min(TSS)))

p <- ggplot(df, aes(x = umap1, y = umap2, color = FRIP_scaled)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_viridis(option = "viridis", direction = -1) +
  #scale_color_gradient(low = "blue", high = "red") +
  theme_classic() +
  labs(color = "TSS (scaled)",
       x = "UMAP1",
       y = "UMAP2")

pdf(paste0(input_output_dir,'/','opt_',input_prefix,'.UMAP.TSS.pdf'),width = 8,height = 8)
print(p)
dev.off()

##assign the df to the ipt_dt
ipt_dt <- df

test_each_cluster_strict <- function(df, sd_fold = 2) {
  
  clusters <- unique(df$LouvainClusters)
  
  #df$TSS <- df$TSS.Tn5/df$Total.Tn5
  
  ## 1. 先计算每个 cluster 的均值
  cluster_means <- tapply(df$TSS,
                          df$LouvainClusters,
                          mean,
                          na.rm = TRUE)
  
  ## 2. cluster 均值之间的 SD
  cluster_sd <- sd(cluster_means, na.rm = TRUE)
  
  results <- data.frame(
    LouvainClusters = clusters,
    p_value = NA_real_,
    mean_FRIP = NA_real_,
    delta = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (cl in clusters) {
    x <- df$TSS[df$LouvainClusters == cl]
    y <- df$TSS[df$LouvainClusters != cl]
    
    res <- wilcox.test(x, y, alternative = "two.sided")
    
    mean_x <- mean(x, na.rm = TRUE)
    mean_y <- mean(y, na.rm = TRUE)
    
    ##we will use the global - cluster
    delta <- mean_y - mean_x
    #delta <- mean_x - mean_y
    
    results[results$LouvainClusters == cl,
            c("p_value", "mean_TSS", "delta")] <-
      c(res$p.value, mean_x, delta)
  }
  
  ## 3. 严格判定：基于 cluster mean SD
  results$significant <- results$p_value < 0.05 &
    results$delta > sd_fold * cluster_sd
  
  ## 可选：保留阈值，方便画图和解释
  results$cluster_sd <- cluster_sd
  results$threshold  <- sd_fold * cluster_sd
  
  results
}

##for the original without filtering the xylem with high doublet scores
res_table <- test_each_cluster_strict(ipt_dt, sd_fold = sd_fold_value)

res_table$mlog10pvalue <- -log10(res_table$p_value)

##replace the Inf to the biggest value
res_table[] <- lapply(res_table, function(x) {
  if (is.numeric(x)) {
    max_val <- max(x[is.finite(x)], na.rm = TRUE)
    x[is.infinite(x)] <- max_val
  }
  x
})

##set the range 
x_range <- range(res_table$mlog10pvalue, na.rm = TRUE)
y_range <- range(res_table$delta, na.rm = TRUE)

cutoff_threshold <-unique(res_table$threshold)

p <- ggplot(res_table, aes(x = mlog10pvalue, y = delta)) +
  geom_point(aes(color = significant), size = 8, alpha = 1) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red")) +
  
  geom_hline(yintercept = cutoff_threshold,
             linetype = "dotted", color = "blue", linewidth = 0.7) +
  geom_vline(xintercept = -log10(0.05),
             linetype = "dotted", color = "darkgreen", linewidth = 0.7) +
  
  
  geom_text_repel(
    aes(label = LouvainClusters),
    size = 6,
    max.overlaps = Inf,      # 确保每个点都有标签
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey70",
    segment.size = 0.3
  ) +
  coord_cartesian(
    xlim = x_range + c(-0.5, 0.5) * diff(x_range),
    ylim = y_range + c(-0.5, 0.5) * diff(y_range)
  ) +
  
  labs(
    x = expression(-log[10](pvalue)),
    y = expression(Delta~"(Global mean - Cluster mean TSS)"),
    color = "Significant",
    title = "Cluster-level TSS Changes"
  ) +
  theme_classic(base_size = 14)

pdf(paste0(input_output_dir,'/','opt_',input_prefix,'.UMAP.TSS.cluster.test.pdf'),width = 8,height = 6 )
p
dev.off()




#######################
##for the doublet score
df <- merged_dt
df <- df %>%
  mutate(doubletscore_scaled = (doubletscore - min(doubletscore)) / (max(doubletscore) - min(doubletscore)))

p <- ggplot(df, aes(x = umap1, y = umap2, color = doubletscore_scaled)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_viridis(option = "turbo", direction = -1) +
  #scale_color_gradient(low = "blue", high = "red") +
  theme_classic() +
  labs(color = "doubletscore (scaled)",
       x = "UMAP1",
       y = "UMAP2")

pdf(paste0(input_output_dir,'/','opt_',input_prefix,'.UMAP.doublet.pdf'),width = 8,height = 8)
print(p)
dev.off()

##assign the df to the ipt_dt
ipt_dt <- df

test_each_cluster_strict <- function(df, sd_fold = 2) {
  
  clusters <- unique(df$LouvainClusters)
  
  ## 1. 先计算每个 cluster 的均值
  cluster_means <- tapply(df$doubletscore,
                          df$LouvainClusters,
                          mean,
                          na.rm = TRUE)
  
  ## 2. cluster 均值之间的 SD
  cluster_sd <- sd(cluster_means, na.rm = TRUE)
  
  results <- data.frame(
    LouvainClusters = clusters,
    p_value = NA_real_,
    mean_FRIP = NA_real_,
    delta = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (cl in clusters) {
    x <- df$doubletscore[df$LouvainClusters == cl]
    y <- df$doubletscore[df$LouvainClusters != cl]
    
    res <- wilcox.test(x, y, alternative = "two.sided")
    
    mean_x <- mean(x, na.rm = TRUE)
    mean_y <- mean(y, na.rm = TRUE)
    
    ##we will use the cluster - global
    #delta <- mean_y - mean_x
    delta <- mean_x - mean_y
    
    results[results$LouvainClusters == cl,
            c("p_value", "mean_FRIP", "delta")] <-
      c(res$p.value, mean_x, delta)
  }
  
  ## 3. 严格判定：基于 cluster mean SD
  results$significant <- results$p_value < 0.05 &
    results$delta > sd_fold * cluster_sd
  
  ## 可选：保留阈值，方便画图和解释
  results$cluster_sd <- cluster_sd
  results$threshold  <- sd_fold * cluster_sd
  
  results
}

##for the original without filtering the xylem with high doublet scores
res_table <- test_each_cluster_strict(ipt_dt, sd_fold = sd_fold_value)

res_table$mlog10pvalue <- -log10(res_table$p_value)

##replace the Inf to the biggest value
res_table[] <- lapply(res_table, function(x) {
  if (is.numeric(x)) {
    max_val <- max(x[is.finite(x)], na.rm = TRUE)
    x[is.infinite(x)] <- max_val
  }
  x
})

##set the range 
x_range <- range(res_table$mlog10pvalue, na.rm = TRUE)
y_range <- range(res_table$delta, na.rm = TRUE)

cutoff_threshold <-unique(res_table$threshold)

p <- ggplot(res_table, aes(x = mlog10pvalue, y = delta)) +
  geom_point(aes(color = significant), size = 8, alpha = 1) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red")) +
  
  geom_hline(yintercept = cutoff_threshold,
             linetype = "dotted", color = "blue", linewidth = 0.7) +
  geom_vline(xintercept = -log10(0.05),
             linetype = "dotted", color = "darkgreen", linewidth = 0.7) +
  
  
  geom_text_repel(
    aes(label = LouvainClusters),
    size = 6,
    max.overlaps = Inf,      # 确保每个点都有标签
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey70",
    segment.size = 0.3
  )+
  
  coord_cartesian(
    xlim = x_range + c(-0.5, 0.5) * diff(x_range),
    ylim = y_range + c(-0.5, 0.5) * diff(y_range)
  ) +
  
  labs(
    x = expression(-log[10](pvalue)),
    y = expression(Delta~"(Cluster mean - Global mean FRIP)"),
    color = "Significant",
    title = "Cluster-level Doublet Score Changes"
  ) +
  theme_classic(base_size = 14)

pdf(paste0(input_output_dir,'/','opt_',input_prefix,'.UMAP.doublet.cluster.test.pdf'),width = 8,height = 6 )
p
dev.off()

























