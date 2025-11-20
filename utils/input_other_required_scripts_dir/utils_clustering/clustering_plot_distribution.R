##this script is to build the quality control plots for the article
library(ggplot2)
library(dplyr)
library(viridis)

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


#p <- ggplot(merged_dt %>% filter(library=="semroot6"),
#            aes(x=umap1, y=umap2, color=library)) +
#  geom_point(alpha=0.7, size=1.5) +
#  scale_color_manual(values=c("semroot6"="skyblue")) +
#  labs(title="UMAP of semroot6") +
#  theme_minimal()+
#  theme_classic()
#pdf('opt_library_semroot6_balanced_UMAP.pdf',height = 8,width = 8)
#p
#dev.off()

# semroot7
#p <- ggplot(merged_dt %>% filter(library=="semroot7"),
#            aes(x=umap1, y=umap2, color=library)) +
# geom_point(alpha=0.7, size=1.5) +
#  scale_color_manual(values=c("semroot7"="salmon")) +
#  labs(title="UMAP of semroot7") +
#  theme_minimal()+
#  theme_classic()
#pdf('opt_library_semroot7_balanced_UMAP.pdf',height = 8,width = 8)
#p
#dev.off()



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

##make a test for the sig test
##min_effect is the difference between the target cluster to average of other clusters 
test_each_cluster_strict <- function(df, min_effect = 0.1) {
  clusters <- unique(df$LouvainClusters)
  results <- data.frame(LouvainClusters = clusters, 
                        p_value = NA, 
                        mean_FRIP = NA, 
                        delta = NA)
  
  for (cl in clusters) {
    x <- df$FRIP[df$LouvainClusters == cl]
    y <- df$FRIP[df$LouvainClusters != cl]
    res <- wilcox.test(x, y, alternative = "less")
    
    mean_x <- mean(x)
    mean_y <- mean(y)
    delta  <- mean_y - mean_x  # how much lower the cluster is than others
    
    results[results$LouvainClusters == cl, ] <- c(cl, res$p.value, mean_x, delta)
  }
  
  results$p_value <- as.numeric(results$p_value)
  results$mean_FRIP <- as.numeric(results$mean_FRIP)
  results$delta <- as.numeric(results$delta)
  
  results$significant <- results$p_value < 0.05 & results$delta >= min_effect
  results
}

#for the df with fake
res_table <- test_each_cluster_strict(df, min_effect = min_effect_value_FRiP)
res_table$mlog10pvalue <- -log10(res_table$p_value)
p <- ggplot(res_table, aes(x = mlog10pvalue, y = delta)) +
  geom_point(aes(color = significant), size = 8, alpha = 1) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red")) +
  
  # Add horizontal and vertical reference lines
  #geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = min_effect_value_FRiP, linetype = "dotted", color = "blue", linewidth = 0.7) +
  geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "darkgreen", linewidth = 0.7) +
  coord_cartesian(xlim = c(0, 150), ylim = c(-0.15, 0.15)) +
  # Add cluster labels
  #geom_text(aes(label = LouvainClusters), vjust = -1, size = 3.2) +
  
  labs(
    x = expression(-log[10](pvalue)),
    y = expression(Delta~"(Cluster mean - Global mean FRIP)"),
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

##perform the testing
test_each_cluster_strict <- function(df, min_effect = 0.1) {
  clusters <- unique(df$LouvainClusters)
  results <- data.frame(LouvainClusters = clusters, 
                        p_value = NA, 
                        mean_FRIP = NA, 
                        delta = NA)
  
  for (cl in clusters) {
    x <- df$TSS[df$LouvainClusters == cl]
    y <- df$TSS[df$LouvainClusters != cl]
    res <- wilcox.test(x, y, alternative = "less")
    
    mean_x <- mean(x)
    mean_y <- mean(y)
    delta  <- mean_y - mean_x  # how much lower the cluster is than others
    
    results[results$LouvainClusters == cl, ] <- c(cl, res$p.value, mean_x, delta)
  }
  
  results$p_value <- as.numeric(results$p_value)
  results$mean_FRIP <- as.numeric(results$mean_FRIP)
  results$delta <- as.numeric(results$delta)
  
  results$significant <- results$p_value < 0.05 & results$delta >= min_effect
  results
}

#for the df with fake
res_table <- test_each_cluster_strict(df, min_effect = min_effect_value_TSS)
res_table$mlog10pvalue <- -log10(res_table$p_value)
p <- ggplot(res_table, aes(x = mlog10pvalue, y = delta)) +
  geom_point(aes(color = significant), size = 8, alpha = 1) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red")) +
  
  # Add horizontal and vertical reference lines
  #geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = min_effect_value_TSS, linetype = "dotted", color = "blue", linewidth = 0.7) +
  geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "darkgreen", linewidth = 0.7) +
  coord_cartesian(xlim = c(0, 150), ylim = c(-0.15, 0.15)) +
  # Add cluster labels
  #geom_text(aes(label = LouvainClusters), vjust = -1, size = 3.2) +
  
  labs(
    x = expression(-log[10](pvalue)),
    y = expression(Delta~"(Cluster mean - Global mean FRIP)"),
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


test_each_cluster_strict <- function(df, max_effect = -3) {
  clusters <- unique(df$LouvainClusters)
  results <- data.frame(LouvainClusters = clusters, 
                        p_value = NA, 
                        mean_FRIP = NA, 
                        delta = NA)
  
  for (cl in clusters) {
    x <- df$doubletscore[df$LouvainClusters == cl]
    y <- df$doubletscore[df$LouvainClusters != cl]
    res <- wilcox.test(x, y, alternative = "greater")
    
    mean_x <- mean(x)
    mean_y <- mean(y)
    delta  <- mean_y - mean_x  # how much lower the cluster is than others
    
    results[results$LouvainClusters == cl, ] <- c(cl, res$p.value, mean_x, delta)
  }
  
  results$p_value <- as.numeric(results$p_value)
  results$mean_FRIP <- as.numeric(results$mean_FRIP)
  results$delta <- as.numeric(results$delta)
  
  results$significant <- results$p_value < 0.05 & results$delta <= max_effect
  results
}

#for the df with fake
res_table <- test_each_cluster_strict(df, max_effect = max_effect_value_doublet)
res_table$mlog10pvalue <- -log10(res_table$p_value)
p <- ggplot(res_table, aes(x = mlog10pvalue, y = delta)) +
  geom_point(aes(color = significant), size = 8, alpha = 1) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red")) +
  
  # Add horizontal and vertical reference lines
  #geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = max_effect_value_doublet, linetype = "dotted", color = "blue", linewidth = 0.7) +
  geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "darkgreen", linewidth = 0.7) +
  coord_cartesian(xlim = c(0, 150), ylim = c(-10, 10)) +
  # Add cluster labels
  #geom_text(aes(label = LouvainClusters), vjust = -1, size = 3.2) +
  
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














