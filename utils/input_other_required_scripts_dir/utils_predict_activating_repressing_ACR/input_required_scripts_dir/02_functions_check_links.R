# #---function to calculate peak center to gene distance.
# cal_dis <- function(x.m){
#   dis <- ifelse(x.m$gene_strand == "-", ifelse(x.m$peak_center < x.m$gene_end, x.m$gene_end - x.m$peak_center,
#                                                ifelse(x.m$peak_center > x.m$gene_start, x.m$peak_center - x.m$gene_start,0)),
#                 ifelse(x.m$peak_center < x.m$gene_start, x.m$gene_start - x.m$peak_center,
#                        ifelse(x.m$peak_center > x.m$gene_end, x.m$peak_center-x.m$gene_end,0)))
#   return(dis)
# }

#---function initialize link gene matrix-----
#link.m have columns: peak_id, gene_id
#gene.bed have columns: "gene_chr","gene_start","gene_end", "gene_id", "gene_strand"
add_distance <- function(link.m, gene.bed){
  link.m <- left_join(link.m, gene.bed, by = "gene_id")
  
  link.m$gene_start <- as.numeric(link.m$gene_start)
  link.m$gene_end <- as.numeric(link.m$gene_end)
  
  link.m$peak_chr <- do.call(rbind,strsplit(link.m$peak_id,split = "_"))[,1]
  link.m$peak_start <- as.numeric(do.call(rbind,strsplit(link.m$peak_id,split = "_"))[,2])
  link.m$peak_end <- as.numeric(do.call(rbind,strsplit(link.m$peak_id,split = "_"))[,3])
  link.m$peak_center <- as.numeric(link.m$peak_start) + 250
  
  # link.m$distance <- ifelse(link.m$gene_strand == "-", ifelse(link.m$peak_center < link.m$gene_start, link.m$gene_start - link.m$peak_center,
  #                                              ifelse(link.m$peak_center > link.m$gene_end, link.m$peak_center - link.m$gene_start,0)),
  #               ifelse(link.m$peak_center < link.m$gene_start, link.m$gene_start - link.m$peak_center,
  #                      ifelse(link.m$peak_center > link.m$gene_end, link.m$peak_center-link.m$gene_end,0)))
  link.m$distance <- ifelse(link.m$peak_center < link.m$gene_start, link.m$gene_start - link.m$peak_center,
                              ifelse(link.m$peak_center > link.m$gene_end, link.m$peak_center-link.m$gene_end,0))
  link.m$distance.type <- ifelse(link.m$distance < 260, "gLink", ifelse(link.m$distance < 2000, "pLink", "dLink"))
  return(link.m)
}

#Function to filter distal link with better proximal link.
filter_distal_links <- function(link.m, min_dis = 2000){
  peaks <- unique(link.m$peak_id)
  link.f <- lapply(peaks, function(x){
    df <- link.m[link.m$peak_id == x,]
    if (max(df$distance) <= min_dis || min(df$distance) > min_dis)
    {
      df$Filter_distal <- "keep"
      return(df)
    }else{
      df.p <- df[df$distance <= min_dis,]
      df.d <- df[df$distance > min_dis,]
      df.p$Filter_distal <- "keep"
      df.d$Filter_distal <- ifelse(abs(df.d$cor.value) < max(abs(df.p$cor.value)), "remove", "keep")
      df <- rbind(df.p, df.d)
      return(df)
    }
  })
  link.f <- do.call(rbind, link.f)
  return(link.f)
}


#Function to filter distal link with better proximal link.
#20240405 
#1. filter link with minor correlation scores if only have distal or proximal link.
#2. when proximal and distal, only keep the distal if cor.distal > cor.proximal
#
filter_distal_links2 <- function(link.m, min_dis = 2000){
  peaks <- unique(link.m$peak_id)
  link.f <- lapply(peaks, function(x){
    df <- link.m[link.m$peak_id == x,]
    if (max(df$distance) <= min_dis || min(df$distance) > min_dis) 
    {
	  df$Filter_distal <- ifelse(abs(df$cor.value) == max(abs(df$cor.value)), "keep", "remove")
      return(df)
    }else{
      df.p <- df[df$distance <= min_dis,]
      df.d <- df[df$distance > min_dis,]
      df.p$Filter_distal <- "keep"
      df.d$Filter_distal <- ifelse(abs(df.d$cor.value) < max(abs(df.p$cor.value)), "remove", "keep")
      df <- rbind(df.p, df.d)
      return(df)
    }
  })
  link.f <- do.call(rbind, link.f)
  return(link.f)
}




#Function to remove all skip link if exist nearby_link
filter_skip_links <- function(link.m, keep_high_skipLink = "FALSE"){
  peaks <- unique(link.m$peak_id)
  link.f <- lapply(peaks, function(x){
    df <- link.m[link.m$peak_id == x,]
    if (sum(df$link.type == "nearby_link") == 0)
    {
      df$Filter_skipLink <- "keep"
      return(df)
    }else if(keep_high_skipLink == "FALSE"){
	  df$Filter_skipLink <- ifelse(df$link.type == "nearby_link","keep", "remove")
      return(df)
    }else{
	  df.n <- df[df$link.type == "nearby_link",]
      df.s <- df[df$link.type == "skip_link",]
      df.n$Filter_skipLink <- "keep"
      df.s$Filter_skipLink <- ifelse(abs(df.s$cor.value) < max(abs(df.n$cor.value)), "remove", "keep")
      df <- rbind(df.n, df.s)
      return(df)
	}
  })
  link.f <- do.call(rbind, link.f)
  return(link.f)
}

#Function to filter distal link with better proximal link.
filter_distal_links1 <- function(link.m, min_dis = 2000){
  peaks <- unique(link.m$peak_id)
  link.f <- lapply(peaks, function(x){
    df <- link.m[link.m$peak_id == x,]
    if (nrow(df) == 1)
    {
      df$Filter_distal <- "keep"
      return(df)
    }else{
	  df$Filter_distal <- ifelse(abs(df$cor.value) == max(abs(df$cor.value)), "keep","remove") #only keep the top core value.
      return(df)
    }
  })
  link.f <- do.call(rbind, link.f)
  return(link.f)
}



#Function to check cre function.
check_cre_type <- function(link.m){
  peaks <- unique(link.m$peak_id)
  link.f <- lapply(peaks, function(x){
    df <- link.m[link.m$peak_id == x,]
    df$cre.type <- ifelse(max(df$cor.value) < 0, "silencer", ifelse(min(df$cor.value) > 0, "enhancer", "bi_function"))
    return(df)
  })
  link.f <- do.call(rbind, link.f)
  return(link.f)
}


summary_link <- function(output_dir,link.m, prefix){
  #---correlation density plot---
  if (length(unique(link.m$distance.type)) >1){
  pdf(
    file = paste0(output_dir,'/',prefix, "_corelation_density.pdf"),
    width = 4,
    height = 4
  )
  den <-
    sm.density.compare(
      link.m$cor.value,
      group = link.m$distance.type,
      model = "equal",
      lwd = 2
    )
  legend(
    "topright",
    den$levels,
    col = den$col,
    lty = den$lty,
    lwd = den$lwd
  )
  abline(
    v = -0.4,
    col = "red",
    lwd = 2,
    lty = 2
  )
  abline(
    v = 0.4,
    col = "red",
    lwd = 2,
    lty = 2
  )
  title(main = prefix,cex.main = 0.6)
  dev.off()
  }else{
  pdf(
    file = paste0(output_dir,'/',prefix, "_corelation_density.pdf"),
    width = 4,
    height = 4
  )
  
  plot(density(link.m$cor.value), lwd = 2)
  abline(
    v = -0.4,
    col = "red",
    lwd = 2,
    lty = 2
  )
  abline(
    v = 0.4,
    col = "red",
    lwd = 2,
    lty = 2
  )
  title(main = prefix,cex.main = 0.6)
  dev.off()
  }
  #---distance density plot---
  # #plot density one by one
  # plot(density(gfg))
  # lines(density(a), col = "red")
  # lines(density(b), col = "green")
  #
  # legend("topright", c("gfg", "a", "b"),
  #        col =c("black","red","green"), lty=1)
  
  #ref modify x-axis: https://stackoverflow.com/questions/25997337/in-r-how-to-set-the-breaks-of-x-axis
  if(length(unique(link.m$cre.type)) > 1)
  {
  pdf(
    file = paste0(output_dir,'/',prefix, "_distance_density.pdf"),
    width = 4,
    height = 4
  )
  den <-
    sm.density.compare(
      log10(link.m$distance + 1),
      group = link.m$cre.type,
      model = "equal",
      lwd = 2
    )
  legend(
    "topright",
    den$levels,
    col = den$col,
    lty = den$lty,
    lwd = den$lwd
  )
  abline(
    v = 2,
    col = "red",
    lwd = 2,
    lty = 2
  )
  title(main = prefix, cex.main = 0.6)
  dev.off()
  }else {
  pdf(
    file = paste0(output_dir,'/',prefix, "_distance_density.pdf"),
    width = 4,
    height = 4
  )
  plot(density(log10(link.m$distance + 1)), lwd = 2)
  abline(
    v = 2,
    col = "red",
    lwd = 2,
    lty = 2
  )
  title(main = prefix, cex.main = 0.6)
  dev.off()
  }
  
  #check silencer enhancer and bi-funcitonal ACRs.
  acr.s <- link.m[link.m$cor.value < 0, ]$peak_id
  acr.e <- link.m[link.m$cor.value > 0, ]$peak_id
  pdf(
    file = paste0(output_dir,'/',prefix, "_Link_type.pdf"),
    width = 6,
    height = 6
  )
  myplot <- ggvenn(
    list(silencer = acr.s, enhancer = acr.e),
    fill_color = c("#0073C2FF", "#EFC000FF"),
    stroke_size = 1,
    set_name_size = 8,
    text_size = 8
  )
  print(myplot)
  dev.off()
  
  #---check acr location type proportion in each cre type---
  #plot type2
  link.m.unque <- link.m %>% distinct(peak_id, .keep_all = TRUE)
  link_type.m <- table(link.m.unque$cre.type, link.m.unque$peak.type2)
  link_type.m <- t(link_type.m)
  link_type.m.r <- link_type.m / rowSums(link_type.m)
  
  pdf(
    file = paste0(output_dir,'/',prefix, "_CRE_type.pdf"),
    width = 6,
    height = 6
  )
  barplot(
    link_type.m.r,
    border = "black",
    col = colorRampPalette(brewer.pal(12, "Paired"))(nrow(link_type.m.r)),
    font.axis = 2,
    beside = T,
    legend = rownames(link_type.m.r),
    font.lab = 2,
    main = prefix,
    cex.main= 0.6
  )
  dev.off()
  
  #---check acr2gene pair type proportion in each cre type---
  #plot type2
  #link.m.unque <- link.m %>% distinct(peak_id, .keep_all = TRUE)
  link_type.m <- table(link.m$link.type, link.m$peak.type2)
  link_type.m <- t(link_type.m)
  link_type.m.r <- link_type.m / rowSums(link_type.m)
  
  pdf(
    file = paste0(output_dir,'/',prefix, "_ACR2gene_pair_type.pdf"),
    width = 6,
    height = 6
  )
  barplot(
    link_type.m.r,
    border = "black",
    col = colorRampPalette(brewer.pal(12, "Paired"))(nrow(link_type.m.r)),
    font.axis = 2,
    beside = T,
    legend = rownames(link_type.m.r),
    font.lab = 2,
    main = prefix,
    cex.main= 0.6
  )
  dev.off()
  
  write.table(link_type.m, file = paste0(output_dir,'/',prefix, "_ACR2gene_pair_type.txt"))
  write.table(link_type.m.r, file = paste0(output_dir,'/',prefix, "_ACR2gene_pair_type_rate.txt"))
}

#Function output the ACRs bed files.
write_bed <- function(df, prefix){
  lapply(unique(df$cre.type), function(x){
    peaks <- unique(df[df$cre.type == x,]$peak_id)
    peaks.m <- as.data.frame(do.call(rbind,(strsplit(peaks, split = "_"))))
    write.table(peaks.m, file = paste0(prefix,"_",x,"_peaks.bed"), quote = F, sep = "\t", col.names = F, row.names = F)
  })
}

#---Motif enrichment function--------------
motif.enricher <- function(peak, ctrl, motif.m){
  
  #peak <- k27.c1.acr
  #ctrl <- c4.acr
  
  peak.motif <- subset(motif.m, row.names(motif.m) %in% peak)
  ctrl.motif <- subset(motif.m, row.names(motif.m) %in% ctrl)
  
  #  peak.motif <- motif.m[peak,]
  #  ctrl.motif <- motif.m[ctrl,]
  
  
  peak.sum <- as.data.frame(colSums(peak.motif))
  colnames(peak.sum) <- "peak.motif"
  peak.sum$peak.non <- c(length(row.names(peak.motif)) - peak.sum$peak.motif)
  
  ctrl.sum <- as.data.frame(colSums(ctrl.motif))
  colnames(ctrl.sum) <- "ctrl.motif"
  ctrl.sum$ctrl.non <- c(length(row.names(ctrl.motif)) - ctrl.sum$ctrl.motif)
  
  motif.all <- cbind(peak.sum, ctrl.sum)
  
  motif.all$p_value <- apply(motif.all, 1, function(x){
    fisher.m <- matrix(0,nrow = 2,ncol = 2)
    fisher.m[1,] <- c(x[1],x[2])
    fisher.m[2,] <- c(x[3],x[4])
    
    if(all(fisher.m == 0))
      p <- 1
    else
      p <- fisher.test(fisher.m, alternative='greater')$p.value
    return(p)
  })
  
  motif.all$p_adjusted <- p.adjust(motif.all$p_value)
  return(motif.all)
}
