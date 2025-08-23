##load the package
library(Matrix)
#library(mclust)
library(ggplot2)
#library(parallel)
library(RColorBrewer)
#library(pheatmap)
#library(dplyr)
#library(VennDiagram)
#install.packages("ggalluvial")
#library(ggalluvial)
library(dplyr)
library(tidyr)
library(rlang)


args <- commandArgs(T)

ipt_syntenic_fl <- as.character(args[1])
output_dir <- as.character(args[2])
spe1 <- as.character(args[3])
spe2 <- as.character(args[4])

########
##step01 build the bar plot for different categories
ipt_dt <- read.delim(ipt_syntenic_fl)
head(ipt_dt)
ipt_dt$all_cate = paste0(ipt_dt[[paste0(spe1,'_ACRcate1')]],'__',ipt_dt[[paste0(spe1,'_ACRcate2')]])
table(ipt_dt$all_cate)
ipt_dt <- ipt_dt[ipt_dt$all_cate != 'na__na',]
ipt_dt$acrpair <- paste0(ipt_dt[[paste0(spe1,'_ACR')]],'__',ipt_dt[[paste0(spe2,'_ACR')]])

ipt_dt <- ipt_dt[c(paste0(spe1,'_ACRcate1'),paste0(spe1,'_ACRcate2'),'all_cate','acrpair',paste0(spe1,'_celltype'),paste0(spe2,'_celltype'))]
dim(ipt_dt)
ipt_dt <- unique(ipt_dt)
ipt_dt$count <- 1

head(ipt_dt)

fml <- as.formula(paste("count ~", paste0(spe1,'_ACRcate1'), "+", paste0(spe1,'_ACRcate2')))
fml2 <- as.formula(paste("count ~", paste0(spe1,'_ACRcate1')))

aggregate_two_cate_dt <- aggregate(fml , data = ipt_dt, sum)
aggregate_one_cate_dt <- aggregate(fml2, data = ipt_dt, sum)

###########
##pie chart
myPalette <- brewer.pal(10, "Set3") 
generate_pie <- function(dt) {
  #df <- read.delim(dt,header = FALSE)
  colnames(dt) <- c('fam','number')
  pct <- round(dt$number/sum(dt$number)*100)
  lbls <- paste(dt$fam,pct)
  lbls <- paste(lbls,"%",sep="") 
  #p_pie <- pie(dt$number,labels = lbls,col=myPalette,cex=4,border="white",edges = 200, radius = 1) ##25 13
  par(mar=c(3,18,3,15)) 
  return(pie(dt$number,labels = lbls,col=myPalette,cex=2,border="white",edges = 200, radius = 1) )
}
pdf(paste0(output_dir,'/opt_synACR_ratio_piechart_',spe1,'_',spe2,'.pdf'),width = 10,height = 10)
generate_pie(aggregate_one_cate_dt)
dev.off()

##########
##bar plot
colnames(aggregate_two_cate_dt) <- c('syn_ACR_category','ACR_category','count')
aggregate_two_cate_dt$syn_ACR_category <- gsub("speciesspec_ACR", "species_specific_ACR", aggregate_two_cate_dt$syn_ACR_category)

p <- ggplot(aggregate_two_cate_dt, aes(x=ACR_category,y=count,fill= syn_ACR_category)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  #scale_fill_manual(values=color_value) + 
  #geom_violin(trim=FALSE) + 
  #geom_boxplot(position = position_dodge2(preserve = "single")) + 
  #geom_boxplot(width=0.2,fill=organct) +
  #geom_jitter(shape=1, position=position_jitter(0.2)) +
  labs(x="\nCateogry", y = 'ACR count\n')+
  theme(
    plot.title = element_text(face="bold.italic",size=35,hjust = 0.5),
    axis.title = element_text(size =35),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=30,colour = "black"),
    
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) +
  #scale_x_discrete(drop=TRUE) +
  facet_wrap(~syn_ACR_category,drop=TRUE)
pdf(paste0(output_dir,'/opt_synACR_count_barplot_',spe1,'_',spe2,'.pdf'),width = 15,height = 7)
p
dev.off()


#################
##updating 082225
##use the bar plot to show
##for each cell type seperated by ',' in the spe1 cell type to show the count of all cell type in spe2 without sperating by ','
ipt_dt_seperate <- ipt_dt %>%
  separate_rows(!!sym(paste0(spe1,'_celltype')), sep = ",")

##########
##bar plot
ipt_dt <- ipt_dt_seperate
ipt_sharedACR_dt <- ipt_dt[ipt_dt[[paste0(spe1,'_ACRcate1')]] == 'shared_ACR',]

ipt_sharedACR_dt$twocatect <- paste0(ipt_sharedACR_dt[[paste0(spe1,'_celltype')]],'__',ipt_sharedACR_dt[[paste0(spe2,'_celltype')]])
ipt_sharedACR_dt <- ipt_sharedACR_dt[ipt_sharedACR_dt$twocatect != 'broadly_accessible__broadly_accessible',]


fml <- as.formula(paste("count ~", paste0(spe1,'_celltype'), "+", paste0(spe2,'_celltype')))
aggregate_dt <- aggregate(fml ,data =ipt_sharedACR_dt, sum)


p <- ggplot(aggregate_dt, aes(x=.data[[paste0(spe2,'_celltype')]],y=count,fill= .data[[paste0(spe1,'_celltype')]])) + 
  geom_bar(stat="identity", position=position_dodge()) +
  #scale_fill_manual(values=color_value) + 
  #geom_violin(trim=FALSE) + 
  #geom_boxplot(position = position_dodge2(preserve = "single")) + 
  #geom_boxplot(width=0.2,fill=organct) +
  #geom_jitter(shape=1, position=position_jitter(0.2)) +
  labs(x=paste0('\n',spe2, "_celltype"), y = 'ACR count\n', fill = paste0(spe1, "_celltype") )+
  theme(
    plot.title = element_text(face="bold.italic",size=35,hjust = 0.5),
    axis.title = element_text(size =35),
    axis.text.x = element_text(colour = "black", size=20,angle = 45, hjust = 1,face = 'bold'),  ##change the text to italic
    axis.text.y = element_text(size=30,colour = "black"),
    
    axis.ticks = element_line(size = rel(2.5)),
    axis.ticks.length = unit(0.5, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),   ##remove the grey background
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    strip.text = element_text(size=20),
    #legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm") ##generate margion
  ) +
  #scale_x_discrete(drop=TRUE) +
  facet_wrap(~.data[[paste0(spe1,'_celltype')]],drop=TRUE, ncol = 2)

pdf(paste0(output_dir,'/opt_synACR_celltype_count_barplot_',spe1,'_',spe2,'.pdf'),width = 15,height = 15)
p
dev.off()



############
##shift plot
#ipt_dt <- ipt_dt_seperate
#ipt_sharedACR_dt <- ipt_dt[ipt_dt[[paste0(spe1,'_ACRcate1')]] == 'shared_ACR',]
#ipt_sharedACR_dt$twocatect <- paste0(ipt_sharedACR_dt[[paste0(spe1,'_celltype')]],'__',ipt_sharedACR_dt[[paste0(spe2,'_celltype')]])
#ipt_sharedACR_dt <- ipt_sharedACR_dt[ipt_sharedACR_dt$twocatect != 'broadly_accessible__broadly_accessible',]

##split the line
#df_clean <- ipt_sharedACR_dt %>%
#  separate_rows(celltype, sep = ",")


#p <- ggplot(ipt_sharedACR_dt,
#            aes( axis1 = .data[[paste0(spe1,'_celltype')]],
#                 axis2 = .data[[paste0(spe2,'_celltype')]], y = 1)) +
#  geom_alluvium(aes(fill = .data[[paste0(spe1,'_celltype')]]), width = 1/12) +
#  geom_stratum(width = 1/12, fill = "grey", color = "black") +
#  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
#  scale_x_discrete(limits = c(paste0(spe1,'_celltype'),paste0(spe2,'_celltype')), expand = c(.05, .05)) +
#  theme_minimal() +
#  labs(title = "Correspondence Between Two Species")

#pdf(paste0(output_dir,'/opt_synACR_shiftplot_no_bACRbACR_',spe1,'_',spe2,'.pdf'),width = 10,height = 20)
#p
#dev.off()

#ipt_dt <- ipt_dt_seperate
#ipt_sharedACR_dt <- ipt_dt[ipt_dt[[paste0(spe1,'_ACRcate1')]] == 'shared_ACR',]
#p <- ggplot(ipt_sharedACR_dt,
#            aes( axis1 = .data[[paste0(spe1,'_celltype')]],
#                 axis2 = .data[[paste0(spe2,'_celltype')]], y = 1)) +
#  geom_alluvium(aes(fill = .data[[paste0(spe1,'_celltype')]]), width = 1/12) +
#  geom_stratum(width = 1/12, fill = "grey", color = "black") +
#  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
#  scale_x_discrete(limits = c(paste0(spe1,'_celltype'),paste0(spe2,'_celltype')), expand = c(.05, .05)) +
#  theme_minimal() +
#  labs(title = "Correspondence Between Two Species")

#pdf(paste0(output_dir,'/opt_synACR_shiftplot_',spe1,'_',spe2,'.pdf'),width = 10,height = 20)
#p
#dev.off()






