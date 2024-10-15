rm(list = ls()) # remove all object in the current enviroment
gc()

# You can restart the R session to avoid any R package conflicts

#load required library
library("dplyr")
library("saveplot")
library("finalfit")
library("ggplot2")
library("survminer")

#please set the working directory
#setwd("/your_path/CRC_Data_Fig7")

load("GSE17536/data/GSE17536_new.rda") # load normalized gene expression data => GSE17536_new
load("GSE17536/data/genes_GSE17536.rda") # load genes 
load("GSE17536/data/GSE17536_celltype_markers.rda")

# read CMS clusters
GSE17536_cms <- read.csv("GSE17536/data/GSE17536_cms.csv")
# remove na values
GSE17536_cms <- na.omit(GSE17536_cms)

#Select only related genes from gene expression data
GSE17536_new <- as.data.frame(t(GSE17536_new))
GSE17536_new <- GSE17536_new[,colnames(GSE17536_new) %in% genes_GSE17536]
GSE17536_new <- GSE17536_new %>% mutate(SampleID=row.names(GSE17536_new))

#load clinical data 
load("GSE17536/data/GSE17536_clinical_tumor.rda")

# Select related survival type (OS or DFS or DSS)
time_ = "DFS.Time" 
ind_ = "DFS"

# create a path for outputs
path <- paste0("GSE17536/Survival/",ind_,"/")
if(!file.exists(path)) dir.create(path, recursive=TRUE, showWarnings = FALSE)

# Get related survival data
clinical_data_tmp <- GSE17536_clinical_tumor %>% select(SampleID, all_of(ind_), all_of(time_))
clinical_data_tmp <- na.omit(clinical_data_tmp)

if(ind_=="DFS"){  # For DFS select only Stage I-II-III
  GSE17536_clinical_tumor_cafs <- GSE17536_clinical_tumor %>% filter((Stage %in% c("Stage I","Stage II","Stage III")) & (SampleID %in% clinical_data_tmp$SampleID)) %>% select(SampleID, all_of(ind_), all_of(time_)) 
}else{
  GSE17536_clinical_tumor_cafs <- GSE17536_clinical_tumor %>% filter((SampleID %in% clinical_data_tmp$SampleID)) %>% select(SampleID, all_of(ind_), all_of(time_)) 
}

GSE17536_clinical_tumor_cafs <- GSE17536_clinical_tumor_cafs %>% arrange(SampleID)
exp_data_cafs <- GSE17536_new %>% filter(SampleID %in% GSE17536_clinical_tumor_cafs$SampleID) %>% arrange(SampleID)

#check the gene expression and clinical data ordered correctly or not
print(identical(as.character(GSE17536_clinical_tumor_cafs$SampleID),as.character(exp_data_cafs$SampleID)))
#combine data
rna_cafs <- cbind(exp_data_cafs[,-ncol(exp_data_cafs)],GSE17536_clinical_tumor_cafs[,-1])

################
################
################

all_hr <- c() # to save HR results

for (i in seq_along(GSE17536_celltype_markers)) {
  if(length(GSE17536_celltype_markers[[i]])<=0) { next}
  print(names(GSE17536_celltype_markers)[i])
  #print( length(GSE17536_celltype_markers[[i]]))
  rna_genes <- as.data.frame(rna_cafs[,colnames(rna_cafs) %in% GSE17536_celltype_markers[[i]]])
  colnames(rna_genes) <- GSE17536_celltype_markers[[i]]
  
  if(length(GSE17536_celltype_markers[[i]])==1) {
    rna_avg <- rna_genes
    colnames(rna_avg) <- 'Avg_Exp'
  } else {
    rna_avg <- as.data.frame(rowMeans(rna_genes))
    names(rna_avg)[names(rna_avg) == 'rowMeans(rna_genes)'] <- 'Avg_Exp'
  }
  
  avg_cafs <- sort(rna_avg$Avg_Exp)[length(as.vector(rna_avg$Avg_Exp))/2]
  rna_avg <- rna_avg %>% mutate(Avg_Exp_lh=factor(ifelse(Avg_Exp>avg_cafs,"High","Low")))
  table(rna_avg$Avg_Exp_lh)
  n_col <- ncol(rna_cafs)
  rna_avg <- cbind(rna_avg,rna_cafs[,c(n_col-1,n_col)])
  rna_avg <- rna_avg %>% mutate(SampleID=row.names(rna_avg))
  
  
  ##############
  #CMS bar plot#
  ##############
  table(GSE17536_cms$prediction)
  
  CMS1 <- GSE17536_cms %>% filter(prediction=="CMS1")
  CMS2 <- GSE17536_cms %>% filter(prediction=="CMS2")
  CMS3 <- GSE17536_cms %>% filter(prediction=="CMS3")
  CMS4 <- GSE17536_cms %>% filter(prediction=="CMS4")

  CMS1_low_cafs <-  rna_avg %>% filter(SampleID %in% CMS1$X & Avg_Exp_lh == "Low")
  CMS1_high_cafs <-  rna_avg %>% filter(SampleID %in% CMS1$X & Avg_Exp_lh == "High")
  
  CMS2_low_cafs <-  rna_avg %>% filter(SampleID %in% CMS2$X & Avg_Exp_lh == "Low")
  CMS2_high_cafs <-  rna_avg %>% filter(SampleID %in% CMS2$X & Avg_Exp_lh == "High")
  
  CMS3_low_cafs <-  rna_avg %>% filter(SampleID %in% CMS3$X & Avg_Exp_lh == "Low")
  CMS3_high_cafs <-  rna_avg %>% filter(SampleID %in% CMS3$X & Avg_Exp_lh == "High")
  
  CMS4_low_cafs <-  rna_avg %>% filter(SampleID %in% CMS4$X & Avg_Exp_lh == "Low")
  CMS4_high_cafs <-  rna_avg %>% filter(SampleID %in% CMS4$X & Avg_Exp_lh == "High")
  
  print(paste0("CMS1-Low:",nrow(CMS1_low_cafs)," CMS1-High:",nrow(CMS1_high_cafs)))
  print(paste0("CMS2-Low:",nrow(CMS2_low_cafs)," CMS2-High:",nrow(CMS2_high_cafs)))
  print(paste0("CMS3-Low:",nrow(CMS3_low_cafs)," CMS3-High:",nrow(CMS3_high_cafs)))
  print(paste0("CMS4-Low:",nrow(CMS4_low_cafs)," CMS4-High:",nrow(CMS4_high_cafs)))
  
  rws_condition <- as.factor(rep(c("Low" , "High"), 4))
  Group <- as.factor(c("CMS1","CMS1","CMS2","CMS2","CMS3","CMS3","CMS4","CMS4"))
  rws_count <- c(nrow(CMS1_low_cafs),nrow(CMS1_high_cafs),nrow(CMS2_low_cafs),nrow(CMS2_high_cafs),nrow(CMS3_low_cafs),nrow(CMS3_high_cafs),nrow(CMS4_low_cafs),nrow(CMS4_high_cafs))
  mypal <- c("#fafd7cff","#b7e4f9ff","#fb6467ff","#5cb85cff")
  rws_all <- data.frame(rws_condition,Group,rws_count)
  
  total <- rws_all %>%
    group_by(rws_condition) %>%
    dplyr::summarise(Sum = sum(rws_count))
  
  bar_p <- ggplot(rws_all, aes(fill=Group, y=rws_count, x=rws_condition)) + 
    geom_bar(position="fill", stat="identity")+ 
    scale_fill_manual(values =mypal)+
    geom_text(aes(rws_condition, 1.05, label = Sum, fill = NULL), data = total)+
    #ggtitle(paste0("CMS Classes"))+
    labs(x=names(GSE17536_celltype_markers)[i], y="Percentage")+
    theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank())+
    theme_classic()
  
  barp_name <- paste0(path,"barplot_GSE17536_",names(GSE17536_celltype_markers)[i],".png")
  barp_name2 <- paste0(path,"barplot_GSE17536_",names(GSE17536_celltype_markers)[i],".svg")
  
  saveplot::save_png(bar_p,barp_name, width=6, height=8, ppi = 600)
  saveplot::save_svg(bar_p,barp_name2, width=6, height=8)
  
  ###################
  #Survival analysis#
  ###################
  
  rna_avg <- rna_avg %>% select(Avg_Exp_lh,all_of(ind_), all_of(time_), Avg_Exp)
  xvar <- as.matrix(rna_avg[,c(1)])
  time <- as.numeric(rna_avg[,time_])
  censor <- as.numeric(rna_avg[,ind_])
  dat <- cbind.data.frame(time, censor, xvar)
  dat <- dat %>% mutate_if(sapply(dat, is.character), as.factor)
  
  exp_names_cafs <- 'xvar'
  dependent_cafs <- "Surv(time, censor)"
  
  p_cafs <- dat %>% mutate(xvar = factor(xvar, levels = c("Low", "High"),
                                         labels = c("Low", "High"))) %>%  surv_plot(dependent_cafs, exp_names_cafs, pval = TRUE,conf.int = TRUE, surv.median.line = "hv",risk.table = TRUE, title=names(GSE17536_celltype_markers)[i])
  
  f_name <- paste0(path,"Surv_GSE17536_",names(GSE17536_celltype_markers)[i],".png")
  saveplot::save_png(p_cafs,f_name, width=6, height=6, ppi = 600)
  ###########
  
  ################
  # Hazard Ratio #
  ################
  xvar_hr <- as.matrix(rna_avg[,c(4)])
  time <- as.numeric(rna_avg[,time_])
  censor <- as.numeric(rna_avg[,ind_])
  dat_hr <- cbind.data.frame(time, censor, xvar_hr)
  
  colnames(dat_hr) <- c("time","censor","xvar")
  
  tab <- dat_hr %>% 
    finalfit(dependent_cafs, c("xvar"), add_dependent_label = FALSE, na_include = FALSE, digits = c(2,2,5)) %>%
    rename("Variable" = label) %>% 
    rename(" " = levels) %>% 
    rename("  " = all)
  
  tab <- tab %>% mutate(Cell_Type = names(GSE17536_celltype_markers)[i])  
  
  all_hr <- rbind(all_hr,tab)
  
  #HR results for individual celltype
  #write.table(tab,file = paste0(path,"HR_GSE17536_",names(GSE17536_celltype_markers)[i],".txt"), row.names = FALSE, sep = "\t")
  #############
}

# Additional file 12: Table S11
write.table(all_hr,file = paste0(path,"HR_GSE17536_ALL.txt"), row.names = FALSE, sep = "\t")
