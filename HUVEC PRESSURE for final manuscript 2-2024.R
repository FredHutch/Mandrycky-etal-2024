################################################################
#Project: scRNAseq analysis for Mandrycky et al, Under Pressure: altered endothelial flow response
#Date: 2-2024
#################################################################

### Load packages
library(monocle3)
library(leidenbase)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)
library(pheatmap)
library(Rcpp)
library(reticulate)
library(VGAM)
library(viridis)
library(magrittr)
library(mygene)
library(ggpubr)
library(scales)
library(garnett)
library(org.Mm.eg.db)
library(ggsci)
library(Rcpp)

rm(list=ls())
writeLines(capture.output(sessionInfo()), "sessionInfo-EC_PRESSURE.txt")

# Set project directories
DIR <- file.path("~/Desktop/Pressure_scrnaseq")
DIR2 <- file.path("~/Desktop/Pressure_scrnaseq/output")

simple_theme <-  theme(text = element_blank(),
                       panel.grid = element_blank(),
                       axis.title = element_blank(),
                       axis.text = element_blank(),
                       axis.ticks = element_blank(),
                       axis.line.x = element_blank(),
                       axis.line.y = element_blank(), legend.position = "none")  ### theme to remove legends and axis/font

simple_themeL <-  theme(text = element_text(size=25),
                        panel.grid = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank(),
                        axis.line.x = element_blank(),
                        axis.line.y = element_blank())   #### theme show the legend with large font

###load demultiplexed scRNAdata sets (demultiplexing via cellranger)
ECstatic_cds<-load_cellranger_data("~/Desktop/Pressure_scrnaseq/Static/")
colData(ECstatic_cds)$Sample="Static"
colData(ECstatic_cds)$flow="no"

EC0mmHg_cds<-load_cellranger_data("~/Desktop/Pressure_scrnaseq/0mmHg/")
colData(EC0mmHg_cds)$Sample="0mmHg"
colData(EC0mmHg_cds)$flow="yes"

EC20mmHg_cds<-load_cellranger_data("~/Desktop/Pressure_scrnaseq/20mmHg/")
colData(EC20mmHg_cds)$Sample="20mmHg"
colData(EC20mmHg_cds)$flow="yes"

EC60mmHg_cds<-load_cellranger_data("~/Desktop/Pressure_scrnaseq/60mmHg/")
colData(EC60mmHg_cds)$Sample="60mmHg"
colData(EC60mmHg_cds)$flow="yes"

### combine to single cds
cds<-combine_cds(list(ECstatic_cds,EC0mmHg_cds,EC20mmHg_cds,EC60mmHg_cds))

### save/read combined cds
#saveRDS(cds, file.path(DIR2, "HUVEC_PRESSURE_COMBINED.RDS"))
cds <- readRDS(file.path(DIR2, "HUVEC_PRESSURE_COMBINED.RDS"))


#filter low quality cells
cds <- detect_genes(cds, min_expr=0.1)     
cds <- cds[,colData(cds)$num_genes_expressed > 1000]
summary(pData(cds)$num_genes_expressed)

colData(cds)$n.umis <- Matrix::colSums(counts(cds))
qplot(colData(cds)$n.umis, geom="density")
qplot(log10(colData(cds)$n.umis), geom="density")
cds <- cds[,Matrix::colSums(counts(cds)) > 5000]
summary(colData(cds)$n.umis)

mito_genes <- rowData(cds)
mito_genes <- mito_genes[str_detect(mito_genes$gene_short_name, "MT-"),]
mito_genes_id <- mito_genes$id
colData(cds)$n.umis <- Matrix::colSums(counts(cds))
num_mito <- Matrix::colSums(counts(cds[mito_genes_id]))
cds$n.mito <- num_mito
perc_mito <- 100*(cds$n.mito / cds$n.umis)
cds$perc_mito_umis <- perc_mito
cds_mito_filter <- cds[,pData(cds)$perc_mito_umis < 10]
cds<- cds_mito_filter

##pre-process
cds <- preprocess_cds(cds, num_dim=10)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds)

#cluster cells and plot by cluster, sample, and static vs flow (Fig 3A)
set.seed(1)
cds<-cluster_cells(cds, random.seed=1)
plot_cells(cds, label_cell_groups = F, cell_size = 0.5, group_label_size = 5, show_trajectory_graph = FALSE)
cell_type_color1 <- c("Static"="darkgrey","0mmHg"= "brown2","20mmHg"="darkgoldenrod1","60mmHg"="cornflowerblue")
plot_cells(cds, color_cells_by = "Sample", label_cell_groups = F, cell_size = 0.5) + scale_color_manual(values = cell_type_color1)
cell_type_color2 <- c("no"="darkgrey","yes"= "red2")
plot_cells(cds, color_cells_by = "flow", label_cell_groups = F, cell_size = 0.5) + scale_color_manual(values = cell_type_color2)

#plot by expression of MKI67 and pan-EC genes (Fig 3A)
plot_cells(cds, genes=c("MKI67", "PECAM1","CDH5","KDR"),label_cell_groups = F, cell_size = 0.5) +simple_theme

#save/load processed cell data set
#saveRDS(cds, file.path(DIR2, "HUVEC_PRESSURE_PROCESSED.RDS")) 
cds <- readRDS(file.path(DIR2, "HUVEC_PRESSURE_PROCESSED.RDS"))


### detect differentially expressed genes between static vs flow (Table 1)
cds_sub <- detect_genes(cds, min_expr=10)
exp_genes <- row.names(subset(fData(cds_sub), fData(cds_sub)$num_cells_expressed>0))
cds_sub <- cds_sub[rowData(cds_sub)$id %in% exp_genes,]
gene_fits <- fit_models(cds_sub, model_formula_str = "~flow")
fit_coefs <- coefficient_table(gene_fits)
identity_terms <- fit_coefs %>% filter(term != "(Intercept)")
identity_DEG <- identity_terms %>% filter (q_value < 0.001) %>%
  dplyr::select(gene_short_name, term, q_value, estimate)
write.csv(identity_DEG, file.path(DIR2, "HUVEC_FLOW_DEG.csv"))


#### Estimate gene set scores for signature gene lists
estimate_score <- function(cds, markers){
  cds_score = cds[fData(cds)$gene_short_name %in% markers,]
  aggregate_score = exprs(cds_score)
  aggregate_score = Matrix::t(Matrix::t(aggregate_score) / pData(cds_score)$Size_Factor)
  aggregate_score = Matrix::colSums(aggregate_score)
  pData(cds)$score = log(aggregate_score +1)
  return(cds)
}

### load gene sets (Fig 3B)
# canonical flow response genes 
flow_genes <- c("KLF4", "KLF2", "SMAD6", "SMAD7", "NOS3", "PTGS2")
# canonical Notch-dependent genes
Notch_genes <- c("HEY1", "HEY2", "JAG1", "JAG2", "HES1", "SOX17", "DLL4")
# GSEA Msig database: GO Response to fluid shear stress
shear_stress_genes <- as.character(c(read.csv(file.path(DIR, "GOBP_RESPONSE_TO_FLUID_SHEAR_STRESS.csv"), header=F))$V1)
# genes identified by differential gene expression: upregulated by flow, top 50 by q-value
Flow_up_genes <- as.character(c(read.csv(file.path(DIR, "HUVEC_FLOW_DEG_UP_NEW.csv"), header=F))$V1)

#estimate scores, plot each separately (Fig 3B)
cds <- estimate_score(cds, markers = flow_genes)
cds <- estimate_score(cds, markers = Notch_genes)  
cds <- estimate_score(cds, markers = shear_stress_genes)
cds <- estimate_score(cds, markers = Flow_up_genes)  
 
### plot scaled gene set scores in UMAP (Fig 3B)
mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
plot_cells(cds, color_cells_by = 'score', cell_size = 0.8, show_trajectory_graph = F)  +
  scale_color_gradientn(colours = mycol)

### violn plot of gene set scores static vs flow (Fig 3B)
dfx <- data.frame(colData(cds)$flow, colData(cds)$score)
names(dfx) <- c("flow", "score")
p <- ggplot(dfx, aes(x= flow, y=score, fill = flow)) + geom_violin(trim=FALSE) + geom_boxplot(fill= "white", width=0.1, outlier.colour=NA)
p <- p + theme_bw() + scale_fill_manual(values = cell_type_color2)
p <- p + xlab(NULL) + ylab(NULL) ##renames axis label
p <- p + theme(plot.margin = unit(c(0,3,0,3), "cm"))
p <- p + theme(panel.grid.major= element_blank(), panel.grid.minor = element_blank())
p <- p+theme(legend.position = "none") + theme(text = element_text(size=18))
p
#Wilcoxon test
p <- p + stat_compare_means() + stat_mean()  ###statistical comparison to calculate p value between 2 samples (Wilcoxon)
p


####PRESSURE ANALYSIS (Fig 3C)
### isolate cluster 1 and remove static sample = cds2
cds2<-cds[,(colData(cds)$flow == "yes"),]
colData(cds2)$cluster=monocle3::clusters(cds2)
cds2<-cds2[,(colData(cds2)$cluster == 1),]
cell_type_color3 <- c("0mmHg"= "brown2","20mmHg"="darkgoldenrod1","60mmHg"="cornflowerblue")
plot_cells(cds2, color_cells_by = "Sample", label_cell_groups = F, cell_size = 0.8) + scale_color_manual(values = cell_type_color3)

### isolate high vs low pressure only for differential gene expression = cds3
cds3<-cds2[,(colData(cds2)$Sample != "20mmHg"),]
cell_type_color4 <- c("0mmHg"= "brown2","60mmHg"="cornflowerblue")
plot_cells(cds3, color_cells_by = "Sample", label_cell_groups = F, cell_size = 0.8) + scale_color_manual(values = cell_type_color4)

### detect DEG between low vs high Pressure (Table 1)
cds_sub <- detect_genes(cds3, min_expr=10)
exp_genes <- row.names(subset(fData(cds_sub), fData(cds_sub)$num_cells_expressed>0))
cds_sub <- cds_sub[rowData(cds_sub)$id %in% exp_genes,]
gene_fits <- fit_models(cds_sub, model_formula_str = "~Sample")
fit_coefs <- coefficient_table(gene_fits)
identity_terms <- fit_coefs %>% filter(term != "(Intercept)")
identity_DEG <- identity_terms %>% filter (q_value < 0.05) %>%
  dplyr::select(gene_short_name, term, q_value, estimate)
write.csv(identity_DEG, file.path(DIR2, "HUVEC_PRESSURE_DEG.csv"))

# load gene sets (Fig 3C)
# Top 50 differentially expressed pressure-response genes by q value (Fig 3C)
Pressure_up_genes <- as.character(c(read.csv(file.path(DIR, "HUVEC_PRES_DEG_UP_NEW.csv"), header=F))$V1)
Pressure_down_genes <- as.character(c(read.csv(file.path(DIR, "HUVEC_PRES_DEG_DOWN_NEW.csv"), header=F))$V1)
# GSEA MSIG Hallmark Myc genesV2
Myc_genes <- as.character(c(read.csv(file.path(DIR, "GSEA_Hallmark_Myc_V2.csv"), header=F))$V1)

### Gene set scores
cds2 <- estimate_score(cds2, markers = Pressure_up_genes)  
cds2 <- estimate_score(cds2, markers = Pressure_down_genes) 
cds2 <- estimate_score(cds2, markers = Myc_genes)   

### plot scaled gene set scores in UMAP (Fig 3B)
mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
plot_cells(cds2, color_cells_by = 'score', cell_size = 0.8, show_trajectory_graph = F)  +
  scale_color_gradientn(colours = mycol) + simple_theme

### violin plot of gene set scores by pressure (Fig 3B)
dfx <- data.frame(colData(cds2)$Sample, colData(cds2)$score)
names(dfx) <- c("Sample", "score")
p <- ggplot(dfx, aes(x= Sample, y=score, fill = Sample)) + geom_violin(trim=FALSE) + geom_boxplot(fill= "white", width=0.1, outlier.colour=NA)
p <- p + theme_bw() + scale_fill_manual(values = cell_type_color1)
p <- p + xlab(NULL) + ylab(NULL) ##renames axis label
p <- p + theme(plot.margin = unit(c(0,3,0,3), "cm"))
p <- p + theme(panel.grid.major= element_blank(), panel.grid.minor = element_blank())
p <- p+theme(legend.position = "none") + theme(text = element_text(size=18))
p



###MURINE EMBRYONIC EC ANALYSIS

### Load cell data set containing E8.5-E11.5 murine embryos
### Cell types classified by Garnette at Venous EC (UEC), Arterial EC/HE (AEC), and Hematopoietic (Hem)
### Combined cds from Hadland et al 2023 Nat Comm, Dignum et al 2022 Cell Reports
#saveRDS(cds, file.path(DIR, "agmE8_11_ECHem.RDS"))  ### save E8-11 EC and Hem only
cdsa <- readRDS(file.path(DIR, "agmE8_11_ECHem.RDS"))

### remove hematopoietic cell types, limit analysis to EC
cdsa<-cdsa[,(colData(cdsa)$cluster_ext_type != "Hem"),]

### plot cells by embryonic stage (Fig 3D)
plot_cells(cdsa, color_cells_by = "Stage", label_cell_groups = F, cell_size = 0.6, show_trajectory_graph = FALSE)
plot_cells(cdsa, color_cells_by = "Sample", label_cell_groups = F, cell_size = 0.6, show_trajectory_graph = FALSE)

### plot cells by cell type (Fig 3D)
cell_type_color5 <- c("UEC" = "wheat", "AEC"="orange")
plot_cells(cdsa,
           color_cells_by = "cluster_ext_type",
           show_trajectory_graph = FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1, cell_size = 0.6) +
  scale_color_manual(values = cell_type_color5)


#### Estimate gene set scores for signature gene sets (Fig 3D)
estimate_score <- function(cds, markers){
  cds_score = cds[fData(cds)$gene_short_name %in% markers,]
  aggregate_score = exprs(cds_score)
  aggregate_score = Matrix::t(Matrix::t(aggregate_score) / pData(cds_score)$Size_Factor)
  aggregate_score = Matrix::colSums(aggregate_score)
  pData(cds)$score = log(aggregate_score +1)
  return(cds)
}

##Code to swap from humnan gene symbols to mouse gene symbols
mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

convert_human_to_mouse <- function(gene_list){
  output = c()
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="human"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      mouse_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="mouse, laboratory"))[,"Symbol"]
      for(mouse_gene in mouse_genes){
        output = append(output,mouse_gene)
      }
    }
  }
  return (output)
}

### gene sets:
# Top 50 differentially expressed pressure-response genes by q value, converted to mouse genes (Fig 3D)
Pressure_up_genes_m <- convert_human_to_mouse(Pressure_up_genes)
Pressure_down_genes_m <- convert_human_to_mouse(Pressure_down_genes)
# AEC genes downregulated in Dll4 mutants (Luo, Benedito Nature 2020)
AEC2_genes <- as.character(c(read.csv(file.path(DIR, "LUO_DLL4_AEC_GENESm.csv"), header=F))$V1)
# GSEA MSig Hallmark Myc genesV2 (mouse)
Myc_genes <- as.character(c(read.csv(file.path(DIR, "GSEA_Hallmark_Myc_V2m.csv"), header=F))$V1)

### repeat for each GSS
cdsa <- estimate_score(cdsa, markers = Pressure_up_genes_m)
cdsa <- estimate_score(cdsa, markers = Pressure_down_genes_m)
cdsa <- estimate_score(cdsa, markers = AEC2_genes)
cdsa <- estimate_score(cdsa, markers = Myc_genes)  

### plot scaled gene set scores in UMAP (Fig 3D)
mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
plot_cells(cdsa, color_cells_by = 'score', cell_size = 0.8, show_trajectory_graph = F)  +
  scale_color_gradientn(colours = mycol) + simple_theme

###makes a dataframe of gene set scores by stage (in arterial EC/HE cell type only)
dfx <- data.frame(colData(cdsa)$cluster_ext_type, colData(cdsa)$Stage, colData(cdsa)$score)
names(dfx) <- c("celltype", "stage", "score")
dfx<-dfx %>% filter(celltype !="UEC")
dfx$stage<-factor(dfx$stage, levels = c("E8-9","E10-11"))
# violin plot of gene set score by stage
p <- ggplot(dfx, aes(x= stage, y=score, fill = stage)) + geom_violin(trim=FALSE) + geom_boxplot(fill= "white", width=0.1, outlier.colour=NA)
p <- p + theme_bw() #+ scale_fill_manual(values = cell_type_color3)
p <- p + xlab(NULL) + ylab("Gene-set Score") ##renames axis label
p <- p + theme(plot.margin = unit(c(0,3,0,3), "cm"))
p <- p + theme(panel.grid.major= element_blank(), panel.grid.minor = element_blank())
p <- p+theme(legend.position = "none") + theme(text = element_text(size=15))
p
#Wilcoxon test
p <- p + stat_compare_means() + stat_mean()  ###statistical comparison to calculate p value between 2 samples (Wilcoxon)
p

###makes a dataframe of gene set scores by EC type (limited to E10-E11 samples)
dfx <- data.frame(colData(cdsa)$cluster_ext_type, colData(cdsa)$Stage, colData(cdsa)$score)
names(dfx) <- c("celltype", "stage", "score")
dfx<-dfx %>% filter(stage =="E10-11")
dfx$celltype<-factor(dfx$celltype, levels = c("UEC","AEC"))
cell_type_color5 <- c("UEC" = "wheat", "AEC"="orange")
p <- ggplot(dfx, aes(x= celltype, y=score, fill = celltype)) + geom_violin(trim=FALSE) + geom_boxplot(fill= "white", width=0.1, outlier.colour=NA)
p <- p + theme_bw() + scale_fill_manual(values = cell_type_color5)
p <- p + xlab(NULL) + ylab("Gene-set Score") ##renames axis label
p <- p + theme(plot.margin = unit(c(0,3,0,3), "cm"))
p <- p + theme(panel.grid.major= element_blank(), panel.grid.minor = element_blank())
p <- p+theme(legend.position = "none") + theme(text = element_text(size=15))
p
#Wilcoxon test
p <- p + stat_compare_means() + stat_mean()  ###statistical comparison to calculate p value between 2 samples (Wilcoxon)
p


