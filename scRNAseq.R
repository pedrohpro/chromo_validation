
rm(list = ls())

########################
### loading library  ###
########################

library(patchwork)
library(devtools)
library(tidyverse) #do it all
library(Seurat) # single cell analysis
library(ComplexHeatmap)
library(pheatmap)
library(data.table)
library(org.Hs.eg.db) #Homo sapiens OrgDb
library(biomaRt) #gene annotations
library(RColorBrewer) # for pretty colors
library(ggrepel)
library(umap)
library(VennDiagram) #venn digrams
library(clusterProfiler) #functional analysis
library(pathview)
library(cowplot)
library(GOplot)
library(enrichplot)
library(gage)
library(vsn)
library(RColorBrewer) # colors pallets
library(scales) # colors pallets
library(scDblFinder) # doublet finder
library(SingleCellExperiment)
library(AUCell) # get set score
library(GSEABase) # build gene set function
library(sessioninfo) # package versions
library(UpSetR)
library(CellChat)
library(circlize)
library(WGCNA)
library(DESeq2)
library(CEMiTool)
library(cluster)
library(metap)
library(readr)
library(writexl)
library(chromo)

set.seed(42)

scdown <- readRDS("data/sc/processed_seurat_obj.rds")

#saveRDS(scdown, file = "data/sc/processed_seurat_obj.rds")

# List of packages with source
packages <- list(
  ComplexHeatmap     = "Bioc",
  pheatmap           = "CRAN",
  data.table         = "CRAN",
  org.Hs.eg.db       = "Bioc",
  biomaRt            = "Bioc",
  RColorBrewer       = "CRAN",
  ggrepel            = "CRAN",
  umap               = "CRAN",
  VennDiagram        = "CRAN",
  clusterProfiler    = "Bioc",
  pathview           = "Bioc",
  cowplot            = "CRAN",
  GOplot             = "CRAN",
  enrichplot         = "Bioc",
  gage               = "Bioc",
  karyoploteR        = "Bioc",
  vsn                = "Bioc",
  scales             = "CRAN",
  scDblFinder        = "Bioc",
  SingleCellExperiment = "Bioc",
  AUCell             = "Bioc",
  GSEABase           = "Bioc",
  sessioninfo        = "CRAN",
  UpSetR             = "CRAN",
  CellChat           = "CRAN", # If not on CRAN, install from GitHub
  circlize           = "CRAN",
  WGCNA              = "CRAN",
  hdWGCNA            = "CRAN", # If not on CRAN, install from GitHub
  DESeq2             = "Bioc",
  CEMiTool           = "Bioc",
  cluster            = "CRAN",
  metap              = "CRAN"
)

# Install packages
for (pkg in names(packages)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    source <- packages[[pkg]]
    message(paste("Installing", pkg, "from", source))
    if (source == "CRAN") {
      install.packages(pkg)
    } else if (source == "Bioc") {
      BiocManager::install(pkg)
    }
  } else {
    message(paste(pkg, "already installed."))
  }
}






#############################
###   Bar plot function   ###
#############################

barfunc <- function(
    df,
    x_axis,
    y_axis,
    fill_groups = F,
    annot_loc = "center", # "center" or "repel" or "top" or F
    fill_col = F,
    title_name = NULL,
    annot_text = y_axis,
    size_annot = 4,
    size_xaxis_text = 15,
    size_xaxis_title = 15,
    size_yaxis_text = 15,
    size_yaxis_title = 15,
    rotate_x = 0,
    x_title = F,
    y_title = F,
    extra_annot = NULL
){
  # Base plot
  bar_plot <- ggplot(df, aes(x = !!sym(x_axis), y = !!sym(y_axis))) +
    geom_bar(stat = "identity") +
    labs(
      title = title_name,
      x = x_axis,
      y = y_axis
    ) +
    theme_minimal() +
    scale_y_continuous(expand = c(0, 0)) + # remove extra padding between bars and x axis
    theme(
      plot.background = element_rect(fill = "white", color = NA),  # Set background color to white
      panel.background = element_rect(fill = "white", color = NA),  # Set panel background color to white
      panel.border = element_blank(),  # Remove panel borders
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.line = element_line(color = "black", linewidth = 0.5),  # Make axes lines bold
      axis.title.x = element_text(face = "bold", size = size_xaxis_title),  # Make axis titles bold
      axis.title.y = element_text(face = "bold", size = size_yaxis_title),  # Make axis titles bold
      axis.text.y = element_text(face = "bold", size = size_yaxis_text),  # Make axis text bold
      axis.text.x = element_text(face = "bold", angle = rotate_x, hjust = ifelse(rotate_x != 0, 1, 0.5), size = size_xaxis_text), # rotate x axis test 45 degrees
      axis.ticks.x = element_blank(),  # Remove x-axis ticks
      axis.ticks.y = element_line(color = "black", linewidth = 0.5)  # Make y-axis ticks bold
      #legend.position = "none" # No legend
    )
  # y axis ticks
  if(size_yaxis_text == 0){
    bar_plot <- bar_plot + theme(axis.ticks.y = element_blank())
  }
  
  # Conditionally add fill aesthetic and scale_fill_manual
  if (fill_groups != F) {
    bar_plot <- bar_plot +
      aes(fill = !!sym(fill_groups)) +
      labs(fill = fill_groups)
  }
  
  # coloring bars
  if(is.list(fill_col)){
    bar_plot <- bar_plot +
      scale_fill_manual(values = fill_col)
  } else if (fill_col != F) {
    bar_plot <- bar_plot +
      geom_bar(stat = "identity", fill = fill_col) +
      theme(legend.position = "none")
  }
  
  # Annotations location
  if(!is.null(extra_annot)){
    bar_plot <- bar_plot +
      annotate(
        geom = "text",
        x = extra_annot[[x_axis]],
        y = extra_annot$y,
        label = extra_annot$annotation_text,
        color = extra_annot$color,
        size = 3,
        angle = 90,
        fontface = "bold"
      )
  }else if (annot_loc == "center") {
    bar_plot <- bar_plot +
      geom_text(
        aes(label = !!sym(annot_text)),
        position = position_stack(vjust = 0.5),
        size = size_annot
      )
  } else if (annot_loc == "repel") {
    bar_plot <- bar_plot +
      geom_text_repel(
        aes(label = !!sym(annot_text)),
        direction = "y",
        size = size_annot
      )
  } else if(annot_loc == "top"){
    bar_plot <- bar_plot +
      scale_y_continuous(expand = c(0,0), limits = c(0, 1.1*max(ifelse(fill_groups != F, df %>% group_by(!!sym(fill_groups)) %>% pull(!!sym(y_axis)), df[[y_axis]])))) +
      geom_text(
        aes(label = !!sym(annot_text)),
        vjust = -0.2,
        size = size_annot
      )
  }
  # x label title
  if(x_title != F){
    bar_plot <- bar_plot +
      xlab(x_title)
  }
  # y label title
  if(y_title != F){
    bar_plot <- bar_plot +
      ylab(y_title)
  }
  
  return(bar_plot)
}


# ############################
# ###         counts       ###
# ############################
# 
# # Runs names lists (maybe better to get this from meta data)
# runs <- sapply(43:71, function(i) paste0("EGAN000033605", i), simplify = T)
# 
# # Raw counts
# counts <- list()
# for (run in runs){ ################# algumas amostras apenas PROVISÓRIO PARA RODAR LOCAL PORQUE ESTÁ MUITO PESADO!!!
#   counts[[run]] <- Read10X(paste0("data/sc/counts/",run,"/"))
# }
# 
# # Separate objects
# seurats <- list()
# for (run in names(counts)) {
#   seurats[[run]] <- CreateSeuratObject(counts = counts[[run]], project = run)
# }
# 
# # merging seurat objects (without batch correction)
# scdown <- merge(
#   x = seurats$EGAN00003360543,
#   y = seurats[names(seurats) != "EGAN00003360543"],
#   add.cell.ids = runs,
#   project = "allRuns"
# )
# 
# # joining layers
# scdown <- JoinLayers(scdown)
# 
# 
# ############################
# ###      meta data       ###
# ############################
# 
# meta <- read.delim("data/sc/meta.csv", sep = ",")
# 
# meta <- meta %>% 
#   dplyr::select(accession_id, biological_sex, subject_id, phenotype) %>% 
#   rename(
#     condition = phenotype,
#     sex = biological_sex,
#     donor = subject_id
#   ) %>% 
#   mutate(
#     condition = case_when(
#       condition == "non-diseased" ~ "CT",
#       condition == "Down Syndrome" ~ "DS"
#     ),
#     age = substr(donor, start = nchar(donor) - 1, stop = nchar(donor) - 1),
#     age = factor(age, levels = c("Y", "M", "O")) 
#   )
# 
# # Including the sample meta data into the merged cell meta data
# scdown@meta.data <- scdown@meta.data %>% 
#   rownames_to_column("cell") %>% 
#   left_join(meta, by = c("orig.ident" = "accession_id")) %>% 
#   column_to_rownames("cell")
#
#
# ############################
# ###   filtering cells    ###
# ############################ 224k
# 
# # Reads per cell
# scdown <- subset(scdown, subset = nCount_RNA > 500 & nCount_RNA < 25000)
# 
# # Features per cell
# scdown <- subset(scdown, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
# 
# # MT genes
# scdown[["percent.mt"]] <- PercentageFeatureSet(scdown, pattern = "^MT-")
# scdown <- subset(scdown, subset = percent.mt < 5)
# 
# # 123k cells
# 
# # doublets detection (removing 39k cells)
# sce <- as.SingleCellExperiment(scdown)
# set.seed(42) # fixes scDblFinder stochastic component
# sce <- scDblFinder(sce)
# scdown$doublet <- sce$scDblFinder.class
# scdown <- subset(scdown, subset = doublet == "singlet")
# 
# 
# ###############################
# ###   filtering features    ###
# ###############################
# 
# chrObj <- chromoInitiate( # filters for chromosome and MT transcripts
#   data.frame(Symbol = rownames(scdown), log2fc = 0, pval = 1), 
#   gene_col = "Symbol",
#   fc_col = "log2fc",
#   p_col = "pval"
# )
# 
# scdown <- scdown[rownames(scdown) %in% chrObj@data$Symbol,]
# 
# # Removing genes with zero expression
# nonzero_genes <- rowSums(scdown[["RNA"]]@layers[["counts"]]) > 0
# scdown <- scdown[nonzero_genes, ]
# 
# 
# ######################
# ###   Processing   ###
# ######################
# 
# scdown <- NormalizeData(scdown)
# 
# scdown <- FindVariableFeatures(scdown) # default 2000 features
# 
# scdown <- ScaleData(scdown, features = rownames(scdown))
# 
# scdown <- RunPCA(scdown, features = VariableFeatures(object = scdown)) 
# 
# elbow <- ElbowPlot(scdown, ndims = 50) + # help define number of PCs 
#   theme(
#     plot.background = element_rect(fill = "white"), 
#     panel.background = element_rect(fill = "white"), 
#     panel.border = element_blank(), 
#     axis.line = element_line(color = "black"),  
#     legend.background = element_rect(fill = "white"),  
#   )
# ggsave("results/sc/sup/elbow.png", plot = elbow, height = 4, width = 4)
# 
# scdown <- RunUMAP(scdown, dims = 1:10)
# 
# scdown <- FindNeighbors(scdown, dims = 1:10, k.param = 30)
# 
# scdown <- FindClusters(scdown, resolution = 0.2)
# 
# 
#########################################
###      Batch effect verification    ###
#########################################

# UMAP
for (i in c("seurat_clusters","donor","sex","condition","age")) {
  clusters <- DimPlot(scdown, reduction = "umap", group.by = i, pt.size = 0.1, label = TRUE, raster=FALSE) +
    ggtitle("") +
    labs(color = i) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    coord_fixed() +
    theme(
      legend.text = element_text(size = 26),  # Adjust legend text size and style
      legend.title = element_text(size = 28, face = "bold"),  # Adjust legend text size and style
      axis.text = element_text(size = 16),    # Adjust axis text size and style
      axis.title = element_text(size = 22, face = "bold"),   # Adjust axis title size and style
      axis.line = element_line(linewidth = 1.5), # Customize axis lines
      axis.ticks = element_line(linewidth = 1.5),                 # Adjust axis ticks size
    )

  ggsave(paste0("results/sc/sup/umap_",i,".png"), plot = clusters, width = 12, height = 12)
}

# Bar plots
for (i in c("donor","sex","condition","age")) {
  aux <- scdown@meta.data %>%
    dplyr::select(seurat_clusters, !!sym(i)) %>%
    group_by(seurat_clusters, !!sym(i)) %>%
    summarise(num_cells = n())%>%
    group_by(seurat_clusters) %>%
    mutate(pct = num_cells / sum(num_cells) * 100) %>%
    ungroup()

  # pct value bar plot
  bar_plot <- barfunc(aux, x_axis = "seurat_clusters", y_axis = "pct", fill_groups = i, annot_loc = "center", annot_text = "num_cells", size_annot = 2)
  ggsave(paste0("results/sc/sup/pctbar_",i,".png"), plot = bar_plot, width = 8, height = 4)
}


####################################
###       Finding markers        ###
####################################

markers_df <- read.delim("data/sc/celltypes_markers.csv", sep = ";")
celltypes_markers <- markers_df$marker

# Verifying assay type. Should be RNA
DefaultAssay(scdown)

# Changing identity to cluster column
Idents(scdown) <- "seurat_clusters"

# hline annotation
auxmark <- markers_df %>% 
  group_by(celltype) %>% 
  summarise(n_markers = n()) %>% 
  mutate(celltype = factor(celltype, levels = unique(markers_df$celltype))) %>% 
  arrange(celltype) %>%  
  mutate(
    y = length(unique(scdown@meta.data$seurat_clusters))*1.05,
    yend = y,
    x = 1, # (just for the loop to work)
    xend = n_markers, # provisory (just for the loop to work)
    colour = "#000000",
    celltype = case_when(
      celltype == "Ast" ~ "Astrocyte",
      celltype == "End" ~ "Endothelial",
      celltype == "Exc" ~ "Excitatory neuron",
      celltype == "Inh" ~ "Inhibitory neuron",
      celltype == "Mic" ~ "Microglia",
      celltype == "Oli" ~ "Oligodendrocyte",
      celltype == "OPC" ~ "Oligodendrocyte Progenitor",
      celltype == "Per" ~ "Pericyte",
      celltype == "NPC" ~ "Neural Progenitor"
    )
  ) %>% 
  as.data.frame()

for (i in 2:nrow(auxmark)) {
  auxmark[i,"x"] <- auxmark[i-1,"xend"] + 1
  auxmark[i,"xend"] <- auxmark[i,"x"] + auxmark[i,"n_markers"] -1
}

# canonical markers
dot_plot <- DotPlot(
  object = scdown, 
  features = celltypes_markers, 
  cols = c("blue","red"),
  cluster.idents = T,
  scale = T 
)+
  theme(
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_blank(),  # Remove panel borders
    axis.line = element_line(color = "black"),  # Set axis line color
    legend.background = element_rect(fill = "white"),  # Set legend background color to white
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "bold"),
    axis.text.y = element_text(size = 17, face = "bold"),
    #legend.position = "none" # no legend
  ) +
  coord_cartesian(clip = "off")

for (i in 1:nrow(auxmark)) {
  dot_plot <- dot_plot + 
    annotate(
      "segment", # type of annotation to be added
      x = auxmark[i,"x"],
      y = auxmark[i,"y"],
      xend = auxmark[i,"xend"],
      yend = auxmark[i,"yend"],
      #colour = auxmark[i,"colour"],
      
    ) +
    annotate(
      "text", # type of annotation to be added
      x = (auxmark[i,"x"] + auxmark[i,"xend"])/2, 
      y = auxmark[i,"y"] * 1.03, 
      label = auxmark[i,"celltype"], 
      #hjust = 1.1, 
      #vjust = -0.5, 
      color = "black",
      fontface = "bold",
      size = 4
    )
}

ggsave("results/sc/sup/canonical_markers.png", plot = dot_plot, height = 5, width = 22)


###############################
###       Annotating        ###
###############################

scdown@meta.data$old_clusters <- scdown@meta.data$seurat_clusters

Idents(scdown) <- "seurat_clusters"

scdown <- RenameIdents(scdown, c(
  `0` = "Oli",
  `1` = "Mic",
  `2` = "Ast",
  `3` = "Exc",
  `4` = "Ast",
  `5` = "Inh",
  `6` = "OPC",
  `7` = "Exc",
  `8` = "End",
  `9` = "Per",
  `10` = "Mic",
  `11` = "Oli"
  )
)

scdown@meta.data <- scdown@meta.data %>%
  mutate(
    seurat_clusters = case_when(
      seurat_clusters == "0" ~ "Oli",
      seurat_clusters == "1" ~ "Mic",
      seurat_clusters == "2" ~ "Ast",
      seurat_clusters == "3" ~ "Exc",
      seurat_clusters == "4" ~ "Ast",
      seurat_clusters == "5" ~ "Inh",
      seurat_clusters == "6" ~ "OPC",
      seurat_clusters == "7" ~ "Exc",
      seurat_clusters == "8" ~ "End",
      seurat_clusters == "9" ~ "Per",
      seurat_clusters == "10" ~ "Mic",
      seurat_clusters == "11" ~ "Oli"
    )
  )

# ordering in decreasing number of cells
scdown@meta.data$seurat_clusters <- factor(scdown@meta.data$seurat_clusters, levels = c("Oli","OPC","Exc","Inh","Ast","Mic","End","Per"))


###################################
###       Annotated UMAP        ###
###################################

umap_plot <- DimPlot(
  scdown, 
  reduction = "umap", 
  label = T, 
  pt.size = 0.1, 
  raster = FALSE, 
  label.size = 12, 
  cols = c("Oli" = brewer.pal(7, "Set3")[1], "Exc" = brewer.pal(8, "Set3")[8], "OPC" = brewer.pal(7, "Set3")[3], 
           "Inh" = brewer.pal(7, "Set3")[4], "Ast" = brewer.pal(7, "Set3")[7], "Mic" = brewer.pal(7, "Set3")[6], 
           "End" = brewer.pal(7, "Set3")[5], "Per" = brewer.pal(8, "Set3")[2])
) +
  NoAxes() + # removes axis
  theme(
    legend.position = "none"
  )
ggsave("results/sc/ANNOTATED_named_umap.png", plot = umap_plot, height = 12, width = 12, dpi = 600)


####################################################
###      Repeating Batch effect verification     ###
####################################################

# Bar plots
for (i in c("donor","sex","condition","age")) {
  aux <- scdown@meta.data %>% 
    select(seurat_clusters, !!sym(i)) %>% 
    group_by(seurat_clusters, !!sym(i)) %>% 
    summarise(num_cells = n())%>%
    group_by(seurat_clusters) %>%
    mutate(pct = num_cells / sum(num_cells) * 100) %>% 
    ungroup()
  
  # pct value bar plot
  bar_plot <- barfunc(aux, x_axis = "seurat_clusters", y_axis = "pct", fill_groups = i, annot_loc = "center", annot_text = "num_cells", size_annot = 2)
  ggsave(paste0("results/sc/sup/ANNOTATED_pctbar_",i,".png"), plot = bar_plot, width = 8, height = 4)
}


##################################################
###          Repeating markers plots           ###
##################################################

Idents(scdown) <- "seurat_clusters"

# hline annotation
auxmark <- markers_df %>% 
  group_by(celltype) %>% 
  summarise(n_markers = n()) %>% 
  mutate(celltype = factor(celltype, levels = unique(markers_df$celltype))) %>% 
  arrange(celltype) %>%  
  mutate(
    y = length(unique(scdown@meta.data$seurat_clusters))*1.05,
    yend = y,
    x = 1, # (just for the loop to work)
    xend = n_markers, # provisory (just for the loop to work)
    colour = "#000000",
    celltype = case_when(
      celltype == "Ast" ~ "Astrocyte",
      celltype == "End" ~ "Endothelial",
      celltype == "Exc" ~ "Excitatory neuron",
      celltype == "Inh" ~ "Inhibitory neuron",
      celltype == "Mic" ~ "Microglia",
      celltype == "Oli" ~ "Oligodendrocyte",
      celltype == "OPC" ~ "Oligodendrocyte Progenitor",
      celltype == "Per" ~ "Pericyte",
      celltype == "NPC" ~ "Neural Progenitor"
    )
  ) %>% 
  as.data.frame()

for (i in 2:nrow(auxmark)) {
  auxmark[i,"x"] <- auxmark[i-1,"xend"] + 1
  auxmark[i,"xend"] <- auxmark[i,"x"] + auxmark[i,"n_markers"] -1
}

# canonical markers
dot_plot <- DotPlot(
  object = scdown, 
  features = celltypes_markers, 
  cols = c("blue","red"),
  cluster.idents = T,
  #cluster.features = T, # não existe!!!
  scale = T # zscore 
)+
  theme(
    plot.background = element_rect(fill = "white"),  # Set background color to white
    panel.background = element_rect(fill = "white"),  # Set panel background color to white
    #panel.grid.major = element_blank(),  # Remove major grid lines
    #panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove panel borders
    axis.line = element_line(color = "black"),  # Set axis line color
    legend.background = element_rect(fill = "white"),  # Set legend background color to white
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    #legend.position = "none" # no legend
  ) +
  coord_cartesian(clip = "off")

for (i in 1:nrow(auxmark)) {
  dot_plot <- dot_plot + 
    annotate(
      "segment", # type of annotation to be added
      x = auxmark[i,"x"],
      y = auxmark[i,"y"],
      xend = auxmark[i,"xend"],
      yend = auxmark[i,"yend"],
      #colour = auxmark[i,"colour"],
      
    ) +
    annotate(
      "text", # type of annotation to be added
      x = (auxmark[i,"x"] + auxmark[i,"xend"])/2, 
      y = auxmark[i,"y"] * 1.03, 
      label = auxmark[i,"celltype"], 
      #hjust = 1.1, 
      #vjust = -0.5, 
      color = "black",
      fontface = "bold",
      size = 4
    )
}

ggsave("results/sc/sup/ANNOTATED_dotplot.png", plot = dot_plot, height = 4, width = 20)


#################################
###      Bias verification    ###
#################################

# Number of cells
aux <- scdown@meta.data %>% 
  mutate(age_sex = paste0(age, "_", sex)) %>% 
  group_by(age_sex, seurat_clusters) %>% 
  summarise(num_cells = n()) %>%
  group_by(age_sex) %>%
  mutate(pct = num_cells / sum(num_cells) * 100) %>%
  ungroup()
bar_plot <- barfunc(aux, x_axis = "age_sex", y_axis = "pct", fill_groups = "seurat_clusters", annot_loc = "center", annot_text = "num_cells")
ggsave("results/sc/sup/bias.png", plot = bar_plot, width = 15, height = 15)

# number of biological replicates
aux <- scdown@meta.data %>% 
  group_by(donor, condition, age, sex) %>% 
  summarise() %>% 
  mutate(age_sex = paste0(age, "_", sex)) %>% 
  group_by(condition, age_sex) %>% 
  summarise(n_donors = n())
bar_plot <- barfunc(aux, x_axis = "condition", y_axis = "n_donors", fill_groups = "age_sex", annot_loc = "center", annot_text = "n_donors")
ggsave("results/sc/sup/bias_donors.png", plot = bar_plot, width = 5, height = 5)


###########################################
###        Cell type composition        ###
###########################################

# Absolute value
aux <- scdown@meta.data %>%
  mutate(samp_id = paste0(orig.ident, "_", sex, "_", age, "_", condition)) %>% 
  group_by(samp_id, orig.ident, condition, sex, age, seurat_clusters) %>%
  summarise(num_cells = n()) %>%
  group_by(orig.ident) %>%
  mutate(pct = num_cells / sum(num_cells) * 100) %>%
  ungroup() %>% 
  arrange(age, condition) %>% 
  mutate(samp_id = factor(samp_id, levels = unique(samp_id)))

# Percentage value
pct_composition <- barfunc(
  aux, 
  x_axis = "samp_id", 
  y_axis = "pct", 
  fill_groups = "seurat_clusters", 
  annot_loc = "center", 
  annot_text = "num_cells",
  #fill_col = c(brewer.pal(n = 10, name = "Set3"), brewer.pal(n = nlevels(aux$seurat_clusters)-10, name = "Paired")),
  size_annot = 3,
  size_yaxis_text = 12,
  rotate_x = 45, 
  size_xaxis_text = 8
) +
  coord_cartesian(clip = "off")

ggsave("results/sc/sup/ANNOTATED_pct_composition.png", plot = pct_composition, height = 5, width = 11)


# amount of cells between conditions for each cell type: original study inh neuron to total neuro. Ver se tem mais microglia em DS. Etc... vários possibilidades!!!

cells_by_sample <- scdown@meta.data %>%
  group_by(orig.ident, seurat_clusters, condition, age) %>%
  summarise(num_cells = n()) %>% 
  group_by(orig.ident) %>%
  mutate(pct = num_cells / sum(num_cells) * 100) %>%
  ungroup()

plot_list <- list()

cell_comp_pvals <- data.frame(celltype = NA, age = NA, pval = NA)

for (i in unique(scdown@meta.data$seurat_clusters)) {
  
  aux <- cells_by_sample %>%
    filter(seurat_clusters == i)
  
  cell_comp_pvals <- rbind(cell_comp_pvals, 
    data.frame(
      celltype = c(i,i), 
      age = c("Y","O"), 
      pval = c(
        wilcox.test(aux %>% filter(age == "Y", condition == "DS") %>% pull(pct), aux %>% filter(age == "Y", condition == "CT") %>% pull(pct))$p.value,
        wilcox.test(aux %>% filter(age == "O", condition == "DS") %>% pull(pct), aux %>% filter(age == "O", condition == "CT") %>% pull(pct))$p.value
      ) # nonparametric Wilcoxon rank-sum test (Mann–Whitney U test)
    )
  )
  
  plot_list[[i]] <- ggplot(aux, aes(x = age, y = pct, fill = condition)) +
    geom_boxplot(color="black", alpha=0.5, position=position_dodge(width=0.8), outlier.shape=NA) +
    geom_point(aes(color=condition), position=position_jitterdodge(dodge.width=0.8, jitter.width=0.1)) +
    scale_fill_manual(values=c("DS"="#2255aa", "CT"="#aa4444")) +
    scale_color_manual(values=c("DS"="#2299ff", "CT"="#ff6666")) +
    theme(
      panel.background = element_rect(fill="white"),
      panel.grid       = element_blank(),
      axis.line        = element_line(color="black"),
      axis.text        = element_text(color="black"),
      axis.ticks       = element_line(color="black"),
      axis.title.y     = element_blank(),
      axis.title.x     = element_blank(),
      legend.position  = "none"
    ) +
    ggtitle(i)
}

combined_plots <- plot_grid(plotlist=plot_list, ncol=2, nrow=4)
ggsave("results/sc/combined_boxplots.png", combined_plots, height = 6, width = 6)


##############################################
###    Differential expression analysis    ###
##############################################

# # All samples: DS vs CT
# scdown@meta.data$type_cond <- paste0(scdown@meta.data$seurat_clusters, "_", scdown@meta.data$condition)
# Idents(scdown) <- "type_cond"
# DElist_all_ages <- list()
# for (i in unique(scdown$seurat_clusters)){
#   DElist_all_ages[[i]] <- FindMarkers(
#     scdown,
#     ident.1 = paste0(i,"_DS"), # find markers for this identity
#     ident.2 = paste0(i,"_CT"), # against this identity
#     logfc.threshold = 0.0, # to analyse every gene
#     min.pct = 0.0, # to analyse every gene
#     min.cells.feature = 0
#   )
# }
# DElist_all_ages <- imap(DElist_all_ages, function(i, ct){
#   i <- i %>%
#     mutate(
#       celltype = ct,
#       gene = rownames(i)
#     )
# })

# Young: DS vs CT
scdown@meta.data$type_cond <- paste0(scdown@meta.data$seurat_clusters, "_", scdown@meta.data$age, "_", scdown@meta.data$condition)
Idents(scdown) <- "type_cond"
DElist_young <- list()
for (i in unique(scdown$seurat_clusters)){
  DElist_young[[i]] <- FindMarkers(
    scdown,
    ident.1 = paste0(i,"_Y_DS"), # find markers for this identity
    ident.2 = paste0(i,"_Y_CT"), # against this identity
    logfc.threshold = 0.0, # to analyse every gene
    min.pct = 0.0, # consider using 5% expression minimum. NOTE: different cell types WILL have different amount of rows!
    min.cells.feature = 0
  )
}
DElist_young <- imap(DElist_young, function(i, ct){
  i <- i %>% 
    mutate(
      celltype = as.character(ct)
    ) %>% 
    rownames_to_column("gene")
})

# Old: DS vs CT
DElist_old <- list()
for (i in unique(scdown$seurat_clusters)){
  DElist_old[[i]] <- FindMarkers(
    scdown,
    ident.1 = paste0(i,"_O_DS"), # find markers for this identity
    ident.2 = paste0(i,"_O_CT"), # against this identity
    logfc.threshold = 0.0, # to analyse every gene
    min.pct = 0.0, # to analyse every gene
    min.cells.feature = 0
  )
}
DElist_old <- imap(DElist_old, function(i, ct){
  i <- i %>% 
    mutate(
      celltype = as.character(ct)
    ) %>% 
    rownames_to_column("gene")
})

# Saving suplementary tables
combined_young <- do.call(rbind, DElist_young)
#write.table(combined_young, file = "data/sc/ST/DE_young.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
combined_old <- do.call(rbind, DElist_old)
#write.table(combined_old, file = "data/sc/ST/DE_old.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write_xlsx(list(DE_young = combined_young, DE_old = combined_old), path = "data/sc/ST/ST1.xlsx")


# # Young men: DS vs CT
# scdown@meta.data$type_cond <- paste0(scdown@meta.data$seurat_clusters, "_", scdown@meta.data$age, "_",  scdown@meta.data$sex, "_", scdown@meta.data$condition)
# Idents(scdown) <- "type_cond"
# DElist_young_men <- list()
# for (i in unique(scdown$seurat_clusters)){
#   DElist_young_men[[i]] <- FindMarkers(
#     scdown,
#     ident.1 = paste0(i,"_Y_male_DS"), # find markers for this identity
#     ident.2 = paste0(i,"_Y_male_CT"), # against this identity
#     logfc.threshold = 0.0, # to analyse every gene
#     min.pct = 0.0, # to analyse every gene
#     min.cells.feature = 0
#   )
# }
# DElist_young_men <- imap(DElist_young_men, function(i, ct){
#   i <- i %>% 
#     mutate(
#       celltype = ct,
#       gene = rownames(i)
#     )
# })
# 
# # Old men: DS vs CT
# DElist_old_men <- list()
# for (i in unique(scdown$seurat_clusters)){
#   DElist_old_men[[i]] <- FindMarkers(
#     scdown,
#     ident.1 = paste0(i,"_O_male_DS"), # find markers for this identity
#     ident.2 = paste0(i,"_O_male_CT"), # against this identity
#     logfc.threshold = 0.0, # to analyse every gene
#     min.pct = 0.0, # to analyse every gene
#     min.cells.feature = 0
#   )
# }
# DElist_old_men <- imap(DElist_old_men, function(i, ct){
#   i <- i %>% 
#     mutate(
#       celltype = ct,
#       gene = rownames(i)
#     )
# })
# 
# # Old women: DS vs CT
# DElist_old_women <- list()
# for (i in unique(scdown$seurat_clusters)){
#   DElist_old_women[[i]] <- FindMarkers(
#     scdown,
#     ident.1 = paste0(i,"_O_female_DS"), # find markers for this identity
#     ident.2 = paste0(i,"_O_female_CT"), # against this identity
#     logfc.threshold = 0.0, # to analyse every gene
#     min.pct = 0.0, # to analyse every gene
#     min.cells.feature = 0
#   )
# }
# DElist_old_women <- imap(DElist_old_women, function(i, ct){
#   i <- i %>% 
#     mutate(
#       celltype = ct,
#       gene = rownames(i)
#     )
# })


#####################################################
###       Chromo composition all cell types       ###
#####################################################

# # All ages: DS vs CT
# all_DEdf <- do.call(rbind, DElist_all_ages)
# all_DEdf$celltype <- factor(all_DEdf$celltype, levels = c("Mic","Ast","End","Oli","Inh","Exc","OPC","Per"))
# chrObj <- chromoInitiate(DEdf = all_DEdf, gene_col = "gene", fc_col = "avg_log2FC", p_col = "p_val_adj", celltype_col = "celltype")
# chrObj <- chromoComposition(chrObj, score_method = "hyp_padj", separate_by = "celltype", only_expr_features = T)
# allPlot <- chromoCompositionPlot(chrObj, size_dot_alt = 0.8, size_dot_no = 0.4) +
#   ggtitle("All ages")

# Young: DS vs CT
all_DEdf <- do.call(rbind, DElist_young)
all_DEdf <- all_DEdf %>% filter(gene != "XIST")
all_DEdf$celltype <- factor(all_DEdf$celltype, levels = c("Mic","Ast","End","Oli","Inh","Exc","OPC","Per"))
chrObj <- chromoInitiate(DEdf = all_DEdf, gene_col = "gene", fc_col = "avg_log2FC", p_col = "p_val_adj", celltype_col = "celltype")
chrObj <- chromoComposition(chrObj, score_method = "hyp_padj", separate_by = "celltype", only_expr_features = T)
youngPLot <- chromoCompositionPlot(chrObj, size_dot_alt = 0.8, size_dot_no = 0.4, title_xaxis = NULL, title_yaxis = NULL) +
  ggtitle("Young")

# Old: DS vs CT
all_DEdf <- do.call(rbind, DElist_old)
all_DEdf$celltype <- factor(all_DEdf$celltype, levels = c("Oli","Mic","Ast","Inh","Exc","End","OPC","Per"))
chrObj <- chromoInitiate(DEdf = all_DEdf, gene_col = "gene", fc_col = "avg_log2FC", p_col = "p_val_adj", celltype_col = "celltype")
chrObj <- chromoComposition(chrObj, score_method = "hyp_padj", separate_by = "celltype", only_expr_features = T)
oldPlot <- chromoCompositionPlot(chrObj, size_dot_alt = 0.8, size_dot_no = 0.4, title_xaxis = NULL, title_yaxis = NULL) +
  ggtitle("Old")

# # Young male: DS vs CT
# all_DEdf <- do.call(rbind, DElist_young_men)
# all_DEdf$celltype <- factor(all_DEdf$celltype, levels = c("Mic","Ast","End","Oli","Inh","Exc","OPC","Per"))
# chrObj <- chromoInitiate(DEdf = all_DEdf, gene_col = "gene", fc_col = "avg_log2FC", p_col = "p_val_adj", celltype_col = "celltype")
# chrObj <- chromoComposition(chrObj, score_method = "hyp_padj", separate_by = "celltype", only_expr_features = T)
# youngMalePLot <- chromoCompositionPlot(chrObj, size_dot_alt = 0.8, size_dot_no = 0.4, title_xaxis = NULL, title_yaxis = NULL) +
#   ggtitle("Young male")
# 
# # Old male: DS vs CT
# all_DEdf <- do.call(rbind, DElist_old_men)
# all_DEdf$celltype <- factor(all_DEdf$celltype, levels = c("Oli","Ast","Mic","Inh","Exc","End","OPC","Per"))
# chrObj <- chromoInitiate(DEdf = all_DEdf, gene_col = "gene", fc_col = "avg_log2FC", p_col = "p_val_adj", celltype_col = "celltype")
# chrObj <- chromoComposition(chrObj, score_method = "hyp_padj", separate_by = "celltype", only_expr_features = T)
# oldMalePlot <- chromoCompositionPlot(chrObj, size_dot_alt = 0.8, size_dot_no = 0.4, title_xaxis = NULL, title_yaxis = NULL) +
#   ggtitle("Old male")
# 
# # Old female: DS vs CT
# all_DEdf <- do.call(rbind, DElist_old_women)
# all_DEdf$celltype <- factor(all_DEdf$celltype, levels = c("Mic","Oli","Inh","Ast","Exc","End","OPC","Per"))
# chrObj <- chromoInitiate(DEdf = all_DEdf, gene_col = "gene", fc_col = "avg_log2FC", p_col = "p_val_adj", celltype_col = "celltype")
# chrObj <- chromoComposition(chrObj, score_method = "hyp_padj", separate_by = "celltype", only_expr_features = T)
# oldFemalePlot <- chromoCompositionPlot(chrObj, size_dot_alt = 0.8, size_dot_no = 0.4, title_xaxis = NULL, title_yaxis = NULL) +
#   ggtitle("Old female")

combined_plots <- plot_grid(plotlist = list(youngPLot, oldPlot), ncol=2, nrow=1)
ggsave("results/sc/chromo/composition_allCellTypes.png", combined_plots, height = 4, width = 8)

# Total num of DEGs
# Young
all_DEdf <- do.call(rbind, DElist_young)
all_DEdf <- all_DEdf %>% filter(gene != "XIST")
all_DEdf$celltype <- factor(all_DEdf$celltype, levels = c("Mic","Ast","End","Oli","Inh","Exc","OPC","Per"))
chrObj <- chromoInitiate(DEdf = all_DEdf, gene_col = "gene", fc_col = "avg_log2FC", p_col = "p_val_adj", celltype_col = "celltype")
chrObj <- chromoComposition(chrObj, score_method = "n_DEGs", separate_by = "celltype", only_expr_features = T)
youngPLot <- chromoCompositionPlot(chrObj, size_dot_alt = 0.8, size_dot_no = 0.4, title_xaxis = NULL, title_yaxis = NULL) +
  ggtitle("Young")

# Old
all_DEdf <- do.call(rbind, DElist_old)
all_DEdf$celltype <- factor(all_DEdf$celltype, levels = c("Oli","Mic","Ast","Inh","Exc","End","OPC","Per"))
chrObj <- chromoInitiate(DEdf = all_DEdf, gene_col = "gene", fc_col = "avg_log2FC", p_col = "p_val_adj", celltype_col = "celltype")
chrObj <- chromoComposition(chrObj, score_method = "n_DEGs", separate_by = "celltype", only_expr_features = T)
oldPlot <- chromoCompositionPlot(chrObj, size_dot_alt = 0.8, size_dot_no = 0.4, title_xaxis = NULL, title_yaxis = NULL) +
  ggtitle("Old")

combined_plots <- plot_grid(plotlist = list(youngPLot, oldPlot), ncol=2, nrow=1)
ggsave("results/sc/chromo/composition_allCellTypes_numDEGs.png", combined_plots, height = 4, width = 8)


#######################
###       GSEA      ###
#######################

DElist_combined <- list(
  #"DElist_all_ages" = DElist_all_ages, 
  "DElist_young" = DElist_young,
  #"DElist_young_men" = DElist_young_men,
  "DElist_old" = DElist_old
  #"DElist_old_men" = DElist_old_men,
  #"DElist_old_women" = DElist_old_women
)

###### GO
gsea_list <- list()

for(i in names(DElist_combined)){
  for(j in names(DElist_combined[[i]])){
    gene_ranks <- DElist_combined[[i]][[j]]$avg_log2FC
    names(gene_ranks) <- DElist_combined[[i]][[j]]$gene
    gene_ranks <- sort(gene_ranks, decreasing = TRUE)
    
    aux <- gseGO(
      geneList      = gene_ranks,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL", 
      ont           = "BP",
      minGSSize     = 10,
      maxGSSize     = 500,
      pvalueCutoff  = 0.05,
      verbose       = FALSE
    )
    gsea_list[[i]][[j]] <- aux
  }
}

# Venn
for(i in names(gsea_list)){
  aux <- lapply(gsea_list[[i]], function(df) df@result$ID)
  file_name <- paste0("results/sc/chromo/upset_", i, ".png")
  
  png(file_name, width = 8000, height = 6400, res = 1200)
  print(upset(fromList(aux), 
              nsets = length(aux), 
              sets = c("Per","End","OPC","Oli","Inh","Exc","Mic","Ast"), 
              keep.order = TRUE))
  dev.off()
}

# test <- pairwise_termsim(gsea_list[["DElist_old"]][["Oli"]])
# test2 <- as.data.frame(test@termsim)

###### KEGG
kegg_list <- list()

for(i in names(DElist_combined)){
  for(j in names(DElist_combined[[i]])){
    gene_ranks <- DElist_combined[[i]][[j]]$avg_log2FC
    names(gene_ranks) <- DElist_combined[[i]][[j]]$gene
    gene_ranks <- sort(gene_ranks, decreasing = TRUE)
    
    # Changing SYMBOL to ENTREZID
    aux <- bitr(names(gene_ranks), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # get entrezid
    aux <- aux[!duplicated(aux$SYMBOL),] # remove one of multiplates
    gene_ranks <- gene_ranks[names(gene_ranks) %in% aux$SYMBOL] # filtering
    aux <- aux %>% arrange(match(SYMBOL, names(gene_ranks))) # reordeing
    names(gene_ranks) <- aux$ENTREZID # renaming
    
    aux <- gseKEGG(
      geneList      = gene_ranks,
      organism     = 'hsa', 
      minGSSize     = 10,
      maxGSSize     = 500,
      pvalueCutoff  = 0.05,
      verbose       = FALSE
    )
    kegg_list[[i]][[j]] <- aux
  }
}

# Venn
for(i in names(kegg_list)){
  aux <- lapply(kegg_list[[i]], function(df) df@result$ID)
  file_name <- paste0("results/sc/chromo/upset_kegg_", i, ".png")
  
  png(file_name, width = 8000, height = 6400, res = 1200)
  print(upset(fromList(aux), 
              nsets = length(aux), 
              sets = c("Per","End","OPC","Oli","Inh","Exc","Mic","Ast"), 
              keep.order = TRUE))
  dev.off()
}

# GSEA
# dummy_chrObj <- chromoInitiate(head(DElist_young[["Oli"]]), gene_col = "gene", fc_col = "avg_log2FC", p_col = "p_val_adj")
# 
# aux <- ORA_result@result %>% mutate(score = -log10(p.adjust), ONTOLOGY = "BP")
# wgcna_results[[gsub("_.*", "", i)]]$ORA[[gsub(".*_", "", i)]] <- aux
# 
# dummy_chrObj@ora[["UP"]][["1"]] <- aux # just a dummy location
# 
# ora_plot <- chromoORAPlot(dummy_chrObj, DEG_type = "UP", cluster = "1") +
#   ggtitle(i) +
#   theme(legend.position = "none")

# Young End
gsea_list[["DElist_young"]][["End"]] <- pairwise_termsim(gsea_list[["DElist_young"]][["End"]])
tree_plot <- treeplot(gsea_list[["DElist_young"]][["End"]], color = "NES", fontsize = 4.2) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "NES") +
  theme(legend.position = "none") +
  coord_cartesian(clip = "off")
tree_plot$data$x <- tree_plot$data$x * 0.08
for (i in c(3,5)) tree_plot$layers[[i]]$data$x <- tree_plot$layers[[i]]$data$x * 0.45 # right text
for (i in c(3,5)) tree_plot$layers[[i]]$aes_params$fontface <- "bold"
tree_plot$layers[[4]]$data$x    <- tree_plot$layers[[4]]$data$x * 0.45 # right lines
tree_plot$layers[[4]]$data$xend <- tree_plot$layers[[4]]$data$xend * 0.45 # right lines
ggsave("results/sc/chromo/gsea_End_young.png", plot = tree_plot, width = 15, height = 5.2)

# Young Mic
gsea_list[["DElist_young"]][["Mic"]] <- pairwise_termsim(gsea_list[["DElist_young"]][["Mic"]])
tree_plot <- treeplot(gsea_list[["DElist_young"]][["Mic"]], color = "NES", fontsize = 4.2) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "NES") +
  theme(legend.position = "none") +
  coord_cartesian(clip = "off")
tree_plot$data$x <- tree_plot$data$x * 0.08
for (i in c(3,5)) tree_plot$layers[[i]]$data$x <- tree_plot$layers[[i]]$data$x * 0.45 # right text
for (i in c(3,5)) tree_plot$layers[[i]]$aes_params$fontface <- "bold"
tree_plot$layers[[4]]$data$x    <- tree_plot$layers[[4]]$data$x * 0.45 # right lines
tree_plot$layers[[4]]$data$xend <- tree_plot$layers[[4]]$data$xend * 0.45 # right lines
ggsave("results/sc/chromo/gsea_Mic_young.png", plot = tree_plot, width = 15, height = 5)

# Young Oli ##### esse ta bugado!!!! why???
gsea_list[["DElist_young"]][["Oli"]] <- pairwise_termsim(gsea_list[["DElist_young"]][["Oli"]])
tree_plot <- treeplot(gsea_list[["DElist_young"]][["Oli"]], color = "NES", fontsize = 4.2) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "NES") +
  theme(legend.position = "none") +
  coord_cartesian(clip = "off")
tree_plot$data$x <- tree_plot$data$x * 0.08
for (i in c(3,5)) tree_plot$layers[[i]]$data$x <- tree_plot$layers[[i]]$data$x * 0.45 # right text
for (i in c(3,5)) tree_plot$layers[[i]]$aes_params$fontface <- "bold"
tree_plot$layers[[4]]$data$x    <- tree_plot$layers[[4]]$data$x * 0.45 # right lines
tree_plot$layers[[4]]$data$xend <- tree_plot$layers[[4]]$data$xend * 0.45 # right lines
ggsave("results/sc/chromo/gsea_Oli_young.png", plot = tree_plot, width = 20, height = 8)


###################################################
###       Chromo composition by cell type       ###
###################################################

for(i in names(DElist_combined)){
  plot_list <- list()
  for(j in names(DElist_combined[[i]])){
    chrObj <- chromoInitiate(DEdf = DElist_combined[[i]][[j]], gene_col = "gene", fc_col = "avg_log2FC", p_col = "p_val_adj")
    chrObj@data$chromosome_name <- factor(chrObj@data$chromosome_name, levels = c(seq(1,22), "X", "Y", "MT"))
    chrObj <- chromoComposition(chrObj, score_method = "hyp") # hyp_padj
    if(nrow(chrObj@composition$compo_df) > 0){
      plot_list[[j]] <- chromoCompositionPlot(chrObj) + # , size_dot_alt = 0.8, size_dot_no = 0.4
        ggtitle(j)
    }
  }
  combined_plots <- plot_grid(plotlist = plot_list, ncol=2, nrow=4)
  ggsave(paste0("results/sc/chromo/composition_", i, ".png"), combined_plots, height = 14, width = 17)
}

# Individual composition plots
chrObj <- chromoInitiate(DEdf = DElist_young[["Oli"]], gene_col = "gene", fc_col = "avg_log2FC", p_col = "p_val_adj")
chrObj@data <- chrObj@data %>% filter(gene != "XIST")
chrObj@data$chromosome_name <- factor(chrObj@data$chromosome_name, levels = c(seq(1,22), "X", "Y", "MT"))
chrObj <- chromoComposition(chrObj, score_method = "hyp")
chromo_plot <- chromoCompositionPlot(chrObj)
ggsave(paste0("results/sc/chromo/composition_Oli_young.png"), chromo_plot, height = 3.5, width = 11)

chrObj <- chromoInitiate(DEdf = DElist_old[["Oli"]], gene_col = "gene", fc_col = "avg_log2FC", p_col = "p_val_adj")
chrObj@data <- chrObj@data %>% filter(gene != "XIST")
chrObj@data$chromosome_name <- factor(chrObj@data$chromosome_name, levels = c(seq(1,22), "X", "Y", "MT"))
chrObj <- chromoComposition(chrObj, score_method = "hyp")
chromo_plot <- chromoCompositionPlot(chrObj)
ggsave(paste0("results/sc/chromo/composition_Oli_old.png"), chromo_plot, height = 3.5, width = 11)

chrObj <- chromoInitiate(DEdf = DElist_young[["End"]], gene_col = "gene", fc_col = "avg_log2FC", p_col = "p_val_adj")
chrObj@data <- chrObj@data %>% filter(gene != "XIST")
chrObj@data$chromosome_name <- factor(chrObj@data$chromosome_name, levels = c(seq(1,22), "X", "Y", "MT"))
chrObj <- chromoComposition(chrObj, score_method = "hyp")
chromo_plot <- chromoCompositionPlot(chrObj)
ggsave(paste0("results/sc/chromo/composition_End_young.png"), chromo_plot, height = 3.5, width = 11)

chrObj <- chromoInitiate(DEdf = DElist_old[["End"]], gene_col = "gene", fc_col = "avg_log2FC", p_col = "p_val_adj")
chrObj@data <- chrObj@data %>% filter(gene != "XIST")
chrObj@data$chromosome_name <- factor(chrObj@data$chromosome_name, levels = c(seq(1,22), "X", "Y", "MT"))
chrObj <- chromoComposition(chrObj, score_method = "hyp")
chromo_plot <- chromoCompositionPlot(chrObj)
ggsave(paste0("results/sc/chromo/composition_End_old.png"), chromo_plot, height = 3.5, width = 11)


#################################
###       Chromo density      ###
#################################

chromo_DElist_young <- list()
for(i in names(DElist_young)){
  aux <- chromoInitiate(DElist_young[[i]], gene_col = "gene", fc_col = "avg_log2FC", p_col = "p_val_adj")
  aux <- chromoDensity(aux, DEG_type = "UP", bandwidth = 3e6, cluster_threshold = 10)
  aux <- chromoDensity(aux, DEG_type = "DOWN", bandwidth = 3e6, cluster_threshold = 10)
  chromo_DElist_young[[i]] <- aux
}

chromo_DElist_old <- list()
for(i in names(DElist_young)){
  aux <- chromoInitiate(DElist_old[[i]], gene_col = "gene", fc_col = "avg_log2FC", p_col = "p_val_adj")
  aux <- chromoDensity(aux, DEG_type = "UP", bandwidth = 3e6, cluster_threshold = 10)
  aux <- chromoDensity(aux, DEG_type = "DOWN", bandwidth = 3e6, cluster_threshold = 10)
  chromo_DElist_old[[i]] <- aux
}

# Saving suplementary table
write_xlsx(
  list(
    End_young = chromo_DElist_young[["End"]]@density[["UP"]][["DEG_clusters"]], 
    Oli_young = chromo_DElist_young[["Oli"]]@density[["UP"]][["DEG_clusters"]]
  ), 
  path = "data/sc/ST/ST2.xlsx"
)

# Plotting density
for(i in names(DElist_young)){
  young_up_plot <- ggplot() + theme_classic()
  if(nrow(chromo_DElist_young[[i]]@density[["UP"]][["DEG_clusters"]]) > 0){
    young_up_plot <- chromoDensityPlot(chromo_DElist_young[[i]], DEG_type = "UP", 
      n_top_clusters = chromo_DElist_young[[i]]@density[["UP"]][["DEG_clusters"]] %>% filter(pval < 0.05) %>% nrow())
  }
  
  young_down_plot <- ggplot() + theme_classic()
  if(nrow(chromo_DElist_young[[i]]@density[["DOWN"]][["DEG_clusters"]]) > 0){
    young_down_plot <- chromoDensityPlot(chromo_DElist_young[[i]], DEG_type = "DOWN",
      n_top_clusters = chromo_DElist_young[[i]]@density[["DOWN"]][["DEG_clusters"]] %>% filter(pval < 0.05) %>% nrow())
  }
  
  old_up_plot <- ggplot() + theme_classic()
  if(nrow(chromo_DElist_old[[i]]@density[["UP"]][["DEG_clusters"]]) > 0){
    old_up_plot <- chromoDensityPlot(chromo_DElist_old[[i]], DEG_type = "UP",
      n_top_clusters = chromo_DElist_old[[i]]@density[["UP"]][["DEG_clusters"]] %>% filter(pval < 0.05) %>% nrow())
  }
  
  old_down_plot <- ggplot() + theme_classic()
  if(nrow(chromo_DElist_old[[i]]@density[["DOWN"]][["DEG_clusters"]]) > 0){
    old_down_plot <- chromoDensityPlot(chromo_DElist_old[[i]], DEG_type = "DOWN",
      n_top_clusters = chromo_DElist_old[[i]]@density[["DOWN"]][["DEG_clusters"]] %>% filter(pval < 0.05) %>% nrow())
  }
  
  combined_plots <- plot_grid(plotlist = list(young_up_plot, young_down_plot, old_up_plot, old_down_plot), 
    ncol=2, nrow=2, labels = c("Young UP", "Young DOWN", "Old UP", "Old DOWN"), 
    label_size = 14, label_fontface = "bold", label_x = 0.3, label_y = 1.005)
  ggsave(paste0("results/sc/chromo/density_", i, ".png"), combined_plots, height = 20, width = 15)
}

# # ORA of specific clusters
# chromo_DElist_young[["End"]] <- chromoORA(chromo_DElist_young[["End"]], DEG_type = "UP", cluster = 2)
# ora_plot <- chromoORAPlot(chromo_DElist_young[["End"]], DEG_type = "UP", cluster = 2, number_of_onto = 20)
# ggsave("results/sc/chromo/ORA_End_young_U2.png", ora_plot, height = 10, width = 5)

# Density plots for specific cell types
density_plot <- chromoDensityPlot(chromo_DElist_young[["End"]], DEG_type = "UP", size_cluster_name = 8,
  n_top_clusters = chromo_DElist_young[["End"]]@density[["UP"]][["DEG_clusters"]] %>% filter(pval < 0.05) %>% nrow())
ggsave("results/sc/chromo/density_End_young_up.png", density_plot, height = 10, width = 7)

density_plot <- chromoDensityPlot(chromo_DElist_young[["Oli"]], DEG_type = "UP", size_cluster_name = 8,
  n_top_clusters = chromo_DElist_young[["Oli"]]@density[["UP"]][["DEG_clusters"]] %>% filter(pval < 0.05) %>% nrow())
ggsave("results/sc/chromo/density_Oli_young_up.png", density_plot, height = 10, width = 7)

density_plot <- chromoDensityPlot(chromo_DElist_young[["Mic"]], DEG_type = "DOWN", size_cluster_name = 8,
  n_top_clusters = chromo_DElist_young[["Mic"]]@density[["DOWN"]][["DEG_clusters"]] %>% filter(pval < 0.05) %>% nrow())
ggsave("results/sc/chromo/density_Mic_young_down.png", density_plot, height = 10, width = 7)


############################
###       cell chat      ###
############################

# Reading saved objects (running cellchat takes long!!!)
cellchat_CT_O <- readRDS(file = "data/sc/cellChat/cellchat_CT_O.rds")
cellchat_CT_Y <- readRDS(file = "data/sc/cellChat/cellchat_CT_Y.rds")
cellchat_DS_O <- readRDS(file = "data/sc/cellChat/cellchat_DS_O.rds")
cellchat_DS_Y <- readRDS(file = "data/sc/cellChat/cellchat_DS_Y.rds")

# saveRDS(cellchat_CT_O, file = "data/sc/cellChat/cellchat_CT_O.rds")
# saveRDS(cellchat_CT_Y, file = "data/sc/cellChat/cellchat_CT_Y.rds")
# saveRDS(cellchat_DS_O, file = "data/sc/cellChat/cellchat_DS_O.rds")
# saveRDS(cellchat_DS_Y, file = "data/sc/cellChat/cellchat_DS_Y.rds")
# 
# Idents(scdown) <- "type_cond"
# 
# # Control old
# scdown_CT_O <- subset(scdown, subset = condition == "CT" & age == "O")
# data.input_CT_O <- GetAssayData(scdown_CT_O, assay = "RNA", slot = "data")
# meta_CT_O <- data.frame(group = Idents(scdown_CT_O), row.names = colnames(scdown_CT_O)) %>% 
#   mutate(group = factor(group, levels = c("Per_O_CT","End_O_CT","OPC_O_CT","Oli_O_CT","Inh_O_CT","Exc_O_CT","Mic_O_CT","Ast_O_CT")))
# cellchat_CT_O <- createCellChat(object = data.input_CT_O, meta = meta_CT_O, group.by = "group")
# cellchat_CT_O@DB <- CellChatDB.human
# cellchat_CT_O <- subsetData(cellchat_CT_O)
# cellchat_CT_O <- identifyOverExpressedGenes(cellchat_CT_O)
# cellchat_CT_O <- identifyOverExpressedInteractions(cellchat_CT_O)
# cellchat_CT_O <- projectData(cellchat_CT_O, PPI.human)
# cellchat_CT_O <- computeCommunProb(cellchat_CT_O, raw.use = TRUE)
# cellchat_CT_O <- filterCommunication(cellchat_CT_O, min.cells = 10)
# cellchat_CT_O <- computeCommunProbPathway(cellchat_CT_O)
# cellchat_CT_O <- aggregateNet(cellchat_CT_O)
# cellchat_CT_O <- netAnalysis_computeCentrality(cellchat_CT_O)
# 
# # Control young
# scdown_CT_Y <- subset(scdown, subset = condition == "CT" & age == "Y")
# data.input_CT_Y <- GetAssayData(scdown_CT_Y, assay = "RNA", slot = "data")
# meta_CT_Y <- data.frame(group = Idents(scdown_CT_Y), row.names = colnames(scdown_CT_Y)) %>% 
#   mutate(group = factor(group, levels = c("Per_Y_CT","End_Y_CT","OPC_Y_CT","Oli_Y_CT","Inh_Y_CT","Exc_Y_CT","Mic_Y_CT","Ast_Y_CT")))
# cellchat_CT_Y <- createCellChat(object = data.input_CT_Y, meta = meta_CT_Y, group.by = "group")
# cellchat_CT_Y@DB <- CellChatDB.human
# cellchat_CT_Y <- subsetData(cellchat_CT_Y)
# cellchat_CT_Y <- identifyOverExpressedGenes(cellchat_CT_Y)
# cellchat_CT_Y <- identifyOverExpressedInteractions(cellchat_CT_Y)
# cellchat_CT_Y <- projectData(cellchat_CT_Y, PPI.human)
# cellchat_CT_Y <- computeCommunProb(cellchat_CT_Y, raw.use = TRUE)
# cellchat_CT_Y <- filterCommunication(cellchat_CT_Y, min.cells = 10)
# cellchat_CT_Y <- computeCommunProbPathway(cellchat_CT_Y)
# cellchat_CT_Y <- aggregateNet(cellchat_CT_Y)
# cellchat_CT_Y <- netAnalysis_computeCentrality(cellchat_CT_Y)
# 
# # DS old
# scdown_DS_O <- subset(scdown, subset = condition == "DS" & age == "O")
# data.input_DS_O <- GetAssayData(scdown_DS_O, assay = "RNA", slot = "data")
# meta_DS_O <- data.frame(group = Idents(scdown_DS_O), row.names = colnames(scdown_DS_O)) %>% 
#   mutate(group = factor(group, levels = c("Per_O_DS","End_O_DS","OPC_O_DS","Oli_O_DS","Inh_O_DS","Exc_O_DS","Mic_O_DS","Ast_O_DS")))
# cellchat_DS_O <- createCellChat(object = data.input_DS_O, meta = meta_DS_O, group.by = "group")
# cellchat_DS_O@DB <- CellChatDB.human
# cellchat_DS_O <- subsetData(cellchat_DS_O)
# cellchat_DS_O <- identifyOverExpressedGenes(cellchat_DS_O)
# cellchat_DS_O <- identifyOverExpressedInteractions(cellchat_DS_O)
# cellchat_DS_O <- projectData(cellchat_DS_O, PPI.human)
# cellchat_DS_O <- computeCommunProb(cellchat_DS_O, raw.use = TRUE)
# cellchat_DS_O <- filterCommunication(cellchat_DS_O, min.cells = 10)
# cellchat_DS_O <- computeCommunProbPathway(cellchat_DS_O)
# cellchat_DS_O <- aggregateNet(cellchat_DS_O)
# cellchat_DS_O <- netAnalysis_computeCentrality(cellchat_DS_O)
# 
# # DS young
# scdown_DS_Y <- subset(scdown, subset = condition == "DS" & age == "Y")
# data.input_DS_Y <- GetAssayData(scdown_DS_Y, assay = "RNA", slot = "data")
# meta_DS_Y <- data.frame(group = Idents(scdown_DS_Y), row.names = colnames(scdown_DS_Y)) %>% 
#   mutate(group = factor(group, levels = c("Per_Y_DS","End_Y_DS","OPC_Y_DS","Oli_Y_DS","Inh_Y_DS","Exc_Y_DS","Mic_Y_DS","Ast_Y_DS")))
# cellchat_DS_Y <- createCellChat(object = data.input_DS_Y, meta = meta_DS_Y, group.by = "group")
# cellchat_DS_Y@DB <- CellChatDB.human
# cellchat_DS_Y <- subsetData(cellchat_DS_Y)
# cellchat_DS_Y <- identifyOverExpressedGenes(cellchat_DS_Y)
# cellchat_DS_Y <- identifyOverExpressedInteractions(cellchat_DS_Y)
# cellchat_DS_Y <- projectData(cellchat_DS_Y, PPI.human)
# cellchat_DS_Y <- computeCommunProb(cellchat_DS_Y, raw.use = TRUE)
# cellchat_DS_Y <- filterCommunication(cellchat_DS_Y, min.cells = 10)
# cellchat_DS_Y <- computeCommunProbPathway(cellchat_DS_Y)
# cellchat_DS_Y <- aggregateNet(cellchat_DS_Y)
# cellchat_DS_Y <- netAnalysis_computeCentrality(cellchat_DS_Y)

# Scatter
scatter_CT_O <- netAnalysis_signalingRole_scatter(cellchat_CT_O) + ggtitle("CT_O")
scatter_CT_Y <- netAnalysis_signalingRole_scatter(cellchat_CT_Y) + ggtitle("CT_Y")
scatter_DS_O <- netAnalysis_signalingRole_scatter(cellchat_DS_O) + ggtitle("DS_O")
scatter_DS_Y <- netAnalysis_signalingRole_scatter(cellchat_DS_Y) + ggtitle("DS_Y")
png("results/sc/cellChat/scatter.png", width = 10000, height = 8000, res = 1200)
scatter_CT_Y + scatter_DS_Y + scatter_CT_O + scatter_DS_O
dev.off()

# scatter combined
data_Y_CT <- ggplot_build(scatter_CT_Y)[["plot"]][["data"]] %>% 
  mutate(x = scale(x), y = scale(y))
data_Y_DS <- ggplot_build(scatter_DS_Y)[["plot"]][["data"]] %>% 
  mutate(x = scale(x), y = scale(y))
data_O_CT <- ggplot_build(scatter_CT_O)[["plot"]][["data"]] %>% 
  mutate(x = scale(x), y = scale(y))
data_O_DS <- ggplot_build(scatter_DS_O)[["plot"]][["data"]] %>% 
  mutate(x = scale(x), y = scale(y))

data_combined <- rbind(data_Y_CT, data_Y_DS, data_O_CT, data_O_DS) %>% 
  mutate(
    celltype = gsub("_.*", "", labels),
    pair = gsub("_[^_]*$", "", labels),
    condition = gsub(".*_", "", labels),
    age = gsub(".*_(.*)_.*", "\\1", labels)
  )

arrow_data <- data_combined %>%
  group_by(pair) %>%
  filter(condition %in% c("CT", "DS")) %>%
  arrange(factor(condition, levels = c("CT", "DS"))) %>%
  summarise(x_start = x[condition == "CT"],
            y_start = y[condition == "CT"],
            x_end = x[condition == "DS"],
            y_end = y[condition == "DS"])

combined_plot <- ggplot(data_combined, aes(x = x, y = y, color = celltype)) +
  geom_point(aes(size = Count, shape = age)) +
  geom_segment(data = arrow_data,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(0.2, "cm")),
               inherit.aes = FALSE) +
  geom_text_repel(aes(label = condition, color = celltype)) +
  labs(
    x = "Sending",
    y = "Receiving"
  ) +
  #scale_color_brewer(palette = "Set1") + # yellow unreadable, but good pallet!
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#777777", "#A65628", "#F781BF")) +
  scale_shape_manual(values = c("O" = 16, "Y" = 17)) +
  theme_classic()

ggsave("results/sc/cellChat/scatter_combined.png", combined_plot, height = 5, width = 6)

# heatmap
ht_CT_O <- netVisual_heatmap(
  cellchat_CT_O, 
  color.heatmap = "Reds", 
  cluster.rows = FALSE, 
  cluster.cols = FALSE#, 
  #title.name = "O_CT"
)

ht_CT_Y <- netVisual_heatmap(
  cellchat_CT_Y, 
  color.heatmap = "Reds", 
  cluster.rows = FALSE, 
  cluster.cols = FALSE#, 
  #title.name = "Y_CT"
)

ht_DS_O <- netVisual_heatmap(
  cellchat_DS_O, 
  color.heatmap = "Reds", 
  cluster.rows = FALSE, 
  cluster.cols = FALSE#, 
  #title.name = "O_DS"
)

ht_DS_Y <- netVisual_heatmap(
  cellchat_DS_Y, 
  color.heatmap = "Reds", 
  cluster.rows = FALSE, 
  cluster.cols = FALSE#, 
  #title.name = "Y_DS"
)

png("results/sc/cellChat/heat.png", width = 12000, height = 10000, res = 1200)
pushViewport(viewport(layout = grid.layout(2, 2)))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(ht_CT_Y, newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(ht_DS_Y, newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(ht_CT_O, newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
draw(ht_DS_O, newpage = FALSE)
popViewport()

popViewport()
dev.off()

##### heatmap combined
# Young
mtx_CT_Y <- as.data.frame(ht_CT_Y@matrix) %>% 
  { rownames(.) <- gsub("_[^_]*$", "", rownames(.)); . } %>% 
  { colnames(.) <- gsub("_[^_]*$", "", colnames(.)); . } %>% 
  mutate(across(everything(), ~ replace_na(., 0)))

mtx_DS_Y <- as.data.frame(ht_DS_Y@matrix) %>% 
  { rownames(.) <- gsub("_[^_]*$", "", rownames(.)); . } %>% 
  { colnames(.) <- gsub("_[^_]*$", "", colnames(.)); . } %>% 
  mutate(across(everything(), ~ replace_na(., 0)))

mtx_young <- mtx_DS_Y - mtx_CT_Y

q_lower <- quantile(mtx_young, 0.05, na.rm = TRUE)
q_upper <- quantile(mtx_young, 0.95, na.rm = TRUE)

col_fun <- colorRamp2(c(min(mtx_young, na.rm = TRUE), 0, max(mtx_young, na.rm = TRUE)), c("#4575b4", "white", "#d73027"))
ht <- Heatmap(mtx_young, name = "Value", col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_names = TRUE, show_column_names = TRUE,
              #top_annotation = HeatmapAnnotation(ColumnSum = anno_barplot(colSums(mtx_young), border = FALSE)),
              #right_annotation = rowAnnotation(RowSum = anno_barplot(rowSums(mtx_young), border = FALSE)),
              cell_fun = function(j, i, x, y, width, height, fill) {
                val <- mtx_young[i, j]
                if(val < q_lower || val > q_upper) {
                  grid.text("*", x = x, y = y, gp = gpar(fontsize = 20, col = "black"))
                }
              },
              heatmap_legend_param = list(title = ""),
              row_names_side = "left"#,
              #column_title = "Young"
              )
ht_young <- draw(ht, heatmap_legend_side = "left")

# Old
mtx_CT_O <- as.data.frame(ht_CT_O@matrix) %>% 
  { rownames(.) <- gsub("_[^_]*$", "", rownames(.)); . } %>% 
  { colnames(.) <- gsub("_[^_]*$", "", colnames(.)); . } %>% 
  mutate(across(everything(), ~ replace_na(., 0)))

mtx_DS_O <- as.data.frame(ht_DS_O@matrix) %>% 
  { rownames(.) <- gsub("_[^_]*$", "", rownames(.)); . } %>% 
  { colnames(.) <- gsub("_[^_]*$", "", colnames(.)); . } %>% 
  mutate(across(everything(), ~ replace_na(., 0)))

mtx_old <- mtx_DS_O - mtx_CT_O

q_lower <- quantile(mtx_old, 0.05, na.rm = TRUE)
q_upper <- quantile(mtx_old, 0.95, na.rm = TRUE)

col_fun <- colorRamp2(c(min(mtx_old, na.rm = TRUE), 0, max(mtx_old, na.rm = TRUE)), c("#4575b4", "white", "#d73027"))
ht <- Heatmap(mtx_old, name = "Value", col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_names = TRUE, show_column_names = TRUE,
              #top_annotation = HeatmapAnnotation(ColumnSum = anno_barplot(colSums(mtx_old), border = FALSE)),
              #right_annotation = rowAnnotation(RowSum = anno_barplot(rowSums(mtx_old), border = FALSE)),
              cell_fun = function(j, i, x, y, width, height, fill) {
                val <- mtx_old[i, j]
                if(val < q_lower || val > q_upper) {
                  grid.text("*", x = x, y = y, gp = gpar(fontsize = 20, col = "black"))
                }
              },
              heatmap_legend_param = list(title = "")#, # "Variation in DS"
              #column_title = "Old"
              )
ht_old <- draw(ht, heatmap_legend_side = "right")

# Plotting
png("results/sc/cellChat/heat_combined.png", width = 12000, height = 5000, res = 1200)
pushViewport(viewport(layout = grid.layout(1, 2)))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(ht_young, newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(ht_old, newpage = FALSE)
popViewport()

popViewport()
dev.off()

# Bubble
bubble_CT_O <- netVisual_bubble(
  cellchat_CT_O, 
  sources.use = c("Per_O_CT", "Mic_O_CT", "End_O_CT"), 
  targets.use = c("Per_O_CT", "Mic_O_CT", "End_O_CT", "OPC_O_CT", "Exc_O_CT"), 
  remove.isolate = FALSE
) + 
  ggtitle("Old control")

bubble_DS_O <- netVisual_bubble(
  cellchat_DS_O, 
  sources.use = c("Per_O_DS", "Mic_O_DS", "End_O_DS"), 
  targets.use = c("Per_O_DS", "Mic_O_DS", "End_O_DS", "OPC_O_DS", "Exc_O_DS"), 
  remove.isolate = FALSE
) + 
  ggtitle("Old DS")

bubble_CT_Y <- netVisual_bubble(
  cellchat_CT_Y, 
  sources.use = c("Per_Y_CT", "Mic_Y_CT", "End_Y_CT", "OPC_Y_CT", "Oli_Y_CT", "Ast_Y_CT"), 
  targets.use = c("Per_Y_CT", "OPC_Y_CT", "Ast_Y_CT"), 
  remove.isolate = FALSE
) + 
  ggtitle("Young control")

bubble_DS_Y <- netVisual_bubble(
  cellchat_DS_Y, 
  sources.use = c("Per_Y_DS", "Mic_Y_DS", "End_Y_DS", "OPC_Y_DS", "Oli_Y_DS", "Ast_Y_DS"), 
  targets.use = c("Per_Y_DS", "OPC_Y_DS", "Ast_Y_DS"), 
  remove.isolate = FALSE
) + 
  ggtitle("Young DS")

combined_plots <- plot_grid(plotlist = list(bubble_CT_O, bubble_DS_O, bubble_CT_Y, bubble_DS_Y), nrow=1, ncol=4)
ggsave("results/sc/cellChat/bubble.png", combined_plots, height = 14, width = 30)

# Bubble combined
# Old
bd_O <- rbind(bubble_CT_O[["data"]], bubble_DS_O[["data"]]) %>% 
  dplyr::select(source, target, ligand, receptor, pathway_name) %>% 
  mutate(
    condition = gsub(".*_", "", source),
    interaction = paste0(ligand, " -> ", receptor),
    source = gsub("_[^_]*$", "", source),
    target = gsub("_[^_]*$", "", target),
    st = paste0(source, " -> ", target)
  ) %>% 
  filter(st %in% c("Per_O -> Per_O", "Per_O -> Mic_O", "End_O -> End_O", "End_O -> OPC_O", "End_O -> Exc_O", "End_O -> Mic_O", "Mic_O -> OPC_O")) %>% 
  drop_na() %>% 
  group_by(source, target, interaction) %>% 
  filter(n() == 1) %>% 
  ungroup() %>%
  complete(interaction, st)

mat <- bd_O %>%
  dplyr::select(interaction, st, condition) %>% 
  mutate(condition = ifelse(condition == "DS", 1, ifelse(condition == "CT", -1, NA))) %>%
  pivot_wider(names_from = st, values_from = condition) %>%
  column_to_rownames("interaction") %>%
  as.matrix()
mat[is.na(mat)] <- 0

d_rows <- dist(mat)
hc_rows <- hclust(d_rows)
row_order <- rownames(mat)[hc_rows$order]
d_cols <- dist(t(mat))
hc_cols <- hclust(d_cols)
col_order <- colnames(mat)[hc_cols$order]
bd_O$interaction <- factor(bd_O$interaction, levels = row_order)
bd_O$st <- factor(bd_O$st, levels = col_order)

old <- ggplot(bd_O, aes(x = st, y = interaction, fill = condition)) +
  geom_tile(color = "black") + 
  scale_fill_manual(
    name = "Changes",
    values = c("DS" = "#d73027", "CT" = "#4575b4"),
    labels = c("DS" = "Gained", "CT" = "Lost"),
    na.value = "white"
  ) +
  theme_classic() +
  theme(
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_blank()
  ) +
  labs(x = NULL, y = NULL)

# Young
bd_Y <- rbind(bubble_CT_Y[["data"]], bubble_DS_Y[["data"]]) %>% 
  dplyr::select(source, target, ligand, receptor, pathway_name) %>% 
  mutate(
    condition = gsub(".*_", "", source),
    interaction = paste0(ligand, " -> ", receptor),
    source = gsub("_[^_]*$", "", source),
    target = gsub("_[^_]*$", "", target),
    st = paste0(source, " -> ", target)
  ) %>% 
  filter(st %in% c("Per_Y -> Ast_Y", "End_Y -> OPC_Y", "OPC_Y -> Per_Y", "OPC_Y -> OPC_Y", 
                   "OPC_Y -> Ast_Y", "Oli_Y -> Ast_Y", "Mic_Y -> OPC_Y", "Ast_Y -> Ast_Y")) %>%
  drop_na() %>% 
  group_by(source, target, interaction) %>% 
  filter(n() == 1) %>% 
  ungroup() %>% 
  complete(interaction, st)

mat <- bd_Y %>%
  dplyr::select(interaction, st, condition) %>% 
  mutate(condition = ifelse(condition == "DS", 1, ifelse(condition == "CT", -1, NA))) %>%
  pivot_wider(names_from = st, values_from = condition) %>%
  column_to_rownames("interaction") %>%
  as.matrix()
mat[is.na(mat)] <- 0

d_rows <- dist(mat)
hc_rows <- hclust(d_rows)
row_order <- rownames(mat)[hc_rows$order]
d_cols <- dist(t(mat))
hc_cols <- hclust(d_cols)
col_order <- colnames(mat)[hc_cols$order]
bd_Y$interaction <- factor(bd_Y$interaction, levels = row_order)
bd_Y$st <- factor(bd_Y$st, levels = col_order)

young <- ggplot(bd_Y, aes(x = st, y = interaction, fill = condition)) +
  geom_tile(color = "black") + 
  scale_fill_manual(
    name = "Changes",
    values = c("DS" = "#d73027", "CT" = "#4575b4"),
    labels = c("DS" = "Gained", "CT" = "Lost"),
    na.value = "white"
  ) +
  theme_classic() +
  theme(
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_blank(),
    legend.position = "none"
  ) +
  labs(x = NULL, y = NULL)

combined_plots <- plot_grid(plotlist = list(young, old), ncol=2, nrow=1)
ggsave("results/sc/cellChat/bubble_combined.png", combined_plots, height = 10, width = 10)

######## falta o violin e depois um dotplot final pra amarrar com o chromo density que nem com o hdWGCNA
aux <- rbind(
  drop_na(bd_O),
  drop_na(bd_Y)
) %>% 
  mutate(
    age = gsub(".*_", "", source),
    source = gsub("_.*", "", source),
    target = gsub("_.*", "", target)
  )

violin_data <- rbind(
  aux %>% group_by(source, ligand, age) %>% summarise(.groups = "drop") %>% rename(celltype = source, gene = ligand),
  aux %>% group_by(target, receptor, age) %>% summarise(.groups = "drop") %>% rename(celltype = target, gene = receptor)
) %>%
  separate_rows(gene, sep = "_") %>% 
  rowwise() %>% 
  mutate(
    DEG_Y = {
      current_gene <- gene
      current_celltype <- celltype
      deg_values <- chromo_DElist_young[[current_celltype]]@data %>% 
        filter(gene == current_gene) %>% 
        pull(DEG)
      deg_values[1]
    },
    DEG_O = {
      current_gene <- gene
      current_celltype <- celltype
      deg_values <- chromo_DElist_old[[current_celltype]]@data %>% 
        filter(gene == current_gene) %>% 
        pull(DEG)
      deg_values[1]
    }
  ) %>%
  ungroup() %>% 
  filter(DEG_O != "NO" | DEG_Y != "NO") %>% 
  filter((age == "Y" & DEG_Y != "NO")|(age == "O" & DEG_O != "NO")) %>% 
  mutate(cell_age = paste0(celltype, "_", age))

# Dot plot
scdown$cell_age <- paste0(scdown$seurat_clusters, "_", scdown$age)
Idents(scdown) <- "cell_age"
sub_scdown <- subset(scdown, idents = unique(violin_data$cell_age))
sub_scdown$cell_age_cond <- paste0(sub_scdown$seurat_clusters, "_", sub_scdown$age, "_", sub_scdown$condition)
#sub_scdown$cell_age_cond <- factor(sub_scdown$cell_age_cond, levels = c(""))
Idents(sub_scdown) <- "cell_age_cond"

dot_plot <- DotPlot(sub_scdown, features = unique(violin_data$gene), group.by = "cell_age_cond") + 
  theme_classic() +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_color_gradient(low = "#4575b4", high = "#d73027")

ggsave("results/sc/cellChat/cellchat_DEGs.png", dot_plot, height = 2.5, width = 4)

# Pathways
pathways <- aux %>% 
  filter(paste0(interaction, "_", st) %in% 
           c("LAMA3 -> ITGAV_ITGB8_End_Y -> OPC_Y", "PDGFB -> PDGFRA_End_Y -> OPC_Y", "ANGPTL4 -> SDC4_Ast_Y -> Ast_Y", 
             "ANGPTL4 -> SDC2_Ast_Y -> Ast_Y", "VSIR -> IGSF11_Ast_Y -> Ast_Y",
             "LAMA5 -> ITGAV_ITGB8_End_Y -> OPC_Y", "TNR -> ITGA8_ITGB1_OPC_Y -> Per_Y", 
             "LAMA4 -> ITGAV_ITGB8_End_Y -> OPC_Y", "MIF -> CD74_CXCR4_Per_O -> Mic_O", "MIF -> CD74_CXCR4_End_O -> Mic_O")) %>% 
  group_by(interaction, pathway_name, st) %>% 
  summarise() %>% 
  arrange(pathway_name) %>% 
  rename(
    "ligand -> receptor" = "interaction",
    "pathway" = "pathway_name",
    "source -> target" = "st"
  )

# Plot
custom_theme <- ttheme_default(
  core = list(bg_params = list(fill = "white", col = "black")),
  colhead = list(bg_params = list(fill = "#cccccc", col = "black"))
)

table_plot <- tableGrob(pathways, rows = NULL, theme = custom_theme)

png("results/sc/cellChat/pathways.png", width = 3200, height = 2000, res = 600)
grid.draw(table_plot)
dev.off()


#######################
###      WGCNA      ###
#######################

scdown@meta.data$cell_sample <- paste0(scdown@meta.data$seurat_clusters, "_", scdown@meta.data$orig.ident)
counts <- scdown@assays[["RNA"]]@layers[["counts"]]

pseudobulk <- do.call(cbind, lapply(unique(scdown@meta.data$cell_sample), function(cs) {
  cell_idx <- which(scdown@meta.data$cell_sample == cs)
  rowSums(counts[, cell_idx])
}))

colnames(pseudobulk) <- unique(scdown@meta.data$cell_sample)
rownames(pseudobulk) <- rownames(scdown)

pseudobulk_cts <- as.data.frame(pseudobulk) #%>% 
  #filter(rownames(.) %in% VariableFeatures(scdown)) #### filter first based on counts and then calculate hvf na mão?

pseudobulk_meta <- data.frame(cell_sample = colnames(pseudobulk_cts)) %>% 
  mutate(
    celltype = gsub("_.*", "", cell_sample),
    sample = gsub(".*_", "", cell_sample)
  ) %>% 
  left_join(
    scdown@meta.data %>% select(`orig.ident`, condition, age) %>% group_by(`orig.ident`, condition, age) %>% summarise(),
    by = c("sample" = "orig.ident")
  ) %>% 
  column_to_rownames("cell_sample")

wgcna_results <- list()
for(i in c("End", "Oli")){ # unique(pseudobulk_meta$celltype)
  
  filt <- pseudobulk_meta$celltype == i
  sub_meta <- pseudobulk_meta[filt, ]
  sub_cts <- pseudobulk_cts[, filt]
  
  dds <- DESeqDataSetFromMatrix(countData = sub_cts, colData = sub_meta, design = ~ 1)
  dds <- DESeq(dds)
  vsd <- vst(dds, blind = FALSE)
  pseudobulk_norm <- as.data.frame(assay(vsd))
  
  cem <- cemitool(
    expr = pseudobulk_norm,
    network_type = "signed",
    tom_type = "signed",
    filter = TRUE,
    plot = FALSE, 
    verbose = FALSE
  )
  
  # Hubs
  wgcna_results[[i]]$hubs <- get_hubs(cem, 1)
  
  # Modules
  wgcna_results[[i]]$modules <- cem@module %>% 
    filter(modules != "Not.Correlated") %>% 
    mutate(
      modules = gsub("M", "", modules),
      modules = as.numeric(modules)
    ) %>% 
    arrange(modules)
  
  # Beta x R2 plot
  p <- cem@beta_r2_plot[["beta_r2_plot"]]
  ggsave(paste0("results/sc/cemitool/beta_R2_", i, ".png"), plot = p, width = 5, height = 5)
  
  # Differential module eigengene expression analysis
  datExpr <- as.data.frame(assay(vsd)) %>% 
    rownames_to_column("gene") %>% 
    filter(gene %in% wgcna_results[[i]]$modules$genes) %>% 
    arrange(match(gene, wgcna_results[[i]]$modules$genes)) %>% 
    column_to_rownames("gene") %>% 
    t() %>% 
    as.data.frame()
  
  MEs <- moduleEigengenes(datExpr, colors = wgcna_results[[i]]$modules$modules)$eigengenes %>% 
    t() %>% 
    as.data.frame()
  
  # Min max scaling (no negative values)
  MEs <- t(apply(MEs, 1, function(x) {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  })) %>% 
    as.data.frame()
  
  # Young DE analysis
  DEdf <- data.frame(module = character(0), pval = numeric(0), log2fc = numeric(0))
  for(j in rownames(MEs)){

    r1 <- as.numeric(MEs[rownames(MEs) == j, sub_meta$condition == "DS" & sub_meta$age == "Y", drop = TRUE])
    r2 <- as.numeric(MEs[rownames(MEs) == j, sub_meta$condition == "CT" & sub_meta$age == "Y", drop = TRUE])

    aux <- data.frame(
      module = j,
      pval = wilcox.test(r1, r2)[["p.value"]],
      log2fc = log2(mean(r1)/mean(r2))
      )

    DEdf <- rbind(DEdf, aux)
  }
  
  wgcna_results[[i]]$young_DEdf <- DEdf %>% mutate(pval = p.adjust(pval, method = "none")) # "BH"
  
  # Old DE analysis
  DEdf <- data.frame(module = character(0), pval = numeric(0), log2fc = numeric(0))
  for(j in rownames(MEs)){
    
    r1 <- as.numeric(MEs[rownames(MEs) == j, sub_meta$condition == "DS" & sub_meta$age == "O", drop = TRUE])
    r2 <- as.numeric(MEs[rownames(MEs) == j, sub_meta$condition == "CT" & sub_meta$age == "O", drop = TRUE])
    
    aux <- data.frame(
      module = j,
      pval = wilcox.test(r1, r2)[["p.value"]],
      log2fc = log2(mean(r1)/mean(r2))
    )
    
    DEdf <- rbind(DEdf, aux)
  }
  
  wgcna_results[[i]]$old_DEdf <- DEdf %>% mutate(pval = p.adjust(pval, method = "none")) # "BH"
  
  # Young plot
  youngPlot <- ggplot(data = wgcna_results[[i]]$young_DEdf, aes(x = log2fc, y = -log10(pval))) +
    geom_point(aes(color = ifelse(pval < 0.05 & log2fc < -1, "#2255dd",
                                  ifelse(pval < 0.05 & log2fc > 1, "#dd1133", "grey"))), 
               shape = 21, size = 10, fill = "#ffffff00", stroke = 2) +
    geom_text(aes(label = gsub("ME", "", module)), vjust = 0.5, hjust = 0.5, size = 5, color = "black") +
    geom_vline(xintercept = c(-1, 1), color = "red", linetype = "dotted") +
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dotted") +
    labs(x = "log2FC", y = "-log10(pval)") +
    theme(legend.position = "none") +
    theme_classic() +
    ggtitle(paste0("Young ", i)) +
    scale_color_identity()
  
  ggsave(paste0("results/sc/cemitool/DME_", i, "_young.png"), plot = youngPlot, width = 5, height = 5)
  
  # Old plot
  oldPlot <- ggplot(data = wgcna_results[[i]]$old_DEdf, aes(x = log2fc, y = -log10(pval))) +
    geom_point(aes(color = ifelse(pval < 0.05 & log2fc < -1, "#2255dd",
                                  ifelse(pval < 0.05 & log2fc > 1, "#dd1133", "grey"))), 
               shape = 21, size = 10, fill = "#ffffff00", stroke = 2) +
    geom_text(aes(label = gsub("ME", "", module)), vjust = 0.5, hjust = 0.5, size = 5, color = "black") +
    geom_vline(xintercept = c(-1, 1), color = "red", linetype = "dotted") +
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dotted") +
    labs(x = "log2FC", y = "-log10(pval)") +
    theme(legend.position = "none") +
    theme_classic() +
    ggtitle(paste0("Old ", i)) +
    scale_color_identity()
  
  ggsave(paste0("results/sc/cemitool/DME_", i, "_old.png"), plot = oldPlot, width = 5, height = 5)
}

# ORA of differentially expressed modules
for(i in c("Ast_10", "Ast_5", "End_18", "Mic_3", "Mic_4", "Oli_5", "Oli_3", "OPC_2", "OPC_5", "Per_12")){ # c("Ast_10", 5) number of ontologies to plot
  ORA_result <- enrichGO(
    keyType = "SYMBOL",
    gene = wgcna_results[[gsub("_.*", "", i)]]$modules %>% filter(modules == gsub(".*_", "", i)) %>% pull(genes),
    universe = rownames(scdown),
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    #pvalueCutoff = 0.05,
    #qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  dummy_chrObj <- chromoInitiate(head(DElist_young[["Oli"]]), gene_col = "gene", fc_col = "avg_log2FC", p_col = "p_val_adj")
  
  aux <- ORA_result@result %>% mutate(score = -log10(p.adjust), ONTOLOGY = "BP")
  wgcna_results[[gsub("_.*", "", i)]]$ORA[[gsub(".*_", "", i)]] <- aux
  
  dummy_chrObj@ora[["UP"]][["1"]] <- aux # just a dummy location
  
  ora_plot <- chromoORAPlot(dummy_chrObj, DEG_type = "UP", cluster = "1") +
    ggtitle(i) +
    theme(legend.position = "none")
  
  ggsave(paste0("results/sc/cemitool/ORA_", i, ".png"), plot = ora_plot, width = 5, height = 3)
}

# chromo enrichment YOUNG
for(i in setdiff(names(wgcna_results), c("Inh"))){
  
  # Modules (WGCNA)
  modules_list <- wgcna_results[[i]][["modules"]] %>% 
    group_by(modules) %>%
    summarise(genes = list(genes)) %>%
    deframe()
  
  # Clusters(chromo)
  clusters_up_list <- chromo_DElist_young[[i]]@density[["UP"]][["DEG_clusters"]] %>% 
    filter(pval < 0.05) %>% 
    dplyr::select(cluster_num, all_features) %>% 
    mutate(all_features = str_split(all_features, ";")) %>%
    deframe()
  names(clusters_up_list) <- paste0("U", names(clusters_up_list))
  
  clusters_down_list <- chromo_DElist_young[[i]]@density[["DOWN"]][["DEG_clusters"]] %>% 
    filter(pval < 0.05) %>% 
    dplyr::select(cluster_num, all_features) %>% 
    mutate(all_features = str_split(all_features, ";")) %>%
    deframe()
  names(clusters_down_list) <- paste0("D", names(clusters_down_list))
  
  clusters_young_list <- c(clusters_up_list, clusters_down_list)
  
  young_data <- do.call(rbind, unlist(lapply(names(modules_list), function(mn) {
    lapply(names(clusters_young_list), function(cn) {
      common <- intersect(modules_list[[mn]], clusters_young_list[[cn]])
      x <- length(common)
      data.frame(
        module = mn, 
        cluster = cn,
        inter = x,
        genes = paste(common, collapse = ";"),
        pval = if(x == 0) 1 else phyper(x - 1, length(modules_list[[mn]]),
                                        length(rownames(scdown)) - length(modules_list[[mn]]), 
                                        length(clusters_young_list[[cn]]), 
                                        lower.tail = FALSE
        )
      )
    })
  }), recursive = FALSE))
  
  # re-ordering columns and padj
  young_data <- young_data %>% 
    group_by(cluster) %>% 
    mutate(pval = p.adjust(pval, method = "BH")) %>% 
    ungroup() %>% 
    mutate(
      module = factor(as.numeric(module))
    )
  
  # saving data
  wgcna_results[[i]]$chromo_enrich[["young"]] <- young_data
  
  young_plot <- ggplot(young_data, aes(x = module, y = cluster)) +
    geom_point(aes(size = ifelse(pval == 1, NA, -log10(pval)), color = pval < 0.05)) +
    scale_color_manual(values = c("FALSE" = "#333333", "TRUE" = "#dd1133")) +
    scale_x_discrete(labels = function(x) paste0("M", x)) +
    theme_classic() +
    theme(
      legend.position = "none",
      panel.grid.major = element_line(color = "lightgrey"),
      panel.grid.minor = element_line(color = "lightgrey"),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_text(face = "bold", size = 12, color = "black")
    ) +
    labs(x = NULL, y = NULL) +
    ggtitle(paste0("Young ", i))
  
  ggsave(paste0("results/sc/cemitool/chromo_", i, "_young.png"), plot = young_plot, 
         width = 3 + 0.25*uniqueN(young_data$module), height = 3 + 0.2*uniqueN(young_data$cluster))
}

# chromo enrichment OLD
for(i in setdiff(names(wgcna_results), c("Inh", "OPC", "Per"))){
  
  # Modules (WGCNA)
  modules_list <- wgcna_results[[i]][["modules"]] %>% 
    group_by(modules) %>%
    summarise(genes = list(genes)) %>%
    deframe()
  
  # Clusters (chromo)
  clusters_up_list <- chromo_DElist_old[[i]]@density[["UP"]][["DEG_clusters"]] %>% 
    filter(pval < 0.05) %>% 
    dplyr::select(cluster_num, all_features) %>% 
    mutate(all_features = str_split(all_features, ";")) %>%
    deframe()
  names(clusters_up_list) <- paste0("U", names(clusters_up_list))
  
  clusters_down_list <- chromo_DElist_old[[i]]@density[["DOWN"]][["DEG_clusters"]] %>% 
    filter(pval < 0.05) %>% 
    dplyr::select(cluster_num, all_features) %>% 
    mutate(all_features = str_split(all_features, ";")) %>%
    deframe()
  names(clusters_down_list) <- paste0("D", names(clusters_down_list))
  
  clusters_old_list <- c(clusters_up_list, clusters_down_list)
  
  # Building df
  old_data <- do.call(rbind, unlist(lapply(names(modules_list), function(mn) {
    lapply(names(clusters_old_list), function(cn) {
      common <- intersect(modules_list[[mn]], clusters_old_list[[cn]])
      x <- length(common)
      data.frame(
        module = mn, 
        cluster = cn,
        inter = x,
        genes = paste(common, collapse = ";"),
        pval = if(x == 0) 1 else phyper(x - 1, length(modules_list[[mn]]),
                                        length(rownames(scdown)) - length(modules_list[[mn]]), 
                                        length(clusters_old_list[[cn]]), 
                                        lower.tail = FALSE
        )
      )
    })
  }), recursive = FALSE))
  
  # re-ordering columns and padj
  old_data <- old_data %>% 
    group_by(cluster) %>% 
    mutate(pval = p.adjust(pval, method = "BH")) %>% 
    ungroup() %>% 
    mutate(
      module = factor(as.numeric(module))
    )
  
  # saving data
  wgcna_results[[i]]$chromo_enrich[["old"]] <- old_data
  
  # Plot
  old_plot <- ggplot(old_data, aes(x = module, y = cluster)) +
    geom_point(aes(size = ifelse(pval == 1, NA, -log10(pval)), color = pval < 0.05)) +
    scale_color_manual(values = c("FALSE" = "#333333", "TRUE" = "#dd1133")) +
    scale_x_discrete(labels = function(x) paste0("M", x)) +
    theme_classic() +
    theme(
      legend.position = "none",
      panel.grid.major = element_line(color = "lightgrey"),
      panel.grid.minor = element_line(color = "lightgrey"),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_text(face = "bold", size = 12, color = "black")
    ) +
    labs(x = NULL, y = NULL) +
    ggtitle(paste0("Old ", i))
  
  ggsave(paste0("results/sc/cemitool/chromo_", i, "_old.png"), plot = old_plot, 
         width = 3 + 0.25*uniqueN(old_data$module), height = 3 + 0.2*uniqueN(old_data$cluster))
}

# Zoom
for(i in c("End_UP_1", "Oli_UP_1", "Mic_DOWN_3")){ # just young #### "Mic_DOWN_2", "Oli_UP_1", "Oli_UP_2", 
  
  zoom_plot <- chromoZoom(
    chromo_DElist_young[[gsub("_.*", "", i)]], 
    DEG_type = gsub("^[^_]*_([^_]+)_.*", "\\1", i),
    cluster = as.numeric(gsub(".*_", "", i)),
    size_gene_name = 4
  )
  
  ggsave(paste0("results/sc/cemitool/zoom_", i, ".png"), plot = zoom_plot, width = 9, height = 3, dpi = 900)
}

############ specific modules ORA
# MODULE 18 Endothelial
ORA_result <- enrichGO(
  keyType = "SYMBOL",
  gene = wgcna_results[["End"]]$modules %>% filter(modules == 18) %>% pull(genes),
  universe = rownames(scdown),
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.05,
  #qvalueCutoff = 0.05,
  readable = TRUE
)

ORA_result <- pairwise_termsim(ORA_result)

tree_plot <- treeplot(ORA_result, fontsize = 4.2, showCategory = nrow(ORA_result@result %>% filter(p.adjust < 0.05)), nCluster = nrow(ORA_result@termsim) %/% 10) +
  #scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + # name = "NES"
  #theme(legend.position = "none") +
  coord_cartesian(clip = "off")

tree_plot$data$x <- tree_plot$data$x * 0.08
for (i in c(3,5)) tree_plot$layers[[i]]$data$x <- tree_plot$layers[[i]]$data$x * 0.45 # right text
for (i in c(3,5)) tree_plot$layers[[i]]$aes_params$fontface <- "bold"
tree_plot$layers[[4]]$data$x    <- tree_plot$layers[[4]]$data$x * 0.45 # right lines
tree_plot$layers[[4]]$data$xend <- tree_plot$layers[[4]]$data$xend * 0.45 # right lines
ggsave("results/sc/cemitool/tree_End_18.png", plot = tree_plot, width = 15, height = 16)

# Circles ORA
aux <- as.data.frame(tree_plot[["data"]]) %>% 
  filter(complete.cases(label)) %>% 
  mutate(group = gsub(".*_", "", group)) %>% 
  rename("pval" = "color") %>% 
  group_by(group) %>% 
  summarise(
    n_onto = n(),
    metap = sumlog(pval)$p,
    color = -log10(metap)
  ) %>% 
  mutate(x = 0, y = 0)

for(i in unique(aux$group)){
  
  circle_plot <- ggplot(data = aux %>% filter(group == i), aes(x = x, y = y)) +
    geom_point(aes(fill = color), color = "#118844", shape = 21, size = 30, stroke = 2) +
    scale_fill_gradient(low = "#ffaaaa", high = "#ff4444", limits = c(min(aux$color), max(aux$color))) +
    geom_text(aes(label = n_onto), vjust = 0.5, hjust = 0.5, size = 10, color = "black") +
    coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1), clip = "off") +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0)
    )
  
  ggsave(paste0("results/sc/cemitool/End_18_group", i, ".png"), plot = circle_plot, width = 4, height = 4, dpi = 900, bg = "transparent")
}

# net plot

# grouped_ont_list <- list(
#   "ion" = list("ont" = c("GO:1990169", "GO:0010273", "GO:0097501", "GO:0061687", "GO:0071294", 
#                          "GO:0071280", "GO:0046688", "GO:0006882", "GO:0071276", "GO:0010035", 
#                          "GO:0010043", "GO:0046686", "GO:0046592", "GO:0009636", "GO:0098754", 
#                          "GO:0071241", "GO:0071248", "GO:0010038"), wid = 10, hei = 10),
#   "oxi" = list("ont" = c("GO:0030182", "GO:0016049", "GO:0051341"), wid = 5, hei = 5),
#   "neu" = list("ont" = c("GO:0021762", "GO:0001504"), wid = 5, hei = 5),
#   "mus" = list("ont" = c("GO:0048660", "GO:0048659", "GO:0033002"), wid = 5, hei = 5),
#   "acd" = list("ont" = c("GO:0015849", "GO:0046942", "GO:1901571"), wid = 5, hei = 5)
# )

aux <- as.data.frame(tree_plot[["data"]]) %>% 
  filter(complete.cases(label))

for(i in unique(aux$group)){
  sub_ora <- ORA_result
  sub_ora@result <- sub_ora@result %>% 
    filter(Description %in% c(aux %>% filter(group == i) %>% pull(label)))
  
  net_plot <- cnetplot(sub_ora, showCategory = nrow(sub_ora@result)) +
    theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    )
  
  ggsave(paste0("results/sc/cemitool/net_End_18_", i, ".png"), plot = net_plot, width = 8, height = 8)
}

# test <- intersect(
#   gsea_list$DElist_young$End@result$ID,
#   wgcna_results[["End"]][["ORA"]][["18"]] %>% filter(p.adjust < 0.05) %>% pull(ID)
# )

# MODULE 5 Oligodentrocyte
ORA_result <- enrichGO(
  keyType = "SYMBOL",
  gene = wgcna_results[["Oli"]]$modules %>% filter(modules == 5) %>% pull(genes),
  universe = rownames(scdown),
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.05,
  #qvalueCutoff = 0.05,
  readable = TRUE
)
ORA_result <- pairwise_termsim(ORA_result)
tree_plot <- treeplot(ORA_result, fontsize = 4.2) + # color = "NES"
  #scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + # name = "NES"
  #theme(legend.position = "none") +
  coord_cartesian(clip = "off")
tree_plot$data$x <- tree_plot$data$x * 0.08
for (i in c(3,5)) tree_plot$layers[[i]]$data$x <- tree_plot$layers[[i]]$data$x * 0.45 # right text
for (i in c(3,5)) tree_plot$layers[[i]]$aes_params$fontface <- "bold"
tree_plot$layers[[4]]$data$x    <- tree_plot$layers[[4]]$data$x * 0.45 # right lines
tree_plot$layers[[4]]$data$xend <- tree_plot$layers[[4]]$data$xend * 0.45 # right lines
ggsave("results/sc/cemitool/tree_Oli_5.png", plot = tree_plot, width = 15, height = 8)

# net plot
grouped_ont_list <- list(
  "ret" = list("ont" = c("GO:0050678", "GO:0050673", "GO:0001936", 
                         "GO:0045820", "GO:0060221"), wid = 5, hei = 5),
  "ske" = list("ont" = c("GO:0007519", "GO:0060538", "GO:0035914"), wid = 5, hei = 5),
  "che" = list("ont" = c("GO:0035767", "GO:2001028", "GO:2001026"), wid = 5, hei = 5),
  "ion" = list("ont" = c("GO:0046686", "GO:0071276", "GO:0071248", 
                         "GO:0046147", "GO:0045926", "GO:1990169", 
                         "GO:0010273", "GO:0097501", "GO:0061687"), wid = 7, hei = 7),
  "imm" = list("ont" = c("GO:0002294", "GO:0042093", "GO:0002287", 
                         "GO:0002293", "GO:0002292", "GO:0046641", 
                         "GO:0035710", "GO:0072538", "GO:0072539", 
                         "GO:0072540"), wid = 8, hei = 8)
)

for(i in names(grouped_ont_list)){
  sub_ora <- ORA_result
  sub_ora@result <- sub_ora@result %>% 
    filter(ID %in% grouped_ont_list[[i]][["ont"]])
  net_plot <- cnetplot(sub_ora, showCategory = nrow(sub_ora@result)) +
    theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    )
  ggsave(paste0("results/sc/cemitool/net_Oli_5_", i, ".png"), plot = net_plot, width = grouped_ont_list[[i]][["wid"]], height = grouped_ont_list[[i]][["hei"]])
}


####################################
###       Expression plots       ###
####################################

##### Understanding cell composition results
scdown$cell_age_cond <- paste0(scdown$seurat_clusters, "_", scdown$age, "_", scdown$condition)
Idents(scdown) <- "cell_age_cond"

list_genes_to_plot <- list(
  "Mic" = c("CSF1R", "CSF1", "IL34", "TGFB1", "SPI1", "TNF", "IL1B"),
  "End" = c("VEGFA", "VEGFB", "VEGFC", "KDR", "FGF2", "ANGPT1", "ANGPT2", "TIE2", "NOTCH1", "DLL4"),
  "Ast" = c("STAT3", "SOX9", "BMP4", "EGFR", "NFIA", "NFIB"),
  "Inh" = c("ASCL1", "DLX1", "DLX2", "LHX6", "NEUROD1", "BDNF", "GDNF")
)
cells_to_plot <- c()

for(i in names(list_genes_to_plot)){
  dot_plot <- DotPlot(scdown, features = list_genes_to_plot[[i]], group.by = "cell_age_cond") + # idents = cells_to_plot,
    theme_classic() +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme(
      axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    scale_color_gradient(low = "#4575b4", high = "#d73027")
  
  ggsave(paste0("results/sc/cell_composition_", i, ".png"), dot_plot, height = 6, width = 0.5*length(list_genes_to_plot[[i]]))
}

# mitochondria genes Dot plot
scdown$cell_age_cond <- paste0(scdown$seurat_clusters, "_", scdown$age, "_", scdown$condition)

mt_genes <- c(
  "MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP8", "MT-ATP6", "MT-CO3",
  "MT-ND3", "MT-ND4L", "MT-ND4", "MT-ND5", "MT-ND6", "MT-CYB"
)

# Young
sub_scdown <- subset(scdown, subset = age == "Y" & seurat_clusters %in% c("Oli", "Mic", "Ast", "Exc", "Inh"))
Idents(sub_scdown) <- "cell_age_cond"

dot_plot <- DotPlot(sub_scdown, features = mt_genes, group.by = "cell_age_cond") + 
  theme_classic() +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_color_gradient(low = "#4575b4", high = "#d73027")

ggsave("results/sc/chromo/mitochondrial_young.png", dot_plot, height = 3, width = 4)
rm(sub_scdown)

# old
sub_scdown <- subset(scdown, subset = age == "O" & seurat_clusters %in% c("Mic", "End", "Inh"))
Idents(sub_scdown) <- "cell_age_cond"

dot_plot <- DotPlot(sub_scdown, features = mt_genes, group.by = "cell_age_cond") + 
  theme_classic() +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_color_gradient(low = "#4575b4", high = "#d73027")

ggsave("results/sc/chromo/mitochondrial_old.png", dot_plot, height = 2, width = 4)
rm(sub_scdown)

# Regulation

mt_regulation_genes <- c(
  "PPARGC1A", "PPARGC1B", "NRF1", "GABPA", "TFAM", "TFB1M", "TFB2M",
  "POLRMT", "SIRT1", "SIRT3", "AMPK", "PRKAA1", "PRKAA2", "PRKAB1", "PRKAB2",
  "PRKAG1", "PRKAG2", "PRKAG3", "ESRRA", "ESRRB", "ESRRG",
  "PPARA", "PPARD", "PPARG", "FOXO1", "FOXO3", "CREB1", "ATF2", "PGC1A",
  "UCP2", "UCP3", "MTOR", "AKT1", "INSR", "LEPR", "HIF1A", "RXRA", "NFE2L2"
)


oxi_stress <- c("SOD1", "SOD2", "CAT", "GPX1", "GPX4", "PRDX1", "PRDX2", "TXN", "TXNRD1",
                "NFE2L2", "KEAP1", "FOXO1", "FOXO3", "ATF4", "HIF1A",
                "GSR", "GCLC", "GCLM", "GSS", "GSTP1", "GSTT1", "GSTK1",
                "NQO1", "IDH1", "IDH2", "G6PD", "PGD", "ME1",
                "HMOX1", "HSPA1A", "HSP90AA1", "UBB", "PSMA1", "SQSTM1")

metallo_set <- c("MT1A", "MT1E", "MT1F", "MT1G", "MT1H", "MT1M", "MT1X",
                 "MT2A", "MT3", "MT4")


##################################
###       Expression UMAP      ###
##################################
expr_matrix <- GetAssayData(scdown, assay = "RNA", layer = "data")

# Mitochondrial gene set
MTset <- GeneSet(mt_genes, setName = "MTset")
cells_rankings <- AUCell_buildRankings(expr_matrix, plotStats = FALSE, verbose = FALSE)
cells_AUC_ast <- AUCell_calcAUC(MTset, cells_rankings, aucMaxRank = nrow(expr_matrix)*0.05)
auc_scores_ast <- as.data.frame(t(getAUC(cells_AUC_ast)))
scdown <- AddMetaData(scdown, auc_scores_ast, col.name = "MTset")

# Metallothionein gene set
metalloSet <- GeneSet(metallo_set, setName = "metalloSet")
cells_rankings <- AUCell_buildRankings(expr_matrix, plotStats = FALSE, verbose = FALSE)
cells_AUC_ast <- AUCell_calcAUC(metalloSet, cells_rankings, aucMaxRank = nrow(expr_matrix)*0.05)
auc_scores_ast <- as.data.frame(t(getAUC(cells_AUC_ast)))
scdown <- AddMetaData(scdown, auc_scores_ast, col.name = "metalloSet")

# End module 18 gene set
EndM18 <- GeneSet(wgcna_results$End$modules %>% filter(modules == 18) %>% pull(genes), setName = "EndM18")
cells_rankings <- AUCell_buildRankings(expr_matrix, plotStats = FALSE, verbose = FALSE)
cells_AUC_ast <- AUCell_calcAUC(EndM18, cells_rankings, aucMaxRank = nrow(expr_matrix)*0.05)
auc_scores_ast <- as.data.frame(t(getAUC(cells_AUC_ast)))
scdown <- AddMetaData(scdown, auc_scores_ast, col.name = "EndM18")

# subseting
for(j in c("MTset", "metalloSet")){
  for(i in c("Y", "O")){
    sub_scdown <- subset(scdown, subset = age == i)
    
    umap_feature <- FeaturePlot(
      sub_scdown,
      features = j,
      pt.size = 0.3
    ) + 
      scale_color_gradient(low = "white", high = "red") + 
      facet_grid(.~sub_scdown$condition) +
      labs(
        title = NULL
      ) +
      theme(
        strip.background = element_rect(fill = "white", color = "white"),  
        strip.text = element_text(color = "black", face = "bold", size = 23),   # Text color and style
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.line = element_blank(),
        # axis.title = element_blank(),
        legend.position = "none" # no legend
      )
    ggsave(paste0("results/sc/feature_plot_", j, "_", i, ".png"), plot = umap_feature, height = 5, width = 7, limitsize = FALSE, dpi = 600)
  }
}
rm(sub_scdown)

# Dot plot

# young and old
scdown$age_cond <- paste0(scdown$age, "_", scdown$condition)
sub_scdown <- subset(scdown, subset = age %in% c("Y", "O"))
Idents(sub_scdown) <- "age_cond"

dot_plot <- DotPlot(sub_scdown, features = c("MTset", "metalloSet", "EndM18"), group.by = "age_cond") + 
  theme_classic() +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_color_gradient(low = "#4575b4", high = "#d73027")

ggsave("results/sc/dotplot_gene_sets.png", dot_plot, height = 2, width = 2)

# young
sub_scdown <- subset(scdown, subset = age %in% c("Y"))
Idents(sub_scdown) <- "age_cond"

dot_plot <- DotPlot(sub_scdown, features = c("MTset", "metalloSet", "EndM18"), group.by = "age_cond") + 
  theme_classic() +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_color_gradient(low = "#4575b4", high = "#d73027")

ggsave("results/sc/dotplot_gene_sets_young.png", dot_plot, height = 1.5, width = 2)

# old
sub_scdown <- subset(scdown, subset = age %in% c("O"))
Idents(sub_scdown) <- "age_cond"

dot_plot <- DotPlot(sub_scdown, features = c("MTset", "metalloSet", "EndM18"), group.by = "age_cond") + 
  theme_classic() +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_color_gradient(low = "#4575b4", high = "#d73027")

ggsave("results/sc/dotplot_gene_sets_old.png", dot_plot, height = 1.5, width = 2)

# End M18 genes
sub_scdown <- subset(scdown, subset = age %in% c("Y", "O"))
dot_plot <- DotPlot(sub_scdown, features = wgcna_results$End$modules %>% filter(modules == 18) %>% arrange(genes) %>% pull(genes), group.by = "age_cond") + 
  theme_classic() +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_color_gradient(low = "#4575b4", high = "#d73027")

ggsave("results/sc/dotplot_End_M18.png", dot_plot, height = 2, width = 10, dpi = 900)

rm(sub_scdown)

EndM18 <- wgcna_results[["End"]]$modules %>% filter(modules == 18) %>% pull(genes)


#########################
###     DEG tile      ###
#########################

# young
all_DEdf <- do.call(rbind, DElist_young)
chrObj <- chromoInitiate(DEdf = all_DEdf, gene_col = "gene", fc_col = "avg_log2FC", p_col = "p_val_adj", celltype_col = "celltype")
ydf <- chrObj@data %>% mutate(age = "Y") %>% dplyr::select(gene, DEG, celltype, age)

# Old
all_DEdf <- do.call(rbind, DElist_old)
chrObj <- chromoInitiate(DEdf = all_DEdf, gene_col = "gene", fc_col = "avg_log2FC", p_col = "p_val_adj", celltype_col = "celltype")
odf <- chrObj@data %>% mutate(age = "O") %>% dplyr::select(gene, DEG, celltype, age)

deg_tile_data <- rbind(ydf, odf) %>% 
  filter(
    !duplicated(interaction(gene, celltype, age))#,
    #DEG != "NO"
  ) %>% 
  mutate(cell_age = paste0(celltype, "_", age))

deg_tile <- function(
    tile_df = deg_tile_data,
    features
){
  tile_df <- tile_df %>% 
    filter(gene %in% features)
  
  tile_plot <- ggplot(tile_df, aes(x = gene, y = cell_age, fill = DEG)) +
    geom_tile(color = "black") + 
    scale_fill_manual(
      name = "Changes",
      values = c("UP"   = "#d73027",
                 "DOWN" = "#4575b4",
                 "NO"   = "white"),
      na.value = "white"
    ) +
    theme_classic() +
    theme(
      axis.ticks.x    = element_blank(),
      axis.ticks.y    = element_blank(),
      axis.text.x     = element_text(angle = 45, hjust = 1),
      axis.line       = element_blank(),
      legend.position = "none"
    ) +
    labs(x = NULL, y = NULL)
  
  col_counts <- tile_df %>%
    filter(DEG %in% c("UP", "DOWN")) %>%
    dplyr::count(gene, DEG) %>%
    tidyr::complete(
      gene = features,
      DEG  = c("UP", "DOWN"),
      fill = list(n = 0)
    )
  
  top_plot <- ggplot(col_counts, aes(x = gene, y = n, fill = DEG)) +
    geom_col() +
    scale_fill_manual(values = c("UP" = "#d73027", "DOWN" = "#4575b4")) +
    theme_classic() +
    theme(
      axis.title      = element_blank(),
      axis.ticks      = element_blank(),
      axis.text       = element_blank(),
      axis.line       = element_blank(),
      legend.position = "none"
    )
    # scale_y_continuous(expand = c(0,0), breaks = pretty_breaks()) +
    # theme(
    #   axis.title      = element_blank(),
    #   axis.ticks.x    = element_blank(),
    #   axis.text.x     = element_blank(),
    #   axis.ticks.y    = element_line(linewidth = 0.2),
    #   axis.text.y     = element_text(size = 6),
    #   axis.line.x     = element_blank(),
    #   axis.line.y     = element_line(linewidth = 0.2, color = "#333333"),
    #   legend.position = "none"
    # )
  
  row_counts <- tile_df %>%
    filter(DEG %in% c("UP", "DOWN")) %>%
    dplyr::count(cell_age, DEG) %>%
    tidyr::complete(
      cell_age = unique(tile_df$cell_age),
      DEG      = c("UP", "DOWN"),
      fill     = list(n = 0)
    )
  
  right_plot <- ggplot(row_counts, aes(y = cell_age, x = n, fill = DEG)) +
    geom_col() +
    scale_fill_manual(values = c("UP" = "#d73027", "DOWN" = "#4575b4")) +
    theme_classic() +
    theme(
      axis.title      = element_blank(),
      axis.ticks      = element_blank(),
      axis.text       = element_blank(),
      axis.line       = element_blank(),
      legend.position = "none"
    )
  
  combined <- top_plot + plot_spacer() +
    tile_plot + right_plot +
    plot_layout(
      ncol    = 2, 
      nrow    = 2,
      heights = c(1, 5),
      widths  = c(10, 1)
    ) & 
    theme(
      plot.margin = margin(0, 0, 0, 0, unit = "cm")
    ) #&
    #coord_cartesian(clip = "off")
  
  return(combined)
}

tile_plot <- deg_tile(features = EndM18)
ggsave("results/sc/tile_deg_EndM18.png", tile_plot, height = 3.5, width = 10)

# Chromatin modifying genes
chromatin_gene_families <- list(
  histone_genes = c(
    "H2AC1", "H2AC4", "H2AC6", "H2AC7", "H2AC8", "H2AC11", "H2AC12", "H2AC13", "H2AC14", "H2AC15", "H2AC16", "H2AC17", "H2AC18", "H2AC19", "H2AC20", "H2AC21", "H2AC25",
    "H2AZ1", "H2AZ2", "H2AX", "H2AJ", "H2AB1", "H2AB2", "H2AB3", "H2AP", "H2AL1Q", "H2AL3", "MACROH2A1", "MACROH2A2",
    "H2BC1", "H2BC3", "H2BC4", "H2BC5", "H2BC6", "H2BC7", "H2BC8", "H2BC9", "H2BC10", "H2BC11", "H2BC12", "H2BC13", "H2BC14", "H2BC15", "H2BC16", "H2BC17", "H2BC18", "H2BC19", "H2BC20", "H2BC21", "H2BC22",
    "H3C1", "H3C2", "H3C3", "H3C4", "H3C5", "H3C6", "H3C7", "H3C8", "H3C9", "H3C10", "H3C11", "H3C12", "H3C13", "H3C14", "H3C15", "H3C16", "H3C17", "H3C18", "H3C19", "H3C20", "H3C21", "H3C22", "H3C23",
    "H3-3A", "H3-3B", "CENPA",
    "H4C1", "H4C2", "H4C3", "H4C4", "H4C5", "H4C6", "H4C8", "H4C9", "H4C11", "H4C12",
    "H1F0", "H1FX", "H1FOO", "HIST1H1A", "HIST1H1B", "HIST1H1C", "HIST1H1D", "HIST1H1E", "HIST2H1A", "HIST2H1C"
  ),
  hat_genes = c("KAT2A", "KAT2B", "EP300", "CREBBP", "KAT5", "KAT6A", "KAT6B", "KAT7", "KAT8"),
  hmt_genes = c("EZH2", "EZH1", "SUV39H1", "SUV39H2", "SETD1A", "SETD1B", "SETD2", "NSD1", "NSD2", "KMT2A", "KMT2B", "KMT2C", "KMT2D", "PRMT1", "PRMT5"),
  hdac_genes = c("HDAC1", "HDAC2", "HDAC3", "HDAC4", "HDAC5", "HDAC6", "HDAC7", "HDAC8", "HDAC9", "HDAC10", "HDAC11", "SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7"),
  kdm_genes = c("KDM1A", "KDM1B", "KDM2A", "KDM2B", "KDM3A", "KDM3B", "KDM4A", "KDM4B", "KDM4C", "KDM5A", "KDM5B", "KDM5C", "KDM5D", "KDM6A", "KDM6B", "JMJD1C", "JMJD6"),
  swi_snf_genes = c("SMARCA4", "SMARCB1", "SMARCC1", "SMARCC2", "ARID1A", "ARID1B", "ARID2"),
  iswi_genes = c("SMARCA1", "SMARCA5", "BAZ1A", "BAZ1B"),
  chd_genes = c("CHD1", "CHD2", "CHD3", "CHD4", "CHD5", "CHD6", "CHD7", "CHD8", "CHD9"),
  ino80_genes = c("INO80", "SRCAP", "ACTL6A", "RUVBL1", "RUVBL2", "EP400"),
  bromo_genes = c("BRD2", "BRD3", "BRD4", "BRDT", "TRIM28", "ATAD2"),
  chromo_genes = c("CBX1", "CBX2", "CBX3", "CBX4", "CBX5", "CBX6", "CBX7", "CBX8"),
  phd_genes = c("KDM5A", "KDM5B", "KDM5C", "ING1", "ING2", "PHF8", "BPTF"),
  tudor_genes = c("L3MBTL1", "TP53BP1", "SPIN1", "JMJD2A"),
  dnmt_genes = c("DNMT1", "DNMT3A", "DNMT3B"),
  tet_genes = c("TET1", "TET2", "TET3")
)

altered_chromatin_genes <- c("BRD3", "CHD1", "CHD5", "CHD7", "CBX7", "HDAC4", "SMARCA1", "KDM6B", "JMJD1C", 
                             "TET1", "H2AX", "H2BC15", "H2BC18", "H2BC4", "H3-3A", "H3-3B", "H4C3")

# for(i in names(chromatin_gene_families)){
#   tile_plot <- deg_tile(features = chromatin_gene_families[[i]])
#   ggsave(paste0("results/sc/chromatin/tile_", i, ".png"), tile_plot, height = 3.5, width = 2 + 0.25*length(chromatin_gene_families[[i]]))
# }

tile_plot <- deg_tile(features = altered_chromatin_genes)
ggsave("results/sc/chromatin/tile_altered_chromatin_genes.png", tile_plot, height = 3.5, width = 4)


##############################
###      session info      ###
##############################

session_info()

