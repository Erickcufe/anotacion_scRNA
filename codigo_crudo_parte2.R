library(Seurat)
library(dplyr)
library(SingleCellExperiment)


# Recuerda descargar primero el dataset del link
so <- readRDS("alsaigh_part2.rds")

# Lo convertimos a sce

sce <- as.SingleCellExperiment(so, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>%
  mutate_if(is.character, as.factor) %>%
  DataFrame(row.names = colnames(sce))

# Definimos resolucion
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)



### #### ### ENFOQUE AUTOMATICO

# Para anotar
# cargar datos de referencia con anotaciones de Ensembl.
library(celldex)
ref.data <- HumanPrimaryCellAtlasData(ensembl=TRUE)

library(biomaRt)
library(SummarizedExperiment)

# Seleccione la base de datos Ensembl
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Extraer los ID de Ensembl del objeto SummarizedExperiment
ensembl_ids <- rownames(ref.data)

# Obtener los símbolos genéticos correspondientes a los ID de Ensembl
annotations <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                     filters = "ensembl_gene_id", values = ensembl_ids, mart = ensembl)

# Fusionar anotaciones con rowData del objeto SummarizedExperiment
rowData(ref.data)$geneSymbol <- annotations$hgnc_symbol[match(ensembl_ids, annotations$ensembl_gene_id)]
rownames(ref.data) <- rowData(ref.data)$geneSymbol


##### #### #### Predicción

library(SingleR)
library(scater)

pred <- SingleR(test = sce, 
                ref = ref.data,
                labels = ref.data$label.main,
                assay.type.test=1,
                BPPARAM= BiocParallel::MulticoreParam(6)) # 6 CPUs.

table(pred$labels)
colnames(pred)

head(sort(table(pred$labels), decreasing=TRUE))


#Primer plot SingleR

SingleR::plotScoreHeatmap(pred)


# Segundo grafico sin podado

table(pred$labels)

sce$ProbLabels <- pred$labels

plotReducedDim(sce, "TSNE", colour_by="ProbLabels")

# Con podado

summary(is.na(pred$pruned.labels))


# Grafico para ver las celulas eliminadas

# plotDeltaDistribution(pred.grun, ncol = 3)

to.remove <- is.na(pred$pruned.labels)
table(Label=pred$labels, Removed=to.remove)

to.remove <- pruneScores(pred, min.diff.med=0.2)
table(Label=pred$labels, Removed=to.remove)

# Modificamos nuestro sce
sce_dummy <- sce[,to.remove]
sce_dummy$new_labels <- pred$labels[to.remove]
table(sce_dummy$new_labels)

plotReducedDim(sce_dummy, "TSNE", colour_by="new_labels")


###### ##### ##### Agregar un cluster

library(Seurat)

# Recuerda que el archivo original era un "so"

so@meta.data$Proplabels <- sce$ProbLabels

# Le transferimos las etiquetas al objeto Seurat
Idents(so) <- so$Proplabels

table(Idents(so))

# Eliminamos los tipos celulares con menos de 350 OJO solo para este ejemplo

`%!in%` <- Negate(`%in%`)

so <- subset(so, (Proplabels %!in% c("Pro-Myelocyte", "MEP", "Pro-B_cell_CD34+", "Pro-B_cell_CD34+", "HSC_-G-CSF", "Hepatocytes", "MSC","Osteoblasts","HSC_CD34+","Erythroblast", "Pre-B_cell_CD34-","GMP","Neutrophils","CMP","DC")))

table(Idents(so))

# Observamos nuestros clusters

DimPlot(so, reduction = "tsne")


####### Marcadores


# Marcadores ya reportados para celulas espumosas
# so <- readRDS( "alsaigh_PC20_part3.rds")

foamy_cells <- c("FABP5", "GPNMB", "TREM2", "APOC1", "CD9", "SPP1")

# so <- ScaleData(so, display.progress = FALSE)

cell_typeA_marker_gene_list <- list(c("FABP5", "GPNMB", "TREM2", "APOC1", "CD9", "SPP1"))
DefaultAssay(so) <- "RNA"
so <- AddModuleScore(object = so,
                     features = cell_typeA_marker_gene_list,
                     # ctrl = 1,
                     name = "espumosas")

FeaturePlot(object = so_new, features = "espumosas1", reduction = "tsne")

select_foamy <- colnames(so)[so$espumosas1 > 1.5]

new_cells <- so$Proplabels

# Hacemos un loop para sustituir el nombre de las etiquetas
for(i in 1:length(colnames(so) %in% select_foamy)){
  if((colnames(so) %in% select_foamy)[i]){
    new_cells[i] <- "Foamy_cells"
  } else{
    new_cells[i] <- so$Proplabels[i]
  }
}

so$new_cells <- new_cells
Idents(so) <- so$new_cells

DimPlot(so, reduction = "tsne")
DimPlot(so, reduction = "umap")
DimPlot(so, reduction = "umap", group.by = "region")
DimPlot(so, reduction = "umap", group.by = "sample_id")



###### ###### Parte final

df_cells_cts <- table(Idents(so), so$region) %>% as.data.frame()
colnames(df_cells_cts) <- c("cell_type", "region", "freq")

df_cells_cts <- df_cells_cts %>%
  group_by(cell_type) %>%
  mutate(prop = freq / sum(freq))

ggplot(df_cells_cts, aes(x = cell_type, y = prop, fill = region)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  theme_classic() +
  labs(x = "Cell Type",
       y = "Proportion",
       fill = "Femoral Artery Plaque") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15))


# Encontrar "nuevos" marcadores
# Cuidado con tu memoria RAM!

f.markers <- FindAllMarkers(so,
                            only.pos = FALSE)

f.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

#so <- ScaleData(so)


DoHeatmap(so, features = top10$gene)