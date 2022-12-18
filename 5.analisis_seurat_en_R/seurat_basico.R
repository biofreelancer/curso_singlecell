# BioFreelancer, 2022
# Cursos con Causa. SingleCell RNA-Seq Básico con cellranger y Seurat

# cargar pacman
library("pacman")

# cargar las otras librerias
p_load( "Seurat" )
p_load( "dplyr", "tidyr", "ggplot2", "ggsci" )
p_load( "ggpubr", "hexbin", "scales" )

p_load( "SingleR", "celldex" )
p_load( "SingleCellExperiment" )

# Definimos la ruta al directorio que contiene los archivos para la matriz de conteos cellranger
# AQUI puedes cambiar las ruta, para tus propios datos
outs_dir <- "experimento_colon/outs/filtered_feature_bc_matrix/"     # este directorio debe contener los archivos: barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz

# Leemos los datos desde 10X cellranger
cellranger_data <- Read10X( data.dir = outs_dir )   

# Creamos objeto seurat
seuobject <- CreateSeuratObject( counts = cellranger_data,
                                 project = "Mi proyecto" )

# ESTE BLOQUE SOLO EJECUTALO SI TIENES BUENA CANTIDAD DE MEMORIA RAM  (al menos 8GB RAM)
# vemos los conteos
conteos <- seuobject@assays$RNA@counts %>% 
  as.data.frame( ) 

# ya lo podemos borrar
rm( conteos )
# liberamos memoria
gc( full = TRUE )
# == FIN DE LA SECCION DE MUCHA MEMORIA

# Vemos numero de celulas
table( seuobject$orig.ident )

# Revisamos el porcentaje de conteos que corresponden a genes mitocondriales
seuobject[[ "percent.mt" ]] <- PercentageFeatureSet( object = seuobject,
                                                     pattern = "^MT-" )

# Visualizamos algunas metricas de control de calidad
VlnPlot( object = seuobject,
         features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
         ncol = 3 )

# creamos un dataframe para 
scatter_df <- data.frame( identity = Idents( seuobject ),
                          nFeature_RNA = FetchData( seuobject,
                                                    "nFeature_RNA" ),
                          nCount_RNA = FetchData( seuobject,
                                                  "nCount_RNA" ),
                          percent.mt = FetchData( seuobject,
                                                  "percent.mt" ) )


# Vemos la distribucion de celulas de acuerdo al % mitocondrial y el numero de genes por celula
ggplot( data = scatter_df,
        mapping = aes( x = nFeature_RNA,
                       y = percent.mt ) ) +
  geom_point( )

# creamos un panel para comparar metricas por celula
mitohex_1 <- ggplot( data = scatter_df,
                     mapping = aes( x = nFeature_RNA,
                                    y = percent.mt ) ) +
  geom_hex( bins = 50 )

# Visualizamos
mitohex_1 

# Marcamos los limites de lo que vamos a filtrar
mitohex_2 <- mitohex_1 +
  geom_hline( yintercept = 5,               # Se recomienda usar celulas con menos del 5% de porcentaje mitocondrial
              linetype = "dashed" ) +
  geom_vline( xintercept = c( 200, 2500 ),  # Se recomienda usar celulas entre 200 y 2,500 genes
              linetype = "dashed" )

# Vis
mitohex_2

# cambiamos color y mejoramos tema
mitohex_3 <- mitohex_2 +
  scale_fill_gradient( low = "gray95", high = "black" ) +
  labs( title = "Celulas segun el %mitocondrial y el numero de genes",
        fill = "cells" ) +
  theme_light( ) 

# Vis
mitohex_3

# Revisamos la correlacion entre el numero de genes y el numero de conteos, en cada celula
ggplot( data = scatter_df,
        mapping = aes( x = nFeature_RNA,
                       y = nCount_RNA ) ) +
  geom_point( alpha = 0.2, size = 0.3 ) +
  labs( title = "Celulas segun el numero de reads (counts) y el numero de genes (features)" ) +
  stat_cor( method = "pearson" ) +
  theme_light( )

# Procedemos a filtrar las celulas, de acuerdo a las metricas de control de calidad evaluadas arriba
seuobject <- subset( x = seuobject,
                     subset = nFeature_RNA > 200 &           # celulas con mas de 200 genes detectados
                       nFeature_RNA < 2500 &                 # celulas con menos de 2,500 genes detectados
                       percent.mt < 5 )                      # y celulas con menos del 5% de genes mitocondriales

# revisamos como queda el QC
# creamos un dataframe para 
scatter_df_2 <- data.frame( identity = Idents( seuobject ),
                            nFeature_RNA = FetchData( seuobject,
                                                      "nFeature_RNA" ),
                            nCount_RNA = FetchData( seuobject,
                                                    "nCount_RNA" ),
                            percent.mt = FetchData( seuobject,
                                                    "percent.mt" ) )

ggplot( data = scatter_df_2,
        mapping = aes( x = nFeature_RNA,
                       y = percent.mt ) ) +
  geom_hex( bins = 50 ) +
  geom_hline( yintercept = 5,               # Se recomienda usar celulas con menos del 5% de porcentaje mitocondrial
              linetype = "dashed" ) +
  geom_vline( xintercept = c( 200, 2500 ),  # Se recomienda usar celulas entre 200 y 2,500 genes
              linetype = "dashed" ) +
  scale_fill_gradient( low = "gray95", high = "black" ) +
  labs( title = "Celulas segun el %mitocondrial y el numero de genes",
        fill = "cells" ) +
  theme_light( ) 

# Revisamos el numero de celulas despues del filtrado
table( seuobject$orig.ident )

## Comenzamos a preparar datos para el clustering
## Normalizamos la data
seuobject <- NormalizeData( object = seuobject,
                            normalization.method = "LogNormalize" )


# Nos quedamos con los 2000 genes mas variables dentro del dataset
seuobject <- FindVariableFeatures( object = seuobject,
                                   selection.method = "vst",
                                   nfeatures = 2000 )

# Sacamos los nombres de los genes
all_genes <- rownames( seuobject )

# Escalamos los datos, para que los genes con mucha expresion no dominen los analisis posteriores
seuobject <- ScaleData( seuobject,
                        features = all_genes )
# Realizamos el PCA
seuobject <- RunPCA( object = seuobject )

# imprimimos los primeros 5 componentes principales
# Y senialamos los 5 primeros genes que influyen en ese componente
print( seuobject[["pca"]],
       dims = 1:5,
       nfeatures = 5 )

# HAcemos un elbow plot para ver la varianza explicada por cada PC
ElbowPlot( seuobject )

# Observamos para N componentes, la influencia de los genes
VizDimLoadings( object = seuobject,
                dims = 1:4, nfeatures = 5,
                reduction = "pca" )

# Observamos el scater de los primeros 2 componentes
DimPlot( object = seuobject,
         reduction = "pca",
         dims = c(1, 2) )

# Cramos un heatmap para el PC1
DimHeatmap( object = seuobject,
            dims = 1,
            cells = 500,
            balanced = TRUE )

# Puedes pedir varios heatmaps al mismo tiempo
DimHeatmap( object = seuobject,
            dims = 1:10,
            cells = 500,
            balanced = TRUE )

## Hacemos la identificacion de clusters por UMAP
# primero encontramos las celulas mas cercanas entre ellas
seuobject <- FindNeighbors( seuobject,
                            dims = 1:10 )

# Agrupamos los vecinos
seuobject <- FindClusters( seuobject,
                           resolution = 0.5 )

# Revisamos cuantos clusters (levels) encontro el UMAP
head( Idents( seuobject ), 5 )

# hacemos el analisis umap
seuobject <- RunUMAP( seuobject,
                      n.neighbors = 30,
                      dims = 1:10 )

# Visualizamos el umap
DimPlot( object = seuobject,
         reduction = "umap",
         pt.size = 1 )

# Visualizamos el PCA, pero ya clusterizado
DimPlot( object = seuobject,
         reduction = "pca",
         pt.size = 1 )

# Se calcula el numero de celulas por cada cluster de umap
labeled_summary1 <- data.frame( cell_type = Idents( object = seuobject ) ) %>%
  group_by( cell_type ) %>%
  summarise( n_cells = n( ) )

# how many cells per label
cell_bars1 <- ggplot( data = labeled_summary1,
                      mapping = aes( x = cell_type,
                                     y = n_cells ) ) +
  geom_col( ) +
  geom_label( mapping = aes( label = n_cells )  ) +
  labs( title = "Numero de celulas por cluster UMAP",
        caption = paste( "Total cells:",
                         sum( labeled_summary1$n_cells ) ) ) +
  theme_classic( )

# Vis
cell_bars1

# Calculamos los marcadores representativos de los clusters encontrados
seuobject.markers <- FindAllMarkers( object = seuobject,
                                     min.pct = 0.25,           # solo genes detectados en al menos el 25% de las celulas en cualquiera de las 2 poblaciones
                                     logfc.threshold = 0.25 )  # le estamos pidiendo que solo analice los genes con un FC minimo de 0.25

# Mostramos los primeros 5 marcadores (ordenados por fold change) para cada cluster
cluster_markers <- seuobject.markers %>%
  group_by(cluster) %>%
  slice_max( n = 5,                            # mostramos los 5 mas
             order_by = abs( avg_log2FC ) )    # ordenados por el valor absoluto de log2 fold change

# vemos el dataframe
cluster_markers

# Hacemos el heatmap para los primeros 5 marcadores de cada cluster
# Primero saca los genes mas sobreexpresados de cada cluter de celulas
tops <- seuobject.markers %>%
  group_by(cluster) %>%
  top_n( n = 5,
         wt = abs( avg_log2FC ) )

# Aqui ya hace el heatmap
DoHeatmap( object = seuobject,
           features = tops$gene )

# Etiquetamos los clusters UMAP con singleR
# usando como referencia celldex::HumanPrimaryCellAtlasData( )$label.main

# creamos el objeto con la referencia
ref <- celldex::HumanPrimaryCellAtlasData( )

# etiquetamos las celulas en un nuevo objeto
etiquetadas <- SingleR( test = as.SingleCellExperiment( x = seuobject ),
                        ref = ref,
                        labels = ref$label.main )

# vemos el etiquetado
etiquetadas

# Ahora le ponemos las etiquetas al objeto seurat
# label the seurat object
seuobject$singler_labels <- etiquetadas$labels

# vemos el UMAP, pero ahora con las etiquetas celulares
DimPlot( seuobject,
         reduction = "umap",
         label = TRUE,
         group.by = "singler_labels" )

# vemos el PCA, pero ahora con las etiquetas celulares
DimPlot( seuobject,
         reduction = "pca",
         label = TRUE,
         group.by = "singler_labels" )

# Contamos el numero de celulas por tipo etiquetado
seuobject$singler_labels %>% 
  table( )

# lo graficamos
labeled_summary2 <- data.frame( cell_type = seuobject$singler_labels ) %>%
  group_by( cell_type ) %>%
  summarise( n_cells = n( ) )

cell_bars2 <- ggplot( data = labeled_summary2,
                      mapping = aes( x = cell_type,
                                     y = n_cells ) ) +
  geom_col( ) +
  geom_label( mapping = aes( label = n_cells )  ) +
  labs( title = "Numero de celulas por etiqueta SingleR",
        caption = paste( "Total cells:",
                         sum( labeled_summary1$n_cells ) ) ) +
  theme_classic( )

# Vis
cell_bars2

## Vamos a comparar el clustering por umap y las identidades por singleR
# Obtenemos las columnas relevantes
cluster_compare_data <- data.frame( umap = seuobject$seurat_clusters,
                                    singler_type = seuobject$singler_labels ) %>%
  mutate( cell_tag = row.names( . ) )

# summarise the combinations
cluster_compare_summarised <- cluster_compare_data %>%
  group_by( umap,
            singler_type ) %>%
  summarise( n_cells = n( ) ) %>%
  ungroup( ) %>%
  group_by( singler_type ) %>%
  mutate( singler_label_total = sum( n_cells ),
          singler_label_percentage = n_cells / singler_label_total )

# clustering concordance plot
clustering_heatmap <- ggplot( data = cluster_compare_summarised,
                              mapping = aes( x = singler_type,
                                             y = umap,
                                             fill = singler_label_percentage,
                                             label = percent( singler_label_percentage,
                                                              accuracy = 0.1 ) ) ) +
  geom_tile( ) +
  geom_text( ) +
  scale_fill_gradient( low = "skyblue",
                       high = "limegreen",
                       limits = c( 0, 1 ),
                       label = percent ) +
  labs( title = "Subgrupos por tipo celular" ) +
  theme_light( ) +
  theme( axis.text.x = element_text( angle = 90, hjust = 1 ),
         panel.grid = element_blank( ),
         legend.position = "top" )

# Vis
clustering_heatmap

## Ahora visualizamos los marcadores para los tipos celulares encontrados
# primero reidentificamos las celulas en el objeto seurat
Idents( object = seuobject ) <- seuobject$singler_labels

# buscamos los marcadores
seuobject.markers <- FindAllMarkers( object = seuobject,
                                     min.pct = 0.25,
                                     logfc.threshold = 0.25 )  # le estamos pidiendo que solo analice los genes con un FC minimo de 0.25

# Mostramos los 5 marcadores de cada tipo celular en el experimento
seuobject.markers %>%
  group_by(cluster) %>%
  slice_max( n = 5,
             order_by = avg_log2FC )

# Heatmap para los 5 principales marcadores de cada PC
# Primero saca los diez genes mas sobreexpresados de cada cluter de celulas
tops <- seuobject.markers %>%
  group_by(cluster) %>%
  top_n( n = 5,
         wt = abs(avg_log2FC) )

# Y hacemos el heatmap para tu paper :)
DoHeatmap( object = seuobject,
           features = tops$gene,
           size = 2, angle = 90 )

# Finalmente extraemos del objeto seurat solo las celulas de un tipo en especial
# Procedemos a filtrar las celulas, de acuerdo a las metricas de control de calidad evaluadas arriba
seuobject_bcells <- subset( x = seuobject,
                            subset = singler_labels == "B_cell" )                      # extraemos los que estan etiquetados con el tipo que nos interesa

# vemos las celulas que se sacaron
# vemos el UMAP, pero ahora con las etiquetas celulares
DimPlot( seuobject_bcells,
         reduction = "umap",
         label = TRUE,
         group.by = "singler_labels" )

# vemos el UMAP, pero ahora con los clusters de umap encontrados originalmente en TODA la muestra
DimPlot( seuobject_bcells,
         reduction = "umap",
         label = TRUE,
         group.by = "seurat_clusters" )

# Para que sirve sacar este grupo celular?
# Por ejemplo, si quisieras construir la red de correlacion de SOLO los b_cells en tus experimentos...
# Podemos sacar estos conteos de solo BCells :)
conteos_bcells <- seuobject_bcells@assays$RNA@counts %>% 
  as.data.frame( ) 

# y guardamos esa tabla de conteos bcells
write.csv( x = conteos_bcells, file = "conteos_bcells_colon.csv", quote = FALSE, row.names = TRUE )

# FELICIDADES, ya sabes procesar datos singlecell RNA en Seurat!
# FIN DE CURSO! Vamonos de vacaciones :)