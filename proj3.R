################################################################################ initial set up
### installing packages
install.packages('hdf5r') # required for Read10X_h5 (called by LOAD10X_Spatial)
### loading libs
library(dplyr)
library(Seurat)
library(patchwork)
library(future)
library(ggplot2)
### checking pwd and changing it if necessary
getwd()
# change if necessary
setwd('/home/augusto/Desktop/UdS/winter_22_23/SCB/project3')
################################################################################ Beginning of Script

################################################################################ Week 1
###                                              Task 0.3 Output of Space Ranger
#

###                                                   Task 1.1: Loading the data

# sample 1
mouseSample1 <- Load10X_Spatial(data.dir = 'Section_1',
                               filename = 'V1_Mouse_Brain_Sagittal_Posterior_filtered_feature_bc_matrix.h5',
                               assay = 'Spatial',
                               slice = 'Mouse_Brain_Sagittal_Posterior_Section1')
# sample 2
mouseSample2 <- Load10X_Spatial(data.dir = 'Section_2',
                                filename = 'V1_Mouse_Brain_Sagittal_Posterior_Section_2_filtered_feature_bc_matrix.h5',
                                assay = 'Spatial',
                                slice = 'MouseBrainSagittalPosteriorSection2_')

###                                       Task 1.2: Inspecting the Seurat-object
View(mouseSample1)
View(mouseSample2)
"
Discussion:
The Gene Expression Data is stored under the metadata attribute and the image 
is stored in the @images attribute.
"

###                                         Task 1.3: visualization of a feature
# plot 
macrophagesVSvasculature <- SpatialFeaturePlot(mouseSample1,  
                                               features = c('Ckb', 'Timp3'))
macrophagesVSvasculature
# saving plot and clearing it from memory
ggsave('macrophagesVSvasculature.jpg')
rm(macrophagesVSvasculature)
"
CKB is associated with microglial cells and Timp3 is associated with smooth 
muscle cells, meaning the vasculature of the brain.
"
cbVSmo <- SpatialFeaturePlot(mouseSample2,  features = c('Hoxb5', 'Crtam')) 
cbVSmo
# saving plot and clearing it from memory
ggsave('cbVSmo.jpg')
rm(cbVSmo)

"
Hoxb5 is a homebox type gene enriched in the medulla oblongata in mice. On the 
other hand Crtam is a gene associated with cell adhesssion and enriched in both
the cerebellum and bone marrow and lymphoid tissue. Furthermore, it is 
particularly associated with the cells belonging to the T-cell lineage.
"

###                                                          Task 2.1: Filtering

## nCount
# sample 1
p1 <- VlnPlot(mouseSample1, features = 'nCount_Spatial') + 
  NoLegend() + 
  ggtitle('Slice 1 Count') + 
  theme(plot.title = element_text(face = 'bold'), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())
p2 <- SpatialFeaturePlot(mouseSample1, features = 'nCount_Spatial') + 
  theme(legend.position = 'right') + 
  ggtitle('Slice 1 Spatial map') + 
  theme(plot.title = element_text(face = 'bold'))

# sample 2
p3 <- VlnPlot(mouseSample2, features = 'nCount_Spatial') + 
  NoLegend() + 
  ggtitle('Slice 2 Count') + 
  theme(plot.title = element_text(face = 'bold'), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())
p4 <- SpatialFeaturePlot(mouseSample2, features = 'nCount_Spatial') + 
  theme(legend.position = 'right') + 
  ggtitle('Slice 2 Spatial map') + 
  theme(plot.title = element_text(face = 'bold')) 

s1plots <- wrap_plots(A = p1, B = p2, design = 'AB') 
s1plots
ggsave('s1plots.png')
s2plots <- wrap_plots(A = p3, B = p4, design = 'AB') 
s2plots
ggsave('s2plots.jpg')

# clear memory
rm(p1, p2, p3, p4)

"
Conclussion:
Count attribute demonstrates many outliers in distribution. Furthermore, tissue 
samples have different levels of accesibility, which naturally implies that some
sort of data filtering must be performed.
"

## checking statistical summary before any filtering is done
summary(mouseSample1@meta.data)
"
         orig.ident   nCount_Spatial  nFeature_Spatial
 SeuratProject:3355   Min.   :   29   Min.   :  28    
                      1st Qu.:10047   1st Qu.:3510    
                      Median :15783   Median :4772    
                      Mean   :18047   Mean   :4723    
                      3rd Qu.:23493   3rd Qu.:5928    
                      Max.   :79686   Max.   :9476  
"
summary(mouseSample2@meta.data)
"
         orig.ident   nCount_Spatial  nFeature_Spatial
 SeuratProject:3289   Min.   :  143   Min.   : 114    
                      1st Qu.: 9129   1st Qu.:3306    
                      Median :14262   Median :4542    
                      Mean   :16508   Mean   :4542    
                      3rd Qu.:21392   3rd Qu.:5665    
                      Max.   :63052   Max.   :8983    
"
## establishing cut off points and subsetting according to them
# sample 1
lower_bound <- quantile(mouseSample1@meta.data$nCount_Spatial, 
                         probs = .05)[[1]]
upper_bound <- quantile(mouseSample1@meta.data$nCount_Spatial,
                         probs = .95)[[1]]
mouseSample1 <- subset(mouseSample1, 
                       nCount_Spatial >= lower_bound 
                       & nCount_Spatial <= upper_bound)
# sample 2
lower_bound2 <- quantile(mouseSample2@meta.data$nCount_Spatial, 
                         probs = .05)[[1]]
upper_bound2 <- quantile(mouseSample2@meta.data$nCount_Spatial,
                         probs = .95)[[1]]
mouseSample2 <- subset(mouseSample2, 
                       nCount_Spatial >= lower_bound2 
                       & nCount_Spatial <= upper_bound2)

# clear memory
rm(lower_bound, lower_bound2, upper_bound, upper_bound2)

## checking statistical summary after outlier exclussion
summary(mouseSample1@meta.data)
"
         orig.ident   nCount_Spatial  nFeature_Spatial
 SeuratProject:2013   Min.   : 8639   Min.   :2143    
                      1st Qu.:12414   1st Qu.:4067    
                      Median :15783   Median :4772    
                      Mean   :16255   Mean   :4741    
                      3rd Qu.:19856   3rd Qu.:5438    
                      Max.   :26086   Max.   :6952  
"
summary(mouseSample2@meta.data)
"
         orig.ident   nCount_Spatial  nFeature_Spatial
 SeuratProject:1973   Min.   : 8237   Min.   :2331    
                      1st Qu.:11108   1st Qu.:3798    
                      Median :14262   Median :4542    
                      Mean   :14717   Mean   :4494    
                      3rd Qu.:17828   3rd Qu.:5158    
                      Max.   :23916   Max.   :6654 
"

## replotting
# sample 1
p1 <- VlnPlot(mouseSample1, features = 'nCount_Spatial') + 
  NoLegend() + 
  ggtitle('Slice 1 Count') + 
  theme(plot.title = element_text(face = 'bold'), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())
p2 <- SpatialFeaturePlot(mouseSample1, features = 'nCount_Spatial') + 
  theme(legend.position = 'right') + 
  ggtitle('Sample 1 Spatial map') + 
  theme(plot.title = element_text(face = 'bold'))

# sample 2
p3 <- VlnPlot(mouseSample2, features = 'nCount_Spatial') + 
  NoLegend() + 
  ggtitle('Slice 2 Count') + 
  theme(plot.title = element_text(face = 'bold'), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())
p4 <- SpatialFeaturePlot(mouseSample2, features = 'nCount_Spatial') + 
  theme(legend.position = 'right') + 
  ggtitle('Sample 2 Spatial map') + 
  theme(plot.title = element_text(face = 'bold')) 

s1plots_2 <- wrap_plots(A = p1, B = p2, design = 'AB') 
s1_filtered <- s1plots / s1plots_2
ggsave('s1_filtered.jpg')
s2plots_2 <- wrap_plots(A = p3, B = p4, design = 'AB') 
s2_filtered <- s2plots / s2plots_2
ggsave('s2_filtered.jpg')
s1_filtered
s2_filtered

# clear memory
rm(s1plots, s1plots_2, s2plots, s2plots_2, p1, p2, p3, p4, s1_filtered, s2_filtered)

"
Where as before we used a cutoff of 20% on both ends, this time a 5% cutoff 
point was used for the sake of retaining information through out the slide.
The biggest reason for the the differences in expression is the anatomy of the 
brain which makes RNA more accessible and abundant at certain location when
compared to others.
"

###                                                  Task 2.2: Apply scTransform

mouseSample1 <- SCTransform(mouseSample1, assay = 'Spatial', verbose = T)
mouseSample2 <- SCTransform(mouseSample2, assay = 'Spatial', verbose = T)

"
SCTransform() replaces the function NormalizeData() as our first choice for a 
data nomralization method.
"

###                                           Task 3.1: Dimensionality reduction

## sample 1
mouseSample1 <- RunPCA(mouseSample1, assay = 'SCT', verbose = T)
# elbow plot to look for point of diminishing return in terms of variance/dim
ePlot1 <- ElbowPlot(mouseSample1)
ePlot1
ggsave('ePlot1.jpg')
"From this graph, it appears that 20 dimensions caputres most of the variation"
# creating UMAP
mouseSample1 <- RunUMAP(mouseSample1, dims = 1:20, reduction = 'pca')

## sample 2
mouseSample2 <- RunPCA(mouseSample2, assay = 'SCT', verbose = T)
# elbow plot to look for point of diminishing return in terms of variance/dim
ePlot2 <- ElbowPlot(mouseSample2)
ePlot2
ggsave('ePlot2.jpg')
"From this graph, it appears that 20 dimensions caputres most of the variation"
# creating UMAP
mouseSample2 <- RunUMAP(mouseSample2, dims = 1:20, reduction = 'pca')

rm(ePlot1, ePlot2)
###                                                         Task 3.2: Clustering

## clustering
# mouse 1
mouseSample1 <- FindNeighbors(mouseSample1, reduction = "pca", dims = 1:20)
mouseSample1 <- FindClusters(mouseSample1, resolution = .5)
# mouse 2
mouseSample2 <- FindNeighbors(mouseSample2, reduction = "pca", dims = 1:20)
mouseSample2 <- FindClusters(mouseSample2, resolution = .5)
"
FindNeighbors() ir run to make a KNN graph which is later used by FindClusters()
to cluster the cells.

A resolution of .05 was used as the recommended vals for resolution for sets of 
~3k cells is [.4 - 1.2], as mentioned in the clustering documentation for 
Seurat.
"
## Plotting results

# mouse 1
s1p1 <- DimPlot(mouseSample1, reduction = 'pca', label = T, label.size = 2) + ggtitle('Sample 1')
s1p2 <- SpatialDimPlot(mouseSample1, label = T) + NoLegend()
s1p3 <- SpatialDimPlot(mouseSample1, 
                     cells.highlight = CellsByIdentities(object = mouseSample1,
                                                         idents = c(1,2,3,4,5,6,7,8,9,10)
                                                         ), 
                     facet.highlight = T, 
                     ncol = 5
                     )
# Plotting PCA and Clusters on Spatial Plot
integratedClusterSample1 <- wrap_plots(s1p1 | s1p2) 
integratedClusterSample1
# save plot and remove from ram
ggsave('integratedClusterSample1.jpg')
rm(integratedClusterSample1, s1p1, s1p2)
# Individual Spatial Plots by Cluster
individualclusters1 <- s1p3
individualclusters1
# save plot and remove from ram
ggsave('individualclusters1.jpg')
rm(individualclusters1, s1p3)

# mouse 2
s2p1 <- DimPlot(mouseSample2, reduction = 'pca', label = T, label.size = 2)
s2p2 <- SpatialDimPlot(mouseSample2, label = T) + NoLegend()
s2p3 <- SpatialDimPlot(mouseSample2, 
                       cells.highlight = CellsByIdentities(object = mouseSample2,
                                                           idents = c(1,2,3,4,5,6,7,8,9,10)
                       ), 
                       facet.highlight = T, 
                       ncol = 5
)
# Plotting PCA and Clusters on Spatial Plot
integratedClusterSample2 <- wrap_plots(s2p1 | s2p2) 
integratedClusterSample2
# save plot and remove from ram
ggsave('integratedClusterSample1.jpg')
rm(integratedClusterSample2, s2p1, s2p2)
# Individual Spatial Plots by Cluster
individualclusters1 <- s2p3
individualclusters1
# save plot and remove from ram
ggsave('individualclusters1.jpg')
rm(individualclusters2, s2p3)


################################################################################ Week 3

###                               Task 4.1: DEG analysis based on the clustering


## mouse 1
# setting up meta data for clustering
cluster_id = 0
sample_id = 1
idCol <- as.data.frame(c(cluster_id, cluster_id, cluster_id))
colnames(idCol) <- 'Cluster_ID'
sampleCol <- as.data.frame(c(sample_id, sample_id, sample_id))
colnames(sampleCol) <- 'Sample_ID'
####################################      Cluster 1
# Find enriched features for cluster 1
mouseSample1.markers1 <- FindMarkers(mouseSample1, 
                                           ident.1 = 1, 
                                           ident.2 = c(2,3,4,5,6,7,8,9,10),
                                           min.pct = .25)
mouseSample1.markers1 <- head(mouseSample1.markers1, 3) %>% as.data.frame()
cluster_id <- cluster_id + 1
idCol$Cluster_ID <- cluster_id
mouseSample1.markers1 <- cbind(mouseSample1.markers1, idCol, sampleCol)
markers <- mouseSample1.markers1

####################################      Cluster 2
# Find enriched features for cluster 2
mouseSample1.markers2 <- FindMarkers(mouseSample1, 
                                    ident.1 = 2, 
                                    ident.2 = c(1,3,4,5,6,7,8,9,10),
                                    min.pct = .25)
mouseSample1.markers2 <- head(mouseSample1.markers2, 3) %>% as.data.frame()
cluster_id <- cluster_id + 1
idCol$Cluster_ID <- cluster_id
mouseSample1.markers2 <- cbind(mouseSample1.markers2, idCol, sampleCol)
markers <- rbind(markers, mouseSample1.markers2)
####################################      Cluster 3
# Find enriched features for cluster 3
mouseSample1.markers3 <- FindMarkers(mouseSample1, 
                                    ident.1 = 3, 
                                    ident.2 = c(1,2,4,5,6,7,8,9,10),
                                    min.pct = .25)
mouseSample1.markers3 <- head(mouseSample1.markers3, 3) %>% as.data.frame()
cluster_id <- cluster_id + 1
idCol$Cluster_ID <- cluster_id
mouseSample1.markers3 <- cbind(mouseSample1.markers3, idCol, sampleCol)
markers <- rbind(markers, mouseSample1.markers3)
####################################      Cluster 4
# Find enriched features for cluster 4
mouseSample1.markers4 <- FindMarkers(mouseSample1, 
                                    ident.1 = 4, 
                                    ident.2 = c(1,2,3,5,6,7,8,9,10),
                                    min.pct = .25)
mouseSample1.markers4 <- head(mouseSample1.markers4, 3) %>% as.data.frame()
cluster_id <- cluster_id + 1
idCol$Cluster_ID <- cluster_id
mouseSample1.markers4 <- cbind(mouseSample1.markers4, idCol, sampleCol)
markers <- rbind(markers, mouseSample1.markers4)

####################################      Cluster 5
# Find enriched features for cluster 5
mouseSample1.markers5 <- FindMarkers(mouseSample1, 
                                    ident.1 = 5, 
                                    ident.2 = c(1,2,3,4,6,7,8,9,10),
                                    min.pct = .25)
mouseSample1.markers5 <- head(mouseSample1.markers5, 3) %>% as.data.frame()
cluster_id <- cluster_id + 1
idCol$Cluster_ID <- cluster_id
mouseSample1.markers5 <- cbind(mouseSample1.markers5, idCol, sampleCol)
markers <- rbind(markers, mouseSample1.markers5)

####################################      Cluster 6
# Find enriched features for cluster 6
mouseSample1.markers6 <- FindMarkers(mouseSample1, 
                                    ident.1 = 1, 
                                    ident.2 = c(1,2,3,4,5,7,8,9,10),
                                    min.pct = .25)
mouseSample1.markers6 <- head(mouseSample1.markers6, 3) %>% as.data.frame()
cluster_id <- cluster_id + 1
idCol$Cluster_ID <- cluster_id
mouseSample1.markers6 <- cbind(mouseSample1.markers6, idCol, sampleCol)
markers <- rbind(markers, mouseSample1.markers6)

####################################      Cluster 7
# Find enriched features for cluster 7
mouseSample1.markers7 <- FindMarkers(mouseSample1, 
                                    ident.1 = 7, 
                                    ident.2 = c(1,2,3,4,5,6,8,9,10),
                                    min.pct = .25)
mouseSample1.markers7 <- head(mouseSample1.markers7, 3) %>% as.data.frame()
cluster_id <- cluster_id + 1
idCol$Cluster_ID <- cluster_id
mouseSample1.markers7 <- cbind(mouseSample1.markers7, idCol, sampleCol)
markers <- rbind(markers, mouseSample1.markers7)

####################################      Cluster 8
# Find enriched features for cluster 8
mouseSample1.markers8 <- FindMarkers(mouseSample1, 
                                    ident.1 = 8, 
                                    ident.2 = c(1,2,3,4,5,6,7,9,10),
                                    min.pct = .25)
mouseSample1.markers8 <- head(mouseSample1.markers8, 3) %>% as.data.frame()
cluster_id <- cluster_id + 1
idCol$Cluster_ID <- cluster_id
mouseSample1.markers8 <- cbind(mouseSample1.markers8, idCol, sampleCol)
markers <- rbind(markers, mouseSample1.markers8)

####################################      Cluster 9
# Find enriched features for cluster 9
mouseSample1.markers9 <- FindMarkers(mouseSample1, 
                                    ident.1 = 9, 
                                    ident.2 = c(1,2,3,4,5,6,7,8,10),
                                    min.pct = .25)
mouseSample1.markers9 <- head(mouseSample1.markers9, 3) %>% as.data.frame()
cluster_id <- cluster_id + 1
idCol$Cluster_ID <- cluster_id
mouseSample1.markers9 <- cbind(mouseSample1.markers9, idCol, sampleCol)
markers <- rbind(markers, mouseSample1.markers9)

####################################      Cluster 10
# Find enriched features for cluster 10
mouseSample1.markers10 <- FindMarkers(mouseSample1, 
                                    ident.1 = 10, 
                                    ident.2 = c(1,2,3,4,5,6,7,8,9),
                                    min.pct = .25)
mouseSample1.markers10 <- head(mouseSample1.markers10, 3) %>% as.data.frame()
cluster_id <- cluster_id + 1
idCol$Cluster_ID <- cluster_id
mouseSample1.markers10 <- cbind(mouseSample1.markers10, idCol, sampleCol)
markers <- rbind(markers, mouseSample1.markers10)

## mouse 2
# setting up meta data for clustering
cluster_id = 0
sample_id = 2
sampleCol$Sample_ID <- sample_id
####################################      Cluster 1
# Find enriched features for cluster 1
mouseSample2.markers1 <- FindMarkers(mouseSample2, 
                                    ident.1 = 1, 
                                    ident.2 = c(2,3,4,5,6,7,8,9,10),
                                    min.pct = .25)
mouseSample2.markers1 <- head(mouseSample2.markers1, 3) %>% as.data.frame()
cluster_id <- cluster_id + 1
idCol$Cluster_ID <- cluster_id
mouseSample2.markers1 <- cbind(mouseSample2.markers1, idCol, sampleCol)
markers2 <- mouseSample2.markers1
####################################      Cluster 2
# Find enriched features for cluster 2
mouseSample2.markers2 <- FindMarkers(mouseSample2, 
                                    ident.1 = 2, 
                                    ident.2 = c(1,3,4,5,6,7,8,9,10),
                                    min.pct = .25)
mouseSample2.markers2 <- head(mouseSample2.markers2, 3) %>% as.data.frame()
cluster_id <- cluster_id + 1
idCol$Cluster_ID <- cluster_id
mouseSample2.markers2 <- cbind(mouseSample2.markers2, idCol, sampleCol)
markers2 <- rbind(markers2, mouseSample2.markers2)
####################################      Cluster 3
# Find enriched features for cluster 3
mouseSample2.markers3 <- FindMarkers(mouseSample2, 
                                    ident.1 = 3, 
                                    ident.2 = c(1,2,4,5,6,7,8,9,10),
                                    min.pct = .25)
mouseSample2.markers3 <- head(mouseSample2.markers3, 3) %>% as.data.frame()
cluster_id <- cluster_id + 1
idCol$Cluster_ID <- cluster_id
mouseSample2.markers3 <- cbind(mouseSample2.markers3, idCol, sampleCol)
markers2 <- rbind(markers2, mouseSample2.markers3)

####################################      Cluster 4
# Find enriched features for cluster 4
mouseSample2.markers4 <- FindMarkers(mouseSample2, 
                                    ident.1 = 4, 
                                    ident.2 = c(1,2,3,5,6,7,8,9,10),
                                    min.pct = .25)
mouseSample2.markers4 <- head(mouseSample2.markers4, 3) %>% as.data.frame()
cluster_id <- cluster_id + 1
idCol$Cluster_ID <- cluster_id
mouseSample2.markers4 <- cbind(mouseSample2.markers4, idCol, sampleCol)
markers2 <- rbind(markers2, mouseSample2.markers4)

####################################      Cluster 5
# Find enriched features for cluster 5
mouseSample2.markers5 <- FindMarkers(mouseSample2, 
                                    ident.1 = 5, 
                                    ident.2 = c(1,2,3,4,6,7,8,9,10),
                                    min.pct = .25)
mouseSample2.markers5 <- head(mouseSample2.markers5, 3) %>% as.data.frame()
cluster_id <- cluster_id + 1
idCol$Cluster_ID <- cluster_id
mouseSample2.markers5 <- cbind(mouseSample2.markers5, idCol, sampleCol)
markers2 <- rbind(markers2, mouseSample2.markers5)

####################################      Cluster 6
# Find enriched features for cluster 6
mouseSample2.markers6 <- FindMarkers(mouseSample2, 
                                    ident.1 = 1, 
                                    ident.2 = c(1,2,3,4,5,7,8,9,10),
                                    min.pct = .25)
mouseSample2.markers6 <- head(mouseSample2.markers6, 3) %>% as.data.frame()
cluster_id <- cluster_id + 1
idCol$Cluster_ID <- cluster_id
mouseSample2.markers6 <- cbind(mouseSample2.markers6, idCol, sampleCol)
markers2 <- rbind(markers2, mouseSample2.markers6)

####################################      Cluster 7
# Find enriched features for cluster 7
mouseSample2.markers7 <- FindMarkers(mouseSample2, 
                                    ident.1 = 7, 
                                    ident.2 = c(1,2,3,4,5,6,8,9,10),
                                    min.pct = .25)
mouseSample2.markers7 <- head(mouseSample2.markers7, 3) %>% as.data.frame()
cluster_id <- cluster_id + 1
idCol$Cluster_ID <- cluster_id
mouseSample2.markers7 <- cbind(mouseSample2.markers7, idCol, sampleCol)
markers2 <- rbind(markers2, mouseSample2.markers7)

####################################      Cluster 8
# Find enriched features for cluster 8
mouseSample2.markers8 <- FindMarkers(mouseSample2, 
                                    ident.1 = 8, 
                                    ident.2 = c(1,2,3,4,5,6,7,9,10),
                                    min.pct = .25)
mouseSample2.markers8 <- head(mouseSample2.markers8, 3) %>% as.data.frame()
cluster_id <- cluster_id + 1
idCol$Cluster_ID <- cluster_id
mouseSample2.markers8 <- cbind(mouseSample2.markers8, idCol, sampleCol)
markers2 <- rbind(markers2, mouseSample2.markers8)

####################################      Cluster 9
# Find enriched features for cluster 9
mouseSample2.markers9 <- FindMarkers(mouseSample2, 
                                    ident.1 = 9, 
                                    ident.2 = c(1,2,3,4,5,6,7,8,10),
                                    min.pct = .25)
mouseSample2.markers9 <- head(mouseSample2.markers9, 3) %>% as.data.frame()
cluster_id <- cluster_id + 1
idCol$Cluster_ID <- cluster_id
mouseSample2.markers9 <- cbind(mouseSample2.markers9, idCol, sampleCol)
markers2 <- rbind(markers2, mouseSample2.markers9)

####################################      Cluster 10
# Find enriched features for cluster 10
mouseSample2.markers10 <- FindMarkers(mouseSample2, 
                                    ident.1 = 10, 
                                    ident.2 = c(1,2,3,4,5,6,7,8,9),
                                    min.pct = .25)
mouseSample2.markers10 <- head(mouseSample2.markers10, 3) %>% as.data.frame()
cluster_id <- cluster_id + 1
idCol$Cluster_ID <- cluster_id
mouseSample2.markers10 <- cbind(mouseSample2.markers10, idCol, sampleCol)
markers2 <- rbind(markers2, mouseSample2.markers10)
# clear namespace
rm(sampleCol, idCol, cluster_id, sample_id)
###                       Task 4.2: DEG analysis based on the spatial patterning

## Pre-Processing
# mouse 1
future('multiprocess', workers = 7) # parallelization aid
# note that we went for top 100 genes on both slices because of hardware limitations
mouseSample1 <- FindSpatiallyVariableFeatures(mouseSample1, 
                                              assay = 'SCT',
                                              features = VariableFeatures(mouseSample1)[1:100],
                                              selection.method = 'markvariogram')
SpatiallyVariableFeatures(mouseSample1, selection.method = 'markvariogram')
top.features <- head(SpatiallyVariableFeatures(mouseSample1, 
                                               selection.method = 'markvariogram'), 
                     3)
# mouse 2
future('multiprocess', workers = 7) # parallelization aid
mouseSample2 <- FindSpatiallyVariableFeatures(mouseSample2, 
                                              assay = 'SCT',
                                              features = VariableFeatures(mouseSample2)[1:100],
                                              selection.method = 'markvariogram')

SpatiallyVariableFeatures(mouseSample2, selection.method = 'markvariogram')
top.features2 <- head(SpatiallyVariableFeatures(mouseSample2, 
                                                selection.method = 'markvariogram'), 
                      3) 

## save, clear memory

# save seurat objs
saveRDS(mouseSample1, 'mouseSample1_605.rds')
saveRDS(mouseSample2, 'mouseSample2_606.rds')
## load if neccessary
# mouseSample1 <- readRDS('mouseSample1_605.rds')
# mouseSample1 <- readRDS('mouseSample2_606.rds')
# plotting
spatiallyVariableFeaturesS1 <- SpatialFeaturePlot(mouseSample1, 
                                                 features = top.features,
                                                 ncol = 3)
spatiallyVariableFeaturesS2 <- SpatialFeaturePlot(mouseSample2, 
                                                 features = top.features2,
                                                 ncol = 3)
spatiallyVariableFeatures <- wrap_plots(spatiallyVariableFeaturesS1 / 
                                         spatiallyVariableFeaturesS2) 

# displaying
spatiallyVariableFeaturesS1
spatiallyVariableFeaturesS2
spatiallyVariableFeatures
# saving
ggsave('spatiallyVariableFeaturesS1.jpg')
ggsave('spatiallyVariableFeaturesS2.jpg')
ggsave('spatiallyVariableFeatures.jpg')
rm(spatiallyVariableFeatureS1, 
   spatiallyVariableFeatureS2, 
   spatiallyVariableFeatures)

## comparing DEG's for individual clusters vs spatially variable ones
allMarkers <- rbind(select(markers, Cluster_ID, Sample_ID), 
                    select(markers2, Cluster_ID, Sample_ID))
allMarkers
saveRDS(allMarkers, 'allMarkers.rds')

# clear memory
rm(mouseSample1.markers1, mouseSample1.markers2, mouseSample1.markers3, mouseSample1.markers4, 
   mouseSample1.markers5, mouseSample1.markers6, mouseSample1.markers7, mouseSample1.markers8, 
   mouseSample1.markers9, mouseSample1.markers10,mouseSample2.markers1, mouseSample2.markers2, 
   mouseSample2.markers3, mouseSample2.markers4, mouseSample2.markers5, mouseSample2.markers6, 
   mouseSample2.markers7, mouseSample2.markers8, mouseSample2.markers9, mouseSample2.markers10,
   top.features, top.features2, allMarkers, markers, markers2)

## load if necessary 
# mouseSample1 <- readRDS('mouseSample1.rds')
# mouseSample2 <- readRDS('mouseSample2.rds')

###                                   Task 5.1: Merging without Batch Correction
mouseSample.merge <- merge(mouseSample1, mouseSample2)
DefaultAssay(mouseSample.merge) <- 'SCT'
VariableFeatures(mouseSample.merge) <- c(VariableFeatures(mouseSample1),
                                         VariableFeatures(mouseSample2)
                                         )
# running PCA
mouseSample.merge <- RunPCA(mouseSample.merge, assay = 'SCT', verbose = T)
# elbow plot to look for point of diminishing return in terms of variance/dim
ePlotMerge <- ElbowPlot(mouseSample.merge)
ePlotMerge
ggsave('ePlotMerge.png')
rm(ePlotMerge)
"From this graph, it appears that 20 dimensions captures most var"
# Clustering
mouseSample.merge <- FindNeighbors(mouseSample.merge, reduction = "pca", dims = 1:20)
mouseSample.merge <- FindClusters(mouseSample.merge)
# creating UMAP
mouseSample.merge <- RunUMAP(mouseSample.merge, dims = 1:20, reduction = 'pca')
# Plotting merged data set with no batch correction
mergedPlot_noBC <- DimPlot(mouseSample.merge, reduction = 'umap')
mergedPlot_noBC
ggsave('mergedPlot_noBC.png')
saveRDS(mergedPlot_noBC, 'mergedPlot_noBC.rds')
rm(mergedPlot_noBC)

# save rds object in case of crash
saveRDS(mouseSample.merge, file = 'mouseSample_merge.rds')
# clear memory
rm(mouseSample.merge)

# TODO: check elbow plots
###                                      Task 5.2: Merging with Batch Correction
mouseSample1 <- SCTransform(mouseSample1,
                            assay = 'Spatial',
                            verbose = T)
mouseSample2 <- SCTransform(mouseSample2,
                            assay = 'Spatial',
                            verbose = T)
# save
saveRDS(mouseSample1, 'mouseSample1_691.rds')
saveRDS(mouseSample2, 'mouseSample2_692.rds')

## load if necessary
# mouseSample1 <- readRDS('mouseSample1_691.rds')
# mouseSample2 <- readRDS('mouseSample2_692.rds')

# finding integration anchors and integrating by them
anchors <- FindIntegrationAnchors(object.list = c(mouseSample1, mouseSample2))
mouseMerged_noBC <- IntegrateData(anchorset = anchors)
# Save original samples and clear memory
rm(anchors)
saveRDS(mouseSample1, 'mouseSample1_691.rds')
rm(mouseSample1)
saveRDS(mouseSample2, 'mouseSample2_692.rds')
rm(mouseSample2)
saveRDS(mouseMerged_noBC, 'mouseMerged_noBC.rds')
## integrated analysis
## load if necessary
# mouseMerged_noBC <- readRDS('mouseMerged_noBC.rds')
# Pre-processing.
mouseMerged <- FindSpatiallyVariableFeatures(mouseMerged, 
                                              assay = 'Spatial',
                                              features = VariableFeatures(mouseMerged)[1:100],
                                              selection.method = 'markvariogram')
# save and load
saveRDS(mouseMerged, 'mouseMerge_BC_716.rds')
mouseMerged_BC <- readRDS('mouseMerge_BC_716.rds')

## Pre-Processing
# normalization
mouseMerged_BC <- SCTransform(mouseMerged_BC, assay = 'Spatial', verbose = F)
# running PCA
mouseMerged_BC <- RunPCA(mouseMerged_BC, assay = 'SCT', verbose = T)
# elbow plot to look for point of diminishing return in terms of variance/dim
ePlotMerge2 <- ElbowPlot(mouseMerged_BC)
ePlotMerge2
ggsave('ePlotMerge2.jpg')
rm(ePlotMerge2)
"From this graph, it appears that 20 dimensions captures most var"
# Clustering
mouseMerged_BC <- FindNeighbors(mouseMerged_BC, reduction = "pca", dims = 1:20)
mouseMerged_BC <- FindClusters(mouseMerged_BC)
# creating UMAP
mouseMerged_BC <- RunUMAP(mouseMerged_BC, dims = 1:20, reduction = 'pca')
# Plotting merged data set with batch correction
mergedPlot_BC_plot <- DimPlot(mouseMerged_BC, reduction = 'umap')
mergedPlot_BC_plot
ggsave('mergedPlot_BC.jpg')
saveRDS(mouseMerged_BC, 'mouseMerged_BC.rds')
rm(mergedPlot_BC)

# NOTES: some cell name are repeated, and this will affect the output of all 
# processes downstream of IntegrateData().

###                                         Task 5.3: Detection of Batch Effects
# TODO: hace un plot de los clusters sobre la imagen, un vln, quizas un volcano para
# comparar bc vs nobc, y obviamente los PCA y UMAps
mergedPlot_noBC <- readRDS('mouseMerged_noBC.rds')

################################################################################ Week 3

###                        Task 6.1: Automatic Annotation Using Data Integration


###                                                  Task 6.2: Manual Annotation

###                                                      Task 7.0: Deconvolution

###                                         Task 7.1: Prepare the Reference Data

###                                     Task 7.2: Select Genes for Deconvolution

###                                     Task 7.3: Create an ExpressionSet-Object

###                                     Task 7.4: Perform Deconvolution

################################################################################  Week 4

###                                              Task 9: Cell-Cell Communication

