#Code for Batf3-dendritic cells and 4-1BB/4-1BB ligand axis are required at the effector phase within the tumor microenvironment for PD-1/PD-L1 blockade efficacy
#Code written by Jason W Shapiro

#The following is intended to allow the reader to reproduce the figures and analyses for the spatial transcriptomics results in the paper.
#Please adjust filepaths accordingly if running your own reanalysis.

#A note on patient/sample ID numbers. For clarity in the paper, some patients are numbered in a different order than their original
#sample labels. In the code, patient samples correspond to identifiers in the format B#, such as B5. The following is a summary of
#how the patient sample IDs in the Seurat objects and data files corresponds to the final identifiers from the paper:
#
#  GEO_Patient    Sample_ID (in Data)
#  Patient 1, sample 1    B5
#  Patient 1, sample 2    B18
#  Patient 2              B9
#  Patient 3              B10
#  Patient 4              B11
#  Patient 5              B12
#  Patient 6              B13 (B12.1 in images)
#  Patient 7              B14
#  Patient 8              B15
#  Patient 9              B16
#  Patient 10, sample 1   B1
#  Patient 10, sample 2   B17


#Load libraries and prepare Giotto paths
library(Seurat)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tibble)
library(data.table)
library(ComplexHeatmap)
library(readxl)
library(Giotto)   #Note: individual Giotto installations might require separate source() commands to specify paths to specific functions.
library(edgeR)
library(DESeq2)
library(cowplot)
library(grid)
library(patchwork)
library(ggrepel)
library(tidyverse)
library(limma)
library(gt)
library(DT) 
library(GSEABase)
library(GSVA)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)

#Note: Giotto requires a set of instructions to run properly. The most important element is the python_path
instrs = createGiottoInstructions(save_dir = getwd(),
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  python_path = '/Users/jason/Software/miniconda3/envs/giotto/bin/pythonw'
)

#Set working directory
workdir="/Users/jason/Documents/Biocore/Gajewski-Lab/Andrea_visium/Data/"
setwd(workdir)

#Load data for each sample. Include original 2
dataB1 = Load10X_Spatial('./spaceranger_AZ_B1/outs', slice='B1')
dataB1$orig.ident='B1'
dataB1 = NormalizeData(dataB1)
dataB1 = ScaleData(dataB1)

dataB5 = Load10X_Spatial('./spaceranger_AZ_slide2_A//outs', slice='B5')
dataB5$orig.ident ='B5'
dataB5 = NormalizeData(dataB5)
dataB5 = ScaleData(dataB5)

dataB9 = Load10X_Spatial('./spaceranger_AZ_B9/outs', slice='B9')
dataB9$orig.ident='B9'
dataB9 = NormalizeData(dataB9)
dataB9 = ScaleData(dataB9)

dataB10 = Load10X_Spatial('./spaceranger_AZ_B10/outs', slice='B10')
dataB10$orig.ident='B10'
dataB10 = NormalizeData(dataB10)
dataB10 = ScaleData(dataB10)

dataB11 = Load10X_Spatial('./spaceranger_AZ_B11/outs', slice='B11')
dataB11$orig.ident='B11'
dataB11 = NormalizeData(dataB11)
dataB11 = ScaleData(dataB11)

dataB12 = Load10X_Spatial('./spaceranger_AZ_B12/outs', slice='B12')
dataB12$orig.ident='B12'
#Don't normalize until after we separate out the tissues

#B12 and B13 are cleanly divided in the capture area and can be separately just using the column coordinates (at position 111)

coords=dataB12@images$B12@coordinates
plot(coords[,c('row','col')])
abline(111,0)   #Spots with col > 111 go to B13 and any with col < 111 go to B12

B12cells = rownames(coords)[coords$col < 111]
B13cells = rownames(coords)[coords$col > 111]
plot(coords[B12cells,c('row','col')])
plot(coords[B13cells,c('row','col')])   #Everything looks good. Add to metadata and subset to get each

dataB12@meta.data$tissue = rep(NA, nrow(dataB12@meta.data))
dataB12$tissue[B12cells] = 'B12'
dataB12$tissue[B13cells] = 'B13'
table(dataB12$tissue)   #sum is 7053, same as ncol(dataB12)

#Subdivide into separate pieces. Do B13 first, since we'll re-use dataB12 for the final B12
dataB13 = subset(dataB12, tissue=='B13')
dataB13$orig.ident = 'B13'
dataB12 = subset(dataB12, tissue=='B12')

dataB12 = NormalizeData(dataB12)
dataB12 = ScaleData(dataB12)
dataB13 = NormalizeData(dataB13)
dataB13 = ScaleData(dataB13)

dataB14 = Load10X_Spatial('./spaceranger_AZ_B14/outs', slice='B14')
dataB14$orig.ident='B14'
dataB14 = NormalizeData(dataB14)
dataB14 = ScaleData(dataB14)

dataB15 = Load10X_Spatial('./spaceranger_AZ_B15/outs', slice='B15')
dataB15$orig.ident='B15'
dataB15 = NormalizeData(dataB15)
dataB15 = ScaleData(dataB15)

dataB16 = Load10X_Spatial('./spaceranger_AZ_B16/outs', slice='B16')
dataB16$orig.ident='B16'
dataB16 = NormalizeData(dataB16)
dataB16 = ScaleData(dataB16)

dataB17 = Load10X_Spatial('./spaceranger_AZ_B17/outs', slice='B17')
dataB17$orig.ident='B17'
dataB17 = NormalizeData(dataB17)
dataB17 = ScaleData(dataB17)

dataB18 = Load10X_Spatial('./spaceranger_AZ_B18/outs', slice='B18')
dataB18$orig.ident='B18'
dataB18 = NormalizeData(dataB18)
dataB18 = ScaleData(dataB18)

#Create merged Seurat Object
data.merge.revision = merge(dataB1, list(dataB5, dataB9, dataB10, dataB11, dataB12, dataB13,
                                         dataB14, dataB15, dataB16, dataB17, dataB18))

#Clean up the workspace to save memory
rm(dataB1, dataB5, dataB9, dataB10, dataB11, dataB12, dataB13, dataB14, dataB15, dataB16, dataB17, dataB18)
gc()

#Define putative DC1+ spots as expressing ITGAX and at least one of BATF3, XCR1, or CLEC9A. Define CD8+ as expressing CD8A (CytAssist does not sequence CD8B)
data.merge.revision$itgax = sign(data.merge.revision$Spatial@counts['ITGAX',])
data.merge.revision$dc1=sign(apply(data.merge.revision$Spatial@counts[c('BATF3','XCR1','CLEC9A'),],2,sum))
data.merge.revision$dc1 = data.merge.revision$itgax*data.merge.revision$dc1
data.merge.revision$cd8 = sign(data.merge.revision$Spatial@counts['CD8A',]) 

#Create Giotto object
slicelist = names(data.merge.revision@images)   #Note: B12 = B12; B12.1 = B13

locations_revision=data.merge.revision@images[[slicelist[1]]]@coordinates[,4:5]
locations_revision$sample=slicelist[1]
for(i in 2:length(slicelist)){
  temp_loc=data.merge.revision@images[[slicelist[i]]]@coordinates[,4:5]
  temp_loc$imagecol=temp_loc$imagecol+max(locations_revision$imagecol)   #need to make sure things stay separate and avoid overlaps with coordinates
  temp_loc$imagerow=temp_loc$imagerow+max(locations_revision$imagerow)
  temp_loc$sample=slicelist[i]
  locations_revision=rbind(locations_revision, temp_loc)
}

counts_revision=data.merge.revision@assays$Spatial@counts
giotto_revision=createGiottoObject(raw_exprs = counts_revision, 
                                   spatial_locs = locations_revision[,1:2], 
                                   norm_expr = data.merge.revision@assays$Spatial@data,
                                   instructions = instrs)

#Used Seurat object normalization to save time, does not affect results
giotto_revision=createSpatialNetwork(gobject=giotto_revision, maximum_distance_delaunay = 500)   
giotto_revision=createSpatialKNNnetwork(gobject=giotto_revision, k=6)   #We will use the KNN version downstream

#Steps for finding exclusions are facilitated with a wrapper function called dc1prep(), found at the end of the file
##First find neighbors to original DC1s and extend inclusion of CD8+ DC1+ to include CD8+ neighboring a DC1+ spot (used in all analyses 
#referring to "Extended" spot definition)

#Visualize the spot exclusions

colscheme=c('DC1+'='#3366FF',
            'DC1-CD8+'='#FFFF00',
            'DC1+CD8+'='#008000',
            'DC1+CD8nb'='#7ECB01',
            'DC1-CD8-'='#B2B2B2',
            'Excluded'='#ECECEC',
            'Excluded CD8+'='#FFCC00')


#Find spot exclusions, using KNN for neighbor distance definition and excluding CD8 within 2 spots of a DC1+ spot.
dc1cells = names(which(data.merge.revision$dc1==1))
data.merge.revision = dc1prep(scobject = data.merge.revision, 
                              giotto_object = giotto_revision,
                              cells = dc1cells, network = 'knn',
                              exclude=2)

# Example of spot exclusion layout using sample B11

#B11
p1<- SpatialDimPlot(data.merge.revision, group.by='Legend', pt.size.factor=1.5,images='B11')+
  scale_fill_manual(values=colscheme)+
  ggtitle('B11')+
  theme(title=element_text(size=18, face='bold'),legend.title=element_text(size=16),
        legend.text = element_text(size=16),legend.key.size = unit(8,'mm'))

pdf('./Results/Spot_type_distances_B11_wExclusions.pdf',height=6,width=15)
print(p1)
dev.off()

#Prepare the gene lists for heatmaps

gene.sets=list(
  stemlike=c('TCF7','SLAMF6','IL7R','BCL6','SELL','TNFRSF25','LEF1','GPR183','LTB','EVL'),
  dysfunctional=c('HAVCR2','TIGIT','ENTPD1','LAG3','CTLA4','PDCD1','TOX','LAYN','PRDM1','BATF','SIRPG','LYST',
                  'IGFLR1','PAG1','IKZF3'),
  effector=c('TNFRSF4','TNFRSF9','CD44','NKG7','GNLY','CST7','LITAF','GZMB','PRF1','ITGB1','LTA',
             'IFNG','BHLHE40','ZEB2','CXCL13'),
  chemokine=c('CCL4','CCL5','FLT3LG','CXCL9','CXCL10','CXCL11','TNFSF9'))

newgenelist=unique(unlist(gene.sets))
setdiff(newgenelist, rownames(data.merge.revision))

#Build objects for each sample with scaling before or after exclusions and normalization 

#simple exclusions, scaled after, excludes a possible CD8 normalization step
nonormlist = list()
samplelist = unique(data.merge.revision$orig.ident)

for(i in 1:length(samplelist)){
  nonormlist[[samplelist[i]]] = subset(data.merge.revision, orig.ident == samplelist[i])
  nonormlist[[samplelist[i]]] = subset(nonormlist[[samplelist[i]]], cd8==1)
  
  #No CD8 NORMALIZATION STEP
  #nonormlist[[samplelist[i]]]$Spatial@counts = t(t(nonormlist[[samplelist[i]]]$Spatial@counts)/nonormlist[[samplelist[i]]]$Spatial@counts['CD8A',])
  
  nonormlist[[samplelist[i]]] = subset(nonormlist[[samplelist[i]]], dc1ex_vis!='Excluded')
  nonormlist[[samplelist[i]]] = NormalizeData(nonormlist[[samplelist[i]]])    
  nonormlist[[samplelist[i]]] = ScaleData(nonormlist[[samplelist[i]]])
  
  nonormlist[[samplelist[i]]]$dc1ex_vis[nonormlist[[samplelist[i]]]$dc1ex_vis=='Negative']='CD8+ DC1-'
  nonormlist[[samplelist[i]]]$dc1ex_vis[nonormlist[[samplelist[i]]]$dc1ex_vis=='Positive']='CD8+ DC1+'
  nonormlist[[samplelist[i]]]$dc1ex_vis[nonormlist[[samplelist[i]]]$dc1ex_vis=='CD8+ DC1-']='DC1-'
  nonormlist[[samplelist[i]]]$dc1ex_vis[nonormlist[[samplelist[i]]]$dc1ex_vis=='CD8+ DC1+']='DC1+'
  nonormlist[[samplelist[i]]]$dc1simp_vis=nonormlist[[samplelist[i]]]$dc1
  nonormlist[[samplelist[i]]]$dc1simp_vis[nonormlist[[samplelist[i]]]$dc1simp_vis==0]='DC1-'
  nonormlist[[samplelist[i]]]$dc1simp_vis[nonormlist[[samplelist[i]]]$dc1simp_vis==1]='DC1+'
} 

temp_simp=array(0,dim=c(length(newgenelist),18))
rownames(temp_simp)=newgenelist
colnames(temp_simp)=c('DC1- B9','DC1+ B9','DC1- B10','DC1+ B10',
                      'DC1- B11','DC1+ B11','DC1- B13','DC1+ B13',
                      'DC1- B14','DC1+ B14','DC1- B15','DC1+ B15',
                      'DC1- B16','DC1+ B16','DC1- B17','DC1+ B17',
                      'DC1- B18','DC1+ B18')

count = 1
for(i in 1:length(samplesubset)){
  temp_simp[,count]=apply(subset(nonormlist[[samplesubset[i]]],dc1simp_vis=='DC1-')$Spatial@scale.data[newgenelist,],1,mean)  
  count = count + 1
  temp_simp[,count]=apply(subset(nonormlist[[samplesubset[i]]],dc1simp_vis=='DC1+')$Spatial@scale.data[newgenelist,],1,mean)
  count = count + 1
}

temp_reorder = temp_simp[,c(seq(1,18,2),seq(2,18,2))]
temp_reorder = temp_reorder[,-c(8,17)]    #Remove B17 due to sample quality issues

coltitles = c('B9<br>DC1-','B10<br>DC1-','B11<br>DC1-','B13<br>DC1-',
              'B14<br>DC1-','B15<br>DC1-','B16<br>DC1-','B18<br>DC1-',
              'B9<br>DC1+','B10<br>DC1+','B11<br>DC1+','B13<br>DC1+',
              'B14<br>DC1+','B15<br>DC1+','B16<br>DC1+','B18<br>DC1+')
colgap = c(rep(0,7),2,rep(0,7))


colrange=brewer.pal(11,'RdBu')
#col_fun=circlize::colorRamp2(c(-2,0,2),c(colrange[9],'#FFFFFF',colrange[2]))
col_fun2=circlize::colorRamp2(c(-1,0,1),c(colrange[9],'#FFFFFF',colrange[2]))

hmean=meanheatmap(object=temp_reorder, geneset='chemokine', colfun=col_fun2, coltitles, colgap)
figsize=calc_ht_size(hmean)
pdf(paste0('./Results/heatmap_simple_mean_scaledafter_chemokine_dist2.pdf'),height=figsize[2],width=figsize[1])
print(hmean)
dev.off()

#Scaled before exclusions, CD8 normalized
prescale_wCD8 = list()

for(i in 1:length(samplelist)){
  prescale_wCD8[[samplelist[i]]] = subset(data.merge.revision, orig.ident==samplelist[i])
  #Subset to only CD8+ and Normalize by CD8 expression
  prescale_wCD8[[samplelist[i]]] = subset(prescale_wCD8[[samplelist[i]]], cd8==1)
  prescale_wCD8[[samplelist[i]]]$Spatial@counts = t(t(prescale_wCD8[[samplelist[i]]]$Spatial@counts)/prescale_wCD8[[samplelist[i]]]$Spatial@counts['CD8A',])
  #log-normalize and scale
  prescale_wCD8[[samplelist[i]]] = NormalizeData(prescale_wCD8[[samplelist[i]]])
  prescale_wCD8[[samplelist[i]]] = ScaleData(prescale_wCD8[[samplelist[i]]])
  #Exclude based on the distane criteria
  prescale_wCD8[[samplelist[i]]] = subset(prescale_wCD8[[samplelist[i]]], dc1ex_vis!='Excluded')
  
  prescale_wCD8[[samplelist[i]]]$dc1ex_vis[prescale_wCD8[[samplelist[i]]]$dc1ex_vis=='Negative'] = 'CD8+ DC1-'
  prescale_wCD8[[samplelist[i]]]$dc1ex_vis[prescale_wCD8[[samplelist[i]]]$dc1ex_vis=='Positive'] = 'CD8+ DC1+'
  prescale_wCD8[[samplelist[i]]]$dc1ex_vis[prescale_wCD8[[samplelist[i]]]$dc1ex_vis=='CD8+ DC1-'] = 'DC1-'
  prescale_wCD8[[samplelist[i]]]$dc1ex_vis[prescale_wCD8[[samplelist[i]]]$dc1ex_vis=='CD8+ DC1+'] = 'DC1+'
  prescale_wCD8[[samplelist[i]]]$dc1simp_vis=prescale_wCD8[[samplelist[i]]]$dc1
  prescale_wCD8[[samplelist[i]]]$dc1simp_vis[prescale_wCD8[[samplelist[i]]]$dc1simp_vis==0] = 'DC1-'
  prescale_wCD8[[samplelist[i]]]$dc1simp_vis[prescale_wCD8[[samplelist[i]]]$dc1simp_vis==1] = 'DC1+'
}

#Means for the extended data
temp=array(0,dim=c(length(newgenelist),24))
rownames(temp)=newgenelist
colnames(temp)=c('DC1- B1','DC1+ B1','DC1- B5','DC1+ B5',
                  'DC1- B9','DC1+ B9','DC1- B10','DC1+ B10',
                  'DC1- B11','DC1+ B11','DC1- B12','DC1+ B12',
                  'DC1- B13','DC1+ B13','DC1- B14','DC1+ B14',
                 'DC1- B15','DC1+ B15','DC1- B16','DC1+ B16',
                 'DC1- B17','DC1+ B17','DC1- B18','DC1+ B18')

count=1
for(i in 1:length(samplelist)){
  temp[,count] = apply(subset(prescale_wCD8[[samplelist[i]]], dc1ex_vis=='DC1-')$Spatial@scale.data[newgenelist,],1,mean)
  count = count + 1
  temp[,count] = apply(subset(prescale_wCD8[[samplelist[i]]], dc1ex_vis=='DC1+')$Spatial@scale.data[newgenelist,],1,mean)
  count = count + 1
}

quantile(temp)

#Reorder with DC1- all on the left, B1, B17 not included due to sample quality
temp_reorder = temp[,-c(1,2,21,22)]
temp_reorder = temp_reorder[,c(seq(1,20,2),seq(2,20,2))]
coltitles = c('B5<br>DC1-','B9<br>DC1-','B10<br>DC1-','B11<br>DC1-','B12<br>DC1-','B13<br>DC1-',
              'B14<br>DC1-','B15<br>DC1-','B16<br>DC1-','B18<br>DC1-',
              'B5<br>DC1+','B9<br>DC1+','B10<br>DC1+','B11<br>DC1+','B12<br>DC1+','B13<br>DC1+',
              'B14<br>DC1+','B15<br>DC1+','B16<br>DC1+','B18<br>DC1+')
colgap = c(rep(0,9),2,rep(0,9))

#Dysfunctional
hmean=meanheatmap(object=temp_reorder, geneset='dysfunctional', colfun=col_fun2, coltitles, colgap)
figsize=calc_ht_size(hmean)
pdf(paste0('./Results/heatmap_extended_mean_scaledbefore_CD8Normalized_dysfunctional_dist2_B5B9toB16B18_DC1sorted_noMYO7B.pdf'),height=figsize[2],width=figsize[1])
print(hmean)
dev.off()

#Stemlike
hmean=meanheatmap(object=temp_reorder, geneset='stemlike', colfun=col_fun2, coltitles, colgap)
figsize=calc_ht_size(hmean)
pdf(paste0('./Results/heatmap_extended_mean_scaledbefore_CD8Normalized_stemlike_dist2_B5B9toB16B18_DC1sorted.pdf'),height=figsize[2],width=figsize[1])
print(hmean)
dev.off()

#Effector
hmean=meanheatmap(object=temp_reorder, geneset='effector', colfun=col_fun2, coltitles, colgap)
figsize=calc_ht_size(hmean)
pdf(paste0('./Results/heatmap_extended_mean_scaledbefore_CD8Normalized_effector_dist2_B5B9toB16B18_DC1sorted.pdf'),height=figsize[2],width=figsize[1])
print(hmean)
dev.off()

#GSEA

#First, build pseudobulk objects

#Subset to only CD8+ spots and to non-excluded spots
data.merge.cd8=subset(data.merge.revision, cd8==1)
data.merge.cd8=subset(data.merge.cd8, dc1ex_vis!='Excluded')  #Exclusion based on exclude=2 from dc1prep function

samplelist = c('B5','B9','B10','B11','B12','B13','B14','B15','B16','B18')  #Excludes B1 and B17 due to sample quality
data.merge.cd8=subset(data.merge.revision, cd8==1)
data.merge.cd8=subset(data.merge.cd8, dc1ex_vis!='Excluded')  #Exclusion based on exclude=2 from dc1prep function

#first make the extended, normalized version
data.merge.cd8$Spatial@data=t(t(data.merge.cd8$Spatial@counts)/data.merge.cd8$Spatial@counts['CD8A',])   #normalize by Actb instead of by total counts. Manually do remaining steps of NormalizeData

pseudolist_extnorm = list()
for(i in 1:length(samplelist)){
  pseudolist_extnorm[[samplelist[i]]] = pseudoprep(data.merge.cd8, samplelist[i], type='extended', metalevel = 'orig.ident', slot='data')
}

pseudoset_extnorm = cbind(pseudolist_extnorm[[samplelist[1]]], pseudolist_extnorm[[samplelist[2]]])
for(i in 3:length(samplelist)){
  pseudoset_extnorm = cbind(pseudoset_extnorm, pseudolist_extnorm[[samplelist[i]]])
}

colnames(pseudoset_extnorm)=c('B5pos','B5neg','B9pos','B9neg','B10pos','B10neg','B11pos','B11neg','B12pos','B12neg','B13pos','B13neg',
                              'B14pos','B14neg','B15pos','B15neg','B16pos','B16neg','B18pos','B18neg')

sampleinfo_extnorm = cbind(str_extract(colnames(pseudoset_extnorm),'B[0-9]+'), str_replace_all(colnames(pseudoset_extnorm),'B[0-9]+',''))
sampleinfo_extnorm = cbind(sampleinfo_extnorm, paste0(sampleinfo_extnorm[,1],sampleinfo_extnorm[,2]))
sampleinfo_extnorm = as.data.frame(sampleinfo_extnorm)
colnames(sampleinfo_extnorm) = c('patient','dc1','sample')
sampleinfo_extnorm$patient = factor(sampleinfo_extnorm$patient, levels=samplelist)
sampleinfo_extnorm$sample = factor(sampleinfo_extnorm$sample, levels=colnames(pseudoset_extnorm))
sampleinfo_extnorm$dc1=factor(sampleinfo_extnorm$dc1, levels=c('pos','neg'))

myDGEList_extnorm <- DGEList(counts = pseudoset_extnorm, 
                             samples = sampleinfo_extnorm,
                             group = sampleinfo_extnorm$patient)

#Repeat but with simple definition (excludes B5 and B12 due to insufficient DC1+ spots)
samplelist_simp = c('B9','B10','B11','B13','B14','B15','B16','B18')

pseudolist_simp = list()
for(i in 1:length(samplelist_simp)){
  pseudolist_simp[[samplelist_simp[i]]] = pseudoprep(data.merge.cd8, samplelist_simp[i], type='simple', metalevel = 'orig.ident', slot='counts')
}

pseudoset_simp = cbind(pseudolist_simp[[samplelist_simp[1]]], pseudolist_simp[[samplelist_simp[2]]])
for(i in 3:length(samplelist_simp)){
  pseudoset_simp = cbind(pseudoset_simp, pseudolist_simp[[samplelist_simp[i]]])
}

colnames(pseudoset_simp)=c('B9pos','B9neg','B10pos','B10neg','B11pos','B11neg','B13pos','B13neg',
                           'B14pos','B14neg','B15pos','B15neg','B16pos','B16neg','B18pos','B18neg')

sampleinfo_simp = cbind(str_extract(colnames(pseudoset_simp),'B[0-9]+'), str_replace_all(colnames(pseudoset_simp),'B[0-9]+',''))
sampleinfo_simp = cbind(sampleinfo_simp, paste0(sampleinfo_simp[,1],sampleinfo_simp[,2]))
sampleinfo_simp = as.data.frame(sampleinfo_simp)
colnames(sampleinfo_simp) = c('patient','dc1','sample')
sampleinfo_simp$patient = factor(sampleinfo_simp$patient, levels=samplelist_simp)
sampleinfo_simp$sample = factor(sampleinfo_simp$sample, levels=colnames(pseudoset_simp))
sampleinfo_simp$dc1=factor(sampleinfo_simp$dc1, levels=c('pos','neg'))

myDGEList_simp <- DGEList(counts = pseudoset_simp, 
                          samples = sampleinfo_simp,
                          group = sampleinfo_simp$patient)

#Next due basic preprocessing and use limma-voom to build logFC results for GSEA
#extnorm results
log2.cpm <- cpm(myDGEList_extnorm,log=TRUE)
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")

keep.exprs <- filterByExpr(myDGEList_extnorm,
                           min.count=5,
                           min.total.count=5,
                           min.prop=0.5)       #Use stringent criteria for keeping genes due to uneven data sampling in Visium
myDGEList.filtered <- myDGEList_extnorm[keep.exprs,]

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "gene_name")

myDGEList.filtered.TMM <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.TMM <- cpm(myDGEList.filtered.TMM, log=TRUE)
log2.cpm.filtered.TMM.df <- as_tibble(log2.cpm.filtered.TMM, rownames = "gene_name")

patient <- factor(sampleinfo_extnorm$patient, levels=samplelist)
dc1 <- factor(sampleinfo_extnorm$dc1)
design <- model.matrix(~0 + patient + dc1)  
colnames(design) = c(samplelist, 'DC1')

v.DEGList <- voomWithQualityWeights(myDGEList.filtered.TMM, design = design, normalize.method = 'none', plot = FALSE)
fit <- lmFit(v.DEGList, design)

contrast.matrix <- makeContrasts(DC1pos.vs.DC1neg = -DC1,
                                 levels=design)

fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)

myTopHits_extnorm=topTable(ebFit, adjust ="BH", coef='DC1pos.vs.DC1neg', number=Inf, sort.by="P")
myTopHits_extnorm=myTopHits_extnorm%>%as_tibble(rownames = 'gene_name')


#Load GSEA reference data
hs_gsea_h <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, gene_symbol)
hs_gsea_go_bp <- msigdbr(species = "Homo sapiens", category="C5", subcategory="GO:BP")%>% dplyr::select(gs_name, gene_symbol)
hs_gsea_go_cc <- msigdbr(species = "Homo sapiens", category="C5", subcategory="GO:CC")%>% dplyr::select(gs_name, gene_symbol)
hs_gsea_go_mf <- msigdbr(species = "Homo sapiens", category="C5", subcategory="GO:MF")%>% dplyr::select(gs_name, gene_symbol)

gsea.set_extnorm=gseaPrep(myTopHits_extnorm)   #estimate mytophits first with B1/B5 ex

#Check number of results for each gene collection
lapply(gsea.set_extnorm, dim)    #no results for GO_MF

#Combine results into one object
allgsea.set_extnorm=rbind(gsea.set_extnorm$H@result, gsea.set_extnorm$GO_BP@result)
allgsea.set_extnorm=rbind(allgsea.set_extnorm, gsea.set_extnorm$GO_CC@result)
allgsea.set_extnorm=rbind(allgsea.set_extnorm, gsea.set_extnorm$GO_MF@result)

write.csv(allgsea.set_extnorm,'./Results/GSEA_extendedDC1_B5B9toB16B18_CD8normalized_notlogscaled_wGOMF.csv',row.names=F,quote=F)

#simple, nonnormalized results
log2.cpm <- cpm(myDGEList_simp,log=TRUE)
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")

keep.exprs <- filterByExpr(myDGEList_simp,
                           min.count=5,
                           min.total.count=5,
                           min.prop=0.5)    #Previously VERY stringent, but had only 2 samples. Relax this?  Still pick up a lot of low expression
myDGEList.filtered <- myDGEList_simp[keep.exprs,]

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "gene_name")

myDGEList.filtered.TMM <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.TMM <- cpm(myDGEList.filtered.TMM, log=TRUE)
log2.cpm.filtered.TMM.df <- as_tibble(log2.cpm.filtered.TMM, rownames = "gene_name")

patient <- factor(sampleinfo_simp$patient, levels=samplelist_simp)
dc1 <- factor(sampleinfo_simp$dc1)
design <- model.matrix(~0 + patient + dc1)  
colnames(design) = c(samplelist_simp, 'DC1')

v.DEGList <- voomWithQualityWeights(myDGEList.filtered.TMM, design = design, normalize.method = 'none', plot = FALSE)
fit <- lmFit(v.DEGList, design)

contrast.matrix <- makeContrasts(DC1pos.vs.DC1neg = -DC1,
                                 levels=design)

fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)

myTopHits_simp=topTable(ebFit, adjust ="BH", coef='DC1pos.vs.DC1neg', number=Inf, sort.by="P")
myTopHits_simp=myTopHits_simp%>%as_tibble(rownames = 'gene_name')

gsea.set_simp=gseaPrep(myTopHits_simp)   #estimate mytophits first with B1/B5 ex

#Check number of results for each gene collection
lapply(gsea.set_simp, dim)    #no results for GO_MF

#Combine results into one object
allgsea.set_simp=rbind(gsea.set_simp$H@result, gsea.set_simp$GO_BP@result)
allgsea.set_simp=rbind(allgsea.set_simp, gsea.set_simp$GO_CC@result)
allgsea.set_simp=rbind(allgsea.set_simp, gsea.set_simp$GO_MF@result)

write.csv(allgsea.set_simp,'./Results/GSEA_simpleDC1_B5B9toB16B18_notnormalized_wGOMF.csv',row.names=F,quote=F)

pathlist = c('HALLMARK_PI3K_AKT_MTOR_SIGNALING',
             'GOBP_REGULATION_OF_ANTIGEN_RECEPTOR_MEDIATED_SIGNALING_PATHWAY',
             'GOBP_POSITIVE_REGULATION_OF_T_CELL_PROLIFERATION',
             'GOBP_ALPHA_BETA_T_CELL_PROLIFERATION',
             'GOBP_NEGATIVE_REGULATION_OF_LEUKOCYTE_APOPTOTIC_PROCESS',
             'GOBP_T_CELL_MIGRATION',
             'GOBP_POSITIVE_REGULATION_OF_ALPHA_BETA_T_CELL_ACTIVATION',
             'GOBP_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE',
             'GOCC_IMMUNOLOGICAL_SYNAPSE',
             'GOBP_RESPONSE_TO_INTERLEUKIN_7',
             'GOBP_INTERLEUKIN_2_PRODUCTION',
             'GOBP_RESPONSE_TO_INTERFERON_GAMMA',
             'GOBP_INTERFERON_GAMMA_PRODUCTION',
             'GOBP_POSITIVE_REGULATION_OF_TUMOR_NECROSIS_FACTOR_SUPERFAMILY_CYTOKINE_PRODUCTION',
             'GOBP_LYMPHOCYTE_COSTIMULATION',
             'GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY',
             'GOBP_LEUKOCYTE_DEGRANULATION',
             'GOBP_CHRONIC_INFLAMMATORY_RESPONSE',
             'GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY',
             'GOBP_POSITIVE_REGULATION_OF_ADAPTIVE_IMMUNE_RESPONSE',
             'GOMF_NF_KAPPAB_BINDING',
             'GOMF_CHEMOKINE_RECEPTOR_BINDING',
             'GOMF_PHOSPHATIDYLINOSITOL_3_KINASE_BINDING')

gsea.set_extnorm_subset = subset(allgsea.set_extnorm, ID%in%pathlist)
pdf('./Results/gsea_extendedDC1_bubbleplot_PathwaySubset_B5B9toB16B18__CD8Normalized.pdf', height=15, width=18)
gsea_bubble(gsea.set_extnorm_subset)+
  ggtitle('GSEA, extended DC1+')
dev.off()


#Pathways for plots with gsea.set_simp
pathlist2 = c(
  'HALLMARK_ALLOGRAFT_REJECTION',
  'HALLMARK_INFLAMMATORY_RESPONSE',
  'HALLMARK_INTERFERON_ALPHA_RESPONSE',
  'GOBP_MYELOID_LEUKOCYTE_ACTIVATION',
  'GOBP_DENDRITIC_CELL_DIFFERENTIATION',
  'GOBP_DENDRITIC_CELL_MIGRATION',
  'GOBP_MYELOID_LEUKOCYTE_MEDIATED_IMMUNITY',
  'GOBP_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY',
  'GOBP_MYELOID_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE',
  'HALLMARK_IL6_JAK_STAT3_SIGNALING',
  'GOMF_CHEMOKINE_BINDING'
)

gsea.set_simp_subset = subset(allgsea.set_simp, ID%in%pathlist2)
pdf('./Results/gsea_simpleDC1_bubbleplot_PathwaySubset_B5B9toB16B18_notnormalized.pdf', height=7.5, width=18)
gsea_bubble(gsea.set_simp_subset)+
  ggtitle('GSEA, simple DC1+')
dev.off()



## Functions

dc1prep<-function(scobject, giotto_object, cells, exclude=2, network = 'knn'){
  if(network=='knn'){
    gnet = 'knn_network'
  }
  if(network=='delaunay'){
    gnet = 'Delaunay_network'
  }
  dc1nb=findNetworkNeighbors(giotto_object, spatial_network_name = gnet, source_cell_ids = cells)
  dc1neighbors=dc1nb$cell_ID[dc1nb$nb_cells!='others']  #keep source + neighbors
  dc1nb_nb=findNetworkNeighbors(giotto_object, spatial_network_name = gnet, source_cell_ids = dc1neighbors)
  dc1nb_neighbors=dc1nb_nb$cell_ID[dc1nb_nb$nb_cells=='neighbor'] #keep new neighbors
  
  scobject$dc1nb=rep(0,nrow(scobject@meta.data))
  scobject$dc1nb[dc1neighbors]=1    #IMPORTANT: the initial NB group is both the original POS spots + Neighbors
  scobject$dc1nb_nb=rep(0,nrow(scobject@meta.data))
  scobject$dc1nb_nb[dc1nb_neighbors]=1  #second group adds the next neighbor "shell"
  
  dc1cd8=which((scobject$dc1==1&scobject$cd8==1)|(scobject$dc1nb==1&scobject$cd8==1))   #define CD8+ DC1+ as any spot that is itself CD8+ DC1+ or is CD8+ and neighbor to DC1+
  scobject$dc1cd8=rep(0,nrow(scobject@meta.data))
  scobject$dc1cd8[dc1cd8]=1   #Expands the possible spots from 11 to 74!
  
  #Next, need to define exclusion zone as any non-positive spot (by dc1cd8) within 3 spaces of a positive spot.
  dc1cd8spots=names(which(scobject$dc1cd8==1))
  dc1cd8nb=findNetworkNeighbors(giotto_object, spatial_network_name = gnet, source_cell_ids = dc1cd8spots)
  dc1cd8nb_1=dc1cd8nb$cell_ID[dc1cd8nb$nb_cells=='neighbor']    #direct neighbor to a CD8+ DC1+
  dc1cd8nb_nb=findNetworkNeighbors(giotto_object, spatial_network_name = gnet, source_cell_ids = dc1cd8nb_1)
  dc1cd8nb_2=dc1cd8nb_nb$cell_ID[dc1cd8nb_nb$nb_cells=='neighbor'&dc1cd8nb$nb_cells!='source'&dc1cd8nb$nb_cells!='both']    #Find neighbors to neighbors but make sure to exclude the original sources
  dc1cd8nb_nb_nb=findNetworkNeighbors(giotto_object, spatial_network_name = gnet, source_cell_ids = dc1cd8nb_2)
  dc1cd8nb_3=dc1cd8nb_nb_nb$cell_ID[dc1cd8nb_nb_nb$nb_cells=='neighbor'&dc1cd8nb$nb_cells!='source'&dc1cd8nb_nb$nb_cells!='source'&dc1cd8nb$nb_cells!='both'&dc1cd8nb_nb$nb_cells!='both']    #Find neighbors to neighbors but make sure to exclude the original sources
  
  #Note: want distance from the original DC1 not from the inferred CD8+ DC1+ spots. So we define this separately using dc1nb_nb and dc1nb
  dc1nb_2 = names(scobject$dc1nb_nb)[scobject$dc1nb_nb==1]
  dc1nb_nb_nb = findNetworkNeighbors(giotto_object, spatial_network_name = gnet, source_cell_ids = dc1nb_2)
  
  #Modified to allow different distance options quickly. Default is exclude=3
  if(exclude==3){
    dc1nb_ex = unique(dc1nb_nb_nb$cell_ID[dc1nb_nb_nb$nb_cells=='neighbor'|dc1nb_nb$nb_cells=='neighbor'|dc1nb_nb$nb_cells=='both'|dc1nb$nb_cells=='neighbor'])
  }
  if(exclude==2){
    dc1nb_ex = unique(dc1nb_nb$cell_ID[dc1nb_nb$nb_cells=='neighbor'|dc1nb_nb$nb_cells=='both'|dc1nb$nb_cells=='neighbor'])
  }
  #This will include the original dc1s and some of the dc1cd8s that we want to exclude due to distances.
  keep=unique(c(cells, dc1cd8spots))
  dc1nb_ex=setdiff(dc1nb_ex, keep)
  if(exclude==3){
    scobject$dc1nb_ex=scobject$dc1nb_nb  
  }
  if(exclude==2){
    scobject$dc1nb_ex=rep(0,nrow(scobject@meta.data))
  }
  
  scobject$dc1nb_ex[dc1nb_ex]=1
  scobject$dc1ex_vis=rep('Negative',nrow(scobject@meta.data))
  scobject$dc1ex_vis[dc1cd8spots]='Positive'
  scobject$dc1ex_vis[scobject$dc1nb_ex==1]='Excluded'
  #scobject$dc1ex_vis[scobject$dc1==1&scobject$cd8==0]='Excluded'
  
  scobject$dc1cd8_vis = rep('Negative',nrow(scobject@meta.data))
  scobject$dc1cd8_vis[dc1cd8spots]='Positive'
  scobject$dc1cd8_vis[dc1cd8nb_1]='Excluded'
  scobject$dc1cd8_vis[dc1cd8nb_2]='Excluded'
  scobject$dc1cd8_vis[dc1cd8nb_3]='Excluded'
  
  scobject$dc1_vis = rep('dc1-',nrow(scobject@meta.data))
  if(exclude==3){
    scobject$dc1_vis[scobject$dc1nb_nb==1]='Excluded'  
  }
  if(exclude==2){
    scobject$dc1_vis[scobject$dc1nb==1]='Excluded'  
  }
  
  scobject$dc1_vis[scobject$dc1nb==1]='DC1 Neighbor'
  scobject$dc1_vis[scobject$dc1==1]='DC1+'
  
  #Visualize the spot exclusions
  newvis = rep(0,ncol(scobject))
  newvis[scobject$dc1==1]='DC1+'
  newvis[scobject$cd8==1&(scobject$dc1==0|scobject$dc1nb==0)]='DC1-CD8+'
  newvis[scobject$cd8==1&scobject$dc1==1]='DC1+CD8+'
  newvis[scobject$cd8==1&scobject$dc1nb==1&scobject$dc1==0]='DC1+CD8nb'
  newvis[scobject$cd8==0&scobject$dc1==0]='DC1-CD8-'
  #newvis[scobject$dc1ex_vis=='Excluded'&scobject$dc1==0]='Excluded'
  newvis[scobject$dc1ex_vis=='Excluded'&scobject$dc1==0]='Excluded'
  newvis[scobject$dc1ex_vis=='Excluded'&scobject$cd8==1]='Excluded CD8+'
  newvis=factor(newvis, levels=c('DC1+','DC1-CD8+','DC1+CD8+','DC1+CD8nb','DC1-CD8-','Excluded','Excluded CD8+'))
  scobject@meta.data$Legend=newvis
  
  return(scobject)
}

meanheatmap<-function(object, geneset, colfun=col_fun, coltitles, colgap){
  Heatmap(object[gene.sets[[geneset]],],
          column_split=factor(colnames(object),levels=colnames(object)),
          border=TRUE,
          cluster_rows=FALSE,
          cluster_columns = FALSE,
          row_names_side='left',
          name='Mean Scale Expr',
          col = colfun,
          row_names_gp=gpar(fontsize = 16, fontface='bold'),
          column_title = gt_render(coltitles),
          column_title_gp = gpar(fontsize=18, fontface='bold'),
          column_gap=unit(colgap,'mm'),
          show_column_names = FALSE,
          heatmap_legend_param = list(title_gp=gpar(fontsize=12,fontface='bold'),labels_gp=gpar(fontsize=12)),
          width=unit(460,'mm'),
          height = nrow(object[gene.sets[[geneset]],])*unit(7, "mm")
  )
}

#calc_ht_size() is a utility function suggested by the authors of ComplexHeatmap to help control cell sizes in pdf outputs across heatmaps
calc_ht_size<-function(ht, unit='inch'){
  pdf(NULL)
  ht = draw(ht)
  w = ComplexHeatmap:::width(ht)
  w = convertX(w, unit, valueOnly = TRUE)
  h = ComplexHeatmap:::height(ht)
  h = convertY(h, unit, valueOnly = TRUE)
  dev.off()
  
  c(w, h)
}

#Pseudobulk functions
pseudoprep <- function(object, sample, type='extended',metalevel='section', slot){
  if(type=='extended'){
    pseudotemp=pseudowrap(subset(data.merge.cd8, subset=(dc1ex_vis=='Positive'&orig.ident==sample)), metalevel, slot)
    pseudotemp_no=pseudowrap(subset(data.merge.cd8, subset=(dc1ex_vis=='Negative'&orig.ident==sample)), metalevel, slot)
  }
  if(type=='simple'){
    pseudotemp=pseudowrap(subset(data.merge.cd8, subset=(dc1==1&orig.ident==sample)), metalevel, slot)
    pseudotemp_no=pseudowrap(subset(data.merge.cd8, subset=(dc1==0&orig.ident==sample)), metalevel, slot)
  }
  colnames(pseudotemp)=paste0(sample,colnames(pseudotemp),'pos')
  colnames(pseudotemp_no)=paste0(sample,colnames(pseudotemp_no),'neg')
  pseudotemp=cbind(pseudotemp, pseudotemp_no)
  
  return(pseudotemp)
}

pseudowrap<-function(scObject, metalevel, slot){
  unilevels=unique(scObject@meta.data[[metalevel]])
  nlevels=length(unilevels)
  pseudocount=array(0,dim=c(nrow(scObject),nlevels))
  rownames(pseudocount)=rownames(scObject)
  colnames(pseudocount)=unilevels
  
  for(i in 1:nlevels){
    keepcols=which(scObject@meta.data[[metalevel]]==unilevels[i])
    if(length(keepcols)>0){
      if(slot=='counts'){
        tempmat=scObject$Spatial@counts[,keepcols]
      }
      if(slot=='data'){
        tempmat=scObject$Spatial@data[,keepcols]
      }
      if(length(keepcols)>1){
        pseudocount[,i]=apply(tempmat,1,sum)
      }
      if(length(keepcols)==1){
        pseudocount[,i]=tempmat
      }
    }
    
    if(length(keepcols)==0){
      pseudocount[,i]=NA
    }
  }
  return(pseudocount)
}

#GSEA Bubbleplots
gsea_bubble<-function(result){
  myGSEA_custom.df <- as_tibble(result)
  myGSEA_custom.df <- myGSEA_custom.df %>%
    mutate(phenotype = 'DC1+ CD8+')
  myGSEA_custom.df <- myGSEA_custom.df %>%arrange(desc(NES))
  myGSEA_custom.df$ID = factor(myGSEA_custom.df$ID, levels=myGSEA_custom.df$ID)
  ggplot(myGSEA_custom.df, aes(x=phenotype, y=ID)) +     
    geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
    theme_bw()+
    scale_color_gradient(low="pink", high="darkred")+
    scale_size_continuous(range=c(5,18))+
    scale_alpha(range=c(0.3,1))+
    theme(text=element_text(face='bold', size=16),
          axis.text.x = element_text(angle=45, hjust=1, size=18),
          axis.text.y = element_text(size=16))
  
}

#GSEA wrapper function
gseaPrep <- function(markers){
  mydata.gsea <- markers[['logFC']]
  names(mydata.gsea) <- as.character(markers$gene_name)
  mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)
  mydata.gsea[mydata.gsea==Inf]=max(mydata.gsea[mydata.gsea!=Inf]+1)
  mydata.gsea[mydata.gsea==-Inf]=min(mydata.gsea[mydata.gsea!=-Inf]-1)
  
  #run GSEA
  
  myGSEA_H <- GSEA(mydata.gsea, TERM2GENE = hs_gsea_h, verbose=FALSE)
  myGSEA_BP <- GSEA(mydata.gsea, TERM2GENE = hs_gsea_go_bp, verbose=FALSE)    
  myGSEA_CC <- GSEA(mydata.gsea, TERM2GENE = hs_gsea_go_cc, verbose=FALSE)
  myGSEA_MF <- GSEA(mydata.gsea, TERM2GENE = hs_gsea_go_mf, verbose=FALSE)
  
  results=list()
  results$H=myGSEA_H
  results$GO_BP = myGSEA_BP
  results$GO_CC = myGSEA_CC
  results$GO_MF = myGSEA_MF
  
  return(results)
}