## dev module for cellular image
library(grImport2)
library(grid)
library(feather)
library(tidyverse)
library(RColorBrewer) 
source("functions/load_data.R")
source("functions/transform.R")
source("functions/utils.R")

source("functions/cellimage_utils.R")
USE_REMOTE_GS = FALSE
USE_REMOTE_BQ = FALSE
panimmune_data <- load_data()

fmx_df <- panimmune_data$fmx_df
sample_group_df <- panimmune_data$sample_group_df
im_expr_df <- panimmune_data$im_expr_df
## feature_df <- panimmune_data$feature_df ## may not be needed

## Data
tt_df <- fmx_df %>% select(ParticipantBarcode,Study)
im_expr_tt_df <- inner_join(im_expr_df,tt_df,by="ParticipantBarcode")

## load("see_inside.RData") yields

## "group_col"    "group_df"
## eg  "Subtype_Curated_Malta_Noushmehr_et_al" and fmx_df filtered to available group annotations 

genes_needed <- c("PDCD1","ICOS")
dft <- build_multi_immunomodulator_expression_df(group_df,genes_needed,group_col)

cells_needed <- c("T_cells_CD8.Aggregate2","Macrophage.Aggregate1")



## Annotations of image objects
variable.annotations <- read_tsv('data/Image_Variables.tsv')
path.annotations <- read.table('data/PathIDs.txt',as.is=T)$V1

tumortype <- "SKCM"

##
##  GENES
##  

## General gene averages and extrema
dfe <- im_expr_tt_df %>% group_by(Symbol,Study) %>% summarize(mean_exp=mean(normalized_count,na.rm=T)) ## averages for all genes and all types
min.per.gene <- dfe %>% group_by(Symbol) %>% summarize(gene_min=min(mean_exp,na.rm=T))
max.per.gene <- dfe %>% group_by(Symbol) %>% summarize(gene_max=max(mean_exp,na.rm=T))

## specific values for the choice
image.ids <- variable.annotations %>% filter(Source=="im_expr_df") %>% select(ImageVariableID) %>% as_vector() %>% as.character()
feature.ids <- variable.annotations %>% filter(Source=="im_expr_df") %>% select(FeatureLabel) %>% as_vector() %>% as.character()
gene.vals <- im_expr_tt_df %>% filter(Study==tumortype,Symbol %in% feature.ids) %>% select(Symbol,normalized_count) %>% 
              group_by(Symbol) %>% summarize(mean_exp=mean(normalized_count))
gene.vals <- dfe %>% filter(Study==tumortype,Symbol %in% feature.ids)

## Compute colors for the specific choice
min.per.gene %>% filter(Symbol %in% feature.ids)
max.per.gene %>% filter(Symbol %in% feature.ids)

gnstep <- 51
allcolors <- colorRampPalette(rev(brewer.pal(n = 7,name="Blues")))(length(breakList))
plotcolors.genes <- character(length=length(feature.ids)) ; names(plotcolors.genes) <- feature.ids

for ( f in feature.ids) {
  gmin <- min.per.gene %>% filter(Symbol==f) %>% .$gene_min
  gmax <- max.per.gene %>% filter(Symbol==f) %>% .$gene_max
  gstep <- (gmax-gmin)/(gnstep-1) ## size of step 
  breakList <- seq(gmin,gmax,gstep) 
  b <- gene.vals %>% filter(Symbol==f) %>% .$mean_exp
  cind <- min(which(!(b-breakList)>0)) ## right turnover point
  plotcolors.genes[f] <- allcolors[cind]
}

##
##  CELLS
##  


# above not sure if we need wrapr::let . See transform.R under functions 

image.ids <- variable.annotations %>% filter(Source=="fmx_df") %>% select(ImageVariableID) %>% as_vector() %>% as.character()
feature.ids <- variable.annotations %>% filter(Source=="fmx_df") %>% select(FeatureLabel) %>% as_vector() %>% as.character()

fmc <- fmx_df %>% select(Study,feature.ids) %>% .[complete.cases(.),] %>% gather(feature,fraction,feature.ids)


## Cell averages - general and extrema
fmc <- fmx_df %>% select(Study,feature.ids) %>% .[complete.cases(.),] %>% gather(feature,fraction,feature.ids)

dfc <- fmc %>% group_by(feature,Study) %>% summarize(mean_cell=mean(fraction,na.rm=T))
min.per.cell <- fmc %>%  group_by(feature) %>% summarize(min_cell=min(fraction,na.rm=T))
max.per.cell <- fmc %>%  group_by(feature) %>% summarize(max_cell=max(fraction,na.rm=T))

df <- fmx_df %>% filter(Study==tumortype) %>% select(feature.ids) ## need to select feature.ids  
cell.vals <- mean(df[[feature.ids]],na.rm = T) ## will need to change when there are more cells

cmin <- min(means.per.cell$cell_mean)
cmax <- max(means.per.cell$cell_mean)
gstep <- (gmax-gmin)/(gnstep-1) ## size of step 
breakList <- seq(gmin,gmax,gstep) 
#b <- gene.vals %>% filter(Symbol==f) %>% .$mean_exp
cind <- min(which(!(cell.vals-breakList)>0)) ## right turnover point
plotcolors[f] <- allcolors[cind]



## Get image and convert to grid object
pic <- grImport2::readPicture("data/tcell-start-cairo-edited.svg")
w <- pictureGrob(pic)
gTree.name <- childNames(w) ## label of overall gTree object
pathlabels <- w$children[[gTree.name]]$childrenOrder ## labels and order of children 

1
fill.color.start <- character(length(pathlabels)) ; names(fill.color.start) <- pathlabels
for (s in pathlabels){
  fill.color.start[s] <- w$children[[gTree.name]]$children[[s]]$gp$fill 
}

fill.color.new <- character(length(pathlabels)) ; names(fill.color.new) <- pathlabels

## Labels of each of the objects. Can occur more than once.
obj.ids <- c("T_cell","ICOS",rep("PD-1",6)) ; names(obj.ids) <- pathlabels

fill.color <- c("#00FFFF","#FF00FF","#FFFF00") 
names(fill.color) <- c("T_cell","ICOS","PD-1")
fill.color.new <- fill.color[obj.ids] ; names(fill.color.new) <- pathlabels

for (s in pathlabels){
  w$children[[gTree.name]]$children[[s]]$gp$fill <- fill.color.new[s]
}

grid.draw(w)

##
## OLDER CODE
##

## Colors - working area
gmin <- -0.25
gmax <- 0.25
gnstep <- 51
gstep <- (gmax-gmin)/(gnstep-1) ## size of step 
##gstep <- 0.01
breakList <- seq(gmin,gmax,gstep) 
allcolors <- colorRampPalette(rev(brewer.pal(n = 7,name="RdBu")))(length(breakList))

b <- 0.1457
cind <- min(which(!(b-breakList)>0)) ## right turnover point

usecolor <- allcolors[cind]

##display.brewer.pal(n=7,name="RdBu")




