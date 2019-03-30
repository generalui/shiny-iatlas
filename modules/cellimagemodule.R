## dev module for cellular image
library(grImport2)
library(grid)
library(feather)

## Data
im_expr_df <- read_feather('data/im_expr_df.feather')
fmx_df <- read_feather('data/fmx_df.feather')
feature_df <- read_feather('data/feature_df.feather') ## may not be needed
tt_df <- fmx_df %>% select(ParticipantBarcode,Study)
im_expr_tt_df <- inner_join(im_expr_df,tt_df,by="ParticipantBarcode")

## Annotations of image objects
variable.annotations <- read_tsv('data/Image_Variables.tsv')
path.annotations <- read.table('data/PathIDs.txt',as.is=T)$V1

tumortype <- "SKCM"

## Gene averages
image.ids <- variable.annotations %>% filter(Source=="im_expr_df") %>% select(ImageVariableID) %>% as_vector() %>% as.character()
feature.ids <- variable.annotations %>% filter(Source=="im_expr_df") %>% select(FeatureLabel) %>% as_vector() %>% as.character()
gene.vals <- im_expr_tt_df %>% filter(Study==tumortype,Symbol %in% feature.ids) %>% select(Symbol,normalized_count) %>% group_by(Symbol) %>% summarize(mean_exp=mean(normalized_count))

## Cell averages
image.ids <- variable.annotations %>% filter(Source=="fmx_df") %>% select(ImageVariableID) %>% as_vector() %>% as.character()
feature.ids <- variable.annotations %>% filter(Source=="fmx_df") %>% select(FeatureLabel) %>% as_vector() %>% as.character()

df <- fmx_df %>% filter(Study==tumortype) %>% select() ## need to select feature.ids  

## Get image and convert to grid object
pic <- grImport2::readPicture("data/tcell-start-cairo-edited.svg")
w <- pictureGrob(pic)
gTree.name <- childNames(w) ## label of overall gTree object
pathlabels <- w$children[[gTree.name]]$childrenOrder ## labels and order of children 


fill.color.start <- character(length(pathlabels)) ; names(fill.color.start) <- pathlabels
for (s in pathlabels){
  fill.color.start[s] <- w$children[[gTree.name]]$children[[s]]$gp$fill 
}

fill.color.new <- character(length(pathlabels)) ; names(fill.color.new) <- pathlabels

## Labels of each of the objects. Can occur more than once.
obj.ids <- c("T_cell","ICOS",rep("PD-1",6)) ; names(obj.ids) <- pathlabels

fill.color <- c("#00FFFFFF","#FF00FFFF","#FFFF00FF") 
names(fill.color) <- c("T_cell","ICOS","PD-1")
fill.color.new <- fill.color[obj.ids] ; names(fill.color.new) <- pathlabels

for (s in pathlabels){
  w$children[[gTree.name]]$children[[s]]$gp$fill <- fill.color.new[s]
}

grid.draw(w)

## numbers will change and we need to capture them
## looks useful :
## https://stat.ethz.ch/R-manual/R-devel/library/grid/doc/grobs.pdf
