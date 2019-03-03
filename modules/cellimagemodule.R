## dev module for cellular image
library(grImport2)
pic <- grImport2::readPicture("../data/tcell-start-cairo-edited.svg")
w <- pictureGrob(pic)

childNames(w)

j <- "GRID.gTree.91"

w$children$GRID.gTree.91$childrenOrder

names(w$children$GRID.gTree.91$children)
# each child is a picComplexPath

w$children$GRID.gTree.91$children$GRID.picComplexPath.70$gp$fill <- "#0000FFFF"
grid.draw(w)

## numbers will change and we need to capture them
## looks useful :
## https://stat.ethz.ch/R-manual/R-devel/library/grid/doc/grobs.pdf
