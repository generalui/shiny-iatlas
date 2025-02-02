# Shiny-iAtlas

Shiny-iAtlas is an interactive web portal that provides multiple analysis modules to visualize and explore immune response characterizations across cancer types. The app is hosted on shinyapps.io at https://isb-cgc.shinyapps.io/shiny-iatlas/ and can also be accessed via the main CRI iAtlas page at http://www.cri-iatlas.org/.

The portal is built entirely in **R** and **Shiny** using the **RStudio** development environment. Layout and interactivity within the portal are achieved by heavy use of the following packages:

+ **`shinydashboard`**
+ **`plotly`**
+ **`crosstalk`**

## Install

Install:
* R: https://www.r-project.org/
* RStudio: https://rstudio.com/products/rstudio/download/

Open RStudio and install dependencies:

```R
source("setup_installs.R")
```

## Data

Input data for the Shiny-iAtlas portal were accessed from multiple remote sources, including **Synapse**, the **ISB Cancer Genomics Cloud**, and **Google Drive**. For convenience, we have created locally cached versions of dataframe objects as **`feather`** files:

+ `fmx_df.feather`
+ `feature_df.feather`
+ `feature_method_df.feather`
+ `im_direct_relationships.feather`
+ `im_potential_factors.feather`
+ `im_expr_df.feather`
+ `sample_group_df.feather`

## Methods

While many of the results presented in tables and plots are taken directly from IRWG data (including the main **feature matrix** and various feature and group annotations), we compute some values internally. Unless otherwise noted, the following methods/tools were used to compute summary statistics:

#### Correlation — Spearman's rank-order correlation:

```R
stats::cor(x, y, method = "spearman", use = "pairwise.complete.obs")
```

#### Concordance Index (CI):

Concordance indexes for survival endpoints with respect to different immune readouts were computed using a custom package developed by Tai-Hsien Ou Yang at Columbia University. The **concordanceIndex** package includes a single synonymous function that can be used as follows:

```R
concordanceIndex::concordanceIndex(predictions, observations)
```

... where `predictions` and `observations` are numerical vectors of the same length.

## Local Shiny-iAtlas Session

To run the app locally, clone this repository and use the following command in the `shiny-iatlas` directory:

```
shiny::runApp()
```

## To deploy, run this line

```
options(repos = BiocInstaller::biocinstallRepos())
```
