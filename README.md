Smmit: Multi-sample single-cell multi-omics integration
====

## Overview
Smmit performs integration both across samples and modalities to produce a single UMAP space. It first uses harmony to integrate across samples and then uses Seurat weighted nearest neighbor function to integrate across modalities.


## Smmit Installation

`Smmit` software can be installed via Github.
Users should have R installed on their computer before installing `Smmit. R` can be downloaded here: http://www.r-project.org/.
To install the latest version of `Smmit` package via Github, run following commands in R:
```{r }
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("zji90/Smmit")
```

## User Manual
Please visit this webpage for the user manual: https://github.com/zji90/Smmit/blob/master/vignettes/Smmit.pdf

## Contact the Author
Author: Changxin Wan, Zhicheng Ji

Report bugs and provide suggestions by sending email to:

Maintainer: Zhicheng Ji (zhicheng.ji@duke.edu)

Or open a new issue on this Github page

