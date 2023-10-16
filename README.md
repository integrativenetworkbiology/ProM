# ProM

Please check the vignette about how to use the package.
To install the R package:
```R
library("devtools")
install_github("integrativenetworkbiology/ProM")
library("ProM")
```
Manuscript is currently under review.  
The R script runReconstruction_CV_pancancer_github.r was used to generate leave-one-out cross-validation results of reconstituted bulk RNAseq in pancancer dataset by ProM and other deconvolution methods. CibersortX was run separately using docker, and thus was not included.
