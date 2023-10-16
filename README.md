# ProM


We developed our novel bulk RNAseq deconvolution algorithm, named ProM, by combining two classical deconvolution strategies, i.e., marker-based and profile-based. Leveraging the complementary nature of the two strategies together with a few other technical improvements, ProM showed attractive advantages over multiple existing deconvolution methods by being more accurate and faster. However, there is no “one-size-fits-all” approach in the field of computational deconvolution. Different methods have been found to excel in different scenarios.  We hope to continue to assess the utility of ProM in other applications in future studies. 

Please check the vignette about how to use the package.
To install the ProM package:
```R
library("devtools")
install_github("integrativenetworkbiology/ProM")
library("ProM")
```
Manuscript is currently under review.  
The R script runReconstruction_CV_pancancer_github.r was used to generate leave-one-out cross-validation results of reconstituted bulk RNAseq in pancancer dataset by ProM and other deconvolution methods. CibersortX was run separately using docker, and thus was not included.
