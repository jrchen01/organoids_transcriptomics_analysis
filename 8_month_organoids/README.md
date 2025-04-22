This folder contains code for the 8‑month organoids scRNA‑seq analysis, including data distribution checks, clustering and annotation, differential expression analysis, and cell–cell communication analysis.

## Software environment
All analyses were performed in R (version 4.2.2) under CentOS Linux 7 (Core) using conda environment. The complete output of `sessionInfo()` is shown below:
This folder contains code for the 8‑month organoids scRNA‑seq analysis, including data distribution checks, clustering and annotation, differential expression analysis, and cell–cell communication analysis.

## Software environment

All analyses were performed in R (version 4.2.2) under CentOS Linux 7 (Core) using a conda-managed environment. For full reproducibility, the complete output of `sessionInfo()` is provided below:

```r
sessionInfo()
R version 4.2.2 (2022-10-31)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /raidixshare_logg01/jchen/miniconda3/envs/r.4.2.2/lib/libopenblasp-r0.3.27.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
 [1] stats4    stats     graphics  grDevices utils     datasets  methods  
 [8] base     

other attached packages:
 [1] ashr_2.2-63                 lubridate_1.9.3            
 [3] forcats_1.0.0               stringr_1.5.1              
 [5] purrr_1.0.4                 readr_2.1.5                
 [7] tidyr_1.3.1                 tibble_3.2.1               
 [9] tidyverse_2.0.0             DESeq2_1.38.0              
[11] SingleCellExperiment_1.20.0 SummarizedExperiment_1.28.0
[13] Biobase_2.58.0              GenomicRanges_1.50.0       
[15] GenomeInfoDb_1.34.9         IRanges_2.32.0             
[17] S4Vectors_0.36.0            BiocGenerics_0.44.0        
[19] MatrixGenerics_1.10.0       matrixStats_1.4.1          
[21] presto_1.0.0                data.table_1.15.4          
[23] Rcpp_1.0.14                 devtools_2.4.5             
[25] usethis_2.2.3               ggrepel_0.9.6              
[27] ggplot2_3.5.1               dplyr_1.1.4                
[29] patchwork_1.3.0             Seurat_5.1.0               
[31] SeuratObject_5.0.2          sp_2.1-4                   

loaded via a namespace (and not attached):
  [1] spatstat.explore_3.2-6 reticulate_1.40.0      tidyselect_1.2.1      
  [4] RSQLite_2.3.4          AnnotationDbi_1.60.0   htmlwidgets_1.6.4     
  [7] grid_4.2.2             BiocParallel_1.32.5    Rtsne_0.17            
 [10] munsell_0.5.1          codetools_0.2-20       ica_1.0-3             
 [13] pbdZMQ_0.3-11          future_1.34.0          miniUI_0.1.1.1        
 [16] withr_3.0.2            spatstat.random_3.2-3  colorspace_2.1-1      
 [19] progressr_0.14.0       uuid_1.2-0             ROCR_1.0-11           
 [22] tensor_1.5             listenv_0.9.1          repr_1.1.7            
 [25] GenomeInfoDbData_1.2.9 mixsqp_0.3-54          polyclip_1.10-6       
 [28] bit64_4.0.5            farver_2.1.2           parallelly_1.42.0     
 [31] vctrs_0.6.5            generics_0.1.3         timechange_0.3.0      
 [34] R6_2.6.1               invgamma_1.1           locfit_1.5-9.9        
 [37] bitops_1.0-7           spatstat.utils_3.0-5   cachem_1.1.0          
 [40] DelayedArray_0.24.0    promises_1.3.0         scales_1.3.0          
 [43] gtable_0.3.6           globals_0.16.3         goftest_1.2-3         
 [46] spam_2.10-0            rlang_1.1.5            genefilter_1.80.0     
 [49] splines_4.2.2          lazyeval_0.2.2         spatstat.geom_3.2-9   
 [52] reshape2_1.4.4         abind_1.4-8            httpuv_1.6.15         
 [55] tools_4.2.2            ellipsis_0.3.2         RColorBrewer_1.1-3    
 [58] sessioninfo_1.2.2      ggridges_0.5.6         plyr_1.8.9            
 [61] base64enc_0.1-3        zlibbioc_1.44.0        RCurl_1.98-1.12       
 [64] deldir_2.0-4           pbapply_1.7-2          cowplot_1.1.3         
 [67] urlchecker_1.0.1       zoo_1.8-12             cluster_2.1.6         
 [70] fs_1.6.4               magrittr_2.0.3         RSpectra_0.16-2       
 [73] scattermore_1.2        lmtest_0.9-40          RANN_2.6.1            
 [76] truncnorm_1.0-9        SQUAREM_2021.1         fitdistrplus_1.1-11   
 [79] pkgload_1.3.4          hms_1.1.3              mime_0.12             
 [82] evaluate_0.24.0        xtable_1.8-4           XML_3.99-0.14         
 [85] fastDummies_1.7.3      gridExtra_2.3          compiler_4.2.2        
 [88] KernSmooth_2.23-24     crayon_1.5.3           htmltools_0.5.8.1    
 [91] later_1.3.2            tzdb_0.4.0             geneplotter_1.76.0   
 [94] DBI_1.2.3              MASS_7.3-60.0.1        Matrix_1.6-5          
 [97] cli_3.6.4              parallel_4.2.2         dotCall64_1.1-1       
 [100] igraph_1.4.2           pkgconfig_2.0.3        IRdisplay_1.1          
 [103] plotly_4.10.4          spatstat.sparse_3.1-0  annotate_1.76.0       
 [106] XVector_0.38.0         digest_0.6.37          sctransform_0.4.1     
 [109] RcppAnnoy_0.0.22       spatstat.data_3.1-2    Biostrings_2.66.0     
 [112] leiden_0.4.3.1         uwot_0.1.16            shiny_1.8.1.1         
 [115] lifecycle_1.0.4        nlme_3.1-165           jsonlite_1.9.0        
 [118] viridisLite_0.4.2      pillar_1.10.2          lattice_0.22-6        
 [121] KEGGREST_1.38.0        fastmap_1.2.0          httr_1.4.7            
 [124] pkgbuild_1.4.4         survival_3.7.0         glue_1.8.0            
 [127] remotes_2.5.0          png_0.1-8              bit_4.0.5            
 [130] stringi_1.8.4          profvis_0.3.8         blob_1.2.4            
 [133] RcppHNSW_0.6.0         memoise_2.0.1          IRkernel_1.3.2        
 [136] irlba_2.3.5.1          future.apply_1.11.3
