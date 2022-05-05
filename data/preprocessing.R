library(dplyr)
library(bigPint)
library(SummarizedExperiment)
library(DelayedArray)

data("se_soybean_cn_sub")


matrix <- seed(assay(se_soybean_cn_sub))
