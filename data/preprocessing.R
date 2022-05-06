library(dplyr)
library(bigPint)
library(SummarizedExperiment)
library(DelayedArray)

data("se_soybean_cn_sub")

# pull normalized counts matrix
matrix <- as.matrix(seed(assay(se_soybean_cn_sub)))
write.table(matrix, file='normalized_counts.csv', sep=',')

# generate sample info matrix
array_id <- c("S1.1", "S1.2", "S1.3", "S2.1", "S2.2", "S2.3", "S3.1", "S3.2", "S3.3")
stage <- c("early", "early", "early", "middle", "middle", "middle", "late", "late", "late")
replicate_num <- c(1 ,2, 3, 1, 2, 3, 1, 2, 3)

sample_info <- data.frame(array_id, stage, replicate_num)
write.table(sample_info, file='sample_info.csv', sep=',')

de_results <- rowData(se_soybean_cn_sub)

write.table(de_results, file="./de_results.csv", sep=",", row.names = TRUE)

