if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("igraph", quietly = TRUE))
  install.packages("igraph")

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")

if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

if (!require("DOSE", quietly = TRUE))
  BiocManager::install("DOSE")

if (!require("STRINGdb", quietly = TRUE))
  BiocManager::install("STRINGdb")

if (!require("linkcomm", quietly = TRUE))
  BiocManager::install("linkcomm")

if (!require("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")


if (!require("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")

