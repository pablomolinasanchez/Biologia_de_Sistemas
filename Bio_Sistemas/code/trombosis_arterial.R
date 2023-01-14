suppressMessages(library(DOSE))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(devtools))
suppressMessages(library(dplyr))
suppressMessages(library(igraph))
suppressMessages(library(ggplot2))
suppressMessages(library(linkcomm))


rm(list=ls())

nodes_prop<-read.csv("string_node_degrees_prop.tsv",sep ="")
links_prop<-read.csv("string_interactions_prop.tsv",sep="")
net_prop <- graph_from_data_frame(d=links_prop,vertices=nodes_prop, directed=T)

# Comunidades
lc<- igraph::as_data_frame(net_prop, what="edges") 
lc<-lc[,c(1,2,13)]
lc

hits.network_lc <- getLinkCommunities(lc , hcmethod = "single")
print(hits.network_lc)
png("../results/Análisis_Comunidades/01_NetworkComunidades.png")
plot(hits.network_lc, type = "graph", layout = layout.fruchterman.reingold, ewidth = 2, vlabel.cex = 0.5)
dev.off()

## Comunidades mas grandes
png("../results/Análisis_Comunidades/02_ComunidadesMayorTamano.png")
par(mfrow = c(2,1))
pie(hits.network_lc$clustsizes[hits.network_lc$clustsizes > 8], radius = 1, main = "Tamaños de las comunidades m?s grandes")
barplot(hits.network_lc$clustsizes[hits.network_lc$clustsizes > 8], xlab = "Comunidades", ylab = "Tamaño (num genes)")
dev.off()
par(mfrow = c(1,1))

## Genes mas conectados
png("../results/Análisis_Comunidades/03_GenesMasCOnectadosOtrasComunidades.png")
plot(hits.network_lc, type = "members")
dev.off()


## Importancia seguncentralidad
cc <- getCommunityCentrality(hits.network_lc)
print(head(sort(cc, decreasing = TRUE)))
png("../results/Análisis_Comunidades/04_GenesMayorCentralidad.png")
barplot(t(head(sort(cc, decreasing= TRUE))), xlab = "genes", ylab = "Community centrality", main = "Genes con mayor CC", ylim = c(0,10))




# Comunidades independientes
getAllNestedComm(hits.network_lc)
png("../results/Análisis_Comunidades/05_Comunidades_Independientes.png")
plot(hits.network_lc, type = "graph", clusterids = c(9, 15))
dev.off()

## Funcion enriquecimiento
data(geneList, package="DOSE")
overrep_enrichment_GO <- function(gene,universe,Db, ont, pAdjustMethod, p,q) {
  ego <- enrichGO(gene          = gene,
                  universe      = universe,
                  OrgDb         = Db,
                  ont           = ont,
                  pAdjustMethod = pAdjustMethod,
                  pvalueCutoff  = p,
                  qvalueCutoff  = q,
                  readable      = TRUE)
  return (ego)
}


## Comunidad 1
community_1<-getNodesIn(hits.network_lc, clusterids = c(1))
community_1<-bitr(community_1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_1<-overrep_enrichment_GO(community_1$ENTREZID,geneList,org.Hs.eg.db,"BP","BH",0.05,0.05)
dir <- paste("../results/Enriquecimientos/Análisis_funcional_Comunidades_1.csv", sep = "")
write.csv(go_1, dir)
head(go_1)


## Comunidad 2
community_2<-getNodesIn(hits.network_lc, clusterids = c(2))
community_2<-bitr(community_2, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_2<-overrep_enrichment_GO(community_2$ENTREZID,geneList,org.Hs.eg.db,"BP","BH",0.05,0.05)
dir <- paste("../results/Enriquecimientos/Análisis_funcional_Comunidades_2.csv", sep = "")
write.csv(go_2, dir)
head(go_2)


## Comunidad 3
community_3<-getNodesIn(hits.network_lc, clusterids = c(3))
community_3<-bitr(community_3, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_3<-overrep_enrichment_GO(community_3$ENTREZID,geneList,org.Hs.eg.db,"BP","BH",0.05,0.05)
dir <- paste("../results/Enriquecimientos/Análisis_funcional_Comunidades_3.csv", sep = "")
write.csv(go_3, dir)
head(go_3)


## Comunidad 4
community_4<-getNodesIn(hits.network_lc, clusterids = c(4))
community_4<-bitr(community_4, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_4<-overrep_enrichment_GO(community_4$ENTREZID,geneList,org.Hs.eg.db,"BP","BH",0.05,0.05)
dir <- paste("../results/Enriquecimientos/Análisis_funcional_Comunidades_4.csv", sep = "")
write.csv(go_4, dir)
head(go_4)


## Comunidad 5
community_5<-getNodesIn(hits.network_lc, clusterids = c(5))
community_5<-bitr(community_5, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_5<-overrep_enrichment_GO(community_5$ENTREZID,geneList,org.Hs.eg.db,"BP","BH",0.05,0.05)
dir <- paste("../results/Enriquecimientos/Análisis_funcional_Comunidades_5.csv", sep = "")
write.csv(go_5, dir)
head(go_5)


## Comunidad 6
community_6<-getNodesIn(hits.network_lc, clusterids = c(6))
community_6<-bitr(community_6, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_6<-overrep_enrichment_GO(community_6$ENTREZID,geneList,org.Hs.eg.db,"BP","BH",0.05,0.05)
dir <- paste("../results/Enriquecimientos/Análisis_funcional_Comunidades_6.csv", sep = "")
write.csv(go_6, dir)
head(go_6)


## Comunidad 7
community_7<-getNodesIn(hits.network_lc, clusterids = c(7))
community_7<-bitr(community_7, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_7<-overrep_enrichment_GO(community_7$ENTREZID,geneList,org.Hs.eg.db,"BP","BH",0.05,0.05)
dir <- paste("../results/Enriquecimientos/Análisis_funcional_Comunidades_7.csv", sep = "")
write.csv(go_7, dir)
head(go_7)


## Comunidad 8
community_8<-getNodesIn(hits.network_lc, clusterids = c(8))
community_8<-bitr(community_8, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_8<-overrep_enrichment_GO(community_8$ENTREZID,geneList,org.Hs.eg.db,"BP","BH",0.05,0.05)
dir <- paste("../results/Enriquecimientos/Análisis_funcional_Comunidades_8.csv", sep = "")
write.csv(go_8, dir)
head(go_8)


## Comunidad 9
community_9<-getNodesIn(hits.network_lc, clusterids = c(9))
community_9<-bitr(community_9, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_9<-overrep_enrichment_GO(community_9$ENTREZID,geneList,org.Hs.eg.db,"BP","BH",0.05,0.05)
dir <- paste("../results/Enriquecimientos/Análisis_funcional_Comunidades_9.csv", sep = "")
write.csv(go_9, dir)
head(go_9)


## Comunidad 10
community_10<-getNodesIn(hits.network_lc, clusterids = c(10))
community_10<-bitr(community_10, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_10<-overrep_enrichment_GO(community_10$ENTREZID,geneList,org.Hs.eg.db,"BP","BH",0.05,0.05)
dir <- paste("../results/Enriquecimientos/Análisis_funcional_Comunidades_10.csv", sep = "")
write.csv(go_10, dir)
head(go_10)


## Comunidad 11
community_11<-getNodesIn(hits.network_lc, clusterids = c(11))
community_11<-bitr(community_11, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_11<-overrep_enrichment_GO(community_11$ENTREZID,geneList,org.Hs.eg.db,"BP","BH",0.05,0.05)
dir <- paste("../results/Enriquecimientos/Análisis_funcional_Comunidades_11.csv", sep = "")
write.csv(go_11, dir)
head(go_11)


## Comunidad 12
community_12<-getNodesIn(hits.network_lc, clusterids = c(12))
community_12<-bitr(community_12, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_12<-overrep_enrichment_GO(community_12$ENTREZID,geneList,org.Hs.eg.db,"BP","BH",0.05,0.05)
dir <- paste("../results/Enriquecimientos/Análisis_funcional_Comunidades_12.csv", sep = "")
write.csv(go_12, dir)
head(go_12)


## Comunidad 13
community_13<-getNodesIn(hits.network_lc, clusterids = c(13))
community_13<-bitr(community_13, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_13<-overrep_enrichment_GO(community_13$ENTREZID,geneList,org.Hs.eg.db,"BP","BH",0.05,0.05)
dir <- paste("../results/Enriquecimientos/Análisis_funcional_Comunidades_13.csv", sep = "")
write.csv(go_13, dir)
head(go_13)


## Comunidad 14
community_14<-getNodesIn(hits.network_lc, clusterids = c(14))
community_14<-bitr(community_14, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_14<-overrep_enrichment_GO(community_14$ENTREZID,geneList,org.Hs.eg.db,"BP","BH",0.05,0.05)
dir <- paste("../results/Enriquecimientos/Análisis_funcional_Comunidades_14.csv", sep = "")
write.csv(go_14, dir)
head(go_14)


## Comunidad 15
community_15<-getNodesIn(hits.network_lc, clusterids = c(15))
community_15<-bitr(community_15, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_15<-overrep_enrichment_GO(community_15$ENTREZID,geneList,org.Hs.eg.db,"BP","BH",0.05,0.05)
dir <- paste("../results/Enriquecimientos/Análisis_funcional_Comunidades_15.csv", sep = "")
write.csv(go_15, dir)
head(go_15)


## Comunidad 16
community_16<-getNodesIn(hits.network_lc, clusterids = c(16))
community_16<-bitr(community_16, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_16<-overrep_enrichment_GO(community_16$ENTREZID,geneList,org.Hs.eg.db,"BP","BH",0.05,0.05)
dir <- paste("../results/Enriquecimientos/Análisis_funcional_Comunidades_16.csv", sep = "")
write.csv(go_16, dir)
head(go_16)


## Comunidad 17
community_17<-getNodesIn(hits.network_lc, clusterids = c(17))
community_17<-bitr(community_17, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_17<-overrep_enrichment_GO(community_17$ENTREZID,geneList,org.Hs.eg.db,"BP","BH",0.05,0.05)
dir <- paste("../results/Enriquecimientos/Análisis_funcional_Comunidades_17.csv", sep = "")
write.csv(go_17, dir)
head(go_17)


## Comunidad 18
community_18<-getNodesIn(hits.network_lc, clusterids = c(18))
community_18<-bitr(community_18, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_18<-overrep_enrichment_GO(community_18$ENTREZID,geneList,org.Hs.eg.db,"BP","BH",0.05,0.05)
dir <- paste("../results/Enriquecimientos/Análisis_funcional_Comunidades_18.csv", sep = "")
write.csv(go_18, dir)
head(go_18)

