#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#setwd("~/Desktop/MiGithub/MapaCalorReloaded") # computadora de la casa
setwd("~/Documentos/MiGitHub/MapaCalorReloaded") # computadora del trabajo
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
library(RColorBrewer)
library(readr)
library(pheatmap)
library(dplyr)
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
NewGenes <- read_csv("Ei_enq_circ.csv")
View(NewGenes)
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
NewGenes$Trophs      <- NewGenes$Trophs-2
NewGenes$`cyst 8h`   <- NewGenes$`cyst 8h`-2;
NewGenes$`cyst 24h`  <- NewGenes$`cyst 24h`-2
NewGenes$`cyst 48h`  <- NewGenes$`cyst 48h`-2
NewGenes$`cyst 72h`  <- NewGenes$`cyst 72h`-2
NewGenes$`excyst 2h` <- NewGenes$`excyst 2h`-2
NewGenes$`excyst 8h` <- NewGenes$`excyst 8h`-2
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
colnames(NewGenes) <- c("GeneID","Trophozoites", "Cyst_8h","Cyst_24h",
                        "Cyst_48h","Cyst_72h","Excyst_2h","Excyst_8h")
rownames(NewGenes) <- NewGenes$GeneID
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
NewGenes_Clean <- NewGenes[,2:8]; View(NewGenes_Clean)
NewGenes_Clean <- as.data.frame(NewGenes_Clean)
rownames(NewGenes_Clean) <- NewGenes$GeneID
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# Generacion de las bdd modificadas
NewGenes_log2   <-  log2(NewGenes_Clean + 1)
NewGenesOrdered <- NewGenes_Clean[order(-NewGenes_Clean$Trophozoites),]
top_genes_NewGenes_Trop <- NewGenesOrdered %>% filter(Trophozoites>0) %>% select(Trophozoites)
k  <- dim(NewGenes_Clean); proporcion <- 1; 
NS <- round(k[1]*proporcion);
random_genes_NewGenes  <- sample(rownames(NewGenes_Clean),NS );
sampledNewGenes_Log2   <- NewGenes_log2[random_genes_NewGenes, ];
NewGenes_log2_filtrado <- NewGenes_log2[rowSums(NewGenes_log2) != 0, ];
sampledNewGenes_Log2   <- NewGenes_log2[random_genes_NewGenes, ];
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# Personalización de la paleta de colores
my_colors = brewer.pal(n = 11, name = "RdBu")
my_colors = colorRampPalette(my_colors)(50)
my_colors = rev(my_colors)
my_colors2 = c("green", "yellow", "pink")
my_colors2 = colorRampPalette(my_colors2)(50)
my_colors3 = brewer.pal(n = 11, name = "RdBu")
my_colors3 = colorRampPalette(my_colors3)(50)
my_colors3 = rev(my_colors3)

my_colors4 = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# Guardado de las bdd
write.csv(NewGenes_log2,"NewGenesLog.csv")
write.csv(sampledNewGenes_Log2,"SampledNewGenesLog.csv")
write.csv(NewGenes_log2_filtrado,"NewGenesLogFiltrado.csv")
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# Creación de gráficas
boxplot(NewGenes_Clean, las = 3)
boxplot(NewGenes_log2, las = 3)

pheatmap(NewGenes_Clean[random_genes_NewGenes, ])
pheatmap(NewGenes_log2[random_genes_NewGenes, ])

pheatmap(NewGenes_log2, scale = "row")
#pheatmap(NewGenes_log2, scale = "column")

pheatmap(NewGenes_log2_filtrado, scale = "row")
#pheatmap(NewGenes_log2_filtrado, scale = "column")

pheatmap(NewGenes_log2_filtrado, scale = "row", color = my_colors)
pheatmap(NewGenes_log2_filtrado, scale = "row", color = my_colors2)
pheatmap(NewGenes_log2_filtrado, scale = "row", color = my_colors3)

pheatmap(NewGenes_log2_filtrado, scale = "row",color = my_colors, border_color = NA, fontsize_row = 6)
pheatmap(NewGenes_log2_filtrado, scale = "row",color = my_colors2,border_color = NA, fontsize_row = 4)
pheatmap(NewGenes_log2_filtrado, scale = "row",color = my_colors3,border_color = NA, fontsize_row = 6)
#pheatmap(NewGenes_log2_filtrado, scale = "column",color = my_colors3,border_color = NA, fontsize_row = 6)


pheatmap(NewGenes_log2, scale = "row", color = my_colors)
pheatmap(NewGenes_log2, scale = "row", color = my_colors2)
pheatmap(NewGenes_log2, scale = "row", color = my_colors3)

pheatmap(NewGenes_log2, scale = "row",color = my_colors3,border_color = NA, fontsize_row = 4)
pheatmap(NewGenes_log2, scale = "row",color = my_colors3,border_color = NA, fontsize_row = 6)
#pheatmap(NewGenes_log2, scale = "column",color = my_colors3,border_color = NA, fontsize_row = 6)
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#pdf("Boxplot.#pdf")
#boxplot(NewGenes_log2, las = 3)
#dev.off()
#pdf("HeathMapNewGenesDef.#pdf")
#pheatmap(NewGenes_log2[random_genes_NewGenes, ])
#dev.off()
#pdf("HeathMapLog2NewGenesDef.#pdf")
#pheatmap(NewGenes_log2_filtrado, scale = "row")
#dev.off()
#pdf("HeathMapLog2NewGenesfiltrado.#pdf")
#pheatmap(NewGenes_log2_filtrado, scale = "row",color = my_colors, border_color = NA, fontsize_row = 6)
#dev.off()
#pdf("HeatMapDatosLimpios.#pdf")
#pheatmap(NewGenes_Clean[random_genes_NewGenes, ])
#dev.off()
#pdf("HeatMapLog2Transformed.#pdf")
#pheatmap(NewGenes_log2[random_genes_NewGenes, ])
#dev.off()
#pdf("HeatMapNormalizadosRenglon.#pdf")
#pheatmap(NewGenes_log2_filtrado, scale = "row")
#dev.off()
#pdf("HeatMapNormalizadosColumna.#pdf")
#pheatmap(NewGenes_log2_filtrado, scale = "column")
#dev.off()
#pdf("HeatMapPaletaPersonal.#pdf")
#pheatmap(NewGenes_log2_filtrado, scale = "row", color = my_colors)
#dev.off()
#pdf("HeatMapPaletaPersonal.#pdf")
#pheatmap(NewGenes_log2_filtrado, scale = "row", color = my_colors2)
#dev.off()
#pdf("HeatMapLog2Transformed.#pdf")
#pheatmap(NewGenes_log2_filtrado, scale = "row", color = my_colors)
#dev.off()
#pdf("HeatMapLog2TransformedLetra4.#pdf")
#pheatmap(NewGenes_log2_filtrado, scale = "row",color = my_colors,border_color = NA, fontsize_row = 4)
#dev.off()
#pdf("MapadeCalorFinal.#pdf")
#pheatmap(NewGenes_log2_filtrado, scale = "row",color = my_colors,border_color = NA, fontsize_row = 6)
#dev.off()
#pheatmap(NewGenes_log2_filtrado, scale = "column",color = my_colors,border_color = NA, fontsize_row = 6)
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
pheatmap(NewGenes_log2_filtrado, scale = "row",
         color = my_colors4, 
         border_color = NA, 
         fontsize_row = 6,
         cluster_cols = TRUE)


# == << >> == << >> ==  == << >> == << >> ==  == << >> == << >> ==  == << >> == << >> ==  == << >> ==
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#my_colors
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 





# == << >> == << >> ==  == << >> == << >> ==  == << >> == << >> ==  == << >> == << >> ==  == << >> ==
# == << >> == << >> ==  == << >> == << >> ==  == << >> == << >> ==  == << >> == << >> ==  == << >> ==
NewGenes_Log2 <- NewGenes_log2_filtrado
hc <- hclust(dist(NewGenes_Log2))
grupos <- cutree(hc, k = 1); colores <- c("purple", "blue");  
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#pdf("DendogramNewGenesJerarquico.#pdf")
plot(hc, main = "Dendograma del análisis de clúster jerárquico", 
     col = colores[grupos], 
     xlab =  "Genes NewGenes",
     ylab = "Distancias",
     hang = -1)
#dev.off()
#
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
hc <- hclust(dist(NewGenes_Log2),"ave")
#pdf("DendogramNewGenesJerarquicoAve.#pdf")
plot(hc, main = "Dendograma del análisis de clúster jerárquico",
     xlab =  "Genes NewGenes",
     ylab = "Distancias",
     hang = -1,
     col=colores[grupos])
#dev.off()
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
hc <- hclust(dist(NewGenes_Log2),"cen")
#pdf("DendogramNewGenesJerarquicoCen.#pdf")
plot(hc, main = "Dendograma del análisis de clúster jerárquico",
     xlab =  "Genes NewGenes",
     ylab = "Distancias",
     hang = -1,
     col=colores[grupos])
#dev.off()
# == << >> == << >> ==  == << >> == << >> ==  == << >> == << >> ==  == << >> == << >> ==  == << >> ==
k <- 3
kmeans_result <- kmeans(NewGenes_Log2, centers = k)
#pdf("KMeansGraph.#pdf")
plot(NewGenes_Log2, col = kmeans_result$cluster)
#dev.off()
variaciones_explicadas <- vector("numeric", length = 10)
for (k in 1:10) {
  kmeans_result <- kmeans(NewGenes_Log2, centers = k)
  variaciones_explicadas[k] <- kmeans_result$tot.withinss}
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#pdf("Kmeansto10.#pdf")
plot(1:10, variaciones_explicadas, type = "b", pch = 19, frame = FALSE,
     xlab = "Número de clústeres", ylab = "Variación explicada",main = "Método del codo")
ss_total <- sum(var(NewGenes_Log2)^2)
variacion_explicada_rel <- 1 - variaciones_explicadas/ss_total
lines(1:10, variacion_explicada_rel, type = "b", pch = 19, col = "red")
abline(v = which.max(variacion_explicada_rel), col = "blue", lty = 2)
k_optimo <- which.max(variacion_explicada_rel); k_optimo
#dev.off()
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
hc <- hclust(dist(NewGenes_Log2))
grupos <- cutree(hc, k = 10); colores <- c("purple", "blue");  
#pdf("DendogramNewGenesJerarquicoJerarq.#pdf")
plot(hc, main = "Dendograma del análisis de clúster jerárquico", 
     col = colores[grupos], xlab =  "Genes NewGenes",ylab = "Distancias",hang = -1)
#dev.off()
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
hc <- hclust(dist(NewGenes_Log2),"ave")
#pdf("DendogramNewGenesJerarquicoAve.#pdf")
plot(hc, main = "Dendograma del análisis de clúster jerárquico",
     xlab =  "Genes NewGenes",ylab = "Distancias",hang = -1,col=colores[grupos])
#dev.off()
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
hc <- hclust(dist(NewGenes_Log2),"cen")
#pdf("DendogramNewGenesJerarquicoCen.#pdf")
plot(hc, main = "Dendograma del análisis de clúster jerárquico",xlab =  "Genes NewGenes",
     ylab = "Distancias",hang = -1,col=colores[grupos])
#dev.off()
# == << >> == << >> ==  == << >> == << >> ==  == << >> == << >> ==  == << >> == << >> ==  == << >> ==
library(mclust)
# Aplicar el clustering basado en el algoritmo de mezcla de Gaussianas
mclust_result <- Mclust(NewGenes_Log2)
#pdf("ClusterGaussClass.#pdf")
plot(mclust_result, what = "classification",
     main = "Cluster basado en mezcla de gaussianos",
     hang = -1, col=my_colors)
#dev.off()
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#pdf("ClusterGaussDensity.#pdf")
plot(mclust_result, what = "density",
     main = "Cluster basado en mezcla de gaussianos",
     hang = -1,col=my_colors)
#dev.off()
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#pdf("ClusterGaussUncertainty.#pdf")
plot(mclust_result, what = "uncertainty",
     main = "Cluster basado en mezcla de gaussianos",
     hang = -1,col=my_colors)
#dev.off()
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#pdf("ClusterGaussBIC.#pdf")
plot(mclust_result, what = "BIC",
     main = "Cluster basado en mezcla de gaussianos",
     hang = -1,col=my_colors)
#dev.off()
# == << >> == << >> ==  == << >> == << >> ==  == << >> == << >> ==  == << >> == << >> ==  == << >> ==
