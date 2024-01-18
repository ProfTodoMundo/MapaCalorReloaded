#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
setwd("~/Desktop/MiGithub/MapaCalorReloaded/DraElisaNewData")
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
library(RColorBrewer)
library(readr)
library(pheatmap)
library(dplyr)
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
NewGenes <- read.csv("ExpresionEiMybs.csv")
colnames(NewGenes) <- c("GeneID","InputOrder","Trophozoites","Encyst_8h",
                        "Encyst_24h","Encyst_48h","Encyst_72h","Excyst_2h","Excyst_8h")
View(NewGenes)
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
rownames(NewGenes) <- NewGenes$GeneID
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
NewGenes_Clean <- NewGenes[,3:9]; View(NewGenes_Clean)
NewGenes_Clean <- as.data.frame(NewGenes_Clean)
rownames(NewGenes_Clean) <- NewGenes$GeneID; View(NewGenes_Clean)
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# Generacion de las bdd modificadas
NewGenes_log2   <-  log2(NewGenes_Clean + 1)
NewGenesOrdered <- NewGenes_Clean[order(-NewGenes_Clean$Trophozoites),]
top_genes_NewGenes_Trop <- NewGenesOrdered %>% filter(Trophozoites>0) %>% select(Trophozoites)
k  <- dim(NewGenes_Clean); proporcion <- 1; 
NS <- round(k[1]*proporcion);
random_genes_NewGenes  <- sample(rownames(NewGenes_Clean),NS );
sampledNewGenes_Log2   <- NewGenes_log2[random_genes_NewGenes, ];
NewGenes_log2_filtrado <- NewGenes_log2[rowSums(NewGenes_log2) != 0, ];View(NewGenes_log2_filtrado)
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
colores <- c("red", "blue","yellow", "green", "grey","orange", "purple")

pdf("NewMapasCalor/BoxPlotNewGenes.pdf")
boxplot(NewGenes_Clean, las = 3, col = colores)
legend("topleft", legend = c("Trophozoites", "Encyst 8h",
                             "Encyst 24h", "Encyst 48h", "Encyst 72h", "Excyst 2h", "Excyst 8h"),
       fill = colores, title = "EiMybs Genes")
dev.off()

pdf("NewMapasCalor/BoxPlotLog2NewGenes.pdf")
boxplot(NewGenes_log2, las = 2, col = colores, ylab="Log2(EiMybs)",
        names = c("Trophozoites", "Encyst 8h","Encyst 24h", "Encyst 48h",
                  "Encyst 72h", "Excyst 2h", "Excyst 8h"),
        cex.axis = 0.8)
legend("topleft", legend = c("Trophozoites", "Encyst 8h",
                             "Encyst 24h", "Encyst 48h", "Encyst 72h", "Excyst 2h", "Excyst 8h"),
       fill = colores, title = "EiMybs Genes")
dev.off()

pdf("NewMapasCalor/BoxPlotFilteredLog2NewGenes.pdf")
boxplot(NewGenes_log2_filtrado, las = 3, col = colores, ylab="Log2(EiMybs)",
        names = c("Trophozoites", "Encyst 8h","Encyst 24h", "Encyst 48h",
                  "Encyst 72h", "Excyst 2h", "Excyst 8h"),
        cex.axis = 0.8)
legend("topleft", legend = c("Trophozoites", "Encyst 8h",
                             "Encyst 24h", "Encyst 48h", "Encyst 72h", "Excyst 2h", "Excyst 8h"),
       fill = colores, title = "EiMybs Genes")
dev.off()
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
LetraSize <- 6
pdf("NewMapasCalor/HeathmapCFNormRenglonSinBordes.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,
         scale = "row",fontsize_row = LetraSize, border_color = NA)
dev.off()
#
pdf("NewMapasCalor/HeathmapCFNormRenglonSinBordesMycolors.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,
         scale = "row",fontsize_row = LetraSize, border_color = NA)
dev.off()
#
pdf("NewMapasCalor/HeathmapCFNormRenglonSinBordesMycolors2.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,
         scale = "row",fontsize_row = LetraSize, border_color = NA)
dev.off()
#
pdf("NewMapasCalor/HeathmapCFNormRenglonSinBordesMycolors3.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,
         scale = "row",fontsize_row = LetraSize, border_color = NA)
dev.off()
#
pdf("NewMapasCalor/HeathmapCFNormRenglonSinBordesMycolors4.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors4,
         scale = "row",fontsize_row = LetraSize, border_color = NA)
dev.off()
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 