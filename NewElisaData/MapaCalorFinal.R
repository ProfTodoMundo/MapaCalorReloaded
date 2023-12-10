#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
setwd("~/Desktop/MiGithub/MapaCalorReloaded") # computadora de la casa
#setwd("~/Documentos/MiGitHub/MapaCalorReloaded") # computadora del trabajo
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
# Creación de gráficas
boxplot(NewGenes_Clean, las = 3)
boxplot(NewGenes_log2, las = 3)
boxplot(NewGenes_log2_filtrado, las = 3)
pheatmap(NewGenes_log2_filtrado, scale = "row",
         color = my_colors4, 
         border_color = NA, 
         fontsize_row = 5,
         display_numbers = FALSE,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         legend_breaks = c(-1,0,1),
         legend_labels = c("Bajo","Medio","Alto"))#,
#         kmeans_k = 50, 
#         cellheight = 10)
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 

pdf("Lab/HeatMapNormalizadosRenglon.pdf")
pheatmap(NewGenes_log2_filtrado, scale = "row")
dev.off()

pdf("Lab/HeatMapNormalizadosColumna.pdf")
pheatmap(NewGenes_log2_filtrado, scale = "column")
dev.off()

pdf("Lab/HeatMapPaletaPersonal.pdf")
pheatmap(NewGenes_log2_filtrado, scale = "row", color = my_colors2)
dev.off()

my_colors = brewer.pal(n = 11, name = "RdBu")
my_colors = colorRampPalette(my_colors)(50)
my_colors = rev(my_colors)
#my_colors

pdf("Lab/HeatMapLog2Transformed.pdf")
pheatmap(NewGenes_log2_filtrado, scale = "row", color = my_colors)
dev.off()

pdf("Lab/HeatMapLog2TransformedLetra4.pdf")
pheatmap(NewGenes_log2_filtrado, scale = "row",color = my_colors, border_color = NA, fontsize_row = 4)
dev.off()

pdf("Lab/MapadeCalorFinal.pdf")
pheatmap(NewGenes_log2_filtrado, scale = "row",color = my_colors, border_color = NA, fontsize_row = 6)
dev.off()
# <<>><<>> <<>><<>> <<>><<>> <<>><<>> <<>><<>> <<>><<>> <<>><<>> <<>><<>> <<>><<>>
# llegamos por fin llegamos
# <<>><<>> <<>><<>> <<>><<>> <<>><<>> <<>><<>> <<>><<>> <<>><<>> <<>><<>> <<>><<>>
pdf("Lab/MapadeCalorFinalSDC.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE)
dev.off()

pdf("Lab/HeathmapCFSinNorm.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE)
dev.off()

#scale = "row",
# A partir de aqui
pdf("Lab/HeathmapCFNormRenglon.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,scale = "row")
dev.off()


#fontsize_row = 6
pdf("Lab/HeathmapCFNormRenglonLetra.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,
         scale = "row",fontsize_row = 2.75)
dev.off()

#border_color = NA, 
pdf("Lab/HeathmapCFNormRenglonLetra8SinBordes.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,
         scale = "row",fontsize_row = 2.75, border_color = NA)
dev.off()


#color = my_colors, 
pdf("Lab/HeathmapCFNormRenglonLetra8SinBordesMycolors.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,
         scale = "row",fontsize_row = 2.75, border_color = NA)
dev.off()


#color = my_colors, 
pdf("Lab/HeathmapCFNormRenglonLetra8conBordesMycolors.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,
         scale = "row",fontsize_row = 2.75, main = "Mi Mapa de Calor")
dev.off()

