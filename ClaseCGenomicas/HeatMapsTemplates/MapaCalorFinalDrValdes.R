# ---- Limpiamos el entorno y la consola ----
rm(list = ls())
cat("\014")
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# ---- Se ajusta el directorio de trabajo ----
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
setwd("~/Documents/GitHub/MapaCalorReloaded/DrJesusValdes")
#setwd("~/Desktop/MiGithub/MapaCalorReloaded") # computadora de la casa
#setwd("~/Documentos/MiGitHub/MapaCalorReloaded") # computadora del trabajo
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# ---- Se cargan las librerias que se requieren ----
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
library(RColorBrewer)
library(readr)
library(pheatmap)
library(dplyr)
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# ---- Lectura de los grenes ----
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
NewGenes <- read_csv("Ei_enq_circ.csv")
View(NewGenes)
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# Reestablecimiento a su valor original, restando 2 a cadsa entrada ----
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# ---- Preprocesamiento de los datos ----
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
# ---- Generacion de las bdd modificadas ---- 
NewGenes_log2   <-  log2(NewGenes_Clean + 1)
NewGenesOrdered <- NewGenes_Clean[order(-NewGenes_Clean$Trophozoites),]
top_genes_NewGenes_Trop <- NewGenesOrdered %>% filter(Trophozoites>0) %>% select(Trophozoites)
k  <- dim(NewGenes_Clean); proporcion <- 1; 
NS <- round(k[1]*proporcion);
random_genes_NewGenes  <- sample(rownames(NewGenes_Clean),NS );
sampledNewGenes_Log2   <- NewGenes_log2[random_genes_NewGenes, ];
NewGenes_log2_filtrado <- NewGenes_log2[rowSums(NewGenes_log2) != 0, ]; View(NewGenes_log2_filtrado)
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# ---- Personalización de la paleta de colores ----
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
# ---- Creación de gráficas ----
## ---- Parametros 60 y 3.75 ----
cell_width <- 60
mi_font_size <- 3.75
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,scale = "row",
         fontsize_row = mi_font_size, cellwidth = cell_width,main ="Amplio de celda 60 y letras de genes 3.25")
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,main ="Amplio de celda 60 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"conBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 60 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 60 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 60 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 60 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
## ---- Parametros 65 y 4.25 ----
cell_width <- 65
mi_font_size <- 4.25
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,scale = "row",
         fontsize_row = mi_font_size, cellwidth = cell_width,main ="Amplio de celda 65 y letras de genes 4.25")
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,main ="Amplio de celda 65 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores1",".pdf")
pdf(nombregrafico)
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 65 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 65 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 65 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 65 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
## ---- Parametros 70 y 3.75 ----
cell_width <- 70
mi_font_size <- 3.75
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,scale = "row",
         fontsize_row = mi_font_size, cellwidth = cell_width,main ="Amplio de celda 70 y letras de genes 3.75")
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,main ="Amplio de celda 70 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 70 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 70 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 70 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 70 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
## ---- Parametros 75 y 3.75 ----
cell_width <- 75
mi_font_size <- 3.75
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,scale = "row",
         fontsize_row = mi_font_size, cellwidth = cell_width,main ="Amplio de celda 75 y letras de genes 3.75")
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,main ="Amplio de celda 75 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 75 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 75 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 75 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 75 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
## ---- Parametros 60 y 4.25 ----
cell_width <- 60
mi_font_size <- 4.25
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,scale = "row",
         fontsize_row = mi_font_size, cellwidth = cell_width,main ="Amplio de celda 60 y letras de genes 4.25")
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,main ="Amplio de celda 60 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 60 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 60 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 60 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/Heathmap","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 60 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/Boxplot",".pdf")
pdf(nombregrafico)
boxplot(NewGenes_log2_filtrado, las = 3)
dev.off()

