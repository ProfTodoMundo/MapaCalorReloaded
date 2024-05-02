# ---- Limpiamos el entorno y la consola ----
rm(list = ls())
cat("\014")
# ---- Se ajusta el directorio de trabajo ----
setwd("~/Documents/GitHub/MapaCalorReloaded/ClaseCGenomicas")
# ---- Se cargan las librerias que se requieren ----
library(RColorBrewer)
library(readr)
library(pheatmap)
library(dplyr)
# ---- Lectura de los grenes ----
ruta_archivo <- "MapaCalor.xlsx"
##---- Lectura de las hojas de excel incluidas en el archivo----
dataset1 <- read_excel(ruta_archivo, sheet = 1)
head(dataset1,5); View(dataset1)
dataset2 <- read_excel(ruta_archivo, sheet = 2)
head(dataset2,5); View(dataset2)
#---- Procesamiento de datos ----
##---- Determinemos la dimension de cada archivo ----
u1 <- dim(dataset1); print(u1)
u2 <- dim(dataset2); print(u2)
##---- Redifinicion de los nombres de las columnas ----
colnames(dataset1) <- c('GeneId','Descripcion',
                        'Ctrl0hr','Ctrl1hr',
                        'Encyst1hr','Ctrl4hr',
                        'Encyst4hr','Ctrl8hr',
                        'Encyst8hr');

print(dataset1); View(dataset1)
colnames(dataset2) <- c('GeneId','Descripcion',
                        'Ctrl0hr','Encyst1hr',
                        'Encyst4hr','Encyst8hr');
print(dataset2); View(dataset2)
##---- Redifinicion de los nombres de los renglones ----
###---- Primero el dataset1 ----
dataset1 <- dataset1[!duplicated(dataset1), ]; View(dataset1)
rownames(dataset1) <- dataset1$GeneId; View(dataset1)
write.csv(dataset1,"Enquistamiento.csv")
#NewGenes <- read_csv("Enquistamiento.csv")
#View(NewGenes)
###---- Ahora el dataset2 ----
dataset2 <- dataset2[!duplicated(dataset2), ]; View(dataset2)
rownames(dataset2) <- dataset2$GeneId; View(dataset2)
write.csv(dataset2,"Enquistamientotoplot.csv")
#NewGenes <- read_csv("Enquistamientotoplot.csv")
#View(NewGenes)
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# Reestablecimiento a su valor original, restando 2 a cadsa entrada ----
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# ---- Preprocesamiento de los datos ----
indices <- c(3,5,7,9)
NewGenes <- dataset1[,indices]; View(NewGenes)
rownames(NewGenes) <- dataset1$GeneId; View(NewGenes) 
NewGenes_Clean <- NewGenes
# ---- Generacion de las bdd modificadas ---- 
NewGenes_log2   <-  log2(NewGenes_Clean + 1)
NewGenesOrdered <- NewGenes_Clean[order(-NewGenes_Clean$Ctrl0hr),]
top_genes_NewGenes_Trop <- NewGenesOrdered %>% filter(Ctrl0hr>0) %>% select(Ctrl0hr)
k  <- dim(NewGenes_Clean); proporcion <- .75; 
NS <- round(k[1]*proporcion);
random_genes_NewGenes  <- sample(rownames(NewGenes_Clean),NS );
sampledNewGenes_Log2   <- NewGenes_log2[random_genes_NewGenes, ];
NewGenes_log2_filtrado <- NewGenes_log2[rowSums(NewGenes_log2) != 0, ]; View(NewGenes_log2_filtrado)
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

