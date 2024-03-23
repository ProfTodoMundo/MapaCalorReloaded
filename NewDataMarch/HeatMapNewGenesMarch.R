# ---- Carga de librerias ----
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
library(readxl)
library(RColorBrewer)
library(readr)
library(pheatmap)
library(dplyr)
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#---- Lectura de los datos ----
setwd("~/Documents/GitHub/MapaCalorReloaded/NewDataMarch")
# Ruta del archivo Excel
archivo_excel <- "~/Documents/GitHub/MapaCalorReloaded/NewDataMarch/EhMybSgenesexpression.xls"
# Obtener los nombres de las hojas del libro Excel
nombres_hojas <- excel_sheets(archivo_excel)
# Lista para almacenar los datos de cada hoja
datos_por_hoja <- list()
# Iterar sobre los nombres de las hojas y leer cada hoja
for (nombre_hoja in nombres_hojas) {
  datos_por_hoja[[nombre_hoja]] <- read_excel(archivo_excel, sheet = nombre_hoja)
}
# Exportar cada hoja a archivos CSV
for (nombre_hoja in nombres_hojas) {
  nombre_archivo_csv <- paste0(nombre_hoja, ".csv")
  write.csv(datos_por_hoja[[nombre_hoja]], file = nombre_archivo_csv, row.names = FALSE)
  cat("Archivo CSV '", nombre_archivo_csv, "' exportado con éxito.\n")
}
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# ---- Personalización de la paleta de colores ----
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
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
# ---- Genes de MicroArray ----
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
NewGenes <- read_csv("Microarray.csv"); View(NewGenes)
u <- dim(NewGenes); n <- u[1];
NewGenes <- NewGenes[2:n,]
colnames(NewGenes)
colnames(NewGenes) <- c("GeneID","Nombre",
                        "HM1IMSSTrophs","RahmanT","MouseAdapted",
                        "Trophs1dP","Trophs29d","HM1IMSSStgeConv",
                        "RahmStgeConv","200NHITrophsTYIStgConv",
                        "200NHITrophsTYILowGlucose","MS753544Ttophscyst1wkrobns",
                        "MS753544Ttophscyst8wkrobns","2592100Ttophscyst3wkrobns"); 
View(NewGenes)
u <- dim(NewGenes)
rownames(NewGenes) <- NewGenes$GeneID
NewGenes_Clean <- NewGenes[,3:u[2]]; View(NewGenes_Clean)
NewGenes_Clean <- as.data.frame(NewGenes_Clean); View(NewGenes_Clean)
rownames(NewGenes_Clean) <- NewGenes$GeneID; View(NewGenes_Clean)

##---- Generacion de las bdd modificadas ----
NewGenes_log2   <-  log2(NewGenes_Clean + 1)
NewGenesOrdered <- NewGenes_Clean[order(-NewGenes_Clean$HM1IMSSTrophs),]; View(NewGenesOrdered)
top_genes_NewGenes_Trop <- NewGenesOrdered %>% filter(HM1IMSSTrophs>0) %>% select(HM1IMSSTrophs)
View(top_genes_NewGenes_Trop)
k  <- dim(NewGenes_Clean); proporcion <- 1; 
NS <- round(k[1]*proporcion);
random_genes_NewGenes  <- sample(rownames(NewGenes_Clean),NS );
sampledNewGenes_Log2   <- NewGenes_log2[random_genes_NewGenes, ];
NewGenes_log2_filtrado <- NewGenes_log2[rowSums(NewGenes_log2) != 0, ];
View(NewGenes_log2_filtrado)
LetraSize <- 6
pdf("MapasCalor/BoxplotEhMybSMicroArray.pdf")
boxplot(NewGenes_log2_filtrado, las = 2)
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonEhMybSMicroArray.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Genes EhMybS - MicroArray")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonEhMybSMicroArrayMycolors.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Genes EhMybS - MicroArray")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonEhMybSMicroArrayMycolors2.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Genes EhMybS - MicroArray")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonEhMybSMicroArrayMycolors3.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Genes EhMybS - MicroArray")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonEhMybSMicroArrayMycolors4.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors4,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Genes EhMybS - MicroArray")
dev.off()
## ---- Creación de gráficas ----
## ---- Parametros 60 y 3.75 ----
cell_width <- 30
mi_font_size <- 3.75
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width,main ="Amplio de celda 30 y letras de genes 3.25")
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,main ="Amplio de celda 30 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"conBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 30 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 30 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 30 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 30 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
## ---- Parametros 65 y 4.25 ----
cell_width <- 25
mi_font_size <- 4.25
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,scale = "row",
         fontsize_row = mi_font_size, cellwidth = cell_width,main ="Amplio de celda 25 y letras de genes 4.25")
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,main ="Amplio de celda 25 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores1",".pdf")
pdf(nombregrafico)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 25 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 25 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 25 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 25 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
## ---- Parametros 70 y 3.75 ----
cell_width <- 20
mi_font_size <- 3.75
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,scale = "row",
         fontsize_row = mi_font_size, cellwidth = cell_width,main ="Amplio de celda 20 y letras de genes 3.75")
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,main ="Amplio de celda 20 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 20 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 20 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 20 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 20 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
## ---- Parametros 75 y 3.75 ----
cell_width <- 15
mi_font_size <- 3.75
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,scale = "row",
         fontsize_row = mi_font_size, cellwidth = cell_width,main ="Amplio de celda 15 y letras de genes 3.75")
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,main ="Amplio de celda 15 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 15 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 15 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 15 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 15 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
## ---- Parametros 60 y 4.25 ----
cell_width <- 20
mi_font_size <- 4.25
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,scale = "row",
         fontsize_row = mi_font_size, cellwidth = cell_width,main ="Amplio de celda 20 y letras de genes 4.25")
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,main ="Amplio de celda 20 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 20 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 20 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 20 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapMICROARRAY","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 20 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/BoxplotMICROARRAY",".pdf")
pdf(nombregrafico)
boxplot(NewGenes_log2_filtrado, las = 3)
dev.off()
# ---- Genes RNASeq ----
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
NewGenes <- read_csv("RNAseq.csv"); View(NewGenes)
colnames(NewGenes)
colnames(NewGenes) <- c("GeneID","Nombre",
                        "TrophozoiteBasal","Normal","AttenuatedCulture",
                        "VirulentColon","VirulentCulture","CtrlStarvation",
                        "SerumStarvated","SerumReplenished","H5_0hr",
                        "H5_2hr","H5_4hr","H5_8hr"); 
View(NewGenes)
u <- dim(NewGenes)
rownames(NewGenes) <- NewGenes$GeneID
NewGenes_Clean <- NewGenes[,3:u[2]]; View(NewGenes_Clean)
NewGenes_Clean <- as.data.frame(NewGenes_Clean); View(NewGenes_Clean)
rownames(NewGenes_Clean) <- NewGenes$GeneID; View(NewGenes_Clean)
##---- Generacion de las bdd modificadas ----
NewGenes_log2   <-  log2(NewGenes_Clean + 1)
NewGenesOrdered <- NewGenes_Clean[order(-NewGenes_Clean$Normal),]; View(NewGenesOrdered)
top_genes_NewGenes <- NewGenesOrdered %>% filter(Normal>0) %>% select(Normal)
View(top_genes_NewGenes)
k  <- dim(NewGenes_Clean); proporcion <- 1; 
NS <- round(k[1]*proporcion);
random_genes_NewGenes  <- sample(rownames(NewGenes_Clean),NS );
sampledNewGenes_Log2   <- NewGenes_log2[random_genes_NewGenes, ];
NewGenes_log2_filtrado <- NewGenes_log2[rowSums(NewGenes_log2) != 0, ];
View(NewGenes_log2_filtrado)
LetraSize <- 6
pdf("MapasCalor/BoxplotRNAseq.pdf")
boxplot(NewGenes_log2_filtrado, las = 2)
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonRNAseq.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Genes EhMybS - RNAseq")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonRNAseqMycolors.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Genes EhMybS - RNAseq")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonRNAseqMycolors2.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Genes EhMybS - RNAseq")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonRNAseqMycolors3.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Genes EhMybS - RNAseq")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonRNAseqMycolors4.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors4,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Genes EhMybS - RNAseq")
dev.off()
## ---- Creación de gráficas ----
## ---- Parametros 60 y 3.75 ----
cell_width <- 30
mi_font_size <- 3.75
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width,main ="Amplio de celda 30 y letras de genes 3.25")
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,main ="Amplio de celda 30 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"conBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 30 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 30 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 30 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 30 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
## ---- Parametros 65 y 4.25 ----
cell_width <- 25
mi_font_size <- 4.25
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,scale = "row",
         fontsize_row = mi_font_size, cellwidth = cell_width,main ="Amplio de celda 25 y letras de genes 4.25")
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,main ="Amplio de celda 25 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores1",".pdf")
pdf(nombregrafico)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 25 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 25 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 25 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 25 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
## ---- Parametros 70 y 3.75 ----
cell_width <- 20
mi_font_size <- 3.75
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,scale = "row",
         fontsize_row = mi_font_size, cellwidth = cell_width,main ="Amplio de celda 20 y letras de genes 3.75")
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,main ="Amplio de celda 20 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 20 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 20 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 20 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 20 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
## ---- Parametros 75 y 3.75 ----
cell_width <- 15
mi_font_size <- 3.75
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,scale = "row",
         fontsize_row = mi_font_size, cellwidth = cell_width,main ="Amplio de celda 15 y letras de genes 3.75")
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,main ="Amplio de celda 15 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 15 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 15 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 15 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 15 y letras de genes 3.75",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
## ---- Parametros 60 y 4.25 ----
cell_width <- 20
mi_font_size <- 4.25
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,scale = "row",
         fontsize_row = mi_font_size, cellwidth = cell_width,main ="Amplio de celda 20 y letras de genes 4.25")
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores1",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,main ="Amplio de celda 20 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 20 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores2",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,main ="Amplio de celda 20 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"ConBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 20 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width)
dev.off()
nombregrafico <-paste("MapasCalor/HeathmapRNAseq","Celda_",cell_width,"y FontSize_",mi_font_size,"SinBordes","Colores3",".pdf")
pdf(nombregrafico)
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,main ="Amplio de celda 20 y letras de genes 4.25",
         scale = "row",fontsize_row = mi_font_size, cellwidth = cell_width, border_color = NA)
dev.off()
nombregrafico <-paste("MapasCalor/BoxplotRNAseq",".pdf")
pdf(nombregrafico)
boxplot(NewGenes_log2_filtrado, las = 3)
dev.off()
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# ---- Genes Blanco S3 ----
GenesBlancoS3 <- read_excel("GenesBlancoS3.xls")
write_csv(GenesBlancoS3,"GenesBlancoS3.csv"); View(GenesBlancoS3)
NewGenes <- read_csv("GenesBlancoS3.csv"); View(NewGenes)
colnames(NewGenes)
colnames(NewGenes) <- c("GeneID",
                        "ProdDescription",
                        "Trophozoite",
                        "SNormalUnique",
                        "SerumStarved",
                        "SerumReplenished",
                        "NormalCulture",
                        "AttCulture",
                        "VirulentColon",
                        "VirulentCulture",
                        "Sense_0hr",
                        "Sense_2hr",
                        "Sense_4hr",
                        "Sense_8hr");
View(NewGenes)
rownames(NewGenes) <- NewGenes$GeneID
u <- dim(NewGenes)
NewGenes_Clean <- NewGenes[,3:u[2]]; View(NewGenes_Clean)
NewGenes_Clean <- as.data.frame(NewGenes_Clean)
rownames(NewGenes_Clean) <- NewGenes$GeneID; View(NewGenes_Clean)
#<< ==  ==  ==  ==  ==  ==  ==  ==  == ==  ==  ==  ==  ==  ==  ==  ==  ==  == >> 
# Generacion de las bdd modificadas
NewGenes_log2   <-  log2(NewGenes_Clean + 1); View(NewGenes_log2)
NewGenesOrdered <- NewGenes_Clean[order(-NewGenes_Clean$NormalCulture),]
top_genes_NewGenes_Trop <- NewGenesOrdered %>% filter(NormalCulture>0) %>% select(NormalCulture)
View(top_genes_NewGenes_Trop)
k  <- dim(NewGenes_Clean); proporcion <- 1; 
NS <- round(k[1]*proporcion);
random_genes_NewGenes  <- sample(rownames(NewGenes_Clean),NS );
sampledNewGenes_Log2   <- NewGenes_log2[random_genes_NewGenes, ];
NewGenes_log2_filtrado <- NewGenes_log2[rowSums(NewGenes_log2) != 0, ];View(NewGenes_log2_filtrado)
LetraSize <- 3.5
pdf("MapasCalor/BoxplotGenesBlancoS3.pdf")
boxplot(NewGenes_log2_filtrado, las = 3)
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonGenesBlancoS3.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Genes Blanco S3")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonGenesBlancoS3Mycolors.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Genes Blanco S3")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonGenesBlancoS3Mycolors2.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Genes Blanco S3")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonGenesBlancoS3Mycolors3.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Genes Blanco S3")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonGenesBlancoS3Mycolors4.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors4,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Genes Blanco S3")
dev.off()
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
