#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# Inicializacion del directorio de trabajo ----
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#setwd("~/Desktop/MiGithub/MapaCalorReloaded/DraElisaNewData") # Mac
setwd("~/Documents/GitHub/MapaCalorReloaded/DraElisaNewData") # Laptop
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# ---- Carga de librerias ----
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
library(RColorBrewer)
library(readr)
library(pheatmap)
library(dplyr)
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
# ---- Genes de Estress-Suero ----
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
NewGenes <- read_csv("ExpresionMybStresSuero.csv"); View(NewGenes)
colnames(NewGenes) <- c("GeneID","Normal","Serum_starved","Serum_Replenished"); View(NewGenes)
rownames(NewGenes) <- NewGenes$GeneID
NewGenes_Clean <- NewGenes[,2:4]; View(NewGenes_Clean)
NewGenes_Clean <- as.data.frame(NewGenes_Clean)
rownames(NewGenes_Clean) <- NewGenes$GeneID; View(NewGenes_Clean)
# Generacion de las bdd modificadas
NewGenes_log2   <-  log2(NewGenes_Clean + 1)
NewGenesOrdered <- NewGenes_Clean[order(-NewGenes_Clean$Normal),]
top_genes_NewGenes_Trop <- NewGenesOrdered %>% filter(Normal>0) %>% select(Normal)
k  <- dim(NewGenes_Clean); proporcion <- 1; 
NS <- round(k[1]*proporcion);
random_genes_NewGenes  <- sample(rownames(NewGenes_Clean),NS );
sampledNewGenes_Log2   <- NewGenes_log2[random_genes_NewGenes, ];
NewGenes_log2_filtrado <- NewGenes_log2[rowSums(NewGenes_log2) != 0, ];View(NewGenes_log2_filtrado)
LetraSize <- 6
pdf("MapasCalor/BoxplotEstresSuero.pdf")
boxplot(NewGenes_log2_filtrado, las = 3)
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonEstresSuero.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Bajo Estrés y Suero")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonEstresSueroMycolors.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Bajo Estrés y Suero")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonEstresSueroMycolors2.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Bajo Estrés y Suero")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonEstresSueroMycolors3.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Bajo Estrés y Suero")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonEstresSueroMycolors4.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors4,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Bajo Estrés y Suero")
dev.off()
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# ---- Genes de Colon Virulento ----
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
NewGenes <- read_csv("ExpresionMybColonVirulento.csv"); View(NewGenes)
colnames(NewGenes) <- c("GeneID",
                        "Attenuated_Clean",
                        "Normal_Culture",
                        "Virulent_Colon",
                        "Virulent_Culture"); View(NewGenes)
rownames(NewGenes) <- NewGenes$GeneID
NewGenes_Clean <- NewGenes[,2:5]; View(NewGenes_Clean)
NewGenes_Clean <- as.data.frame(NewGenes_Clean)
rownames(NewGenes_Clean) <- NewGenes$GeneID; View(NewGenes_Clean)
#<< ==  ==  ==  ==  ==  ==  ==  ==  == ==  ==  ==  ==  ==  ==  ==  ==  ==  == >> 
# Generacion de las bdd modificadas
NewGenes_log2   <-  log2(NewGenes_Clean + 1); View(NewGenes_log2)
NewGenesOrdered <- NewGenes_Clean[order(-NewGenes_Clean$Normal_Culture),]
top_genes_NewGenes_Trop <- NewGenesOrdered %>% filter(Normal_Culture>0) %>% select(Normal_Culture)
k  <- dim(NewGenes_Clean); proporcion <- 1; 
NS <- round(k[1]*proporcion);
random_genes_NewGenes  <- sample(rownames(NewGenes_Clean),NS );
sampledNewGenes_Log2   <- NewGenes_log2[random_genes_NewGenes, ];
NewGenes_log2_filtrado <- NewGenes_log2[rowSums(NewGenes_log2) != 0, ];View(NewGenes_log2_filtrado)
LetraSize <- 6
pdf("MapasCalor/BoxplotColonVirulento.pdf")
boxplot(NewGenes_log2_filtrado, las = 3)
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonColonVirulento.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "ColonVirulento")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonColonVirulentoMycolors.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "ColonVirulento")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonColonVirulentoMycolors2.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "ColonVirulento")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonColonVirulentoMycolors3.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "ColonVirulento")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonColonVirulentoMycolors4.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors4,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "ColonVirulento")
dev.off()
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# ---- Genes Estres-Calor ----
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
NewGenes <- read_csv("ExpresionMybStressCalor.csv"); View(NewGenes)
colnames(NewGenes)
## ---- Reetiquetado de nombres de columnas y se extraen las columnas a trabajar  ---
colnames(NewGenes) <- c("GeneID",
                        "HM1IMSS_Rahman_RNA",
                        "Sense_0hr",
                        "Sense_2hr",
                        "Sense_4hr",
                        "Sense_8hr","Empty"); View(NewGenes)
rownames(NewGenes) <- NewGenes$GeneID
NewGenes_Clean <- NewGenes[,2:6]; View(NewGenes_Clean)
NewGenes_Clean <- as.data.frame(NewGenes_Clean)
rownames(NewGenes_Clean) <- NewGenes$GeneID; View(NewGenes_Clean)
## ---- Generacion de las bdd modificadas ---- 
NewGenes_log2   <-  log2(NewGenes_Clean + 1); View(NewGenes_log2)
#  ----  Orden con respecto a una de las columnas ---- 
NewGenesOrdered <- NewGenes_Clean[order(-NewGenes_Clean$HM1IMSS_Rahman_RNA),]
top_genes_NewGenes_Trop <- NewGenesOrdered %>% filter(HM1IMSS_Rahman_RNA>0) %>% select(HM1IMSS_Rahman_RNA)
#  ---- Se continua con el procesamiento de los datos ---- 
k  <- dim(NewGenes_Clean); proporcion <- 1; 
NS <- round(k[1]*proporcion);
random_genes_NewGenes  <- sample(rownames(NewGenes_Clean),NS );
sampledNewGenes_Log2   <- NewGenes_log2[random_genes_NewGenes, ];
NewGenes_log2_filtrado <- NewGenes_log2[rowSums(NewGenes_log2) != 0, ];View(NewGenes_log2_filtrado)
LetraSize <- 6
## ---- Generacion de graficos ----
pdf("MapasCalor/BoxplotEstresCalor.pdf")
boxplot(NewGenes_log2_filtrado, las = 3)
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonEstresCalor.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Estres-Calor")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonEstresCalorMycolors.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Estres-Calor")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonEstresCalorMycolors2.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Estres-Calor")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonEstresCalorMycolors3.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Estres-Calor")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonEstresCalorMycolors4.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors4,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Estres-Calor")
dev.off()
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# ---- Genes Enquistamiento ----
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# ---- Limpiamos el entorno y la consola ----
rm(list = ls())
cat("\014")
# ---- Lectura de los genes ----
NewGenes <- read_csv("ExpresionMybEnquistamiento.csv"); View(NewGenes)
colnames(NewGenes)
## ---- Reetiquetado de nombres de columnas y se extraen las columnas a trabajar  ---
colnames(NewGenes) <- c("GeneID",
                        "Trophs",
                        "Trophs_Low_Glucose",
                        "Cysts_1wk",
                        "Cysts_3wk",
                        "Cysts_8wk"); View(NewGenes)
rownames(NewGenes) <- NewGenes$GeneID
NewGenes_Clean <- NewGenes[,2:6]; View(NewGenes_Clean)
NewGenes_Clean <- as.data.frame(NewGenes_Clean)
rownames(NewGenes_Clean) <- NewGenes$GeneID; View(NewGenes_Clean)
#NewGenesClean <- na.omit(NewGenes_Clean); View(NewGenesClean)
NewGenesClean <- NewGenes_Clean[1:32,]; View(NewGenesClean)
NewGenes_Clean <- NewGenesClean; View(NewGenes_Clean)
## ---- Generacion de las bdd modificadas ---- 
NewGenes_Clean <- as.data.frame(lapply(NewGenes_Clean, as.numeric))
NewGenes_log2   <-  log2(NewGenes_Clean + 1); View(NewGenes_log2)
#  ----  Orden con respecto a una de las columnas ---- 
NewGenesOrdered <- NewGenes_Clean[order(-NewGenes_Clean$Trophs),]
top_genes_NewGenes_Trop <- NewGenesOrdered %>% filter(Trophs>0) %>% select(Trophs)
#  ---- Se continua con el procesamiento de los datos ---- 
k  <- dim(NewGenes_Clean); proporcion <- 1; 
NS <- round(k[1]*proporcion);
random_genes_NewGenes  <- sample(rownames(NewGenes_Clean),NS );
sampledNewGenes_Log2   <- NewGenes_log2[random_genes_NewGenes, ];
NewGenes_log2_filtrado <- NewGenes_log2[rowSums(NewGenes_log2) != 0, ];View(NewGenes_log2_filtrado)
LetraSize <- 6
## ---- Generacion de graficos ----
pdf("MapasCalor/Enquistamiento.pdf")
boxplot(NewGenes_log2_filtrado, las = 3)
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonEnquistamiento.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Enquistamiento")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonEnquistamientoMycolors.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Enquistamiento")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonEnquistamientoMycolors2.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Enquistamiento")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonEnquistamientoMycolors3.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Enquistamiento")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonEnquistamientoMycolors4.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors4,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Enquistamiento")
dev.off()
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# ---- Genes Modelo Intestinal ----
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# ---- Limpiamos el entorno y la consola ----
rm(list = ls())
cat("\014")
## ---- Personalización de la paleta de colores ----
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
# ---- Lectura de los genes ----
NewGenes <- read_csv("ExpresionMybModIntestinal.csv"); View(NewGenes)
colnames(NewGenes)
## ---- Reetiquetado de nombres de columnas y se extraen las columnas a trabajar  ---
colnames(NewGenes) <- c("GeneID",
                        "HM1_Trophs_TYI",
                        "Rahman",
                        "Mouse_Adptd",
                        "HTrophs1d",
                        "Trophs29d"); View(NewGenes)
rownames(NewGenes) <- NewGenes$GeneID
NewGenes_Clean <- NewGenes[,2:6]; View(NewGenes_Clean)
NewGenes_Clean <- as.data.frame(NewGenes_Clean)
rownames(NewGenes_Clean) <- NewGenes$GeneID; View(NewGenes_Clean)
#NewGenesClean <- na.omit(NewGenes_Clean); View(NewGenesClean)
#NewGenesClean <- NewGenes_Clean[1:32,]; View(NewGenesClean)
#NewGenes_Clean <- NewGenesClean; View(NewGenes_Clean)
## ---- Generacion de las bdd modificadas ---- 
NewGenes_Clean <- as.data.frame(lapply(NewGenes_Clean, as.numeric))
NewGenes_log2   <-  log2(NewGenes_Clean + 1); View(NewGenes_log2)
#  ----  Orden con respecto a una de las columnas ---- 
NewGenesOrdered <- NewGenes_Clean[order(-NewGenes_Clean$HM1_Trophs_TYI),]
top_genes_NewGenes_Trop <- NewGenesOrdered %>% filter(HM1_Trophs_TYI>0) %>% select(HM1_Trophs_TYI)
#  ---- Se continua con el procesamiento de los datos ---- 
k  <- dim(NewGenes_Clean); proporcion <- 1; 
NS <- round(k[1]*proporcion);
random_genes_NewGenes  <- sample(rownames(NewGenes_Clean),NS );
sampledNewGenes_Log2   <- NewGenes_log2[random_genes_NewGenes, ];
NewGenes_log2_filtrado <- NewGenes_log2[rowSums(NewGenes_log2) != 0, ];View(NewGenes_log2_filtrado)
LetraSize <- 6
                                        
## ---- Generacion de graficos ----
pdf("MapasCalor/ModIntestinal.pdf")
boxplot(NewGenes_log2_filtrado, las = 3)
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonModIntestinal.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Modelo-Intestinal")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonModIntestinalMycolors.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Modelo-Intestinal")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonModIntestinalMycolors2.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Modelo-Intestinal")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonModIntestinalMycolors3.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Modelo-Intestinal")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonModIntestinalMycolors4.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors4,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Modelo-Intestinal")
dev.off()
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
