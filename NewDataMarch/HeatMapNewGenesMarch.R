#---- Lectura de los datos ----
setwd("~/Documents/GitHub/MapaCalorReloaded/NewDataMarch")
EhMybSgenesexpression <- read_excel("EhMybSgenesexpression.xls")
write_csv(EhMybSgenesexpression,"EhMybSGenes.csv")
View(EhMybSgenesexpression)
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# ---- Carga de librerias ----
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
library(readxl)
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
NewGenes <- read_csv("EhMybSGenes.csv"); View(NewGenes)
colnames(NewGenes)
colnames(NewGenes) <- c("GeneID","Nombre",
                        "TrophozoiteBasal","Normal",
                        "AttCulture","VirulentColon",
                        "VirulentCulture","CtrlStarvation",
                        "SerumStarved","SerumReplenished",
                        "HS_0hr","HS_2hr",
                        "HS_4hr","HS_8hr"); 
View(NewGenes)
u <- dim(NewGenes)
rownames(NewGenes) <- NewGenes$GeneID
NewGenes_Clean <- NewGenes[,3:u[2]]; View(NewGenes_Clean)
NewGenes_Clean <- as.data.frame(NewGenes_Clean); View(NewGenes_Clean)
rownames(NewGenes_Clean) <- NewGenes$GeneID; View(NewGenes_Clean)
#---- Generacion de las bdd modificadas ----
NewGenes_log2   <-  log2(NewGenes_Clean + 1)
NewGenesOrdered <- NewGenes_Clean[order(-NewGenes_Clean$Normal),]; View(NewGenesOrdered)
top_genes_NewGenes_Trop <- NewGenesOrdered %>% filter(Normal>0) %>% select(Normal)
View(top_genes_NewGenes_Trop)
k  <- dim(NewGenes_Clean); proporcion <- 1; 
NS <- round(k[1]*proporcion);
random_genes_NewGenes  <- sample(rownames(NewGenes_Clean),NS );
sampledNewGenes_Log2   <- NewGenes_log2[random_genes_NewGenes, ];
NewGenes_log2_filtrado <- NewGenes_log2[rowSums(NewGenes_log2) != 0, ];
View(NewGenes_log2_filtrado)
LetraSize <- 6
pdf("MapasCalor/BoxplotEhMybS.pdf")
boxplot(NewGenes_log2_filtrado, las = 2)
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonEhMybS.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Genes EhMybS")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonEhMybSMycolors.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Genes EhMybS")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonEhMybSMycolors2.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Genes EhMybS")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonEhMybSMycolors3.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Genes EhMybS")
dev.off()
pdf("MapasCalor/HeathmapCFNormRenglonEhMybSMycolors4.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors4,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Genes EhMybS")
dev.off()
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# ---- Genes de Colon Virulento ----
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
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
LetraSize <- 6
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
pdf("MapasCalor/HeathmapCFNormRenglonGenesBlancoS3oMycolors4.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors4,
         scale = "row",fontsize_row = LetraSize, border_color = NA,
         main = "Genes Blanco S3")
dev.off()
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
