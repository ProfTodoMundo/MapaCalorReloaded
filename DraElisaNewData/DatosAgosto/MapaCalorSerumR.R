#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#setwd("~/Desktop/MiGithub/MapaCalorReloaded") # computadora de la casa
#setwd("~/Desktop/MiGithub/MapaCalorReloaded/DraElisaNewData")
setwd("~/Documents/GitHub/MapaCalorReloaded/DraElisaNewData/DatosAgosto")
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
library(RColorBrewer)
library(readr)
library(pheatmap)
library(dplyr)
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
#NewGenes <- read_csv("ExpresiónMybStresSuero.csv")
NewGenes <- read_csv("Transcriptoma_estres_suero.csv")
colnames(NewGenes)
colnames(NewGenes) <- c("GeneID","DescripcionProd",
                        "Normal","Serum_starved","Serum_Replenished")
View(NewGenes)
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
MisMaterias <- NewGenes$DescripcionProd; 
n <- length(MisMaterias); TodasMaterias <- as.character(MisMaterias)
MateriasCorrected <- c(); MateriasTemporal <- TodasMaterias
for(i in 1:n){
  MateriasCorrected <- toupper(TodasMaterias[i])
  MateriasCorrected <- gsub("á","a",as.character(MateriasCorrected,fixed=TRUE));
  MateriasCorrected <- gsub("é","e",as.character(MateriasCorrected,fixed=TRUE));
  MateriasCorrected <- gsub("í","i",as.character(MateriasCorrected,fixed=TRUE));
  MateriasCorrected <- gsub("ó","o",as.character(MateriasCorrected,fixed=TRUE));
  MateriasCorrected <- gsub("ú","u",as.character(MateriasCorrected,fixed=TRUE));
  MateriasCorrected <- gsub("Á","A",as.character(MateriasCorrected,fixed=TRUE));
  MateriasCorrected <- gsub("É","E",as.character(MateriasCorrected,fixed=TRUE));
  MateriasCorrected <- gsub("Í","I",as.character(MateriasCorrected,fixed=TRUE));
  MateriasCorrected <- gsub("Ó","O",as.character(MateriasCorrected,fixed=TRUE));
  MateriasCorrected <- gsub("Ú","U",as.character(MateriasCorrected,fixed=TRUE));
  MateriasCorrected <- gsub("Ü","U",as.character(MateriasCorrected,fixed=TRUE));
  MateriasCorrected <- gsub("Ö","OE",as.character(MateriasCorrected,fixed=TRUE));
  MateriasCorrected <- gsub("Ñ","n",as.character(MateriasCorrected,fixed=TRUE));
  MateriasCorrected <- gsub("'","",as.character(MateriasCorrected,fixed=TRUE));
  MateriasCorrected <- gsub(",","",as.character(MateriasCorrected,fixed=TRUE));
  MateriasCorrected <- gsub("/","",as.character(MateriasCorrected,fixed=TRUE));
  TodasMaterias[i] <- MateriasCorrected
}
NewGenes$DescripcionProd<- TodasMaterias;
summary(NewGenes$DescripcionProd)
TemasSelectos <- c();TemasSelectos <- NewGenes$DescripcionProd; 
n<- length(TemasSelectos)
todasmismaterias <- sort(TemasSelectos)
TemasSelectos <- todasmismaterias;
materia <- TemasSelectos[1]; n <- length(TemasSelectos); k <- 2;
ListaTemasSelectos <- c(); ListaTemasSelectos[1] <- materia; Materia2 <- c()
for(i in 2:n){
  Materia2 <- TemasSelectos[i]
  if(materia==Materia2){materia <- Materia2;
  }else{
    ListaTemasSelectos[k]<- Materia2;materia <- Materia2; k <- k+1}
}
typeof(ListaTemasSelectos); length(ListaTemasSelectos)
temp <- factor(ListaTemasSelectos); niveles <- levels(temp)
NewGenes$DescripcionProd <- factor(NewGenes$DescripcionProd,
                          levels = niveles,
                          labels = ListaTemasSelectos)
summary(NewGenes$DescripcionProd);
rownames(NewGenes) <- NewGenes$GeneID
colnames(NewGenes)
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
NewGenes_Clean <- NewGenes[,3:5]; View(NewGenes_Clean)
NewGenes_Clean <- as.data.frame(NewGenes_Clean)
rownames(NewGenes_Clean) <- NewGenes$GeneID; View(NewGenes_Clean)
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
# Generacion de las bdd modificadas
NewGenes_log2   <-  log2(NewGenes_Clean + 1)
NewGenesOrdered <- NewGenes_Clean[order(-NewGenes_Clean$Normal),]
top_genes_NewGenes_Trop <- NewGenesOrdered %>% filter(Normal>0) %>% select(Normal)
k  <- dim(NewGenes_Clean); proporcion <- 1; 
NS <- round(k[1]*proporcion);
random_genes_NewGenes  <- sample(rownames(NewGenes_Clean),NS );
sampledNewGenes_Log2   <- NewGenes_log2[random_genes_NewGenes, ];
NewGenes_log2_filtrado <- NewGenes_log2[rowSums(NewGenes_log2) != 0, ];View(NewGenes_log2_filtrado)
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
boxplot(NewGenes_Clean, las = 3)
boxplot(NewGenes_log2, las = 3)
boxplot(NewGenes_log2_filtrado, las = 3)
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
pdf("MapasCalor/HeathmapCFNormRenglonSinBordes.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,
         scale = "row",fontsize_row = LetraSize, border_color = NA)
dev.off()
#
pdf("MapasCalor/HeathmapCFNormRenglonSinBordesMycolors.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors,
         scale = "row",fontsize_row = LetraSize, border_color = NA)
dev.off()
#
pdf("MapasCalor/HeathmapCFNormRenglonSinBordesMycolors2.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors2,
         scale = "row",fontsize_row = LetraSize, border_color = NA)
dev.off()
#
pdf("MapasCalor/HeathmapCFNormRenglonSinBordesMycolors3.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors3,
         scale = "row",fontsize_row = LetraSize, border_color = NA)
dev.off()
#
pdf("MapasCalor/HeathmapCFNormRenglonSinBordesMycolors4.pdf") # sin dendogramas columnas
pheatmap(NewGenes_log2_filtrado, cluster_cols = FALSE,color = my_colors4,
         scale = "row",fontsize_row = LetraSize, border_color = NA)
dev.off()
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
colores <- c("red", "blue","yellow", "green", "grey","orange", "purple")

pdf("MapasCalor/BoxPlotNewGenes.pdf")
boxplot(NewGenes_Clean, las = 3, col = colores)
legend("topleft", legend = c("Normal","Serum_Starved","Serum_Replenised"),
       fill = colores, title = "Transcriptoma Estres Suero")
dev.off()

pdf("MapasCalor/BoxPlotLog2NewGenes.pdf")
boxplot(NewGenes_log2, las = 2, col = colores, ylab="Log2(EiMybs)",
        names = c("Normal","Serum_Starved","Serum_Replenised"),
        cex.axis = 0.8)
legend("topleft", legend = c("Normal","Serum_Starved","Serum_Replenised"),
       fill = colores, title = "Transcriptoma Estres Suero")
dev.off()

pdf("MapasCalor/BoxPlotFilteredLog2NewGenes.pdf")
boxplot(NewGenes_log2_filtrado, las = 3, col = colores, ylab="Log2(EiMybs)",
        names = c("Normal","Serum_Starved","Serum_Replenised"),
        cex.axis = 0.8)
legend("topleft", legend = c("Normal","Serum_Starved","Serum_Replenised"),
       fill = colores, title = "Transcriptoma Estres Suero")
dev.off()