#----Instalacion y carga de librerias a utilizar ----
# Instala y carga la librería readxl si aún no lo has hecho
# install.packages("readxl")
library(readxl)
#---- Definicion del directorio de trabajo ----
setwd("~/Documents/GitHub/MapaCalorReloaded/ClaseCGenomicas")
#setwd("~/Desktop/MiGithub/MapaCalorReloaded/ClaseCGenomicas")
##---- Limpieza del ambiente de RStudio ----
rm(list = ls())
##----  Especifica la ruta de tu archivo Excel----
ruta_archivo <- "Histograma.xlsx"
##---- Lectura de las hojas de excel incluidas en el archivo----
dataset1 <- read_excel(ruta_archivo, sheet = 1)
head(dataset1,5)
dataset2 <- read_excel(ruta_archivo, sheet = 2)
head(dataset2,5)
#---- Procesamiento de datos ----
##---- Determinemos la dimension de cada archivo ----
u1 <- dim(dataset1); print(u1)
u2 <- dim(dataset2); print(u2)
##---- Obtenemos solo las columnas no vacias ----
indices1 <- c(1,3,6); dataset1 <- dataset1[,indices1]; print(dataset1)
indices2 <- c(1,3,4); dataset2 <- dataset2[,indices2]; print(dataset2)
##---- Redifinicion de los nombres de las columnas ----
colnames(dataset1) <- c('GeneId','PDescription','TPM_Basal'); print(dataset1)
colnames(dataset2) <- c('GeneId','PDescription','TPM_Basal'); print(dataset2)
##---- Redifinicion de los nombres de los renglones ----
###---- Primero el dataset1 ----
dataset1 <- dataset1[!duplicated(dataset1), ]; View(dataset1)
rownames(dataset1) <- dataset1$GeneId; View(dataset1)
###---- Ahora el dataset2 ----
dataset2 <- dataset2[!duplicated(dataset2), ]; View(dataset2)
rownames(dataset2) <- dataset2$GeneId; View(dataset2)
#----Grafica de datos ----
## recordemos como se ven los datasets
View(dataset1);View(dataset2);
summary(dataset1); summary(dataset2)
##---- Cargemos las librerias para graficar ----
#install.packages('ggplot2')
##---- Graficando con la libreria precargada ----
nbreaks <- 10
tBE <- hist(dataset1$TPM_Basal, breaks = nbreaks, col= rainbow(1,0.7), main = 'BasalExpresion')
BE <- dataset1$TPM_Basal
Log2BE <- log2(BE)
nBE    <- length(Log2BE)
hist(Log2BE, breaks = nbreaks, col= rainbow(25,0.3), 
     main = ' Log2 Basal Expresion')
dataset1$Log2BE <- Log2BE
##---- Graficando con ggplot2 ----
###---- En caso de ser necesario instalamos las librerias que se requieren ----
if (!requireNamespace("ggplot2", quietly = TRUE)) {install.packages("ggplot2")}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {install.packages("RColorBrewer")}
library(ggplot2)
library(RColorBrewer)
###---- Template 1, dodge, jitter, identity, stack ----
numbins = 60;
anchobins = 0.35;
ggplot(dataset1, aes(x = dataset1$Log2BE)) +
  geom_histogram(bins=numbins,
                 binwidth = anchobins, 
                 fill = "green", 
                 color = "red",
                 position = "dodge") +
  theme_classic()+ 
  theme(legend.position="left")+  
  labs(title = "TPM Basal Expresion", 
       x = "Valores", 
       y = "Frecuencia")
###---- Template 2, jitter + position legend:none, top, bottom, left, right, topright, ----
# topleft, bottomtright, bottomleft
ggplot(dataset1, aes(x = dataset1$Log2BE)) +
  geom_histogram(bins=numbins,
                 binwidth = anchobins, 
                 fill = "green", 
                 color = "red",
                 position = "jitter") +
  theme_classic()+ 
  theme(legend.position="right")+  
  labs(title = "TPM Basal Expresion", 
       x = "Valores", 
       y = "Frecuencia")
###---- Template 3, dodge + uso de la paleta de colores (no corre) ----
ggplot(dataset1, aes(x = Log2BE)) +
  geom_histogram(bins=numbins,
                 binwidth = anchobins, 
                 #fill = "green", 
                 color = "gray",
                 position = "dodge") +
  theme_classic() + 
  theme(legend.position="bottomright") +  
  scale_fill_brewer(7,palette="GnBu")+
  labs(title = "TPM Basal Expresion", 
       x = "Valores", 
       y = "Frecuencia")
# - - - - - - - -  - - - - - - - -  - - - - - - - -  - - - - - - - -  - - - - - - - -  - - - - - - - - 
ggplot(dataset1, aes(x = Log2BE)) +
  geom_histogram(bins=numbins,
                 binwidth = anchobins, 
                 fill = "green", 
                 color = "gray",
                 position = "dodge") +
  theme_classic() + 
  theme(legend.position="bottomright") +  
  scale_fill_brewer(7,palette="YlOrRd")+
  labs(title = "TPM Basal Expresion", 
       x = "Valores", 
       y = "Frecuencia")
###---- Template 4, cambio de tonalidadades ----
ggplot(dataset1, aes(x = Log2BE)) +
  geom_histogram(bins=numbins,
                 binwidth = anchobins, 
                 alpha = 0.25,
                 fill = "green", 
                 color = "gray",
                 position = "dodge") +
  theme_classic() + 
  theme(legend.position="bottomright") +  
  labs(title = "TPM Basal Expresion", 
       x = "Valores", 
       y = "Frecuencia")
###---- Template 5, tipo de lineas ----
ggplot(dataset1, aes(x = Log2BE)) +
  geom_histogram(bins=numbins,
                 binwidth = anchobins, 
                 alpha = 0.95,
                 fill = "yellow", 
                 color = "green",
                 position = "dodge",
                 linetype = "dashed") +
  theme_classic() + 
  theme(legend.position="bottomright") +  
  labs(title = "TPM Basal Expresion", 
       x = "Valores", 
       y = "Frecuencia")