
#===========================================================================
setwd("~/Desktop/Datos Elisa")
#setwd("~/Escritorio/Datos Elisa")
#===========================================================================
library(ggplot2);   library(dplyr);         library(readxl);
library(pastecs);   library(sciplot);       library(MASS);
library(gridExtra); library("gplots");      library("lattice");
#library(car);       library(gridExtra);     library(lattice);
library(corrplot);  library(readr);         library(readxl);   
#library(rvest);    library(RSQLite); 
library(DBI);    
#library(xml2);    
#library(RCurl);     
#library(devtools);
library(ggplot2);   library(datasets);      library(dplyr);
library(sciplot);   library(scatterplot3d); #library("car")
library(psych);     library("gplots");      library("plotrix")
library(gplots);    library(moments);       #library(univariateML)
#===========================================================================
library("fitdistrplus"); library("MASS"); library("survival")
#===========================================================================
setwd("~/Desktop/Datos Elisa")
url1 <- 'GenesERMDef.xlsx'
dataset1 <- read_excel(url1)
dataset <- na.omit(dataset1)
dataset$Secuencia = factor(dataset$Secuencia,
                           levels = c('CAACGG','CAACTG','TAACGG','TAACTG'),
                           labels = c('CAACGG','CAACTG','TAACGG','TAACTG'))
colnames(dataset) <- c( 'GenId1', 'ProdDescription','ExpBasal', 'GeneId2',
                        'Sequence','Position','DigitVer')
datosn <- dataset$Position;d <- as.numeric(datosn);dataset$Position <- d
dataset$DigitVerif <- dataset$Sequence
dataset$DigitVerif = factor(dataset$DigitVerif,
                            levels = c('CAACGG','CAACTG','TAACGG','TAACTG'),
                            labels = c(1,2,3,4))
CtoDatos <- dataset[,c('Sequence','DigitVerif','Position','ExpBasal','GenId1')]
colnames(CtoDatos) <- c('Sequence','CodedSeq','Position','BasalExp','GenId')
bdd <- CtoDatos
summary(bdd)
nbreaks <- 10
tBE <- hist(bdd$BasalExp, breaks = nbreaks, col= rainbow(1,0.7), main = 'BasalExpresion')
BE <- bdd$BasalExp
#par(mfrow=c(2,1))
Log2BE <- log2(BE)
nBE    <- length(Log2BE)
hist(Log2BE, breaks = nbreaks, col= rainbow(25,0.3), 
     main = ' Log2 Basal Expresion')
meanL2BE <- mean(Log2BE)
StdDevL2BE <- sd(Log2BE)
NormLog2BE <- (Log2BE-meanL2BE)/StdDevL2BE
tst<- NormLog2BE
hist(tst, breaks = nbreaks, col= 1:5, 
     main = 'Normalized Log2 Basal Expresion',
     xlab='Basal Expresion',
     ylab= 'Frequency Basal Expresion')

fw1<-fitdist(tst, "norm")
plotdist(tst, histo = TRUE, demp = TRUE)
nnorm.f <- fitdist(tst,"norm")
summary(nnorm.f)
#par(mfrow=c(2,2))
denscomp(nnorm.f,legendtext = 'Dist Normal')
qqcomp(nnorm.f,legendtext = 'Dist Normal')
cdfcomp(nnorm.f,legendtext = 'Dist Normal')
ppcomp(nnorm.f,legendtext = 'Dist Normal')

probs <- 0;
probs[8] = 0.175;  probs[9] = 0.825; 
probs[7] = 0.15;   probs[10] = 0.85;   
probs[6] = 0.125;  probs[11] = 0.875; 
probs[5] = 0.1;    probs[12] = 0.9;    
probs[4] = 0.075;  probs[13] = 0.925; 
probs[3] = 0.05;   probs[14] = 0.95;   
probs[2] = 0.025;  probs[15] = 0.975; 
probs[1] = 0.005;  probs[16] = 0.995;  
CuantilesData <- quantile(tst,prob = probs)
CuantilesModel <- qnorm(probs, mean=0, sd=1)
Cuantilillos <- t(CuantilesModel)
colnames(Cuantilillos) <- c('0.5%','2.5%','5%','7.5%',
                            '10%','12.5%','15%','17.5%',
                            '82.5%','85%','87.5%','90%',
                            '92.5%','95%','97.5%','99.5%')
Cuantilillos <- t(Cuantilillos)
colnames(Cuantilillos) <- c('Cuantiles Ajuste')
CuantilesA <- matrix(0,8,2)
colnames(CuantilesA) <- c('LimInf','LimSup')
rownames(CuantilesA) <- c('65','70','75','80','85','90','95','99')
CuantilesD <- matrix(0,8,2)
colnames(CuantilesD) <- c('LimInf','LimSup')
rownames(CuantilesD) <- c('65','70','75','80','85','90','95','99')
CuantilesA[1,1] <-CuantilesData[8]; CuantilesA[1,2] <-CuantilesData[9]
CuantilesA[2,1] <-CuantilesData[7]; CuantilesA[2,2] <-CuantilesData[10]
CuantilesA[3,1] <-CuantilesData[6]; CuantilesA[3,2] <-CuantilesData[11]
CuantilesA[4,1] <-CuantilesData[5]; CuantilesA[4,2] <-CuantilesData[12]
CuantilesA[5,1] <-CuantilesData[4]; CuantilesA[5,2] <-CuantilesData[13]
CuantilesA[6,1] <-CuantilesData[3]; CuantilesA[6,2] <-CuantilesData[14]
CuantilesA[7,1] <-CuantilesData[2]; CuantilesA[7,2] <-CuantilesData[15]
CuantilesA[8,1] <-CuantilesData[1]; CuantilesA[8,2] <-CuantilesData[16]
CuantilesD[1,1] <-Cuantilillos[8];  CuantilesD[1,2] <-Cuantilillos[9]
CuantilesD[2,1] <-Cuantilillos[7];  CuantilesD[2,2] <-Cuantilillos[10]
CuantilesD[3,1] <-Cuantilillos[6];  CuantilesD[3,2] <-Cuantilillos[11]
CuantilesD[4,1] <-Cuantilillos[5];  CuantilesD[4,2] <-Cuantilillos[12]
CuantilesD[5,1] <-Cuantilillos[4];  CuantilesD[5,2] <-Cuantilillos[13]
CuantilesD[6,1] <-Cuantilillos[3];  CuantilesD[6,2] <-Cuantilillos[14]
CuantilesD[7,1] <-Cuantilillos[2];  CuantilesD[7,2] <-Cuantilillos[15]
CuantilesD[8,1] <-Cuantilillos[1];  CuantilesD[8,2] <-Cuantilillos[16]

CuantilesData <- CuantilesA
CuantilesA <- CuantilesD
CuantilesD <- CuantilesData

dataset <- cbind(bdd,tst);
summary(dataset[,c('Sequence','BasalExp','tst')])
#write.csv(dataset,"ExpBasalDataset.csv")

umbral <- 0.0000001

tt1 <- min(tst)
VLI <- tt1
VLS <- CuantilesA[4,1] - umbral;
MLI <- CuantilesA[4,1]
MLS <- CuantilesA[1,1] - umbral;
MI  <- CuantilesA[1,1]
MS  <- CuantilesA[1,2] - umbral;
MHI <- CuantilesA[1,2]
MHS <- CuantilesA[4,2] - umbral;
VHI <- CuantilesA[4,2]
VHS <- max(tst);                
Limites <- matrix(0,1,10)
Limites <- c(VLI,VLS,MLI,MLS,MI,MS,MHI,MHS,VHI,VHS);
#print( tst)
print( Limites)

EBVL   <- dataset %>% filter(dataset$tst>=VLI & dataset$tst<VLS); 
EBML   <- dataset %>% filter(dataset$tst>=MLI & dataset$tst<MLS); 
EBM    <- dataset %>% filter(dataset$tst>=MI & dataset$tst<MS); 
EBMH   <- dataset %>% filter(dataset$tst>=MHI & dataset$tst<MHS); 
EBVH   <- dataset %>% filter(dataset$tst>=VHI & dataset$tst<=VHS); 
#print(EBVL)
#print(EBML)

n1 <- length(EBVL$tst);
n2 <- length(EBML$tst); 
n3 <- length(EBM$tst); 
n4 <- length(EBMH$tst);
n5 <- length(EBVH$tst); 

categoria <- rep('Muy Baja',n1);
EBVLCateg <- mutate(EBVL,categoria);#print(EBVLCateg)
categoria <- rep('Moderada Baja',n2);
EBMLCateg <- mutate(EBML,categoria);#print(EBMLCateg)
categoria <- rep('Moderada',n3);
EBMCateg  <- mutate(EBM,categoria); #print(EBMCateg)
categoria <- rep('Moderada Alta',n4);
EBMHCateg <- mutate(EBMH,categoria); #print(EBMHCateg)
categoria <- rep('Muy Alta',n5);
EBVHCateg <- mutate(EBVH,categoria); #print(EBVHCateg)

write.csv(EBVLCateg,"EBVLCategoria.csv")
write.csv(EBMLCateg,"EBMLCategoria.csv")
write.csv(EBMCateg,"EBMCategoria.csv")
write.csv(EBMHCateg,"EBMHCategoria.csv")
write.csv(EBVHCateg,"EBVHCategoria.csv")
#==================================================================


NormBasalExp <- rbind(EBVLCateg,EBMLCateg,
                      EBMCateg,EBMHCateg,EBVHCateg)
colnames(NormBasalExp)
ExpBasal <- NormBasalExp; #head(ExpBasal)
colnames(ExpBasal) <- c('Sequence','CodedSeq',
                        'Position','BasalExp',
                        'GenId','Log2Basal',
                        'Categoria')
head(ExpBasal,5)
#=======================================================
ExpBasal$Categoria = factor(ExpBasal$Categoria,
                           levels = c('Muy Baja',
                                      'Moderada Baja',
                                      'Moderada',
                                      'Moderada Alta',
                                      'Muy Alta'),
                           labels = c('Muy Baja',
                                      'Moderada Baja',
                                      'Moderada',
                                      'Moderada Alta',
                                      'Muy Alta'))
#=======================================================
summary(ExpBasal)

#=================================================================
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>

bexp <- ExpBasal$Log2Basal
Categorias <- ExpBasal$Categoria 

n <- length(bexp)
mirango <- c()
for (i in 1:n) {
  if(bexp[i]>=VHI){
    mirango[i] <- 'Very High'
  }else if (bexp[i]<=MHS & bexp[i]>=MHI){
    mirango[i] <- 'Moderate High'
  }else if (bexp[i]<=MS & bexp[i]>=MI){
    mirango[i] <- 'Moderate'
  }else if (bexp[i]<= MLS & bexp[i]>=MLI){
    mirango[i] <- 'Moderate Low'
  }else{
    mirango[i] <- 'Very Low'
  }
  }
mirango
plot_dat = data.frame(bexp, mirango)
plot_dat

plot_dat$mirango <- factor(plot_dat$mirango, levels = c("Very Low", 
                                                    "Moderate Low",
                                                    "Moderate", 
                                                    "Moderate High", 
                                                    "Very High"))
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
hist(plot_dat$bexp, )
nbreaks = 30
VLS
MLS
MS
MHS

# son 5 bins por cada unidad, es decir, cada bin es de aprox 0.2
vl  <-  rep("#104E8B",8)
ml  <-  rep("#8B2323",2); mivec <- c(vl,ml)
mmd <-  rep("#8B7355",8); mivec <- c(vl,ml,mmd)
mh  <-  rep("#8B2323",2); mivec <- c(vl,ml,mmd,mh)
vh  <-  rep("#104E8B",10); mivec <- c(vl,ml,mmd,mh,vh)
hist(bexp, breaks = nbreaks, col= mivec, 
     ylim = c(0,100),
     main = 'Normalized Log2 Basal Expresion',
     xlab='Basal Expresion',
     ylab= 'Frequency Basal Expresion')

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
LIMINF <- VLI-0.5
LIMSUP <- VHS + 0.5
ggplot(plot_dat, 
       aes(x=bexp, colour = factor(mirango))) + 
  geom_histogram(binwidth = 0.15) +
  ggtitle("Sample Histogram")+
  scale_color_manual(values = c("#104E8B",
                                "#8B2323", 
                                "#8B7355",
                                "#8B2323",
                                "#104E8B"))+
  scale_fill_manual(values=c("#104E8B",
                             "#8B2323", 
                             "#8B7355",
                             "#8B2323",
                             "#104E8B"))+
  lims(colour=c(mirango))+
  xlim(LIMINF,LIMSUP)+
  ylim(0,70)
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
LIMINF <- VLI-0.5
LIMSUP <- VHS + 0.5
ggplot(plot_dat, 
       aes(x=bexp, fill=mirango, colour = factor(mirango))) + 
  geom_histogram(binwidth = 0.125) +
  ggtitle("Sample Histogram")+
  scale_color_manual(values = c("#104E8B",
                                "#8B2323", 
                                "#8B7355",
                                "#8B2323",
                                "#104E8B"))+
  scale_fill_manual(values=c("#104E8B",
                             "#8B2323", 
                             "#8B7355",
                             "#8B2323",
                             "#104E8B"))+
  lims(colour=c(mirango))+
  xlim(LIMINF,LIMSUP)+
  ylim(0,70)
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>

EBVL   <- dataset %>% filter(dataset$tst>=VLI & dataset$tst<VLS); 
EBML   <- dataset %>% filter(dataset$tst>=MLI & dataset$tst<MLS); 
EBM    <- dataset %>% filter(dataset$tst>=MI & dataset$tst<MS); 
EBMH   <- dataset %>% filter(dataset$tst>=MHI & dataset$tst<MHS); 
EBVH   <- dataset %>% filter(dataset$tst>=VHI & dataset$tst<=VHS); 
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
ggplot(data=ExpBasal,
       aes(Log2Basal,
           fill=Categoria)) + 
  geom_histogram(bins = 60,
                 binwidth = 0.35) +
  theme_classic()+ 
  theme(legend.position="left")+
  scale_color_manual(values = c("#104E8B",
                                "#8B2323", 
                                "#8B7355",
                                "#8B2323",
                                "#104E8B"))+
  scale_fill_manual(values=c("#104E8B",
                             "#8B2323", 
                             "#8B7355",
                             "#8B2323",
                             "#104E8B"))+
  labs(title="Expresion Basal Normalizada",
       x="Categoria por Expresion Basal (Log2)",
       y = "Frecuencia")
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>



library(RColorBrewer)
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
tambinwidth <- 69
numbins <- 25
nbreaks = c(VLI,MLI,MI,MHI,VHI)

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
ggplot(data=ExpBasal,
       aes(Log2Basal)) + 
  geom_histogram(colour = 4, 
                 fill = "white", 
                 bins = 15,
                 position = "dodge") +
  theme(legend.position="right")+
  labs(title="Expresion Basal Normalizada",
       x="Categoria por Expresion Basal (Log2)",
       y = "Frecuencia")
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
ggplot(data=ExpBasal,
       aes(Log2Basal,
           fill=Categoria,
           color=Categoria))+
  geom_histogram(bins=30,
                 binwidth = 0.5,
                 position = "dodge")+
  theme(legend.position="right")+
  labs(title="Expresion Basal Normalizada",
       x="Categoria por Expresion Basal (Log2)",
       y = "Frecuencia")
# - - - - - - - -  - - - - - - - -  - - - - - - - -  - - - - - - - -  - - - - - - - -  - - - - - - - - 

#<> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
ggplot(data=ExpBasal,
       aes(Log2Basal,fill=Categoria)) + 
  geom_histogram(colour = 5, 
                 binwidth = 0.07,
                 bins = 25,
                 position = "dodge") +
  theme(legend.position="left")+
  theme_classic()+
  scale_fill_brewer(7,palette="Dark2")+
  labs(title="Expresion Basal Normalizada",
       x="Categoria por Expresion Basal (Log2)",
       y = "Frecuencia")
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>

ggplot(data=ExpBasal,
       aes(Log2Basal,fill=Categoria))+
  geom_histogram(binwidth = 0.5,
                 bins = 15,
                 position = "dodge",
                 alpha=0.85)+
  theme_classic()+ 
  theme(legend.position="left")+
  scale_color_manual(values = c("#0066FF",
                                "#33FF33", 
                                "#999999",
                                "#33FF33",
                                "#0066FF"))+
  scale_fill_manual(values=c("#0066FF",
                             "#33FF33", 
                             "#999999",
                             "#33FF33",
                             "#0066FF"))+
  labs(title="Expresion Basal Normalizada",
       x="Categoria por Expresion Basal (Log2)", 
       y = "Frecuencia")

ggplot(data=ExpBasal,
       aes(Log2Basal, fill = Categoria)) + 
  geom_histogram(colour = 4, fill = "white", 
                 bins = 25,
                 position = "dodge") +
  theme(legend.position="right")+
  labs(title="Expresion Basal Normalizada",
       x="Categoria por Expresion Basal (Log2)",
       y = "Frecuencia")

ggplot(data=ExpBasal,aes(Log2Basal,fill=Categoria))+
  geom_histogram(position = "dodge",alpha=0.85,bins=50,binwidth = 0.05)+
  theme_classic()+ theme(legend.position="left")+
  scale_color_manual(values = c("#104E8B","#8B2323", "#8B7355","#8B2323","#104E8B"))+
  scale_fill_manual(values=c("#104E8B","#8B2323", "#8B7355","#8B2323","#104E8B"))+
  labs(title="Expresion Basal Normalizada",x="Categoria por Expresion Basal (Log2)", y = "Frecuencia")

#==================================================================



#==================================================================
ggplot(data=ExpBasal,
       aes(Log2Basal,
           fill= Categoria))+
  geom_histogram(binwidth = 0.35, 
                 bins = 30)+ 
  theme_classic()+ 
  theme(legend.position="left")+
  scale_color_manual(values = c("#104E8B",
                                "#8B2323",
                                "#8B7355",
                                "#8B2323",
                                "#104E8B"))+
  scale_fill_manual(values=c("#104E8B",
                             "#8B2323", 
                             "#8B7355",
                             "#8B2323",
                             "#104E8B"))+
  labs(title="Expresion Basal Normalizada",
       x="Categoria por Expresion Basal (Log2)",
       y = "Frecuencia")
#=======================================================


# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
nbreaks = 30
# son 5 bins por cada unidad, es decir, cada bin es de aprox 0.2
vl  <-  rep("#104E8B",8)
ml  <-  rep("#8B2323",2); mivec <- c(vl,ml)
mmd <-  rep("#8B7355",8); mivec <- c(vl,ml,mmd)
mh  <-  rep("#8B2323",2); mivec <- c(vl,ml,mmd,mh)
vh  <-  rep("#104E8B",10); mivec <- c(vl,ml,mmd,mh,vh)
hist(bexp, breaks = nbreaks, col= mivec, 
     ylim = c(0,100),
     xlim = c(-4.5,4.5),
     main = 'Normalized Log2 Basal Expresion',
     xlab='Basal Expresion',
     ylab= 'Frequency Basal Expresion')
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
nbreaks = 30
# son 5 bins por cada unidad, es decir, cada bin es de aprox 0.2
vl  <-  rep("#104E8B",8)
ml  <-  rep("#8B2323",2); mivec <- c(vl,ml)
mmd <-  rep("#8B7355",8); mivec <- c(vl,ml,mmd)
mh  <-  rep("#8B2323",2); mivec <- c(vl,ml,mmd,mh)
vh  <-  rep("#104E8B",10); mivec <- c(vl,ml,mmd,mh,vh)
hist(bexp, breaks = nbreaks, col= mivec, 
     ylim = c(0,100),
     xlim = c(-4,4),
     main = 'Normalized Log2 Basal Expresion',
     xlab='Basal Expresion',
     ylab= 'Frequency Basal Expresion')
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
nbreaks = 30
# son 5 bins por cada unidad, es decir, cada bin es de aprox 0.2
vl  <-  rep("#104E8B",8)
ml  <-  rep("#8B2323",2); mivec <- c(vl,ml)
mmd <-  rep("#8B7355",8); mivec <- c(vl,ml,mmd)
mh  <-  rep("#8B2323",2); mivec <- c(vl,ml,mmd,mh)
vh  <-  rep("#104E8B",10); mivec <- c(vl,ml,mmd,mh,vh)
hist(bexp, breaks = nbreaks, col= mivec, 
     ylim = c(0,100),
     xlim = c(-3,3),
     main = 'Normalized Log2 Basal Expresion',
     xlab='Basal Expresion',
     ylab= 'Frequency Basal Expresion')
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>


# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
nbreaks = 30
# son 5 bins por cada unidad, es decir, cada bin es de aprox 0.2
vl  <-  rep("#104E8B",8)
ml  <-  rep("#8B2323",2); mivec <- c(vl,ml)
mmd <-  rep("#8B7355",8); mivec <- c(vl,ml,mmd)
mh  <-  rep("#8B2323",2); mivec <- c(vl,ml,mmd,mh)
vh  <-  rep("#104E8B",10); mivec <- c(vl,ml,mmd,mh,vh)
hist(bexp, breaks = nbreaks, col= mivec, 
     ylim = c(0,100),
     xlim = c(-3,3.9),
     main = 'Normalized Log2 Basal Expresion',
     xlab='Basal Expresion',
     ylab= 'Frequency Basal Expresion')
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>



# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
nbreaks = 30
# son 5 bins por cada unidad, es decir, cada bin es de aprox 0.2
vl  <-  rep("aliceblue",8)
ml  <-  rep("powderblue",2); mivec <- c(vl,ml)
mmd <-  rep("lightskyblue",8); mivec <- c(vl,ml,mmd)
mh  <-  rep("powderblue",2); mivec <- c(vl,ml,mmd,mh)
vh  <-  rep("aliceblue",10); mivec <- c(vl,ml,mmd,mh,vh)
hist(bexp, breaks = nbreaks, col= mivec, 
     ylim = c(0,100),
     xlim = c(-3,3.9),
     main = 'Normalized Log2 Basal Expresion',
     xlab='Basal Expresion',
     ylab= 'Frequency Basal Expresion')
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>


# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
nbreaks = 30
# son 5 bins por cada unidad, es decir, cada bin es de aprox 0.2
vl  <-  rep("gray100",8)
ml  <-  rep("gray83",2); mivec <- c(vl,ml)
mmd <-  rep("gray65",8); mivec <- c(vl,ml,mmd)
mh  <-  rep("gray83",2); mivec <- c(vl,ml,mmd,mh)
vh  <-  rep("gray100",10); mivec <- c(vl,ml,mmd,mh,vh)
hist(bexp, breaks = nbreaks, col= mivec, 
     ylim = c(0,100),
     xlim = c(-3,3.9),
     main = 'Normalized Log2 Basal Expresion',
     xlab='Basal Expresion',
     ylab= 'Frequency Basal Expresion')
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
nbreaks = 30
# son 5 bins por cada unidad, es decir, cada bin es de aprox 0.2
vl  <-  rep("lightcyan",8)
ml  <-  rep("honeydew",2); mivec <- c(vl,ml)
mmd <-  rep("lightblue",8); mivec <- c(vl,ml,mmd)
mh  <-  rep("honeydew",2); mivec <- c(vl,ml,mmd,mh)
vh  <-  rep("lightcyan",10); mivec <- c(vl,ml,mmd,mh,vh)
hist(bexp, breaks = nbreaks, col= mivec, 
     ylim = c(0,100),
     xlim = c(-3,3.9),
     main = 'Normalized Log2 Basal Expresion',
     xlab='Basal Expresion',
     ylab= 'Frequency Basal Expresion')
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
nbreaks = 30
# son 5 bins por cada unidad, es decir, cada bin es de aprox 0.2
vl  <-  rep("lightyellow",8)
ml  <-  rep("lightcyan",2); mivec <- c(vl,ml)
mmd <-  rep("lightblue",8); mivec <- c(vl,ml,mmd)
mh  <-  rep("lightcyan",2); mivec <- c(vl,ml,mmd,mh)
vh  <-  rep("honeydew",10); mivec <- c(vl,ml,mmd,mh,vh)
hist(bexp, breaks = nbreaks, col= mivec, 
     ylim = c(0,100),
     xlim = c(-3,3.9),
     main = 'Normalized Log2 Basal Expresion',
     xlab='Basal Expresion',
     ylab= 'Frequency Basal Expresion')
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>


# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
nbreaks = 30
# son 5 bins por cada unidad, es decir, cada bin es de aprox 0.2
vl  <-  rep("#FFFF99",8)
ml  <-  rep("#FFCC00",2); mivec <- c(vl,ml)
mmd <-  rep("#FF9933",8); mivec <- c(vl,ml,mmd)
mh  <-  rep("#FF3300",2); mivec <- c(vl,ml,mmd,mh)
vh  <-  rep("#990000",12); mivec <- c(vl,ml,mmd,mh,vh)
hist(bexp, breaks = nbreaks, col= mivec, 
     ylim = c(0,100),
     xlim = c(-3,3.9),
     main = 'Normalized Log2 Basal Expresion',
     xlab='Basal Expresion',
     ylab= 'Frequency Basal Expresion')
legend("topright",
       legend=c(" Very Low","Moderate Low","Moderate",
                "Moderated High","Very High"),
       pch = 15,
       col=c("#FFFF99","#FFCC00",
                       "#FF9933","#FF3300",
                       "#990000"),
       ncol = 1, 
       cex = 0.65)

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
nbreaks = 30
# son 5 bins por cada unidad, es decir, cada bin es de aprox 0.2
vl  <-  rep("#FFFF99",8)
ml  <-  rep("#FFCC00",2); mivec <- c(vl,ml)
mmd <-  rep("#FF9933",8); mivec <- c(vl,ml,mmd)
mh  <-  rep("#FF3300",2); mivec <- c(vl,ml,mmd,mh)
vh  <-  rep("#990000",12); mivec <- c(vl,ml,mmd,mh,vh)
m <- hist(bexp, breaks = nbreaks, col= mivec, 
     ylim = c(0,100),
     xlim = c(-3,3.9),
     main = 'Normalized Log2 Basal Expresion',
     xlab='Basal Expresion',
     ylab= 'Frequency Basal Expresion')
legend("topright",
       legend=c(" Very Low","Moderate Low","Moderate",
                "Moderated High","Very High"),
       pch = 15,
       col=c("#FFFF99","#FFCC00",
             "#FF9933","#FF3300",
             "#990000"),
       ncol = 1, 
       cex = 0.65)
print(m)
text(m$mids, m$counts, labels = m$counts,
     adj = c(0.5, -0.5))   
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
EBVL   <- dataset %>% filter(dataset$tst>=VLI & dataset$tst<VLS); 
EBML   <- dataset %>% filter(dataset$tst>=MLI & dataset$tst<MLS); 
EBM    <- dataset %>% filter(dataset$tst>=MI & dataset$tst<MS); 
EBMH   <- dataset %>% filter(dataset$tst>=MHI & dataset$tst<MHS); 
EBVH   <- dataset %>% filter(dataset$tst>=VHI & dataset$tst<=VHS); 
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
r1 <- seq(VLI,VLS,by= 0.1)
r2 <- seq(MLI,MLS,by= 0.1)
r3 <- seq(MI,MS,by= 0.1)
r4 <- seq(MHI,MHS,by= 0.1)
r5 <- seq(VHI,VHS+0.1,by= 0.1)
nbreaks <- c(r1,r2,r3,r4,r5)
Limites
n <- length(nbreaks)

vl  <-  rep("#FFFF99",16)
ml  <-  rep("#FFCC00",5); mivec <- c(vl,ml)
mmd <-  rep("#FF9933",18); mivec <- c(vl,ml,mmd)
mh  <-  rep("#FF3300",4); mivec <- c(vl,ml,mmd,mh)
vh  <-  rep("#990000",22); mivec <- c(vl,ml,mmd,mh,vh)
hist(bexp, breaks = nbreaks, col= mivec, 
     main = 'Normalized-Log2 Basal Expresion',
     xlim = c(-3,3.9),
     xlab='Basal Expresion',
     ylab= 'Frequency Basal Expresion')
legend("topright",
       legend=c(" Very Low","Moderate Low","Moderate",
                "Moderated High","Very High"),
       pch = 15,
       col=c("#FFFF99","#FFCC00",
             "#FF9933","#FF3300",
             "#990000"),
       ncol = 1, 
       cex = 0.65)
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
r1 <- seq(VLI,VLS,by= 0.1)
r2 <- seq(MLI,MLS,by= 0.1)
r3 <- seq(MI,MS,by= 0.1)
r4 <- seq(MHI,MHS,by= 0.1)
r5 <- seq(VHI,VHS+0.1,by= 0.1)
nbreaks <- c(r1,r2,r3,r4,r5)
Limites
n <- length(nbreaks)

vl  <-  rep("#FFFF99",16)
ml  <-  rep("#FFCC00",5); mivec <- c(vl,ml)
mmd <-  rep("#FF9933",18); mivec <- c(vl,ml,mmd)
mh  <-  rep("#FF3300",4); mivec <- c(vl,ml,mmd,mh)
vh  <-  rep("#990000",22); mivec <- c(vl,ml,mmd,mh,vh)
m <- hist(bexp, breaks = nbreaks, col= mivec, 
     main = 'Normalized-Log2 Basal Expresion',
     xlim = c(-3,3.9),
     xlab='Basal Expresion',
     ylab= 'Frequency Basal Expresion')
legend("topright",
       legend=c(" Very Low","Moderate Low","Moderate",
                "Moderated High","Very High"),
       pch = 15,
       col=c("#FFFF99","#FFCC00",
             "#FF9933","#FF3300",
             "#990000"),
       ncol = 1, 
       cex = 0.65)
print(m)
text(m$mids, m$counts, labels = m$counts,
     adj = c(0.5, -0.5))   
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>

