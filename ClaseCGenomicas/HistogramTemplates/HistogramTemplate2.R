setwd("~/Desktop/Datos Elisa")
library(ggplot2);   library(dplyr);         library(readxl);
library(pastecs);   library(sciplot);       library(MASS);
library(gridExtra); library("gplots");      library("lattice");
library(corrplot);  library(readr);         library(readxl);   
library(DBI);    
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
View(dataset)
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



nbreaks <- 30

tBE <- hist(bdd$BasalExp, breaks = nbreaks, col= rainbow(1,0.7), main = 'BasalExpresion')


BE <- bdd$BasalExp
par(mfrow=c(2,1))
Log2BE <- log2(BE)
nBE    <- length(Log2BE)
hist(Log2BE, breaks = nbreaks, col= rainbow(25,0.3), 
     main = ' Log2 Basal Expresion')
meanL2BE <- mean(Log2BE)
StdDevL2BE <- sd(Log2BE)
NormLog2BE <- (Log2BE-meanL2BE)/StdDevL2BE
tst<- NormLog2BE
hist(tst, breaks = nbreaks, col= 1:5, main = 'Normalized Log2 Basal Expresion',xlab='Basal Expresion',ylab= 'Frequency Basal Expresion')

fw1<-fitdist(tst, "norm")
plotdist(tst, histo = TRUE, demp = TRUE)



nnorm.f <- fitdist(tst,"norm")
summary(nnorm.f)
par(mfrow=c(2,2))
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



print(CuantilesD)



par(mfrow=c(2,1))
hist(tst, breaks = nbreaks, col= rainbow(1,0.7),
     main = 'Normalized Log2 Basal Expresion - DATA', lty=9)
#===========================================================================
abline(v=CuantilesA[1,1], lty=2, col="darkgoldenrod4"); # 65% INFERIOR
abline(v=CuantilesA[1,2], lty=2, col="darkgoldenrod4"); # 65% SUPERIOR
abline(v=CuantilesA[2,1], lty=2, col="darkblue"); # 70% INFERIOR
abline(v=CuantilesA[2,2], lty=2, col="darkblue") # 70% SUPERIOR
abline(v=CuantilesA[3,1], lty=2, col="aquamarine4"); # 75% INFERIOR
abline(v=CuantilesA[3,2], lty=2, col="aquamarine4"); # 75% SUPERIOR
abline(v=CuantilesA[4,1], lty=2, col="green");  # 80% INFERIOR
abline(v=CuantilesA[4,2], lty=2, col="green"); # 80% SUPERIOR
abline(v=CuantilesA[5,1], lty=2, col="brown"); # 85% INFERIOR
abline(v=CuantilesA[5,2], lty=2, col="brown"); # 85% SUPERIOR
abline(v=CuantilesA[6,1], lty=2, col="red");  # 90% INFERIOR
abline(v=CuantilesA[6,2], lty=2, col="red"); # 90% SUPERIOR
abline(v=CuantilesA[7,1], lty=2, col="blue");  # 95% INFERIOR
abline(v=CuantilesA[7,2], lty=2, col="blue"); # 95% SUPERIOR
abline(v=CuantilesA[8,1], lty=2, col="orange");  # 99% INFERIOR
abline(v=CuantilesA[8,2], lty=2, col="orange"); # 99% SUPERIOR
legend("topright",
       legend=c("65%","70%","75%","80%","85%","90%","95%","99%"), 
       pch=c(1,2,3,4,5,6,7,8),
       col=c("darkgoldenrod4","darkblue","aquamarine4",
             "green", "brown","red","blue","orange"))
hist(tst, breaks = nbreaks, col= rainbow(1,0.7),
     main = 'Normalized Log2   Basal Expresion - ADJUSTED', lty=9)
#===========================================================================
abline(v=CuantilesD[1,1], lty=2, col="darkgoldenrod4"); # 65% INFERIOR
abline(v=CuantilesD[1,2], lty=2, col="darkgoldenrod4"); # 65% SUPERIOR
abline(v=CuantilesD[2,1], lty=2, col="darkblue");  # 70% INFERIOR
abline(v=CuantilesD[2,2], lty=2, col="darkblue"); # 70% SUPERIOR
abline(v=CuantilesD[3,1], lty=2, col="aquamarine4");  # 75% INFERIOR
abline(v=CuantilesD[3,2], lty=2, col="aquamarine4"); # 75% SUPERIOR
abline(v=CuantilesD[4,1], lty=2, col="green");  # 80% INFERIOR
abline(v=CuantilesD[4,2], lty=2, col="green"); # 80% SUPERIOR
abline(v=CuantilesD[5,1], lty=2, col="brown");  # 85% INFERIOR
abline(v=CuantilesD[5,2], lty=2, col="brown"); # 85% SUPERIOR
abline(v=CuantilesD[6,1], lty=2, col="red");  # 90% INFERIOR
abline(v=CuantilesD[6,2], lty=2, col="red"); # 90% SUPERIOR
abline(v=CuantilesD[7,1], lty=2, col="blue");  # 95% INFERIOR
abline(v=CuantilesD[7,2], lty=2, col="blue"); # 95% SUPERIOR
abline(v=CuantilesD[8,1], lty=2, col="orange");  # 99% INFERIOR
abline(v=CuantilesD[8,2], lty=2, col="orange"); # 99% SUPERIOR
legend("topright",
       legend=c("65%","70%","75%","80%","85%","90%","95%","99%"), 
       pch=c(1,2,3,4,5,6,7,8),
       col=c("darkgoldenrod4","darkblue","aquamarine4",
             "green", "brown","red","blue","orange"))


par(mfrow=c(2,1))
hist(tst, breaks = nbreaks, col= rainbow(1,0.7),
     main = 'Normalized Log2  Basal Expresion - DATA', lty=9)
#===========================================================================
abline(v=CuantilesA[1,1], lty=2, col="darkgoldenrod4"); # 65% INFERIOR
abline(v=CuantilesA[1,2], lty=2, col="darkgoldenrod4"); # 65% SUPERIOR
abline(v=CuantilesA[4,1], lty=2, col="green");  # 80% INFERIOR
abline(v=CuantilesA[4,2], lty=2, col="green"); # 80% SUPERIOR
legend("topright",legend=c("65%","80%"), 
       pch=c(1,2),col=c("darkgoldenrod4","green"))
hist(tst, breaks = nbreaks, col= rainbow(1,0.7),
     main = 'Normalized Log2  Basal Expresion - ADJUSTED', lty=9)
#===========================================================================
abline(v=CuantilesD[1,1], lty=2, col="darkgoldenrod4"); # 65% INFERIOR
abline(v=CuantilesD[1,2], lty=2, col="darkgoldenrod4"); # 65% SUPERIOR
abline(v=CuantilesD[4,1], lty=2, col="green");  # 80% INFERIOR
abline(v=CuantilesD[4,2], lty=2, col="green"); # 80% SUPERIOR
legend("topright",legend=c("65%","80%"), 
       pch=c(1,2),col=c("darkgoldenrod4","green"))


par(mfrow=c(2,1))
hist(tst, breaks = nbreaks, col= rainbow(1,0.7), 
     main = 'Normalized Log2  Basal Expresion-DATA', lty=9)
#===========================================================================
abline(v=CuantilesA[2,1], lty=2, col="darkblue"); abline(v=CuantilesA[2,2], lty=2, col="darkblue")
abline(v=CuantilesA[5,1], lty=2, col="brown"); abline(v=CuantilesA[5,2], lty=2, col="brown")
legend("topright",legend=c("70%","85%"),
       pch=c(1,2),#3,4,5,6,7,8),
       col=c("brown"))
hist(tst, breaks = nbreaks, col= rainbow(1,0.7),
     main = 'Normalized Log2  Basal Expresion-ADJUSTED', lty=9)
#===========================================================================
abline(v=CuantilesD[2,1], lty=2, col="darkblue"); 
abline(v=CuantilesD[2,2], lty=2, col="darkblue")
abline(v=CuantilesD[5,1], lty=2, col="brown"); 
abline(v=CuantilesD[5,2], lty=2, col="brown")
legend("topright",legend=c("70%","85%"),
       pch=c(1,2),#3,4,5,6,7,8),
       col=c("brown"))


CuantilesData <- CuantilesA
CuantilesA <- CuantilesD
CuantilesD <- CuantilesData

counts<- table(bdd$Sequence)
PropSeq <- table(bdd$Sequence)
prop.table(PropSeq)
PC <- prop.table(PropSeq)
#print(PropSeq)


dataset <- cbind(bdd,tst);
summary(dataset[,c('Sequence','BasalExp','tst')])
#write.csv(dataset,"ExpBasalDataset.csv")





tt1 <- min(tst)
VLI <- tt1;             VLS <- CuantilesA[4,1] - 0.0000001;
MLI <- CuantilesA[4,1]; MLS <- CuantilesA[1,1] - 0.0000001;
MI  <- CuantilesA[1,1]; MS  <- CuantilesA[1,2] - 0.0000001; 
MHI <- CuantilesA[1,2]; MHS <- CuantilesA[4,2] - 0.0000001; 
VHI <- CuantilesA[4,2]; VHS <- max(tst);                
Limites <- matrix(0,1,10)
Limites <- c(VLI,VLS,MLI,MLS,MI,MS,MHI,MHS,VHI,VHS);
N <- length(tst)
ContVL<- 0; ContML<- 0; ContM <- 0; ContMH<- 0; ContVH<- 0;
for(i in 1:N){
  if((tst[i]>=VLI) & (tst[i]<=VLS)){ContVL<- ContVL+1;}
  if((tst[i]>=MLI) & (tst[i]<=MLS)){ContML<- ContML+1;}
  if((tst[i]>=MI)  & (tst[i]<=MS)){ContM <- ContM+1;}
  if((tst[i]>=MHI) & (tst[i]<=MHS)){ContMH<- ContMH+1;}
  if((tst[i]>=VHI) & (tst[i]<=VHS)){ContVH<- ContVH+1;}
}
Conteo <- matrix(0,2,6);
Conteo[1,1] <- ContVL;   Conteo[1,2] <- ContML
Conteo[1,3] <- ContM;    Conteo[1,4] <- ContMH
Conteo[1,5] <- ContVH;   Conteo[1,6] <- sum(Conteo[1,])
Conteo[2,1] <- ContVL/N; Conteo[2,2] <- ContML/N;
Conteo[2,3] <- ContM/N;  Conteo[2,4] <- ContMH/N;
Conteo[2,5] <- ContVH/N; Conteo[2,6] <- sum(Conteo[2,])
colnames(Conteo) <- c('VL','ML','M','MH','VH','Ttl')
rownames(Conteo) <- c('fr','Prob')
ProbClEB <- Conteo

Limites <- matrix(0,5,2)
Limites[1,1] <- VLI; Limites[1,2] <- VLS
Limites[2,1] <- MLI; Limites[2,2] <- MLS
Limites[3,1] <- MI;  Limites[3,2] <- MS
Limites[4,1] <- MHI; Limites[4,2] <- MHS
Limites[5,1] <- VHI; Limites[5,2] <- VHS
colnames(Limites) <- c('LimInf','Limsup')
rownames(Limites) <- c('VL','ML','M','MH','VH')
Limites
summary(dataset)



EBVL   <- dataset %>% filter(dataset$tst>=VLI & dataset$tst<=VLS); 
EBML   <- dataset %>% filter(dataset$tst>=MLI & dataset$tst<=MLS); 
EBM    <- dataset %>% filter(dataset$tst>=MI & dataset$tst<=MS); 
EBMH   <- dataset %>% filter(dataset$tst>=MHI & dataset$tst<=MHS); 
EBVH   <- dataset %>% filter(dataset$tst>=VHI & dataset$tst<=VHS); 

n1 <- length(EBVL$tst); categoria <- rep('Muy Baja',n1); 
EBVLCateg <- mutate(EBVL,categoria)
n2 <- length(EBML$tst); categoria <- rep('Moderada Baja',n2);
EBMLCateg <- mutate(EBML,categoria)
n3 <- length(EBM$tst);  categoria <- rep('Moderada',n3);
EBMCateg  <- mutate(EBM,categoria)
n4 <- length(EBMH$tst); categoria <- rep('Moderada Alta',n4);
EBMHCateg <- mutate(EBMH,categoria)
n5 <- length(EBVH$tst); categoria <- rep('Muy Alta',n5);
EBVHCateg <- mutate(EBVH,categoria)
#===========================================================================
NormBasalExp <- rbind(EBVLCateg,EBMLCateg,EBMCateg,EBMHCateg,EBVHCateg)
ExpBasal <- NormBasalExp[,c('GenId','BasalExp','tst','categoria')]
colnames(ExpBasal) <- c('GenId','BasalExp','Log2Basal','Categoria')
#===========================================================================
ExpBasal$Categoria = factor(ExpBasal$Categoria,
                           levels = c('Muy Baja','Moderada Baja','Moderada','Moderada Alta','Muy Alta'),
                           labels = c('Muy Baja','Moderada Baja','Moderada','Moderada Alta','Muy Alta'))
#===========================================================================
summary(ExpBasal)
#===========================================================================
bexp <- ExpBasal$Log2Basal
Categorias <- ExpBasal$Categoria 
library(RColorBrewer)
anchobin <- 15
numbins  <- 9
# - - - - - - - -  - - - - - - - -  - - - - - - - -  - - - - - - - -  - - - - - - - -  - - - - - - - - 
ggplot(data=ExpBasal,aes(Log2Basal,fill=Categoria,color=Categoria))+
  geom_histogram(binwidth = anchobin, bins = numbins,position = "dodge")+
  theme(legend.position="right")+
  labs(title="Expresion Basal Normalizada",x="Categoria por Expresion Basal (Log2)", y = "Frecuencia")
# - - - - - - - -  - - - - - - - -  - - - - - - - -  - - - - - - - -  - - - - - - - -  - - - - - - - - 
ggplot(data=ExpBasal,aes(Log2Basal,fill=Categoria,color=Categoria))+
  geom_histogram(binwidth = anchobin, bins = numbins,position = "dodge")+ theme(legend.position="right")+
  theme(legend.position="right")+theme_classic()
  labs(title="Expresion Basal Normalizada",x="Categoria por Expresion Basal (Log2)", y = "Frecuencia")
# - - - - - - - -  - - - - - - - -  - - - - - - - -  - - - - - - - -  - - - - - - - -  - - - - - - - - 
ggplot(data=ExpBasal,aes(Log2Basal,fill=Categoria))+geom_histogram(binwidth = anchobin, bins = numbins,position = "dodge")+
  theme(legend.position="left")+scale_fill_brewer(7,palette="GnBu")
ggplot(data=ExpBasal,aes(Log2Basal,fill=Categoria))+geom_histogram(binwidth = anchobin, bins = numbins,position = "dodge")+
  theme(legend.position="left")+scale_fill_brewer(7,palette="YlOrRd")
#===========================================================================

ggplot(data=ExpBasal,aes(Log2Basal,fill=Categoria))+geom_histogram(binwidth = anchobin, bins = numbins,position = "dodge")+
  theme_classic()+scale_fill_brewer(7,palette="Dark2")+ theme(legend.position="left")+
  labs(title="Expresion Basal Normalizada",x="Categoria por Expresion Basal (Log2)", y = "Frecuencia")
#===========================================================================

ggplot(data=ExpBasal,aes(Log2Basal,fill=Categoria))+
  geom_histogram(binwidth = anchobin, bins = numbins,position = "dodge",alpha=0.85)+
  theme_classic()+ theme(legend.position="left")+
  scale_color_manual(values = c("#104E8B","#8B2323", "#8B7355","#8B2323","#104E8B"))+
  scale_fill_manual(values=c("#104E8B","#8B2323", "#8B7355","#8B2323","#104E8B"))+
  labs(title="Expresion Basal Normalizada",x="Categoria por Expresion Basal (Log2)", y = "Frecuencia")
#===========================================================================

ggplot(data=ExpBasal,aes(Log2Basal,fill=Categoria))+
  geom_histogram(position = "dodge",alpha=0.85,bins = numbins,binwidth = anchobin)+
  theme_classic()+ theme(legend.position="left")+
  scale_color_manual(values = c("#104E8B","#8B2323", "#8B7355","#8B2323","#104E8B"))+
  scale_fill_manual(values=c("#104E8B","#8B2323", "#8B7355","#8B2323","#104E8B"))+
  labs(title="Expresion Basal Normalizada",x="Categoria por Expresion Basal (Log2)", y = "Frecuencia")

#===========================================================================

ggplot(data=ExpBasal,aes(Log2Basal,fill=Categoria))+geom_histogram(binwidth = anchobin, bins = numbins,
                                                                   position = "dodge")+
  theme(legend.position="left")+scale_fill_brewer(7,palette="GnBu")+
  theme_classic()+ theme(legend.position="left")+
  labs(title="Expresion Basal Normalizada",x="Categoria por Expresion Basal (Log2)", y = "Frecuencia")

#===========================================================================

ggplot(data=ExpBasal,aes(Log2Basal,fill=Categoria))+geom_histogram(bins = numbins,
                                                                   position = "dodge")+
  theme(legend.position="left")+scale_fill_brewer(7,palette="YlOrRd")

#===========================================================================
ggplot(data=ExpBasal,aes(Log2Basal,color=Categoria))+geom_histogram(binwidth = anchobin, bins = numbins,fill='white',
                                                                    position = "dodge",
                                                                    linetype = "dashed")

#===========================================================================
ggplot(data=ExpBasal,aes(Log2Basal,fill= Categoria))+
  geom_histogram(binwidth = anchobin, bins = numbins)+ 
  theme_classic()+ theme(legend.position="left")+
  scale_color_manual(values = c("#104E8B","#8B2323", "#8B7355","#8B2323","#104E8B"))+
  scale_fill_manual(values=c("#104E8B","#8B2323", "#8B7355","#8B2323","#104E8B"))+
  labs(title="Expresion Basal Normalizada",x="Categoria por Expresion Basal (Log2)", y = "Frecuencia")
#===========================================================================








