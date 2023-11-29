#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
if (!requireNamespace("magrittr"))install.packages("magritt");           library(magrittr)
if (!requireNamespace("pheatmap"))install.packages("pheatmap");          library(pheatmap)
if (!requireNamespace("RColorBrewer"))install.packages("RcolorBrewer");  library(RColorBrewer)
if (!requireNamespace("Rio"))install.packages("rio");                    library(rio)
if (!requireNamespace("readr"))install.packages("readr");                library(readr)
if (!requireNamespace("tidyverse"))install.packages("");                 library(tidyverse)
if (!requireNamespace("mclust"))install.packages("mclust");              library(mclust)
if (!requireNamespace("venn"))install.packages("venn");                  library(venn)
if (!requireNamespace("dplyr"))install.packages("dplyr");                library(dplyr)
if (!requireNamespace("ggplot2"))install.packages("ggplot2");            library(ggplot2)
if (!requireNamespace("cowplot"))install.packages("cowplot");            library(cowplot)
if (!requireNamespace("RColorBrewer"))install.packages("RColorBrewer");  library(RColorBrewer)
if (!requireNamespace("ggVennDiagram"))install.packages("ggVennDiagram");library(ggVennDiagram)
if (!requireNamespace("VennDiagram"))install.packages("VennDiagram");    library(VennDiagram)
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
setwd("~/Desktop/MiGithub/MapaCalorReloaded")
library(readr)
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
NewGenes <- read_csv("Ei_enq_circ.csv")
View(NewGenes)
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
NewGenes$Trophs <- NewGenes$Trophs-2
NewGenes$`cyst 8h` <- NewGenes$`cyst 8h`-2;
NewGenes$`cyst 24h` <- NewGenes$`cyst 24h`-2
NewGenes$`cyst 48h` <- NewGenes$`cyst 48h`-2
NewGenes$`cyst 72h` <- NewGenes$`cyst 72h`-2
NewGenes$`excyst 2h` <- NewGenes$`excyst 2h`-2
NewGenes$`excyst 8h` <- NewGenes$`excyst 8h`-2
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
colnames(NewGenes) <- c("GeneID","Trophozoites", "Cyst_8h","Cyst_24h",
                        "Cyst_48h","Cyst_72h","Excyst_2h","Excyst_8h")
rownames(NewGenes) <- NewGenes$GeneID
#<< == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> << == >> 
GenesClean <- NewGenes[,2:8]; View(GenesClean)
rownames(GenesClean) <- rownames(NewGenes)
GenesClean <- as.data.frame(GenesClean)
rownames(GenesClean) <- NewGenes$GeneID