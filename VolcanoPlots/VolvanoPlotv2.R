#---- Librerias ----
library(ggplot2)
library(dplyr)
library(ggrepel)
#---- Directorio de trabajo ----
#setwd("~/Documentos/VolcanoPlots")
setwd("~/Documents/GitHub/MapaCalorReloaded/VolcanoPlots")
#---- Lectura de datos ----
datos <- read.table("pEhExvsEhMyb10.txt", 
                    header = TRUE,    # Si la primera fila tiene nombres de columna
                    sep = "\t",       # \t = tabulador; usa "," si es CSV
                    stringsAsFactors = FALSE)
summary(datos)
colnames(datos)<- c("Genes","baseMean","log2FC",
                    "lfcSE","stat","pvalue","padj",
                    "EhMyb10_1","EhMyb10_2","EhMyb10_3",
                    "pEhEx_1","pEhEx_2","pEhEx_3",
                    "EhMyb10_1.1","EhMyb10_2.1","EhMyb10_3.1",
                    "pEhEx_1.1","pEhEx_2.1","pEhEx_3.1",
                    "Certeza","DE" )
# Parametro central
top_n_labels <- 0
#---- Procesamiento de datos ----
alpha_cut    <- 0.05; lfc_cut      <- 1.8; 
##---- dataframe ----
datos <- datos %>%
  mutate(
    fuente_pvalue = padj,
    p_selecc = pvalue,
    p_selecc = pmax(p_selecc, .Machine$double.xmin),        # evitar -Inf
    menoslog10p = -log10(p_selecc),
    signif = case_when(
      p_selecc < alpha_cut & abs(log2FC) > lfc_cut ~ "Significativo", TRUE ~ "No Significativo"),
    direccion = case_when(
      p_selecc < alpha_cut & log2FC >=  lfc_cut ~ "Sobreexpresado",
      p_selecc < alpha_cut & log2FC <= -lfc_cut ~ "Subexpresado",
      TRUE ~ "No Significativo")
  )
##---- Etiquetar los top por nivel de significancia ----
to_label <- datos %>%
  filter(signif == "Significativo") %>% arrange(desc(menoslog10p)) %>% slice_head(n = top_n_labels)
print(to_label)
##---- Generacion de la grafica de Volcan ----
p <- ggplot(datos, aes(x = log2FC, y = menoslog10p)) +
  geom_point(aes(color = direccion), alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Subexpresado" = "#1f77b4", "No Significativo" = "grey70", "Sobreexpresado" = "#d62728")) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
  geom_hline(yintercept = -log10(alpha_cut), linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = to_label,aes(label = Genes),
    size = 1.5, max.overlaps = 50, segment.size=0.15
  ) +
  labs(
    title = "Gráfico para ______________",
    x = "log2(Fold Change)",y = "-log10(p)", color = "Dirección"
  ) +  theme_minimal(base_size = 12)
print(p)

nombre0 <- paste0("VolcanoPlot_pvalue0",top_n_labels,".pdf")
pdf(nombre0)
print(p)
dev.off()
#---- Contando ----
conteo_colores <- datos %>%
  filter(!is.na(log2FC), !is.na(menoslog10p),
         direccion %in% c("Sobreexpresado", "Subexpresado")) %>%
  count(direccion, name = "n_genes")
conteo_colores <- as.data.frame(conteo_colores); 
nombre1 <- paste0("cuantos_color_pvalue0",top_n_labels,".csv")
write.csv(conteo_colores,nombre1)
##----- Data frame para genes sobreexpresados -----
df_sobreexpresado <- datos %>%
  filter(!is.na(log2FC), !is.na(menoslog10p),
         direccion == "Sobreexpresado") %>%
  select(Genes, log2FC, menoslog10p, pvalue)
nombre3 <- paste0("Genes_sobreexpresados_pvalue0",top_n_labels,".csv")
write.csv(df_sobreexpresado,nombre3)
##----- Data frame para genes subexpresados -----
df_subexpresado <- datos %>%
  filter(!is.na(log2FC), !is.na(menoslog10p),
         direccion == "Subexpresado") %>%
  select(Genes, log2FC, menoslog10p, pvalue)
nombre4 <- paste0("Genes_subexpresados_pvalue0",top_n_labels,".csv")
write.csv(df_subexpresado,nombre4)