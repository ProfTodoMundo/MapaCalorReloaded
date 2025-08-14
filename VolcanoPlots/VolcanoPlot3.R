#---- Librerias ----
library(ggplot2)
library(dplyr)
library(ggrepel)

#---- (opcional) setwd ----
# setwd("~/Documents/GitHub/MapaCalorReloaded/VolcanoPlots")

#---- Lectura de datos ----
datos <- read.table("pEhExvsEhMyb10.txt",
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)

colnames(datos) <- c("Genes","baseMean","log2FC","lfcSE","stat","pvalue","padj",
                     "EhMyb10_1","EhMyb10_2","EhMyb10_3",
                     "pEhEx_1","pEhEx_2","pEhEx_3",
                     "EhMyb10_1.1","EhMyb10_2.1","EhMyb10_3.1",
                     "pEhEx_1.1","pEhEx_2.1","pEhEx_3.1",
                     "Certeza","DE")

# Parámetros
alpha_cut    <- 0.05
lfc_cut      <- 1.8
top_n_labels <- 50

#---- Procesamiento: SOLO pvalue ----
datos <- datos %>%
  mutate(
    p_selecc    = pvalue,
    p_selecc    = pmax(p_selecc, .Machine$double.xmin),
    menoslog10p = -log10(p_selecc),
    signif = ifelse(p_selecc < alpha_cut & abs(log2FC) > lfc_cut,
                    "Significativo", "No significativo"),
    direccion = case_when(
      p_selecc < alpha_cut & log2FC >=  lfc_cut ~ "Sobreexpresado",
      p_selecc < alpha_cut & log2FC <= -lfc_cut ~ "Subexpresado",
      TRUE ~ "No significativo"
    )
  )

#---- Etiquetas top ----
to_label <- datos %>%
  filter(signif == "Significativo") %>%
  arrange(desc(menoslog10p)) %>%
  slice_head(n = top_n_labels)

#---- Plot ----
p <- ggplot(datos, aes(x = log2FC, y = menoslog10p)) +
  geom_point(aes(color = direccion), alpha = 0.8, size = 2) +
  scale_color_manual(values = c(
    "Subexpresado"     = "#1f77b4",
    "No significativo" = "grey70",
    "Sobreexpresado"   = "#d62728"
  )) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
  geom_hline(yintercept = -log10(alpha_cut), linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = to_label, aes(label = Genes),
    size = 1.5, max.overlaps = 50, segment.size = 0.15
  ) +
  labs(
    title = "Volcano (p-value)",
    x = "log2(Fold Change)", y = "-log10(p)",
    color = "Dirección"
  ) +
  theme_minimal(base_size = 12)

# Mostrar
print(p)

# Guardar PDF
nombre0 <- paste0("VolcanoPlot_pvalue", top_n_labels, ".pdf")
pdf(nombre0); print(p); dev.off()

#---- Conteos & exports ----
conteo_colores <- datos %>%
  filter(!is.na(log2FC), !is.na(menoslog10p),
         direccion %in% c("Sobreexpresado", "Subexpresado")) %>%
  count(direccion, name = "n_genes")
write.csv(conteo_colores, paste0("cuantos_color_pvalue", top_n_labels, ".csv"), row.names = FALSE)

df_sobreexpresado <- datos %>%
  filter(!is.na(log2FC), !is.na(menoslog10p), direccion == "Sobreexpresado") %>%
  select(Genes, log2FC, menoslog10p, pvalue)
write.csv(df_sobreexpresado, paste0("Genes_sobreexpresados_pvalue", top_n_labels, ".csv"), row.names = FALSE)

df_subexpresado <- datos %>%
  filter(!is.na(log2FC), !is.na(menoslog10p), direccion == "Subexpresado") %>%
  select(Genes, log2FC, menoslog10p, pvalue)
write.csv(df_subexpresado, paste0("Genes_subexpresados_pvalue", top_n_labels, ".csv"), row.names = FALSE)
