#---- Primer prueba ----
library(ggplot2)
setwd("~/Documentos/VolcanoPlots")
datos <- read.table("pEhExvsEhMyb10.txt", 
                    header = TRUE,    # Si la primera fila tiene nombres de columna
                    sep = "\t",       # \t = tabulador; usa "," si es CSV
                    stringsAsFactors = FALSE)
View(datos)
# Create a basic volcano plot
library(dplyr)
library(ggplot2)
library(ggrepel)

alpha_cut    <- 0.05
lfc_cut      <- 1
top_n_labels <- 10

df2 <- datos %>%
  mutate(
    fuente_p = ifelse(is.na(padj), "pvalue", "padj"),
    p_chosen = ifelse(is.na(padj), pvalue, padj),           # fallback a pvalue
    p_chosen = pmax(p_chosen, .Machine$double.xmin),        # evitar -Inf
    neglog10p = -log10(p_chosen),
    signif = case_when(
      p_chosen < alpha_cut & abs(log2FoldChange) > lfc_cut ~ "Significativo",
      TRUE ~ "No significativo"
    ),
    direccion = case_when(
      p_chosen < alpha_cut & log2FoldChange >=  lfc_cut ~ "Up",
      p_chosen < alpha_cut & log2FoldChange <= -lfc_cut ~ "Down",
      TRUE ~ "NS"
    )
  )

# Etiquetar los top por significancia
to_label <- df2 %>%
  filter(signif == "Significativo") %>%
  arrange(desc(neglog10p)) %>%
  slice_head(n = top_n_labels)

# Volcano plot (colores por dirección; forma por fuente de p-valor)
ggplot(df2, aes(x = log2FoldChange, y = neglog10p)) +
  geom_point(aes(color = direccion, shape = fuente_p), alpha = 0.8, size = 2) +
  scale_color_manual(values = c(Down = "#1f77b4", NS = "grey70", Up = "#d62728")) +
  scale_shape_manual(values = c(padj = 16, pvalue = 17)) +  # círculo = padj, triángulo = pvalue
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
  geom_hline(yintercept = -log10(alpha_cut), linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = to_label,
    aes(label = X),
    size = 3, max.overlaps = 50
  ) +
  labs(
    title = "Volcano plot con fallback a pvalue cuando padj es NA",
    x = "log2(Fold Change)",
    y = "-log10(p)",
    color = "Dirección",
    shape = "Fuente p"
  ) +
  theme_minimal(base_size = 12)

# (Opcional) resumen rápido de cuántos usaron padj vs pvalue
table(df2$fuente_p)


