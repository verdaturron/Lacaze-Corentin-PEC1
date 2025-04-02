# Script en R con comentarios explicativos


# Cargamos las librerías necesarias
library(SummarizedExperiment)  # Para estructurar los datos
library(readxl)                # Para leer archivos Excel
library(ggplot2)               # Para gráficos
library(pheatmap)              # Para mapas de calor
library(FactoMineR)            # Para análisis PCA
library(factoextra)            # Para visualizar PCA
library(reshape2)              # Para reformatear data.frames


## Función para crear un objeto SummarizedExperiment

create_se_object <- function(df, mode) {
  # Eliminamos columnas innecesarias
  df <- df[, !(colnames(df) %in% c("METABOLITE_ID", "Units"))]
  
  # rowData: nombres de metabolitos
  rowData <- DataFrame(MetaboliteName = df$`Metabolite Name`)
  
  # Extraemos matriz de datos numéricos
  assay_data <- as.matrix(df[, -1])                  # quitamos columna de nombres
  assay_data <- apply(assay_data, 2, as.numeric)     # convertimos a numérico
  
  # Creamos objeto SummarizedExperiment
  SummarizedExperiment(
    assays = list(counts = assay_data),               # matriz de intensidades
    rowData = rowData,                                # nombres de filas
    colData = DataFrame(Sample = colnames(assay_data)), # nombres de columnas
    metadata = list(mode = mode)                      # metadatos con tipo de modo
  )}


## Función para realizar un análisis PCA

run_pca <- function(se, mode) {
  mat <- t(assay(se))         # transponemos: muestras en filas
  mat <- log2(mat + 1)        # transformación log2 para estabilizar varianza
  pca <- PCA(mat, graph = FALSE)  # ejecutamos PCA sin graficar automáticamente
  fviz_pca_ind(pca, title = paste("PCA - Modo", mode))  # visualizamos con ggplot
}


## Función para generar un heatmap

plot_heatmap <- function(se, mode) {
  mat <- assay(se)
  mat <- log2(mat + 1)         # log-transformación
  mat <- t(scale(t(mat)))      # normalizamos por metabolito (z-score)
  
  # Generamos mapa de calor
  pheatmap(mat, main = paste("Heatmap - Modo", mode), show_rownames = FALSE)}


## Función para analizar variación individual de un sujeto

analyze_subject_variation <- function(df, subject_id) {
  # Creamos nombres de columnas para el sujeto (b: basal, a: manzana, c: arándano)
  cols <- paste0(c("b", "a", "c"), subject_id)
  
  # Verificamos que todas las columnas existan
  if (!all(cols %in% colnames(df))) {
    stop("Columnas faltantes para el sujeto ", subject_id)}
  
  # Extraemos solo las columnas del sujeto
  df_subj <- df[, cols]
  
  # Definimos los nombres de las filas
  rownames(df_subj) <- make.unique(as.character(df$`Metabolite Name`))
  
  # Convertimos los valores a numéricos
  df_subj <- as.data.frame(lapply(df_subj, as.numeric), row.names = rownames(df_subj))
  
  # Calculamos log2 fold change respecto al basal
  delta_apple <- log2(df_subj[[paste0("a", subject_id)]] + 1) -
    log2(df_subj[[paste0("b", subject_id)]] + 1)
  delta_cran  <- log2(df_subj[[paste0("c", subject_id)]] + 1) -
    log2(df_subj[[paste0("b", subject_id)]] + 1)
  
  # Unimos en un data.frame para graficar
  variations <- data.frame(
    Metabolito     = rownames(df_subj),
    Delta_Manzana  = delta_apple,
    Delta_Arandano = delta_cran)
  
  # Reformateamos para ggplot
  df_melt <- melt(variations, id.vars = "Metabolito")
  
  # Graficamos boxplot con ggplot
  ggplot(df_melt, aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    labs(
      title = paste("Variación - Sujeto", subject_id),
      x = "Condición",
      y = "Log2 Fold Change vs basal") + theme_minimal()}


## Ejecución del análisis completo

# Cargamos los datos desde archivos Excel
file_pos <- "Gu_pos-metabolites_PA-WB.xlsx"
file_neg <- "Gu_neg-metabolites_PA-WB.xlsx"

df_pos <- read_excel(file_pos)
df_neg <- read_excel(file_neg)

# Creamos objetos SummarizedExperiment
se_pos <- create_se_object(df_pos, "POS")
se_neg <- create_se_object(df_neg, "NEG")

# Ejecutamos PCA para ambos modos
run_pca(se_pos, "POS")
run_pca(se_neg, "NEG")

# Generamos mapas de calor
plot_heatmap(se_pos, "POS")
plot_heatmap(se_neg, "NEG")

# Analizamos variación individual para el sujeto 1
analyze_subject_variation(df_pos, 1)
