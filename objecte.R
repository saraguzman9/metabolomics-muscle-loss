# Instal·lar els paquets 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("SummarizedExperiment", "limma", "ggplot2", "pheatmap", "dplyr", "tidyverse"))

# Carregar els paquets
library(SummarizedExperiment)
library(limma)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyverse)


# URL del dataset
url <- "https://rest.xialab.ca/api/download/metaboanalyst/human_cachexia.csv"

# Llegeir el dataset directament des de l'URL
cachexia_data <- read.csv(url, header = TRUE, row.names = 1)

# Mostrar les primeres files
head(cachexia_data)


# Verificar el tipus de dades
str(cachexia_data)

# Comprovar si hi ha valors NA
sum(is.na(cachexia_data))


cachexia_data$Muscle.loss <- as.factor(cachexia_data$Muscle.loss)
levels(cachexia_data$Muscle.loss)  # Comprovar nivells

library(dplyr)

cachexia_data_norm <- cachexia_data %>%
  mutate(across(where(is.numeric), log1p))

cachexia_data_norm <- cachexia_data %>%
  mutate(across(where(is.numeric), log1p))  # Log(1+x) per evitar log(0)


head(col_data)  # Mostrar les primeres files
colnames(assay_data)  # Comprovar si aquests noms coincideixen amb rownames(col_data)
rownames(col_data) <- colnames(assay_data)  # Si cal, reassignar els noms correctament


dim(assay_data)  # Ha de ser (77, 63)
dim(col_data)  # Ha de ser (63, X), on X és el nombre de variables associades a cada mostra

# Assegurem que col_data només conté informació sobre les mostres
col_data <- DataFrame(Muscle.loss = cachexia_data_norm$Muscle.loss)

# Reassignem perquè només tingui 63 files (una per mostra)
col_data <- col_data[1:ncol(assay_data), , drop = FALSE]  

# Assignem els noms de les files per assegurar la correspondència
rownames(col_data) <- colnames(assay_data)

# Creem el SummarizedExperiment
se <- SummarizedExperiment(assays = list(counts = assay_data),
                           colData = col_data)

# Comprovem que tot sigui correcte
dim(se)
se

saveRDS(se, file = "summarized_experiment.rds")
saveRDS(se, file = "C:/Users/sarag/Desktop/PEC1_OMIQUES/summarized_experiment.rds")

