# 1. Carreguem l'objecte
se <- readRDS("C:/Users/sarag/Desktop/PEC1_OMIQUES/summarized_experiment.rds")


# 2. Comprobació bàsica de l'objecte
library(SummarizedExperiment)
dim(se) # Dimensions de l'objecte
summary(colData(se)) # Resum de les dades de les mostres
summary(assay(se, "counts")) # Resum de les dades de l'expressió


# 3. Visualització de les dades
library(ggplot2)
df_long <- as.data.frame(t(assay(se, "counts")))
df_long <- stack(df_long) # Convertim a format llarg per visualitzar millor
ggplot(df_long, aes(x = values)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.5) +
  labs(title = "Distribució dels valors d'expressió",
       x = "Valor log-transformat",
       y = "Freqüència") +
  theme_minimal()


# 4. PCA per veure patrons globals
install.packages("FactoMineR")
install.packages("factoextra")
library(FactoMineR)
library(factoextra)
pca_res <- PCA(t(assay(se, "counts")), graph = FALSE) # PCA amb les dades d'expressió
fviz_pca_ind(pca_res, # Gràfic PCA
             col.ind = colData(se)$Muscle.loss, # Color segons Muscle.loss
             palette = c("red", "blue"),
             addEllipses = TRUE, # Afegir el·lipses per categoria
             title = "PCA de les mostres")


# 5. Boxplot per comparar grups
df_long$Muscle.loss <- rep(colData(se)$Muscle.loss, each = nrow(assay(se))) 
ggplot(df_long, aes(x = ind, y = values, fill = Muscle.loss)) +
  geom_boxplot() +
  labs(title = "Distribució dels valors d'expressió per grup",
       x = "Metabòlits",
       y = "Expressió log-transformat") +
  theme(axis.text.x = element_blank(), # Amagar noms de metabòlits
        legend.title = element_blank())


# 6. Comprovació de valors perduts
sum(is.na(assay(se, "counts")))  # Nombre total de valors NA


# 7. Visualització de valors perduts amb un mapa de calor
install.packages("heatmaply")
library(heatmaply)
heatmaply(
  is.na(assay(se, "counts")) * 1,  # Converteix TRUE/FALSE en 1/0
  main = "Mapa de calor de valors perduts")


# 8. Comprovació de la normalitat de les dades
# Histogrames per veure la distribució dels metabòlits
ggplot(df_long, aes(x = values)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.5) +
  facet_wrap(~ind, scales = "free_x") +  # Crear histogrames separats per metabòlits
  labs(title = "Distribució dels metabòlits",
       x = "Valor log-transformat",
       y = "Freqüència") +
  theme_minimal()

# Prova de normalitat de Shapiro-Wilk per a cada metabolit
shapiro_test <- apply(assay(se, "counts"), 1, function(x) shapiro.test(x)$p.value)
shapiro_test[shapiro_test < 0.05]  # Mostra els metabòlits amb p-value < 0.05 (no normal)


# 9. Correlació entre metabòlits
library(corrplot)
# Calcular la matriu de correlació entre les mostres
cor_matrix <- cor(assay(se, "counts"), use = "pairwise.complete.obs")
# Visualitzar la matriu de correlació
plot.new(); dev.off()
corrplot(cor_matrix, method = "circle", type = "upper", order = "hclust")
corrplot(cor_matrix, method = "circle", type = "upper", order = "hclust", tl.cex = 0.6)
dev.new(width = 10, height = 10)  # Obre una nova finestra gràfica més gran
corrplot(cor_matrix, method = "circle", type = "upper", order = "hclust")



# 10. ANOVA per comparar entre grups Muscle.loss
anova_res <- apply(assay(se, "counts"), 1, function(x) aov(x ~ colData(se)$Muscle.loss))
summary(anova_res[[1]])  # Comprovar resultats de la primera fila (metabolit)

