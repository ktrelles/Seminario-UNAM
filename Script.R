# ============================================================
# CASO PRÁCTICO: Intensidad turística en las provincias de Italia (2023)
# 1) Mapa coroplético  2) I de Moran  3) LISA  4) Modelos espaciales
#
# Carpeta de trabajo con: final_dataset.xlsx y ProvCM01012026_g_WGS84.(shp/dbf/shx/prj)
# Guardar este script en codificación UTF-8 (hay acentos en nombres de columnas).
# ============================================================
library(sf)
library(dplyr)
library(readxl)
library(ggplot2)
library(spdep)
library(spatialreg)

set.ZeroPolicyOption(TRUE)

setwd("C:/Users/Usuario/OneDrive - unibs.it/ADSEM/Caso práctico")



# ------------------------------------------------------------
# 0. DATOS (corte transversal 2023) Y SHAPEFILE
# ------------------------------------------------------------
panel <- read_excel("final_dataset.xlsx", sheet = "panel_regioni", na = "")
panel$ID_provincia <- as.character(panel$ID_provincia)  # "NA" = Napoli, no es dato faltante

datos <- panel %>%
  filter(anno == 2023) %>%
  transmute(
    ID        = ID_provincia,
    provincia = nome_provincia,
    tasso     = `Tasso turisticità`,
    prod      = `Produttività del lavoro`,
    partec    = `Partecipazione in mercato del lavoro`,
    ricavi    = `Ricavi imprese culturali`,
    inc_va    = `Incidenza valore aggiunto`,
    pop       = `Popolazione residente`
  ) %>%
  mutate(ln_tasso  = log(tasso),
         ln_prod   = log(prod),
         ln_ricavi = log(ricavi),
         ln_pop    = log(pop))

mapa <- st_read("ProvCM01012026_g_WGS84.shp", quiet = TRUE)

# Unir SOLO por ID (el nombre difiere entre fuentes: Aosta, Reggio Calabria, Forlì-Cesena)
map_data <- inner_join(mapa, datos, by = "ID")
stopifnot(nrow(map_data) == nrow(datos))   # deben ser 106 provincias
cat("Provincias con datos y geometría:", nrow(map_data), "\n")

# Para los modelos: misma tabla sin geometría, en el MISMO orden que map_data
datos_m <- st_drop_geometry(map_data)

# ------------------------------------------------------------
# 1. MATRICES W (se construyen desde map_data, así el orden coincide)
# ------------------------------------------------------------
nb_queen <- poly2nb(map_data, queen = TRUE)
W_queen  <- nb2listw(nb_queen, style = "W")
summary(nb_queen)   # revisar: nº de vecinos promedio y subgrafos desconectados (islas grandes)

coords <- st_coordinates(st_centroid(st_geometry(map_data)))  # UTM 32N (metros): distancia euclidiana
nb_knn <- knn2nb(knearneigh(coords, k = 4))
W_knn  <- nb2listw(nb_knn, style = "W")

# ------------------------------------------------------------
# 2. MAPA COROPLÉTICO
# ------------------------------------------------------------
p_mapa <- ggplot(map_data) +
  geom_sf(aes(fill = tasso), color = "white", linewidth = 0.15) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, trans = "log10",
                       name = "Tasso\nturisticità\n(escala log)") +
  labs(title = "Intensidad turística por provincia, 2023",
       caption = "Fuente: final_dataset.xlsx; ISTAT (límites 2026). Elaboración propia.") +
  theme_void()
p_mapa
ggsave("fig1_mapa.png", p_mapa, width = 7, height = 8, dpi = 300)

# ------------------------------------------------------------
# 3. I DE MORAN GLOBAL (variable dependiente del modelo: ln del tasso)
# ------------------------------------------------------------
moran_q <- moran.test(map_data$ln_tasso, W_queen)
moran_k <- moran.test(map_data$ln_tasso, W_knn)    # sensibilidad a W
set.seed(123)
moran_mc <- moran.mc(map_data$ln_tasso, W_queen, nsim = 999)
print(moran_q); print(moran_k); print(moran_mc)

# Diagrama de dispersión de Moran (la pendiente de la recta = I de Moran)
z  <- as.numeric(scale(map_data$ln_tasso))
wz <- lag.listw(W_queen, z)
p_scatter <- ggplot(data.frame(z, wz), aes(z, wz)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.7, color = "#7A1F3D") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(x = "Valor estandarizado de la provincia",
       y = "Promedio estandarizado de sus vecinas",
       title = "Diagrama de dispersión de Moran") +
  theme_minimal()
p_scatter
ggsave("fig2_moran.png", p_scatter, width = 6, height = 5, dpi = 300)

# ------------------------------------------------------------
# 4. LISA
# ------------------------------------------------------------
lisa   <- localmoran(map_data$ln_tasso, W_queen)
p_lisa <- lisa[, 5]   # p-valor (última columna)

map_data$cluster <- case_when(
  p_lisa >= 0.05  ~ "No significativo",
  z > 0 & wz > 0  ~ "Alto-Alto",
  z < 0 & wz < 0  ~ "Bajo-Bajo",
  z > 0 & wz < 0  ~ "Alto-Bajo",
  z < 0 & wz > 0  ~ "Bajo-Alto",
  TRUE            ~ "No significativo"
)
map_data$cluster <- factor(map_data$cluster,
                           levels = c("Alto-Alto", "Bajo-Bajo", "Alto-Bajo", "Bajo-Alto", "No significativo"))
print(table(map_data$cluster))

# OJO: con 106 pruebas simultáneas conviene corregir por pruebas múltiples
cat("Significativas sin corrección:", sum(p_lisa < 0.05),
    "| con FDR:", sum(p.adjust(p_lisa, method = "fdr") < 0.05), "\n")

cols_lisa <- c("Alto-Alto" = "#B2182B", "Bajo-Bajo" = "#4393C3",
               "Alto-Bajo" = "#F4A582", "Bajo-Alto" = "#92C5DE",
               "No significativo" = "#EEEEEE")
p_lisa_map <- ggplot(map_data) +
  geom_sf(aes(fill = cluster), color = "white", linewidth = 0.15) +
  scale_fill_manual(values = cols_lisa, name = "LISA") +
  labs(title = "Agrupamientos locales de intensidad turística (LISA)",
       caption = "p < 0.05 sin corrección por pruebas múltiples (exploratorio).") +
  theme_void()
p_lisa_map
ggsave("fig3_lisa.png", p_lisa_map, width = 7, height = 8, dpi = 300)

# ------------------------------------------------------------
# 5. MODELOS
# ------------------------------------------------------------
f <- ln_tasso ~ ln_prod + partec + ln_ricavi + inc_va + ln_pop

# 5.1 MCO y diagnósticos
ols <- lm(f, data = datos_m)
summary(ols)
lm.morantest(ols, W_queen)                 # Moran de los residuos
lm.RStests(ols, W_queen, test = "all")     # LM-error, LM-lag, versiones robustas, SARMA

# 5.2 Modelos espaciales
sar <- lagsarlm(f,   data = datos_m, listw = W_queen)
sem <- errorsarlm(f, data = datos_m, listw = W_queen)
sdm <- lagsarlm(f,   data = datos_m, listw = W_queen, type = "mixed")   # SDM
slx <- lmSLX(f,      data = datos_m, listw = W_queen)

# 5.3 Comparación (menor AIC = mejor)
aics <- c(OLS = AIC(ols), SAR = AIC(sar), SEM = AIC(sem), SDM = AIC(sdm), SLX = AIC(slx))
print(round(sort(aics), 1))

# Prueba de razón de verosimilitud: ¿SDM mejora a SAR? (H0: theta = 0)
LR.Sarlm(sdm, sar)

# 5.4 Modelo elegido: SDM -> efectos directos, indirectos y totales
summary(sdm)
set.seed(123)
imp_sdm <- impacts(sdm, listw = W_queen, R = 999)
summary(imp_sdm, zstats = TRUE, short = TRUE)

# 5.5 Sensibilidad a W (misma especificación con k vecinos)
sdm_knn <- lagsarlm(f, data = datos_m, listw = W_knn, type = "mixed")
cat("SDM queen: rho =", round(sdm$rho, 3), "| AIC =", round(AIC(sdm), 1), "\n")
cat("SDM k-NN : rho =", round(sdm_knn$rho, 3), "| AIC =", round(AIC(sdm_knn), 1), "\n")
set.seed(123)
summary(impacts(sdm_knn, listw = W_knn, R = 999), zstats = TRUE, short = TRUE)

# ============================================================
# VALORES DE REFERENCIA (validados en Python con estos mismos datos;
# en R pueden diferir un poco por el método de estimación y los p-valores):
#  - 106 provincias; W queen sin islas (3 componentes: península, Sicilia, Cerdeña)
#  - I de Moran de ln(tasso): ~0.18, p ~ 0.006 (queen y k-NN)
#  - MCO: R2 ~ 0.085; Moran de residuos ~ 0.12 (p ~ 0.04)
#  - LM-lag p ~ 0.04; LM-error p ~ 0.07 (queen); las versiones robustas NO son significativas
#  - AIC aprox.: OLS 343.8 | SAR 341.5 | SEM 339.7 | SDM 334.4 | SLX 334.4
#  - LR SDM vs SAR: p ~ 0.004
#  - LISA: ~21 provincias con p < 0.05; con FDR, ninguna
# ============================================================