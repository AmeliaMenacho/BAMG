############################################################
# 0. RUTA DE TRABAJO
############################################################
setwd("C:/Users/Amelia/Documents/Biotech/BAMG/Coding-02")
getwd()
list.files()
# Ver las primeras líneas del PED (solo texto)
readLines("genotype.ped", n = 3)

# Ver las primeras líneas del MAP
readLines("genotype.map", n = 5)

# Ver las primeras líneas del PHENO
readLines("phenotypes.pheno", n = 5)


# convertir PED + MAP a formato PLINK binario usando plink2.exe
system('./plink2 --ped genotype.ped --map genotype.map --make-bed --out arabidopsis')


# Probar PLINK2 (sin ./)
system("plink2.exe")

list.files()

system("plink2.exe --ped genotype.ped --map genotype.map --make-bed --out arabidopsis")

# Leer archivo de fenotipos
pheno <- read.table("phenotypes.pheno", header = TRUE)

# Ver las primeras columnas
head(pheno)[,1:10]

# Ver cuántas columnas tiene
ncol(pheno)

# Mostrar nombres de todas las columnas
colnames(pheno)

sort(colMeans(is.na(pheno)))[1:30]

############################################################
# 2. PREPARAR FENOTIPO PARA rMVP (FT10)
############################################################
library(dplyr)
library(readr)

# Cargar fenotipo
pheno <- read.table("phenotypes.pheno", header = TRUE)

# Seleccionar FID, IID y FT10
pheno_mvp <- pheno %>% 
  select(FID, IID, FT10)

# Guardar en CSV
write_csv(pheno_mvp, "pheno_ft10.csv")

head(pheno_mvp)

pheno <- read.table("phenotypes.pheno", header = TRUE)

# Crear fenotipo en formato rMVP: Taxa (ID) + FT10
pheno_mvp <- pheno %>% 
  transmute(Taxa = IID,
            FT10 = FT10)

# Ver cuántos datos NO NA tiene FT10
sum(!is.na(pheno_mvp$FT10))

# Guardar en CSV
write_csv(pheno_mvp, "pheno_ft10.csv")

head(pheno_mvp)


list.files()


###########
setwd("C:/Users/Amelia/Documents/BAMNG/BAMG/Coding-02")

library(dplyr)
library(readr)

# 1. Leer fenotipos originales
pheno_raw <- read.table("phenotypes.pheno", header = TRUE)

# 2. Armar fenotipo para rMVP: Taxa + FT10
pheno_mvp <- pheno_raw %>% 
  transmute(Taxa = IID,
            FT10  = FT10)

cat("N no-NA en FT10:", sum(!is.na(pheno_mvp$FT10)), "\n")


list.files(pattern = "\\.bed$")

bed_file <- list.files(pattern = "\\.bed$")

if (length(bed_file) == 0) stop("No hay ningún .bed en esta carpeta")

bed_prefix <- sub("\\.bed$", "", bed_file[1])
bed_prefix

############################################################
# 3. PREPARAR DATOS PARA rMVP (GENOTIPO + FENOTIPO)
############################################################
setwd("C:/Users/Amelia/Documents/Biotech/BAMG/Coding-02")


library(rMVP)

setwd("C:/Users/Amelia/Documents/Biotech/BAMG/Coding-02")

MVP.Data(
  fileBed  = "arabidopsis",     # usa arabidopsis.bed/bim/fam
  filePhe  = "pheno_FT10.csv",  # ahora SIN subcarpetas
  sep.phe  = ",",
  priority = "speed",
  fileKin  = FALSE,
  filePC   = FALSE,
  out      = "gwas_arabidopsis"
)

list.files(pattern = "gwas_arabidopsis")


############################################################
# 4. CALCULAR MATRIZ DE KINSHIP
############################################################
MVP.Data.Kin(
  fileKin    = TRUE,
  mvp_prefix = "gwas_arabidopsis",
  out        = "gwas_arabidopsis"
)

############################################################
# 5. LEER DATOS DE MVP Y ARMAR TRAIT
############################################################
library(bigmemory)
library(dplyr)

phenotype <- read.table("gwas_arabidopsis.phe", header = TRUE)
genotype  <- attach.big.matrix("gwas_arabidopsis.geno.desc")
mapping   <- read.table("gwas_arabidopsis.geno.map", header = TRUE)
kinship   <- attach.big.matrix("gwas_arabidopsis.kin.desc")

names(phenotype)[1:10]

# Trait: Taxa + FT10 (igual que en el script original con Seed.volume)
trait <- phenotype %>% select(Taxa, FT10)

head(trait)

############################################################
# 6. CREAR CARPETA DE SALIDA PARA EL GWAS
############################################################
setwd("C:/Users/Amelia/Documents/Biotech/BAMG/Coding-02")

dir.create("gwas_arabidopsis_output", showWarnings = FALSE)

sum(!is.na(trait$FT10))
############################################################
# 7. CORRER EL GWAS (GLM, MLM, FarmCPU)
############################################################
#GWAS
imMVP <- MVP(
  phe         = trait,
  geno        = genotype,
  map         = mapping,
  K           = kinship,
  method      = c("GLM", "MLM", "FarmCPU"),
  outpath     = "gwas_arabidopsis_output",
  file.output = TRUE,
  file.type   = "csv",
)
############################################################
# 8. GRÁFICOS: DENSIDAD, HISTOGRAMA, MANHATTAN, QQ
############################################################

# Densidad
MVP.Report(
  MVP       = imMVP,
  plot.type = "d",
  col       = c("darkgreen", "yellow", "red"),
  file.type = "jpg",
  dpi       = 300,
  outpath   = "gwas_arabidopsis_output"
)

# Histograma
MVP.Hist(
  phe       = trait,
  outpath   = "gwas_arabidopsis_output",
  file.type = "jpg",
  breakNum  = 18,
  dpi       = 300
)

# Manhattan
MVP.Report(
  imMVP,
  plot.type     = "m",
  multracks     = TRUE,
  threshold     = c(1e-6, 1e-4),
  threshold.lty = c(1,2),
  threshold.lwd = c(1,1),
  threshold.col = c("black","grey"),
  amplify       = TRUE,
  bin.size      = 1e6,
  chr.den.col   = c("darkgreen","yellow","red"),
  signal.col    = c("red","green"),
  signal.cex    = c(1,1),
  file.type     = "jpg",
  memo          = "",
  dpi           = 300,
  outpath       = "gwas_arabidopsis_output"
)

# Q-Q
MVP.Report(
  imMVP,
  plot.type  = "q",
  col        = c("dodgerblue1","olivedrab3","darkgoldenrod1"),
  threshold  = 1e6,
  signal.pch = 19,
  signal.cex = 1.5,
  signal.col = "red",
  box        = FALSE,
  multracks  = TRUE,
  file.type  = "jpg",
  memo       = "",
  dpi        = 300,
  outpath    = "gwas_arabidopsis_output"
)
