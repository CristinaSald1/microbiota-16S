# Microbiota bacteriana y potencial funcional asociado a virulencia en fauna en cautiverio

Este repositorio contiene los scripts, metadatos y resultados del análisis bioinformático, estadístico y funcional orientado a caracterizar la microbiota bacteriana presente en fauna mantenida en cautiverio y determinar su potencial funcional, con énfasis en la detección de rutas metabólicas asociadas a virulencia y factores de patogenicidad.


## Hipótesis

La composición bacteriana y las funciones metabólicas difieren entre especies hospederas, y algunas comunidades microbianas presentan rutas asociadas a virulencia que podrían representar riesgos zoonóticos.


## Objetivos

### Objetivo General

Caracterizar la microbiota bacteriana asociada a fauna mantenida en cautiverio y determinar su potencial funcional, enfocándose en rutas metabólicas relacionadas con virulencia y factores de patogenicidad.


### Objetivo Específicos

1. Procesar las secuencias 16S rRNA mediante un flujo estandarizado de QIIME2 para obtener tablas de ASVs y clasificación taxonómica.
2. Evaluar la diversidad bacteriana e identificar taxones diferencialmente abundantes entre especies hospederas.
3. Inferir el potencial funcional microbiano y detectar rutas metabólicas relacionadas con virulencia.


## Muestra
Este repositorio se encuentra inicialmente estructurado con fines de prueba y estandarización utilizando un dataset público.

- **Número total de muestras:** 37
- **Organismos hospederos:** ganado vacuno de cuatro regiones de Kazajistán (Western, Southern, Northern and Southeast)
- **Tipo de muestra:** heces
- **Tecnología:** 16S rRNA (Illumina)
- **Bioproject ID:** PRJNA847733

**Referencia:** 
Daugaliyeva, A., Daugaliyeva, S., Ashanin, A. et al. (2022). Study of cattle microbiota in different regions of Kazakhstan using 16S metabarcoding analysis. Scientific Reports, 12, 16410. https://doi.org/10.1038/s41598-022-20732-4


---

## Métodos

## 1. Análisis bioinformático
Las secuencias crudas fueron evaluadas con **FastQC v0.12.1** (Andrews, 2010).  
El procesamiento de datos se realizó con **QIIME2 v2025.7.0** (Bolyen et al., 2019) bajo **Linux/Ubuntu**:

1. **Corte de adaptadores y cebadores:** `q2-cutadapt` (Cutadapt v5.1; Martin, 2011)
2. **Filtrado y denoising:** `q2-dada2` (DADA2 v1.22.0; Callahan et al., 2016, 2017)
3. **Clasificación taxonómica:** `q2-feature-classifier` con clasificador **Naïve Bayes (SILVA v138)** y `--p-min-confidence 0.8`
4. **Filtrado posterior:** exclusión de muestras < 5000 lecturas y ASVs no clasificadas a nivel de filo o con < 3 lecturas.

--

## 2. Análisis estadísticos
Los archivos de salida de QIIME2 se importaron a **R v4.5.0** mediante `qiime2R`.

### Diversidad alfa y beta
Usando:
- `phyloseq v1.44.0` (McMurdie & Holmes, 2013)
- `microbiome v1.30.0` (Shetty, 2017)
- `vegan v2.8-0` (Oksanen et al., 2025)

1. Índices de **Shannon**, **Simpson** y **Chao1**  
2. **Kruskal–Wallis/Wilcoxon** para diferencias entre especies hospederas  
3. **Bray–Curtis** y **UniFrac** para diversidad beta (PCoA visualizations)  
4. **PERMANOVA** (`adonis2`) para diferencias multivariadas  
5. Transformación **CLR (Centered Log-Ratio)** para evitar sesgos composicionales

### Abundancia diferencial
Los taxones diferencialmente abundantes se identificaron con **DESeq2 v3.21** (Love et al., 2014) usando un **modelo binomial negativo**.  
Los resultados se reportan como **log₂ fold-change (LFC)** con **FDR < 0.05**.

--

## 3. Análisis funcional y ecológico

### Predicción funcional
Predicción de metagenomas con **q2-PICRUSt2** (Douglas et al., 2020), validada con **Tax4Fun2** (Wemheuer et al., 2020).  
La anotación funcional se basa en **KEGG Orthology (KO)** (Kanehisa & Goto, 2000).

### Interacción ecológica y redes
Los perfiles funcionales se integraron con datos taxonómicos para definir comunidades bacterianas funcionales por especie hospedera.  
Se realizaron:
- **PERMANOVA** y **agrupamiento jerárquico**  
- Construcción de redes de co-ocurrencia con **microbiomeSeq**  
- Validación de géneros zoonóticos mediante **VFDB (Zhou et al., 2024)** y **PATRIC (Wattam et al., 2017)**

---

## Ejecución del pipeline

Descar secuencias desde SRA
```
mkdir -p datasets
cd datasets
cat SraAccList.txt | xargs -n 1 -I{} fastq-dump --split-files --gzip -O raw-data {}
```

Análisis de calidad de secuencias
```
mkdir -p fasqc_reports
fastqc raw-data/*.fastq.gz -o fasqc_reports/
multiqc fasqc_reports/ -o fasqc_reports/
```


Crear entorno QIIME2

```
# Actualizar conda
conda update conda
```
```
# Crear un entorno conda para qiime2 v25.7
conda env create \
  --name qiime2-amplicon \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.7/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml
```
- Consultar guía oficial para obtener otras versiones de qiime2: https://library.qiime2.org/quickstart/amplicon

```
# Verificar la creación del entorno
conda deactivate
conda activate qiime2-amplicon
qiime info
```


Ejecutar scripts del repositorio

```
# Generar tablas de ASVs y clasificación taxonómica
bash scripts/qiime2.sh
```
```
# Analizar la diversidad bacteriana
bash scripts/stats.R
```
