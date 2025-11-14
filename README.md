# Microbiota Bacteriana de Fauna Silvestre en Cautiverio
Este repositorio contiene los scripts, metadatos y resultados del flujo de análisis bioinformático, estadístico y funcional para caracterizar la diversidad bacteriana asociada a la fauna silvestre mantenida en cautiverio.

```
mkdir -p datasets
cd datasets
cat SraAccList.txt | xargs -n 1 -I{} fastq-dump --split-files --gzip -O raw-data {}
```

---

## 1. Análisis bioinformático

Las secuencias crudas fueron evaluadas con **FastQC v0.12.1** (Andrews, 2010).  
El procesamiento de datos se realizó con **QIIME2 v2025.7.0** (Bolyen et al., 2019) bajo **Linux/Ubuntu**:

1. **Corte de adaptadores y cebadores:** `q2-cutadapt` (Cutadapt v5.1; Martin, 2011)
2. **Filtrado y denoising:** `q2-dada2` (DADA2 v1.22.0; Callahan et al., 2016, 2017)
3. **Clasificación taxonómica:** `q2-feature-classifier` con clasificador **Naïve Bayes (SILVA v138)** y `--p-min-confidence 0.8`
4. **Filtrado posterior:** exclusión de muestras < 5000 lecturas y ASVs no clasificadas a nivel de filo o con < 3 lecturas.


```
conda create -env qiime2-amplicon
```

---

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

---

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

## Dependencias
Consulta el archivo [`environment.yml`](./environment.yml) para replicar el entorno con Conda:

```bash
conda env create -f environment.yml
conda activate microbiota_env

