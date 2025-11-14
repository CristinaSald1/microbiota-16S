#!/bin/bash

# ===========================================================
#  Script: qiime2.sh
#  Description: Preparation of single-end or paired-end reads
#  and taxonomic classification of ASVs using QIIME 2.
#  Author: Cristina Saldaña
#  Date: 2025-11-08
# ===========================================================

set -e  # Stop execution on error

echo ""
echo "===================================================="
echo "             STARTING QIIME2 PIPELINE               "
echo "===================================================="
echo ""

# ===============================================
# ---------------- CONFIGURATION ----------------
# ===============================================

DATA_IN="datasets/clean-data"         # Cleaned fastq.gz files
METADATA="datatasets/metadata.tsv"    # Metadata file

QIIME_IN="qiime2/data"                # QIIME2 inputs
CUTAD_OUT="qiime2/cutadapt"           # Cutadapt outputs
DADA_OUT="qiime2/dada2"               # DADA2 outputs
TAXA_OUT="qiime2/taxonomy"            # Taxonomic classification results

ENV_NAME="qiime2-amplicon"            # QIIME2 environment name

# Primers/Adapters
Forward="TCTCGGCGAGGTCAGATGTATGAGAGAGAGAGAGAGCGGGGGWGCAGGWGCAG"       
Reverse="GTCTGGGGGCTGAGATGTAGTAGAGAGAGAGCAGACTCHVGTATCTATCCATC"           

# Trained classifier
CLASSIFIER="qiime2/classifier/silva-V3-V4-classifier.qza"

# Create output folders if they don't exist
mkdir -p "$QIIME_IN" "$CUTAD_OUT" "$DADA_OUT" "$TAXA_OUT"


# ===============================================
# --------------- QIIME2 ENVIRONMENT ------------
# ===============================================

echo "Activating environment $ENV_NAME..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME" || { echo "Failed to activate environment $ENV_NAME"; exit 1; }

# Check that QIIME2 is active in the correct environment
if [[ $(which qiime) != *"envs/$ENV_NAME/bin/qiime" ]]; then
    echo "Warning: QIIME2 does not appear to be running from the $ENV_NAME environment"
    echo "Detected path: $(which qiime)"
    exit 1
else
    echo "QIIME2 environment successfully activated"
fi


# ===============================================
# -------------- DETECT DATA TYPE ---------------
# ===============================================

if [ $(ls "$DATA_IN"/*_R2*.clean.fastq.gz 2>/dev/null | wc -l) -gt 0 ]; then
    TYPE="paired-end"
else
    TYPE="single-end"
fi
echo "Detected data type: $TYPE"
echo


# ===============================================
# ---------------- MANIFEST FILE ----------------
# ===============================================

MANIFEST="$QIIME_IN/manifest.csv"

if [ "$TYPE" == "single-end" ]; then
    echo -e "sample-id\tabsolute-filepath\tdirection" > $MANIFEST
    for f in $DATA_IN/*.clean.fastq.gz; do
        id=$(basename "$f" .clean.fastq.gz)
        abs_path=$(readlink -f "$f")
        echo -e "${id}\t${abs_path}\tforward" >> $MANIFEST
    done

else
    echo -e "sample-id\tabsolute-filepath-forward\tabsolute-filepath-reverse" > $MANIFEST
    for f in $DATA_IN/*_R1*.fastq.gz; do
        id=$(basename "$f" _R1.clean.fastq.gz)
        forward=$(readlink -f "$f")
        reverse=$(readlink -f "${DATA_IN}/${id}_R2.clean.fastq.gz")
        if [ ! -f "$reverse" ]; then
            echo "No se encontró el archivo reverse para $id"; exit 1
        fi
        echo "${id}\t${forward}\t${reverse}" >> $MANIFEST
    done

fi

echo "Manifest file generated at $MANIFEST"
echo

# ===============================================
# --------------- IMPORT TO QIIME2 --------------
# ===============================================

if [ "$TYPE" == "single-end" ]; then
    qiime tools import \
        --type 'SampleData[SequencesWithQuality]' \
        --input-path "$MANIFEST" \
        --output-path "$QIIME_IN/seqs.qza" \
        --input-format SingleEndFastqManifestPhred33V2

else
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "$MANIFEST" \
        --output-path "$QIIME_IN/seqs.qza" \
        --input-format PairedEndFastqManifestPhred33V2

fi

echo "Import data completed."
echo

# ===============================================
# ----------- PRIMER/ADAPTER TRIMMING -----------
# ===============================================

echo "Running Cutadapt to remove primers/adapters..."
if [ "$TYPE" == "single-end" ]; then
    qiime cutadapt trim-single \
        --i-demultiplexed-sequences "$QIIME_IN/seqs.qza" \
        --p-front "$Forward" \
        --p-error-rate 0.1 \
        --o-trimmed-sequences "$CUTAD_OUT/trimmed-seqs.qza" \
        
else
    qiime cutadapt trim-paired \
        --i-demultiplexed-sequences "$QIIME_IN/seqs.qza" \
        --p-front-f "$Forward" \
        --p-front-r "$Reverse" \
        --p-error-rate 0.1 \
        --o-trimmed-sequences "$CUTAD_OUT/trimmed-seqs.qza" \

fi

    qiime demux summarize \
        --i-data "$CUTAD_OUT/trimmed-seqs.qza" \
        --o-visualization "$CUTAD_OUT/trimmed-seqs.qzv"

echo "Trimming completed."
echo
   
# ===============================================
# ---------------- DADA2 ANALYSIS ---------------
# ===============================================

while true; do

    if [ "$TYPE" == "single-end" ]; then
        read -p "Enter truncation length (default 240): " TRUNC_LEN
        TRUNC_LEN=${TRUNC_LEN:-240}
        echo "Using truncation length: $TRUNC_LEN"

        qiime dada2 denoise-single \
            --i-demultiplexed-seqs "$CUTAD_OUT/trimmed-seqs.qza" \
            --p-trunc-len "$TRUNC_LEN" \
            --p-max-ee 2 \
            --o-representative-sequences "$DADA_OUT/denoising-seqs.qza" \
            --o-table "$DADA_OUT/table.qza" \
            --o-denoising-stats "$DADA_OUT/stats.qza"

    else
        read -p "Enter forward truncation length (default 240): " TRUNC_F
        read -p "Enter reverse truncation length (default 240): " TRUNC_R
        TRUNC_F=${TRUNC_F:-240}
        TRUNC_R=${TRUNC_R:-240}
        echo "Forward truncation: $TRUNC_F, Reverse truncation: $TRUNC_R"

        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs "$CUTAD_OUT/trimmed-sequences.qza" \
            --p-trunc-len-f "$TRUNC_F" \
            --p-trunc-len-r "$TRUNC_R" \
            --p-max-ee 2 \
            --o-representative-sequences "$DADA_OUT/denoising-seqs.qza" \
            --o-table "$DADA_OUT/table.qza" \
            --o-denoising-stats "$DADA_OUT/stats.qza"
    fi

    echo "DADA2 denoising completed."

    # DADA2 statistics visualization
    qiime metadata tabulate \
        --m-input-file "$DADA_OUT/stats.qza" \
        --o-visualization "$DADA_OUT/stats.qzv"

    # Summary table visualization
    qiime feature-table summarize \
        --i-table "$DADA_OUT/table.qza" \
        --o-visualization "$DADA_OUT/table-summary.qzv"

    # Representative sequences visualization
    qiime feature-table tabulate-seqs \
        --i-data "$DADA_OUT/denoising-seqs.qza" \
        --o-visualization "$DADA_OUT/denoising-seqs.qzv"

    echo "Results stored in $DADA_OUT"
    echo "Check the exported stats.qzv before proceeding."
    echo

    read -p "Continue with taxonomic classification? (Y/N): " ANSW
    if [[ "$ANSW" == "Y" || "$ANSW" == "y" ]]; then
        echo "Proceeding to taxonomic classification..."
        break
    else
        echo "Restarting DADA2 with new parameters..."
    fi
done

echo "Import and DADA2 steps completed successfully."
echo


# ===============================================
# ---------- TAXONOMIC CLASSIFICATION -----------
# ===============================================

echo
echo "===================================================="
echo "          STARTING TAXONOMIC CLASSIFICATION         "
echo "===================================================="
echo

# 1. Classify representative sequences
qiime feature-classifier classify-sklearn \
    --i-classifier "$CLASSIFIER" \
    --i-reads "$DADA_OUT/denoising-seqs.qza" \
    --p-confidence 0.8 \
    --o-classification "$TAXA_OUT/taxonomy.qza"

# 2. Visualize taxonomy table
qiime metadata tabulate \
    --m-input-file "$TAXA_OUT/taxonomy.qza" \
    --o-visualization "$TAXA_OUT/taxonomy.qzv"

# 3. Create taxonomy barplots
qiime taxa barplot \
    --i-table "$DADA_OUT/table.qza" \
    --i-taxonomy "$TAXA_OUT/taxonomy.qza" \
    --m-metadata-file "$METADATA" \
    --o-visualization "$TAXA_OUT/taxa-bar-plots.qzv"


echo
echo "===================================================="
echo "   TAXONOMIC CLASSIFICATION COMPLETED SUCCESSFULLY  "
echo "===================================================="
echo

echo "===================================================="
echo "       QIIME2 PIPELINE COMPLETED SUCCESSFULLY       "
echo "===================================================="
















