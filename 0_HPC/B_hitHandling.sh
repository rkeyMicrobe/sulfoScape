#!/bin/bash
#SBATCH --job-name=hitHandling
#SBATCH --output=logs/hitHandling.%j.out
#SBATCH --error=logs/hitHandling.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32GB
#SBATCH --time=02:00:00
#SBATCH --account=b.durham
#SBATCH --partition=hpg-default

echo "Started: $(date)"
echo "Hostname: $(hostname)"
echo "Working directory: $(pwd)"

echo "Defining variables and paths"
NAME="G2PA"
MODE="PA"

PROJ="/blue/b.durham/rebeccakey/sGrad_metaT"
HIT_DIR="${PROJ}/results/hits/${NAME}"
CNT_DIR="${PROJ}/data/countFiles"

OUTCSV="${HIT_DIR}/all_${NAME}_hmmHits.csv"
COUNTS="${CNT_DIR}/${NAME}_counts.csv"
TAX="${PROJ}/data/taxonomy/${NAME}_taxonomy.csv"
CNT_FILTERED="${HIT_DIR}/${NAME}_hitCounts2.csv"
TAX_FILTERED="${HIT_DIR}/${NAME}_hitTaxonomy2.csv"

#---------------------------------------------------------------
# Step 1
#echo "Combining HMMER CSVs..."
#hit_files=$(ls "${HIT_DIR}"/*_hmmHits_final.csv | grep -v "all_")

#first_file=$(echo "$hit_files" | head -n1)
#head -n 1 "$first_file" > "$OUTCSV"

#for f in "${HIT_DIR}"/*_hmmHits_final.csv; do
    tail -n +2 "$f" >> "$OUTCSV"
#done

#---------------------------------------------------------------
# Step 2
#echo "SubSet the countFiles to hmmerHits only..."
module purge; module load python/3.10

#python "/blue/b.durham/rebeccakey/sGrad_metaT/tools/filter_CountSubset.py" \
#  --mode "$MODE" \
#  --hmmer "$OUTCSV" \
#  --counts "$COUNTS" \
#  --output "$CNT_FILTERED" \
#  --chunksize 50000

#---------------------------------------------------------------
# Step 3
echo "Subset the taxonomyFiles to hmmerHits only..."
python "/blue/b.durham/rebeccakey/sGrad_metaT/tools/filter_TaxSubset.py" \
  --mode "$MODE" \
  --hmmer "$OUTCSV" \
  --taxonomy "$TAX" \
  --output "$TAX_FILTERED" \
  --chunksize 50000 \
  --tax_contig_col nt_id

echo "Finished processing script: $(date)"

