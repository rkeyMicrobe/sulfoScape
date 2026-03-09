#!/bin/bash
#SBATCH --job-name=hmmer_array
#SBATCH --output=logs/hmmer_array.%A_%a.out
#SBATCH --error=logs/hmmer_array.%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32GB
#SBATCH --time=00:45:00
#SBATCH --account=b.durham
#SBATCH --qos=b.durham-b
#SBATCH --partition=hpg-default
#SBATCH --array=1-15%3

echo "Started: $(date)"
echo "Hostname: $(hostname)"
echo "Working directory: $(pwd)"

#---------------------------------------------------------------
# SETUP

PROJ="/blue/b.durham/rebeccakey/sGrad_metaT"
NAME="G2NS"
GENE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${PROJ}/geneList.txt")
QUERY="${PROJ}/genes/${GENE}.trimal2.fasta"
OUTDIR="${PROJ}/results/${NAME}_outputs/${GENE}"
CONTIG="${PROJ}/data/${NAME}_aaComb.fasta"

mkdir -p "$OUTDIR"

echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "Processing gene: $GENE"

#---------------------------------------------------------------
# Step 0: trimal

echo "Running trimAl QC on MSA..."
module purge; module load trimal
TRIMMED="${OUTDIR}/${GENE}.trimd.fasta"

trimal -in "$QUERY" \
       -out "$TRIMMED" \
       -gt 0.6 \
       -resoverlap 0.6 \
       -seqoverlap 50

ORIG_COUNT=$(grep -c '^>' "$QUERY")
TRIM_COUNT=$(grep -c '^>' "$TRIMMED")

echo "Original: $ORIG_COUNT seqs"
echo "Trimmed:  $TRIM_COUNT seqs"
echo "$GENE: $TRIM_COUNT seqs after trimming" \
>> "${PROJ}/results/hits/${NAME}_QC_summary.txt"

if [[ "$TRIM_COUNT" -eq 0 ]]; then
  echo "⚠ WARNING: No sequences remained after trimming for $GENE. Using untrimmed MSA instead."
  FINAL_MSA="$QUERY"
else
  FINAL_MSA="$TRIMMED"
fi

#---------------------------------------------------------------
# Step 1: seqmagick

echo "Converting to Stockholm format..."
module purge; module load seqmagick
seqmagick convert "$FINAL_MSA" "$OUTDIR/${GENE}.sto"

#---------------------------------------------------------------
# Step 2: hmmbuild

echo "Building HMM profile..."
module purge; module load hmmer
hmmbuild "${OUTDIR}/${GENE}.hmm" "${OUTDIR}/${GENE}.sto"

#---------------------------------------------------------------
# Step 3: hmmsearch

echo "Searching against $CONTIG..."
hmmsearch --cpu 2 \
  --tblout "${OUTDIR}/${GENE}.tbl" \
  "${OUTDIR}/${GENE}.hmm" "$CONTIG" \
  > "${OUTDIR}/${GENE}.out"

echo "Filtering hits by E-value..."
grep -v '^#' "${OUTDIR}/${GENE}.tbl" | awk '$5 <= 1e-5' > "${OUTDIR}/${GENE}_filtered.tbl"

echo "Total hits:"
wc -l "${OUTDIR}/${GENE}_filtered.tbl"

#---------------------------------------------------------------
# Step 4: CSV output

echo "Generating CSV for R..."
OUTCSV="${OUTDIR}/${GENE}_${NAME}_hmmHits.csv"
ORIGINAL="${OUTDIR}/${GENE}_${NAME}_hmmHits_noFilter.csv"


echo "TargetName,QueryName,Evalue,BitScore,Bias,ExpNumDomains,BestDomEvalue,BestDomScore,BestDomBias,Description" \
> "$OUTCSV"
grep -v '^#' "${OUTDIR}/${GENE}_filtered.tbl" | \
awk 'BEGIN {OFS=","} {print $1, $3, $5, $6, $7, $12, $9, $10, $11, $13}' \
>> "$OUTCSV"

echo "TargetName,QueryName,Evalue,BitScore,Bias,ExpNumDomains,BestDomEvalue,BestDomScore,BestDomBias,Description" \
> "$ORIGINAL"
grep -v '^#' "${OUTDIR}/${GENE}.tbl" | \
awk 'BEGIN {OFS=","} {print $1, $3, $5, $6, $7, $12, $9, $10, $11, $13}' \
>> "$ORIGINAL"


#---------------------------------------------------------------
# Step 5: Append sequences to the hits

module purge; module load gcc/5.2.0
module load bioawk/1.0

bioawk -c fastx '{print $name "\t" $seq}' "$CONTIG" > "$OUTDIR/seqs.tsv"

awk -F'\t' '
  NR==FNR { seq[$1]=$2; next }                 # read TSV
  FNR==1  { FS=","; OFS=","; print $0 ",Sequence"; next }  # header
  { split($1,a,/ /); key=a[1]; print $0 "," ((key in seq)?seq[key]:"") }
' "$OUTDIR/seqs.tsv" "$OUTCSV" \
> "${OUTCSV%.csv}_final.csv"

#---------------------------------------------------------------
# Step 5: push to hits dir

cp "${OUTCSV%.csv}_final.csv" "${PROJ}/results/hits/${NAME}"
echo "Finished processing $GENE: $(date)"




