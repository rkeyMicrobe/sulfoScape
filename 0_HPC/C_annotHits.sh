#!/bin/bash
#SBATCH --job-name=diamonTax
#SBATCH --output=logs/diamond.%A_%a.out
#SBATCH --error=logs/diamond.%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=12:00:00
#SBATCH --account=b.durham
#SBATCH --qos=b.durham-b
#SBATCH --partition=hpg-default

echo "Started: $(date)"
echo "Hostname: $(hostname)"
echo "Working directory: $(pwd)"

PROJ="/blue/b.durham/rebeccakey/sGrad_metaT"
TOOLS="${PROJ}/tools"

CRUISE="G2NS"

HIT_DIR="${PROJ}/results/hits/${CRUISE}"
HITS_CSV="${HIT_DIR}/all_${CRUISE}_hmmHits.csv"
DB_FASTA="${TOOLS}/marmicDB.faa"
DB_TAXID="${TOOLS}/marmicDB.uid2tax.tab"
DIAMOND="${TOOLS}/marmicrodb.dmnd"

# -----------------------------
# Modules
module purge
module load diamond

# -----------------------------
# 1) Build DIAMOND DB
echo "Build a diamond db"


#if [[ ! -f "${DIAMOND}" ]]; then
#  echo "Building DIAMOND DB at ${DIAMOND}"
#  diamond makedb --in "${DB_FASTA}" -d "${TOOLS}/diamond_MarFMARM"
#fi

# -----------------------------
# 2) Extract query FASTA from HMM hits
echo "Extract contig queries from HMM hits"

awk -F',' 'NR>1{
  seq=$11
  gsub(/\r/,"",seq)
  gsub(/[^A-Za-z\*\-]/,"",seq)   # keep AA,*,-
  gsub(/\*/,"",seq)              # drop stop codons
  print ">"$1"\n"seq
}' "${HITS_CSV}" > "${HIT_DIR}/dbQuery_seqs.faa"

# -----------------------------
# 3) Run DIAMOND blastp
echo "Running DIAMOND blastp..."

ANNOT="${HIT_DIR}/${CRUISE}_annot_MM2_n10.tsv"
diamond blastp \
  -q "${HIT_DIR}/dbQuery_seqs.faa" \
  -d "${DIAMOND}" \
  -o "${ANNOT}" \
  -p "${SLURM_CPUS_PER_TASK:-8}" --sensitive -k 10 --evalue 1e-10 --query-cover 70 --id 45 \
  -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle

# -----------------------------
# 4) Best hit per contig (lowest evalue, then highest bitscore)
echo "Select blastp seq with the lowest eval"

BEST_TSV="${HIT_DIR}/${CRUISE}_topHit_MM2_n10.tsv"
RAW_TSV="${HIT_DIR}/${CRUISE}_annot_MM2_n10.tsv"

#sort -k1,1 -k11,11g -k12,12r "${RAW_TSV}" | awk '!seen[$1]++' > "${BEST_TSV}"
cp "${RAW_TSV}" "${BEST_TSV}"

# -----------------------------
# 5) Build Taxid map
echo "Building UID → taxid map from decompressed DB.tab file"

MAP="${TOOLS}/uid2tax_MM2.map"
awk 'FNR>1 {print $2 "\t" $3}' "${DB_TAXID}" > "${MAP}"

# -----------------------------
# 6) Attach taxid
echo "Append taxID"

awk 'BEGIN{OFS="\t"}
     NR==FNR {tax[$1]=$2; next}
     {
       qseqid=$1; sseqid=$2; evalue=$11; stitle=$13;
       taxid=(sseqid in tax ? tax[sseqid] : "NA");
       print qseqid, sseqid, evalue, stitle, taxid
     }' "${MAP}" "${BEST_TSV}" > "${HIT_DIR}/${CRUISE}_topHit_wTax_MM2_n10.tsv"

{
  echo "qseqid,ref_uid,evalue,stitle,taxid"
  awk 'BEGIN{OFS=","} {gsub(/\t/,","); print}' "${HIT_DIR}/${CRUISE}_topHit_wTax_MM2_n10.tsv"
} > "${HIT_DIR}/${CRUISE}_topHit_wTax_MM2_n10.csv"

# -----------------------------
# 7) Clean

rm "${HIT_DIR}/${CRUISE}_topHit_wTax_MM2_n10.tsv"
rm "${HIT_DIR}/${CRUISE}_topHit_MM2_n10.tsv"
rm "${HIT_DIR}/${CRUISE}_annot_MM2_n10.tsv"
rm "${HIT_DIR}/dbQuery_seqs.faa"

# -----------------------------
# END--------------------------
# -----------------------------

echo "Done: $(date)"
