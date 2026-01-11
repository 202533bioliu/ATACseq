#!/bin/bash
set -e
set -o pipefail

# === Configuration ===
BASE_DATA_DIR="/home/icac_yl/ATAC/fastq"
DATASETS=("GSE279911" "GSE103588")
BOWTIE2_INDEX="/home/icac_yl/ATAC/bt2_index/hg38"
BASE_OUTPUT_DIR="/home/icac_yl/ATAC/fastq"
SAMPLE_SUFFIX="_1.fastq"
THREADS=16
GENOME_SIZE="hs"
FILTER_CHRS="chrM|MT"

# === Process each dataset ===
for DATASET in "${DATASETS[@]}"; do
RAW_DATA_DIR="${BASE_DATA_DIR}/${DATASET}"
OUTPUT_DIR="${BASE_OUTPUT_DIR}/${DATASET}/atac_analysis_results2"

# === Create output directories ===
mkdir -p ${OUTPUT_DIR}/{trimmed,aligned,bam_filtered,peaks,annotation,tss_heatmap,merged_peaks,bigwig,trimmed_fastqc}

# === Get sample list ===
SAMPLE_LIST=$(ls ${RAW_DATA_DIR}/*${SAMPLE_SUFFIX} | sed -e "s|${RAW_DATA_DIR}/||g" -e "s|${SAMPLE_SUFFIX}||g")

# === Trimming: TrimGalore ===
for sample in ${SAMPLE_LIST}; do
R1=${RAW_DATA_DIR}/${sample}_1.fastq
R2=${RAW_DATA_DIR}/${sample}_2.fastq
trim_galore --paired --quality 20 --phred33 --length 20 --gzip \
--output_dir ${OUTPUT_DIR}/trimmed \
--cores ${THREADS} ${R1} ${R2}
done

# === FastQC on trimmed files ===
for sample in ${SAMPLE_LIST}; do
R1_TRIM=${OUTPUT_DIR}/trimmed/${sample}_1_val_1.fq.gz
R2_TRIM=${OUTPUT_DIR}/trimmed/${sample}_2_val_2.fq.gz
fastqc -o ${OUTPUT_DIR}/trimmed_fastqc -t ${THREADS} ${R1_TRIM} ${R2_TRIM}
done

# === Alignment and BAM processing ===
for sample in ${SAMPLE_LIST}; do
R1_TRIM=${OUTPUT_DIR}/trimmed/${sample}_1_val_1.fq.gz
R2_TRIM=${OUTPUT_DIR}/trimmed/${sample}_2_val_2.fq.gz
SAM=${OUTPUT_DIR}/aligned/${sample}.sam
BAM_SORTED=${OUTPUT_DIR}/aligned/${sample}.sorted.bam
BAM_Q30=${OUTPUT_DIR}/bam_filtered/${sample}.q30.bam
BAM_NAME_SORT=${OUTPUT_DIR}/bam_filtered/${sample}.name_sorted.bam
BAM_FIXMATE=${OUTPUT_DIR}/bam_filtered/${sample}.fixmate.bam
BAM_DEDUP=${OUTPUT_DIR}/bam_filtered/${sample}.rmdup.bam
BAM_CLEAN=${OUTPUT_DIR}/bam_filtered/${sample}.clean.bam

echo "Now ${sample} is aligning !"
bowtie2 -x ${BOWTIE2_INDEX} -1 ${R1_TRIM} -2 ${R2_TRIM} \
-X 2000 --very-sensitive -p ${THREADS} --mm -S ${SAM}

samtools view -@${THREADS} -bS ${SAM} | samtools sort -@${THREADS} -o ${BAM_SORTED}
rm ${SAM}

samtools view -@${THREADS} -b -q 30 ${BAM_SORTED} -o ${BAM_Q30}
samtools sort -n -@${THREADS} -o ${BAM_NAME_SORT} ${BAM_Q30}
samtools fixmate -@${THREADS} -m ${BAM_NAME_SORT} ${BAM_FIXMATE}
samtools sort -@${THREADS} -o ${BAM_Q30} ${BAM_FIXMATE}
samtools markdup -@${THREADS} -r ${BAM_Q30} ${BAM_DEDUP}

rm ${BAM_NAME_SORT} ${BAM_FIXMATE} ${BAM_Q30}
samtools index -@${THREADS} "${BAM_DEDUP}"

samtools idxstats "${BAM_DEDUP}" | cut -f1 | grep -v -E "${FILTER_CHRS}" | \
xargs samtools view -@${THREADS} -b "${BAM_DEDUP}" > "${BAM_CLEAN}"

samtools index -@${THREADS} "${BAM_CLEAN}"
rm "${BAM_DEDUP}" "${BAM_DEDUP}.bai"
done

# === Peak calling with MACS2 ===
for sample in ${SAMPLE_LIST}; do
BAM=${OUTPUT_DIR}/bam_filtered/${sample}.clean.bam
macs2 callpeak -t ${BAM} -f BAMPE -g ${GENOME_SIZE} \
--outdir ${OUTPUT_DIR}/peaks -n ${sample} \
-B --call-summits
done

# === Peak annotation with HOMER ===
for sample in ${SAMPLE_LIST}; do
PEAK=${OUTPUT_DIR}/peaks/${sample}_peaks.narrowPeak
annotatePeaks.pl ${PEAK} hg38 > ${OUTPUT_DIR}/annotation/${sample}.annotation.txt
done

# === NEW: Merge peaks across samples (using pipeline_step4_merge_counts.sh logic) ===
PEAK_DIR="${OUTPUT_DIR}/peaks"
ALIGN_DIR="${OUTPUT_DIR}/bam_filtered"
OUT_DIR="${OUTPUT_DIR}/merged_peaks"
mkdir -p "${OUT_DIR}"

echo "===== [1/3] Build consensus peak set ====="
cat "${PEAK_DIR}"/*_peaks.narrowPeak \
| awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' \
| sort -k1,1 -k2,2n \
| bedtools merge -i - \
> "${OUT_DIR}/consensus_peaks.bed"
echo "Consensus peaks: $(wc -l < "${OUT_DIR}/consensus_peaks.bed")"

echo "===== [2/3] Compute read counts per peak (bedtools multicov) ====="
BAM_LIST=()
while IFS= read -r -d '' bam; do
BAM_LIST+=("$bam")
done < <(find "${ALIGN_DIR}" -name "*.clean.bam" -print0 | sort -z)

if [ "${#BAM_LIST[@]}" -eq 0 ]; then
echo "ERROR: No .clean.bam files found!" >&2
exit 1
fi

bedtools multicov -bed "${OUT_DIR}/consensus_peaks.bed" -bams "${BAM_LIST[@]}" \
> "${OUT_DIR}/peak_counts.bed"

echo "===== [3/3] Create count matrix TSV with header ====="
{
  printf "peak_id"
  for bam in "${BAM_LIST[@]}"; do
  s=$(basename "$bam" .clean.bam)
  printf "\t%s" "$s"
  done
  printf "\n"
} > "${OUT_DIR}/peak_counts.tsv"

# body£ºpeak_id(chr_start_end) + counts
awk 'BEGIN{OFS="\t"} {peak=$1"_"$2"_"$3; printf "%s", peak; for(i=4;i<=NF;i++) printf "\t%s",$i; printf "\n"}' \
"${OUT_DIR}/peak_counts.bed" >> "${OUT_DIR}/peak_counts.tsv"

echo "Done!"
echo " - ${OUT_DIR}/consensus_peaks.bed"
echo " - ${OUT_DIR}/peak_counts.bed"
echo " - ${OUT_DIR}/peak_counts.tsv"

# === Annotate consensus peaks ===
CONSENSUS_PEAK="${OUTPUT_DIR}/merged_peaks/consensus_peaks.bed"
CONSENSUS_ANNOT="${OUTPUT_DIR}/merged_peaks/consensus_peaks.annotation.txt"
annotatePeaks.pl "${CONSENSUS_PEAK}" hg38 > "${CONSENSUS_ANNOT}"

# === TSS enrichment heatmap ===
HOMER_TSS_FILE="/home/icac_yl/anaconda3/envs/ATAC/share/homer/data/genomes/hg38/tss.txt"
awk '
    BEGIN { FS = "\t"; OFS = "\t" }
    NF == 0 || /^#/ { next }
    NF >= 3 && $1 ~ /^chr[0-9XYM]+$/ && ($2 == "+" || $2 == "-") && $3 ~ /^[0-9]+$/ {
        chr = $1
        strand = $2
        tss_1based = $3 + 0
        tss_0based = tss_1based - 1
        if (strand == "+") {
            start = tss_0based
            end = tss_0based + 1
        } else {
            start = tss_0based - 1
            end = tss_0based
            if (start < 0) start = 0
        }
        name = sprintf("TSS_%s_%d_%s", chr, tss_1based, strand)
        score = 0
        print chr, start, end, name, score, strand
    }' "${HOMER_TSS_FILE}" | \
sort -k1,1V -k2,2n | \
tr -d '\r' > "${OUTPUT_DIR}/tss_heatmap/hg38_TSS_point.bed"

for sample in ${SAMPLE_LIST}; do
BAM="${OUTPUT_DIR}/bam_filtered/${sample}.clean.bam"
BIGWIG="${OUTPUT_DIR}/bigwig/${sample}.bw"
if [ ! -f "${BIGWIG}" ]; then
bamCoverage -b "${BAM}" -o "${BIGWIG}" \
--normalizeUsing RPKM \
--binSize 10 \
--numberOfProcessors ${THREADS}
fi
done

computeMatrix reference-point \
--referencePoint TSS \
-b 1000 -a 1000 \
-R "${OUTPUT_DIR}/tss_heatmap/hg38_TSS_point.bed" \
-S ${OUTPUT_DIR}/bigwig/*.bw \
--numberOfProcessors ${THREADS} \
-o ${OUTPUT_DIR}/tss_heatmap/tss_matrix.mat.gz

plotProfile \
-m ${OUTPUT_DIR}/tss_heatmap/tss_matrix.mat.gz \
-out ${OUTPUT_DIR}/tss_heatmap/TSS_profile.png \
--plotTitle "Average ATAC-seq Signal at TSS (1kb)" \
--legendLocation "upper-right"

plotHeatmap \
-m ${OUTPUT_DIR}/tss_heatmap/tss_matrix.mat.gz \
-out ${OUTPUT_DIR}/tss_heatmap/TSS_heatmap.png \
--colorMap viridis \
--plotTitle "ATAC-seq Signal at TSS (1kb, All Samples)" \
--zMin 0 \
--whatToShow 'heatmap and colorbar'
done