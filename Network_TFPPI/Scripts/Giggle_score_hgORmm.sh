#!/usr/bin/env bash

# Usage: bash classify_and_process.sh input_file annotation.gtf
# input_file = peaks.bed or genes.txt
# annotation.gtf = gene annotation file (e.g. from Ensembl)

file="$1"
species=$2
if [[ $species == "mouse" ]]; then
  ref_tss="GenCodeVM25.gene.tss100K.filtered.bed" ### overlap with gene, get the genelist TSS region
  ref_gene="Gencode.vM25.gene.sorted.bed"  ### get the closest gene for the input bed file###
  giggle_index="/workspace/rsrch1/tmp/DataBase_datatable0618/Network_TFGiggle1010/202501_CirTrom/Top1k_Giggle_index" ## need change
else
   ref_tss="GRCh38.94.tss100kb.bed"
   ref_gene="GRCh38.94.gene.sorted.bed"
   giggle_index="/workspace/rsrch1/tmp/DataBase_datatable0618/Network_TFGiggle1010/202501_CirTrom/Top1k_hg38Giggle_index" ## need change
fi

# Get first non-empty, non-comment line
cols=$(awk 'NF>0{print NF; exit}' "$file")

# Convert GTF to BED-like gene file (chr, start, end, gene_name)
gene_bed="genes_from_gtf.bed"


if [[ $cols -ge 3 ]]; then
    echo "🔹 Detected BED file ($file) — finding nearest genes..."
    type="bed"
    bedtools closest -a ${file} -b ${ref_gene} -t first | awk '{print $NF}' | sort | uniq | grep -v "\\." > geneList.txt
    grep '^[cC]' ${file} | sort -k1,1 -k2,2n > peakfile.bed
    echo "✅ Output: ${file%.bed}_nearest_genes.bed"
else
    echo "🔹 Detected gene list ($file) — generating ±100 kb gene region..."
    # Extract TSS coordinates from GTF for each gene in list
    echo "Generating TSS ±100kb for listed genes..."
    grep -f ${file} ${ref_tss} > peakfile.bed
    grep -v '^\s' ${file} > geneList.txt
    type="genes"
fi


###### go giggle analysis 

######as giggle input must be gz file, so need to bgzip -c $i > $i.gz
bgzip -f peakfile.bed >peakfile.bed.gz
/workspace/rsrch1/tmp/DataBase_datatable0618/Network_TFGiggle1010/giggle/bin/giggle search -i ${giggle_index} -q peakfile.bed.gz -s > peakfile_GIGGLE_res_Top1K.xls

/usr/bin/Rscript /workspace/rsrch1/tmp/DataBase_datatable0618/Network_TFGiggle1010/Scripts/Giggle_score.network_hg38ORmm.R /workspace/rsrch1/tmp/DataBase_datatable0618/Network_TFGiggle1010/hg38_mm10_TestInput ${species} ${type}_network.pdf


