---
output: html_document
---

------------------------------------------------------------------------

## Process BAM files obtained from cellranger-atac tool

**Inputs** (Required):
1. -open_procBAM yes
2. -script_dir (a script directory downloaded from GitHub: Socrates2/utils/input_other_required_scripts_dir)
3. -BAM_fl (a bam file built from cellranger-atac output)
4.  -o output_dir (mkdir output_dir)

**Parameters:**
**-core 1** number of cores used for the running (Default: 1)
**-mpq 10** a mapping score for filtering out multiple mapped reads (Default: 10)

**Outputs:**
1. opt_tn5_mq10.bed

**Running:**
python Socrates2/utils/prepare_tn5_gene_accessibiliy.py \
-BAM_fl /Path/to/bam_fl.bam \
-core 1 \
-mpq 10 \
-o output_dir