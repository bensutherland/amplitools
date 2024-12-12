#!/bin/bash

# Global variables
DATAFOLDER="14_extract_mhap"
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
LOG_FOLDER="10_log_files"
SCRIPT=$0
NAME=$(basename $0)


# User set variables
NUM_CPU=36
INPUT_BCF="mpileup_calls.bcf"
INDEL_DIST=5      # indel surrounding window size around which to remove SNPs
MISSING_PPN=0.15   # missing data proportion (default: filter if > 10%)
QUAL=99           # minimum quality for SNP at least one indiv
#MIN_TOTAL_DP=10   # minimum total depth for site
MIN_AVG_DP=10     # minimum depth when averaged across all samples
MIN_GENO_DP=10     # minimum geno depth or change to missing
MAX_GENO_DP=10000   # maximum geno depth or change to missing
LOW_GQ_FILT=20

# Copy script as it was run
cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

# Store output in log file
RESULT_LOG=$LOG_FOLDER/"$TIMESTAMP"_filter_bcf.log
rm $RESULT_LOG 2> /dev/null


## Filter
# 1. Filter near indels
echo "Filtering out variants within $INDEL_DIST bp of an indel" | tee -a $RESULT_LOG
OUTPUT_FN_1="${INPUT_BCF%.bcf}"_noindel"$INDEL_DIST".bcf
bcftools filter --threads $NUM_CPU -g $INDEL_DIST $DATAFOLDER/$INPUT_BCF \
        -Ob -o $DATAFOLDER/$OUTPUT_FN_1
echo "Completed, output as $DATAFOLDER/$OUTPUT_FN_1" | tee -a $RESULT_LOG

# Reporting
echo "Output BCF now has the following number of variants: " | tee -a $RESULT_LOG
bcftools view $DATAFOLDER/$OUTPUT_FN_1 --threads $NUM_CPU | grep -vE '^#' - | wc -l | tee -a $RESULT_LOG


# 2. Filter based on missing
echo "Filtering out variants missing in more than $MISSING_PPN (proportion) individuals" | tee -a $RESULT_LOG
OUTPUT_FN_2="${OUTPUT_FN_1%.bcf}"_miss"$MISSING_PPN".bcf
bcftools view -i "F_missing < $MISSING_PPN" --threads $NUM_CPU $DATAFOLDER/$OUTPUT_FN_1 \
        -Ob -o $DATAFOLDER/$OUTPUT_FN_2
echo "Completed, output as $DATAFOLDER/$OUTPUT_FN_2" | tee -a $RESULT_LOG

# Reporting
echo "Output BCF now has the following number of variants: " | tee -a $RESULT_LOG
bcftools view $DATAFOLDER/$OUTPUT_FN_2 --threads $NUM_CPU | grep -vE '^#' - | wc -l | tee -a $RESULT_LOG


# 3. Filter based on type
echo "Keeping only variants of type 'SNP'" | tee -a $RESULT_LOG
OUTPUT_FN_3="${OUTPUT_FN_2%.bcf}"_SNP.bcf
bcftools view -i 'TYPE="snp"' --threads $NUM_CPU $DATAFOLDER/$OUTPUT_FN_2 \
        -Ob -o $DATAFOLDER/$OUTPUT_FN_3
echo "Completed, output as $DATAFOLDER/$OUTPUT_FN_3" | tee -a $RESULT_LOG

# Reporting
echo "Output BCF now has the following number of variants: " | tee -a $RESULT_LOG
bcftools view $DATAFOLDER/$OUTPUT_FN_3 --threads $NUM_CPU | grep -vE '^#' - | wc -l | tee -a $RESULT_LOG


# 4. Filter based on quality
echo "Keeping only variants with quality > $QUAL in at least one indiv" | tee -a $RESULT_LOG
OUTPUT_FN_4="${OUTPUT_FN_3%.bcf}"_q$QUAL.bcf
bcftools view -i "QUAL>=$QUAL" --threads $NUM_CPU $DATAFOLDER/$OUTPUT_FN_3 \
        -Ob -o $DATAFOLDER/$OUTPUT_FN_4
echo "Completed, output as $DATAFOLDER/$OUTPUT_FN_4" | tee -a $RESULT_LOG

# Reporting
echo "Output BCF now has the following number of variants: " | tee -a $RESULT_LOG
bcftools view $DATAFOLDER/$OUTPUT_FN_4 --threads $NUM_CPU | grep -vE '^#' - | wc -l | tee -a $RESULT_LOG


# 5. Filter based on average depth
echo "Keeping only variants with average depth across samples > $MIN_AVG_DP" | tee -a $RESULT_LOG
OUTPUT_FN_5="${OUTPUT_FN_4%.bcf}"_avgDP"$MIN_AVG_DP".bcf
bcftools view -i "AVG(FMT/DP)>$MIN_AVG_DP" --threads $NUM_CPU $DATAFOLDER/$OUTPUT_FN_4 \
        -Ob -o $DATAFOLDER/$OUTPUT_FN_5
echo "Completed, output as $DATAFOLDER/$OUTPUT_FN_5" | tee -a $RESULT_LOG

# Reporting
echo "Output BCF now has the following number of variants: " | tee -a $RESULT_LOG
bcftools view $DATAFOLDER/$OUTPUT_FN_5 --threads $NUM_CPU | grep -vE '^#' - | wc -l | tee -a $RESULT_LOG


# 6. Filter based on biallelic only
echo "Keeping only biallelic SNPs" | tee -a $RESULT_LOG
OUTPUT_FN_6="${OUTPUT_FN_5%.bcf}"_biallele.bcf
bcftools view --min-alleles 2 --max-alleles 2 --threads $NUM_CPU $DATAFOLDER/$OUTPUT_FN_5 \
        -Ob -o $DATAFOLDER/$OUTPUT_FN_6
echo "Completed, output as $DATAFOLDER/$OUTPUT_FN_6" | tee -a $RESULT_LOG

# Reporting
echo "Output BCF now has the following number of variants: " | tee -a $RESULT_LOG
bcftools view $DATAFOLDER/$OUTPUT_FN_6 --threads $NUM_CPU | grep -vE '^#' - | wc -l | tee -a $RESULT_LOG


# 7. Remove low or high coverage genotypes
echo "Set to missing low (<$MIN_GENO_DP) or high (>$MAX_GENO_DP) coverage genotypes" | tee -a $RESULT_LOG
OUTPUT_FN_7="${OUTPUT_FN_6%.bcf}"_minDP"$MIN_GENO_DP"_maxDP"$MAX_GENO_DP".bcf
bcftools filter -S . -e "(FORMAT/AD[*:0] + FORMAT/AD[*:1]) < $MIN_GENO_DP | (FORMAT/AD[*:0] + FORMAT/AD[*:1]) > $MAX_GENO_DP" $DATAFOLDER/$OUTPUT_FN_6 -Ob -o $DATAFOLDER/$OUTPUT_FN_7 --threads $NUM_CPU
echo "Completed, output as $DATAFOLDER/$OUTPUT_FN_7" | tee -a $RESULT_LOG

# Reporting
echo "Output BCF now has the following number of variants: " | tee -a $RESULT_LOG
bcftools view $DATAFOLDER/$OUTPUT_FN_7 --threads $NUM_CPU | grep -vE '^#' - | wc -l | tee -a $RESULT_LOG


# 8. Remove low GQ genotypes
echo "Set to missing low quality (GQ < $LOW_GQ_FILT) genotypes" | tee -a $RESULT_LOG
OUTPUT_FN_8="${OUTPUT_FN_7%.bcf}"_minGQ"$LOW_GQ_FILT".bcf
bcftools +setGT $DATAFOLDER/$OUTPUT_FN_7 -- -t q -n . -i "FMT/GQ<$LOW_GQ_FILT" > $DATAFOLDER/$OUTPUT_FN_8
echo "Completed, output as $DATAFOLDER/$OUTPUT_FN_8" | tee -a $RESULT_LOG

# Reporting
echo "Output BCF now has the following number of variants: " | tee -a $RESULT_LOG
bcftools view $DATAFOLDER/$OUTPUT_FN_8 --threads $NUM_CPU | grep -vE '^#' - | wc -l | tee -a $RESULT_LOG


# 9. Filter again for missing data after the previous
echo "Filtering out variants missing in more than $MISSING_PPN (proportion) individuals, following previous substitution" | tee -a $RESULT_LOG
OUTPUT_FN_9="${OUTPUT_FN_8%.bcf}"_miss$MISSING_PPN.bcf
bcftools view -i "F_missing < $MISSING_PPN" --threads $NUM_CPU $DATAFOLDER/$OUTPUT_FN_8 \
        -Ob -o $DATAFOLDER/$OUTPUT_FN_9
echo "Completed, output as $DATAFOLDER/$OUTPUT_FN_9" | tee -a $RESULT_LOG

# Reporting
echo "Output BCF now has the following number of variants: "
bcftools view $DATAFOLDER/$OUTPUT_FN_9 --threads $NUM_CPU | grep -vE '^#' - | wc -l | tee -a $RESULT_LOG

echo "Filtering complete" | tee -a $RESULT_LOG
