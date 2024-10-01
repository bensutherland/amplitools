### Convert VCF file positions to different genome ###
It might be necessary to convert the positions of variants in a VCF file to a different genome, in which case the program SNPlift provides a useful solution.     

Brief instructions are given here, but please see full instructions in the source repo: [SNPlift](https://github.com/enormandeau/snplift).    

Prepare inputs   
```
# Move out of amplitools into the parent directory, then clone snplift
cd ..
git clone https://github.com/enormandeau/snplift.git

# Copy the target VCF into SNPlift
cp ./amplitools/03_prepped_data/all_sample_renamed.vcf ./snplift/04_input_vcf/
```

Prepare the provided default SNPlift config file (`02_infos/snplift_config.sh`) by editing the following parameters:     
- full path to the original decompressed genome (bwa indexed)
- full path to the target decompressed genome (bwa indexed)
- relative path to the original VCF filename
- relative path to the new VCF filename
- `CORRECT_ALLELES=1` to convert the ref/alt alleles when alignments are reverse complemented
Note: setting the `CORRECT_ID` to 0 above prevents the ID column from being recalculated, so that your original IDs are carried through to the new VCF.

Run SNPlift and fix output
```
# Run using config
./snplift 02_infos/snplift_config.sh

# Add header (adjust w/ fai-indexed genome location)
bcftools reheader all_sample_renamed_snplift.vcf --fai ~/Documents/genomes/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna.fai --output ./all_sample_renamed_snplift_rehead.vcf

# Compress
bgzip all_sample_renamed_snplift_rehead.vcf

# Convert to BCF file
bcftools view ./all_sample_renamed_snplift_rehead.vcf.gz -Ob -o all_sample_renamed_snplift_rehead.bcf

# Index
bcftools index all_sample_renamed_snplift_rehead.bcf

# Save output
cp all_sample_renamed_snplift_rehead.bcf* ./../amplitools/03_prepped_data/
```

Additional resources including excluding novel SNPs are linked from the [main README](https://github.com/bensutherland/amplitools/tree/main).

