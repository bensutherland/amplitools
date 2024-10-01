### Combining VCF files produced by panel ###
In some cases, the amplicon panel results may be provided as Torrent VariantCaller VCF files.      
To load this into amplitools, or other applications, it may be necessary to merge the VCF files into a single file.     

##### Merge multiple VCF files #####
```
# Copy all VCF files in 02_input_data

# If the files are compressed, decompress
gunzip 02_input_data/*.vcf.gz

# Compress with bgzip
ls 02_input_data/*.vcf | xargs -n 1 bgzip

# Index the files with bcftools
ls 02_input_data/*.vcf.gz | xargs -n 1 bcftools index

# Create filelist for merging VCF files
ls -1 02_input_data/*.vcf.gz > 02_input_data/sample_list.txt

# Merge VCF files
bcftools merge --file-list 02_input_data/sample_list.txt -Ov -o 03_prepped_data/all_sample.vcf
```

##### Rename samples #####
At this stage, you can rename your sample names as follows:    
```
# Identify samples present in the VCF file
bcftools query -l 03_prepped_data/all_sample.vcf > 03_prepped_data/sample_list_for_rename.txt

# Manually edit the above text file as a space separated file with two columns, adding the second column with desired renamed

# Rename samples
bcftools reheader --samples 03_prepped_data/sample_list_for_rename.txt -o 03_prepped_data/all_sample_renamed.vcf 03_prepped_data/all_sample.vcf
```

Additional resources including (1) converting to a different reference genome; or (2) excluding novel SNPs; are linked from the [main README](https://github.com/bensutherland/amplitools/tree/main).         

