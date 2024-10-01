### Retain only hotspot (target) SNPs ###
In some cases, it will be necessary to remove all but hotspot SNPs, for example, if using TVC VCF file output. The following assumes the data is in a single BCF file, and the marker name is retained for hotspot SNPs. In this example, the input VCF has been converted by SNPlift.      

Remove novel variants from BCF file:     
```
# Identify hotspots (i.e., field 3 of VCF file has a marker name, and is not '.')
bcftools view 03_prepped_data/all_sample_renamed_snplift_rehead.bcf | grep -vE '^#' - | awk '$3 != "." {print $1 "\t" $2 }' - > 03_prepped_data/include_snps.txt

# Keep only the include SNPs in VCF
bcftools view --targets-file ./03_prepped_data/include_snps.txt 03_prepped_data/all_sample_renamed_snplift_rehead.bcf -Ov -o 03_prepped_data/all_sample_renamed_snplift_rehead_hotspot.vcf

```

Next: the hotspot-only VCF file is ready for analysis.   

Back to [main README](https://github.com/bensutherland/amplitools/tree/main).

