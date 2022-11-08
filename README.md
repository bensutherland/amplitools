# amplitools
Tools for working with amplicon sequencing data.       
Currently designed for S5XL variantCaller output data.     

## 00. Setup ##
Put your input torrent data in folder `02_input_data`.      

## 01. Conversion from proton to genepop (R steps) ##
Open `01_scripts/proton_to_genepop.R`.        
Set the variable `proton.FN` to the first input torrent data filename.      
Set the variable `hotspot_only` true or false to only include hotspot SNPs (exclude novel)        

This script will:      
**Load the data, reformat**         
Although many columns will be present (~51), the first step will reduce this to only 14 columns:      
"Chrom", "Position", "Ref", "Variant", "Allele.Call", "Type", "Allele.Source", "Allele.Name", "Region.Name", "Coverage", "Strand.Bias", "Sample.Name", "Barcode", "Run.Name"        

Based on user input, novel SNPs will be dropped, keeping hotspot SNPs only.      
A new column, 'identifier' will be created, comprised of `<RunName>__<Barcode>__<SampleName>`.      
e.g., `R_2022_08_04_09_19_56_user_S5XL-00533-1089-OYR-20220729_7__IonCode_0501__F2-03`        


**Per sample, per marker, convert from Allele.Call to the actual genotypes (nucleotides)**        
This will also provide the column, 'genepop', which will give numeric genotypes.        

Note: in variantCaller Allele.Call data,        
'Absent' means homozygous reference (0101)       
'No Call' means missing data (0000)         
'Heterozygous' means heterozygous (0102)        
'Homozygous' means homozygous variant (0202)        

Note: this assumes that per marker, the identity of the reference and variant alleles in the line item is always the same, regardless of the specific sample.      

Your output will be written out as: `03_results/<run_name>_proton_data_converted.txt`, but this is purely for troubleshooting, this file will not be used again in the pipeline.      

**Collect data into a genetic block**         
The genepop genotypes, connected to the alleles, will be generated into a large dataframe, and written out as `03_results/<run_name>_genetic_data_only.txt`. (temp file only)         
The final genepop genotypes will be named `03_results/<run_name>_genetic_data_only_final.txt`.       

## 02. Finalization of proton to genepop (bash steps) ##
Use the following script to format the genepop genotypes to a genepop file
`01_scripts/format_genepop.sh 03_results/R_2022_08_04_09_19_56_user_S5XL-00533-1089-OYR-20220729_7_2_.xls_genetic_data_only_final.txt`        
Note: there cannot be any spaces or special characters in the name.    



## 03. Other functions ##
Align a FASTA against a reference genome to see where your markers are:       
`./01_scripts/bwa_align.sh`

`samtools view -q 30 05_genome_plot/query_seqs.sam -o 05_genome_plot/query_seqs_above_30_MAPQ.bam`      

`samtools view 05_genome_plot/query_seqs_above_30_MAPQ.bam | awk '{ print $1 }' - | sort -n | uniq -c | sort -nk1 | awk '{ print $1 }' - | uniq -c | sort -nk1 | less`       

`samtools view 05_genome_plot/query_seqs_above_30_MAPQ.bam | sed 's/:://g' | grep 'LR761' - | awk '{ print $1 "\t" $3 "\t" $4 }' - > 05_genome_plot/mnames_to_alt_genome_pos.txt`      

Note: here is where you could technically reduce your list based on markers passing filters (#TODO#).       
Get `retained_loci_keep_targets.txt` from `simple_pop_stats`.    
`awk '{ print "mname_"$1"JH"}' 02_input_data/retained_loci.txt > 02_input_data/retained_loci_keep_targets.txt`
`grep -f 02_input_data/retained_loci_keep_targets.txt 05_genome_plot/mnames_to_alt_genome_pos.txt > 05_genome_plot/mnames_to_alt_genome_pos_QCd_loci.txt`    

`awk '!seen[$1]++' 05_genome_plot/mnames_to_alt_genome_pos.txt > 05_genome_plot/mnames_to_alt_genome_pos_no_dups.txt` 
or
`awk '!seen[$1]++' 05_genome_plot/mnames_to_alt_genome_pos_QCd_loci.txt > 05_genome_plot/mnames_to_alt_genome_pos_QCd_loci_no_dups.txt`


`awk '{ print $2 }' 05_genome_plot/mnames_to_alt_genome_pos_QCd_loci_no_dups.txt | sort -n | uniq -c > 05_genome_plot/markers_no_dups_per_chr.txt`




o integrate:      
Markers included in proton file:     
`awk -F '\t' '{ print $12 }' 02_input_data/R_2022_08_04_09_19_56_user_S5XL-00533-1089-OYR-20220729_7\ \(2\).xls | grep -v 'tvc' | grep -v 'Allele' - | sort | uniq > markers_included_in_proton_file.txt`

