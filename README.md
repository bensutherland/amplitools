# amplitools
Tools for working with amplicon sequencing data

Put your input torrent data in folder `02_input_data`.      
Open `proton_to_genepop.R`, and set the user-set variable to point to the input torrent data.     
Choose if you want to only include hotspot SNPs, or also include novel SNPs, using the true or false variable `hotspot_only`.        

This script will:      
1. Load the data, reformat; 
2. Convert from Allele Call to the actual alleles; 

Note: in variantCaller data,        
'Absent' means homozygous reference       
'No Call' means missing data         
'Heterozygous' means heterozygous        
'Homozygous' means homozygous variant         

3. Convert from Allele Call to genepop 0101 format
4. Create a block of data that can be used to format a genepop, currently `genetic_data_only_final.txt`     


Back to bash:    
```
head -n 1 genetic_data_only_final.txt | sed 's/\t/\n/g' > header.txt 
echo "POP" > pop.txt
tail -n+2 genetic_data_only_final.txt > tail.txt
sed 's/\t/,\t/' ./tail.txt > tail_w_comma.txt
cat header.txt pop.txt tail_w_comma.txt > my_data.gen
```

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

