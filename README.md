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

