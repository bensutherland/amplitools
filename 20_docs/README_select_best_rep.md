# Select the best replicate    
## Input is in fastq format
To determine the deepest sequenced replicate, determine the number of reads in each fastq file:
`01_scripts/count_reads.sh`
...note: adjust variable SUFFIX if your files do not have the suffix '.fastq.gz'.

Prepare an interpretation file, `00_archive/filename_to_sample_map.txt`, with column headers `filename` (fastq name per file, no path) and `sample_id` (name of sample, identical for replicates).
To create an interpretation file from scratch, use the following:
`basename -a 12_input_mhap/*.fastq.gz > 00_archive/filename_to_sample_map.txt`
 ...then use spreadsheet editor to complete the info.

In RStudio, source amplitools initiator, then run the following:
```
select_best_rep_fastq(input_folder = "12_input_mhap", metadata.FN = "00_archive/filename_to_sample_map.txt", counts.FN = "12_input_mhap/reads_per_sample_table.txt")
```
This will output a histogram of reads per sample, and a reads per sample table with only the best sample retained ('best' is based on number of reads).


Then in terminal, run the following to move any file that was not set to be retained into the removed files folder:
```
for file in $(cat 12_input_mhap/remove_files.txt); do mv "$file" 12_input_mhap/removed_files/; done
```

Finally, add any manually-selected files that you would like to keep in the analysis into the `12_input_mhap` folder. If these files are paired-end, put the R2 files into the `removed_files` subfolder, and only use the R1 file to match the amplicon panel output data.



