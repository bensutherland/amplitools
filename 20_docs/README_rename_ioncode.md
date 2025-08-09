If need to rename files to add a custom filename, for example if you have multiple projects that all use the same IonCodes, you can use the following simplified approach for renaming:    

note: backup data first, use at own risk (mv can be destructive!)    


In this case, the custom project ID is `G1124-06-VIUN`
```
cd 12_input_mhap/ 

ls -1 *.fastq.gz | awk '{ print "mv " $0 " G1124-06-VIUN_" $0 }' - > ../01_scripts/rename_fastq_add_G1124-06-VIUN.sh

vi ../01_scripts/rename_fastq_add_G1124-06-VIUN.sh # add shebang

# will probably need to make .sh script executable using chmod

# run script from within the 12_input_mhap folder
./../01_scripts/rename_fastq_add_G1124-06-VIUN.sh
```

