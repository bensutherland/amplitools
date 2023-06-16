#!/bin/bash
# Use tab-delim file with marker names to extract required info from a VCF
# Note: will depend on the VCF mnames being the same as the selected mnames
# Expects that your marker names file has the first column as mname, in format of
#  mname_pos_nucleotide; e.g., 171077_68_G

# User-set variables
MNAMES_FN=$(basename $1)
INPUT_VCF=$(basename $2)

# Global variables
INPUT_FOLDER="02_input_data"
OUTPUT_FOLDER="10_designer"
OUTPUT="vcf_selection.csv"

# Remove any existing output
echo "Removing $OUTPUT to make new file" 
rm $OUTPUT_FOLDER/$OUTPUT

# Report on presence of duplicate markers
echo "*** If present, duplicate markers reported here: ***"
awk '{ print $1 }' $INPUT_FOLDER/$MNAMES_FN | 
    grep -vE '^mname' |
    sort | uniq -c | sort -n | 
    awk '$1 >1 { print  }' - 

# Report on the number of records to collect
echo "*** Expect the output to have this many markers: ***"
awk '{ print $1 }' $INPUT_FOLDER/$MNAMES_FN | 
    grep -vE '^mname' |
    sort | uniq | 
    wc -l  

# Prepare bed file
cat $INPUT_FOLDER/$MNAMES_FN | 
    
    # Keep only the mname and remove duplicates
    awk '{ print $1 }' - | 
    grep -vE '^mname' | 
    sort |
    uniq |

    # Remove positional and variant information from the mname. 
    # Add colon to match VCF mname. Add variant position, plus 1 due to 0-based count.   
    # note: assumes that the positional data does not need +1 due to 0- vs. 1-based count
    awk -F_ '{ print $1 ":" $2 ":"}' - | 

    # For each mname, extract the relevant section of the vcf
    while read i
    do
        # Reporting
        echo $i

        # Extract the necessary information from the VCF, matching i after the tab
        grep -P "\t"$i $INPUT_FOLDER/$INPUT_VCF |
        awk '{ print $1 "," $2 "," $3 "," $4 "," $5 }' - >> $OUTPUT_FOLDER/$OUTPUT
 
    done    

echo "Completed, the relevant sections of the VCF are in $OUTPUT_FOLDER/$OUTPUT"

echo "Removing any records that have fewer than 201 bp upstream of the variant"     
awk -F, '( $2 > 201 )' $OUTPUT_FOLDER/$OUTPUT > $OUTPUT_FOLDER/temp.txt 
mv $OUTPUT_FOLDER/temp.txt $OUTPUT_FOLDER/$OUTPUT 

echo "The output file now has the following number of records available" 
wc -l $OUTPUT_FOLDER/$OUTPUT 

echo "Process finished."
