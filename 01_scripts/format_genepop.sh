#!/bin/bash
# Add the header to a genepop

# Usage: 
#   ./01_scripts/format_genepop.sh 03_results/your_file_genetic_data_only_final.txt

# Uses the output of the function proton_to_genepop() 
# Has not been generalized yet
# NOTE: the input filename must not have any special characters or spaces

# The input file is in the following format:
# Must have the 'indiv' as first row 
# indiv marker1 marker2 marker3
# F2-03 0101    0000    0102
# etc. 

# The output file is to be in the following format:
# FILENAME
# marker1
# marker2
# marker3
# POP
# F2-03, 0101   0000    0102 

# Global variables
INPUT=$1

# Reporting
echo "formatting $INPUT"

# Create pathless FILENAME (see above)
shortname=$(basename $INPUT)

# Create output filename
OUT_FN=$(echo $INPUT | sed 's/.txt/.gen/g')

# Create temp header (i.e., filename in top row followed by each marker name per row
head -n 1 $INPUT | sed 's/\t/\n/g' | 
       
        sed "s/indiv/$shortname/g" > header_temp.txt

# Create pop line
echo "POP" > pop_temp.txt

# Create genetic data line
tail -n+2 $INPUT |  sed 's/\t/,\t/' > tail_temp.txt

# Assemble all pieces
cat header_temp.txt pop_temp.txt tail_temp.txt > $OUT_FN

# Remove temp files
rm header_temp.txt pop_temp.txt tail_temp.txt

