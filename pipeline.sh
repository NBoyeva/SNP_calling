#!/bin/bash

:'
 -i input directoty with fastq files
 -o output directory for pipeline work results
 -n number of cores
 -rd directory with reference file
 -rf reference filename
 -ib index base for reference file indexes
 
'

n_cores=4 # by default



######## Arguments #######

while getopts 'i:o:n:d:f:b:' flag; do
  case "${flag}" in
    i) INPUT_DIR=${OPTARG} ;;
    o) OUTPUT_DIR=${OPTARG} ;;
    n) n_cores=${OPTARG} ;;
    d) REF_DIR=${OPTARG} ;;
    f) ref_fa=${OPTARG} ;;
    b) index_base=${OPTARG} ;;
    *) print_usage
       exit 1 ;;
  esac
done 


# Arbitrary decision about the maximal length of sequences to be trimmed of
length=145

# Make OUTPUT_DIR if it doesn't exist
mkdir -p $OUTPUT_DIR



######## Retrieve list of sample numbers ########

# Define an empty array to store the sample numbers
sample_numbers=()

# Define the regular expression pattern to extract sample numbers
pattern="[0-9]+_S([0-9]+)_"

# Loop through the files in the folder
for file in "$INPUT_DIR"/*; do
    if [[ $file =~ $pattern ]]; then
        sample_number="${BASH_REMATCH[1]}"
        sample_numbers+=("$sample_number")
    fi
done


# Remove duplicates and sort the array naturally
unique_sorted_numbers=($(echo "${sample_numbers[@]}" | tr ' ' '\n' | sort -n -k1,1 | uniq))

######## Run pipeline by samples ########

for n in "${unique_sorted_numbers[@]}" # Iterate over samples
do
	for file_R1 in $INPUT_DIR/*S${n}_*_R1_001.fastq.gz # Iterate over files
	do
		
		base=$(basename $file_R1 _R1_001.fastq.gz) # Basename, e.g. 48_S1_L001

		# Make output directory for the lane
		OUTPUT_SUBDIR=$OUTPUT_DIR/$base 
		mkdir -p $OUTPUT_SUBDIR
	
		# File with reverse reads
		file_R2=$INPUT_DIR/${base}_R2_001.fastq.gz
		#html_file=$OUTPUT_SUBDIR/${base}.html
		
		# Paths for trimmed files
		fastq_trim1=$OUTPUT_SUBDIR/${base}_R1_001.fastq
		fastq_trim2=$OUTPUT_SUBDIR/${base}_R2_001.fastq
		
		echo $fastq_trim1 $fastq_trim2

		# Trimming
		fastp -i $file_R1 -I $file_R2 -l $length -y 50 -n 0 -o $fastq_trim1 -O $fastq_trim2

		# if ! [ -f  $REF_DIR/${ref_fa}.pac ]; then #Check if index file exists
	    		# bowtie2-build $REF_DIR/$ref_fa $index_base #Make index files (I believe It does)
		# fi
		
		# Path for alignment
		sam_output=${base}.sam
		
		# Alignment
		bowtie2 -p $n_cores \
		-x $REF_DIR/$index_base \
		-1 $fastq_trim1 \
		-2 $fastq_trim2 \
		-S $OUTPUT_SUBDIR/$sam_output
		
		# Convert .sam to .bam
		samtools sort $OUTPUT_SUBDIR/$sam_output > $OUTPUT_SUBDIR/${base}.bam
		
	done
done







