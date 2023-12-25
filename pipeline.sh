#!/bin/bash


n_cores=4 # by default
SAM_DIR=sam_folder
BAM_DIR=bam_folder
VCF_NAME=vcf_folder 


######## Arguments #######

while getopts 'i:o:n:d:f:b:g:' flag; do
  case "${flag}" in
    i) INPUT_DIR=${OPTARG} ;;
    o) OUTPUT_DIR=${OPTARG} ;;
    n) n_cores=${OPTARG} ;;
    d) REF_DIR=${OPTARG} ;;
    f) ref_fa=${OPTARG} ;;
    b) index_base=${OPTARG} ;;
    g) PATH_TO_GATK_PYFILE=${OPTARG} ;; 
    *) print_usage
       exit 1 ;;
  esac
done 

SAM_DIR=$OUTPUT_DIR/$SAM_DIR
BAM_DIR=$OUTPUT_DIR/$BAM_DIR
VCF_DIR=$OUTPUT_DIR/$VCF_NAME


# Arbitrary decision about the maximal length of sequences to be trimmed of
length=140

# Make OUTPUT_DIR if it doesn't exist
mkdir -p $OUTPUT_DIR
mkdir -p $BAM_DIR
mkdir -p $VCF_DIR

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

echo "Running samtools indexing..."

samtools faidx $REF_DIR/$ref_fa -o $REF_DIR/$index_base.fasta.fai
samtools dict $REF_DIR/$ref_fa -o $REF_DIR/$index_base.dict

echo "Analyzing samples..."	

######## Run pipeline by samples ########

for n in "${unique_sorted_numbers[@]}" # Iterate over samples
do
	echo "Analyzing sample $n..."
	
	for file_R1 in $INPUT_DIR/*S${n}_*_R1_001.fastq.gz # Iterate over files
	do
		
		echo "Working with file $file_R1..."
		
		base=$(basename $file_R1 _R1_001.fastq.gz) # Basename, e.g. 48_S1_L001

		# Make output directory for the lane
		OUTPUT_SUBDIR=$OUTPUT_DIR/$base 
		mkdir -p $OUTPUT_SUBDIR
	
		# File with reverse reads
		file_R2=$INPUT_DIR/${base}_R2_001.fastq.gz
		
		# Paths for trimmed files
		fastq_trim1=$OUTPUT_SUBDIR/${base}_R1_001.fastq
		fastq_trim2=$OUTPUT_SUBDIR/${base}_R2_001.fastq

		# Trimming
		echo "Launching fastp..."
		
		fastp -i $file_R1 -I $file_R2 -l $length -y 50 -n 0 -o $fastq_trim1 -O $fastq_trim2

		if ! [ -f  $REF_DIR/${index_base}.1.bt2 ]; then #Check if index file exists			
	    		bowtie2-build $REF_DIR/$ref_fa $index_base #Make index files 
		fi

		# Path for alignment
		sam_output=${base}.sam
		
		# Alignment
		echo "Making alignment..."
		
		bowtie2 -p $n_cores \
			-x $REF_DIR/$index_base \
			-1 $fastq_trim1 \
			-2 $fastq_trim2 \
			-S $OUTPUT_SUBDIR/$sam_output
		
		echo "Launching GATK..."
		
		# sort .sam files
		
		echo "Sorting alignments..."
		
		python3 $PATH_TO_GATK_PYFILE SortSam \
			-I $OUTPUT_SUBDIR/$sam_output \
			-O $OUTPUT_SUBDIR/${base}.sorted.sam \
			--SORT_ORDER coordinate \
			--VALIDATION_STRINGENCY LENIENT

		# .sam file deduplication
		
		echo "Marking duplicates..."
		
		python3 $PATH_TO_GATK_PYFILE MarkDuplicates \
			-I $OUTPUT_SUBDIR/${base}.sorted.sam \
		       	-O $OUTPUT_SUBDIR/${base}.dedup.sam \
			--METRICS_FILE $OUTPUT_SUBDIR/metrix.txt \
			--ASSUME_SORTED true \
			--VALIDATION_STRINGENCY LENIENT
				
		# Convert .sam to .bam
		
		echo "Converting to .bam..."
		
		samtools sort $OUTPUT_SUBDIR/${base}.dedup.sam > $OUTPUT_SUBDIR/${base}.bam
		
		# Add name for read group .bam file
		python3 $PATH_TO_GATK_PYFILE AddOrReplaceReadGroups \
			I=$OUTPUT_SUBDIR/${base}.bam \
			O=$BAM_DIR/${base}.grouped.bam \
			RGID=4 \
			RGLB=lib1 \
			RGPL=ILLUMINA \
			RGPU=unit1 \
			RGSM=20
		
	done
	
	samtools merge -f $BAM_DIR/S${n}_merged.bam $BAM_DIR/*_S${n}_L00[1-4].grouped.bam

	samtools index $BAM_DIR/S${n}_merged.bam

	python3 $PATH_TO_GATK_PYFILE --java-options "-Xmx11g" HaplotypeCaller  \
			-R $REF_DIR/$ref_fa \
			-I $BAM_DIR/S${n}_merged.bam \
			-O $VCF_DIR/${n}_sample.vcf.gz \
			-ERC GVCF
	
	for item in "$OUTPUT_DIR"/*; do
  		if [[ -d $item && $(basename "$item") != $VCF_NAME ]]; then
    
    			rm -f "$OUTPUT_DIR/$item"/*
  		fi
	done
done





