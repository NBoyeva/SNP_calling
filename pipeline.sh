#!/bin/bash


n_cores=4 # by default
temp_dir=temp_folder
VCF_dir=vcf_folder 


######## Arguments #######

while getopts 'i:o:n:d:f:b:g:s:' flag; do
  case "${flag}" in
    i) INPUT_DIR=${OPTARG} ;;
    o) OUTPUT_DIR=${OPTARG} ;;
    n) n_cores=${OPTARG} ;;
    d) REF_DIR=${OPTARG} ;;
    f) ref_fa=${OPTARG} ;;
    b) index_base=${OPTARG} ;;
    g) PATH_TO_GATK_PYFILE=${OPTARG} ;;
    s) PATH_TO_SNPEFF_FOLDER=${OPTARG} ;; 
    *) print_usage
       exit 1 ;;
  esac
done 

TEMP_DIR=$OUTPUT_DIR/$temp_dir
VCF_DIR=$OUTPUT_DIR/$VCF_dir

# Arbitrary decision about the maximal length of sequences to be trimmed of

length=140

# Make OUTPUT_DIR if it doesn't exist

mkdir -p $OUTPUT_DIR
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

for sample_number in "${unique_sorted_numbers[@]}" # Iterate over samples
do
	echo "Analyzing sample $sample_number..."
	
	# Create temporary directories
	
	mkdir -p $TEMP_DIR
	
	SUBDIR_SAMPLE=$TEMP_DIR/S${sample_number}_subdir
	mkdir -p $SUBDIR_SAMPLE
	
	for file_R1 in $INPUT_DIR/*S${sample_number}_*_R1_001.fastq.gz # Iterate over files
	do

		echo "Working with file $file_R1..."
		
		base=$(basename $file_R1 _R1_001.fastq.gz) # Basename, e.g. 48_S1_L001
		# Make output directory for the lane
		
		SUBDIR_LANE=$TEMP_DIR/$base 
		mkdir -p $SUBDIR_LANE
	
		# File with reverse reads
		
		file_R2=$INPUT_DIR/${base}_R2_001.fastq.gz
		
		# Paths for trimmed files
		
		fastq_trim1=$SUBDIR_LANE/${base}_R1_001.fastq
		fastq_trim2=$SUBDIR_LANE/${base}_R2_001.fastq

		# Trimming

		echo "Launching fastp..."
		
		fastp -i $file_R1 -I $file_R2 -l $length -y 50 -n 0 -o $fastq_trim1 -O $fastq_trim2

		if ! [ -f  $REF_DIR/${index_base}.1.bt2 ]; then #Check if index file exists			
	    		bowtie2-build $REF_DIR/$ref_fa $REF_DIR/$index_base #Make index files 
		fi

		# Path for alignment
		
		sam_output=${base}.sam
		
		# Alignment
		
		echo "Making alignment..."
		
		bowtie2 -p $n_cores \
			-x $REF_DIR/$index_base \
			-1 $fastq_trim1 \
			-2 $fastq_trim2 \
			-S $SUBDIR_LANE/$sam_output

		
		echo "Launching GATK..."
		
		# sort .sam files
		
		echo "Sorting alignments..."
		
		python3 $PATH_TO_GATK_PYFILE SortSam \
			-I $SUBDIR_LANE/$sam_output \
			-O $SUBDIR_LANE/${base}.sorted.sam \
			--SORT_ORDER coordinate \
			--VALIDATION_STRINGENCY LENIENT

		# .sam file deduplication
		
		echo "Marking duplicates..."
		
		python3 $PATH_TO_GATK_PYFILE MarkDuplicates \
			-I $SUBDIR_LANE/${base}.sorted.sam \
		       	-O $SUBDIR_LANE/${base}.dedup.sam \
			--METRICS_FILE $SUBDIR_LANE/metrix.txt \
			--ASSUME_SORTED true \
			--VALIDATION_STRINGENCY LENIENT
				
		# Convert .sam to .bam
		
		echo "Converting to .bam..."
		
		samtools sort $SUBDIR_LANE/${base}.dedup.sam > $SUBDIR_LANE/${base}.bam
		
		# Add name for read group in .bam file
		
		python3 $PATH_TO_GATK_PYFILE AddOrReplaceReadGroups \
			I=$SUBDIR_LANE/${base}.bam \
			O=$SUBDIR_SAMPLE/${base}.grouped.bam \
			RGID=4 \
			RGLB=lib1 \
			RGPL=ILLUMINA \
			RGPU=unit1 \
			RGSM=20
	
	done
	
	# Merge .bam files from each lane of the sample
	
	echo "Merging .bam files..."

	samtools merge -f $SUBDIR_SAMPLE/S${sample_number}_merged.bam $SUBDIR_SAMPLE/*_S${sample_number}_L00[1-4].grouped.bam
	
	# Index .bam files
	
	echo "Indexing merged .bam file..."
	
	samtools index $SUBDIR_SAMPLE/S${sample_number}_merged.bam
	
	# Use HaplotypeCaller to create .vcf files
	
	echo "Create .vcf files by HaplotypeCaller"
	
	python3 $PATH_TO_GATK_PYFILE --java-options "-Xmx11g" HaplotypeCaller  \
			-R $REF_DIR/$ref_fa \
			-I $SUBDIR_SAMPLE/S${sample_number}_merged.bam \
			-O $SUBDIR_SAMPLE/S${sample_number}.vcf.gz \
			-ERC GVCF 

	# SnpEff annotation
	
	echo "Perform SnpEff annotation..."
	
	java -Xmx8g -jar $PATH_TO_SNPEFF_FOLDER/snpEff.jar -v \
	     -lof \
	     GRCh37.75 \
	     $SUBDIR_SAMPLE/S${sample_number}.vcf.gz > $SUBDIR_SAMPLE/S${sample_number}.ann.vcf 
	     
	# vcf quality filtration
	
	echo "Filter .vcf file by quality..."
	
	cat $SUBDIR_SAMPLE/S${sample_number}.ann.vcf | java -jar $PATH_TO_SNPEFF_FOLDER/SnpSift.jar filter " ( QUAL >= 50 )" > $SUBDIR_SAMPLE/S${sample_number}.ann.filtered.vcf
	
	# .vcf annotation with clinvar database
	
	echo "Perform .vcf file annotation with clinvar database..."
	
	java -Xmx11g -jar $PATH_TO_SNPEFF_FOLDER/SnpSift.jar annotate -id \
	-v $PATH_TO_SNPEFF_FOLDER/data/GRCh37.75/clinvar/clinvar.vcf.gz \
	$SUBDIR_SAMPLE/S${sample_number}.ann.filtered.vcf > $SUBDIR_SAMPLE/S${sample_number}.clinvar.ann.filtered.vcf

	# Filtration by clinvar database
	
	echo "Filter .vcf files basing on clinvar data..."
	
	java -Xmx11g -jar $PATH_TO_SNPEFF_FOLDER/SnpSift.jar filter \
	-v " ( (ANN[0].IMPACT has 'HIGH') | (ANN[0].IMPACT has 'MODERATE') | (exists CLNSIGINCL) ) " \
	$SUBDIR_SAMPLE/S${sample_number}.clinvar.ann.filtered.vcf > $SUBDIR_SAMPLE/S${sample_number}.filtered_clinvar.ann.filtered.vcf

	echo "WTF..."
	# 
	java -Xmx11g -jar $PATH_TO_SNPEFF_FOLDER/SnpSift.jar varType -v \
	$SUBDIR_SAMPLE/S${sample_number}.filtered_clinvar.ann.filtered.vcf > $VCF_DIR/S${sample_number}.vcf
	
	echo "Obtain .tsv table from .vcf file..."

	# Obtain .tsv table from vcf files
	java -Xmx8g -jar $PATH_TO_SNPEFF_FOLDER/SnpSift.jar extractFields -s "," \
	-v $VCF_DIR/S${sample_number}.vcf \
	CHROM POS ID REF ALT QUAL VARTYPE SNP MNP INS DEL MIXED HOM HET ANN[*].ALLELE ANN[*].EFFECT ANN[*].IMPACT ANN[*].GENE ANN[*].GENEID ANN[*].FEATURE ANN[*].FEATUREID ANN[*].BIOTYPE ANN[*].RANK ANN[*].HGVS_C ANN[*].HGVS_P ANN[*].CDNA_POS ANN[*].CDNA_LEN ANN[*].CDS_POS ANN[*].CDS_LEN ANN[*].AA_POS ANN[*].AA_LEN ANN[*].DISTANCE ANN[*].ERRORS LOF[*].NUMTR LOF[*].PERC NMD[*].NUMTR NMD[*].PERC ANN DBVARID ALLELEID CLNSIG CLNVCSO CLNREVSTAT RS CLNDNINCL ORIGIN MC CLNDN CLNVC CLNVI AF_EXAC AF_ESP CLNSIGINCL CLNDISDB GENEINFO CLNDISDBINCL AF_TGP CLNSIGCONF CLNHGVS \
	> $VCF_DIR/S${sample_number}.tsv
	
	# Remove subfolders and files in temporal directory
	
	echo "Done with sample $sample_number"	
	rm -r $TEMP_DIR
done


