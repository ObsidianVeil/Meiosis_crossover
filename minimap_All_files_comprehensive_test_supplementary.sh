#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300g
#SBATCH --time=48:00:00
#SBATCH --job-name=test_interchrom_8_files
#SBATCH --output=/share/BioinfMSc/ml_projects/FH_files/logs/%x-%j.out
#SBATCH --error=/share/BioinfMSc/ml_projects/FH_files/logs/%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbxsh12@nottingham.ac.uk

source $HOME/.bash_profile

#################################################################################################################################################

output=/home/mbxsh12/ind_proj/minimaptest_supplementary #set to wherever you want all your files

input=/home/mbxsh12/ind_proj/barcodes_small #set to the source directory of all the files

ref_input=/home/mbxsh12/ind_proj/edited_chr_refs.fasta #set to location of reference fasta file

#################################################################################################################################################

conda activate meiosis

mkdir -p "$output"

cross_chrom_python="Cross_chrom_mapping_correct.py"

Crossover_csv_script="Crossover_points.py"

matplotlib_script="Matplotlibgraphs.py"

retag_script="Retag_supplements.py"

################################################################################################################################################

bigwigfolder="bigwigs"

pyexpdir="interchrom_BAM_files"

train="train_output"

splitfiles="maf_files"

bamfiles1="bam_files"

sortbamdir="bamfiles_sorted"

sortbamdir2="bamfiles_filtered_sorted"

headerbamdir="bamfilewheader"

coveragedir="coverageplots"

graphoutput="mergedbarplots"

csvoutput="csvfiles"

aligndir="minimap_files"

################################################################################################################################################

checkfile () {
	if [ ! -f "$1" ]; then
		echo "$1 does not exist"
		exit 1
	else
		echo "$1 does exist"
	fi
}

mapreadcount () {
	mapped_reads=$(samtools view -c -F 4 "$1")
	if [ "$mapped_reads" -eq 0 ]; then
		echo "No mapped reads in $1."
		
	fi
	
	total_reads=$(samtools view -c "$1")
	echo "$1: $mapped_reads mapped out of $total_reads total reads"
}

cp "$cross_chrom_python" "$output"

cp "$matplotlib_script" "$output"

if [ -f "$Crossover_csv_script" ]; then
	scp "$Crossover_csv_script" "$output"
fi

cd "$output" 

mkdir -p barcodes "$train" "$splitfiles" dotplots "$bamfiles1" "$sortbamdir" "$coveragedir" "$bigwigfolder" "$graphoutput" "$pyexpdir" "$headerbamdir" "$sortbamdir2" "$csvoutput" tmp "$aligndir"

cp -r "$input"/* ./barcodes

cp "$ref_input" ./                                  

reference=$(basename "$ref_input")

samtools faidx "$reference"

ref_fai="$reference".fai

checkfile "$cross_chrom_python"

checkfile "$ref_input"

checkfile "$ref_fai"

for barcode in barcodes/*; do
	filename=$(basename "$barcode" | cut -c1-9)
	alignoutput="${aligndir}/${filename}"
	samtoolssort_output="${sortbamdir}/${filename}_sorted.bam"
	pyoutput1="${pyexpdir}/${filename}_interchrom_filtered.bam"
	tmp="tmp/${filename}_tmp.bam"
	
	checkfile "$barcode"
	
	#echo "Starting last-train on $filename"                  #does this really need running?
	#last-train -P32 -Q0 refdb "$barcode" > "$train_output"
	
	#checkfile "$train_output"
	
	echo "Starting minimap for $filename"
	minimap2 -ax map-ont "$reference" "$barcode" > "$alignoutput"
	
	checkfile "$alignoutput"
	
	echo "Sorting $filename"
	samtools sort "$alignoutput" -o "$samtoolssort_output"
	samtools index "$samtoolssort_output"
	
	checkfile "$samtoolssort_output"
	
	#readlength=$(samtools view "$samtoolssort_output" | awk '{total += length($10); count++} END {print total/count}')

	#echo "$samtoolssort_output has an average read length of $readlength"
	
	#samtools view "$samtoolssort_output" | cut -f1,2,3 | grep -E 'SA:|2048|256'

	samtools view -c -f 2048 "$samtoolssort_output"  # supplementary
	samtools view -c -f 256 "$samtoolssort_output"  # secondary
	samtools view -c -F 256 -F 2048 "$samtoolssort_output"  # primary (neither secondary nor supplementary)

	
	#mapreadcount "$samtoolssort_output"
	
	if [ -f "$cross_chrom_python" ]; then
		echo "Beginning Cross Chromosomal mapping for $filename"
		echo "exporting $samtoolssort_output"
		export samtoolssort_output
		echo "exporting $pyoutput1"
		export pyoutput1
		python "$cross_chrom_python"
	
		checkfile "$pyoutput1"
	
		echo "Sorting BAM file after Python filtering"
		tmp="${pyoutput1}_temp_sorted.bam"
		samtools sort "$pyoutput1" -o "$tmp"
		mv "$tmp" "$pyoutput1"
	else
		echo "Chromosomal crossover filtering disabled"
	fi
	
	#temporarily hashed to test above if statement
	
	#echo "Beginning Cross Chromosomal mapping for $filename"
	#echo "exporting $samtoolssort_output"
	#export samtoolssort_output
	#echo "exporting $pyoutput1"
	#export pyoutput1
	#python "$cross_chrom_python"
	
	#checkfile "$pyoutput1"
	
	#echo "Sorting BAM file after Python filtering"
	#tmp="${pyoutput1}_temp_sorted.bam"
	#samtools sort "$pyoutput1" -o "$tmp"
	#mv "$tmp" "$pyoutput1"
	
done

#silenced for testing
#if cross chromosomal filtering doesn't occur, redirects to sorted bam file.
#if [ -z "$( ls -A '$pyoutput1' )" ]; then
#   "$pyexpdir"="$sortbamdir"
#fi


for bamfile in "$pyexpdir"/*.bam; do
	filename=$(basename "$bamfile" | cut -c1-9)
	echo "filename is $filename"
	Pysamoutput="${pyexpdir}/${filename}_interchrom_filtered.bam"
	bamoutput="${headerbamdir}/$filename.bam"
	sortedoutput="${sortbamdir2}/${filename}_sorted.bam"
	bigwigoutput="${bigwigfolder}/${filename}_coverage.bw"
	
	checkfile "$Pysamoutput"
	
	mapreadcount "$Pysamoutput"
	
	mapped_reads=$(samtools view -c -F 4 "$bamfile")
	if [ "$mapped_reads" -eq 0 ]; then
		echo "No mapped reads in $bamfile."
		continue
	fi
	
	echo "Indexing $bamfile"
	samtools index "$bamfile"
	
	#Add header to file
	echo "Adding header to $filename bam file using $ref_fai"
	samtools view -bt  "$ref_fai" -o "$bamoutput" "$Pysamoutput"
	
	checkfile "$bamoutput"

	#make bai index for $filename to go with bam file
	echo "Sorting $filename and making bai index file"
	samtools sort -o "$sortedoutput" "$bamoutput"
	samtools index "$sortedoutput"
	
	checkfile "$sortedoutput"
	
	echo "Generating bigwig file for $filename"
	echo "bigwig file is $bigwigoutput"
	bamCoverage -b "$sortedoutput" -o "$bigwigoutput" -of bigwig
	
	checkfile "$bigwigoutput"
	
	echo "Calculating coverage for $filename"
	plotCoverage -b "$sortedoutput" -o "${coveragedir}/${filename}_coverage.png" --numberOfProcessors max
done

if [ -f "$matplotlib_script" ]; then
	echo "Beginning python graphing run"
	export bigwigfolder
	export graphoutput
	python "$matplotlib_script"
else
	echo "Matplotlib Python script does not exist. Graphs will not be generated"
	rmdir "$graphoutput"
fi

if [ -f "$Crossover_csv_script" ]; then
	echo "Generating CSV files of crossover points"
	export bigwigfolder
	export csvoutput
	python "$Crossover_csv_script"
else
	echo "Crossover script does not exist"
	rmdir "$csvoutput"
fi

#silence for debugging
echo "Cleaning up"
rm -r barcodes "$train" "$splitfiles" "$bamfiles1" "$sortbamdir" "$ref_fai" "$Crossover_csv_script" tmp "$sortbamdir2" "$headerbamdir" "$pyexpdir" "$aligndir"

echo "Finished. Deactivating Conda and ending run"

conda deactivate
exit