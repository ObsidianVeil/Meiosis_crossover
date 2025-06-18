#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300g
#SBATCH --time=72:00:00
#SBATCH --job-name=Allfiles_all_analysis
#SBATCH --output=/share/BioinfMSc/ml_projects/FH_files/logs/%x-%j.out
#SBATCH --error=/share/BioinfMSc/ml_projects/FH_files/logs/%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbxsh12@nottingham.ac.uk

source $HOME/.bash_profile

#to do:
#make one conda env with bigwig, LAST and python in
#Python script needs a for loop to process each chromosome #should be done now? 18/06/25



conda activate LASTaligner

cd /share/BioinfMSc/ml_projects/FH_files

mkdir -p barcodes train_files maf_files dotplots sam_files bam_files bamfiles_sorted coverageplots bigwigs brokenbarplots

scp /share/BioinfMSc/ml_projects/Yeasts/SGray/*.pass.fastq.gz ./barcodes

samtools faidx edited_chr_refs.fasta


for barcode in barcodes/*.pass.fastq.gz; do
	filename=$(basename "$barcode" | cut -c1-9)
	echo "Starting lastdb on $barcode"
	lastdb -P32 refdb edited_chr_refs.fasta
	
	echo "Starting last-train on $barcode"
	last-train -P32 -Q0 refdb $barcode > train_files/"$filename".train
	
	echo "Starting lastal for $barcode"
	lastal -P32 -Q0 --split refdb $barcode | last-split > maf_files/"$filename"_lastsplit.maf
	
	echo "Making dotplot for $barcode"
	last-dotplot dotplots/"$filename"_lastsplit.maf > dotplots/"$filename".png #make dotplot
	
	echo "converting $barcode maf file to sam file"
	#convert maf to sam file
	maf-convert sam maf_files/"$filename"_lastsplit.maf > sam_files/"$filename".sam

	#convert sam file to bam file
	echo "Converting $barcode sam file to bam file"
	samtools view -bt edited_chr_refs.fasta.fai -o bam_files/"$filename".bam sam_files/"$filename".sam

	#make bai index to go with bam file
	samtools sort -o bamfiles_sorted/"$filename"_sorted.bam bam_files/"$filename".bam
	samtools index bamfiles_sorted/"$filename"_sorted.bam
done

conda deactivate

conda activate bigwig

for barcode in bamfiles_sorted/*.bam; do
	[[ "$barcode" == *.bam.bai ]] && continue #skips the .bam.bai files
	filename=$(basename "barcode" | cut -c1-9) #creates filename with just barcodeXX
	echo "Generating bigwig file for $barcode"
	bamCoverage -b "$barcode" -o bigwigs/"$filename"_coverage.bw -of bigwig
	
	echo "Calculating coverage for $barcode"
	plotCoverage -b "$barcode" -o coverageplots/"$filename"_coverage.png --numberOfProcessors max --verbose
done

conda deactivate

conda activate lectures
for file in bigwigs/*.bw; do
	filename=$(basename "$file")
	fullpath="$file"
	barcode="${filename%%_coverage.bw}"
	echo "$fullpath"
	export barcode
	export fullpath
	mkdir -p brokenbarplots/"$barcode"
	python ./placeholderpath
	
done

conda deactivate
exit