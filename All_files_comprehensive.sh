#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300g
#SBATCH --time=48:00:00
#SBATCH --job-name=Allfiles_all_analysis
#SBATCH --output=/share/BioinfMSc/ml_projects/FH_files/logs/%x-%j.out
#SBATCH --error=/share/BioinfMSc/ml_projects/FH_files/logs/%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbxsh12@nottingham.ac.uk

source $HOME/.bash_profile


conda activate meiosis

cd /share/BioinfMSc/ml_projects/FH_files #set to wherever you want all your files

mkdir -p barcodes train_files maf_files dotplots sam_files bam_files bamfiles_sorted coverageplots bigwigs mergedbarplots

scp /share/BioinfMSc/ml_projects/Yeasts/SGray/*.pass.fastq.gz ./barcodes #amend source directory to wherever you saved your fastq files

samtools faidx edited_chr_refs.fasta


for barcode in barcodes/*.pass.fastq.gz; do
	filename=$(basename "$barcode" | cut -c1-9)
	echo "Starting lastdb on $barcode"
	lastdb -P32 refdb edited_chr_refs.fasta #insert reference genome here
	
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

for barcode in bamfiles_sorted/*.bam; do
	filename=$(basename "$barcode" | cut -c1-9) #creates filename with just barcodeXX
	echo "Generating bigwig file for $barcode"
	bamCoverage -b "$barcode" -o bigwigs/"$filename"_coverage.bw -of bigwig
	
	echo "Calculating coverage for $barcode"
	plotCoverage -b "$barcode" -o coverageplots/"$filename"_coverage.png --numberOfProcessors max
done

python Matplotlibgraphs.py

conda deactivate
exit