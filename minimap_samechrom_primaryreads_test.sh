#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300g
#SBATCH --time=48:00:00
#SBATCH --job-name=test_samechrom_primary
#SBATCH --output=/share/BioinfMSc/ml_projects/FH_files/logs/%x-%j.out
#SBATCH --error=/share/BioinfMSc/ml_projects/FH_files/logs/%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbxsh12@nottingham.ac.uk

source $HOME/.bash_profile

#################################################################################################################################################

output=/share/BioinfMSc/ml_projects/FH_files/output_samechrom_primary #set to wherever you want all your files

input=/share/BioinfMSc/ml_projects/FH_files/barcodes_small #set to the source directory of all the files

ref_input=/share/BioinfMSc/ml_projects/Yeasts/edited_chr_refs.fasta #set to location of reference fasta file

#################################################################################################################################################
#Define scripts

conda activate meiosis

#cross_chrom_python="Cross_chrom_mapping_correct.py"

#Crossover_csv_script="Crossover_points.py"

matplotlib_script="/share/BioinfMSc/ml_projects/FH_files/matplotlibgraphs_with_transition_csv.py"

################################################################################################################################################
#Define directories
bigwigfolder="bigwigs"

pyexpdir="interchrom_BAM_files"

train="train_output"

splitfiles="maf_files"

bamfiles1="bam_files"

sortbamdir="bamfiles_sorted"

sortbamdir2="bamfiles_filtered_sorted"

headerbamdir="bamfilewheader"

graphoutput="mergedbarplots"

csvoutput="csvfiles"

################################################################################################################################################
#Define functions

echoerr() { echo "$@" 1>&2; }

cleanup () {
	echo "Cleaning up"
	echo "Removing: ${removedirs[*]}"
	for item in "${removedirs[@]}"; do
		if [[ -n "$item" && -e "$item" ]]; then
			echo "Deleting $item"
			rm -r -- "$item"
		else
			echo "Skipping $item (does not exist)"
		fi
	done
}

endrun () {
	echo "Ending run"
	cleanup
	exit 1
}

checkfile () {
	if [ ! -f "$1" ]; then
		echoerr "$1 does not exist"
		endrun
	else
		echo "$1 does exist"
	fi
}

mapreadcount () {
	mapped_reads=$(samtools view -c -F 4 "$1")
	if [ "$mapped_reads" -eq 0 ]; then
		echoerr "No mapped reads in $1."
	fi
	
	total_reads=$(samtools view -c "$1")
	echo "$1: $mapped_reads mapped out of $total_reads total reads"
}

################################################################################################################################################

mkdir -p "$output"

#cp "$cross_chrom_python" "$output"

cp "$matplotlib_script" "$output"

if [ -f "$Crossover_csv_script" ]; then
	scp "$Crossover_csv_script" "$output"
fi

cd "$output" 

mkdir -p barcodes "$train" "$splitfiles"  "$bamfiles1" "$sortbamdir" "$bigwigfolder" "$graphoutput" "$pyexpdir" "$headerbamdir" "$sortbamdir2" "$csvoutput" tmp

removedirs=(barcodes "$train" "$splitfiles" "$bamfiles1" "$sortbamdir" "$ref_fai" "$Crossover_csv_script" tmp "$sortbamdir2" "$headerbamdir" "$pyexpdir")

cp -r "$input"/* ./barcodes

cp "$ref_input" ./                                  

reference=$(basename "$ref_input")

samtools faidx "$reference"

ref_fai="$reference".fai

#checkfile "$cross_chrom_python"

checkfile "$ref_input"

checkfile "$ref_fai"

for barcode in barcodes/*; do
	filename=$(basename "$barcode" | cut -c1-9)
	alignoutput="${aligndir}/${filename}"
	samtoolssort_output="${sortbamdir}/${filename}_sorted.bam"
	pyoutput1="${pyexpdir}/${filename}_interchrom_filtered.bam"
	tmp="tmp/${filename}_tmp.bam"
	
	checkfile "$barcode"
	
	#echo "Starting minimap for $filename"
	#minimap2 -ax map-ont "$reference" "$barcode" > "$alignoutput"
	
	#checkfile "$alignoutput"
	
	echo "Running minimap on $filename, followed by removing secondary reads with samtools view, and then sorting by name with samtools sort"
	minimap2 -t 3 -ax map-ont "$reference" "$barcode" | samtools view -@ "$SLURM_CPUS_PER_TASK" -b -h -F 260  | samtools sort -n -@ "$SLURM_CPUS_PER_TASK" -o "$samtoolssort_output" #Does this keep supplementary as well as primary while cutting secondary?
	
	checkfile "$samtoolssort_output"

	suppnum=$(samtools view -@ 32 -c -f 2048 "$samtoolssort_output")  # supplementary
	secnum=$(samtools view -@ 32 -c -f 256 "$samtoolssort_output")  # secondary
	primnum=$(samtools view -@ 32 -c -F 2308 "$samtoolssort_output")  # primary (neither secondary nor supplementary)
	echo "Primary reads: $primnum"
	echo "Secondary reads: $secnum"
	echo "supplementary reads: $suppnum"
	
	if (($suppnum + $secnum + $primnum <1 )); then
		echoerr "No reads in $filename"
		continue
	fi
		
	#mapreadcount "$samtoolssort_output"	
done


for bamfile in "$sortbamdir"/*.bam; do
	filename=$(basename "$bamfile" | cut -c1-9)
	echo "filename is $filename"
	bamoutput="${headerbamdir}/$filename.bam"
	sortedoutput="${sortbamdir2}/${filename}_sorted.bam"
	bigwigoutput="${bigwigfolder}/${filename}_coverage.bw"
	
	checkfile "$bamfile"
	
	mapreadcount "$bamfile"
	
	#Add header to file
	echo "Adding header to $filename bam file using $ref_fai"
	samtools view -bh -@ "$SLURM_CPUS_PER_TASK" -t "$ref_fai" -o "$bamoutput" "$bamfile"
	
	checkfile "$bamoutput"

	#make bai index for $filename to go with bam file
	echo "Sorting $filename and making bai index file"
	samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$sortedoutput" "$bamoutput"
	samtools index "$sortedoutput"
	
	checkfile "$sortedoutput"
	
	echo "Generating bigwig file for $filename"
	echo "bigwig file is $bigwigoutput"
	bamCoverage -b "$sortedoutput" -o "$bigwigoutput" -of bigwig
	
	checkfile "$bigwigoutput"
	
	#echo "Plotting coverage for $filename"
	#plotCoverage -b "$sortedoutput" -o "${coveragedir}/${filename}_coverage.png" --numberOfProcessors max
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
cleanup

echo "Finished. Deactivating Conda and ending run"

conda deactivate
exit