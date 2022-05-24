#!/bin/bash

set -e

###################
## Assign inputs ##
###################

# Define usage of the script
function print_usage {
  printf """
Usage: mlstpipeline.sh    [-h or --help]
                          [-f or --fastqfolder]
                          [-o or --outname]
                          [-t or --threads]
                          [-c or --conversion]

"""
}
# Describe usage of the tool and provide help
function print_help {
  print_usage
  printf """
Optional arguments:
    -h, --help:
                Show this help message and exit.
    -o, --outname:
                Name of your analysis.
                It will be used to name the output files.
                Default: mymlst.
    -t, --threads:
                Number of threads that will be used.
                It must be an integer.
                Default: 8.
    -c, --conversion:
                If the argument is included, convert the prokka.gff to prokka.gtf.
                Default: no conversion.
Required arguments:
    -f, --fastqfolder:
                Path to the folder that contains ALL your FASTQ files.
                Only FASTQ files should be placed in it.
                You need forward and reverse paired-end reads.
"""
}

# Define inputs
for ARGS in "$@"; do
  shift
        case "$ARGS" in
                "--fastqfolder") set -- "$@" "-f" ;;
                "--outname") set -- "$@" "-o" ;;
                "--threads") set -- "$@" "-t" ;;
                "--conversion") set -- "$@" "-c" ;;
                "--help") set -- "$@" "-h" ;;
                *) set - "$@" "$ARGS"
        esac
done

# Define defaults
outn="mymlst"; threads=8; conversion=0

# Define all parameters
while getopts 'f:o::t::ch' flag; do
        case "${flag}" in
                f) fastqfolder=${OPTARG} ;;
                o) outn=${OPTARG} ;;
                t) threads=${OPTARG} ;;
                c) conversion=${OPTARG} ;;
                h) print_help
                   exit 1;;
                *) print_usage
                    exit 1;;
        esac
done

##############################
## Identify Software Errors ##
##############################

printf "\nChecking if required software is installed...\n"
# Check installation of the required software.
# If something is missing, show how to install it.
if ! [ -x "$(command -v mlst)" ]; then
  echo "Missing: MLST not found"
  echo "Information on the installation:"
  echo "https://github.com/tseemann/mlst"
  exit 127
fi
if ! [ -x "$(command -v spades.py)" ]; then
  echo "Missing: SPAdes not found"
  echo "Information on the installation:"
  echo "https://github.com/ablab/spades"
  exit 127
fi
if ! [ -x "$(command -v prokka)" ]; then
  echo "Missing: Prokka not found"
  echo "Information on the installation:"
  echo "https://github.com/tseemann/prokka"
  exit 127
fi
if ! [ -x "$(command -v agat_convert_sp_gff2gtf.pl)" ]; then
  echo "Missing: agat_convert_sp_gff2gtf.pl not found"
  echo "Information on the installation:"
  echo "https://github.com/NBISweden/AGAT"
  echo "If agat is installed, but this error persists, reactivate the conda environment."
  exit 127
fi
echo "Required software is properly installed."


########################################
## Identify Errors in Inputs Required ##
########################################

printf "\nChecking if required inputs are correct...\n"
# Check if directory containing FASTQ files exist
if [ ! -d ${fastqfolder} ]; then
  echo "Error: --fastqfolder doesn't exist."
  echo "Solution: check if the path to this directory is correct."
  exit 1
fi

# Check if only FASTQ files are provided in the fastqfolder
for seqfile in ${fastqfolder}/*; do
  # Get the extension of the file, if file is compressed
  if file --mime-type "${seqfile}" | grep -q gzip; then
    filename=${seqfile%.*}
    extension=${filename##*.}
  else # Get the extension of the file, if file is NOT compressed
    extension=${seqfile##*.}
  fi
  # Check if extension is fastq or fq
  extension=$(tr "[:upper:]" "[:lower:]" <<< ${extension})
  if [[ ${extension} != "fq" && ${extension} != "fastq" ]]; then
    echo "Error: --fastqfolder should only contain FASTQ files."
    echo "Solution: remove any other file from this directory."
    exit 1
  fi
done
echo "Required inputs seem correct."


########################################
## Identify Errors in Optional Inputs ##
########################################

printf "\nChecking if optional inputs are correct...\n"
# Check if the number of threads is an integer
if ! [[ ${threads} =~ ^[0-9]+$ ]]; then
  echo "Error: --threads is not an integer."
  echo "Solution: remove this optional parameter or use an integer."
  exit 1
fi
echo "Optional inputs seem correct."


############################
## Perform Reads Assembly ##
############################

printf "\nPerforming genome assemblies...\n"
# Count number of files that are found in the folder
end=$(ls -1q ${fastqfolder} | wc -l)

for i in $(seq 2  2 ${end}); do
  # Get the 1st file of a pair
  filefirst=$(ls -1q ${fastqfolder} | head -n$i | tail -n2 | head -n1)
  # Get the 2nd file of a pair
  filesecond=$(ls -1q ${fastqfolder} | head -n$i | tail -n1)
  # Perform assembly with spades
  spades.py -1 ${fastqfolder}/${filefirst} \
    -2 ${fastqfolder}/${filesecond} \
    --careful --threads ${threads} --cov-cutoff auto -o ${outn}_${i}/assembly
  # Concatenate contigs into one long sequence that could be used as reference
  echo ">Assembly" > ${outn}_${i}/assembly/ref_assembly.fasta
  sed '/^>/d' ${outn}_${i}/assembly/contigs.fasta >> \
    ${outn}_${i}/assembly/ref_assembly.fasta
done
echo "Genome assemblies finished successfully."


#####################################
## Copy results to temporal folder ##
#####################################

# Remove temporal directory if existing and create new one with all permissions
rm -rf tmp_mlst && \
  mkdir tmp_mlst && \
  chmod +xwr tmp_mlst

for i in $(seq 2  2 ${end}); do
  # Get name of the file
  filefirst=$(ls -1q ${fastqfolder} | head -n$i | tail -n1)
  # Check if file is compressed and get file name without the extension
  if file --mime-type ${fastqfolder}/${filefirst} | grep -q gzip; then
    nameext=${filefirst%.*}
    fname=${nameext%.*}
  else
    fname=${filefirst%.*}
  fi
  # Copy the contig results to the temporal folder
  cp ${outn}_${i}/assembly/contigs.fasta tmp_mlst/${fname}.fasta
done


###########################
## Perform MLST analysis ##
###########################

printf "\nPerforming MLST analysis...\n"
mlst -t ${threads} -q tmp_mlst/* > ${outn}_MLST.tsv
echo "MLST analysis finished successfully."


###############################
## Compute Prokka annotation ##
###############################

printf "\nPerforming gene annotations...\n"
for i in $(seq 2  2 ${end}); do
  # Compute annotation
  prokka --outdir ${outn}_${i}/annotation --prefix prokka --cpus ${threads} \
    ${outn}_${i}/assembly/contigs.fasta
done
echo "Gene annotations finished successfully."


###########################
## Remove temporal files ##
###########################

rm -r tmp_mlst

######################################
## Convert prokka.gff to prokka.gtf ##
######################################
if [[ $conversion != 0 ]]; then
	printf "\nConverting gff to gtf...\n"
	for foldername in $outn_*/; do   # */ --> only directories
		agat_convert_sp_gff2gtf.pl --gff $foldername/annotation/prokka.gff -o $foldername/annotation/prokka.gtf
	done	
	echo "Conversion finished successfully."
fi

################################################
## Rename folders for an easier understanding ##
################################################

printf "\nOrdering output files...\n"
for i in $(seq 2  2 ${end}); do
  # Get name of the file
  filefirst=$(ls -1q ${fastqfolder} | head -n$i | tail -n1)
  # Check if file is compressed and get file name without the extension
  if file --mime-type ${fastqfolder}/${filefirst} | grep -q gzip; then
    nameext=${filefirst%.*}
    fname=${nameext%.*}
  else
    fname=${filefirst%.*}
  fi
  # Move the folder to a new name matching the files names
  mv ${outn}_${i} ${outn}_${fname}
done

# Remove final directory if existing and create new one with all permissions
rm -rf ${outn} && \
  mkdir ${outn} && \
  chmod +xwr ${outn}

# Move outputs inside final directory
mv ${outn}_* ${outn}
echo "Done!"
echo "All analyses finished successfully. Good luck with the results!"
