###Trimmomatic, HumanN2, and PICRUSt2 script###
#Performed using the University of Connecticut's HPC cluster 'Xanadu'

###Trimmomatic###

#!/bin/bash

date +"%T"

#Load trimmomatic
module load Trimmomatic/0.36

#Verify trimmomatic
module list

#Run trimmomatic: change initial directory to the appropriate directory containing the RNA and DNA files
for f in $(ls $file_directoy/.fastq.gz | sed 's/?.fastq.gz//' | sort -u)
do
java -jar $trimmomatic_file_directory/trimmomatic-0.36.jar \
 PE\
 -threads 35 \
 -phred33\
 ${f}\
 ${f}\
 ${f}_paired.fq.gz\
 ${f}_paired.fq.gz\
 ${f}_unpaired.fq.gz\
 ${f}_unpaired.fq.gz\
 LEADING:3\
 TRAILING:3\
 SLIDINGWINDOW:4:20\
 ILLUMINACLIP:/TruSeq3-PE.fa:2:30:10 \
 MINLEN:35
done

#Make new directory for the data splitting them by paired and unpaired 
mkdir $paired_directory/paired/
mkdir $unpaired_directory/unpaired/

#Move trimmed files to their appropriate directories, use paired reads for future analyses 
mv $trimmomatic_output_directory/*_unpaired.fq.gz $unpaired_directory/unpaired/
mv $trimmomatic_output_diretory/*_paired.fq.gz $paired_directory/paired/

date +"%T"


###HumanN2###

#Load Necssary Modules
module load humann2/0.11.1
module load bowtie2/2.3.3.1
module load gnu-parallel/20160622


parallel -j 1 --eta \
'humann2 \
        --threads 30 \
        --memory-use maximum \
        --input-format  fastq \
        --input {} \
        --output $HumanN2_output_directory/ \
        --nucleotide-database $chocophlan_database_directory/ \
        --protein-database $uniref50_database_directory \
        --metaphlan-option "--mpa_pkl $MetaPhlAn_database_directory/ --bowtie2db $MetaPhlAn_directory/db_v20/"' \
        ::: $directory_containing_concatenated_merged_DNA_or_RNA_files/ZR-* () #Change to either "ZR-*" or "ZD-*" for RNA and DNA samples respectively

#HumanN2 output files contain uniprot ID's and MetaCyc pathways. They can be regrouped as desribed in the HumanN2 tutorial from the Huttenhower Lab: https://huttenhower.sph.harvard.edu/humann2/


###PICRUSt2###
#PICRUSt2 was exucted using the basic pipline: available at https://github.com/picrust/picrust2/wiki/Full-pipeline-script

picrust2_pipeline.py -s $path_to_fna_file -i $path_to_biom_file -o $picrust2_output_directory/ -p 1
