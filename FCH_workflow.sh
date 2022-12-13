#!/bin/bash
set -e
# setting to break on newline only
IFS=$'\n'

## Set project variables ##

read -r -p "What directory are the raw sequencing reads in? Please give the absolute path to 00_Raw_Reads without a trailing '/'  " reads_dir || exit 100
read -r -p "Please provide the absolute path to the project working directory:     " project_dir || exit 100
read -r -p "What is your project name? No spaces or special characters please.   " project || exit 100
read -r -p "How much RAM do you have available in gigabytes (Gb)? Use digits only.   " RAM || exit 100
read -r -p "How many cores do you want to use? Use digits only   " CPU || exit 100
read -r -p "Create a samples.txt file with unique sample ids without spaces, saved in the same directory as the raw sequencing reads. Do you have the samples.txt file made? yes or no   " samples || exit 100

if [[ "$samples" == no ]]
then
  echo Change into the 00_Raw_Reads directory and in the terminal run: "ls ./*_R1.fastq.gz | cut -f1 -d '.' | rev | cut -f2- -d '_' | rev | sort -u >> samples.txt"
fi

# Make and Change into new project directory
mkdir -p $project_dir
cd ${project_dir}
###################
## PREPROCESSING ##
###################

mkdir -p ./00_Raw_Reads/fastqc
mkdir -p ./01_HQ_Reads/merged_reads/fastqc
mkdir -p ./01_HQ_Reads/error_corrected/fastqc
mkdir -p ./01_HQ_Reads/fastqc
mkdir -p ./01_HQ_Reads/logs
mkdir -p ./01_HQ_Reads/nonpareil
mkdir -p ./01_HQ_Reads/phyloFlash

# Copy file of unique sample ids into the project directory
cp ${reads_dir}/samples.txt .

# for file in `ls ${reads_dir}/*.f*q.gz`; do
#   echo Counting all the sequences in $file
#   counts=$(gunzip -c $file | grep '^@' | wc -l)
#   echo $file,$counts >> ${reads_dir}/00_Raw_Reads_counts.csv
# done

# ##########
# # CONDA ##
# ##########
# Install bbmap, fastqc, and multiqc in a conda environment:
# conda create -y -n pp bbmap fastqc multiqc
read -r -p "Do you need to preprocess your raw sequencing reads? yes or no   " reads_QC || exit 100


#The following workflow, creates high-quality trimmed, error corrected, and merged reads.

if [[ "${reads_QC}" == yes ]]
then

  conda activate pp

  while read -r sample; do



    #Input raw reads from 00_Raw_Reads directory

    R1=( ${reads_dir}/${sample}*1.f*q.gz )
    R2=( ${reads_dir}/${sample}*2.f*q.gz )

    echo Assessing intial read quality on ${sample} raw reads in ${reads_dir}
    fastqc $R1 $R2 -q -o ./00_Raw_Reads/fastqc
    ###########################
    #### Format Conversion ####
    ###########################
    #Create interleaved file and verify paired end reads.
    reformat.sh in=$R1 in2=$R2 out=./01_HQ_Reads/iReads.fq.gz verifypaired=t ow  >> ./01_HQ_Reads/logs/${sample}.log 2>&1

    echo Starting to quality filter ${sample} raw reads

    ###########################
    #### Adapter Trimming  ####
    ###########################
    bbduk.sh in=./01_HQ_Reads/iReads.fq.gz out=./01_HQ_Reads/atrimmed.fq.gz ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo ow=t >> ./01_HQ_Reads/logs/${sample}.log 2>&1

    ###########################
    #### Quality Trimming  ####
    ###########################
    #removes Ns from sequences
    # 'entropy' means to filter out reads with low complexity
    # 'maq' is 'mininum average quality' to filter out overall poor reads
    #trimpolygright=0 will remove poly G's from the 3' end that get added erroneously because 2-dye chemistry cannot call 'undefined' bases. Undefined bases get called as "G".
    #ftm=5 removes the extra base sometimes added by the sequencer. i.e.: if a read is 151bp it will trim it to 150bp, but leave a 150bp read alone.
    #bbduk.sh in=./01_HQ_Reads/atrimmed.fq.gz out=./01_HQ_Reads/qtrimmed.fq.gz qtrim=rl trimq=6 ftm=5 trimpolygright=0 minlen=70 ordered maxns=0 maq=8 entropy=.95 ow=t >> ./01_HQ_Reads/logs/${sample}.log 2>&1

    #ftl=16 to remove the 16bp of barcode left on the 5' read. Removes Nextera 5' noise
    bbduk.sh in=./01_HQ_Reads/atrimmed.fq.gz out=./01_HQ_Reads/qtrimmed.fq.gz ftl=16 qtrim=rl trimq=6 ftm=5 trimpolygright=0 minlen=70 ordered maxns=0 maq=8 entropy=.95 ow=t >> ./01_HQ_Reads/logs/${sample}.log 2>&1

    ###########################
    ## Contaminant filtering ##
    ###########################
    #Remove synthetic artifacts and spike-ins by kmer-matching
    # 'cardinality' will generate an accurate estimation of the number of unique kmers in the dataset using the LogLog algorithm
    bbduk.sh in=./01_HQ_Reads/qtrimmed.fq.gz out=./01_HQ_Reads/${sample}_HQ.fq.gz ref=artifacts,phix k=31 hdist=1 stats=./01_HQ_Reads/logs/${sample}_contamination_stats.txt ordered cardinality ow=t >> ./01_HQ_Reads/logs/${sample}.log 2>&1

    echo finished bbduk
    # remove extra files
    rm ./01_HQ_Reads/iReads.fq.gz ./01_HQ_Reads/atrimmed.fq.gz ./01_HQ_Reads/qtrimmed.fq.gz >> ./01_HQ_Reads/logs/${sample}.log 2>&1
    #Merge reads
    echo merging reads
    bbmerge-auto.sh in=./01_HQ_Reads/${sample}_HQ.fq.gz out=./01_HQ_Reads/merged_reads/${sample}_HQ_merged.fq.gz outu=./01_HQ_Reads/merged_reads/${sample}_HQ_unmerged.fq.gz rem k=62 extend2=50 ecct vstrict ow=t >> ./01_HQ_Reads/logs/${sample}.log 2>&1
    echo starting error correction
    bbcms.sh in=./01_HQ_Reads/${sample}_HQ.fq.gz out=./01_HQ_Reads/error_corrected/${sample}_HQ_ecc.fq.gz outb=./01_HQ_Reads/error_corrected/low_${sample}.fq.gz bits=4 hashes=3 k=31 mincount=2 hcf=0.4 tossjunk=t ow=t >> ./01_HQ_Reads/logs/${sample}.log 2>&1

    #Create R1.fa for nonpareil and phyloFlash input
    reformat.sh in=./01_HQ_Reads/error_corrected/${sample}_HQ_ecc.fq.gz out=./01_HQ_Reads/error_corrected/${sample}_HQ_ecc_R1.fa out2=./01_HQ_Reads/error_corrected/${sample}_HQ_ecc_R2.fa verifypaired >> ./01_HQ_Reads/logs/${sample}.log 2>&1

    echo finished ${sample}


  done < samples.txt


  #Run multiqc on all clean Paired End reads.
  echo Assessing quality of processed ${project} reads.
  find ./01_HQ_Reads/ -name "*.fq.gz" -exec fastqc '{}' -q -o ./01_HQ_Reads/fastqc \;
  multiqc ./01_HQ_Reads/fastqc -o ./01_HQ_Reads/

  conda deactivate

  echo Finished read quality control. Starting analysis of sequencing effort.
  echo

  conda activate nonpareil

  while read -r sample; do
    echo Estimating sequencing effort for ${sample}
    nonpareil -s ./01_HQ_Reads/error_corrected/${sample}_HQ_ecc_R1.fa -T alignment -L 35 -t ${CPU} -R ${RAM}000 -b ./01_HQ_Reads/nonpareil/${sample}
  done < samples.txt

  conda deactivate

  echo Starting assembly of rRNA SSU and classification

  conda activate phyloflash

  while read -r sample; do
    cd ${project_dir}/01_HQ_Reads/phyloFlash
    echo running phyloFlash on ${sample}
    phyloFlash.pl -lib ${sample} -read1 ${project_dir}/01_HQ_Reads/error_corrected/${sample}_HQ_ecc_R1.fa -read2 ${project_dir}/01_HQ_Reads/error_corrected/${sample}_HQ_ecc_R2.fa -CPU ${CPU} -almosteverything -taxlevel 20

    cd ${project_dir}

  done < samples.txt

  conda deactivate

elif [[ "${reads_QC}" == no ]]; then :
  #statements
fi


################
### Assembly ###
################

###############
### CONDA #####
###############

# conda create -n assembly -y -c bioconda spades megahit idba
# conda activate assembly
# pip install biopython

#Create directories

mkdir -p ./02_Assembly/metaspades/logs
mkdir -p ./02_Assembly/megahit-large/logs
mkdir -p ./02_Assembly/idba_ud/logs
mkdir -p ./02_Assembly/stats


for sample in `cat samples.txt`; do

  mkdir -p ./02_Assembly/quast/${sample}

  conda activate assembly

  read -r -p "Which HQ reads dataset do you want to use to assemble $sample? HQ HQ_ecc HQ_merged Skip:  " HQ_Reads || exit 100


  #Assemble reads per sample using metaspades

## Do not make a directory for megahit, it will make it itself and will fail if the output directory already exists, as to not overwrite results. ##

  if [[ "${HQ_Reads}" == HQ ]]
  then
    echo Using metaspades to assemble ${HQ_Reads} ${sample} reads.
    mkdir -p ./02_Assembly/metaspades/${sample}
    metaspades.py --only-assembler -k 21,33,55,77,99,127 --12 ./01_HQ_Reads/${sample}_HQ.fq.gz -t ${CPU} -m ${RAM} -o ./02_Assembly/metaspades/${sample} >> ./02_Assembly/metaspades/logs/${sample}.log 2>&1
    cp ./02_Assembly/metaspades/${sample}/scaffolds.fasta ./02_Assembly/quast/${sample}/metaspades.fasta
    #Assemble reads per site using megahit
    echo Using megahit to assemble ${HQ_Reads} ${sample} reads.
#    megahit --12 ./01_HQ_Reads/${sample}_HQ.fq.gz --presets meta-large -m ${RAM} -t ${CPU} -o ~/${sample} >> ./02_Assembly/megahit-large/logs/${sample}.log 2>&1
    megahit --12 ./01_HQ_Reads/${sample}_HQ.fq.gz --k-min 21 --k-max 91 --k-step 10 -m ${RAM} -t ${CPU} -o ~/${sample} >> ./02_Assembly/megahit-large/logs/${sample}.log 2>&1
    cp  ~/${sample}/final.contigs.fa ./02_Assembly/quast/${sample}/megahit.fasta
    mkdir -p ./02_Assembly/megahit-large/${sample}
    mv ~/${sample} ./02_Assembly/megahit-large
    #IDBA can only take fasta files as input; therefore, fastq.gz files must be reformated to fasta files
    reformat.sh in=./01_HQ_Reads/${sample}_HQ.fq.gz out=./01_HQ_Reads/${sample}_HQ.fa ow=t >> ./02_Assembly/idba_ud/logs/${sample}.log 2>&1
    #fq2fa --paired ./01_HQ_Reads/${sample}_HQ.fq.gz ./01_HQ_Reads/${sample}_HQ.fa
    # #Assemble reads per site using idba_ud
    echo Using idba to assemble ${HQ_Reads} ${sample} reads.
    mkdir -p ./02_Assembly/idba_ud/${sample}
    idba_ud -r ./01_HQ_Reads/${sample}_HQ.fa --num_threads=${CPU} -o ./02_Assembly/idba_ud/${sample} >> ./02_Assembly/idba_ud/logs/${sample}.log 2>&1
    cp ./02_Assembly/idba_ud/${sample}/scaffold.fa ./02_Assembly/quast/${sample}/idba.fasta

    # --- Evaluation ---
    echo metagenome_id,reads,assembler,n_contigs,contig_bp,gap_pct,ctg_L50,ctg_max,gene_counts,mapped_pct,avg_cov,assembly_performance >> ./02_Assembly/stats/${sample}_assembly_stats.csv
    echo Evaluating ${sample} assemblies with AssemblyStats
    ##statswrapper.sh ./02_Assembly/quast/${sample}/*.fasta minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/stats/${sample}.txt
    stats.sh ./02_Assembly/metaspades/${sample}/scaffolds.fasta minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/metaspades/${sample}/stats.txt
    spades_stats=$(cat ./02_Assembly/metaspades/${sample}/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d' | sed 's/\t/,/g')
    stats.sh ./02_Assembly/megahit-large/${sample}/final.contigs.fa minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/megahit-large/${sample}/stats.txt
    megahit_stats=$(cat ./02_Assembly/megahit-large/${sample}/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d' | sed 's/\t/,/g')
    stats.sh ./02_Assembly/idba_ud/${sample}/scaffold.fa minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/idba_ud/${sample}/stats.txt
    idba_stats=$(cat ./02_Assembly/idba_ud/${sample}/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d' | sed 's/\t/,/g')
    #Evaluate assemblies with Quast (leave out "-R ref.fa if you don't have a reference)

    echo Predicting prokaryotic genes in each ${sample} assembly
    #Count gene predictions for prokaryotes using prodigal for each assembly
    python ~/bin/LengthFilter.py ./02_Assembly/metaspades/${sample}/scaffolds.fasta -m 1000 > ./02_Assembly/metaspades/${sample}/${sample}_c1k.fa
    prodigal -i ./02_Assembly/metaspades/${sample}/${sample}_c1k.fa -a ./02_Assembly/metaspades/${sample}/${sample}_c1k_proteins.faa -f gff -p meta -o ./02_Assembly/metaspades/${sample}/${sample}_c1k.gff >> ./02_Assembly/metaspades/logs/${sample}.log 2>&1
    spades_counts=$(grep '^>' ./02_Assembly/metaspades/${sample}/${sample}_c1k_proteins.faa | wc -l )

    python ~/bin/LengthFilter.py ./02_Assembly/megahit-large/${sample}/final.contigs.fa -m 1000 > ./02_Assembly/megahit-large/${sample}/${sample}_c1k.fa
    prodigal -i ./02_Assembly/megahit-large/${sample}/${sample}_c1k.fa -a ./02_Assembly/megahit-large/${sample}/${sample}_c1k_proteins.faa -f gff -p meta -o ./02_Assembly/megahit-large/${sample}/${sample}_c1k.gff >> ./02_Assembly/megahit-large/logs/${sample}.log 2>&1
    megahit_counts=$(grep '^>' ./02_Assembly/megahit-large/${sample}/${sample}_c1k_proteins.faa | wc -l )

    python ~/bin/LengthFilter.py ./02_Assembly/idba_ud/${sample}/scaffold.fa -m 1000 > ./02_Assembly/idba_ud/${sample}/${sample}_c1k.fa
    prodigal -i ./02_Assembly/idba_ud/${sample}/${sample}_c1k.fa -a ./02_Assembly/idba_ud/${sample}/${sample}_c1k_proteins.faa -f gff -p meta -o ./02_Assembly/idba_ud/${sample}/${sample}_c1k.gff >> ./02_Assembly/idba_ud/logs/${sample}.log 2>&1
    idba_counts=$(grep '^>' ./02_Assembly/idba_ud/${sample}/${sample}_c1k_proteins.faa | wc -l )

    echo Determining metagenome content within ${sample} assemblies
    # Note that these are reads mapped to contigs > 2,000 bp.
    #(i.e. % mapped) by mapping reads by to assembly
    bbmap.sh in=./01_HQ_Reads/${sample}_HQ.fq.gz ref=./02_Assembly/metaspades/${sample}/${sample}_c1k.fa nodisk covhist=./02_Assembly/metaspades/${sample}/covhist.txt covstats=./02_Assembly/metaspades/${sample}/covstats.txt out=./02_Assembly/metaspades/${sample}/${sample}_reads_assembled.fq.gz outu=./02_Assembly/metaspades/${sample}/${sample}_reads_unassembled.fq.gz fast=t ambig=best >> ./02_Assembly/metaspades/${sample}/map.log 2>&1
    spades_match=$(grep Percent\ mapped: ./02_Assembly/metaspades/${sample}/map.log | cut -f2 -d$'\t' | tail -1 | tail -1)
    spades_cov=$(grep Average\ coverage: ./02_Assembly/metaspades/${sample}/map.log | cut -f2 -d$'\t' | tail -1 | tail -1)
    spades_Pa=$(echo "$spades_match*$spades_cov" | bc)

    bbmap.sh in=./01_HQ_Reads/${sample}_HQ.fq.gz ref=./02_Assembly/megahit-large/${sample}/${sample}_c1k.fa nodisk covhist=./02_Assembly/megahit-large/${sample}/covhist.txt covstats=./02_Assembly/megahit-large/${sample}/covstats.txt outm=./02_Assembly/megahit-large/${sample}/${sample}_reads_assembled.fq.gz outu=./02_Assembly/megahit-large/${sample}/${sample}_reads_unassembled.fq.gz fast=t ambig=best >> ./02_Assembly/megahit-large/${sample}/map.log 2>&1
    megahit_match=$(grep Percent\ mapped: ./02_Assembly/megahit-large/${sample}/map.log | cut -f2 -d$'\t' | tail -1 | tail -1)
    megahit_cov=$(grep Average\ coverage: ./02_Assembly/megahit-large/${sample}/map.log | cut -f2 -d$'\t' | tail -1 | tail -1)
    megahit_Pa=$(echo "$megahit_match*$megahit_cov" | bc)

    bbmap.sh in=./01_HQ_Reads/${sample}_HQ.fq.gz ref=./02_Assembly/idba_ud/${sample}/${sample}_c1k.fa nodisk covhist=./02_Assembly/idba_ud/${sample}/covhist.txt covstats=./02_Assembly/idba_ud/${sample}/covstats.txt outm=./02_Assembly/idba_ud/${sample}/${sample}_reads_assembled.fq.gz outu=./02_Assembly/idba_ud/${sample}/${sample}_reads_unassembled.fq.gz fast=t ambig=best >> ./02_Assembly/idba_ud/${sample}/map.log 2>&1
    idba_match=$(grep Percent\ mapped: ./02_Assembly/idba_ud/${sample}/map.log | cut -f2 -d$'\t' | tail -1 | tail -1)
    idba_cov=$(grep Average\ coverage: ./02_Assembly/idba_ud/${sample}/map.log | cut -f2 -d$'\t' | tail -1 | tail -1)
    idba_Pa=$(echo "$idba_match*$idba_cov" | bc)


    echo Creating assembly comparision for ${sample}
    echo ${sample},${HQ_Reads},metaspades,${spades_stats},${spades_counts},${spades_match},${spades_cov},${spades_Pa} >> ./02_Assembly/stats/${sample}_assembly_stats.csv
    echo ${sample},${HQ_Reads},megahit,${megahit_stats},${megahit_counts},${megahit_match},${megahit_cov},${megahit_Pa} >> ./02_Assembly/stats/${sample}_assembly_stats.csv
    echo ${sample},${HQ_Reads},idba_ud,${idba_stats},${idba_counts},${idba_match},${idba_cov},${idba_Pa} >> ./02_Assembly/stats/${sample}_assembly_stats.csv

  elif [[ "${HQ_Reads}" == HQ_ecc ]]
  then
    echo Using metaspades to assemble ${HQ_Reads} ${sample} reads.
    mkdir -p ./02_Assembly/metaspades/${sample}_ecc
    metaspades.py -k 21,33,55,77,99,127 --12 ./01_HQ_Reads/error_corrected/${sample}_HQ_ecc.fq.gz -t ${CPU} -m ${RAM} -o ./02_Assembly/metaspades/${sample}_ecc >> ./02_Assembly/metaspades/logs/${sample}.log 2>&1
    cp ./02_Assembly/metaspades/${sample}_ecc/scaffolds.fasta ./02_Assembly/quast/${sample}/metaspades_ecc.fasta
    #Assemble reads per site using megahit
    echo Using megahit to assemble ${HQ_Reads} ${sample} reads.

  #  megahit --12 ./01_HQ_Reads/error_corrected/${sample}_HQ_ecc.fq.gz --presets meta-large -m ${RAM} -t ${CPU} -o ./02_Assembly/megahit-large/${sample} >> ./02_Assembly/megahit-large/logs/${sample}.log 2>&1
    megahit --12 ./01_HQ_Reads/error_corrected/${sample}_HQ_ecc.fq.gz --k-min 21 --k-max 91 --k-step 10 -m ${RAM} -t ${CPU} -o ~/${sample}_ecc >> ./02_Assembly/megahit-large/logs/${sample}.log 2>&1
    cp ~/${sample}_ecc/final.contigs.fa ./02_Assembly/quast/${sample}/megahit_ecc.fasta
    mkdir -p ./02_Assembly/megahit-large/${sample}_ecc
    mv ~/${sample}_ecc/ ./02_Assembly/megahit-large

    #IDBA can only take fasta files as input; therefore, fastq.gz files must be reformated to fasta files
    reformat.sh in=./01_HQ_Reads/error_corrected/${sample}_HQ_ecc.fq.gz out=./01_HQ_Reads/error_corrected/${sample}_HQ_ecc.fa ow=t >> ./02_Assembly/idba_ud/logs/${sample}.log 2>&1
    #fq2fa --paired ./01_HQ_Reads/error_corrected/${sample}_HQ_ecc.fq.gz ./01_HQ_Reads/error_corrected/${sample}_HQ_ecc.fa
    # #Assemble reads per site using idba_ud
    echo Using idba to assemble ${HQ_Reads} ${sample} reads.

    mkdir -p ./02_Assembly/idba_ud/${sample}_ecc
    idba_ud -r ./01_HQ_Reads/error_corrected/${sample}_HQ_ecc.fa --num_threads=${CPU} --pre_correction -o ./02_Assembly/idba_ud/${sample}_ecc >> ./02_Assembly/idba_ud/logs/${sample}.log 2>&1
    cp ./02_Assembly/idba_ud/${sample}_ecc/scaffold.fa ./02_Assembly/quast/${sample}/idba_ecc.fasta

    # --- Evaluation ---
    echo metagenome_id,reads,assembler,n_contigs,contig_bp,gap_pct,ctg_L50,ctg_max,gene_counts,mapped_pct,avg_cov,assembly_performance >> ./02_Assembly/stats/${sample}_assembly_stats.csv
    echo Evaluating ${sample} assemblies with AssemblyStats
    stats.sh ./02_Assembly/metaspades/${sample}_ecc/scaffolds.fasta minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/metaspades/${sample}_ecc/stats.txt
    spades_stats=$(cat ./02_Assembly/metaspades/${sample}_ecc/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d' | sed 's/\t/,/g')
    stats.sh ./02_Assembly/megahit-large/${sample}_ecc/final.contigs.fa minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/megahit-large/${sample}_ecc/stats.txt
    megahit_stats=$(cat ./02_Assembly/megahit-large/${sample}_ecc/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d' | sed 's/\t/,/g')
    stats.sh ./02_Assembly/idba_ud/${sample}_ecc/scaffold.fa minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/idba_ud/${sample}_ecc/stats.txt
    idba_stats=$(cat ./02_Assembly/idba_ud/${sample}_ecc/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d' | sed 's/\t/,/g')
    #Evaluate assemblies with Quast (leave out "-R ref.fa if you don't have a reference)

    echo Predicting prokaryotic genes in each ${sample} assembly
    #Count gene predictions for prokaryotes using prodigal for each assembly
    python ~/bin/LengthFilter.py ./02_Assembly/metaspades/${sample}_ecc/scaffolds.fasta -m 1000 > ./02_Assembly/metaspades/${sample}_ecc/${sample}_c1k.fa
    prodigal -i ./02_Assembly/metaspades/${sample}_ecc/${sample}_c1k.fa -a ./02_Assembly/metaspades/${sample}_ecc/${sample}_c1k_proteins.faa -f gff -p meta -o ./02_Assembly/metaspades/${sample}_ecc/${sample}_c1k.gff >> ./02_Assembly/metaspades/logs/${sample}.log 2>&1
    spades_counts=$(grep '^>' ./02_Assembly/metaspades/${sample}_ecc/${sample}_c1k_proteins.faa | wc -l )

    python ~/bin/LengthFilter.py ./02_Assembly/megahit-large/${sample}_ecc/final.contigs.fa -m 1000 > ./02_Assembly/megahit-large/${sample}_ecc/${sample}_c1k.fa
    prodigal -i ./02_Assembly/megahit-large/${sample}_ecc/${sample}_c1k.fa -a ./02_Assembly/megahit-large/${sample}_ecc/${sample}_c1k_proteins.faa -f gff -p meta -o ./02_Assembly/megahit-large/${sample}_ecc/${sample}_c1k.gff >> ./02_Assembly/megahit-large/logs/${sample}.log 2>&1
    megahit_counts=$(grep '^>' ./02_Assembly/megahit-large/${sample}_ecc/${sample}_c1k_proteins.faa | wc -l )

    python ~/bin/LengthFilter.py ./02_Assembly/idba_ud/${sample}_ecc/scaffold.fa -m 1000 > ./02_Assembly/idba_ud/${sample}_ecc/${sample}_c1k.fa
    prodigal -i ./02_Assembly/idba_ud/${sample}_ecc/${sample}_c1k.fa -a ./02_Assembly/idba_ud/${sample}_ecc/${sample}_c1k_proteins.faa -f gff -p meta -o ./02_Assembly/idba_ud/${sample}_ecc/${sample}_c1k.gff >> ./02_Assembly/idba_ud/logs/${sample}.log 2>&1
    idba_counts=$(grep '^>' ./02_Assembly/idba_ud/${sample}_ecc/${sample}_c1k_proteins.faa | wc -l )

    echo Determining metagenome content within ${sample} assemblies
    # Note that these are reads mapped to contigs > 2,000 bp.
    #(i.e. % mapped) by mapping reads by to assembly
    bbmap.sh in=./01_HQ_Reads/${sample}_HQ.fq.gz ref=./02_Assembly/metaspades/${sample}_ecc/${sample}_c1k.fa nodisk covhist=./02_Assembly/metaspades/${sample}_ecc/covhist.txt covstats=./02_Assembly/metaspades/${sample}_ecc/covstats.txt out=./02_Assembly/metaspades/${sample}_ecc/${sample}_reads_assembled.fq.gz outu=./02_Assembly/metaspades/${sample}_ecc/${sample}_reads_unassembled.fq.gz fast=t ambig=best >> ./02_Assembly/metaspades/${sample}_ecc/map.log 2>&1
    spades_match=$(grep Percent\ mapped: ./02_Assembly/metaspades/${sample}_ecc/map.log | cut -f2 -d$'\t' | tail -1 | tail -1)
    spades_cov=$(grep Average\ coverage: ./02_Assembly/metaspades/${sample}_ecc/map.log | cut -f2 -d$'\t' | tail -1 | tail -1)
    spades_Pa=$(echo "$spades_match*$spades_cov" | bc)

    bbmap.sh in=./01_HQ_Reads/${sample}_HQ.fq.gz ref=./02_Assembly/megahit-large/${sample}_ecc/${sample}_c1k.fa nodisk covhist=./02_Assembly/megahit-large/${sample}_ecc/covhist.txt covstats=./02_Assembly/megahit-large/${sample}_ecc/covstats.txt outm=./02_Assembly/megahit-large/${sample}_ecc/${sample}_reads_assembled.fq.gz outu=./02_Assembly/megahit-large/${sample}_ecc/${sample}_reads_unassembled.fq.gz fast=t ambig=best >> ./02_Assembly/megahit-large/${sample}_ecc/map.log 2>&1
    megahit_match=$(grep Percent\ mapped: ./02_Assembly/megahit-large/${sample}_ecc/map.log | cut -f2 -d$'\t' | tail -1)
    megahit_cov=$(grep Average\ coverage: ./02_Assembly/megahit-large/${sample}_ecc/map.log | cut -f2 -d$'\t' | tail -1)
    megahit_Pa=$(echo "$megahit_match*$megahit_cov" | bc)

    bbmap.sh in=./01_HQ_Reads/${sample}_HQ.fq.gz ref=./02_Assembly/idba_ud/${sample}_ecc/${sample}_c1k.fa nodisk covhist=./02_Assembly/idba_ud/${sample}_ecc/covhist.txt covstats=./02_Assembly/idba_ud/${sample}_ecc/covstats.txt outm=./02_Assembly/idba_ud/${sample}_ecc/${sample}_reads_assembled.fq.gz outu=./02_Assembly/idba_ud/${sample}_ecc/${sample}_reads_unassembled.fq.gz fast=t ambig=best >> ./02_Assembly/idba_ud/${sample}_ecc/map.log 2>&1
    idba_match=$(grep Percent\ mapped: ./02_Assembly/idba_ud/${sample}_ecc/map.log | cut -f2 -d$'\t' | tail -1)
    idba_cov=$(grep Average\ coverage: ./02_Assembly/idba_ud/${sample}_ecc/map.log | cut -f2 -d$'\t' | tail -1)
    idba_Pa=$(echo "$idba_match*$idba_cov" | bc)

    echo Creating assembly comparision for ${sample}
    echo ${sample},${HQ_Reads},metaspades,${spades_stats},${spades_counts},${spades_match},${spades_cov},${spades_Pa} >> ./02_Assembly/stats/${sample}_assembly_stats.csv
    echo ${sample},${HQ_Reads},megahit,${megahit_stats},${megahit_counts},${megahit_match},${megahit_cov},${megahit_Pa} >> ./02_Assembly/stats/${sample}_assembly_stats.csv
    echo ${sample},${HQ_Reads},idba_ud,${idba_stats},${idba_counts},${idba_match},${idba_cov},${idba_Pa} >> ./02_Assembly/stats/${sample}_assembly_stats.csv

  elif [[ "${HQ_Reads}" == HQ_merged ]]
  then
    echo Using metaspades to assemble ${HQ_Reads} ${sample} reads.
    # Assemble using merged/unmerged reads
    mkdir -p ./02_Assembly/metaspades/${sample}_merged
    metaspades.py  --only-assembler -k 21,33,55,77,99,127 --merged ./01_HQ_Reads/merged_reads/${sample}_HQ_merged.fq.gz --12 ./01_HQ_Reads/merged_reads/${sample}_HQ_unmerged.fq.gz -t ${CPU} -m ${RAM} -o ./02_Assembly/metaspades/${sample}_merged >> ./02_Assembly/metaspades/logs/${sample}.log 2>&1
    cp ./02_Assembly/metaspades/${sample}_merged/scaffolds.fa ./02_Assembly/quast/${sample}/metaspades_merged.fasta
    echo Using megahit to assemble ${HQ_Reads} ${sample} reads.
#    megahit -r ./01_HQ_Reads/merged_reads/${sample}_HQ_merged.fq.gz --12 ./01_HQ_Reads/merged_reads/${sample}_HQ_unmerged.fq.gz --presets meta-large -m ${RAM} -t ${CPU} -o ./02_Assembly/megahit-large/${sample} >> ./02_Assembly/megahit-large/logs/${sample}.log 2>&1
    megahit -r ./01_HQ_Reads/merged_reads/${sample}_HQ_merged.fq.gz --12 ./01_HQ_Reads/merged_reads/${sample}_HQ_unmerged.fq.gz --k-min 21 --k-max 91 --k-step 10 -m ${RAM} -t ${CPU} -o ~/${sample}_merged >> ./02_Assembly/megahit-large/logs/${sample}.log 2>&1
    cp ~/${sample}_merged/final.contigs.fa ./02_Assembly/quast/${sample}/megahit_merged.fasta
    mkdir -p ./02_Assembly/megahit-large/${sample}_merged/
    mv ~/${sample}_merged/ ./02_Assembly/megahit-large/

    echo idba_ud does not assemble merged reads and will not be used for comparison.

    echo metagenome_id,reads,assembler,n_contigs,contig_bp,gap_pct,ctg_L50,ctg_max,gene_counts,mapped_pct,avg_cov,assembly_performance >> ./02_Assembly/stats/${sample}_assembly_stats.csv
    echo Evaluating ${sample} assemblies with AssemblyStats
    stats.sh ./02_Assembly/metaspades/${sample}_merged/scaffolds.fasta minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/metaspades/${sample}_merged/stats.txt
    spades_stats=$(cat ./02_Assembly/metaspades/${sample}_merged/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d' | sed 's/\t/,/g')
    stats.sh ./02_Assembly/megahit-large/${sample}_merged/final.contigs.fa minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/megahit-large/${sample}_merged/stats.txt
    megahit_stats=$(cat ./02_Assembly/megahit-large/${sample}_merged/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d' | sed 's/\t/,/g')
    #Evaluate assemblies with Quast (leave out "-R ref.fa if you don't have a reference)

    echo Predicting prokaryotic genes in each ${sample} assembly
    #Count gene predictions for prokaryotes using prodigal for each assembly
    python ~/bin/LengthFilter.py ./02_Assembly/metaspades/${sample}_merged/scaffolds.fasta -m 1000 > ./02_Assembly/metaspades/${sample}_merged/${sample}_c1k.fa
    prodigal -i ./02_Assembly/metaspades/${sample}_merged/${sample}_c1k.fa -a ./02_Assembly/metaspades/${sample}_merged/${sample}_c1k_proteins.faa -f gff -p meta -o ./02_Assembly/metaspades/${sample}_merged/${sample}_c1k.gff >> ./02_Assembly/metaspades/logs/${sample}.log 2>&1
    spades_counts=$(grep '^>' ./02_Assembly/metaspades/${sample}_merged/${sample}_c1k_proteins.faa | wc -l )

    python ~/bin/LengthFilter.py ./02_Assembly/megahit-large/${sample}_merged/final.contigs.fa -m 1000 > ./02_Assembly/megahit-large/${sample}_merged/${sample}_c1k.fa
    prodigal -i ./02_Assembly/megahit-large/${sample}_merged/${sample}_c1k.fa -a ./02_Assembly/megahit-large/${sample}_merged/${sample}_c1k_proteins.faa -f gff -p meta -o ./02_Assembly/megahit-large/${sample}_merged/${sample}_c1k.gff >> ./02_Assembly/megahit-large/logs/${sample}.log 2>&1
    megahit_counts=$(grep '^>' ./02_Assembly/megahit-large/${sample}_merged/${sample}_c1k_proteins.faa | wc -l )

    echo Determining metagenome content within ${sample} assemblies
    # Note that these are reads mapped to contigs > 2,000 bp.
    #(i.e. % mapped) by mapping reads by to assembly
    bbmap.sh in=./01_HQ_Reads/${sample}_HQ.fq.gz ref=./02_Assembly/metaspades/${sample}_merged/${sample}_c1k.fa nodisk covhist=./02_Assembly/metaspades/${sample}_merged/covhist.txt covstats=./02_Assembly/metaspades/${sample}_merged/covstats.txt out=./02_Assembly/metaspades/${sample}_merged/${sample}_reads_assembled.fq.gz outu=./02_Assembly/metaspades/${sample}_merged/${sample}_reads_unassembled.fq.gz fast=t ambig=best >> ./02_Assembly/metaspades/${sample}_merged/map.log 2>&1
    spades_match=$(grep Percent\ mapped: ./02_Assembly/metaspades/${sample}_merged/map.log | cut -f2 -d$'\t' | tail -1)
    spades_cov=$(grep Average\ coverage: ./02_Assembly/metaspades/${sample}_merged/map.log | cut -f2 -d$'\t' | tail -1)
    spades_Pa=$(echo "$spades_match*$spades_cov" | bc)

    bbmap.sh in=./01_HQ_Reads/${sample}_HQ.fq.gz ref=./02_Assembly/megahit-large/${sample}_merged/${sample}_c1k.fa nodisk covhist=./02_Assembly/megahit-large/${sample}_merged/covhist.txt covstats=./02_Assembly/megahit-large/${sample}_merged/covstats.txt outm=./02_Assembly/megahit-large/${sample}_merged/${sample}_reads_assembled.fq.gz outu=./02_Assembly/megahit-large/${sample}_merged/${sample}_reads_unassembled.fq.gz fast=t ambig=best >> ./02_Assembly/megahit-large/${sample}_merged/map.log 2>&1
    megahit_match=$(grep Percent\ mapped: ./02_Assembly/megahit-large/${sample}_merged/map.log | cut -f2 -d$'\t' | tail -1)
    megahit_cov=$(grep Average\ coverage: ./02_Assembly/megahit-large/${sample}_merged/map.log | cut -f2 -d$'\t' | tail -1)
    metahit_Pa=$(echo "$megahit_match*$megahit_cov" | bc)

    echo Creating assembly comparision for ${sample}
    echo ${sample},${HQ_Reads},metaspades,${spades_stats},${spades_counts},${spades_match},${spades_cov},${spades_Pa} >> ./02_Assembly/stats/${sample}_assembly_stats.csv
    echo ${sample},${HQ_Reads},megahit,${megahit_stats},${megahit_counts},${megahit_match},${megahit_cov},${megahit_Pa} >> ./02_Assembly/stats/${sample}_assembly_stats.csv
    echo ${sample},${HQ_Reads},idba_ud,${idba_stats},${idba_counts},${idba_match},${idba_cov},NA >> ./02_Assembly/stats/${sample}_assembly_stats.csv

  elif [[ "${HQ_Reads}" == Skip ]]
  then :
  fi

# --- Evaluation ---

  echo Evaluating ${sample} assemblies

  #statswrapper.sh ./02_Assembly/quast/${sample}/*.fasta minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/stats/${sample}.txt
  ~/GitHub/quast/metaquast.py -f -m 1000 -o ./02_Assembly/quast/${sample} ./02_Assembly/quast/${sample}/*.fasta >> ./02_Assembly/quast/${sample}.log 2>&1



  read -r -p "Analyze the ./02_Assembly/stats/${sample}_assembly_stats.csv and ./02_Assembly/quast/${sample}/report.html and determine which ${sample} assembly is best. Please enter the absoute path to the assembly.fasta file:  " best_assembly || exit 100
  # #Pick which assembly you like best by comparing assembly_stats.csv
  # Copy best assembly to 02_Assembly directory
  # For my analysis, 58 of the best assemblies were megahit and 2 were metaspades.
  mkdir -p ./02_Assembly/taxonomy
  echo Determining the overall taxonomic makeup of the ${sample} ${best_assembly} assembly
  sendsketch.sh in=${best_assembly} >> ./02_Assembly/taxonomy/${sample}_RefSeq.tsv 2>&1
  sendsketch.sh in=${best_assembly} nt >> ./02_Assembly/taxonomy/${sample}_nt.tsv 2>&1

  conda deactivate
  #Or, try to determine taxonomy on a per-contig basis.  If this is not sensitive enough, try BLAST instead.
  # sendsketch.sh in=${best_assembly} persequence minhits=1 records=4 >> ./02_Assembly/taxonomy/${sample}.log 2>&1
  echo Starting anvi\'o on ${sample}

  conda activate anvio-7.1

  prefix=$(echo $sample | sed 's/-/_/g')
  echo running anvio on ${sample}
  mkdir -p ./03_Anvio/${sample}/COG
  mkdir -p ./03_Anvio/${sample}/KOfam
  mkdir -p ./03_Anvio/${sample}/Pfam
  mkdir -p ./03_Anvio/${sample}/prodigal
  mkdir -p ./03_Anvio/${sample}/bowtie2
  mkdir -p ./03_Anvio/logs
  #Simplify headers as to not cause future anger; make a key so that contigs can be matched up later if needed; length cut off of 500bp; the prefix argument doesn't take hyphens
  echo Reformating ${sample} fasta file and removing contigs \< 1000bp
  anvi-script-reformat-fasta --seq-type NT --simplify-names --prefix ${prefix} -r ./03_Anvio/${sample}/${sample}_rename_key.txt -l 1000 -o ./03_Anvio/${sample}/${sample}_renamed_c1k.fa ${best_assembly} >> ./03_Anvio/logs/${sample}.log 2>&1
  #Generate a contigs database for each assembly
  echo making contigs database for ${sample} assembly
  anvi-gen-contigs-database -T ${CPU} -f ./03_Anvio/${sample}/${sample}_renamed_c1k.fa -n ${sample} -o ./03_Anvio/${sample}/${sample}.db >> ./03_Anvio/logs/${sample}.log 2>&1
  echo running hmm
  anvi-run-hmms -c ./03_Anvio/${sample}/${sample}.db -T ${CPU} --also-scan-trnas >> ./03_Anvio/logs/${sample}.log 2>&1
  echo running ncbi COGS
  anvi-run-ncbi-cogs -c ./03_Anvio/${sample}/${sample}.db -T ${CPU} --sensitive >> ./03_Anvio/logs/${sample}.log 2>&1
  echo running KEGG
  anvi-run-kegg-kofams -c ./03_Anvio/${sample}/${sample}.db -T ${CPU} >> ./03_Anvio/logs/${sample}.log 2>&1
  echo running pfams
  anvi-run-pfams -c ./03_Anvio/${sample}/${sample}.db -T ${CPU} >> ./03_Anvio/logs/${sample}.log 2>&1
  echo exporting functions to contigs database
  anvi-export-functions -c ./03_Anvio/${sample}/${sample}.db -o ./03_Anvio/${sample}/COG/${sample}_COG20Category.txt --annotation-sources COG20_CATEGORY >> ./03_Anvio/logs/${sample}.log 2>&1
  anvi-export-functions -c ./03_Anvio/${sample}/${sample}.db -o ./03_Anvio/${sample}/COG/${sample}_COG20Function.txt --annotation-sources COG20_FUNCTION >> ./03_Anvio/logs/${sample}.log 2>&1
  anvi-export-functions -c ./03_Anvio/${sample}/${sample}.db -o ./03_Anvio/${sample}/KOfam/${sample}_KOfam.txt --annotation-sources KOfam >> ./03_Anvio/logs/${sample}.log 2>&1
  anvi-export-functions -c ./03_Anvio/${sample}/${sample}.db -o ./03_Anvio/${sample}/Pfam/${sample}_Pfam.txt --annotation-sources Pfam >> ./03_Anvio/logs/${sample}.log 2>&1
  echo exporing gene calls to contig database
  anvi-export-gene-calls -c ./03_Anvio/${sample}/${sample}.db --gene-caller prodigal -o ./03_Anvio/${sample}/prodigal/${sample}_AllGeneCalls.txt >> ./03_Anvio/logs/${sample}.log 2>&1
  echo running SCG taxonomy
  anvi-run-scg-taxonomy -c ./03_Anvio/${sample}/${sample}.db -T ${CPU} >> ./03_Anvio/logs/${sample}.log 2>&1

  #Remove zeros from column 5 e-value and only write significant hits
  echo creating functional profile on ${sample}
  awk -F "\t" '{ if(($5 != 0) && ($5 <= .0001)) { print } }' ./03_Anvio/${sample}/KOfam/${sample}_KOfam.txt >> ./03_Anvio/${sample}/KOfam/${sample}_KOfamEvalue_filtered.txt
  awk -F "\t" '{ if(($5 != 0) && ($5 <= .0001)) { print } }' ./03_Anvio/${sample}/Pfam/${sample}_Pfam.txt >> ./03_Anvio/${sample}/Pfam/${sample}_PfamEvalue_filtered.txt
  awk -F "\t" '{ if(($5 != 0) && ($5 <= .0001)) { print } }' ./03_Anvio/${sample}/COG/${sample}_COG20Function.txt >> ./03_Anvio/${sample}/COG/${sample}_COG20FunctionEvalue_filtered.txt
  awk -F "\t" '{ if(($5 != 0) && ($5 <= .0001)) { print } }' ./03_Anvio/${sample}/COG/${sample}_COG20Category.txt >> ./03_Anvio/${sample}/COG/${sample}_COG20CategoryEvalue_filtered.txt
  #Cut unique values from column 4 function, sort them and then count them
  cut -f4 ./03_Anvio/${sample}/KOfam/${sample}_KOfamEvalue_filtered.txt | sort | uniq -c >> ./03_Anvio/${sample}/KOfam/${sample}_KOfamCounts.txt
  cut -f4 ./03_Anvio/${sample}/Pfam/${sample}_PfamEvalue_filtered.txt | sort | uniq -c >> ./03_Anvio/${sample}/Pfam/${sample}_PfamCounts.txt
  cut -f4 ./03_Anvio/${sample}/COG/${sample}_COG20FunctionEvalue_filtered.txt | sort | uniq -c >> ./03_Anvio/${sample}/COG/${sample}_COG20Function.txt
  cut -f4 ./03_Anvio/${sample}/COG/${sample}_COG20CategoryEvalue_filtered.txt | sort | uniq -c >> ./03_Anvio/${sample}/COG/${sample}_COG20Category.txt

  conda deactivate

  read -r -p "What minimum contig length do you want to analyze for biosynthetic gene clusters?    " min_BGC || exit 100

  mkdir -p ./04_Annotations/antismash/${sample}_c${min_BGC}
  echo Running antismash to annotate ${sample} secondary metabolites

  ~/bin/run_antismash \
  ./03_Anvio/${sample}/${sample}_renamed_c1k.fa \
  ./04_Annotations/antismash/${sample}_c${min_BGC} \
  --hmmdetection-strictness relaxed \
  --logfile ${sample}_c${min_BGC}.log \
  --cb-general \
  --cb-knownclusters \
  --cb-subclusters \
  --asf --rre --pfam2go --smcog-trees \
  --genefinding-tool prodigal-m \
  --cc-mibig \
  --minlength ${min_BGC} \
  --logfile ./04_Annotations/antismash/${sample}_c${min_BGC}/${sample}.log \
  -c ${CPU} >> ./04_Annotations/antismash/${sample}.log 2>&1

  mkdir -p ./04_Annotations/bigscape/c${min_BGC}

  cd ./04_Annotations/antismash/${sample}_c${min_BGC}/${sample}_renamed_c1k

  echo Running bigscape on ${sample} antismash genbank files

  if test -f "c*.gbk"; then
    for genbank_file in $(ls c*.gbk); do
      cluster=$(echo ${genbank_file} | cut -f1 -d '_')
      echo copying and renaming ${genbank_file} to cluster_${cluster}_${sample}.gbk
      cp ${genbank_file} ${project_dir}/04_Annotations/bigscape/c${min_BGC}/cluster_${cluster}_${sample}.gbk
    done
  else echo No BGC found in ${sample}
  fi
  cd ${project_dir}
done

cd ${project_dir}/04_Annotations/bigscape

~/bin/run_bigscape.py c${min_BGC} c${min_BGC}_results --mix --mibig -c ${CPU} >> ./04_Annotations/bigscape/${project}.log 2&>1



# echo starting to assemble HQ reads into one co-assembly
#
# #Concatenate HQ reads into one file.
# cat ./01_HQ_Reads/*.fq.gz >> ./01_HQ_Reads/${project}_HQ.fq.gz
# cat ./01_HQ_Reads/merged_reads/*_merged.fq.gz >> ./01_HQ_Reads/merged_reads/${project}_HQ_merged.fq.gz
# cat ./01_HQ_Reads/merged_reads/*_unmerged.fq.gz >> ./01_HQ_Reads/merged_reads/${project}_HQ_unmerged.fq.gz

# mkdir -p ./02_Assembly/quast/co-assembly
# mkdir -p ./02_Assembly/stats/co-assembly
# mkdir -p ./02_Assembly/taxonomy/co-assembly

# echo Co-assembling ${project}
# #Assemble reads per sample using metaspades
# echo start metaspades on ${project}
# mkdir -p ./02_Co-Assembly/metaspades/merged_reads
# # Assemble using HQ reads
# metaspades.py  --only-assembler -k 21,33,55,77,99,127 --12 ./01_HQ_Reads/${project}_HQ.fq.gz -t 90 -m 900 -o ./02_Co-Assembly/metaspades >> ./02_Co-Assembly/metaspades/logs/${project}.log 2>&1
# mv ./02_Co-Assembly/metaspades/scaffolds.fasta ./02_Co-Assembly/metaspades/metaspades.fasta
# # Assemble using merged/unmerged reads
# metaspades.py  --only-assembler -k 21,33,55,77,99,127 --merged ./01_HQ_Reads/merged_reads/${project}_HQ_merged.fq.gz --12 ./01_HQ_Reads/merged_reads/${project}_HQ_unmerged.fq.gz -t 90 -m 900 -o ./02_Co-Assembly/metaspades/merged_reads >> ./02_Co-Assembly/metaspades/logs/${project}.log 2>&1
# mv ./02_Co-Assembly/metaspades/merged_reads/scaffolds.fasta ./02_Co-Assembly/metaspades/merged_reads/metaspades_merged.fasta
#
# #Assemble reads per site using megahit
#
# echo start megahit on ${project}
#
# # Assemble using normalized reads on meta-large pre-setting
# megahit --12 ./01_HQ_Reads/${project}_HQ.fq.gz --presets meta-large -m ${RAM} -t 90 -o ./02_Co-Assembly/megahit-large >> ./02_Co-Assembly/megahit-large/logs/${project}.log 2>&1
# mv  ./02_Co-Assembly/megahit-large/final.contigs.fa ./02_Co-Assembly/megahit-large/megahit_lg.fasta
# # Assemble using normalized and merged/unmerged reads on meta-large pre-setting
# megahit -r ./01_HQ_Reads/merged_reads/${project}_HQ_merged.fq.gz --12 ./01_HQ_Reads/merged_reads/${project}_HQ_unmerged.fq.gz --presets meta-large -m ${RAM} -t 90 -o ./02_Co-Assembly/megahit-large/merged_reads >> ./02_Co-Assembly/megahit-large/logs/${project}.log 2>&1
# mv  ./02_Co-Assembly/megahit-large/merged_reads/final.contigs.fa ./02_Co-Assembly/megahit-large/merged_reads/megahit_lg_merged.fasta
#
#
# #IDBA can only take fasta files as input; therefore, fastq.gz files must be reformated to fasta files
# reformat.sh in=./01_HQ_Reads/${project}_HQ.fq.gz out=./01_HQ_Reads/${project}_HQ.fa -Xmx${RAM}G ow=t >> ./02_Co-Assembly/logs/${project}.log 2>&1
#
# # #Assemble reads per site using idba_ud
# echo start idba_ud on ${project}
# mkdir ./02_Co-Assembly/idba_ud
# idba_ud -r ./01_HQ_Reads/${project}_HQ.fa --num_threads=${CPU} --pre_correction -o ./02_Co-Assembly/idba_ud >> ./02_Co-Assembly/idba_ud/logs/${project}.log 2>&1
# mv ./02_Co-Assembly/idba_ud/scaffold.fa ./02_Co-Assembly/idba_ud/idba_ud.fasta
#
#
# # --- Evaluation ---
#
# #Evaluate assemblies with Quast (leave out "-R ref.fa if you don't have a reference)
# metaquast.py -o ./02_Co-Assembly/quast/co-assembly ./02_Co-Assembly/idba_ud/idba_ud.fasta ./02_Co-Assembly/megahit-large/megahit_lg.fasta ./02_Co-Assembly/metaspades/metaspades.fasta ./02_Co-Assembly/megahit-large/merged_reads/megahit_lg_merged.fasta ./02_Co-Assembly/metaspades/merged_reads/metaspades_merged.fasta  >> ./02_Co-Assembly/quast/${project}.log 2>&1
#
# echo metagenome_id,assembler,n_contigs,contig_bp,gap_pct,ctg_L50,ctg_max,gene_counts,mapped_pct,avg_cov,assembly_performance >> ./02_Co-Assembly/stats/${project}_co-assembly_stats.csv
# echo Evaluating ${project} assemblies with AssemblyStats
# # #statswrapper.sh ./02_Co-Assembly/quast/metaspades.fasta ./02_Co-Assembly/quast/megahit.fasta ./02_Co-Assembly/quast/idba.fasta minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Co-Assembly/stats
# stats.sh ./02_Co-Assembly/metaspades/metaspades.fasta minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Co-Assembly/metaspades/stats.txt
# spades_stats=$(cat ./02_Co-Assembly/metaspades/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d' | sed 's/\t/,/g')
# stats.sh ./02_Co-Assembly/megahit-large/megahit.fasta minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Co-Assembly/megahit-large/stats.txt
# megahit_stats=$(cat ./02_Co-Assembly/megahit-large/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d' | sed 's/\t/,/g')
# stats.sh ./02_Co-Assembly/idba_ud/idba.fasta minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Co-Assembly/idba_ud/stats.txt
# idba_stats=$(cat ./02_Co-Assembly/idba_ud/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d' | sed 's/\t/,/g')
# # Evaluate stats for merged assemblies
# stats.sh ./02_Co-Assembly/metaspades/merged_reads/metaspades_merged.fasta minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Co-Assembly/metaspades/merged_reads/stats.txt
# spades_merged_stats=$(cat ./02_Co-Assembly/metaspades/merged_reads/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d' | sed 's/\t/,/g')
# stats.sh ./02_Co-Assembly/megahit-large/merged_reads/megahit_lg_merged.fasta minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Co-Assembly/megahit-large/merged_reads/stats.txt
# megahit_merged_stats=$(cat ./02_Co-Assembly/megahit-large/merged_reads/tats.txt| cut -f2,4,5,9,15 -d$'\t' | sed '2q;d' | sed 's/\t/,/g')
#
# echo Determining metagenome content within ${project} assemblies
# #(i.e. % mapped) by mapping reads by to assembly
# bbmap.sh in=./01_HQ_Reads/${project}_HQ.fq.gz ref=./02_Co-Assembly/metaspades/metaspades.fasta nodisk covhist=./02_Co-Assembly/metaspades/covhist.txt covstats=./02_Co-Assembly/metaspades/covstats.txt outm=./02_Co-Assembly/metaspades/${sample}_reads_assembled.fq.gz outu=./02_Co-Assembly/metaspades/${sample}_reads_unassembled.fq.gz fast=t ambig=best >> ./02_Co-Assembly/metaspades/map.log 2>&1
# spades_match=$(grep Percent\ mapped: ./02_Co-Assembly/metaspades/map.log | cut -f2 -d$'\t' | tail -1)
# spades_cov=$(grep Average\ coverage: ./02_Co-Assembly/metaspades/map.log | cut -f2 -d$'\t' | tail -1)
#
# bbmap.sh in=./01_HQ_Reads/${project}_HQ.fq.gz ref=./02_Co-Assembly/megahit-large/megahit.fasta nodisk covhist=./02_Co-Assembly/megahit-large/covhist.txt covstats=./02_Co-Assembly/megahit-large/covstats.txt outm=./02_Co-Assembly/megahit-large/${sample}_reads_assembled.fq.gz outu=./02_Co-Assembly/megahit-large/${sample}_reads_unassembled.fq.gz fast=t ambig=best >> ./02_Co-Assembly/megahit-large/map.log 2>&1
# megahit_match=$(grep Percent\ mapped: ./02_Co-Assembly/megahit-large/map.log | cut -f2 -d$'\t' | tail -1)
# megahit_cov=$(grep Average\ coverage: ./02_Co-Assembly/megahit-large/map.log | cut -f2 -d$'\t' | tail -1)
#
# bbmap.sh in=./01_HQ_Reads/${project}_HQ.fq.gz ref=./02_Co-Assembly/idba_ud/idba.fasta nodisk covhist=./02_Co-Assembly/idba_ud/covhist.txt covstats=./02_Co-Assembly/idba_ud/covstats.txt outm=./02_Co-Assembly/idba_ud/${sample}_reads_assembled.fq.gz outu=./02_Co-Assembly/idba_ud/${sample}_reads_unassembled.fq.gz fast=t ambig=best >> ./02_Co-Assembly/idba_ud/map.log 2>&1
# idba_match=$(grep Percent\ mapped: ./02_Co-Assembly/idba_ud/map.log | cut -f2 -d$'\t' | tail -1)
# idba_cov=$(grep Average\ coverage: ./02_Co-Assembly/idba_ud/map.log | cut -f2 -d$'\t' | tail -1)
#
# #Evaluate genomic content of merged assemblies
# bbmap.sh in=./01_HQ_Reads/${project}_HQ.fq.gz ref=./02_Co-Assembly/metaspades/merged_reads/metaspades_merged.fasta nodisk covhist=./02_Co-Assembly/metaspades/merged_reads/covhist.txt covstats=./02_Co-Assembly/metaspades/merged_reads/covstats.txt outm=./02_Co-Assembly/metaspades/merged_reads/${sample}_reads_assembled.fq.gz outu=./02_Co-Assembly/metaspades/merged_reads/${sample}_reads_unassembled.fq.gz fast=t ambig=best >> ./02_Co-Assembly/metaspades/merged_reads/map.log 2>&1
# spades_merged_match=$(grep Percent\ mapped: ./02_Co-Assembly/metaspades/merged_reads/map.log | cut -f2 -d$'\t' | tail -1)
# spades__merged_cov=$(grep Average\ coverage: ./02_Co-Assembly/metaspades/merged_reads/map.log | cut -f2 -d$'\t' | tail -1)
#
# bbmap.sh in=./01_HQ_Reads/${project}_HQ.fq.gz ref=./02_Co-Assembly/megahit-large/merged_reads/megahit_lg_merged.fasta nodisk covhist=./02_Co-Assembly/megahit-large/merged_reads/covhist.txt covstats=./02_Co-Assembly/megahit-large/merged_reads/covstats.txt outm=./02_Co-Assembly/megahit-large/merged_reads/${sample}_reads_assembled.fq.gz outu=./02_Co-Assembly/megahit-large/merged_reads/${sample}_reads_unassembled.fq.gz fast=t ambig=best >> ./02_Co-Assembly/megahit-large/merged_reads/map.log 2>&1
# megahit_merged_match=$(grep Percent\ mapped: ./02_Co-Assembly/megahit-large/merged_reads/map.log | cut -f2 -d$'\t' | tail -1)
# megahit_merged_cov=$(grep Average\ coverage: ./02_Co-Assembly/megahit-large/merged_reads/map.log | cut -f2 -d$'\t' | tail -1)
#
# echo Predicting prokaryotic genes in each ${project} assembly
# #Count gene predictions for prokaryotes using prodigal for each assembly
# python ~/bin/LengthFilter.py ./02_Co-Assembly/metaspades/metaspades.fasta -m 1000 > ./02_Co-Assembly/metaspades/${project}_c1k.fa
# prodigal -i ./02_Co-Assembly/metaspades/${project}_c1k.fa -a ./02_Co-Assembly/metaspades/${project}_c1k_proteins.faa -f gff -p meta -o ./02_Co-Assembly/metaspades/${project}_c1k.gff >> ./02_Co-Assembly/metaspades/logs/${project}.log 2>&1
# spades_counts=$(grep '^>' ./02_Co-Assembly/metaspades/${project}_c1k_proteins.faa | wc -l )
#
# python ~/bin/LengthFilter.py ./02_Co-Assembly/megahit-large/megahit.fasta -m 1000 > ./02_Co-Assembly/megahit-large/${project}_c1k.fa
# prodigal -i ./02_Co-Assembly/megahit-large/${project}_c1k.fa -a ./02_Co-Assembly/megahit-large/${project}_c1k_proteins.faa -f gff -p meta -o ./02_Co-Assembly/megahit-large/${project}_c1k.gff >> ./02_Co-Assembly/megahit-large/logs/${project}.log 2>&1
# megahit_counts=$(grep '^>' ./02_Co-Assembly/megahit-large/${project}_c1k_proteins.faa | wc -l )
#
# python ~/bin/LengthFilter.py ./02_Co-Assembly/idba_ud/idba.fasta -m 1000 > ./02_Co-Assembly/idba_ud/${project}_c1k.fa
# prodigal -i ./02_Co-Assembly/idba_ud/${project}_c1k.fa -a ./02_Co-Assembly/idba_ud/${project}_c1k_proteins.faa -f gff -p meta -o ./02_Co-Assembly/idba_ud/${project}_c1k.gff >> ./02_Co-Assembly/idba_ud/logs/${project}.log 2>&1
# idba_counts=$(grep '^>' ./02_Co-Assembly/idba_ud/${project}_c1k_proteins.faa | wc -l )
# # Evaluage the merged assemblies
# python ~/bin/LengthFilter.py ./02_Co-Assembly/metaspades/merged_reads/metaspades_merged.fasta -m 1000 > ./02_Co-Assembly/metaspades/merged_reads/${project}_c1k.fa
# prodigal -i ./02_Co-Assembly/metaspades/merged_reads/${project}_c1k.fa -a ./02_Co-Assembly/metaspades/merged_reads/${project}_c1k_proteins.faa -f gff -p meta -o ./02_Co-Assembly/metaspades/merged_reads/${project}_c1k.gff >> ./02_Co-Assembly/metaspades/merged_reads/logs/${project}.log 2>&1
# spades_merged_counts=$(grep '^>' ./02_Co-Assembly/metaspades/merged_reads/${project}_c1k_proteins.faa | wc -l )
#
# python ~/bin/LengthFilter.py ./02_Co-Assembly/megahit-large/merged_reads/megahit_lg_merged.fasta -m 1000 > ./02_Co-Assembly/megahit-large/merged_reads/${project}_c1k.fa
# prodigal -i ./02_Co-Assembly/megahit-large/merged_reads/${project}_c1k.fa -a ./02_Co-Assembly/megahit-large/merged_reads/${project}_c1k_proteins.faa -f gff -p meta -o ./02_Co-Assembly/megahit-large/merged_reads/${project}_c1k.gff >> ./02_Co-Assembly/megahit-large/merged_reads/logs/${project}.log 2>&1
# megahit_merged_counts=$(grep '^>' ./02_Co-Assembly/megahit-large/merged_reads/${project}_c1k_proteins.faa | wc -l )
#
# echo Creating assembly comparision for ${project}
# echo ${project},metaspades,${spades_stats},${spades_counts},${spades_match},${spades_cov} >> ./02_Co-Assembly/stats/${project}_co-assembly_stats.csv
# echo ${project},megahit,${megahit_stats},${megahit_counts},${megahit_match},${megahit_cov} >> ./02_Co-Assembly/stats/${project}_co-assembly_stats.csv
# echo ${project},idba_ud,${idba_stats},${idba_counts},${idba_match},${idba_cov} >> ./02_Co-Assembly/stats/${project}_co-assembly_stats.csv
# # Add merged evaluation
# echo ${project},metaspades_merged,${spades_merged_stats},${spades_merged_counts},${spades_merged_match},${spades_merged_cov} >> ./02_Co-Assembly/stats/${project}_co-assembly_stats.csv
# echo ${project},megahit_merged,${megahit_merged_stats},${megahit_merged_counts},${megahit_merged_match},${megahit_merged_cov} >> ./02_Co-Assembly/stats/${project}_co-assembly_stats.csv
#
# #Pick which assembly you like best by comparing ${project}_co-assembly_stats.csv
# read -r -p "Analyze co-assembly results at ./02_Co-Assembly/quast/co-assembly and ./02_Co-Assembly/stats/${project}_co-assembly_stats.csv to determine the best co-assembly. Please enter the absoute path to the co-assembly.fasta file:    " best_coassembly || exit 100
# echo Determining the overall taxonomic makeup of the ${sample} ${best_coassembly} assembly
# sendsketch.sh in=${best_coassembly} >> ./02_Co-Assembly/taxonomy/${project}.log 2>&1
#
# #Or, try to determine taxonomy on a per-contig basis.  If this is not sensitive enough, try BLAST instead.
# # sendsketch.sh in=./02_Co-Assembly/${project}.*.fa persequence minhits=1 records=4 >> ./02_Co-Assembly/taxonomy/${project}.log 2>&1
#
# echo Running anvio workflow on ${project} co-assembly
# mkdir -p ./03_anvio/COG
# mkdir -p ./03_anvio/KOfam
# mkdir -p ./03_anvio/Pfam
# mkdir -p ./03_anvio/prodigal
# mkdir -p ./03_anvio/logs
# mkdir -p ./03_anvio/profiles
# # Simplity headers as to not cause future anger; make a key so that contigs can be matched up later if needed; length cut off of 500bp; the prefix argument doesn't take hyphens
# echo Reformating ${project} co-assembly and removing short contigs
# anvi-script-reformat-fasta --seq-type NT --simplify-names --prefix ${project} -r ./02_Co-Assembly/${project}_rename_key.txt  -l 2000 -o ./03_anvio/${project}_reformated_assembly_c1k.fa ${best_coassembly} >> ./03_anvio/logs/${project}.log 2>&1
# #cp ./03_anvio/${project}_reformated_assembly_c1k.fa ./02_Co-Assembly/
# # Generate a contigs database for each assembly
# echo making contigs database for ${project} assembly
# anvi-gen-contigs-database -T 90 -f ./03_anvio/${project}_reformated_assembly_c1k.fa -n ${project} -o ./03_anvio/${project}.db >> ./03_anvio/logs/${project}.log 2>&1
# echo running hmm
# anvi-run-hmms -c ./03_anvio/${project}.db -T 90 --also-scan-trnas >> ./03_anvio/logs/${project}.log 2>&1
# # anvi-run-hmms -c ./03_anvio/${project}.db -T 90 >> ./03_anvio/logs/${project}.log 2>&1
#
# echo running ncbi COGS
# anvi-run-ncbi-cogs -c ./03_anvio/${project}.db -T 90 --sensitive >> ./03_anvio/logs/${project}.log 2>&1
# echo running KEGG
# anvi-run-kegg-kofams -c ./03_anvio/${project}.db -T 90 >> ./03_anvio/logs/${project}.log 2>&1
# echo running pfams
# anvi-run-pfams -c ./03_anvio/${project}.db -T 90 >> ./03_anvio/logs/${project}.log 2>&1
# echo exporting functions to contigs database
# anvi-export-functions -c ./03_anvio/${project}.db -o ./03_anvio/COG/${project}_COG20Category.txt --annotation-sources COG20_CATEGORY >> ./03_anvio/logs/${project}.log 2>&1
# anvi-export-functions -c ./03_anvio/${project}.db -o ./03_anvio/COG/${project}_COG20Function.txt --annotation-sources COG20_FUNCTION >> ./03_anvio/logs/${project}.log 2>&1
# anvi-export-functions -c ./03_anvio/${project}.db -o ./03_anvio/KOfam/${project}_KOfam.txt --annotation-sources KOfam >> ./03_anvio/logs/${project}.log 2>&1
# anvi-export-functions -c ./03_anvio/${project}.db -o ./03_anvio/Pfam/${project}_Pfam.txt --annotation-sources Pfam >> ./03_anvio/logs/${project}.log 2>&1
# echo exporing gene calls to contig database
# anvi-export-gene-calls -c ./03_anvio/${project}.db --gene-caller prodigal -o ./03_anvio/prodigal/${project}_AllGeneCalls.txt >> ./03_anvio/logs/${project}.log 2>&1
# echo running SCG taxonomy
# anvi-run-scg-taxonomy -c ./03_anvio/${project}.db -T 90 >> ./03_anvio/logs/${project}.log 2>&1
#
# # Index co-assembly in preparation for mapping
# echo Indexing ${project} assembly for mapping
# bbmap.sh ref=./03_anvio/${project}_reformated_assembly_c1k.fa
# for sample in `cat samples.txt`; do
#   echo creating bam file for ${sample}
#   sample_rev=$( echo ${sample} | sed 's/-/_/g' )
#   mkdir -p ./03_anvio/${sample}
#   bbmap.sh in=./01_HQ_Reads/${sample}_HQ_reads.fq.gz out=./03_anvio/${sample}.sam bamscript=bs.sh outm=./03_anvio/${sample}/${sample}_mapped.fq.gz outu=./03_anvio/${sample}/${sample}_unmapped.fq.gz ihist=./03_anvio/${sample}/${sample}_insertsize.txt interleaved=t pigz=t threads=12 -Xmx${RAM}G ow=t >> ./03_anvio/logs/${project}.log 2>&1
#   sh bs.sh >> ./03_anvio/logs/${project}.log 2>&1
#   #
#   echo creating anvi profile for sample ${sample}
#
#   anvi-profile -i ./03_anvio/${sample}_sorted.bam -c ./03_anvio/${project}.db -T 90 -M 900 --sample-name ${sample_rev} --output-dir ./03_anvio/profiles/${sample} >> ./03_anvio/logs/${project}.log 2>&1
# done
# echo merging sample profiles
# anvi-merge ./03_anvio/profiles/*/PROFILE.db -o ./03_anvio/${project}_merged -c ./03_anvio/${project}.db --enforce-hierarchical-clustering >> ./03_anvio/logs/${project}.log 2>&1
#
# ## Automatic Binning ##
# mkdir -p ./03_Binning/
# # Use anvio's built-in function to bin contigs
# for binner in `cat binners.txt`; do
#   anvi-cluster-contigs -c ./03_anvio/${project}.db -p /03_anvio/${project}_merged/PROFILE.db -C ${binner} --driver ${binner} -T 90 --log-file ./04_Mapping/${binner}.log --just-do-it
# done
#
# for binner in `cat binners.txt`; do
#   echo ${binner}
#   cd ./${binner}
#   for bin_fasta in `ls *.fasta`; do
#     # echo ${bin_fasta}
#     bin_number=$(echo ${bin_fasta} | cut -f2 -d '.')
#     # echo ${bin_number}
#     bin=( bin_${bin_number} )
#     # echo ${bin}
#     anvi-script-reformat-fasta --prefix ${binner} -l 1000 -o ${bin}_c1k.fasta ${bin_fasta}
#
#     for contig in `grep '^>' ${bin}_c1k.fasta | cut -f2 -d '>'`; do
#     #
#       echo ${contig},${bin} >> temp.txt
#     done
#   done
#   sed 's/\,/\t/g' temp.txt > ${binner}_binning_results.txt
#   rm temp.txt
#   #Import bins to contigs database
#   anvi-import-collection ${binner}_binning_results.txt -p ../../../03_anvio_c1k/FCHc02/FCHc02_merged/PROFILE.db -c ../../../03_anvio_c1k/FCHc02/FCHc02_c1k.db -C ${binner} --contigs-mode >> ../../../03_anvio_c1k/logs/FCHc02.log 2>&1
#
#   cd ..
# done
#
#
# conda deactivate
#
# #####
# ##Antismash###
# ####
#
# conda activate antismash
#
# mkdir -p ./04_Annotations/antismash
# antismash \
# --cb-general \
# --cb-knownclusters \
# --cb-subclusters \
# --asf --pfam2go --smcog-trees \
# ./03_anvio/${project}_reformated_assembly_c1k.fa \
# --genefinding-tool prodigal-m \
# --output-dir ./04_Annotations/antismash >> ./04_Annotations/antismash/${project}.log 2>&1
#
# conda deactivate
#
# mkdir ./04_Annotations/bigscape
#
# for sample in `cat samples.txt`; do
#   # echo ${sample}
#   cd ./04_Annotations/antismash/${sample}
#   # echo ${sample}
#     for genbank_file in $(ls c*.gbk); do
#       cluster=$(echo ${genbank_file} | cut -f1 -d '_')
#       echo copying and renaming ${genbank_file} to cluster_${cluster}_${sample}.gbk
#       cp ${genbank_file} ${project_dir}/04_Annotations/bigscape/cluster_${cluster}_${sample}.gbk
#     done
#   cd ${project_dir}
# done
#
# conda activate bigscape
#
# cd ${project_dir}
#
# bigscape.py -i ./04_Annotations/bigscape/ -o ./04_Annotations/bigscape/results >> ./04_Annotations/bigscape/${project}.log 2&>1
#
# conda deactivate

echo finished meow
