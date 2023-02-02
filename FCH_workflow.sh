#!/bin/bash
# Copyright (c) 2022, René KM Xavier
set -e
# setting to break on newline only
IFS=$'\n'


set -o errtrace # Enable the err trap, code will get called when an error is detected
trap "echo ERROR: There was an error check error log" ERR
###################
#### VARIABLES ####
###################
read -r -p "What directory are the raw sequencing reads in? Please give the absolute path to 00_Raw_Reads without a trailing '/'  " reads_dir || exit 100
read -r -p "Please provide the absolute path to the project working directory:     " project_dir || exit 100
read -r -p "What is your project name? No spaces or special characters please.   " project || exit 100
read -r -p "How much RAM do you have available in gigabytes (Gb)? Use digits only.   " RAM || exit 100
read -r -p "How many cores do you want to use? Use digits only   " CPU || exit 100
read -r -p "Create a samples.txt file with unique sample ids without spaces, saved in the same directory as the raw sequencing reads. Do you have the samples.txt file made? yes or no:   " samples || exit 100
read -r -p "Do you want to count the raw sequencing reads?  yes or no:   " count || exit 100
read -r -p "Do you want to preprocess your raw sequencing reads? yes or no:   " reads_QC || exit 100
read -r -p "Do you want to analyze your HQ reads for Eukaryotes? yes or no:  " taxatarget || exit 100
read -r -p "Do you want to analyze your HQ reads for trans-kingdom taxonomy? yes or no :" metaphlan || exit 100
# Be sure samples.txt is formatted correctly.
if [[ "$samples" == no ]]
then
  echo Change into the 00_Raw_Reads directory and in the terminal run: "ls ./*_R1.f*q.gz | cut -f1 -d '.' | rev | cut -f2- -d '_' | rev | sort -u >> samples.txt"
  read -r -p "Create a samples.txt file with unique sample ids without spaces, saved in the same directory as the raw sequencing reads. Do you have the samples.txt file made? yes or no:   " samples || exit 100
fi
# Make and Change into new project directory
mkdir -p $project_dir
cd ${project_dir}
###################
## PREPROCESSING ##
###################
mkdir -p ./00_Raw_Reads/fastqc
mkdir -p ./01_HQ_Reads/merged_reads
mkdir -p ./01_HQ_Reads/error_corrected
mkdir -p ./01_HQ_Reads/fastqc
mkdir -p ./01_HQ_Reads/logs
# Copy file of unique sample ids into the project directory
cp ${reads_dir}/samples.txt .
# Get a quick look at the size of each sample.
#read -r -p "Do you want to count the raw sequencing reads?  yes or no:   " count || exit 100
if [[ "${count}" == yes ]]; then
  while read -r sample; do
    echo Counting all the sequences in $sample
    reformat.sh in=${reads_dir}/${sample}_R1.fq.gz in2=${reads_dir}/${sample}_R2.fq.gz > ${reads_dir}/${sample}_counts.txt 2>&1
    counts=$(grep 'Output:' ${reads_dir}/${sample}_counts.txt | cut -f2,3)
    echo -e "${sample}\t${counts}" >> ${project_dir}/00_Raw_Reads_counts.tsv
  done < samples.txt
elif [[ "${count}" == no ]]; then :
fi
############################
### READ QUALITY CONTROL ###
############################
#The following preprocesing workflow, creates high-quality trimmed, error corrected, and merged reads.
###############
### CONDA #####
###############
# Install bbmap, fastqc, and multiqc in a conda environment:
# conda create -y -n pp bbmap fastqc multiqc
#read -r -p "Do you want to preprocess your raw sequencing reads? yes or no:   " reads_QC || exit 100
if [[ "${reads_QC}" == yes ]]
then

  conda activate pp
  while read -r sample; do
    if test -f "./01_HQ_Reads/${sample}_HQ.fq"; then :
    else
      #Input raw reads from 00_Raw_Reads directory

      R1=( ${reads_dir}/${sample}*1.f*q.gz )
      R2=( ${reads_dir}/${sample}*2.f*q.gz )

      # echo Assessing intial read quality on ${sample} raw reads in ${reads_dir}
      # fastqc $R1 $R2 -q -o ./00_Raw_Reads/fastqc
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
      bbduk.sh in=./01_HQ_Reads/atrimmed.fq.gz out=./01_HQ_Reads/qtrimmed.fq.gz qtrim=rl trimq=6 ftm=5 trimpolygright=0 minlen=70 ordered maxns=0 maq=8 entropy=.95 ow=t >> ./01_HQ_Reads/logs/${sample}.log 2>&1

      ########## NOT RECOMMENDED ##############
      #ftl=16 to remove the 16bp of barcode left on the 5' read. Removes Nextera 5' noise
      #bbduk.sh in=./01_HQ_Reads/atrimmed.fq.gz out=./01_HQ_Reads/qtrimmed.fq.gz ftl=16 qtrim=rl trimq=6 ftm=5 trimpolygright=0 minlen=70 ordered maxns=0 maq=8 entropy=.95 ow=t >> ./01_HQ_Reads/logs/${sample}.log 2>&1
      #########################################

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

      #Create R1.fa for nonpareil and phyloflash input
      reformat.sh in=./01_HQ_Reads/error_corrected/${sample}_HQ_ecc.fq.gz out=./01_HQ_Reads/error_corrected/${sample}_HQ_ecc_R1.fa out2=./01_HQ_Reads/error_corrected/${sample}_HQ_ecc_R2.fa verifypaired ow=t >> ./01_HQ_Reads/logs/${sample}.log 2>&1
      #Create fastq R1 and R2 input files for TaxaTarget
      reformat.sh in=./01_HQ_Reads/${sample}_HQ.fq.gz out=./01_HQ_Reads/${sample}_HQ_R1.fastq out2=./01_HQ_Reads/${sample}_HQ_R2.fastq verifypaired ow=t >> ./01_HQ_Reads/logs/${sample}.log 2>&1
      # Create interleaved fastq file for MetaPhlAn
      gunzip -k ./01_HQ_Reads/${sample}_HQ.fq.gz
      echo finished ${sample}

    fi
  done < samples.txt

  if test -f "./01_HQ_Reads/multiqc_report.html"; then :
  else
  #Run multiqc on all clean Paired End reads.
    echo Assessing quality of processed ${project} reads.
    find ./01_HQ_Reads/ -name "*_HQ.fq.gz" -exec fastqc '{}' -t 6 -q -o ./01_HQ_Reads/fastqc \;
    multiqc ./01_HQ_Reads/fastqc -o ./01_HQ_Reads/
  fi

  conda deactivate

  #########################################################################################################################################
  ##############################                READ-BASED ANALYSIS                                           ##############################
  #########################################################################################################################################

  echo Finished ${sample} read quality control. Starting analysis of sequencing effort.

  conda activate nonpareil
  mkdir -p ./01_HQ_Reads/nonpareil
  while read -r sample; do
    if test -f "./01_HQ_Reads/nonpareil/${sample}.log"; then :
    else
      echo Estimating sequencing effort for ${sample}
      nonpareil -s ./01_HQ_Reads/error_corrected/${sample}_HQ_ecc_R1.fa -T alignment -f fasta -t 12 -b ./01_HQ_Reads/nonpareil/${sample} >> ./01_HQ_Reads/nonpareil/${sample}.log 2>&1
      #nonpareil -s ./01_HQ_Reads/error_corrected/${sample}_HQ_ecc_R1.fa -T alignment -L 35 -t ${CPU} -R ${RAM}000 -b ./01_HQ_Reads/nonpareil/${sample} >> ./01_HQ_Reads/nonpareil/${sample}.log 2>&1
    fi
  done < samples.txt

  conda deactivate
  ###########################
  ##  READ-BASED TAXONOMY  ##
  ###########################
  echo Starting assembly of rRNA SSU and classification

  conda activate phyloflash
  mkdir -p ./02_Taxonomy/phyloflash
  while read -r sample; do
    if test -f "${project_dir}/02_Taxonomy/phyloflash/${sample}.phyloFlash.log"; then :
    else
      cd ${project_dir}/02_Taxonomy/phyloflash
      echo Mining ${sample} metagenome for rRNA genes
      phyloFlash.pl -lib ${sample} -read1 ${project_dir}/01_HQ_Reads/error_corrected/${sample}_HQ_ecc_R1.fa -read2 ${project_dir}/01_HQ_Reads/error_corrected/${sample}_HQ_ecc_R2.fa -CPU ${CPU} -almosteverything -taxlevel 20 -log
    fi
  done < samples.txt
  conda deactivate
elif [[ "${reads_QC}" == no ]]; then :
fi

# Install taxaTarget and dependencies following instructions at: https://github.com/SethCommichaux/taxaTarget
#read -r -p "Do you want to analyze your HQ reads for Eukaryotes? yes or no:  " taxatarget || exit 100
if [[ "${taxatarget}" == no ]]; then :
elif [[ "${taxatarget}" == yes ]]; then
  #statements
  mkdir -p ./02_Taxonomy/taxatarget

  source activate taxatarget

  # Run taxaTarget on sample reads
  for sample in `cat samples.txt`; do
    if test -f "./02_Taxonomy/taxatarget/${sample}_results/Taxonomic_report.txt"; then :
    else

      echo Analyzing ${sample} sample for Eukaryotes
      HQ1=./01_HQ_Reads/${sample}_HQ_R1.fastq
      HQ2=./01_HQ_Reads/${sample}_HQ_R2.fastq

      python /home/fch/GitHub/taxaTarget/run_pipeline_scripts/run_protist_pipeline_fda.py -r ${HQ1} -r2 ${HQ2} -e /home/fch/GitHub/taxaTarget/run_pipeline_scripts/environment.txt -o ./02_Taxonomy/taxatarget/${sample}_results -t ${CPU} >> ./02_Taxonomy/taxatarget/${sample}.log 2>&1
      python /home/fch/bin/modify_taxonomic_report.py ./02_Taxonomy/taxatarget/${sample}_results/Taxonomic_report.txt
    fi
  done

  if test -f "${project_dir}/02_Taxonomy/taxatarget/files.txt"; then :
  else
    # Data Wrangling step. Create input file for aggregate_results.py
    ls ${project_dir}/02_Taxonomy/taxatarget/*_results/Taxonomic_report.txt.modified > tmp.txt
    for line in `cat tmp.txt`; do
      sample=$(echo ${line} | rev | cut -f2 -d '/' | rev | cut -f1 -d '_')
      echo -e "${sample}\t${line}" >> ${project_dir}/02_Taxonomy/taxatarget/files.txt
    done
    rm tmp.txt
    # Aggregate TaxaTarget results into one matrix for statistical analysis
    cd ${project_dir}/02_Taxonomy/taxatarget
    python /home/fch/bin/aggregate_results.py files.txt
    cd ${project_dir}
  fi
  conda deactivate
fi

cd ${project_dir}
# Install metaphlan into conda environment
# conda create -n metaphlan -y -c bioconda python=3.7 metaphlan
#activate metaphlan conda environment
# conda activate metaphlan
# Install metaphlan databases
# metaphlan --install --bowtie2db <database folder>
#read -r -p "Do you want to analyze your HQ reads for trans-kingdom taxonomy? yes or no :" metaphlan || exit 100
# MetaPhlAn introduces the capability of characterizing organisms at the strain level using non
# aggregated marker information. Such capability comes with several slightly different flavours and
# are a way to perform strain tracking and comparison across multiple samples.
# Usually, MetaPhlAn is first ran with the default -t to profile the species present in
# the community, and then a strain-level profiling can be performed to zoom-in into specific species
# of interest. This operation can be performed quickly as it exploits the --bowtie2out intermediate
# file saved during the execution of the default analysis type.
if [[ "${metaphlan}" == no ]]; then :
elif [[ "${metaphlan}" == yes ]]; then
  conda activate metaphlan
  while read -r sample; do
    if test -f "./02_Taxonomy/metaphlan/${sample}/${sample}_viral_profile.txt"; then :
    else
      mkdir -p ./02_Taxonomy/metaphlan/${sample}
      echo Creating trans-kingdom taxonomic profile for ${sample}
      #metaphlan ./01_HQ_Reads/${sample}_HQ.fq --input_type fastq --add_viruses --bowtie2out ./02_Taxonomy/metaphlan/${sample}/${sample}_viral_bowtie2.bz2 --nproc ${CPU} -o ./02_Taxonomy/metaphlan/${sample}/${sample}_viral_profile.txt
      metaphlan ./01_HQ_Reads/${sample}_HQ.fq --input_type fastq --unknown_estimation --add_viruses --bowtie2out ./02_Taxonomy/metaphlan/${sample}/${sample}_viral_ukn_bowtie2.bz2 --nproc ${CPU} -o ./02_Taxonomy/metaphlan/${sample}/${sample}_viral_ukn_profile.txt
      #metaphlan ./01_HQ_Reads/${sample}_HQ.fq --input_type fastq --bowtie2out ./02_Taxonomy/metaphlan/${sample}/${sample}_bowtie2.bz2 --nproc ${CPU} -o ./02_Taxonomy/metaphlan/${sample}/${sample}_profile.txt
      #metaphlan ./01_HQ_Reads/${sample}_HQ.fq --input_type fastq --unknown_estimation --bowtie2out ./02_Taxonomy/metaphlan/${sample}/${sample}_ukn_bowtie2.bz2 --nproc ${CPU} -o ./02_Taxonomy/metaphlan/${sample}/${sample}_ukn_profile.txt
    fi
  done < samples.txt

  if test -f "./02_Taxonomy/metaphlan/profiles_merged/${project}_viral_profile.txt"; then :
  else
    #Make a matrix of taxonomy data
    mkdir -p ./02_Taxonomy/metaphlan/profiles_merged
    #merge_metaphlan_tables.py ./02_Taxonomy/metaphlan/*/*_viral_profile.txt > ./02_Taxonomy/metaphlan/profiles_merged/${project}_viral_profile.txt
    merge_metaphlan_tables.py ./02_Taxonomy/metaphlan/*/*_viral_ukn_profile.txt > ./02_Taxonomy/metaphlan/profiles_merged/${project}_viral_ukn_profile.txt
    #merge_metaphlan_tables.py ./02_Taxonomy/metaphlan/*/*_profile.txt > ./02_Taxonomy/metaphlan/profiles_merged/${project}_profile.txt
    #merge_metaphlan_tables.py ./02_Taxonomy/metaphlan/*/*_ukn_profile.txt > ./02_Taxonomy/metaphlan/profiles_merged/${project}_ukn_profile.txt
  fi
  conda deactivate
fi
#########################################################################################################################################
##############################                ASSEMBLY-BASED ANALYSIS                                           #########################
#########################################################################################################################################

# The following assembly workflow creates assemblies for three debruijn graph assemblers and compares them for each sample.
###############
### CONDA #####
###############
# conda create -n assembly -y -c bioconda spades megahit idba
# conda activate assembly
# pip install biopython
read -r -p "Would you like to assemble your preprocessed sample reads?  yes or no:  " assemble || exit 100
if [[ "${assemble}" == yes ]]; then
  #Create directories
  mkdir -p ./02_Assembly/metaspades/logs
  mkdir -p ./02_Assembly/megahit-large/logs
  mkdir -p ./02_Assembly/idba_ud/logs
  mkdir -p ./02_Assembly/stats
  ## Do not make a directory for megahit, it will make it itself and will fail if the output directory already exists, as to not overwrite results. ##
  ## Since megahit will not write to external harddrives, the code is written to write data to home directory then move it to the working directory.
  ################
  ### Assembly ###
  ################
  read -r -p "Which HQ reads dataset do you want to use to assemble? HQ HQ_ecc HQ_merged Skip:  " HQ_Reads || exit 100
  read -r -p "Would you like to evaluate the assemblies? yes or no:  " quast || exit 100

  for sample in `cat samples.txt`; do

    mkdir -p ./02_Assembly/quast/${sample}
    conda activate assembly
    if [[ "${HQ_Reads}" == HQ ]]
    then

      if test -f "./02_Assembly/metaspades/${sample}/scaffolds.fasta"; then :
      else

        echo Using metaspades to assemble ${HQ_Reads} ${sample} reads.
        mkdir -p ./02_Assembly/metaspades/${sample}
        metaspades.py --only-assembler -k 21,33,55,77,99,127 --12 ./01_HQ_Reads/${sample}_HQ.fq.gz -t ${CPU} -m ${RAM} -o ./02_Assembly/metaspades/${sample} >> ./02_Assembly/metaspades/logs/${sample}.log 2>&1
        cp ./02_Assembly/metaspades/${sample}/scaffolds.fasta ./02_Assembly/quast/${sample}/metaspades.fasta
      fi

      if test -f "./02_Assembly/quast/${sample}/megahit.fasta"; then :
      else
        #Assemble reads per site using megahit
        echo Using megahit to assemble ${HQ_Reads} ${sample} reads.
        megahit --12 ./01_HQ_Reads/${sample}_HQ.fq.gz --presets meta-large -m 0.5 -t ${CPU} -o ~/${sample} >> ./02_Assembly/megahit-large/logs/${sample}.log 2>&1
        cp  ~/${sample}/final.contigs.fa ./02_Assembly/quast/${sample}/megahit.fasta
        mkdir -p ./02_Assembly/megahit-large/${sample}
        mv ~/${sample} ./02_Assembly/megahit-large
      fi

      if test -f "./01_HQ_Reads/${sample}_HQ.fa"; then :
      else
        #IDBA can only take fasta files as input; therefore, fastq.gz files must be reformated to fasta files
        reformat.sh in=./01_HQ_Reads/${sample}_HQ.fq.gz out=./01_HQ_Reads/${sample}_HQ.fa ow=t >> ./02_Assembly/idba_ud/logs/${sample}.log 2>&1
        #fq2fa --paired ./01_HQ_Reads/${sample}_HQ.fq.gz ./01_HQ_Reads/${sample}_HQ.fa
      fi

      if test -f "./02_Assembly/idba_ud/${sample}/scaffold.fa"; then :
      else
        # #Assemble reads per site using idba_ud
        echo Using idba to assemble ${HQ_Reads} ${sample} reads.
        mkdir -p ./02_Assembly/idba_ud/${sample}
        idba_ud -r ./01_HQ_Reads/${sample}_HQ.fa --num_threads=${CPU} -o ./02_Assembly/idba_ud/${sample} >> ./02_Assembly/idba_ud/logs/${sample}.log 2>&1
        cp ./02_Assembly/idba_ud/${sample}/scaffold.fa ./02_Assembly/quast/${sample}/idba.fasta
      fi

      if grep -m 1 -q "${HQ_Reads}" ./02_Assembly/stats/${sample}_assembly_summary.tsv; then :
      else
        # --- Evaluation ---
        echo -e "metagenome_id\treads\tassembler\tn_contigs\tcontig_bp\tgap_pct\tctg_L50\tctg_max\tgene_counts\tmapped_pct\tavg_cov\tassembly_performance" >> ./02_Assembly/stats/${sample}_assembly_summary.tsv
        echo Evaluating ${sample} assemblies with AssemblyStats
        ##statswrapper.sh ./02_Assembly/quast/${sample}/*.fasta minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/stats/${sample}.txt
        stats.sh ./02_Assembly/metaspades/${sample}/scaffolds.fasta minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/metaspades/${sample}/stats.txt
        spades_stats=$(cat ./02_Assembly/metaspades/${sample}/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d')
        stats.sh ./02_Assembly/megahit-large/${sample}/final.contigs.fa minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/megahit-large/${sample}/stats.txt
        megahit_stats=$(cat ./02_Assembly/megahit-large/${sample}/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d')
        stats.sh ./02_Assembly/idba_ud/${sample}/scaffold.fa minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/idba_ud/${sample}/stats.txt
        idba_stats=$(cat ./02_Assembly/idba_ud/${sample}/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d')

        spades_L50=$(cat ./02_Assembly/metaspades/${sample}/stats.txt | cut -f9 -d$'\t' | sed '2q;d')
        megahit_L50=$(cat ./02_Assembly/megahit-large/${sample}/stats.txt | cut -f9 -d$'\t' | sed '2q;d')
        idba_L50=$(cat ./02_Assembly/idba_ud/${sample}/stats.txt | cut -f9 -d$'\t' | sed '2q;d')

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
        spades_match=$(grep Percent\ mapped: ./02_Assembly/metaspades/${sample}/map.log | cut -f2 -d$'\t' | tail -1)
        spades_cov=$(grep Average\ coverage: ./02_Assembly/metaspades/${sample}/map.log | cut -f2 -d$'\t' | tail -1)
        spades_Pa=$(echo "$spades_match*$spades_L50" | bc)

        bbmap.sh in=./01_HQ_Reads/${sample}_HQ.fq.gz ref=./02_Assembly/megahit-large/${sample}/${sample}_c1k.fa nodisk covhist=./02_Assembly/megahit-large/${sample}/covhist.txt covstats=./02_Assembly/megahit-large/${sample}/covstats.txt outm=./02_Assembly/megahit-large/${sample}/${sample}_reads_assembled.fq.gz outu=./02_Assembly/megahit-large/${sample}/${sample}_reads_unassembled.fq.gz fast=t ambig=best >> ./02_Assembly/megahit-large/${sample}/map.log 2>&1
        megahit_match=$(grep Percent\ mapped: ./02_Assembly/megahit-large/${sample}/map.log | cut -f2 -d$'\t' | tail -1)
        megahit_cov=$(grep Average\ coverage: ./02_Assembly/megahit-large/${sample}/map.log | cut -f2 -d$'\t' | tail -1)
        megahit_Pa=$(echo "$megahit_match*$megahit_L50" | bc)

        bbmap.sh in=./01_HQ_Reads/${sample}_HQ.fq.gz ref=./02_Assembly/idba_ud/${sample}/${sample}_c1k.fa nodisk covhist=./02_Assembly/idba_ud/${sample}/covhist.txt covstats=./02_Assembly/idba_ud/${sample}/covstats.txt outm=./02_Assembly/idba_ud/${sample}/${sample}_reads_assembled.fq.gz outu=./02_Assembly/idba_ud/${sample}/${sample}_reads_unassembled.fq.gz fast=t ambig=best >> ./02_Assembly/idba_ud/${sample}/map.log 2>&1
        idba_match=$(grep Percent\ mapped: ./02_Assembly/idba_ud/${sample}/map.log | cut -f2 -d$'\t' | tail -1)
        idba_cov=$(grep Average\ coverage: ./02_Assembly/idba_ud/${sample}/map.log | cut -f2 -d$'\t' | tail -1)
        idba_Pa=$(echo "$idba_match*$idba_L50" | bc)


        echo Creating assembly comparision for ${sample}
        echo -e "${sample}\t${HQ_Reads}\tmetaspades\t${spades_stats}\t${spades_counts}\t${spades_match}\t${spades_cov}\t${spades_Pa}" >> ./02_Assembly/stats/${sample}_assembly_summary.tsv
        echo -e "${sample}\t${HQ_Reads}\tmegahit\t${megahit_stats}\t${megahit_counts}\t${megahit_match}\t${megahit_cov}\t${megahit_Pa}" >> ./02_Assembly/stats/${sample}_assembly_summary.tsv
        echo -e "${sample}\t${HQ_Reads}\tidba_ud\t${idba_stats}\t${idba_counts}\t${idba_match}\t${idba_cov}\t${idba_Pa}" >> ./02_Assembly/stats/${sample}_assembly_summary.tsv

      fi

    elif [[ "${HQ_Reads}" == HQ_ecc ]]
    then

      if test -f "./02_Assembly/metaspades/${sample}_ecc/scaffolds.fasta"; then :
      else
        echo Using metaspades to assemble ${HQ_Reads} ${sample} reads.
        mkdir -p ./02_Assembly/metaspades/${sample}_ecc
        metaspades.py -k 21,33,55,77,99,127 --only-assembler --12 ./01_HQ_Reads/error_corrected/${sample}_HQ_ecc.fq.gz -t ${CPU} -m ${RAM} -o ./02_Assembly/metaspades/${sample}_ecc >> ./02_Assembly/metaspades/logs/${sample}.log 2>&1
        cp ./02_Assembly/metaspades/${sample}_ecc/scaffolds.fasta ./02_Assembly/quast/${sample}/metaspades_ecc.fasta
      fi

      if test -f "./02_Assembly/quast/${sample}/megahit_ecc.fasta"; then :
      else

        #Assemble reads per site using megahit
        echo Using megahit to assemble ${HQ_Reads} ${sample} reads.

        megahit --12 ./01_HQ_Reads/error_corrected/${sample}_HQ_ecc.fq.gz --presets meta-large -m 0.5 -t ${CPU} -o ./02_Assembly/megahit-large/${sample} >> ./02_Assembly/megahit-large/logs/${sample}.log 2>&1
        cp ~/${sample}_ecc/final.contigs.fa ./02_Assembly/quast/${sample}/megahit_ecc.fasta
        mkdir -p ./02_Assembly/megahit-large/${sample}_ecc
        mv ~/${sample}_ecc/ ./02_Assembly/megahit-large
      fi

      if test -f "./01_HQ_Reads/error_corrected/${sample}_HQ_ecc.fa"; then :
      else

        #IDBA can only take fasta files as input; therefore, fastq.gz files must be reformated to fasta files
        reformat.sh in=./01_HQ_Reads/error_corrected/${sample}_HQ_ecc.fq.gz out=./01_HQ_Reads/error_corrected/${sample}_HQ_ecc.fa ow=t >> ./02_Assembly/idba_ud/logs/${sample}.log 2>&1
        #fq2fa --paired ./01_HQ_Reads/error_corrected/${sample}_HQ_ecc.fq.gz ./01_HQ_Reads/error_corrected/${sample}_HQ_ecc.fa

      fi
      if test -f "./02_Assembly/idba_ud/${sample}_ecc/scaffold.fa"; then :
      else

        # #Assemble reads per site using idba_ud
        echo Using idba to assemble ${HQ_Reads} ${sample} reads.

        mkdir -p ./02_Assembly/idba_ud/${sample}_ecc
        idba_ud -r ./01_HQ_Reads/error_corrected/${sample}_HQ_ecc.fa --num_threads=${CPU} --pre_correction -o ./02_Assembly/idba_ud/${sample}_ecc >> ./02_Assembly/idba_ud/logs/${sample}.log 2>&1
        cp ./02_Assembly/idba_ud/${sample}_ecc/scaffold.fa ./02_Assembly/quast/${sample}/idba_ecc.fasta
      fi

      if grep -m 1 -q "${HQ_Reads}" ./02_Assembly/stats/${sample}_assembly_summary.tsv; then :
      else

        # --- Evaluation ---
        echo -e "metagenome_id\treads\tassembler\tn_contigs\tcontig_bp\tgap_pct\tctg_L50\tctg_max\tgene_counts\tmapped_pct\tavg_cov\tassembly_performance" >> ./02_Assembly/stats/${sample}_assembly_summary.tsv
        echo Evaluating ${sample} assemblies with AssemblyStats
        stats.sh ./02_Assembly/metaspades/${sample}_ecc/scaffolds.fasta minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/metaspades/${sample}_ecc/stats.txt
        spades_stats=$(cat ./02_Assembly/metaspades/${sample}_ecc/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d')
        stats.sh ./02_Assembly/megahit-large/${sample}_ecc/final.contigs.fa minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/megahit-large/${sample}_ecc/stats.txt
        megahit_stats=$(cat ./02_Assembly/megahit-large/${sample}_ecc/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d')
        stats.sh ./02_Assembly/idba_ud/${sample}_ecc/scaffold.fa minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/idba_ud/${sample}_ecc/stats.txt
        idba_stats=$(cat ./02_Assembly/idba_ud/${sample}_ecc/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d')

        spades_L50=$(cat ./02_Assembly/metaspades/${sample}_ecc/stats.txt | cut -f9 -d$'\t' | sed '2q;d')
        megahit_L50=$(cat ./02_Assembly/megahit-large/${sample}_ecc/stats.txt | cut -f9 -d$'\t' | sed '2q;d')
        idba_L50=$(cat ./02_Assembly/idba_ud/${sample}_ecc/stats.txt | cut -f9 -d$'\t' | sed '2q;d')

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
        # Note that these are reads mapped to contigs > 1,000 bp.
        #(i.e. % mapped) by mapping reads by to assembly
        bbmap.sh in=./01_HQ_Reads/${sample}_HQ.fq.gz ref=./02_Assembly/metaspades/${sample}_ecc/${sample}_c1k.fa nodisk covhist=./02_Assembly/metaspades/${sample}_ecc/covhist.txt covstats=./02_Assembly/metaspades/${sample}_ecc/covstats.txt out=./02_Assembly/metaspades/${sample}_ecc/${sample}_reads_assembled.fq.gz outu=./02_Assembly/metaspades/${sample}_ecc/${sample}_reads_unassembled.fq.gz fast=t ambig=best >> ./02_Assembly/metaspades/${sample}_ecc/map.log 2>&1
        spades_match=$(grep Percent\ mapped: ./02_Assembly/metaspades/${sample}_ecc/map.log | cut -f2 -d$'\t' | tail -1)
        spades_cov=$(grep Average\ coverage: ./02_Assembly/metaspades/${sample}_ecc/map.log | cut -f2 -d$'\t' | tail -1)
        spades_Pa=$(echo "$spades_match*$spades_L50" | bc)

        bbmap.sh in=./01_HQ_Reads/${sample}_HQ.fq.gz ref=./02_Assembly/megahit-large/${sample}_ecc/${sample}_c1k.fa nodisk covhist=./02_Assembly/megahit-large/${sample}_ecc/covhist.txt covstats=./02_Assembly/megahit-large/${sample}_ecc/covstats.txt outm=./02_Assembly/megahit-large/${sample}_ecc/${sample}_reads_assembled.fq.gz outu=./02_Assembly/megahit-large/${sample}_ecc/${sample}_reads_unassembled.fq.gz fast=t ambig=best >> ./02_Assembly/megahit-large/${sample}_ecc/map.log 2>&1
        megahit_match=$(grep Percent\ mapped: ./02_Assembly/megahit-large/${sample}_ecc/map.log | cut -f2 -d$'\t' | tail -1)
        megahit_cov=$(grep Average\ coverage: ./02_Assembly/megahit-large/${sample}_ecc/map.log | cut -f2 -d$'\t' | tail -1)
        megahit_Pa=$(echo "$megahit_match*$megahit_L50" | bc)

        bbmap.sh in=./01_HQ_Reads/${sample}_HQ.fq.gz ref=./02_Assembly/idba_ud/${sample}_ecc/${sample}_c1k.fa nodisk covhist=./02_Assembly/idba_ud/${sample}_ecc/covhist.txt covstats=./02_Assembly/idba_ud/${sample}_ecc/covstats.txt outm=./02_Assembly/idba_ud/${sample}_ecc/${sample}_reads_assembled.fq.gz outu=./02_Assembly/idba_ud/${sample}_ecc/${sample}_reads_unassembled.fq.gz fast=t ambig=best >> ./02_Assembly/idba_ud/${sample}_ecc/map.log 2>&1
        idba_match=$(grep Percent\ mapped: ./02_Assembly/idba_ud/${sample}_ecc/map.log | cut -f2 -d$'\t' | tail -1)
        idba_cov=$(grep Average\ coverage: ./02_Assembly/idba_ud/${sample}_ecc/map.log | cut -f2 -d$'\t' | tail -1)
        idba_Pa=$(echo "$idba_match*$idba_L50" | bc)

        echo Creating assembly summary for ${sample}
        echo -e "${sample}\t${HQ_Reads}\tmetaspades\t${spades_stats}\t${spades_counts}\t${spades_match}\t${spades_cov}\t${spades_Pa}" >> ./02_Assembly/stats/${sample}_assembly_summary.tsv
        echo -e "${sample}\t${HQ_Reads}\tmegahit\t${megahit_stats}\t${megahit_counts}\t${megahit_match}\t${megahit_cov}\t${megahit_Pa}" >> ./02_Assembly/stats/${sample}_assembly_summary.tsv
        echo -e "${sample}\t${HQ_Reads}\tidba_ud\t${idba_stats}\t${idba_counts}\t${idba_match}\t${idba_cov}\t${idba_Pa}" >> ./02_Assembly/stats/${sample}_assembly_summary.tsv
      fi

    elif [[ "${HQ_Reads}" == HQ_merged ]]
    then

      if test -f "./02_Assembly/metaspades/${sample}_merged/scaffolds.fasta"; then :
      else

        echo Using metaspades to assemble ${HQ_Reads} ${sample} reads.
        # Assemble using merged/unmerged reads
        mkdir -p ./02_Assembly/metaspades/${sample}_merged
        metaspades.py  --only-assembler -k 21,33,55,77,99,127 --merged ./01_HQ_Reads/merged_reads/${sample}_HQ_merged.fq.gz --12 ./01_HQ_Reads/merged_reads/${sample}_HQ_unmerged.fq.gz -t ${CPU} -m ${RAM} -o ./02_Assembly/metaspades/${sample}_merged >> ./02_Assembly/metaspades/logs/${sample}.log 2>&1
        cp ./02_Assembly/metaspades/${sample}_merged/scaffolds.fasta ./02_Assembly/quast/${sample}/metaspades_merged.fasta
      fi

      if test -f "./02_Assembly/quast/${sample}/megahit_merged.fasta"; then :
      else

        echo Using megahit to assemble ${HQ_Reads} ${sample} reads.
        megahit -r ./01_HQ_Reads/merged_reads/${sample}_HQ_merged.fq.gz --12 ./01_HQ_Reads/merged_reads/${sample}_HQ_unmerged.fq.gz --presets meta-large -m 0.5 -t ${CPU} -o ./02_Assembly/megahit-large/${sample} >> ./02_Assembly/megahit-large/logs/${sample}.log 2>&1
        cp ~/${sample}_merged/final.contigs.fa ./02_Assembly/quast/${sample}/megahit_merged.fasta
        mkdir -p ./02_Assembly/megahit-large/${sample}_merged/
        mv ~/${sample}_merged/ ./02_Assembly/megahit-large/
        echo idba_ud does not assemble merged reads and will not be used for comparison.
      fi
      #########
      ## FIX ##
      #########
      ## Fix code so 'No such file or directory' Error does not show
      if grep -m 1 -q "${HQ_Reads}" ./02_Assembly/stats/${sample}_assembly_summary.tsv; then :
      else

        echo -e "metagenome_id\treads\tassembler\tn_contigs\tcontig_bp\tgap_pct\tctg_L50\tctg_max\tgene_counts\tmapped_pct\tavg_cov\tassembly_performance" >> ./02_Assembly/stats/${sample}_assembly_summary.tsv
        echo Evaluating ${sample} assemblies with AssemblyStats
        stats.sh ./02_Assembly/metaspades/${sample}_merged/scaffolds.fasta minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/metaspades/${sample}_merged/stats.txt
        spades_stats=$(cat ./02_Assembly/metaspades/${sample}_merged/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d')
        stats.sh ./02_Assembly/megahit-large/${sample}_merged/final.contigs.fa minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/megahit-large/${sample}_merged/stats.txt
        megahit_stats=$(cat ./02_Assembly/megahit-large/${sample}_merged/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d')

        spades_L50=$(cat ./02_Assembly/metaspades/${sample}_merged/stats.txt | cut -f9 -d$'\t' | sed '2q;d')
        megahit_L50=$(cat ./02_Assembly/megahit-large/${sample}_meged/stats.txt | cut -f9 -d$'\t' | sed '2q;d')

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

        spades_Pa=$(echo "$spades_match*$spades_L50" | bc)

        bbmap.sh in=./01_HQ_Reads/${sample}_HQ.fq.gz ref=./02_Assembly/megahit-large/${sample}_merged/${sample}_c1k.fa nodisk covhist=./02_Assembly/megahit-large/${sample}_merged/covhist.txt covstats=./02_Assembly/megahit-large/${sample}_merged/covstats.txt outm=./02_Assembly/megahit-large/${sample}_merged/${sample}_reads_assembled.fq.gz outu=./02_Assembly/megahit-large/${sample}_merged/${sample}_reads_unassembled.fq.gz fast=t ambig=best >> ./02_Assembly/megahit-large/${sample}_merged/map.log 2>&1
        megahit_match=$(grep Percent\ mapped: ./02_Assembly/megahit-large/${sample}_merged/map.log | cut -f2 -d$'\t' | tail -1)
        megahit_cov=$(grep Average\ coverage: ./02_Assembly/megahit-large/${sample}_merged/map.log | cut -f2 -d$'\t' | tail -1)
        metahit_Pa=$(echo "$megahit_match*$megahit_L50" | bc)

        echo Creating assembly comparision for ${sample}
        echo -e "${sample}\t${HQ_Reads}\tmetaspades\t${spades_stats}\t${spades_counts}\t${spades_match}\t${spades_cov}\t${spades_Pa}" >> ./02_Assembly/stats/${sample}_assembly_summary.tsv
        echo -e "${sample}\t${HQ_Reads}\tmegahit\t${megahit_stats}\t${megahit_counts}\t${megahit_match}\t${megahit_cov}\t${megahit_Pa}" >> ./02_Assembly/stats/${sample}_assembly_summary.tsv
      fi

    elif [[ "${HQ_Reads}" == Skip ]]; then :
    fi
    ###################
    ### Assembly QC ###
    ###################
    #read -r -p "Would you like to evaluate the ${sample} assemblies? yes or no:  " quast || exit 100
    if [[ "${quast}" == yes ]]; then

      if test -f "./02_Assembly/quast/${sample}.log"; then :
      else
      # --- Evaluation ---
        echo Evaluating ${sample} assemblies
        #statswrapper.sh ./02_Assembly/quast/${sample}/*.fasta minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/stats/${sample}.txt
        ~/GitHub/quast/metaquast.py -f -m 1000 -o ./02_Assembly/quast/${sample} ./02_Assembly/quast/${sample}/*.fasta >> ./02_Assembly/quast/${sample}.log 2>&1
      fi

    elif [[ "${quast}" == no ]]; then :
    fi

    cut -f2,3,5,12 -d$'\t' ./02_Assembly/stats/${sample}_assembly_summary.tsv
    ls `pwd`/02_Assembly/quast/${sample}/*.fasta
    ## Future, would be nice to also include the paths to each. So don't have to leave the command line.
    read -r -p "Analyze the ./02_Assembly/stats/${sample}_assembly_summary.tsv and ./02_Assembly/quast/${sample}/report.html and determine which ${sample} assembly is best. Please enter the absoute path to the assembly.fasta file:  " best_assembly || exit 100
    read -r -p "Would you like to estimate the taxonomy for the ${sample} metagenome? yes or no: " taxa_1 || exit 100
    read -r -p "Would you like to analyze the ${sample} assembly using anvi'o? yes or no: " anvio_1 || exit 100
    read -r -p "Would you like to annotate your ${sample} metagenome for biosynthetic gene clusters? yes or no: " smash || exit 100
    read -r -p "What minimum contig length do you want to analyze for biosynthetic gene clusters? or type NA   " min_BGC || exit 100
    read -r -p "Would you like to analyze annotated biosynthetic gene clusters for each sample against the MiBIG database?  yes or no: " bigscape || exit 100
    ###################
    #### TAXONOMY #####
    ###################
    if [[ "${taxa_1}" == yes ]]; then
      # #Pick which assembly you like best by comparing assembly_stats.csv
      # Copy best assembly to 02_Assembly directory
      # For my analysis, 58 of the best assemblies were megahit and 2 were metaspades.
      mkdir -p ./02_Assembly/taxonomy
      echo Determining the overall taxonomic makeup of the ${sample} ${best_assembly} assembly
      sendsketch.sh in=${best_assembly} >> ./02_Assembly/taxonomy/${sample}_RefSeq.tsv 2>&1
      sendsketch.sh in=${best_assembly} nt >> ./02_Assembly/taxonomy/${sample}_nt.tsv 2>&1
    elif [[ "${taxa_1}" ==  no ]]; then :
    fi
    # Deactivate assembly conda environment
    conda deactivate
    ###################
    ##### ANVI'O ######
    ###################
    ###############
    ### CONDA #####
    ###############
    # Please follow the conda installation instructions at: https://anvio.org/install
    #read -r -p "Would you like to analyze the ${sample} assembly using anvi'o? yes or no: " anvio_1 || exit 100
    if [[ "${anvio_1}" == yes ]]; then
      #Or, try to determine taxonomy on a per-contig basis.  If this is not sensitive enough, try BLAST instead.
      # sendsketch.sh in=${best_assembly} persequence minhits=1 records=4 >> ./02_Assembly/taxonomy/${sample}.log 2>&1
      echo Starting anvi\'o on ${sample}

      conda activate anvio-7.1

      prefix=$(echo $sample | sed 's/-/_/g')
      echo running anvi\'o on ${sample}
      mkdir -p ./03_Anvio/${sample}/COG
      mkdir -p ./03_Anvio/${sample}/KOfam
      mkdir -p ./03_Anvio/${sample}/Pfam
      mkdir -p ./03_Anvio/${sample}/prodigal
      mkdir -p ./03_Anvio/logs
      #Simplify headers as to not cause future anger; make a key so that contigs can be matched up later if needed; length cut off of 500bp; the prefix argument doesn't take hyphens
      echo Reformating ${sample} fasta file and removing contigs \< 1000bp
      anvi-script-reformat-fasta --seq-type NT --simplify-names --prefix ${prefix} -r ./03_Anvio/${sample}/${sample}_rename_key.txt -l 1000 -o ./03_Anvio/${sample}/${sample}_renamed_c1k.fa ${best_assembly} >> ./03_Anvio/logs/${sample}.log 2>&1
      #Generate a contigs database for each assembly
      echo making contigs database for ${sample} assembly
      anvi-gen-contigs-database -T ${CPU} -f ./03_Anvio/${sample}/${sample}_renamed_c1k.fa -n ${sample} -o ./03_Anvio/${sample}/${sample}.db >> ./03_Anvio/logs/${sample}.log 2>&1
      echo running hmm
      anvi-run-hmms -c ./03_Anvio/${sample}/${sample}.db -T ${CPU} --also-scan-trnas >> ./03_Anvio/logs/${sample}.log 2>&1
      echo running ncbi COGS on ${sample}
      anvi-run-ncbi-cogs -c ./03_Anvio/${sample}/${sample}.db -T ${CPU} --sensitive >> ./03_Anvio/logs/${sample}.log 2>&1
      echo running KEGG on ${sample}
      anvi-run-kegg-kofams -c ./03_Anvio/${sample}/${sample}.db -T ${CPU} >> ./03_Anvio/logs/${sample}.log 2>&1
      echo running pfams on ${sample}
      anvi-run-pfams -c ./03_Anvio/${sample}/${sample}.db -T ${CPU} >> ./03_Anvio/logs/${sample}.log 2>&1
      echo exporting functions to ${sample} contigs database
      anvi-export-functions -c ./03_Anvio/${sample}/${sample}.db -o ./03_Anvio/${sample}/COG/${sample}_COG20Category.txt --annotation-sources COG20_CATEGORY >> ./03_Anvio/logs/${sample}.log 2>&1
      anvi-export-functions -c ./03_Anvio/${sample}/${sample}.db -o ./03_Anvio/${sample}/COG/${sample}_COG20Function.txt --annotation-sources COG20_FUNCTION >> ./03_Anvio/logs/${sample}.log 2>&1
      anvi-export-functions -c ./03_Anvio/${sample}/${sample}.db -o ./03_Anvio/${sample}/KOfam/${sample}_KOfam.txt --annotation-sources KOfam >> ./03_Anvio/logs/${sample}.log 2>&1
      anvi-export-functions -c ./03_Anvio/${sample}/${sample}.db -o ./03_Anvio/${sample}/Pfam/${sample}_Pfam.txt --annotation-sources Pfam >> ./03_Anvio/logs/${sample}.log 2>&1
      echo exporing gene calls to ${sample} contig database
      anvi-export-gene-calls -c ./03_Anvio/${sample}/${sample}.db --gene-caller prodigal -o ./03_Anvio/${sample}/prodigal/${sample}_AllGeneCalls.txt >> ./03_Anvio/logs/${sample}.log 2>&1
      echo running SCG taxonomy on ${sample}
      anvi-run-scg-taxonomy -c ./03_Anvio/${sample}/${sample}.db -T ${CPU} --all-hits-output-file ./03_Anvio/${sample}/${sample}_hits.txt >> ./03_Anvio/logs/${sample}.log 2>&1
      anvi-estimate-scg-taxonomy -c ./03_Anvio/${sample}/${sample}.db --metagenome-mode --per-scg-output-file ./03_Anvio/${sample}/${sample}-scg.txt

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

    elif [[ "${anvio_1}" == no ]]; then :
    fi
    ######################
    ### BGC Annotation ###
    ######################
    ################
    ### DOCKER #####
    ################
    # Both antismash and bigscape are run in Docker containers. Please visit their websites for installation instructions:
    #https://docs.antismash.secondarymetabolites.org/install/
    #https://nselem.github.io/bigscape-corason/installation/
    #read -r -p "Would you like to annotate your ${sample} metagenome for biosynthetic gene clusters? yes or no: " smash || exit 100
    if [[ "${smash}" == yes ]]; then
      #read -r -p "What minimum contig length do you want to analyze for biosynthetic gene clusters?    " min_BGC || exit 100

      mkdir -p ./04_Annotations/antismash/${sample}_c${min_BGC}
      echo Running antismash to annotate ${sample} secondary metabolites

      ~/bin/run_antismash \
      ./03_Anvio/${sample}/${sample}_renamed_c1k.fa \
      ./04_Annotations/antismash/${sample}_c${min_BGC} \
      --hmmdetection-strictness relaxed \
      --cb-general \
      --cb-knownclusters \
      --cb-subclusters \
      --asf --rre --pfam2go --smcog-trees \
      --genefinding-tool prodigal-m \
      --cc-mibig \
      --minlength ${min_BGC} \
      -c ${CPU} >> ./04_Annotations/antismash/${sample}.log 2>&1 || "No BGCs found in ${sample}" && exit 100

      cd ${project_dir}
    elif [[ "${smash}" == no ]]; then :
    fi
  done

  #read -r -p "Would you like to analyze annotated biosynthetic gene clusters for each sample against the MiBIG database?  yes or no: " bigscape || exit 100
  if [[ "${bigscape}" == yes ]]; then
    echo Running bigscape on ${project} antismash genbank files

    ~/bin/run_bigscape ./04_Annotations/antismash ./04_Annotations/bigscape --mix --mibig -c ${CPU} >> ./04_Annotations/bigscape/${project}.log 2&>1

    cd ${project_dir}

  elif [[ "${bigscape}" == no ]]; then :
  fi

elif [[ "${assemble}" == no ]]; then :
fi
#########################################################################################################################################
##############################                BEGIN CO-ASSEMBLY                                            ##############################
#########################################################################################################################################
read -r -p "Would you like to co-assemble reads from your entire project?  yes or no:  " coassemble || exit 100
if [[ "${coassemble}" == yes ]]; then
  echo starting to assemble ${project} reads into one co-assembly
  mkdir -p ./02_Assembly/metaspades/logs
  mkdir -p ./02_Assembly/megahit-large/logs
  mkdir -p ./02_Assembly/idba_ud/logs
  mkdir -p ./02_Assembly/stats
  mkdir -p ./02_Assembly/quast/${project}
  ## Do not make a directory for megahit, it will make it itself and will fail if the output directory already exists, as to not overwrite results. ##
  ## Since megahit will not write to external harddrives, the code is written to write data to home directory then move it to the working directory.
  ################
  ### Assembly ###
  ################
  conda activate assembly
  read -r -p "Which HQ reads dataset do you want to use to co-assemble $project? HQ HQ_ecc HQ_merged Skip:  " HQ_Reads || exit 100
  read -r -p "Would you like to evaluate the quality of all the assemblies? yes or no:  " quast || exit 100



  if [[ "${HQ_Reads}" == HQ ]]
  then
    # Check to see if the project reads have already been created. If not, then create them.
    if test -f "./01_HQ_Reads/${project}_HQ.fq.gz";
    then :
    else
      cat ./01_HQ_Reads/*_HQ.fq.gz > ./01_HQ_Reads/${project}_HQ.fq.gz
    fi
    # Check to see if metaspades was already run
    if test -f "./02_Assembly/metaspades/${project}/scaffolds.fasta";
    then :
    else
      echo Using metaspades to assemble all ${HQ_Reads} ${project} reads.
      mkdir -p ./02_Assembly/metaspades/${project}
      metaspades.py --only-assembler -k 21,33,55,77,99,127 --12 ./01_HQ_Reads/${project}_HQ.fq.gz -t ${CPU} -m ${RAM} -o ./02_Assembly/metaspades/${project} >> ./02_Assembly/metaspades/logs/${project}.log 2>&1
      cp ./02_Assembly/metaspades/${project}/scaffolds.fasta ./02_Assembly/quast/${project}/metaspades.fasta
    fi
    if test -f "./02_Assembly/quast/${project}/megahit.fasta";
    then :
    else
      #Assemble reads per site using megahit
      echo Using megahit to assemble ${HQ_Reads} ${project} reads.
      megahit --12 ./01_HQ_Reads/${project}_HQ.fq.gz --presets meta-large -m 0.5 -t ${CPU} -o ~/${project} >> ./02_Assembly/megahit-large/logs/${project}.log 2>&1
      cp  ~/${project}/final.contigs.fa ./02_Assembly/quast/${project}/megahit.fasta
      mkdir -p ./02_Assembly/megahit-large/${project}
      mv ~/${project} ./02_Assembly/megahit-large
    fi

    if test -f "./01_HQ_Reads/${project}_HQ.fa";
    then :
    else
      #IDBA can only take fasta files as input; therefore, fastq.gz files must be reformated to fasta files
      reformat.sh in=./01_HQ_Reads/${project}_HQ.fq.gz out=./01_HQ_Reads/${project}_HQ.fa ow=t >> ./02_Assembly/idba_ud/logs/${project}.log 2>&1
      #fq2fa --paired ./01_HQ_Reads/${project}_HQ.fq.gz ./01_HQ_Reads/${project}_HQ.fa
    fi

    if test -f "./02_Assembly/idba_ud/${project}/scaffold.fa";
    then :
    else
      # #Assemble reads per site using idba_ud
      echo Using idba to assemble ${HQ_Reads} ${project} reads.
      mkdir -p ./02_Assembly/idba_ud/${project}
      idba_ud -r ./01_HQ_Reads/${project}_HQ.fa --num_threads=${CPU} -o ./02_Assembly/idba_ud/${project} >> ./02_Assembly/idba_ud/logs/${project}.log 2>&1
      cp ./02_Assembly/idba_ud/${project}/scaffold.fa ./02_Assembly/quast/${project}/idba.fasta
    fi

    if grep -m 1 -q "${HQ_Reads}" ./02_Assembly/stats/${project}_assembly_summary.tsv;
    then echo ${HQ_Reads} assemblies have already had stats run.
    else
      # --- Evaluation ---
      echo -e "metagenome_id\treads\tassembler\tn_contigs\tcontig_bp\tgap_pct\tctg_L50\tctg_max\tgene_counts\tmapped_pct\tavg_cov\tassembly_performance" >> ./02_Assembly/stats/${project}_assembly_summary.tsv
      echo Evaluating ${project} assemblies with AssemblyStats
      ##statswrapper.sh ./02_Assembly/quast/${project}/*.fasta minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/stats/${project}.txt
      stats.sh ./02_Assembly/metaspades/${project}/scaffolds.fasta minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/metaspades/${project}/stats.txt
      spades_stats=$(cat ./02_Assembly/metaspades/${project}/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d')
      spades_L50=$(cat ./02_Assembly/metaspades/${project}/stats.txt | cut -f9 -d$'\t' | sed '2q;d')
      stats.sh ./02_Assembly/megahit-large/${project}/final.contigs.fa minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/megahit-large/${project}/stats.txt
      megahit_stats=$(cat ./02_Assembly/megahit-large/${project}/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d')
      megahit_L50=$(cat ./02_Assembly/megahit-large/${project}/stats.txt | cut -f9 -d$'\t' | sed '2q;d')
      stats.sh ./02_Assembly/idba_ud/${project}/scaffold.fa minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/idba_ud/${project}/stats.txt
      idba_stats=$(cat ./02_Assembly/idba_ud/${project}/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d')
      idba_L50=$(cat ./02_Assembly/idba_ud/${project}/stats.txt | cut -f9 -d$'\t' | sed '2q;d')

      #Evaluate assemblies with Quast (leave out "-R ref.fa if you don't have a reference)

      echo Predicting prokaryotic genes in each ${project} assembly
      #Count gene predictions for prokaryotes using prodigal for each assembly
      python ~/bin/LengthFilter.py ./02_Assembly/metaspades/${project}/scaffolds.fasta -m 1000 > ./02_Assembly/metaspades/${project}/${project}_c1k.fa
      prodigal -i ./02_Assembly/metaspades/${project}/${project}_c1k.fa -a ./02_Assembly/metaspades/${project}/${project}_c1k_proteins.faa -f gff -p meta -o ./02_Assembly/metaspades/${project}/${project}_c1k.gff >> ./02_Assembly/metaspades/logs/${project}.log 2>&1
      spades_counts=$(grep '^>' ./02_Assembly/metaspades/${project}/${project}_c1k_proteins.faa | wc -l )

      python ~/bin/LengthFilter.py ./02_Assembly/megahit-large/${project}/final.contigs.fa -m 1000 > ./02_Assembly/megahit-large/${project}/${project}_c1k.fa
      prodigal -i ./02_Assembly/megahit-large/${project}/${project}_c1k.fa -a ./02_Assembly/megahit-large/${project}/${project}_c1k_proteins.faa -f gff -p meta -o ./02_Assembly/megahit-large/${project}/${project}_c1k.gff >> ./02_Assembly/megahit-large/logs/${project}.log 2>&1
      megahit_counts=$(grep '^>' ./02_Assembly/megahit-large/${project}/${project}_c1k_proteins.faa | wc -l )

      python ~/bin/LengthFilter.py ./02_Assembly/idba_ud/${project}/scaffold.fa -m 1000 > ./02_Assembly/idba_ud/${project}/${project}_c1k.fa
      prodigal -i ./02_Assembly/idba_ud/${project}/${project}_c1k.fa -a ./02_Assembly/idba_ud/${project}/${project}_c1k_proteins.faa -f gff -p meta -o ./02_Assembly/idba_ud/${project}/${project}_c1k.gff >> ./02_Assembly/idba_ud/logs/${project}.log 2>&1
      idba_counts=$(grep '^>' ./02_Assembly/idba_ud/${project}/${project}_c1k_proteins.faa | wc -l )

      echo Determining metagenome content within ${project} assemblies
      # If project reads have not already been concatenated, do so now.
      if test -f "./01_HQ_Reads/${project}_HQ.fq.gz"; then :
      else
        echo Concatenating ${project} HQ reads
        cat ./01_HQ_Reads/*_HQ.fq.gz > ./01_HQ_Reads/${project}_HQ.fq.gz
      fi
      # Note that these are reads mapped to contigs > 2,000 bp.
      #(i.e. % mapped) by mapping reads by to assembly
      bbmap.sh in=./01_HQ_Reads/${project}_HQ.fq.gz ref=./02_Assembly/metaspades/${project}/${project}_c1k.fa nodisk covhist=./02_Assembly/metaspades/${project}/covhist.txt covstats=./02_Assembly/metaspades/${project}/covstats.txt out=./02_Assembly/metaspades/${project}/${project}_reads_assembled.fq.gz outu=./02_Assembly/metaspades/${project}/${project}_reads_unassembled.fq.gz fast=t ambig=best ow >> ./02_Assembly/metaspades/${project}/map.log 2>&1
      spades_match=$(grep Percent\ mapped: ./02_Assembly/metaspades/${project}/map.log | cut -f2 -d$'\t' | tail -1)
      spades_cov=$(grep Average\ coverage: ./02_Assembly/metaspades/${project}/map.log | cut -f2 -d$'\t' | tail -1)
      spades_Pa=$(echo "$spades_match*$spades_L50" | bc)

      bbmap.sh in=./01_HQ_Reads/${project}_HQ.fq.gz ref=./02_Assembly/megahit-large/${project}/${project}_c1k.fa nodisk covhist=./02_Assembly/megahit-large/${project}/covhist.txt covstats=./02_Assembly/megahit-large/${project}/covstats.txt outm=./02_Assembly/megahit-large/${project}/${project}_reads_assembled.fq.gz outu=./02_Assembly/megahit-large/${project}/${project}_reads_unassembled.fq.gz fast=t ambig=best ow >> ./02_Assembly/megahit-large/${project}/map.log 2>&1
      megahit_match=$(grep Percent\ mapped: ./02_Assembly/megahit-large/${project}/map.log | cut -f2 -d$'\t' | tail -1)
      megahit_cov=$(grep Average\ coverage: ./02_Assembly/megahit-large/${project}/map.log | cut -f2 -d$'\t' | tail -1)
      megahit_Pa=$(echo "$megahit_match*$megahit_L50" | bc)

      bbmap.sh in=./01_HQ_Reads/${project}_HQ.fq.gz ref=./02_Assembly/idba_ud/${project}/${project}_c1k.fa nodisk covhist=./02_Assembly/idba_ud/${project}/covhist.txt covstats=./02_Assembly/idba_ud/${project}/covstats.txt outm=./02_Assembly/idba_ud/${project}/${project}_reads_assembled.fq.gz outu=./02_Assembly/idba_ud/${project}/${project}_reads_unassembled.fq.gz fast=t ambig=best ow >> ./02_Assembly/idba_ud/${project}/map.log 2>&1
      idba_match=$(grep Percent\ mapped: ./02_Assembly/idba_ud/${project}/map.log | cut -f2 -d$'\t' | tail -1)
      idba_cov=$(grep Average\ coverage: ./02_Assembly/idba_ud/${project}/map.log | cut -f2 -d$'\t' | tail -1)
      idba_Pa=$(echo "$idba_match*$idba_L50" | bc)


      echo Creating assembly comparision for ${project}
      echo -e "${project}\t${HQ_Reads}\tmetaspades\t${spades_stats}\t${spades_counts}\t${spades_match}\t${spades_cov}\t${spades_Pa}" >> ./02_Assembly/stats/${project}_assembly_summary.tsv
      echo -e "${project}\t${HQ_Reads}\tmegahit\t${megahit_stats}\t${megahit_counts}\t${megahit_match}\t${megahit_cov}\t${megahit_Pa}" >> ./02_Assembly/stats/${project}_assembly_summary.tsv
      echo -e "${project}\t${HQ_Reads}\tidba_ud\t${idba_stats}\t${idba_counts}\t${idba_match}\t${idba_cov}\t${idba_Pa}" >> ./02_Assembly/stats/${project}_assembly_summary.tsv
    fi
  elif [[ "${HQ_Reads}" == HQ_ecc ]]
  then
    # First test if any output files were already created, if not then create them.
    if test -f "./01_HQ_Reads/error_corrected/${project}_HQ_ecc.fq.gz";
    then :
    else
      echo Using metaspades to assemble ${HQ_Reads} ${project} reads.
      cat ./01_HQ_Reads/error_corrected/*_HQ_ecc.fq.gz > ./01_HQ_Reads/error_corrected/${project}_HQ_ecc.fq.gz
    fi

    if test -f "./02_Assembly/metaspades/${project}_ecc/scaffolds.fasta";
    then :
    else
      echo Using metaspades to assemble ${HQ_Reads} ${project} reads.
      mkdir -p ./02_Assembly/metaspades/${project}_ecc
      metaspades.py -k 21,33,55,77,99,127 --12 ./01_HQ_Reads/error_corrected/${project}_HQ_ecc.fq.gz -t ${CPU} -m ${RAM} -o ./02_Assembly/metaspades/${project}_ecc >> ./02_Assembly/metaspades/logs/${project}.log 2>&1
      cp ./02_Assembly/metaspades/${project}_ecc/scaffolds.fasta ./02_Assembly/quast/${project}/metaspades_ecc.fasta
    fi

    if test -f "./02_Assembly/quast/${project}/megahit_ecc.fasta";
    then :
    else
    #Assemble reads per site using megahit
      echo Using megahit to assemble ${HQ_Reads} ${project} reads.
      # Becaue megahit won't write to an external hard drive, the working directory is being made in the home directory then copied over to the project working directory.
      megahit --12 ./01_HQ_Reads/error_corrected/${project}_HQ_ecc.fq.gz --presets meta-large -m 0.5 -t ${CPU} -o ~/${project}_ecc >> ./02_Assembly/megahit-large/logs/${project}.log 2>&1
      cp ~/${project}_ecc/final.contigs.fa ./02_Assembly/quast/${project}/megahit_ecc.fasta
      mv ~/${project}_ecc/ ./02_Assembly/megahit-large/
    fi

    if test -f "./01_HQ_Reads/error_corrected/${project}_HQ_ecc.fa";
    then :
    else
      #IDBA can only take fasta files as input; therefore, fastq.gz files must be reformated to fasta files
      reformat.sh in=./01_HQ_Reads/error_corrected/${project}_HQ_ecc.fq.gz out=./01_HQ_Reads/error_corrected/${project}_HQ_ecc.fa ow=t >> ./02_Assembly/idba_ud/logs/${project}.log 2>&1
      #fq2fa --paired ./01_HQ_Reads/error_corrected/${project}_HQ_ecc.fq.gz ./01_HQ_Reads/error_corrected/${project}_HQ_ecc.fa
    fi

    if test -f "./02_Assembly/idba_ud/${project}_ecc/scaffold.fa";
    then :
    else
      # #Assemble reads per site using idba_ud
      echo Using idba to assemble ${HQ_Reads} ${project} reads.

      mkdir -p ./02_Assembly/idba_ud/${project}_ecc
      idba_ud -r ./01_HQ_Reads/error_corrected/${project}_HQ_ecc.fa --num_threads=${CPU} --pre_correction -o ./02_Assembly/idba_ud/${project}_ecc >> ./02_Assembly/idba_ud/logs/${project}.log 2>&1
      cp ./02_Assembly/idba_ud/${project}_ecc/scaffold.fa ./02_Assembly/quast/${project}/idba_ecc.fasta
    fi
    # --- Evaluation ---
    if grep -m 1 -q "${HQ_Reads}" ./02_Assembly/stats/${project}_assembly_summary.tsv;
    then echo ${HQ_Reads} assemblies have already had stats run.
    else
      echo -e "metagenome_id\treads\tassembler\tn_contigs\tcontig_bp\tgap_pct\tctg_L50\tctg_max\tgene_counts\tmapped_pct\tavg_cov\tassembly_performance" >> ./02_Assembly/stats/${project}_assembly_summary.tsv
      echo Evaluating ${project} assemblies with AssemblyStats
      stats.sh ./02_Assembly/metaspades/${project}_ecc/scaffolds.fasta minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/metaspades/${project}_ecc/stats.txt
      spades_stats=$(cat ./02_Assembly/metaspades/${project}_ecc/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d')
      stats.sh ./02_Assembly/megahit-large/${project}_ecc/final.contigs.fa minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/megahit-large/${project}_ecc/stats.txt
      megahit_stats=$(cat ./02_Assembly/megahit-large/${project}_ecc/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d')
      stats.sh ./02_Assembly/idba_ud/${project}_ecc/scaffold.fa minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/idba_ud/${project}_ecc/stats.txt
      idba_stats=$(cat ./02_Assembly/idba_ud/${project}_ecc/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d')

      spades_L50=$(cat ./02_Assembly/metaspades/${project}_ecc/stats.txt | cut -f9 -d$'\t' | sed '2q;d')
      megahit_L50=$(cat ./02_Assembly/megahit-large/${project}_ecc/stats.txt | cut -f9 -d$'\t' | sed '2q;d')
      idba_L50=$(cat ./02_Assembly/idba_ud/${project}_ecc/stats.txt | cut -f9 -d$'\t' | sed '2q;d')

      #Evaluate assemblies with Quast (leave out "-R ref.fa if you don't have a reference)

      echo Predicting prokaryotic genes in each ${project} assembly
      #Count gene predictions for prokaryotes using prodigal for each assembly
      python ~/bin/LengthFilter.py ./02_Assembly/metaspades/${project}_ecc/scaffolds.fasta -m 1000 > ./02_Assembly/metaspades/${project}_ecc/${project}_c1k.fa
      prodigal -i ./02_Assembly/metaspades/${project}_ecc/${project}_c1k.fa -a ./02_Assembly/metaspades/${project}_ecc/${project}_c1k_proteins.faa -f gff -p meta -o ./02_Assembly/metaspades/${project}_ecc/${project}_c1k.gff >> ./02_Assembly/metaspades/logs/${project}.log 2>&1
      spades_counts=$(grep '^>' ./02_Assembly/metaspades/${project}_ecc/${project}_c1k_proteins.faa | wc -l )

      python ~/bin/LengthFilter.py ./02_Assembly/megahit-large/${project}_ecc/final.contigs.fa -m 1000 > ./02_Assembly/megahit-large/${project}_ecc/${project}_c1k.fa
      prodigal -i ./02_Assembly/megahit-large/${project}_ecc/${project}_c1k.fa -a ./02_Assembly/megahit-large/${project}_ecc/${project}_c1k_proteins.faa -f gff -p meta -o ./02_Assembly/megahit-large/${project}_ecc/${project}_c1k.gff >> ./02_Assembly/megahit-large/logs/${project}.log 2>&1
      megahit_counts=$(grep '^>' ./02_Assembly/megahit-large/${project}_ecc/${project}_c1k_proteins.faa | wc -l )

      python ~/bin/LengthFilter.py ./02_Assembly/idba_ud/${project}_ecc/scaffold.fa -m 1000 > ./02_Assembly/idba_ud/${project}_ecc/${project}_c1k.fa
      prodigal -i ./02_Assembly/idba_ud/${project}_ecc/${project}_c1k.fa -a ./02_Assembly/idba_ud/${project}_ecc/${project}_c1k_proteins.faa -f gff -p meta -o ./02_Assembly/idba_ud/${project}_ecc/${project}_c1k.gff >> ./02_Assembly/idba_ud/logs/${project}.log 2>&1
      idba_counts=$(grep '^>' ./02_Assembly/idba_ud/${project}_ecc/${project}_c1k_proteins.faa | wc -l )

      echo Determining metagenome content within ${project} assemblies
      # If project reads have not already been concatenated, do so now.
      if test -f "./01_HQ_Reads/${project}_HQ.fq.gz"; then :
      else
        echo Concatenating ${project} HQ reads
        cat ./01_HQ_Reads/*_HQ.fq.gz > ./01_HQ_Reads/${project}_HQ.fq.gz
      fi
      # Note that these are reads mapped to contigs > 2,000 bp.
      #(i.e. % mapped) by mapping reads by to assembly
      bbmap.sh in=./01_HQ_Reads/${project}_HQ.fq.gz ref=./02_Assembly/metaspades/${project}_ecc/${project}_c1k.fa nodisk covhist=./02_Assembly/metaspades/${project}_ecc/covhist.txt covstats=./02_Assembly/metaspades/${project}_ecc/covstats.txt out=./02_Assembly/metaspades/${project}_ecc/${project}_reads_assembled.fq.gz outu=./02_Assembly/metaspades/${project}_ecc/${project}_reads_unassembled.fq.gz fast=t ambig=best >> ./02_Assembly/metaspades/${project}_ecc/map.log 2>&1
      spades_match=$(grep Percent\ mapped: ./02_Assembly/metaspades/${project}_ecc/map.log | cut -f2 -d$'\t' | tail -1)
      spades_cov=$(grep Average\ coverage: ./02_Assembly/metaspades/${project}_ecc/map.log | cut -f2 -d$'\t' | tail -1)
      spades_Pa=$(echo "$spades_match*$spades_L50" | bc)

      bbmap.sh in=./01_HQ_Reads/${project}_HQ.fq.gz ref=./02_Assembly/megahit-large/${project}_ecc/${project}_c1k.fa nodisk covhist=./02_Assembly/megahit-large/${project}_ecc/covhist.txt covstats=./02_Assembly/megahit-large/${project}_ecc/covstats.txt outm=./02_Assembly/megahit-large/${project}_ecc/${project}_reads_assembled.fq.gz outu=./02_Assembly/megahit-large/${project}_ecc/${project}_reads_unassembled.fq.gz fast=t ambig=best >> ./02_Assembly/megahit-large/${project}_ecc/map.log 2>&1
      megahit_match=$(grep Percent\ mapped: ./02_Assembly/megahit-large/${project}_ecc/map.log | cut -f2 -d$'\t' | tail -1)
      megahit_cov=$(grep Average\ coverage: ./02_Assembly/megahit-large/${project}_ecc/map.log | cut -f2 -d$'\t' | tail -1)
      megahit_Pa=$(echo "$megahit_match*$megahit_L50" | bc)

      bbmap.sh in=./01_HQ_Reads/${project}_HQ.fq.gz ref=./02_Assembly/idba_ud/${project}_ecc/${project}_c1k.fa nodisk covhist=./02_Assembly/idba_ud/${project}_ecc/covhist.txt covstats=./02_Assembly/idba_ud/${project}_ecc/covstats.txt outm=./02_Assembly/idba_ud/${project}_ecc/${project}_reads_assembled.fq.gz outu=./02_Assembly/idba_ud/${project}_ecc/${project}_reads_unassembled.fq.gz fast=t ambig=best >> ./02_Assembly/idba_ud/${project}_ecc/map.log 2>&1
      idba_match=$(grep Percent\ mapped: ./02_Assembly/idba_ud/${project}_ecc/map.log | cut -f2 -d$'\t' | tail -1)
      idba_cov=$(grep Average\ coverage: ./02_Assembly/idba_ud/${project}_ecc/map.log | cut -f2 -d$'\t' | tail -1)
      idba_Pa=$(echo "$idba_match*$idba_L50" | bc)

      echo Creating assembly comparision for ${project}
      echo -e "${project}\t${HQ_Reads}\tmetaspades\t${spades_stats}\t${spades_counts}\t${spades_match}\t${spades_cov}\t${spades_Pa}" >> ./02_Assembly/stats/${project}_assembly_summary.tsv
      echo -e "${project}\t${HQ_Reads}\tmegahit\t${megahit_stats}\t${megahit_counts}\t${megahit_match}\t${megahit_cov}\t${megahit_Pa}" >> ./02_Assembly/stats/${project}_assembly_summary.tsv
      echo -e "${project}\t${HQ_Reads}\tidba_ud\t${idba_stats}\t${idba_counts}\t${idba_match}\t${idba_cov}\t${idba_Pa}" >> ./02_Assembly/stats/${project}_assembly_summary.tsv
    fi
  elif [[ "${HQ_Reads}" == HQ_merged ]]
  then

    if test -f "./01_HQ_Reads/merged_reads/${project}_HQ_merged.fq.gz";
    then :
    else
      cat ./01_HQ_Reads/merged_reads/*_HQ_merged.fq.gz > ./01_HQ_Reads/merged_reads/${project}_HQ_merged.fq.gz
      cat ./01_HQ_Reads/merged_reads/*_HQ_unmerged.fq.gz > ./01_HQ_Reads/merged_reads/${project}_HQ_unmerged.fq.gz
    fi

    if test -f "./02_Assembly/metaspades/${project}_merged/scaffolds.fasta";
    then :
    else
      echo Using metaspades to assemble ${HQ_Reads} ${project} reads.
      # Assemble using merged/unmerged reads
      mkdir -p ./02_Assembly/metaspades/${project}_merged
      metaspades.py  --only-assembler -k 21,33,55,77,99,127 --merged ./01_HQ_Reads/merged_reads/${project}_HQ_merged.fq.gz --12 ./01_HQ_Reads/merged_reads/${project}_HQ_unmerged.fq.gz -t ${CPU} -m ${RAM} -o ./02_Assembly/metaspades/${project}_merged >> ./02_Assembly/metaspades/logs/${project}.log 2>&1
      cp ./02_Assembly/metaspades/${project}_merged/scaffolds.fasta ./02_Assembly/quast/${project}/metaspades_merged.fasta
    fi

    if test -f "./02_Assembly/quast/${project}/megahit_merged.fasta";
    then :
    else
      echo Using megahit to assemble ${HQ_Reads} ${project} reads.
      megahit -r ./01_HQ_Reads/merged_reads/${project}_HQ_merged.fq.gz --12 ./01_HQ_Reads/merged_reads/${project}_HQ_unmerged.fq.gz --presets meta-large -m 0.5 -t ${CPU} -o ~/${project}_merged >> ./02_Assembly/megahit-large/logs/${project}.log 2>&1
      cp ~/${project}_merged/final.contigs.fa ./02_Assembly/quast/${project}/megahit_merged.fasta
      mkdir -p ./02_Assembly/megahit-large/${project}_merged/
      mv ~/${project}_merged/ ./02_Assembly/megahit-large/
      echo idba_ud does not assemble merged reads and will not be used for comparison.
    fi

    if grep -m 1 -q "${HQ_Reads}" ./02_Assembly/stats/${project}_assembly_summary.tsv;
    then echo ${HQ_Reads} assemblies have already had stats run.
    else
      echo -e "metagenome_id\treads\tassembler\tn_contigs\tcontig_bp\tgap_pct\tctg_L50\tctg_max\tgene_counts\tmapped_pct\tavg_cov\tassembly_performance" >> ./02_Assembly/stats/${project}_assembly_summary.tsv
      echo Evaluating ${project} assemblies with AssemblyStats
      stats.sh ./02_Assembly/metaspades/${project}_merged/scaffolds.fasta minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/metaspades/${project}_merged/stats.txt
      spades_stats=$(cat ./02_Assembly/metaspades/${project}_merged/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d')
      stats.sh ./02_Assembly/megahit-large/${project}_merged/final.contigs.fa minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/megahit-large/${project}_merged/stats.txt
      megahit_stats=$(cat ./02_Assembly/megahit-large/${project}_merged/stats.txt | cut -f2,4,5,9,15 -d$'\t' | sed '2q;d')

      spades_L50=$(cat ./02_Assembly/metaspades/${project}_merged/stats.txt | cut -f9 -d$'\t' | sed '2q;d')
      megahit_L50=$(cat ./02_Assembly/megahit-large/${project}_merged/stats.txt | cut -f9 -d$'\t' | sed '2q;d')

      #Evaluate assemblies with Quast (leave out "-R ref.fa if you don't have a reference)

      echo Predicting prokaryotic genes in each ${project} assembly
      #Count gene predictions for prokaryotes using prodigal for each assembly
      python ~/bin/LengthFilter.py ./02_Assembly/metaspades/${project}_merged/scaffolds.fasta -m 1000 > ./02_Assembly/metaspades/${project}_merged/${project}_c1k.fa
      prodigal -i ./02_Assembly/metaspades/${project}_merged/${project}_c1k.fa -a ./02_Assembly/metaspades/${project}_merged/${project}_c1k_proteins.faa -f gff -p meta -o ./02_Assembly/metaspades/${project}_merged/${project}_c1k.gff >> ./02_Assembly/metaspades/logs/${project}.log 2>&1
      spades_counts=$(grep '^>' ./02_Assembly/metaspades/${project}_merged/${project}_c1k_proteins.faa | wc -l )

      python ~/bin/LengthFilter.py ./02_Assembly/megahit-large/${project}_merged/final.contigs.fa -m 1000 > ./02_Assembly/megahit-large/${project}_merged/${project}_c1k.fa
      prodigal -i ./02_Assembly/megahit-large/${project}_merged/${project}_c1k.fa -a ./02_Assembly/megahit-large/${project}_merged/${project}_c1k_proteins.faa -f gff -p meta -o ./02_Assembly/megahit-large/${project}_merged/${project}_c1k.gff >> ./02_Assembly/megahit-large/logs/${project}.log 2>&1
      megahit_counts=$(grep '^>' ./02_Assembly/megahit-large/${project}_merged/${project}_c1k_proteins.faa | wc -l )

      echo Determining metagenome content within ${project} assemblies
      # If project reads have not already been concatenated, do so now.
      if test -f "./01_HQ_Reads/${project}_HQ.fq.gz"; then :
      else
        echo Concatenating ${project} HQ reads
        cat ./01_HQ_Reads/*_HQ.fq.gz > ./01_HQ_Reads/${project}_HQ.fq.gz
      fi
      # Note that these are reads mapped to contigs > 2,000 bp.
      #(i.e. % mapped) by mapping reads by to assembly

      bbmap.sh in=./01_HQ_Reads/${project}_HQ.fq.gz ref=./02_Assembly/metaspades/${project}_merged/${project}_c1k.fa nodisk covhist=./02_Assembly/metaspades/${project}_merged/covhist.txt covstats=./02_Assembly/metaspades/${project}_merged/covstats.txt out=./02_Assembly/metaspades/${project}_merged/${project}_reads_assembled.fq.gz outu=./02_Assembly/metaspades/${project}_merged/${project}_reads_unassembled.fq.gz fast=t ambig=best >> ./02_Assembly/metaspades/${project}_merged/map.log 2>&1
      spades_match=$(grep Percent\ mapped: ./02_Assembly/metaspades/${project}_merged/map.log | cut -f2 -d$'\t' | tail -1)
      spades_cov=$(grep Average\ coverage: ./02_Assembly/metaspades/${project}_merged/map.log | cut -f2 -d$'\t' | tail -1)
      spades_Pa=$(echo "$spades_match*$spades_L50" | bc)

      bbmap.sh in=./01_HQ_Reads/${project}_HQ.fq.gz ref=./02_Assembly/megahit-large/${project}_merged/${project}_c1k.fa nodisk covhist=./02_Assembly/megahit-large/${project}_merged/covhist.txt covstats=./02_Assembly/megahit-large/${project}_merged/covstats.txt outm=./02_Assembly/megahit-large/${project}_merged/${project}_reads_assembled.fq.gz outu=./02_Assembly/megahit-large/${project}_merged/${project}_reads_unassembled.fq.gz fast=t ambig=best >> ./02_Assembly/megahit-large/${project}_merged/map.log 2>&1
      megahit_match=$(grep Percent\ mapped: ./02_Assembly/megahit-large/${project}_merged/map.log | cut -f2 -d$'\t' | tail -1)
      megahit_cov=$(grep Average\ coverage: ./02_Assembly/megahit-large/${project}_merged/map.log | cut -f2 -d$'\t' | tail -1)
      metahit_Pa=$(echo "$megahit_match*$megahit_L50" | bc)

      echo Creating assembly comparision for ${project}
      echo -e "${project}\t${HQ_Reads}\tmetaspades\t${spades_stats}\t${spades_counts}\t${spades_match}\t${spades_cov}\t${spades_Pa}" >> ./02_Assembly/stats/${project}_assembly_summary.tsv
      echo -e "${project}\t${HQ_Reads}\tmegahit\t${megahit_stats}\t${megahit_counts}\t${megahit_match}\t${megahit_cov}\t${megahit_Pa}" >> ./02_Assembly/stats/${project}_assembly_summary.tsv
    fi
  elif [[ "${HQ_Reads}" == Skip ]]
  then :
  fi
  ######################
  #### ASSEMBLY QC #####
  ######################
  #read -r -p "Would you like to evaluate the quality of all the assemblies? yes or no:  " quast || exit 100
  if [[ "${quast}" == yes ]]; then
    echo Evaluating ${project} assemblies

    #statswrapper.sh ./02_Assembly/quast/${project}/*.fasta minscaf=1000 format=3 -Xmx${RAM}G ow=t out=./02_Assembly/stats/${project}.txt
    ~/GitHub/quast/metaquast.py -f -m 1000 -o ./02_Assembly/quast/${project} ./02_Assembly/quast/${project}/*.fasta >> ./02_Assembly/quast/${project}.log 2>&1
  elif [[ "${quast}" == no ]]; then :
  fi
  cut -f2,3,5,12 -d$'\t' ./02_Assembly/stats/${project}_assembly_summary.tsv
  ls `pwd`/02_Assembly/quast/${project}/*.fasta
  # Pick which assembly you like best by comparing assembly_stats.csv

  read -r -p "Analyze the ./02_Assembly/stats/${project}_assembly_summary.tsv and ./02_Assembly/quast/${project}/report.html and determine which ${project} assembly is best. Please enter the absoute path to the assembly.fasta file:  " best_assembly || exit 100
  read -r -p "Would you like to estimate the taxonomy of your metagenome? yes or no : " taxa || exit 100
  read -r -p "Would you like to analyze your metagenome using anvi'o? yes or no: " anvio || exit 100
  read -r -p "What minimum contig length do you want to analyze for biosynthetic gene clusters?    " min_BGC || exit 100
  read -r -p "Would you like to analyze annotated biosynthetic gene clusters for ${project} against the MiBIG database?  yes or no: " bigscape2 || exit 100

  ###################
  #### TAXONOMY #####
  ###################
  #read -r -p "Would you like to estimate the taxonomy of your metagenome? yes or no : " taxa || exit 100
  if [[ "${taxa}" == yes ]]; then
    # For my analysis, 58 of the best assemblies were megahit and 2 were metaspades.
    mkdir -p ./02_Assembly/taxonomy
    echo Determining the overall taxonomic makeup of the ${project} ${best_assembly} assembly
    sendsketch.sh in=${best_assembly} >> ./02_Assembly/taxonomy/${project}_RefSeq.tsv 2>&1
    sendsketch.sh in=${best_assembly} nt >> ./02_Assembly/taxonomy/${project}_nt.tsv 2>&1
  elif [[ "${taxa}" == no ]]; then :
  fi
  # Deactivate the assembly conda environment
  conda deactivate
  ###################
  ##### ANVI'O ######
  ###################
  #read -r -p "Would you like to analyze your metagenome using anvi'o? yes or no: " anvio || exit 100
  if [[ "${anvio}" == yes ]]; then
    #Or, try to determine taxonomy on a per-contig basis.  If this is not sensitive enough, try BLAST instead.
    # sendsketch.sh in=${best_assembly} persequence minhits=1 records=4 >> ./02_Assembly/taxonomy/${project}.log 2>&1
    echo Starting anvi\'o on ${project}

    conda activate anvio-7.1
    # anvi'o hates dashes.
    prefix=$(echo $project | sed 's/-/_/g')
    echo running anvi\'o on ${project}
    mkdir -p ./03_Anvio/${project}/COG
    mkdir -p ./03_Anvio/${project}/KOfam
    mkdir -p ./03_Anvio/${project}/Pfam
    mkdir -p ./03_Anvio/${project}/prodigal
    mkdir -p ./03_Anvio/

    if test -f "./03_Anvio/${project}/${project}_renamed_c1k.fa"; then :
    else
      #Simplify headers as to not cause future anger; make a key so that contigs can be matched up later if needed; length cut off of 500bp; the prefix argument doesn't take hyphens
      echo Reformating ${project} fasta file and removing contigs \< 1000bp
      anvi-script-reformat-fasta --seq-type NT --simplify-names --prefix ${prefix} -r ./03_Anvio/${project}/${project}_rename_key.txt -l 1000 -o ./03_Anvio/${project}/${project}_renamed_c1k.fa ${best_assembly} >> ./03_Anvio/logs/${project}.log 2>&1
      #Generate a contigs database for each assembly
      echo making contigs database for ${project} assembly
      anvi-gen-contigs-database -T ${CPU} -f ./03_Anvio/${project}/${project}_renamed_c1k.fa -n ${project} -o ./03_Anvio/${project}/${project}.db >> ./03_Anvio/logs/${project}.log 2>&1
      echo running hmm
      anvi-run-hmms -c ./03_Anvio/${project}/${project}.db -T ${CPU} --also-scan-trnas >> ./03_Anvio/logs/${project}.log 2>&1
      echo running ncbi COGS
      anvi-run-ncbi-cogs -c ./03_Anvio/${project}/${project}.db -T ${CPU} --sensitive >> ./03_Anvio/logs/${project}.log 2>&1
      echo running KEGG
      anvi-run-kegg-kofams -c ./03_Anvio/${project}/${project}.db -T ${CPU} >> ./03_Anvio/logs/${project}.log 2>&1
      echo running pfams
      anvi-run-pfams -c ./03_Anvio/${project}/${project}.db -T ${CPU} >> ./03_Anvio/logs/${project}.log 2>&1
      echo exporting functions to contigs database
      anvi-export-functions -c ./03_Anvio/${project}/${project}.db -o ./03_Anvio/${project}/COG/${project}_COG20Category.txt --annotation-sources COG20_CATEGORY >> ./03_Anvio/logs/${project}.log 2>&1
      anvi-export-functions -c ./03_Anvio/${project}/${project}.db -o ./03_Anvio/${project}/COG/${project}_COG20Function.txt --annotation-sources COG20_FUNCTION >> ./03_Anvio/logs/${project}.log 2>&1
      anvi-export-functions -c ./03_Anvio/${project}/${project}.db -o ./03_Anvio/${project}/KOfam/${project}_KOfam.txt --annotation-sources KOfam >> ./03_Anvio/logs/${project}.log 2>&1
      anvi-export-functions -c ./03_Anvio/${project}/${project}.db -o ./03_Anvio/${project}/Pfam/${project}_Pfam.txt --annotation-sources Pfam >> ./03_Anvio/logs/${project}.log 2>&1
      echo exporing gene calls to contig database
      anvi-export-gene-calls -c ./03_Anvio/${project}/${project}.db --gene-caller prodigal -o ./03_Anvio/${project}/prodigal/${project}_AllGeneCalls.txt >> ./03_Anvio/logs/${project}.log 2>&1
      echo running SCG taxonomy
      anvi-run-scg-taxonomy -c ./03_Anvio/${project}/${project}.db -T ${CPU} --all-hits-output-file ./03_Anvio/${project}/${project}_hits.txt >> ./03_Anvio/logs/${project}.log 2>&1
      anvi-estimate-scg-taxonomy -c ./03_Anvio/${project}/${project}.db --metagenome-mode --per-scg-output-file ./03_Anvio/${project}/${project}-scg.txt

      #Remove zeros from column 5 e-value and only write significant hits
      echo creating functional profile on ${project}
      awk -F "\t" '{ if(($5 != 0) && ($5 <= .0001)) { print } }' ./03_Anvio/${project}/KOfam/${project}_KOfam.txt >> ./03_Anvio/${project}/KOfam/${project}_KOfamEvalue_filtered.txt
      awk -F "\t" '{ if(($5 != 0) && ($5 <= .0001)) { print } }' ./03_Anvio/${project}/Pfam/${project}_Pfam.txt >> ./03_Anvio/${project}/Pfam/${project}_PfamEvalue_filtered.txt
      awk -F "\t" '{ if(($5 != 0) && ($5 <= .0001)) { print } }' ./03_Anvio/${project}/COG/${project}_COG20Function.txt >> ./03_Anvio/${project}/COG/${project}_COG20FunctionEvalue_filtered.txt
      awk -F "\t" '{ if(($5 != 0) && ($5 <= .0001)) { print } }' ./03_Anvio/${project}/COG/${project}_COG20Category.txt >> ./03_Anvio/${project}/COG/${project}_COG20CategoryEvalue_filtered.txt
      #Cut unique values from column 4 function, sort them and then count them
      cut -f4 ./03_Anvio/${project}/KOfam/${project}_KOfamEvalue_filtered.txt | sort | uniq -c >> ./03_Anvio/${project}/KOfam/${project}_KOfamCounts.txt
      cut -f4 ./03_Anvio/${project}/Pfam/${project}_PfamEvalue_filtered.txt | sort | uniq -c >> ./03_Anvio/${project}/Pfam/${project}_PfamCounts.txt
      cut -f4 ./03_Anvio/${project}/COG/${project}_COG20FunctionEvalue_filtered.txt | sort | uniq -c >> ./03_Anvio/${project}/COG/${project}_COG20Function.txt
      cut -f4 ./03_Anvio/${project}/COG/${project}_COG20CategoryEvalue_filtered.txt | sort | uniq -c >> ./03_Anvio/${project}/COG/${project}_COG20Category.txt
    fi
    # Anvi'o is only effienct when analyzing < 20,000 splits. So choose a bp length cutoff so that only <20,000 splits are included
    read -r -p "Run :  sqlite3 ${project_dir}/03_Anvio/${project}/${project}.db 'select count(*) from contigs_basic_info where length > ####;   and enter the length cutoff that results in <20,000 splits. Numbers only:  " len_cutoff || exit 100

    for sample in `cat samples.txt`; do
      if test -f "./03_Anvio/${project}/profiles_c${len_cutoff}/${sample}/PROFILE.db"; then :
      else
        echo creating bam file for ${sample} HQ reads against the ${project} best assembly.
        prefix2=$(echo $sample | sed 's/-/_/g')
        mkdir -p ./04_Mapping/${project}
        bbmap.sh in=./01_HQ_Reads/${sample}_HQ.fq.gz ref=./03_Anvio/${project}/${project}_renamed_c1k.fa out=./04_Mapping/${project}/${sample}.sam bamscript=bs.sh outm=./04_Mapping/${project}/${project}_mapped.fq.gz outu=./04_Mapping/${project}/${project}_unmapped.fq.gz ihist=./04_Mapping/${project}/${project}_insertsize.txt interleaved=t pigz=t threads=${CPU} -Xmx${RAM}g ow >> ./04_Mapping/${project}/${project}.log 2>&1
        sh bs.sh >> ./04_Mapping/${project}/${project}.log 2>&1
        rm bs.sh
        echo creating anvid\'o profile for ${sample} sample in project ${project}
        anvi-profile -i ./04_Mapping/${project}/${sample}_sorted.bam -c ./03_Anvio/${project}/${project}.db -T ${CPU} -M ${len_cutoff} --sample-name ${prefix2} --output-dir ./03_Anvio/${project}/profiles_c${len_cutoff}/${sample} >> ./03_Anvio/logs/${project}.log 2>&1
      fi

    done

    if test -f "./03_Anvio/${project}/profiles_c${len_cutoff}/merged/PROFILE.db"; then :
    else
      echo merging ${project} profiles
      anvi-merge --enforce-hierarchical-clustering ./03_Anvio/${project}/profiles_c${len_cutoff}/*/PROFILE.db -o ./03_Anvio/${project}/profiles_c${len_cutoff}/merged -c ./03_Anvio/${project}/${project}.db >> ./03_Anvio/logs/${project}.log 2>&1
    fi

    read -r -p "Which binner would you like to use to create Metagenomic Assembled Genomes? concoct metabat2 maxbin2 dastool binsanity :  " binner || exit 100

    echo Creating MAGs using ${binner}
    anvi-cluster-contigs -c ./03_Anvio/${project}/${project}.db -p ./03_Anvio/${project}/profiles_c${len_cutoff}/merged/PROFILE.db -C ${binner} --driver ${binner} -T ${CPU} --log-file ./03_Anvio/${project}/${project}.log --just-do-it >> ./03_Anvio/logs/${project}.log 2>&1
    echo Create an anvi\'o summary
    anvi-summarize -p ./03_Anvio/${project}/profiles_c${len_cutoff}/merged/PROFILE.db -c ./03_Anvio/${project}/${project}.db -C ${binner} -o ./03_Anvio/${project}/SUMMARY_${binner}

    conda deactivate
  elif [[ "${anvio}" == no ]]; then :
  fi
  ######################
  ### BGC Annotation ###
  ######################
  if test -f "./04_Annotations/antismash/${project}_c${min_BGC}/${project}_renamed_c1k/index.html"; then :
  else
    #read -r -p "What minimum contig length do you want to analyze for biosynthetic gene clusters?    " min_BGC || exit 100

    mkdir -p ./04_Annotations/antismash/${project}_c${min_BGC}
    echo Running antismash to annotate ${project} secondary metabolites

    ~/bin/run_antismash \
    ./03_Anvio/${project}/${project}_renamed_c1k.fa \
    ./04_Annotations/antismash/${project}_c${min_BGC} \
    --hmmdetection-strictness relaxed \
    --cb-general \
    --cb-knownclusters \
    --cb-subclusters \
    --asf --rre --pfam2go --smcog-trees \
    --genefinding-tool prodigal-m \
    --cc-mibig \
    --minlength ${min_BGC} \
    -c ${CPU} >> ./${project}.log 2>&1 || "No BGCs found in ${project}" && exit 100

    cd ${project_dir}
  fi

  #read -r -p "Would you like to analyze annotated biosynthetic gene clusters for ${project} against the MiBIG database?  yes or no: " bigscape2 || exit 100
  if [[ "${bigscape2}" == yes ]]; then
    echo Running bigscape on ${project} antismash genbank files
    mkdir -p ./04_Annotations/bigscape
    ~/bin/run_bigscape ./04_Annotations/antismash ./04_Annotations/bigscape --mix --mibig -c ${CPU} >> ./${project}.log 2&>1
  elif [[ "${bigscape2}" == no ]]; then :
  fi

  cd ${project_dir}
elif [[ "${coassemble}" == no ]]; then
  echo Co-Assemblies have not been generated nor analyzed for ${project}
fi

echo finished meow
