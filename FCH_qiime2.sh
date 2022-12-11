#!/bin/bash

# https://github.com/biovcnet/topic-amplicons/blob/master/Lesson03a/analysis.md
# https://rachaellappan.github.io/VL-QIIME2-analysis/index.html

#Navigate to the project directory
#Which should have two subdirectories : 16S and 18S , each having the subdirectory 00_Raw_Reads


#Before running this script, be sure your raw sequencing reads are in the 00_Raw_Reads folder
#Read file names must have the format: sampleid_00_L001_R#_001.fastq.gz
#The first column of metadata file should have header 'sampleid'

cd ./16S
pwd
#Make directory structure
mkdir -p work/alpha_rarefaction
mkdir -p work/beta_rarefaction
echo Starting qiime2
#Verify RAW amplicon sequences are in 00_Raw_Reads directory with no other files.


echo Importing fastq files into a qiime2 qza file format

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path 00_Raw_Reads --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path work/demux-paired-end.qza

echo Trimming primers from amplicons

# https://docs.qiime2.org/2022.8/plugins/available/cutadapt/trim-paired/
# Please see those docs at https://cutadapt.readthedocs.io for complete details.
# Cutadapt finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads.
# Inputs:
#   --i-demultiplexed-sequences ARTIFACT
#     SampleData[PairedEndSequencesWithQuality]
#                           The paired-end sequences to be trimmed.   [required]
# Parameters:
#   --p-cores INTEGER       Number of CPU cores to use.
#     Range(1, None)                                                [default: 1]
#   --p-adapter-f TEXT...   Sequence of an adapter ligated to the 3' end. The
#     List[Str]             adapter and any subsequent bases are trimmed. If a
#                           `$` is appended, the adapter is only found if it is
#                           at the end of the read. Search in forward read. If
#                           your sequence of interest is "framed" by a 5' and a
#                           3' adapter, use this parameter to define a "linked"
#                           primer - see https://cutadapt.readthedocs.io for
#                           complete details.                         [optional]
#   --p-front-f TEXT...     Sequence of an adapter ligated to the 5' end. The
#     List[Str]             adapter and any preceding bases are trimmed. Partial
#                           matches at the 5' end are allowed. If a `^`
#                           character is prepended, the adapter is only found if
#                           it is at the beginning of the read. Search in
#                           forward read.                             [optional]
#   --p-anywhere-f TEXT...  Sequence of an adapter that may be ligated to the
#     List[Str]             5' or 3' end. Both types of matches as described
#                           under `adapter` and `front` are allowed. If the
#                           first base of the read is part of the match, the
#                           behavior is as with `front`, otherwise as with
#                           `adapter`. This option is mostly for rescuing failed
#                           library preparations - do not use if you know which
#                           end your adapter was ligated to. Search in forward
#                           read.                                     [optional]
#   --p-adapter-r TEXT...   Sequence of an adapter ligated to the 3' end. The
#     List[Str]             adapter and any subsequent bases are trimmed. If a
#                           `$` is appended, the adapter is only found if it is
#                           at the end of the read. Search in reverse read. If
#                           your sequence of interest is "framed" by a 5' and a
#                           3' adapter, use this parameter to define a "linked"
#                           primer - see https://cutadapt.readthedocs.io for
#                           complete details.                         [optional]
#   --p-front-r TEXT...     Sequence of an adapter ligated to the 5' end. The
#     List[Str]             adapter and any preceding bases are trimmed. Partial
#                           matches at the 5' end are allowed. If a `^`
#                           character is prepended, the adapter is only found if
#                           it is at the beginning of the read. Search in
#                           reverse read.                             [optional]


# Trim amplicon 16S V4 515F–806R primers amplicon size 300-350bp
# The primer sequences without linker, pad, barcode, or adapter are as follows:
#
# Updated sequences: 515F (Parada)–806R (Apprill), forward-barcoded:
# FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT
# Original sequences: 515F (Caporaso)–806R (Caporaso), reverse-barcoded:
# FWD:GTGCCAGCMGCCGCGGTAA; REV:GGACTACHVGGGTWTCTAAT

# 515F forward primer, barcoded
# Field descriptions (space-delimited):
#
# 5′ Illumina adapter
# Golay barcode
# Forward primer pad
# Forward primer linker
# Forward primer (515F)
# AATGATACGGCGACCACCGAGATCTACACGCT XXXXXXXXXXXX TATGGTAATT GT GTGYCAGCMGCCGCGGTAA
#
# 806R reverse primer
# Field descriptions (space-delimited):
#
# Reverse complement of 3′ Illumina adapter
# Reverse primer pad
# Reverse primer linker
# Reverse primer (806R)
# CAAGCAGAAGACGGCATACGAGAT AGTCAGCCAG CC GGACTACNVGGGTWTCTAAT

qiime cutadapt trim-paired --i-demultiplexed-sequences work/demux-paired-end.qza --p-cores 12 --p-front-f GTGYCAGCMGCCGCGGTAA --p-front-r GGACTACNVGGGTWTCTAAT --o-trimmed-sequences work/primer-trimmed-demux-PE.qza --verbose &> 16S_primer_trimming.log

qiime demux summarize --i-data work/primer-trimmed-demux-PE.qza --o-visualization work/primer-trimmed-demux-PE.qzv

echo Visualize work/primer-trimmed-demux-PE.qzv in https://view.qiime2.org/ look for the basepair position where the median quality score is consistently less than 30.
read -r -p "At what bp position does the quality score on the FORWARD read consistently drop below 30? Enter numbers only:   " fwd_trim || exit 100
read -r -p "At what bp position does the quality score on the REVERSE read consistently drop below 30? Enter numbers only:   " rvs_trim || exit 100


echo Truncating reads for quality using DADA2
#Run DADA2 denoising and filtering
#Remove low quality bases from 3 prime side by truncating the reads
qiime dada2 denoise-paired --i-demultiplexed-seqs work/primer-trimmed-demux-PE.qza --p-trunc-len-f ${fwd_trim} --p-trunc-len-r ${rvs_trim} --output-dir work/DADA2_denoising_output --verbose &> 16S_DADA2_denoising.log

echo Visualizing DADA2 outputs

# Denoising stats
qiime metadata tabulate --m-input-file work/DADA2_denoising_output/denoising_stats.qza --o-visualization work/DADA2_denoising_output/denoising_stats.qzv

Representative sequences
qiime feature-table tabulate-seqs --i-data work/DADA2_denoising_output/representative_sequences.qza --o-visualization work/DADA2_denoising_output/rep_seqs.qzv

# Feature table
qiime feature-table summarize --i-table work/DADA2_denoising_output/table.qza --o-visualization work/DADA2_denoising_output/table.qzv

echo Classify the representative sequences

# Classify the 16S representative sequences
qiime feature-classifier classify-sklearn --i-classifier /mnt/d/HBOI/FCCHH/qiime2/16S_qiime2_wd/silva-138-99-515-806-nb-classifier.qza --i-reads work/DADA2_denoising_output/representative_sequences.qza --output-dir work/classified_sequences --verbose   &> classify_16S_rep_seqs.log

# Tabulate the features, their taxonomy and the confidence of taxonomy assignment
qiime metadata tabulate --m-input-file work/classified_sequences/classification.qza --o-visualization work/classified_sequences/taxonomy.qzv

# The pipeline requires only the rep seqs, not the classified rep seqs
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences work/DADA2_denoising_output/representative_sequences.qza --output-dir work/phylogeny --p-n-threads 16 --verbose   &> phylogenetic_tree_generation.log

echo making initial barplot
qiime taxa barplot --i-table work/DADA2_denoising_output/table.qza --i-taxonomy work/classified_sequences/classification.qza --m-metadata-file ../FCCHHc2.5.map.txt --o-visualization work/phylogeny/barplots_initial.qzv

echo remove low abundant features less than 10
qiime feature-table filter-samples --i-table work/DADA2_denoising_output/table.qza --p-min-frequency 10 --o-filtered-table work/phylogeny/frequency-filtered-table.qza

echo Removing chloroplast, mitochondria, and eukaryotes from 16S data
qiime taxa filter-table --i-table work/phylogeny/frequency-filtered-table.qza --i-taxonomy work/classified_sequences/classification.qza --p-exclude chloroplast,mitochondria,eukaryota --o-filtered-table work/phylogeny/final_table_no_mce.qza

qiime taxa filter-seqs --i-sequences work/DADA2_denoising_output/representative_sequences.qza --i-taxonomy work/classified_sequences/classification.qza --p-exclude chloroplast,mitochondria,eukaryota --o-filtered-sequences work/phylogeny/rep_seqs_no_mce.qza

echo Generating final barplots
qiime taxa barplot --i-table work/phylogeny/final_table_no_mce.qza --i-taxonomy work/classified_sequences/classification.qza --m-metadata-file ../FCCHHc2.5.map.txt --o-visualization work/phylogeny/barplots_no_mce_final.qzv

  # Also check remaining read counts
qiime feature-table summarize --i-table work/phylogeny/final_table_no_mce.qza --o-visualization work/phylogeny/feature_table_no_mce_final.qzv --m-sample-metadata-file ../FCCHHc2.5.map.txt

echo qiime2 is finished with 16S taxonomy
echo starting statistics

# Alpha and beta diversity
# The commands I describe show how these steps were generally carried out for the 16S rRNA and 18S rRNA datasets.

# Alpha rarefaction was done as follows:
echo The maximum depth parameter on line 160 can be changed by visualizing the work/DADA2_denoising_output/table.qva
echo and move the slide bar to determine ideal depth for retaining most features and samples.
# 16S rRNA dataset
qiime diversity alpha-rarefaction --i-table work/phylogeny/final_table_no_mce.qza --i-phylogeny work/phylogeny/rooted_tree.qza --p-max-depth 16386 --m-metadata-file ../FCCHHc2.5.map.txt --o-visualization work/alpha_rarefaction/rarefaction_3138.qzv

# The core metrics pipeline produces results of standard alpha and beta diversity metrics:

qiime diversity core-metrics-phylogenetic --i-phylogeny work/phylogeny/rooted_tree.qza --i-table work/phylogeny/final_table_no_mce.qza --m-metadata-file ../FCCHHc2.5.map.txt --p-sampling-depth 16386 --output-dir work/core_metrics --p-n-jobs-or-threads 12 --verbose   &> 16S_core_metrics_samples.log
# I used the alpha group significance plugin to test for differences in alpha diversity:

# Run plugin for each alpha diversity result
for result in work/core_metrics/*vector.qza; do
  outname=${result/_vector.qza/_group_significance.qzv};
  qiime diversity alpha-group-significance --i-alpha-diversity $result --m-metadata-file ../FCCHHc2.5.map.txt --o-visualization $outname;
done
# Beta rarefaction plots were made to view the effect of multiple rarefactions on the data:

qiime diversity beta-rarefaction --i-table work/phylogeny/final_table_no_mce.qza --p-metric weighted_unifrac --p-clustering-method nj --m-metadata-file ../FCCHHc2.5.map.txt --p-sampling-depth 16386 --i-phylogeny work/phylogeny/rooted_tree.qza --o-visualization work/beta_rarefaction/weighted_unifrac.qzv
# Beta diversity plots came from the core metrics plugin.
#
# Any time I needed to add metadata to a beta diversity plot, I used the same PCoA matrix each time (rather than rarefying again) and just recreated the Emperor plot with the new metadata, for example:

# qiime emperor plot # --i-pcoa jaccard_pcoa_results.qza # --m-metadata-file ../metadata_for_qiime2_with_blasto_ANCOM_families.txt # --o-visualization jaccard_emperor_with_blasto_ANCOM_families.qzv

cd ../18S
pwd
echo Starting qiime2 on 18S
#Verify RAW amplicon sequences are in 00_Raw_Reads directory with no other files.

mkdir work
mkdir work/alpha_rarefaction
mkdir work/beta_rarefaction
echo Importing fastq files into a qiime2 qza file format

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path 00_Raw_Reads --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path work/demux-paired-end.qza

echo Trimming primers from amplicons

# For EMP 18S 1391f-EukBr primers. amplicon size 210-
# 1391f: GTACACACCGCCCGTC
# EukBr: TGATCCTTCTGCAGGTTCACCTAC
qiime cutadapt trim-paired --i-demultiplexed-sequences work/demux-paired-end.qza --p-cores 12 --p-front-f GTACACACCGCCCGTC --p-front-r TGATCCTTCTGCAGGTTCACCTAC --p-discard-untrimmed --o-trimmed-sequences work/primer-trimmed-demux-PE.qza --verbose &> 18S_primer_trimming.log

qiime demux summarize --i-data work/primer-trimmed-demux-PE.qza --o-visualization work/primer-trimmed-demux-PE.qzv

echo Visualize work/primer-trimmed-demux-PE.qzv in https://view.qiime2.org/ look for the basepair position where the median quality score is consistently less than 30.
echo enter basepair position for forward read on line 244 and bp position for reverse read on line 245.
echo ... or not, and just keep the default parameters.
echo ...
echo ...
echo ...
echo Truncating reads for quality using DADA2
# Run DADA2 denoising and filtering

#18S:
qiime dada2 denoise-paired --i-demultiplexed-seqs work/primer-trimmed-demux-PE.qza --p-trunc-len-f 87 --p-trunc-len-r 87 --p-n-threads 0 --output-dir work/DADA2_denoising_output --verbose &> 18S_DADA2_denoising.log

echo Visualizing DADA2 outputs

# Denoising stats
qiime metadata tabulate --m-input-file work/DADA2_denoising_output/denoising_stats.qza --o-visualization work/DADA2_denoising_output/denoising_stats.qzv

# Representative sequences
qiime feature-table tabulate-seqs --i-data work/DADA2_denoising_output/representative_sequences.qza --o-visualization work/DADA2_denoising_output/rep_seqs.qzv

# Feature table
qiime feature-table summarize --i-table work/DADA2_denoising_output/table.qza --o-visualization work/DADA2_denoising_output/table.qzv

echo Classify the representative sequences

## Classify the 18S representative sequences
qiime feature-classifier classify-sklearn --i-classifier /mnt/d/HBOI/FCCHH/qiime2/18S_qiime2_wd/silva-138-NR99-18Sclassifier.qza --i-reads work/DADA2_denoising_output/representative_sequences.qza --output-dir work/classified_sequences --verbose &> classify_18S_rep_seqs.log

# Tabulate the features, their taxonomy and the confidence of taxonomy assignment
qiime metadata tabulate --m-input-file work/classified_sequences/classification.qza --o-visualization work/classified_sequences/taxonomy.qzv

# The pipeline requires only the rep seqs, not the classified rep seqs
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences work/DADA2_denoising_output/representative_sequences.qza --output-dir work/phylogeny --p-n-threads 12 --verbose &> phylogenetic_tree_generation.log

echo making initial barplot
qiime taxa barplot --i-table work/DADA2_denoising_output/table.qza --i-taxonomy work/classified_sequences/classification.qza --m-metadata-file ../FCCHHc2.5.map.txt --o-visualization work/phylogeny/barplots_initial.qzv

echo remove low abundant features less than 10
qiime feature-table filter-samples --i-table work/DADA2_denoising_output/table.qza --p-min-frequency 10 --o-filtered-table work/phylogeny/frequency-filtered-table.qza

echo Removing chloroplast, mitochondria, and bacteria from 18S data
  qiime taxa filter-table --i-table work/phylogeny/frequency-filtered-table.qza --i-taxonomy work/classified_sequences/classification.qza --p-exclude chloroplast,mitochondria,bacteria --o-filtered-table work/phylogeny/final_table_no_mcb.qza

  qiime taxa filter-seqs --i-sequences work/DADA2_denoising_output/representative_sequences.qza --i-taxonomy work/classified_sequences/classification.qza --p-exclude chloroplast,mitochondria,bacteria --o-filtered-sequences work/phylogeny/rep_seqs_no_mcb.qza

echo Generating final barplots
  qiime taxa barplot --i-table work/phylogeny/final_table_no_mcb.qza --i-taxonomy work/classified_sequences/classification.qza --m-metadata-file ../FCCHHc2.5.map.txt --o-visualization work/phylogeny/barplots_no_mcb_final.qzv

  # Also check remaining read counts
  qiime feature-table summarize --i-table work/phylogeny/final_table_no_mcb.qza --o-visualization work/phylogeny/feature_table_no_mcb_final.qzv --m-sample-metadata-file ../FCCHHc2.5.map.txt

echo qiime2 is finished with 18S taxonomy
echo starting stats
# Alpha and beta diversity
# The commands I describe show how these steps were generally carried out for the 16S rRNA and 18S rRNA datasets.

# Alpha rarefaction was done as follows:
echo The maximum depth parameter on line 344 can be changed by visualizing the work/DADA2_denoising_output/table.qva
echo and move the slide bar to determine ideal depth for retaining most features and samples.

# 18S rRNA dataset
qiime diversity alpha-rarefaction --i-table work/phylogeny/final_table_no_mcb.qza --i-phylogeny work/phylogeny/rooted_tree.qza --p-max-depth 18821 --m-metadata-file ../FCCHHc2.5.map.txt --o-visualization work/alpha_rarefaction/rarefaction_18821.qzv
# The core metrics pipeline produces results of standard alpha and beta diversity metrics:

qiime diversity core-metrics-phylogenetic --i-phylogeny work/phylogeny/rooted_tree.qza --i-table work/phylogeny/final_table_no_mcb.qza --p-sampling-depth 18821 --m-metadata-file ../FCCHHc2.5.map.txt --output-dir work/core_metrics --p-n-jobs-or-threads 12 --verbose &> 18S_core_metrics_samples.log
# I used the alpha group significance plugin to test for differences in alpha diversity:

# Run plugin for each alpha diversity result
for result in work/core_metrics/*vector.qza; do
  outname=${result/_vector.qza/_group_significance.qzv};
  qiime diversity alpha-group-significance --i-alpha-diversity $result --m-metadata-file ../FCCHHc2.5.map.txt --o-visualization $outname;
done
# Beta rarefaction plots were made to view the effect of multiple rarefactions on the data:

qiime diversity beta-rarefaction --i-table work/phylogeny/final_table_no_mcb.qza --p-metric weighted_unifrac --p-clustering-method nj --m-metadata-file ../FCCHHc2.5.map.txt --p-sampling-depth 18821 --i-phylogeny work/phylogeny/rooted_tree.qza --o-visualization work/beta_rarefaction/weighted_unifrac_18821.qzv
# Beta diversity plots came from the core metrics plugin.

# Any time I needed to add metadata to a beta diversity plot, I used the same PCoA matrix each time (rather than rarefying again) and just recreated the Emperor plot with the new metadata, for example:

# qiime emperor plot # --i-pcoa work/core_metrics/jaccard_pcoa_results.qza # --m-metadata-file ../metadata_for_qiime2_with_blasto_ANCOM_families.txt # --o-visualization jaccard_emperor_with_blasto_ANCOM_families.qzv
#

echo qiime2 is done.
