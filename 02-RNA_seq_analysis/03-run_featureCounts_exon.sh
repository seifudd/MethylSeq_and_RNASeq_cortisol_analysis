#!/bin/sh -l
#$ -cwd

#SBATCH --job-name="FeatureCounts"
#SBATCH --partition=largemem
#SBATCH --time=12:00:00
#SBATCH --mem=64g
#SBATCH --cpus-per-task=4

. /data/NHLBI_BCB/Mehdi_P/RLee/RLee-RNA-seq-rlee-v3/richardlee_136891/04-featurecount/config.sh

SAMPLEID=$1
OUTPUTDIR=$2
OUTPUTDIR_HISAT2=$3

featureCounts -T $NUMCPUS \
	-t exon \
	-g gene_id \
	-a $HUMAN_GTF \
	-s 0 -p -O -f \
	-o "$OUTPUTDIR/$SAMPLEID.genefeatureCounts.txt" \
	"$OUTPUTDIR_HISAT2/$SAMPLEID/$SAMPLEID.out.sorted.bam.unique.bam"

: <<'END'
Version 1.4.6-p3

Usage: featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2] ... 

    Required parameters:

    -a <input>	Give the name of the annotation file. The program assumes
              	that the provided annotation file is in GTF format. Use -F
              	option to specify other annotation formats.
    
    -o <input>	Give the name of the output file. The output file contains
              	the number of reads assigned to each meta-feature (or each
              	feature if -f is specified). A meta-feature is the aggregation
              	of features, grouped by using gene identifiers. Please refer
              	to the users guide for more details.
    
   input_files	Give the names of input read files that include the read
              	mapping results. Format of input files is automatically
              	determined (SAM or BAM). Paired-end reads will be
              	automatically re-ordered if it is found that reads from the
              	same pair are not adjacent to each other. Multiple files can
              	be provided at the same time.
    
    Optional parameters:
    
    -A <input>	Specify the name of a file including aliases of chromosome
              	names. The file should be a comma delimited text file that
              	includes two columns. The first column gives the chromosome
              	names used in the annotation and the second column gives the
              	chromosome names used by reads. This file should not contain
              	header lines. Names included in this file are case sensitive.
    
    -F <input>	Specify the format of the annotation file. Acceptable formats
              	include `GTF' and `SAF'. `GTF' by default. Please refer to the
              	users guide for SAF annotation format.
    
    -t <input>	Specify the feature type. Only rows which have the matched
              	matched feature type in the provided GTF annotation file
              	will be included for read counting. `exon' by default.
    
    -g <input>	Specify the attribute type used to group features (eg. exons)
              	into meta-features (eg. genes), when GTF annotation is provided.
              	`gene_id' by default. This attribute type is usually the gene
              	identifier. This argument is useful for the meta-feature level
              	summarization.
    
    -f        	If specified, read summarization will be performed at the 
              	feature level (eg. exon level). Otherwise, it is performed at
              	meta-feature level (eg. gene level).
    
    -O        	If specified, reads (or fragments if -p is specified) will
              	be allowed to be assigned to more than one matched meta-
              	feature (or feature if -f is specified). 
    
    -s <int>  	Indicate if strand-specific read counting should be performed.
              	It has three possible values:  0 (unstranded), 1 (stranded) and
              	2 (reversely stranded). 0 by default.
    
    -M        	If specified, multi-mapping reads/fragments will be counted (ie.
              	a multi-mapping read will be counted up to N times if it has N
              	reported mapping locations). The program uses the `NH' tag to
              	find multi-mapping reads.
    
    -Q <int>  	The minimum mapping quality score a read must satisfy in order
              	to be counted. For paired-end reads, at least one end should
              	satisfy this criteria. 0 by default.
    
    -T <int>  	Number of the threads. 1 by default.
    
    -R        	Output read counting result for each read/fragment. For each
              	input read file, read counting results for reads/fragments will
              	be saved to a tab-delimited file that contains four columns
              	including read name, status(assigned or the reason if not
              	assigned), name of target feature/meta-feature and number of
              	hits if the read/fragment is counted multiple times. Name of
              	the file is the same as name of the input read file except a
              	suffix `.featureCounts' is added.
    
    --minReadOverlap <int>      Specify the minimum number of overlapped bases
              	required to assign a read to a feature. 1 by default. Negative
              	values are permitted, indicating a gap being allowed between a
              	read and a feature.
    
    --readExtension5 <int>      Reads are extended upstream by <int> bases from
              	their 5' end.
    
    --readExtension3 <int>      Reads are extended upstream by <int> bases from
              	their 3' end.
    
    --read2pos <5:3>            The read is reduced to its 5' most base or 3'
              	most base. Read summarization is then performed based on the
              	single base which the read is reduced to.
    
    --fraction	If specified, a fractional count 1/n will be generated for each
              	multi-mapping read, where n is the number of alignments (indica-
              	ted by 'NH' tag) reported for the read. This option must be used
              	together with the '-M' option.
    
    --primary 	If specified, only primary alignments will be counted. Primary
              	and secondary alignments are identified using bit 0x100 in the
              	Flag field of SAM/BAM files. All primary alignments in a dataset
              	will be counted no matter they are from multi-mapping reads or
              	not ('-M' is ignored). 
    
    --ignoreDup                 If specified, reads that were marked as
              	duplicates will be ignored. Bit Ox400 in FLAG field of SAM/BAM
              	file is used for identifying duplicate reads. In paired end
              	data, the entire read pair will be ignored if at least one end
              	is found to be a duplicate read.
    
    --countSplitAlignmentsOnly  If specified, only split alignments (CIGAR
              	strings containing letter `N') will be counted. All the other
              	alignments will be ignored. An example of split alignments is
              	the exon-spanning reads in RNA-seq data.
    
    Optional paired-end parameters:
    
    -p        	If specified, fragments (or templates) will be counted instead
              	of reads. This option is only applicable for paired-end reads.
              	The two reads from the same fragment must be adjacent to each
              	other in the provided SAM/BAM file.
    
    -P        	If specified, paired-end distance will be checked when assigning
              	fragments to meta-features or features. This option is only
              	applicable when -p is specified. The distance thresholds should
              	be specified using -d and -D options.
    
    -d <int>  	Minimum fragment/template length, 50 by default.
    
    -D <int>  	Maximum fragment/template length, 600 by default.
    
    -B        	If specified, only fragments that have both ends 
              	successfully aligned will be considered for summarization.
              	This option is only applicable for paired-end reads.
    
    -C        	If specified, the chimeric fragments (those fragments that 
              	have their two ends aligned to different chromosomes) will
              	NOT be included for summarization. This option is only 
              	applicable for paired-end read data.
    
    -v        	Output version of the program.
    
    --donotsort   If specified, paired end reads will not be reordered even if
              	reads from the same pair were found not to be next to each other
              	in the input. 
END

