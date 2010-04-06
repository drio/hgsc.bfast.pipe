------------------------------BAMStats.jar-------------------------------------

This application can be used to see the statistical information about a BAM/SAM
file. 

Typical information shown for a BAM/SAM file is
    a) Number of reads
    b) Number of mapped reads
    c) Percentage of mapped reads
    d) Number of unmapped reads
    e) Number of duplicate reads
    f) Number of reads with exactly one alignment
    g) Percentage of reads with exactly one alignment
    h) Throughput (Product of mapped read count and read length)
    i) Uniqueness percentage
    j) Effective throughput (Product of uniquely mapped read count and read length)
    
This application computes the statistical information for the specified read
types only, i.e., for a paired-read files, it can compute statistics for either
read 1 or read 2 (which must be specified in the input). 

Usage:
java -jar BAMStats.jar FileName ReadType SequencerType
	FileName             - BAM/SAM file to be analyzed
	ReadType             - 1:Read1, 2:Read2, 3:Fragment
	Sequencer {Optional} - Please specify "solid" for
	                     - solid reads, omit for other types

  * For BAMs generated using solid reads, this tool reports paired reads as R3 and F3.
    In this case, read1 is reported as R3 and read2 as F3.
  * For other BAMs, this tool reports the paired reads as read1 and read2.
   
Pre-requisites:

Before calculating statistics on a BAM/SAM file, it is essential that the BAM/
SAM file is sorted and duplicates entries are marked.

1) Sorting
It is required that the input BAM/SAM file be coordinate sorted. If this file is
not already coordinate sorted, the following command can be used to achieve that.

java -jar -Xmx4g /stornext/snfs1/next-gen/software/picard-tools/current/SortSam.jar \
I=$InputFile O=$outputFile SO=coordinate TMP_DIR="/space1/tmp" \
VALIDATION_STRINGENCY=SILENT

SortSam tool sorts the file specified by the input (I) parameter and writes the
sorted file as specified by the output (O). It needs a temporary directory to
write intermediate data. Thus, an appropriate directory should be specified as
the value for "TMP_DIR". If this is not specified, it attempts to write to the 
root directory.

2) Marking duplicates
Duplicate reads should be appropriately marked using Picard or similar tool.
This is necessary to identify the number of uniquely mapped reads and hence, 
determine uniqueness percentage and effective throughput. Use the following
command to mark the duplicates

java -jar Xmx4g /stornext/snfs1/next-gen/software/picard-tools/current/MarkDuplicates.jar \
I=$InputFile O=$OutputFile M=$MetricsFile TMP_DIR="/space1/tmp" \
VERBOSITY=ERROR VALIDATION_STRINGENCY=SILENT

At a high-level, the algorithm used to mark duplicate entries works as explained
below:

To consider the reads as candidates for marking as duplicates, their reference
index and the start (or end) position in the reference, orientation and the
library are compared. If these are all equal, the read with the highest score
is considered the original read and all the subsequent matching reads having
lower score are marked as duplicates.

Sample Output and Explanation of Results:

BAM/SAM File : input.bam
Read Type    : Fragment

Total Reads in file        : 98,546,410
Number of Reads considered : 98,546,410

Mapped Reads               : 55,011,114
Mapped Reads Percentage    : 55.82%
Throughput                 : 2,750,555,700 bp

Num Reads With 1 Alignment : 40,621,051
% Reads With 1 Alignment   : 41.22%

Duplicate Reads            : 6,604,430
Uniquely Mapped Reads      : 48,406,684
Uniqueness Percentage      : 87.99%
Effective Throughput       : 2,420,334,200 bp

Logging Mapping Quality Distribution to : input_log.txt

Logging Mapping Quality Distribution to : input_nonunique_log.txt

Computation Time      : 1076.837 sec

1) "Read Type" shows what type of reads are considered for the given input. Reads
could be either of type "Read 1", "Read 2", or "Fragment".

2) "Total Reads in file" represents the sum of all the reads in the input,
irrespective of their type.

3) "Number of Reads Considered" is the total number of reads of the specified type.

4) "Mapped Reads" is the number of mapped reads of the specified type, i.e., in the
example above, it represents the total number of "fragment" mapped reads.

5) "Mapped Reads Percentage" is the percentage of mapped reads with the the total number
of the reads of the specified type, i.e., with "Number of Reads Considered".

6) "Throughput" is the sum of read lengths of all mapped reads. It represents
the size (in base pairs) of the sequences produced.

7) "Num Reads With 1 Alignment" is the number of mapped reads having exactly one matching
alignment in the reference, i.e., number of mapped reads where value of NH tag is 1.

8) "% Reads With 1 Alignment" is the percentage of number of reads with 1 alignments with
the total number of reads considered.

9) "Duplicate Reads" represents the number of duplicate reads in the file. How duplicates
are marked is described in the earlier section.

10) "Uniquely Mapped Reads" is the number of mapped reads (of the specified type) that are
not duplicates. Note that this number also contains the number of reads with 1 alignment.

11) "Uniqueness Percentage" is the percentage of uniquely mapped reads with the total
number of reads of the specified type.

12) "Effective Throughput" is the sum of read lengths for all uniquely mapped reads. It
represents the size (in base pairs) of unique sequences produced. 
 