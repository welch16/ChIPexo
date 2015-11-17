
## the E Coli genome is in:

/p/keles/genome_data/EColi_U00096.2/fasta/E.coli_K-12_MG1655_GENOME.fa


## the old ChIP-exo alignment is

@PG	ID:Bowtie	VN:0.12.5	

CL:"bowtie /home/GLBRCORG/iong/MG1655_data/bowtie_index/ecoli_NC_000913.2 /home/GLBRCORG/iong/mds/931__1002_Ecoli_RLEcoli2493_MOPSminGlu_37C_Sigma70_ChIPexo__ChIPseq/run299.456-166-1002-Sig-70-Lib_ACAGTG_L007_R1.fastq -m 1 -l 55 -k 1 -5 3 -3 40 --best -S /home/GLBRCORG/iong/mdsp/931__1002_Ecoli_RLEcoli2493_MOPSminGlu_37C_Sigma70_ChIPexo__ChIPseq/run299.456-166-1002-Sig-70-Lib_ACAGTG_L007_R1_bowtie_trim-5_3_-3_40.sam -X 500"

## the new alignment was done with bwa

bwa samse E.coli_K-12_MG1655_GENOME.fa edsn931_042814_Left.sai p2.run299.456-166-1002-Sig-70-Lib_ACAGTG_L007_R1.fastq

or

bwa samse E.coli_K-12_MG1655_GENOME.fa edsn1318_042814_Left.sai p2.run423.edsn-1318-cult-1202-ChIPexo-0minRif-BetaPrimeFLAG_GATCAG_L001_R1.fastq

#####

The complete command to align with bwa was something of the type:

# 1 - create index from the same genome we have

bwa index E.coli_K-12_MG1655_GENOME.fa

# 2 - Align reads to index

bwa aln E.coli_K-12_MG1655_GENOME.fa p2.run423.edsn-1313-cult-1197-ChIPexo-20minRif-Beta_TGACCA_L001_R1.fastq > edsn1313_042814_Left.sai

# 3 - this is what appears in bam file header

bwa samse E.coli_K-12_MG1655_GENOME.fa edsn1313_042814_Left.sai p2.run423.edsn-1313-cult-1197-ChIPexo-20minRif-Beta_TGACCA_L001_R1.fastq > edsn1313_042814_qc.unsorted.sam 

# 4 - sort and create index

samtools view -uS -t E.coli_K-12_MG1655_GENOME.fa edsn1313_042814_qc.unsorted.sam | samtools sort - edsn1313_042814_qc.sorted

samtools index edsn1313_042814_qc.sorted.bam

more on bwa is [this
slides](https://www.broadinstitute.org/files/shared/mpg/nextgen2010/nextgen_li.pdf)

#####

The Bowtie source and binary packages come with a pre-built index of the E. coli genome, and a set of 1,000 35-bp reads simulated from that genome. To use Bowtie to align those reads, issue the following command. If you get an error message "command not found", try adding a ./ before the bowtie.

bowtie e_coli reads/e_coli_1000.fq

### Building a new index

The pre-built E. coli index included with Bowtie is built from the
sequence for strain 536, known to cause urinary tract infections. We
will create a new index from the sequence of E. coli strain O157:H7, a
strain known to cause food poisoning. Download the sequence file by
right-clicking this link and selecting "Save Link As..." or "Save
Target As...". The sequence file is named NC_002127.fna. When the
sequence file is finished downloading, move it to the Bowtie install
directory and issue this command:

bowtie-build NC_002127.fna e_coli_O157_H7

The command should finish quickly, and print several lines of status
messages. When the command has completed, note that the current
directory contains four new files named e_coli_O157_H7.1.ebwt,
e_coli_O157_H7.2.ebwt, e_coli_O157_H7.rev.1.ebwt, and
e_coli_O157_H7.rev.2.ebwt. These files constitute the index. Move
these files to the indexes subdirectory to install it.  To test that
the index is properly installed, issue this command:

e_coli_O157_H7 <--- name

bowtie -c e_coli_O157_H7 GCGTGAGCTATGAGAAAGCGCCACGCTTCC

If the index is installed properly, this command should print a single
alignment and then exit.


##### CTCF alignment options from Pugh and Rhen ChIP-exo 2011 paper

Alignment to genome, peak calling, and data sharing. The Saccharomyces
reference genome was obtained from www.yeastgenome.org (build:
19-Jan-2007). The entire length of the sequenced tags were aligned to
the reference genome using Corona Lite software provided by the SOLiD
system, allowing up to 3 mismatches. This process was repeated for the
remaining tags, after removal of the 3’ most 6 bp, which tend to have
higher error rates. Raw sequencing data are available at NCBI Sequence
Read Archive (accession number: SRA044886). The resulting sequence
read distribution was used to identify peaks on the forward (W) and
reverse (C) strand separately using the peak calling algorithm in
GeneTrack . Processed data (the set of bound locations) can be found
in Supplementary Data File 1. For purposes of certain browser
displays, tag coordinates of the 5’ end were shifted by half the modal
towards the 3’ end. For ChIP-seq, the 5’ ends of all tags on the
forward and reverse strands were shifted 47 bp (half of the average
DNA fragment length) in the 3’ direction to maximize 5’-end overlap on
opposite strands, then plotted to be indistinguishable from their
source strand.
