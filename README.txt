--------------------------------------------------------------------------------
FAVR - Filtering and Annotation of Variants that are Rare
--------------------------------------------------------------------------------

Version: 1.0

Authors: Tu Nguyen-Dumont (1), Fabrice Odefrey (1),
         Melissa C Southey (1), Daniel J Park (1),
         Bernard J Pope (2).

         (1) Genetic Epidemiology Laboratory, Department of Pathology.
         (2) Victorian Life Sciences Computation Initiative.

         The University of Melbourne, Australia.

Web:     https://github.com/bjpop/favr

License: BSD

Requirements: Python 2.x (where x is >= 6), and the PySam library.

--------------------------------------------------------------------------------
General description
--------------------------------------------------------------------------------

Characterizing genetic diversity through the analysis of massively parallel
sequence (MPS) data offers enormous potential in terms of our understanding of
predisposition to complex human disease. Great challenges remain, however,
regarding our ability to resolve those genetic variants that are genuinely
associated with disease from the millions of "bystanders" and artefactual
signals. FAVR is designed to assist in the resolution of some of these issues
in the context of rare germline variants by facilitating "platform-steered"
artefact filtering. FAVR is a suite of tools that allow:

   (i) favr_35s_filter.py:

       Filtering of artefacts derived from "imbalanced" paired end sequencing,
       such that variants are filtered out if they only occur on 35 base reads
       (from SOLiD4 50-35 paired end chemistry).

  (ii) favr_nonfamily_filter.py:

       Filtering of variants based on their presence/absence and
       abundance in samples from non-family members.

 (iii) favr_family_annotate.py:

       Annotation of variants based on their presence/absence in
       samples from family members.

 (iv)  favr_refgene_annotate.py:

       Annotation based on RefGene co-ordinates of genetic features.

A note on command-line arguments and file names:

Each of the programs accepts command-line arguments, some of which are file
names. In all instances where a file name is expected you may use either a
relative or absolute path. Relative paths are interpreted "relative" to the
directory from which the program was executed, whereas absolute paths are
not, and must be specified in full. If the file you wish to specify resides
in the same directory from which the program is run then a relative path
would be just the name of the file on its own. This is likely to be the
simplest way of referring to the file.

When the input is a list of variants, we assume it is stored in tab separated
format (TSV).

--------------------------------------------------------------------------------
favr_35s_filter
--------------------------------------------------------------------------------

Filter variants based on whether they appear exclusively on 35 base reads
(from SOLID paired-end reads). Variants which appear only on 35 base reads are
binned. Other variants are kept.

Command line usage:

   ./favr_35s_filter.py
      [-h | --help]
      --variants=<variant list as TSV file>
      --bam=<bam file of reads for the same sample as variants>
      --bin=<bin filename>
      --keep=<keep filename>
      --log=<log filename>

Explanation of the arguments:

   --variants=<variant list>
      same as the favr_family_annotate.py tool (described above).

   --bam=<bam file of reads for the same sample as variants>

      This is the bam file of reads for the same sample that generated the
      list of variants.

   --bin=<bin filename>

      The bin file is an output of the program that contains the variants
      which were filtered out (binned).

   --keep=<keep filename>

      The keep file is an output of the program that contains the variants
      which were kept (not binned).

   --log=<log filename>

      The logfile records the reasons why each variant was either binned
      or kept.

--------------------------------------------------------------------------------
favr_nonfamily_filter
--------------------------------------------------------------------------------

Filter rare variants by comparing to samples from different families. The
basic idea is that, for each variant in the input, we are looking for evidence
of the variant in samples (bam files) of other non-family members. If there is
enough evidence of the variant in non-family members then we filter out (bin)
the variant because it is not sufficiently rare. The presence (or absence) of
evidence for a variant is based on two parameters: varLikeThresh and
samplesPercent.

The rule for finding evidence for a particular variant is based on what
percentage of total samples are "variant-like". A sample is considered to be
variant-like if it contains the same base (in the same position) as the
variant in at least "varLikeThresh" number of reads.

A variant is binned if "samplesPercent" percentage of the input samples are
variant-like.

Command line usage:

   ./favr_nonfamily_filter.py
      [-h | --help]
      --variants=<variant list as TSV file>
      --bin=<bin filename>
      --keep=<keep filename>
      --log=<log filename>
      --varLikeThresh=<variant read threshold>
      --samplesPercent=<percent of total samples which pass the threshold>
      reads1.bam reads2.bam ...

Explanation of the arguments:

   --variants=<variant list>
      same as the favr_family_annotate.py tool (described above).

   --bin=<bin filename>

      The bin file is an output of the program that contains the variants
      which were filtered out (binned).

      The format of the bin file is zero or more groups of lines of the form:

      chr1:1897857
         <vars/coverage: 1/9>
         <vars/coverage: 0/4>

      The group starts with the coordinate of the variant. It is followed
      by N lines of the form: <vars/coverage: X/Y>. Each one of these
      corresponds to one of the sample bam files,
      in the same order they were specified on the command line. X equals
      the number of reads in the sample whose base was the same as the
      variant (at the same position as the variant). Y equals the total
      number of reads in the sample which covered the same position as the
      variant. So a ratio of 1/9 means: at the same position as the
      variant, 1 base in the sample was the same as the variant base, and
      9 reads in the sample covered that position.

   --keep=<keep filename>

      The keep file is an output of the program that contains the variants
      which were kept (not binned).

      The keepfile is the same format as the input variants file, and
      is sorted on chromosome coordinates.

   --log=<log filename>

      The logfile records the reasons why each variant was either binned
      or kept, based on the output of the classify function.

   --varLikeThresh=<variant read threshold>

      Threshold for number of reads in a sample which have the same base
      as the variant (in the same position). Samples with a number of
      "variant" reads greater than or equal to varLikeThresh are considered
      variant-like.

   --samplesPercent=<percent of total samples which pass the threshold>

      A variant is binned (filtered out) if this many percent of the samples
      (bam files) are variant-like.

   reads1.bam reads2.bam ...
      same as the favr_family_annotate.py tool (described above).

--------------------------------------------------------------------------------
favr_family_annotate
--------------------------------------------------------------------------------

Annotate rare variants by comparing to samples from the same family.

Command line usage:

   ./favr_family_annotate.py
      [-h | --help]
      --variants=<variant list>
      --annotations=<output TSV file with annotations added>
      reads1.bam reads2.bam ...

Explanation of the arguments:

   --variants=<variant list>

      List of variants, one per line. Lines that start with a
      "Residue Based Coordinate System" (COMMA SEPARATED) are considered
      variants, other lines are ignored. The Residue Based Coordinate System has
      the form:

          chromosome,coordinate,orientation,alleles

      for example:

          22,30163533,1,A/C

      that is: chromosome 22, position 30163533, orientation 1 (forward),
      reference base 'A', variant base 'C'.

      This format is supported by various genomics tools including SIFT.

      We assume the position is 1 based, which means that the first coordinate
      in a chromosome is 1 (not zero).

      The rest of the line after the coordinates can be any text, which
      will be preserved in the output.

   --annotations=<output TSV file with annotations added>

      Save the annotated variants to this file. Annotations are appended to
      each variant line. Variants which are found
      in at least one family member sample are annotated with 'IN RELATIVE',
      whereas variants which are not found in at least one family member
      are annotated with 'NOT IN RELATIVE'.

   reads1.bam reads2.bam ...

      A list of bam files containing aligned sequence reads for
      relatives. Each bam file must be accompanied by
      an index (.bai) file (but you don't mention those index
      files on the command line).

--------------------------------------------------------------------------------
favr_refgene_annotate
--------------------------------------------------------------------------------

Annotate variants based on their position relative to exon features in the
genome reference.

   ./favr_refgene_annotate.py
      [-h | --help]
      --variants=<variant list as TSV file>
      --startslack=<distance from start of coding region>
      --spliceslack=<distance from exon start/end sites>
      --refGene=<refGene.txt file>
      --output=<output file name>

 Below is a diagram of a typical gene:

 5' -- forward strand --> 3'

      B      C   D           E     F    G    H          I 
      --------   -------------     ------    ------------
 A    |      |   |           |     |    |    |          |  J   K        L
 ------      |   |           |     |    |    |          ----   ----------
 ------  e1  |   |     e2    |     | e3 |    |    e4    ----   --- e5 ---
      |      |   |           |     |    |    |          |
      --------   -------------     ------    ------------

    On the forward strand, the symbols mean:
    ----------------------------------------
    e1..e5    5 exons
    A         transcription start site AND exon e1 start coordinate
    B         coding region start site
    C         exon e1 end coordinate
    D         exon e2 start coordinate
    E         exon e2 end coordinate
    F         exon e3 start coordinate
    G         exon e3 end coordinate
    H         exon e4 start coordinate
    I         coding region end site
    J         exon e4 end coordinate
    K         exon e5 start coordinate
    L         transcritpion end site AND exon e5 end coordinate

    Exons e1 and e4 are considered "partial coding" because only part of their
    extent lies within the coding region. Exon e5 is "non coding" because none
    of its extent lies within the coding region. Whereas exons e2 and e3 are
    "coding" exons because their entire extent lies within the coding region.

    If the gene was on the reverse (-) strand then the start and end positions
    would swap over. That is, the 5' end would be on the right hand side of
    the diagram, L would be the transcription start site and A would be the
    transcription end site, and so forth.

Explanation of the arguments:

   --variants=<variant list>
      same as the favr_family_annotate.py tool (described above).

   --startslack=<distance from start of coding region>

      Annotate variants which are within startslack base positions
      before the start of a coding region.

   --spliceslack=<distance from exon start/end sites>

      Annotate variants which are within +/- base positions
      of an exon start or end boundary.

   --refGene=<refGene.txt file>

      The refGene.txt coordinates file from UCSC.

   --output=<output file name>

      The annotated list of variants as output. The possible annotations are:
         - Within <startslack> before coding region start (towards the 5' end
           of the strand).
         - Within +/- <spliceslack> of coding exon start boundary
         - Within +/- <spliceslack> of coding exon end boundary
         - Within +/- <spliceslack> of NON-coding exon start boundary
         - Within +/- <spliceslack> of NON-coding exon end boundary
         - Within +/- <spliceslack> of PARTIAL-coding exon start boundary
         - Within +/- <spliceslack> of PARTIAL-coding exon end boundary

      Note "start" and "end" are determined by the strand on which the gene
      is located.
