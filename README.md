# ExtendedSplitORFPipeline
This pipeline is an extension of the Split-ORF pipeline available at https://github.com/SchulzLab/SplitOrfs.

The pipeline will invoke the Split-ORF pipeline and use the resulting files to determine unique regions in the DNA and protein sequences of the predicted split-ORFs.

## Dependencies
In order for the pipeline to work correctly you need:


1) *NCBI blast installation*. The workflow uses a local NCBI Blast installation. It uses the makedb and blastp executables. The latest version can be downloaded [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi).
2) *bedtools2* (version v2.27.1 or higher). Only needed if the predicted protein regions should be overlapped with functional annotations.
3) *Python 2.7* All python scripts of the Split-ORF pipeline are written for Python 2.7.
4) *Anaconda* with environment py2 (Python 2.7) and py3_7 (Python 3.7). For more information see [here](https://docs.anaconda.com/anaconda/install/)
5) *Biopython* (version 1.77 or higher) Used to read in fasta files into python. The latest version can be downloaded [here](https://biopython.org/wiki/Download.)
6) *pybedtools* (version 0.8.1 or higher). Needed for uniqueness determination.
7) *MUMmer* (version 3.23). Needed for exact matching of the predicted ORF sequences against the reference files. The latest version can be downloaded [here](http://mummer.sourceforge.net/)
8) *R with Rmarkdown* (R version 3.6). Used to create end reports for both the Split-ORF pipeline and the extended pipeline. For more information see [here] (https://www.r-project.org/)


Blast, bedtools and mummer should be added to the unix PATH variable
Biopython and pybedtools only need to be available in the py3_7 environment.

## Usage
The run_pipeline.sh is the masterscript that executes the Split-ORF pipeline and the additional scripts that are used to determine the unique regions.
For further information on the split-ORF pipeline refer to https://github.com/SchulzLab/SplitOrfs.
For more information on the extended pipeline refer to run_pipeline.sh.

As input you need to give three fasta files and an annotation (bed) file:

1.The first is a Fasta file of the proteins that you want to use for the analysis (all protein-coding transcripts). It has the following format for the header: 

**>ENSG00000001626|ENST00000003084**

Where the first entry is the gene identifier and the second one the protein ID/transcript ID from which the protein was made. They need to be separated by a | character

2.The second is a Fasta file that contains transcript sequences. From these transcript sequences possible open reading frames (ORFs) are generated and translated to peptides and then used for alignment with BlastP against the proteins.

The format looks similar to above:

**>ENSG00000100767|ENST00000216658**

The only difference is that the second ID is the transcript ID that represents the transcript to be subjected to ORF generation.

3.The **annotation bed file** is used for checking with overlap of known protein domains. It has the following format (header only shown for illustration should not be in the file):

|Protein/Transcript ID   | Start   | End   | Identifier  |
|---|---|---|---|
|ENST00000308027| 21  |    274  |   PF07690|
|ENST00000574588|104 |    414   |  PF00038|

The first column denotes the Protein or Transcript ID representing the protein (here the Ensembl human transcipt ID of the protein).The second and third denote the start and end of the domain annotation in the protein. The last column is the identifier of the domain type (here PFAM domain). The annotation of human and mouse proteins can be found in the folder *annotations* in the Split-ORF directory in the repo.

4.The last is a Fasta file that contains the DNA sequence that you want to use as reference (all protein-coding transcipts). The predicted split-ORF sequences will be aligned to these sequences to determine uniqueness. Each region that does not match will be annotated as unique.

The format looks similar to above:

**>ENSG00000005801|ENST00000005082**

The default function call is:

```javascript
Usage: ./run_pipeline.sh [-h] proteins.fa transcripts.fa annotation.bed proteincodingtranscripts.fa
```

proteins.fa 			should be a multi fasta file containing the amino acid sequences of the proteins that are used as reference (whole transcriptome).
transcripts.fa 			should be a multi fasta file containing the DNA sequences of the reads/transcripts that shall be analyzed.
annotation.bed 			should be a bedfile containing the annotations for the used genome build.
				The standard annotation files for human and mouse (ENSEMBL 95) can be found in the annotations directory in the SplitORF directory.
proteinCodingTranscripts.fa 	should be a multi fasta file containing the DNA sequences of the protein coding transcripts that are used as reference.

## Output
For easier usability a folder named "run_day.month.year_hour.minute.second" is created in the "Output" folder everytime the pipeline is envoked and all output files produced by the pipeline are saved within.

The Split-ORF pipeline produces a number of files as output, some are just intermediates not of relevance. The relevant ones are:

1. OrfProteins.fa -  a fasta file of all generated proteins (in all three reading frames) from the transcripts in the transcripts.fa file supplied to runSplitOrfs.sh.
2. ProteinDatabase.phr, ProteinDatabase.psq, ProteinDatabase.pin - database files from Blast
3. UniqueProteinORFPairs.txt - the final set of transcripts, that have at least 2 ORFs matching to one of the proteins supplied in proteins.fa. Format explained below.
4. UniqueProteinORFPairs_annotated.txt - an extended file from above, when you also ad an annotation bed file to the pipeline.

The different columns of the UniqueProteinORFPairs_annotated.txt file are explained below.

|column name|explanation|
|---|---|
|**geneID** | Gene identifier|
|**targetTransID** | Identifier of the target protein/transcript that the ORFs have been aligned to.|
|**OrfTransID** | Identifier from which the ORFs have been generated.|
|**NumOrfs** | Number of ORFs matching to targetTransID.|
|**OrfIDs** | The unique ORF identifiers representing the ORFs that aligned (comma separated). These ORFs can be found in the file OrfProteins.fa in the output folder.|
|**OrfPos**| The nucleotide start-stop positions from which the ORF was generated (comma separated list for all matching ORFs).|
|**OrfLengths** | The nucleotide length of the matching ORFs (comma separated).|
|**OrfSeqIdents** | The sequence identity values of the ORF-protein alignments as reported by BlastP (comma separated).|
|**MinSeqIdent** | Minimal observed sequence identity of all the ORF-protein matches.|
|**MaxSeqIdent** | Maximal observed sequence identity of all the ORF-protein matches.|
|**protAlignPos** | Alignment start-stop positions of the ORF in the protein.|
|**ProtCoverage** | Number of amino acid positions covered of the original protein by alignment from all ORFs.|
|**ORF-DomainAnnot** | Identifiers of annotations that overlap with an ORF (comma separated list in order of the ORFs). NA means *not available*, when no ORF annotation existed.|
|**NumOrfAnnot** | Number of ORFs that have at least one overlapping annotation.|
|**AnnotPercent** | The ratio of ORFs that have at least one overlapping annotation.|

The extended pipeline also produces several output files. The relevant ones are:
1. Unique_DNA_regions.bed	-	a bedfile annotating the unique regions of the predicted split-ORF proteins (as DNA sequence)
2. Unique_DNA_regions.fa	-	a Fasta file containing the unique regions of the predicted split-ORF proteins (as DNA sequence)
3. Unique_Protein_regions.bed	-	a bedfile annotating the unique regions of the predicted split-ORF proteins (as AA sequence)
4. Unique_Protein_regions.fa	-	a Fasta file containing the unique regions of the predicted split-ORF proteins (as AA sequence)
5. Split-ORF_Report.html	-	A report with basic statistics of the split-ORF pipeline
6. Uniqueness_Report.html	-	A report with basic statistics of the uniqueness scripts
