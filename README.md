# NPP

npp is a java script intended for the screening of KEPs. It will cleave all input protein sequences at the residues of KK,KR,RK and RR.
Cleaved products are compared using blastp from blast+ command line toolkit to look for similar cleavage products (blast identity of 70%, resp 40% for very short sequences)
For more detailed Methods, please see: xxx


## Needed input file

Fasta file with protein sequences (Copci_AmutBmut1_GeneModels_FrozenGeneCatalog_20160912_aa.fasta)
NOTE: Header sequence to split fasta file into protein sequences has been hardcoded and would need to be adapted for other input files.


## Dependencies

Apache Commons Text library commons-text-1.0.jar [Apache](https://archive.apache.org/dist/commons/text/binaries/)
blastp from blast-2.7.1+ [blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)


## Software Versions

Script was originally run with blast-2.7.1+ and commons-text-1.0.jar.

