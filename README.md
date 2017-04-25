# UCE TO PROTEIN
A workflow for extracting protein-coding sequences from UCE data used in phylogenetics.

If you are using this program, please cite @@@:


@@@


This script uses Biopython:


Cock PA, Antao T, Chang JT, Bradman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B and de Hoon MJL (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 25:1422-1423.



And NCBI BLASTX:


Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL (2008) BLAST+: architecture and applications. BMC Bioinformatics 10:421.

The core of this workflow consists of three steps:
1) Using NCBI BLASTX for all sequences (taxa) within each UCE locus to find matches in a protein database,
2) Identifying best matches for each sequence for downstream analysis, and
3) Extraction of protein queries and their nucleotide equivalents from those hits.

More specifically, the initial requirements are unaligned FASTA files for each UCE locus produced by Phyluce and one or more protein sets from taxa closely related or conspecific with those present in the UCE data set. 

As the first step, a BLAST database needs to be constructed from the reference proteins. Following this, BLASTX is ran against this database for all sequences within each unaligned UCE locus file. This results in one BLASTX output per UCE locus, containing multiple matches (hits) for each UCE sequence (taxon). Each hit, in turn, may be composed of one or more ranges that correspond to protein fragments (exons). BLAST scores for those ranges are used for identification of best hits for each taxon and UCE. This is done with custom Python code using Biopython's module for parsing BLAST XML output and the following logic:

For each sequence (taxon) and hit within, both total and maximum scores are tallied. If a hit's total score is equal to the maximum score of its ranges and it corresponds to the maximum score for the taxon, such hit is considered best and is kept. This means that this hit was composed of a single range and its score was not exceeded by any other hits, composed of one or multiple ranges. If a hit's total score is higher than the maximum score of any one of its ranges, and that hit's total score is the best hit score for a species, this hit is kept unless it contains overlapping ranges. Finally, if the total hit score is equal to its maximum score but not to the best hit score for a species, its total score is checked to see if it corresponds to the highest individual range score for a species. If this is true, the hit is kept. Such hit would be composed of single range and considered best even if hits with higher total scores but overlapping ranges are present. If composed of multiple ranges, a best hit is concatenated into a query in the order based on its coordinates and its presence on either forward or reverse strand. These protein queries are then matched to corresponding input nucleotide sequence or its reverse complement.

If introns that do not change reading frame are present, translations of the query sequence may span across them. Because of this additional trimming is done if long (4 sites or more) gaps in the subject protein sequence are found. All sites corresponding to those long gaps are trimmed from the protein query and its nucleotide equivalent. If at this point there is still a stop codon in the protein query, such record is discarded.

For each record a protein query and its corresponding nucleotides are now considered ready for downstream analyses. For convenience, UCE locus name, taxon name, protein database hit name, trimmed nucleotide query, trimmed protein query, untrimmed protein query, and protein subject are written to an SQLite database.

## Installation and requirements

You can download this repository zipped (button on the right-hand side of the screen) and use `uce_to_protein.py` in this directory as a stand-alone program or clone it if you have [git installed](http://git-scm.com/book/en/v2/Getting-Started-Installing-Git) on your system.

This script requires you have Python version 3 or newer, Biopython, NCBI BLAST+, and SQLite3 installed.

Python 3 is available from [here](http://www.python.org/downloads/). On Ubuntu you can install it from the command line using

```
sudo apt-get install python3
```
See [Biopython project documentation](http://biopython.org/wiki/Documentation) for installing instructions.
BLAST can be downloaded from [here](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) or using command line on Ubuntu
```
sudo apt-get install ncbi-blast+
```
BLAST needs to be your system path such that you need to be able to call `blastx` and `makeblastdb` from the terminal.
Sqlite3 can be [downloaded here](http://www.sqlite.org/download.html) or, on Ubuntu, installed with 
```
sudo apt-get install sqlite3 libsqlite3-dev
```

# Interface
`uce_to_protein` can be run from the command line. Here is the general usage (you can view this in your command line with `uce_to_protein.py -h`):

```
usage: uce_to_protein <command> [<args>]
The commands are:
  blastdb        Create BLAST database
  queryblast     Query existing BLAST database and output BLASTX xml
  parse          Parse BLASTX xml output and create or update database of best hits
  queryprot      Query existing database of best hits
Use uce_to_protein <command> -h for help with arguments of the command of interest

positional arguments:
  command     Subcommand to run

optional arguments:
  -h, --help  show this help message and exit
```

To show help for individual commands, use `uce_to_protein.py <command> -h` or `uce_to_protein.py <command> --help`.

## Examples
This workflow for obtaining protein-coding sequences consists of four major steps, each of which is wrapped by the script:  

1) Creation of protein database to be used by BLAST, which can be done with `blastdb` command.

2) Querying of the BLAST database is done with `queryblast` command.

3) Parsing BLAST XML output and creation of SQLite database is done with `parse`.

4) Writing FASTA output of matched amino acid and corresponding nucleotide sequences can be done with `queryprot`.

### Creation of BLAST database
You will need some reference protein sequences to match your UCE data to. You can use one or concatenate more protein FASTA files for species closely related to the taxa for which you have sequenced UCEs.
```
uce_to_protein.py blastdb -i proteins.fasta 
```
Will create a database, by default called `prot_db` (you can change the default with `-o` or `-output`). This is a simple wrapper for MAKEBLASTDB. You can supply any additional arguments this program supports with `-a`, for example: `-a '-log makeblastdb.log'`. Not that the additional arguments must be provided as a literal string, quoted and complete with preceding dash(es).

### BLAST search
Once you have a database with proteins you wish to match against your UCEs, you need unaligned FASTA files of your UCE loci of interest. There are two considerations when preparing input FASTA files. First, they should have "clean" taxon names, that is without corresponding locus description. You can use `sed` or Phyluce's `phyluce_align_remove_locus_name_from_nexus_lines` script for that. Second, your files should follow `uce-name.unaligned.fasta` naming scheme or something similar. This script will chop off anything after the first period to determine your locus name and you will run into trouble if this means duplicates. You can match them against the protein database with the following:
```
uce_to_protein.py queryblast -i *unaligned.fasta
```
This assumes that you kept the default name `prot_db` from the first step, but you can change this with `-d`. This is the most time-consuming part of the process. You can use multiple cores for this step and/or freely customize the BLASTX search, which by default assumes E-value cutoff of 10e-5. See `uce_to_protein.py queryblast -h` for more options. This command will produce two kinds of output: one XML file per each FASTA input file and a configuration file, by default called `fasta_to_xml.conf`, needed for the next step. This configuration file has the following structure:
```
[fasta_to_xml]
uce-12001.unaligned.fasta = uce-12001.unaligned.fasta.xml
uce-12002.unaligned.fasta = uce-12002.unaligned.fasta.xml
uce-12003.unaligned.fasta = uce-12003.unaligned.fasta.xml
...
```
It will be used to determine match protein queries to nucleotide data in your input and assumes that both FASTA and XML files are in the same directory. You will need to modify this file if you move those files before proceeding to the next step.

### Parsing of BLAST results and retrieving best hits
In this step the script searches BLAST XML output for best hits, trims protein queries, matches them to nucleotide sequences from the original FASTA input, and puts all the data into a SQLite database.

The basic usage is this:

```
uce_to_protein.py parse -i fasta_to_xml.conf -o my_uce_hits.sqlite
```
As mentioned above, this command uses a configuration file to match nucleotide FASTA files and proteins in BLAST XML output. The script finds one best hit for each taxon in each locus and constructs a simple SQLite database where the information on locus name, taxon name, name of annotated protein matched in your reference, nucleotide query trimmed from introns, trimmed and untrimmed protein queries, and protein subject that was matched in your reference set. See `uce_to_protein.py parse -h` for description of options.

You can mine this database directly or use the `queryprot` command to write unaligned FASTA files for downstream use. If you specify an already existing database file with `-o`, it should be updated without introducing duplicates. This means that you can add data to it when you sequence additional taxa or loci. The locus and taxon names should be consistent between the data sets, however. Furthermore, if you are adding taxa, you may want to perform this and the previous step on FASTA files that contain only new taxa to speed up the process.   

This step is not parallelized and it takes several hours for ~2,500 loci and ~150 taxa on an example Hymeoptera dataset. 

### Writing FASTA output from SQLite database
Finally, the script provides a function that will quickly write out FASTA files from the SQLite database:
```
uce_to_protein.py queryprot -d my_uce_hits.sqlite -c taxon_config.conf -g taxon_group 
```
This command takes the name of your SQLite database and taxon configuration file that is essentially the same as the file used by Phyluce, where:
```
[set1]
Genus_species_1
Genus_species_2

[extended_set]
Genus_species_1
Genus_species_2
Genus_species_3
Genus_species_4

[other_taxa]
Genus_unkn
Genus_species_UG03
```
The output will be written for all loci in the database, as three FASTA files for each locus: 1) trimmed protein sequence found, 2) its corresponding nucleotide sequence, and 3) the nucleotide sequence with 3rd codon position removed. These files should be ready for alignment.

It is probably a good idea to further trim the alignments for both ambiguously aligned sites and taxa. This workflow is experimental. Use caution and examine your trimmed alignments visually for outliers using a fast and light-weight tool such as [Aliview](http://www.ormbunkar.se/aliview/).