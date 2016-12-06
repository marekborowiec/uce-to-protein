#! /usr/bin/env python3

"""
uce_to_protein.py is a script that provides an automated workflow 
to extract protein-coding sequences from Phyluce-derived datasets 
of ultraconserved elements.

Copyright (C) 2016 Marek Borowiec 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import argparse, re, sqlite3, subprocess, sys
from collections import defaultdict
from operator import itemgetter
from multiprocessing.dummy import Pool as dummyPool
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from xml.parsers.expat import ExpatError

__version__ = "0.1"


class ParsedArgs:

    def __init__(self):
        parser = argparse.ArgumentParser(
            usage='''uce_to_protein <command> [<args>]
The commands are:
  db      Create BLAST database
  query   Query existing database
  parse   Print out formatted output of xml BLASTX output
Use uce_to_protein <command> -h for help with arguments of the command of interest
'''
        )

        parser.add_argument(
            "command", 
            help="Subcommand to run"
        )

        # parse_args defaults to [1:] for args, but you need to
        # exclude the rest of the args too, or validation will fail
        self.args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, self.args.command):
            print('Unrecognized command')
            parser.print_help()
            exit(1)
        # use dispatch pattern to invoke method with same name
        getattr(self, self.args.command)()

    def db(self):
        """ command that wraps BLAST for creating protein database """
        parser = argparse.ArgumentParser(
            description="Create a BLAST database",
        )
        parser.add_argument(
            "-i",
            "--input",
            dest = "prot_input",
            help = "Input filename with reference protein sequences"
        )
        parser.add_argument(
            "-o",
            "--output",
            dest = "output_db_name",
            default = "prot_db",
            help = "File name for BLAST database created from reference sequences. Default: 'prot_db'"
        )

        args = parser.parse_args(sys.argv[2:])
        return args

    def query(self):
        """ command that wraps BLAST to query UCE alignment against protein database """
        parser = argparse.ArgumentParser(
            description="Query existing protein database with unaligned UCE sequences",
        )
        parser.add_argument(
            "-i",
            "--input-seqs",
            nargs = "+",
            type = str,
            required = True,
            dest = "query",
            help = "Input filenames with nucleotide sequences to be compared against reference"
        )
        parser.add_argument(
            "-d",
            "--database",
            dest = "input_db_name",
            default = "prot_db",
            help = "File name for BLAST database created from reference sequences. Default: 'prot_db'"
        )
        parser.add_argument(
        # parallelization can be used for BLAST
            "-c",
            "--cores",
            dest = "cores",
            default = 1,
            help = "Number of cores used. Default: 1"
        )

        args = parser.parse_args(sys.argv[2:])
        return args

    def parse(self):
        """ command to parse BLASTX output produced by query command """
        parser = argparse.ArgumentParser(
            description="Parse xml output of BLAST",
        )
        parser.add_argument(
            "-i",
            "--input-xml",
            nargs = "+",
            type = str,
            required = True,
            dest = "xml_files",
            help = "Input filenames with xml BLASTX output"
        )
        parser.add_argument(
            "-n",
            "--input-fasta",
            nargs = "+",
            type = str,
            required = True,
            dest = "fasta_files",
            help = "Input filenames with unaligned FASTA sequences"
        )

        args = parser.parse_args(sys.argv[2:])
        return args

    def get_args_dict(self):
        """ store arguments in a dictionary """
        command = self.args.__dict__
        arguments = getattr(self, self.args.command)().__dict__
        argument_dictionary = command.copy()
        argument_dictionary.update(arguments)
        
        return argument_dictionary

class DbCreator():
    """ given arguments, create a protein database using 'makeblastdb' """
    def __init__(self, **kwargs):
        self.command = kwargs.get("command")
        
def find_intron(string):
    """ find long gaps in string """
    matches = re.finditer("-{4,}", string)
    if matches:
        spans = []
        for match in matches:
            spans.append(match.span())
        return spans

def get_bits(alignment):
    """ get bits from all matches in alignment """
    return [(hsp.bits) for hsp in alignment.hsps]

def get_frames(alignment):
    """ get frames from all matches in alignment """
    return [hsp.frame for hsp in alignment.hsps]

def get_nt_extracted(sign, nt_string, starts, ends):
    start_end = zip(starts, ends)
    increasing_start_end = sorted(start_end, key=itemgetter(0)) # this is needed to extract from nt query
    seq = Seq(nt_string)
    new_nt = []
    for start, end in increasing_start_end:
        if sign == "plus":
            fragment = str(seq[(start - 1):end])
            new_nt.append(fragment)
        elif sign == "minus":
            fragment = str(seq[(start - 1):end])
            fragment = Seq(fragment)
            fragment = str(fragment.reverse_complement())
            new_nt.insert(0, fragment)

    return "".join(new_nt)

def get_sorted_queries(sign, nucleotides, alignment):
    """ get hsps queries and subject sorted by query start,
    either ascending or descenging, depending on strand """
    queries = [hsp.query for hsp in alignment.hsps]
    subjects = [hsp.sbjct for hsp in alignment.hsps]
    starts = [int(hsp.query_start) for hsp in alignment.hsps]
    ends = [int(hsp.query_end) for hsp in alignment.hsps]

    zipped = zip(starts, ends, queries, subjects)
    if sign == "plus":
        sorted_zipped = sorted(zipped, key=itemgetter(0))
    elif sign == "minus":
        sorted_zipped = sorted(zipped, key=itemgetter(0), reverse=True)
    sorted_starts, sorted_ends, sorted_queries, sorted_subjects = zip(*sorted_zipped)

    nt_query_string = get_nt_extracted(sign, nucleotides, starts, ends)
    query_string = ''.join(sorted_queries)
    subject_string = ''.join(sorted_subjects)

    return (nt_query_string, query_string, subject_string, sorted_starts, sorted_ends)

def add_to_highest_scoring_dict(highest_score_dict, nucleotides, key, alignment, check_length="no", query_length=0):
    """ assign items to dictionary of highest scoring hits """
    subject_gene_name = alignment.title
    bits = get_bits(alignment)
    frames = get_frames(alignment)
    if frames[0][0] > 0: # checking sign of first frame 
        sign = "plus"
    elif frames[0][0] < 0:
        sign = "minus"
    nt_query, query, subject, sorted_starts, sorted_ends = get_sorted_queries(sign, nucleotides, alignment)
    trimmed_queries = trim_introns_and_seqs_w_stop(nt_query, query, subject)
    try:
        (trimmed_nt_query, trimmed_query) = trimmed_queries
    except TypeError:
        trimmed_query = False
    if trimmed_query:
        if check_length == "no":
            highest_score_dict[key] = (trimmed_nt_query, trimmed_query, query, subject, 
             bits, frames, sorted_starts, sorted_ends, subject_gene_name)
        elif check_length == "yes":
            if len(query) >= 30 and len(query) <= query_length / 3:
                highest_score_dict[key] = (trimmed_nt_query, trimmed_query, query, subject, 
                 bits, frames, sorted_starts, sorted_ends, subject_gene_name)
            
def get_highest_scoring(records, nucleotide):
    """ get a dictionary of { species : 'best' hit } from parsed BLASTX output """
    # dictionary of highest cumulative scores for each hit ("alignment")
    nt_dict = {record.id : str(record.seq) for record in nucleotide} # { species : nt sequence } made from parsed FASTA file
    highest_scoring_dict = {}
    try:
        for item in records:
            query_len = item.query_letters
            species_dict = defaultdict(list)

            for alignment in item.alignments: # alignment here refers to BLAST hit
                species_dict[item.query].append(alignment)

            for species, alignments in species_dict.items():
                total_species_scores = []
                max_species_scores = []
                for alignment in alignments:
                    scores = [hsp.score for hsp in alignment.hsps]
                    if scores:
                        total_alignment_score = sum(scores)
                        max_alignment_score = max(scores)
                        total_species_scores.append(total_alignment_score)
                        max_species_scores.append(max_alignment_score)

                for alignment in alignments:
                    scores = [hsp.score for hsp in alignment.hsps]
                    nt = nt_dict[species]
                    if scores:
                        total_alignment_score = sum(scores)
                        max_alignment_score = max(scores)
                        # if total score equals max score and is also best total score, take it
                        if total_alignment_score == max_alignment_score and total_alignment_score == max(total_species_scores):
                            add_to_highest_scoring_dict(highest_scoring_dict, nt, species, alignment)
                        # if total score > max score, check if this is the best total score
                        elif total_alignment_score > max_alignment_score and total_alignment_score == max(total_species_scores):
                            # if the query is not too long, keep it
                            add_to_highest_scoring_dict(highest_scoring_dict, nt, species, alignment, check_length="yes", query_length=query_len)
                        # if total score equals max score but is not the best total score, check if it is best max score, then take it
                        # else skip it
                        elif total_alignment_score == max_alignment_score and max_alignment_score == max(max_species_scores):
                            add_to_highest_scoring_dict(highest_scoring_dict, nt, species, alignment)

    except ExpatError:
        print("Unexpected end of xml file...")

    return highest_scoring_dict

def trim_introns_and_seqs_w_stop(nt_query, query, subject):    
    """ given highest scoring concatenated query and subject,
    trim bases that correspond to long gaps in subject (likely introns)
    and delete sequences that still have stop codons """
    trimm_nt = []
    trimm_aa = []
    if find_intron(subject):
        introns = list(find_intron(subject))
        if len(introns) == 1:
            intron_start = introns[0][0]
            intron_end = introns[0][1]
            trimm_nt.append(nt_query[:(intron_start * 3)])
            trimm_nt.append(nt_query[(intron_end * 3):])
            trimm_aa.append(query[:intron_start])
            trimm_aa.append(query[intron_end:])
        elif len(introns) == 2:
            intron1_start = introns[0][0]
            intron1_end = introns[0][1]
            intron2_start = introns[1][0]
            intron2_end = introns[1][1]
            trimm_nt.append(nt_query[:(intron1_start * 3)])
            trimm_nt.append(nt_query[(intron1_end * 3):(intron2_start * 3)])
            trimm_nt.append(nt_query[(intron2_end * 3):])
            trimm_aa.append(query[:intron1_start])
            trimm_aa.append(query[intron1_end:intron2_start])
            trimm_aa.append(query[intron2_end:])
        elif len(introns) > 2:
            for index, intron in enumerate(introns[:-1]):
                if index == 0:
                    intron1_start = intron[0]
                    intron1_end = intron[1]
                    trimm_nt.append(nt_query[:(intron1_start * 3)]) # first exon
                    trimm_aa.append(query[:intron1_start])
                if index > 0:
                    current_intron_start = intron[0]
                    current_intron_end = intron[1]
                    next_intron_start = introns[index+1][0]
                    trimm_nt.append(nt_query[(current_intron_end * 3):(next_intron_start * 3)]) # intermediate exons
                    trimm_aa.append(query[current_intron_end:next_intron_start])
            last_intron_start = introns[-1][1]
            last_intron_end = introns[-1][1]
            trimm_nt.append(nt_query[(last_intron_end * 3):]) # last exon
            trimm_aa.append(query[last_intron_end:])

    else:
        trimm_nt.append(nt_query)
        trimm_aa.append(query)

    trimmed_nt = "".join(trimm_nt)
    trimmed_aa = "".join(trimm_aa)
        
    if "*" not in trimmed_aa:
        #print(coordinates)
    # exclude sequences that still have stop codon
        return (trimmed_nt, trimmed_aa)

def create_sqlite_db():
    """ create sqlite database """
    conn = sqlite3.connect("test.db")

    with conn:
        cur = conn.cursor()
        cur.execute("DROP TABLE IF EXISTS Uces")        
        cur.execute("DROP TABLE IF EXISTS Taxa")
        cur.execute("DROP TABLE IF EXISTS Hitnames")
        cur.execute("DROP TABLE IF EXISTS Sequences")

        cur.execute("CREATE TABLE {tn} ({uid} INTEGER PRIMARY KEY, {un} UNIQUE)".format(
         tn="Uces", uid="uce_name_ID", un="uce_name"))

        cur.execute("CREATE TABLE {tn} ({tid} INTEGER PRIMARY KEY, {txn} UNIQUE)".format(
         tn="Taxa", tid="taxon_name_ID", txn="taxon_name"))

        cur.execute("CREATE TABLE {tn} ({tid} INTEGER PRIMARY KEY, {txn} UNIQUE)".format(
         tn="Hitnames", tid="hit_name_ID", txn="hit_name"))

        cur.execute("""CREATE TABLE {tn} ({sid} INTEGER PRIMARY KEY, {uid} INT, {tid} INT, {hid} INT, {tq} TEXT, {q} TEXT, {s} TEXT,
         FOREIGN KEY({uid}) REFERENCES Uces({uid}),
         FOREIGN KEY({tid}) REFERENCES Taxa({tid}),
         FOREIGN KEY({hid}) REFERENCES Hitnames({hid}))
         """.format(
         tn="Sequences", sid="seq_ID", uid="uce_name_ID",
         tid="taxon_name_ID", hid="hit_name_ID",
         tq="trimmed_query", q="query", s="subject"))

def populate_sqlite_db(uce_name, seq_dict):
    """ add data from highest scoring dictionary to sqlite database """
    conn = sqlite3.connect("test.db")

    with conn:
        cur = conn.cursor()
        cur.execute("INSERT OR IGNORE INTO Uces(uce_name) VALUES (?)", (uce_name, ))
        cur.execute("SELECT uce_name_ID FROM Uces WHERE uce_name = ?", (uce_name, ))
        uce_ID = int(cur.fetchone()[0])

        for species, seq in seq_dict.items():
            (nt_trimmed_query, trimmed_query, query, subject, bits, frames, query_s, query_e, gene) = seq


            #print('File: {}\nSpecies: {}\nGene: {}\nNucleotides: {}\nQuery: {}\nTrimmed: {}\nSubject: {}\nFrames: {}\nQuery start: {}\nQuery end: {}'.format(
            # uce_name, species, gene, nt_trimmed_query, query, trimmed_query, subject, frames, query_s, query_e))

            cur.execute("INSERT OR IGNORE INTO Taxa(taxon_name) VALUES (?)", (species, ))
            cur.execute("INSERT OR IGNORE INTO Hitnames(hit_name) VALUES (?)", (gene, ))

            cur.execute("SELECT taxon_name_ID FROM Taxa WHERE taxon_name = ?", (species, ))
            taxon_ID = int(cur.fetchone()[0])
            cur.execute("SELECT hit_name_ID FROM Hitnames WHERE hit_name = ?", (gene, ))
            hit_ID = int(cur.fetchone()[0])

            cur.execute("INSERT INTO Sequences(uce_name_ID, taxon_name_ID, hit_name_ID, trimmed_query, query, subject) VALUES (?, ?, ?, ?, ?, ?)", 
             (uce_ID, taxon_ID, hit_ID, trimmed_query, query, subject))

def output_fasta(highest_scoring_dict, n=80):
    """ produce a FASTA string from dictionary of best hits """
    # each sequence line will have 80 characters 
    # for each element of species : seq dictionary,
    # split sequence into list of string, each n chars long
    # then join everything with newline 
    fasta_string = '\n'.join(['>{}\n{}'.format(species, '\n'.join([seq[1][i:i+n] for i in range(0, len(seq[1]), n)])) \
     for species, seq in sorted(highest_scoring_dict.items())])
    if fasta_string:
        return "{}\n".format(fasta_string) # add last end of line

def write_fasta(blast_file_name, fasta_string):
    """ write FASTA file """
    new_output_name = re.sub(".xml", "", blast_file_name)
    out_file_name = "protein-{}".format(new_output_name)
    with open(out_file_name, "w") as f:
        f.write(fasta_string)

def input_parse(xml_file_name, nt_file_name):
    """ parse XML BLASTX output file """
    with open(xml_file_name, "r") as f, open(nt_file_name, "r") as nt:
        records = NCBIXML.parse(f)
        nucleotide = SeqIO.parse(nt, "fasta")
        return get_highest_scoring(records, nucleotide)

def common_entries(*dcts):
    """ from StackOverflow """
    for i in set(dcts[0]).intersection(*dcts[1:]):
        yield (i,) + tuple(d[i] for d in dcts)

def output_parsed(blast_file_name, fasta_file_name):
    """ write output of query command """
    base_uce_name = fasta_file_name.split(".")[0] # move this somewhere else
    new_output_name = "{}.unaligned.fasta".format(base_uce_name)
    highest_scoring_dict = input_parse(blast_file_name, fasta_file_name)

    fasta_string = output_fasta(highest_scoring_dict)
    if fasta_string:
        print("Writing FASTA file protein-{}...".format(new_output_name))
        write_fasta(blast_file_name, fasta_string)
        populate_sqlite_db(base_uce_name, highest_scoring_dict)
    else:
        print("File {} contains no protein matches".format(blast_file_name))

def get_blast_call_string(in_file, db_name):
    """ make command line call string for BLASTX """
    call_string = "blastx -query {0} -db {1} -outfmt 5 -evalue 10e-5 > {0}.xml".format(in_file, db_name)
    return call_string

def call_blast(call_string):
    """ make command line call for BLASTX """
    a = re.search("-query (\S+)", call_string)
    b = re.search("-db (\S+)", call_string)
    file_name = a.group(1)
    database = b.group(1)
    print("Blasting file {} against protein database {}...".format(file_name, database))
    try:
        p = subprocess.call(call_string, shell=True)
    except KeyboardInterrupt:
        print("\nYou killed the query")
        sys.exit()

def main():

    # initialize parsed arguments and batabase creator object
    kwargs = run()
    db_creator = DbCreator(**kwargs)
       
    if db_creator.command == "db":
        call_string = "makeblastdb -in {} -dbtype prot -out {}".format(kwargs["prot_input"], kwargs["output_db_name"])
        subprocess.call(call_string, shell=True)

    if db_creator.command == "query":
        if int(kwargs["cores"]) == 1:
            for file in kwargs["query"]:
                blast_command = get_blast_call_string(file, kwargs["input_db_name"])
                call_blast(blast_command)

        if int(kwargs["cores"]) > 1:
            commands = [get_blast_call_string(file, kwargs["input_db_name"]) for file in kwargs["query"]]
            pool = dummyPool(int(kwargs["cores"]))
            try:
                pool.map(call_blast, commands)
            except KeyboardInterrupt:
                pool.close()
                pool.terminate()
                print("\nYou killed the query")
                sys.exit()

    if db_creator.command == "parse":
        create_sqlite_db()
        blast_fasta_pairs = zip(kwargs["xml_files"], kwargs["fasta_files"]) 
        for pair in blast_fasta_pairs:
            (blast_fn, fasta_fn) = pair 
            output_parsed(blast_fn, fasta_fn)

def run():

    # initialize parsed arguments
    config = ParsedArgs()
    # get arguments
    config_dict = config.get_args_dict()
    return config_dict
    
if __name__ == '__main__':
        
        main()