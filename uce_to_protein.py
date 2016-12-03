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
from multiprocessing.dummy import Pool as dummyPool
from Bio.Blast import NCBIXML
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

def get_query_string(alignment):
    """ get query string from all matches in alignment """
    return ''.join([hsp.query for hsp in alignment.hsps])

def get_match_string(alignment):
    """ get match string from all matches in alignment """
    return ''.join([hsp.match for hsp in alignment.hsps])

def get_subject_string(alignment):
    """ get subject string from all matches in alignment """
    return ''.join([hsp.sbjct for hsp in alignment.hsps])

def get_bits_string(alignment):
    """ get bits string from all matches in alignment """
    return ''.join([str(hsp.bits) for hsp in alignment.hsps])

def get_frames(alignment):
    """ get frames from all matches in alignment """
    return [str(hsp.frame) for hsp in alignment.hsps]

def get_query_starts(alignment):
    """ get query starts from all matches in alignment """
    return [str(hsp.query_start) for hsp in alignment.hsps]

def add_to_highest_scoring_dict(dictionary, key, alignment, check_length="no", query_length=0):
    """ assign items to dictionary of highest scoring hits """
    subject_gene_name = alignment.title
    query = get_query_string(alignment)
    subject = get_subject_string(alignment)
    bits = get_bits_string(alignment)
    frames = get_frames(alignment)
    query_starts = get_query_starts(alignment)
    if check_length == "yes":
            if len(query) >= 30 and len(query) <= query_length / 3:
                dictionary[key] = (query, subject, bits, frames, query_starts, subject_gene_name)
    elif check_length == "no":
        dictionary[key] = (query, subject, bits, frames, query_starts, subject_gene_name)

def get_highest_scoring(records):
    """ get a dictionary of { species : 'best' hit } from parsed BLASTX output """
    # dictionary of highest cumulative scores for each hit ("alignment")
    highest_scoring_dict = {}
    try:
        for item in records:
            query_len = item.query_letters
            species_dict = defaultdict(list)

            for alignment in item.alignments:
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
                    if scores:
                        total_alignment_score = sum(scores)
                        max_alignment_score = max(scores)
                        # if total score equals max score and is also best total score, take it
                        if total_alignment_score == max_alignment_score and total_alignment_score == max(total_species_scores):
                            add_to_highest_scoring_dict(highest_scoring_dict, species, alignment)
                        # if total score > max score, check if this is the best total score
                        elif total_alignment_score > max_alignment_score and total_alignment_score == max(total_species_scores):
                            # if the query is not too long, keep it
                            add_to_highest_scoring_dict(highest_scoring_dict, species, alignment, check_length="yes", query_length=query_len)
                        # if total score equals max score but is not the best total score, check if it is best max score, then take it
                        # else skip it
                        elif total_alignment_score == max_alignment_score and max_alignment_score == max(max_species_scores):
                            add_to_highest_scoring_dict(highest_scoring_dict, species, alignment)

    except ExpatError:
        print("Unexpected end of xml file...")

    return highest_scoring_dict

def trim_introns_and_seqs_w_stop(highest_scoring_dict):    
    """ given highest scoring hits dictrionary, trim bases 
    that correspond to long gaps in subject (likely introns)
    and delete remaining sequences that still have stop codons """
    trimmed_seq_dict = {}

    for species, seqs in highest_scoring_dict.items():
        (query, subject, bits, frames, query_s, gene) = seqs

        trimm = []
        if find_intron(seqs[1]):
            introns = list(find_intron(seqs[1]))
            if len(introns) == 1:
                intron_start = introns[0][0]
                intron_end = introns[0][1]
                trimm.append(query[:intron_start])
                trimm.append(query[intron_end:])
            elif len(introns) == 2:
                intron1_start = introns[0][0]
                intron1_end = introns[0][1]
                intron2_start = introns[1][0]
                intron2_end = introns[1][1]
                trimm.append(query[:intron1_start])
                trimm.append(query[intron1_end:intron2_start])
                trimm.append(query[intron2_end:])
            elif len(introns) > 2:
                for index, intron in enumerate(introns[:-1]):
                    if index == 0:
                        intron1_start = intron[0]
                        trimm.append(query[:intron1_start]) # first exon
                    if index > 0:
                        current_intron_end = intron[1]
                        next_intron_start = introns[index+1][0]
                        trimm.append(query[current_intron_end:next_intron_start]) # intermediate exons
                last_intron_end = introns[-1][1]
                trimm.append(query[last_intron_end:]) # last exon
        else:
            trimm.append(query)

        trimmed = "".join(trimm)
        
        if "*" not in trimmed:
        # exclude sequences that still have stop codon
            trimmed_no_stop = trimmed
            trimmed_seq_dict[species] = (trimmed_no_stop, query, subject, bits, frames, query_s, gene)


        print('Species: {}\nGene: {}\nQuery: {}\nTrimmed: {}\nSubject: {}\nFrames: {}\nQuery start: {}\n'.format(
            species, gene, query, trimmed, subject, frames, query_s))

    return trimmed_seq_dict

def create_sqlite_db():
    """ create sqlite database """
    conn = sqlite3.connect("test.db")

    with conn:
        cur = conn.cursor()
        cur.execute("DROP TABLE IF EXISTS Filename")        
        cur.execute("DROP TABLE IF EXISTS InputSequence")
        cur.execute("DROP TABLE IF EXISTS Hitname")
        cur.execute("DROP TABLE IF EXISTS Sequences")

        cur.execute("CREATE TABLE {tn} ({fid} INTEGER PRIMARY KEY, {fn} UNIQUE)".format(
         tn="Filename", fid="file_name_ID", fn="file_name"))

        cur.execute("""CREATE TABLE {tn} ({sid} INTEGER PRIMARY KEY, {fid} TEXT, {txn} TEXT,
         FOREIGN KEY({fid}) REFERENCES Filename({fid}))""".format(
         tn="InputSequence", sid="input_seq_ID", fid="file_ID", txn="tax_name"))

        cur.execute("""CREATE TABLE {tn} ({hid} INTEGER PRIMARY KEY, {hn} TEXT,
         FOREIGN KEY({hid}) REFERENCES InputSequence({sid}))""".format(
         tn="Hitname", hid="hit_ID", sid="input_seq_ID", hn="hit_name"))

        cur.execute("""CREATE TABLE {tn} ({sid} INTEGER PRIMARY KEY, {tq} TEXT, {q} TEXT, {s} TEXT,
         FOREIGN KEY({sid}) REFERENCES Hitname({hid}))""".format(
         tn="Sequences", sid="seq_ID", hid="hit_ID", tq="trimmed_query", q="query", s="subject"))

def populate_sqlite_db(file_name, seq_dict):
    """ add data from highest scoring dictionary to sqlite database """
    conn = sqlite3.connect("test.db")

    with conn:
        cur = conn.cursor()
        cur.execute("INSERT INTO Filename(file_name) VALUES (?)", (file_name, ))
        file_ID = cur.lastrowid
        for species, seq in seq_dict.items():
            (trimmed_query, query, subject, bits, frames, query_s, gene) = seq
            cur.execute("INSERT INTO InputSequence(file_ID, tax_name) VALUES (?, ?)", (file_ID, species))
            cur.execute("INSERT INTO Hitname(hit_name) VALUES (?)", (gene, ))
            hit_ID = cur.lastrowid
            cur.execute("INSERT INTO Sequences(seq_ID, trimmed_query, query, subject) VALUES (?, ?, ?, ?)", 
             (hit_ID, trimmed_query, query, subject))

def output_fasta(trimmed_seq_dict, n=80):
    """ produce a FASTA string from dictionary of best hits """
    # each sequence line will have 80 characters 
    # for each element of species : seq dictionary,
    # split sequence into list of string, each n chars long
    # then join everything with newline 
    fasta_string = '\n'.join(['>{}\n{}'.format(species, '\n'.join([seq[0][i:i+n] for i in range(0, len(seq[0]), n)])) \
     for species, seq in sorted(trimmed_seq_dict.items())])
    return fasta_string

def write_fasta(in_file_name, fasta_string):
    """ write FASTA file """
    new_output_name = re.sub(".xml", "", in_file_name)
    out_file_name = "protein-{}".format(new_output_name)
    with open(out_file_name, "w") as f:
        f.write(fasta_string)

def xml_parse(xml_file_name):
    """ parse XML BLASTX output file """
    with open(xml_file_name, "r") as f:
        records = NCBIXML.parse(f)
        return get_highest_scoring(records)

def output_parsed(file_name):
    """ write output of query command """
    new_output_name = re.sub(".xml", "", file_name)
    highest_scoring_dict = xml_parse(file_name)
    trimmed_seq_dict = trim_introns_and_seqs_w_stop(highest_scoring_dict)
    fasta_string = output_fasta(trimmed_seq_dict)
    if fasta_string:
        print("Writing FASTA file protein-{}...".format(new_output_name))
        write_fasta(file_name, fasta_string)
        populate_sqlite_db(file_name, trimmed_seq_dict)
    else:
        print("File {} contains no protein matches".format(file_name))

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
        for file_name in kwargs["xml_files"]:
            output_parsed(file_name)

def run():

    # initialize parsed arguments
    config = ParsedArgs()
    # get arguments
    config_dict = config.get_args_dict()
    return config_dict
    
if __name__ == '__main__':
        
        main()
