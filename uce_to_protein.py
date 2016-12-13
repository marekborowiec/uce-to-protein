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

__version__ = "0.1"

import argparse
import configparser as cfgp
import re
import sqlite3
import subprocess
import os.path
import sys
from collections import defaultdict
from multiprocessing.dummy import Pool as dummyPool
from operator import itemgetter
from xml.parsers.expat import ExpatError

from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class ParsedArgs:

    def __init__(self):
        parser = argparse.ArgumentParser(
            usage="""uce_to_protein <command> [<args>]
The commands are:
  blastdb        Create BLAST database
  queryblast     Query existing BLAST database and output BLASTX xml
  parse          Parse BLASTX xml output and create or update database of best hits
  queryprot      Query existing database of best hits
Use uce_to_protein <command> -h for help with arguments of the command of interest
""")

        parser.add_argument(
         "command",
         help="Subcommand to run")

        # parse_args defaults to [1:] for args, but you need to
        # exclude the rest of the args too, or validation will fail
        self.args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, self.args.command):
            print('Unrecognized command')
            parser.print_help()
            exit(1)
        # use dispatch pattern to invoke method with same name
        getattr(self, self.args.command)()

    def blastdb(self):
        """ command that wraps BLAST for creating protein database """
        parser = argparse.ArgumentParser(
         description="Create a BLAST database", )
        parser.add_argument(
         "-i",
         "--input",
         dest = "prot_input",
         help = "Input filename with reference protein sequences")
        parser.add_argument(
         "-o",
         "--output",
         dest = "output_db_name",
         default = "prot_db",
         help = "File name for BLAST database created from reference sequences. "
          "Default: 'prot_db'")
        parser.add_argument(
         "-a",
         "--other-args",
         dest = "other_args",
         default = "",
         type=str,
         help = "A string of any additional commands to MAKEBLASTDB. "
          "Must be quoted. Example: --other-args '-logfile blast-db.log'")

        args = parser.parse_args(sys.argv[2:])
        return args

    def queryblast(self):
        """ command that wraps BLASTX to query UCE alignment against protein database """
        parser = argparse.ArgumentParser(
         description="Query existing protein database with unaligned UCE sequences", )
        parser.add_argument(
         "-i",
         "--input-seqs",
         nargs = "+",
         type = str,
         required = True,
         dest = "query",
         help = "Input FASTA files with nucleotide sequences "
          "to be compared against reference")
        parser.add_argument(
         "-d",
         "--database",
         dest = "input_db_name",
         default = "prot_db",
         help = "File name for BLAST database created from reference sequences. "
          "Default: 'prot_db'")
        parser.add_argument(
         "-e",
         "--e-value",
         dest = "e_value",
         default = "10e-5",
         help = "E-value cutoff for BLAST hits to be retained. Default: '10e-5'")
        parser.add_argument(
         "-c",
         "--cores",
         dest = "cores",
         default = 1,
         help = "Number of cores used for parallel BLASTX searches. Default: 1")
        parser.add_argument(
         "-f",
         "--out-config",
         dest = "out_config",
         default = "fasta_to_xml.conf",
         help = "Name of output config file used for parsing. "
          "Default: 'fasta_to_xml.conf'")
        parser.add_argument(
         "-a",
         "--other-args",
         dest = "other_args",
         type=str,
         default = "",
         help = "A string of any additional commands to BLASTX. Must be quoted. "
          "Example: --other-args '-seg no'")

        args = parser.parse_args(sys.argv[2:])
        return args

    def parse(self):
        """ command to parse BLASTX output produced by query command """
        parser = argparse.ArgumentParser(
         description="Parse XML output of BLASTX", )
        parser.add_argument(
         "-i",
         "--input-config",
         type = str,
         required = True,
         dest = "fasta_to_xml_config",
         help = "Name of configuration file that lists corresponding nucleotide "
          "FASTA and BLASTX XML files")
        parser.add_argument(
         "-o",
         "--output",
         dest = "best_hits",
         default = "best_hits.sqlite",
         help = "File name for best hits database created from BLASTX XML files. "
          "Default: 'best_hits.sqlite'")

        args = parser.parse_args(sys.argv[2:])
        return args

    def queryprot(self):
        """ command to extract sequences from database of best hits """
        parser = argparse.ArgumentParser(
         description="Get sequences from best hits database", )
        parser.add_argument(
         "-d",
         "--database",
         type = str,
         required = True,
         dest = "best_hits_db",
         help = "Name of sqlite database containing best hits")
        parser.add_argument(
         "-c",
         "--config",
         dest = "taxon_config",
         help = "Name of configuration file containing all taxa to be queried")
        parser.add_argument(
         "-g",
         "--taxon-group",
         dest = "taxon_group",
         help = "Name of group in configuration file with taxa to be extracted")

        args = parser.parse_args(sys.argv[2:])
        return args

    def get_args_dict(self):
        """ store arguments in a dictionary """
        command = self.args.__dict__
        arguments = getattr(self, self.args.command)().__dict__
        argument_dictionary = command.copy()
        argument_dictionary.update(arguments)
        return argument_dictionary


class ArgCreator():
    """ arguments class """
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


def get_frames(alignment):
    """ get frames from all matches in alignment """
    return [hsp.frame for hsp in alignment.hsps]


def get_nt_extracted(sign, nt_string, increasing_start_end):
    """ match and extract sequence from nucleotide string using best hit coordinates """
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


def check_overlap(tpl):
    """ check if numbers overlap in tuple """
    for index, pos in enumerate(tpl):
        if index > 0:
            if pos <= tpl[index - 1]:
                return True


def queries_overlap(alignment):
    """ check if ranges of hsps queries overlap """
    starts = [int(hsp.query_start) for hsp in alignment.hsps]
    ends = [int(hsp.query_end) for hsp in alignment.hsps]
    start_end = zip(starts, ends)
    increasing_start_end = sorted(start_end, key=itemgetter(0))
    flat_incr_start_end = tuple(sum(increasing_start_end, ()))

    if check_overlap(flat_incr_start_end):
        return True
    else:
        return False


def get_sorted_queries(sign, nucleotides, alignment):
    """ get hsps queries and subject sorted by query start,
    either ascending or descenging, depending on strand """
    starts = [int(hsp.query_start) for hsp in alignment.hsps]
    ends = [int(hsp.query_end) for hsp in alignment.hsps]
    queries = [hsp.query for hsp in alignment.hsps]
    subjects = [hsp.sbjct for hsp in alignment.hsps]
    start_end = zip(starts, ends)
    increasing_start_end = sorted(start_end, key=itemgetter(0))   # this is needed to check for overlap and extract from nt query
    zipped = zip(starts, ends, queries, subjects)                 # this is needed to sort all according to strand

    if sign == "plus":
        sorted_zipped = sorted(zipped, key=itemgetter(0))
    elif sign == "minus":
        sorted_zipped = sorted(zipped, key=itemgetter(0), reverse=True)

    (sorted_starts, sorted_ends, sorted_queries, sorted_subjects) = zip(*sorted_zipped)

    nt_query_string = get_nt_extracted(sign, nucleotides, increasing_start_end)
    query_string = ''.join(sorted_queries)
    subject_string = ''.join(sorted_subjects)

    return (nt_query_string, query_string, subject_string, sorted_starts, sorted_ends)


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
                    trimm_nt.append(nt_query[:(intron1_start * 3)])                               # first exon
                    trimm_aa.append(query[:intron1_start])

                if index > 0:
                    current_intron_start = intron[0]
                    current_intron_end = intron[1]
                    next_intron_start = introns[index+1][0]
                    trimm_nt.append(nt_query[(current_intron_end * 3):(next_intron_start * 3)])   # intermediate exons
                    trimm_aa.append(query[current_intron_end:next_intron_start])

            last_intron_start = introns[-1][1]
            last_intron_end = introns[-1][1]
            trimm_nt.append(nt_query[(last_intron_end * 3):])                                     # last exon
            trimm_aa.append(query[last_intron_end:])

    else:
        trimm_nt.append(nt_query)
        trimm_aa.append(query)

    trimmed_nt = "".join(trimm_nt)
    trimmed_aa = "".join(trimm_aa)

    if "*" not in trimmed_aa:    # exclude sequences that still have stop codon
        return (trimmed_nt, trimmed_aa)


def add_to_highest_scoring_dict(highest_score_dict, nucleotides, key, alignment):
    """ assign items to dictionary of highest scoring hits """

    if queries_overlap(alignment) is False:    # do not add if ranges overlap
        frames = get_frames(alignment)
        if frames[0][0] > 0:                   # checking sign of first frame
            sign = "plus"
        elif frames[0][0] < 0:
            sign = "minus"

        (nt_query, query, subject, sorted_starts, sorted_ends) = get_sorted_queries(sign, nucleotides, alignment)
        trimmed_queries = trim_introns_and_seqs_w_stop(nt_query, query, subject)
        subject_gene_name = alignment.title

        try:
            (trimmed_nt_query, trimmed_query) = trimmed_queries
        except TypeError:
            trimmed_query = False

        if trimmed_query:    # in case trimmed query was excluded because it still had stop codons
                highest_score_dict[key] = (
                 trimmed_nt_query, trimmed_query, query, subject,
                 frames, sorted_starts, sorted_ends, subject_gene_name)


def get_highest_scoring(records, nucleotide):
    """ get a dictionary of { species : 'best' hit } from parsed BLASTX output """
    nt_dict = {record.id : str(record.seq) for record in nucleotide}    # { species : nt sequence } made from parsed FASTA file
    highest_scoring_dict = {}
    try:
        for item in records:

            species_dict = defaultdict(list)

            for alignment in item.alignments:                           # alignment here refers to BLAST hit
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
                    try:
                        nt = nt_dict[species]
                    except KeyError:
                        print("ERROR: It appears taxa names do not match between your FASTA and XML files")
                        sys.exit()
                    if scores:
                        total_alignment_score = sum(scores)
                        max_alignment_score = max(scores)
                        # if total score equals max score and is also best total score, take it
                        if total_alignment_score == max_alignment_score \
                                and total_alignment_score == max(total_species_scores):
                            add_to_highest_scoring_dict(highest_scoring_dict, nt, species, alignment)

                        # if total score > max score, check if this is the best total score
                        # at this step alignment will be rejected if ranges overlap
                        elif total_alignment_score > max_alignment_score \
                                and total_alignment_score == max(total_species_scores):
                            add_to_highest_scoring_dict(highest_scoring_dict, nt, species, alignment)

                        # if total score equals max score but is not the best total score,
                        # check if it is best max score, then take it
                        # else skip it
                        elif total_alignment_score == max_alignment_score \
                                and max_alignment_score == max(max_species_scores):
                            add_to_highest_scoring_dict(highest_scoring_dict, nt, species, alignment)

    except ExpatError:
        print("Unexpected end of xml file...")
        sys.stdout.flush()

    return highest_scoring_dict


def create_sqlite_db(db_name):
    """ create sqlite database """
    conn = sqlite3.connect(db_name)

    with conn:
        cur = conn.cursor()

        cur.execute(
         "CREATE TABLE {tn} ({uid} INTEGER PRIMARY KEY, {un} UNIQUE)".format(
         tn="Uces", uid="uce_name_ID", un="uce_name"))

        cur.execute(
         "CREATE TABLE {tn} ({tid} INTEGER PRIMARY KEY, {txn} UNIQUE)".format(
         tn="Taxa", tid="taxon_name_ID", txn="taxon_name"))

        cur.execute(
         "CREATE TABLE {tn} ({tid} INTEGER PRIMARY KEY, {txn} UNIQUE)".format(
         tn="Hitnames", tid="hit_name_ID", txn="hit_name"))

        cur.execute(
         """CREATE TABLE {tn} 
         ({sid} INTEGER PRIMARY KEY, {uid} INT, {tid} INT, {hid} INT, {tnq} TEXT, {tpq} TEXT, {upq} TEXT, {s} TEXT,
         FOREIGN KEY({uid}) REFERENCES Uces({uid}),
         FOREIGN KEY({tid}) REFERENCES Taxa({tid}),
         FOREIGN KEY({hid}) REFERENCES Hitnames({hid}))
         """.format(
         tn="Sequences", sid="seq_ID", uid="uce_name_ID",
         tid="taxon_name_ID", hid="hit_name_ID",
         tnq="trimmed_nuc_query", tpq="trimmed_prot_query", upq="untrimmed_prot_query", s="subject"))


def add_to_sqlite_db(uce_name, seq_dict, db_name):
    """ add data from highest scoring dictionary to sqlite database """
    conn = sqlite3.connect(db_name)

    with conn:
        cur = conn.cursor()
        cur.execute("INSERT OR IGNORE INTO Uces(uce_name) VALUES (?)", (uce_name, ))
        cur.execute("SELECT uce_name_ID FROM Uces WHERE uce_name = ?", (uce_name, ))
        uce_ID = int(cur.fetchone()[0])

        for species, seq in seq_dict.items():
            (trimmed_nuc_query, trimmed_prot_query, untrimmed_prot_query, subject, frames, query_s, query_e, gene) = seq

            cur.execute("INSERT OR IGNORE INTO Taxa(taxon_name) VALUES (?)", (species, ))
            cur.execute("INSERT OR IGNORE INTO Hitnames(hit_name) VALUES (?)", (gene, ))

            cur.execute("SELECT taxon_name_ID FROM Taxa WHERE taxon_name = ?", (species, ))
            taxon_ID = int(cur.fetchone()[0])
            cur.execute("SELECT hit_name_ID FROM Hitnames WHERE hit_name = ?", (gene, ))
            hit_ID = int(cur.fetchone()[0])

            # this insert statement ensures that no duplicate records are added when updating existing database
            # duplicates are defined as entries that have both uce name ID and taxon name ID that already exist in database
            cur.execute(
             """INSERT INTO 
             Sequences(uce_name_ID, taxon_name_ID, hit_name_ID, trimmed_nuc_query, trimmed_prot_query, untrimmed_prot_query, subject)
             SELECT ?, ?, ?, ?, ?, ?, ?
             WHERE NOT EXISTS(SELECT 1 FROM Sequences WHERE uce_name_ID = ? AND taxon_name_ID = ?); """, 
             (uce_ID, taxon_ID, hit_ID, trimmed_nuc_query, trimmed_prot_query, untrimmed_prot_query, subject, uce_ID, taxon_ID))


def input_parse(xml_file_name, nt_file_name):
    """ parse XML BLASTX output file """
    with open(xml_file_name, "r") as f, open(nt_file_name, "r") as nt:
        records = NCBIXML.parse(f)
        nucleotide = SeqIO.parse(nt, "fasta")
        return get_highest_scoring(records, nucleotide)


def populate_sqlite_db(fasta_file_name, blast_file_name, db_name):
    """ write output of query command """
    base_uce_name = fasta_file_name.split(".")[0]
    base_uce_name_xml = blast_file_name.split(".")[0]
    if base_uce_name != base_uce_name_xml:
        print("ERROR: Locus names of your files '{}' and '{}' do no match.\n"
         "Are you sure you are comparing the right files?".format(
         blast_file_name, fasta_file_name))
        sys.stdout.flush()
        sys.exit()
    highest_scoring_dict = input_parse(blast_file_name, fasta_file_name)
    if highest_scoring_dict:
        print("Got best hit for {} ...".format(base_uce_name))
        sys.stdout.flush()
        add_to_sqlite_db(base_uce_name, highest_scoring_dict, db_name)
    else:
        print("Locus {} contains no protein matches".format(base_uce_name))
        sys.stdout.flush()


def get_blast_call_string(in_file, db_name, e_value, other_args):
    """ make command line call string for BLASTX """
    call_string = "blastx -query {0} -db {1} -outfmt 5 -evalue {2} {3} > {0}.xml".format(
     in_file, db_name, e_value, other_args)
    return call_string


def call_blast(call_string):
    """ make command line call for BLASTX """
    a = re.search("-query (\S+)", call_string)
    b = re.search("-db (\S+)", call_string)
    file_name = a.group(1)
    database = b.group(1)
    print("Blasting file '{}' against protein database '{}' ...".format(file_name, database))
    sys.stdout.flush()
    try:
        p = subprocess.call(call_string, shell=True)
    except KeyboardInterrupt:
        print("\nYou killed the query")
        sys.stdout.flush()
        sys.exit()


def query_taxa_from_sqlite(db_name, taxa):
    """ query sqlite database for all records and extract those in config """
    conn = sqlite3.connect(db_name)

    with conn:
        cur = conn.cursor()
        # get data for all taxa
        cur.execute(
         """SELECT uce_name, taxon_name, trimmed_nuc_query, trimmed_prot_query FROM Sequences 
         INNER JOIN Uces ON Sequences.uce_name_ID = Uces.uce_name_ID 
         INNER JOIN Taxa ON Sequences.taxon_name_ID = Taxa.taxon_name_ID""")
        rows = cur.fetchall()
        group_rows = [row for row in rows if row[1] in taxa]
        sorted_group_rows = sorted(group_rows, key=itemgetter(1))
        uces_dict = defaultdict(list)

        for row in sorted_group_rows:
            row_uce, taxon, nt_seq, aa_seq = row
            uces_dict[row_uce].append((taxon, nt_seq, aa_seq))

        return uces_dict


def write_fasta_to_xml_config(**kwargs):
    """ write config file with fasta and xml file name pairs """
    config = cfgp.ConfigParser()
    fasta_to_xml_dict = {file : "{}.xml".format(file) for file in kwargs["query"]}
    config["fasta_to_xml"] = fasta_to_xml_dict
    with open(kwargs["out_config"], "w") as outconfig:
        config.write(outconfig)
    print("Wrote config file {}".format(kwargs["out_config"]))
    sys.stdout.flush()


def parse_fasta_to_xml_config(config_name):
    """ parse fasta, xml file name pairs from config """
    config = cfgp.ConfigParser()
    config.optionxform = lambda option: option
    config.read(config_name)
    file_pairs = sorted(dict(config["fasta_to_xml"]).items())
    for pair in file_pairs:
        (fasta, xml) = pair
    return file_pairs


def parse_taxon_config(config_name, group):
    """ parse taxon group from taxon config file """
    config = cfgp.ConfigParser(allow_no_value=True)
    config.optionxform = lambda option: option
    config.read(config_name)
    try:
        return [taxon for taxon in config[group]]
    except KeyError:
        print("ERROR: No group named {} in {}".format(group, config_name))
        sys.stdout.flush()
        sys.exit()


def drop_3rd_codon_pos(seq):
    """ drop every third character from sequence """
    return "".join([char for index, char in enumerate(seq) if (index + 1) % 3 != 0])


def write_fasta(uces_dict, group):
    """ write fasta output given dictionary { uce : (taxon, nuc, aa) } """
    for uce, record in sorted(uces_dict.items()):
        nt_file_name = "{}-coding-nuc-{}.unaligned.fasta".format(group, uce)
        no3rd_file_name = "{}-coding-no3rd-{}.unaligned.fasta".format(group, uce)
        prot_file_name = "{}-protein-{}.unaligned.fasta".format(group, uce)

        bio_nt_seqs = []
        bio_no3rd_seqs = []
        bio_prot_seqs = []

        for species in record:
            (taxon_name, nt_seq, aa_seq) = species

            nuc_seq = SeqRecord(Seq(nt_seq), id=taxon_name, description="")
            bio_nt_seqs.append(nuc_seq)

            no3rd_seq = SeqRecord(Seq(drop_3rd_codon_pos(nt_seq)), id=taxon_name, description="")
            bio_no3rd_seqs.append(no3rd_seq)

            prot_seq = SeqRecord(Seq(aa_seq), id=taxon_name, description="")
            bio_prot_seqs.append(prot_seq)

        print("Writing coding nucleotide fasta file '{}'".format(nt_file_name))
        sys.stdout.flush()
        SeqIO.write(bio_nt_seqs, nt_file_name, "fasta")
        print("Writing 3rd codon positions removed fasta file '{}'".format(no3rd_file_name))
        sys.stdout.flush()
        SeqIO.write(bio_no3rd_seqs, no3rd_file_name, "fasta")
        print("Writing protein fasta file '{}'".format(prot_file_name))
        sys.stdout.flush()
        SeqIO.write(bio_prot_seqs, prot_file_name, "fasta")


def main():

    kwargs = run()
    arg_creator = ArgCreator(**kwargs)  # initialize parsed arguments and argument creator object

    if arg_creator.command == "blastdb":
        call_string = "makeblastdb -in {} -dbtype prot -out {} {}".format(
         kwargs["prot_input"], kwargs["output_db_name"], kwargs["other_args"])
        subprocess.call(call_string, shell=True)

    if arg_creator.command == "queryblast":
        if int(kwargs["cores"]) == 1:
            for file in kwargs["query"]:
                blast_command = get_blast_call_string(
                 file, kwargs["input_db_name"], kwargs["e_value"], kwargs["other_args"])
                call_blast(blast_command)

        if int(kwargs["cores"]) > 1:
            commands = [get_blast_call_string(
             file, kwargs["input_db_name"], kwargs["e_value"],
             kwargs["other_args"]) for file in kwargs["query"]]
            pool = dummyPool(int(kwargs["cores"]))
            try:
                pool.map(call_blast, commands)
            except KeyboardInterrupt:
                pool.close()
                pool.terminate()
                print("\nYou killed the query")
                sys.stdout.flush()
                sys.exit()

        write_fasta_to_xml_config(**kwargs)   # write config file

    if arg_creator.command == "parse":
        db_name = kwargs["best_hits"]
        fasta_xml_pairs = parse_fasta_to_xml_config(kwargs["fasta_to_xml_config"])

        if not os.path.isfile(db_name) or os.stat(db_name).st_size == 0:
            create_sqlite_db(db_name)

        for pair in fasta_xml_pairs:
            (fasta_fn, blast_fn) = pair
            populate_sqlite_db(fasta_fn, blast_fn, db_name)

        print("Wrote results to {} database".format(db_name))
        sys.stdout.flush()

    if arg_creator.command == "queryprot":
        db_name = kwargs["best_hits_db"]
        config_name = kwargs["taxon_config"]
        group = kwargs["taxon_group"]
        if os.path.exists(config_name):
            taxa = parse_taxon_config(config_name, group)
            uces_dict = query_taxa_from_sqlite(db_name, taxa)
            if uces_dict:
                write_fasta(uces_dict, group)
            else:
                print("ERROR: Taxa from '{}' do not match any in your '{}' database.".format(
                 group, db_name)) 
                sys.stdout.flush()
                sys.exit()
        else:
            print("ERROR: File '{}' does not exist.".format(config_name))
            sys.stdout.flush()
            sys.exit()


def run():

    config = ParsedArgs()                   # initialize parsed arguments
    config_dict = config.get_args_dict()    # get arguments
    return config_dict

if __name__ == '__main__':

        main()
