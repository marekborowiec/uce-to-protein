#! /usr/bin/env python3

import argparse, re, subprocess, sys
from collections import defaultdict
from multiprocessing.dummy import Pool as dummyPool
from Bio.Blast import NCBIXML
from xml.parsers.expat import ExpatError


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
       # database creation command
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
       # database query command
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
       # parse xml output of BLASTX
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
    # store arguments in a dictionary
        command = self.args.__dict__
        arguments = getattr(self, self.args.command)().__dict__
        argument_dictionary = command.copy()
        argument_dictionary.update(arguments)
        
        return argument_dictionary

class DbCreator():
    """Given arguments, creates a protein database using 'makeblastdb'"""
    def __init__(self, **kwargs):
        self.command = kwargs.get("command")
        

def main():

    def get_highest_scoring(records):
        # list of missing characters in protein sequences
        missing = ["X",".","*","-","?"]
        # dictionary of highest cumulative scores for each hit ("alignment")
        highest_scoring_dict = {}
        try:
            for item in records:
                #print(dir(item))
                query_len = item.query_letters
                #print(query_len)
                species_dict = defaultdict(list)

                for alignment in item.alignments:
                    species_dict[item.query].append(alignment)
                #print(species_dict.items())
                for species, alignments in species_dict.items():
                    #print(species)
                    total_species_scores = []
                    max_species_scores = []
                    for alignment in alignments:
                        scores = [hsp.score for hsp in alignment.hsps]
                        total_alignment_score = sum(scores)
                        max_alignment_score = max(scores)
                        total_species_scores.append(total_alignment_score)
                        max_species_scores.append(max_alignment_score)

                    for alignment in alignments:
                        scores = [hsp.score for hsp in alignment.hsps]
                        total_alignment_score = sum(scores)
                        max_alignment_score = max(scores)
                        query = ''.join([hsp.query for hsp in alignment.hsps])
                        # if total score equals max score and is also best total score, take it
                        if total_alignment_score == max_alignment_score and total_alignment_score == max(total_species_scores):
                            highest_scoring_dict[species] = query
                        # if total score > max score, check if this is the best total score
                        elif total_alignment_score > max_alignment_score and total_alignment_score == max(total_species_scores):
                            # if the query is not too long, keep it
                            if len(query) >= 30 and len(query) <= query_len / 3:
                                highest_scoring_dict[species] = query
                        # if total score equals max score but is not the best total score, check if it is best max score, then take it
                        elif total_alignment_score == max_alignment_score and max_alignment_score == max(max_species_scores):
                            highest_scoring_dict[species] = query
                         # else skip it

                            #print(">{}".format(species))
                        #print('sequence:', alignment.title) 
                        #print('length:', alignment.length)
                                #print(hsp.score)   
                        #print('gaps:', hsp.gaps)
                        #print('e value:', hsp.expect)
                            #print(hsp.query)
                        #print(hsp.match[0:90] + '...')
                        #print(hsp.sbjct[0:90] + '...')
                        #for species, hsp in best_scoring_records.items():
                            #print(species)
                    #print(tpl[1].query)
                    #print(tpl[1].score)
        except ExpatError:
            print("Unexpected end of xml file...")

        return highest_scoring_dict

    def output_fasta(highest_scoring_dict):
        fasta_string = '\n'.join(['>{}\n{}'.format(species, seq) for species, seq in highest_scoring_dict.items()])
        return fasta_string

    def write_fasta(in_file_name, fasta_string):
        new_output_name = re.sub(".xml", "", file_name)
        out_file_name = "protein-{}".format(new_output_name)
        f = open(out_file_name, "w")
        f.write(fasta_string)
        f.close()

    def xml_parser(xml_file):
        f = open(xml_file)
        records = NCBIXML.parse(f)
        return records

    def output_parsed(file_name, records):
        new_output_name = re.sub(".xml", "", file_name)
        f = open(file_name)
        highest_scoring_dict = get_highest_scoring(records)
        fasta_string = output_fasta(highest_scoring_dict)
        if fasta_string:
            print("Writing FASTA file protein-{}...".format(new_output_name))
            write_fasta(file_name, fasta_string)
        else:
            print("File {} contains no protein matches".format(file_name))
        f.close()

    def get_blast_call_string(in_file, db_name):
        call_string = "blastx -query {0} -db {1} -outfmt 5 -evalue 10e-5 > {0}.xml".format(in_file, db_name)
        return call_string

    def call_blast(call_string):
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
        for file_name in kwargs["xml_files"]:
            records = xml_parser(file_name)
            output_parsed(file_name, records)

def run():

    # initialize parsed arguments
    config = ParsedArgs()
    # get arguments
    config_dict = config.get_args_dict()
    return config_dict
    
if __name__ == '__main__':
        
        main()
