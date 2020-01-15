from Bio import SeqIO, AlignIO
from Bio import Entrez
import pandas as pd
import random
import os
import pickle
import ete3
import numpy

Entrez.email = "gabriel.foley@uqconnect.edu.au"


def load_sequences(*args, split_char=""):
    """
    Join multiple sequence files together
    :param args: The sequence files to open
    :param split_char: A character which comes after the sequence ID
    :return:
    """
    full_dict = {}
    for file in args:
        current_handle = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
        full_dict.update(current_handle)

    for name in full_dict.keys():

        if name.startswith(("XP", "XM", "XR", "WP", "NP", "NC", "NG", "NM", "NR")):
            full_dict[name].annotations["Database"] = "NCBI"

            # Split on the given symbol if it is in the header
            if split_char:
                if split_char in name:
                    full_dict[name].id = name.split(split_char)[0]

        elif name.startswith(("sp", "tr", "gi")):
            full_dict[name].id = name.split("|")[1]
            full_dict[name].annotations["Database"] = "UniProt"

        elif name.startswith("pdb"):
            full_dict[name].annotations["Database"] = "PDB"

        else:
            full_dict[name].annotations["Database"] = "Unknown"

            # Split on the given symbol if it is in the header
            if split_char:
                if split_char in name:
                    full_dict[name].id = name.split(split_char)[0]

    return full_dict


def load_alignment(filepath, file_type="fasta"):
    handle = list(AlignIO.parse(filepath, file_type))
    alignment = handle[0]
    return alignment

def write_alignment(alignment, filepath, file_type="fasta"):
    AlignIO.write(alignment, filepath, file_type)

def load_tree(filepath, tree_format=1):
    tree = ete3.Tree(filepath, format=tree_format)
    return tree


def save_ids(*args, percent_identity="", output_dir="", concatenate=False):
    """
    Load in a list of IDs and retreive the sequences for them
    :param args: The lists to open
    :param percent_identity: The minimum percent identity to filter by
    :param output_dir: Optional filepath to write file to
    :param concatenate: Whether we should concatenate the lists and remove duplicates
    :return: SeqRecord object if no file path specified, otherwise nothing
    """

    id_list = []
    output_id = random_string(5)

    for file in args:
        outpath = file[file.rindex('/') + 1:].replace(".", "_output.")

        # df = pd.read_csv(file, delimiter='\t', header=None)
        df = pd.read_csv(file, delimiter='\t', comment='#', names=["query acc.ver", "subject acc.ver", "% identity",
                                                                   "alignment length", "mismatches", "gap opens",
                                                                   "q. start", "q. end", "s. start", "s. end",
                                                                   "evalue", "bit score"])

        # If user has set a minimum percent identity then filter by that
        if percent_identity:
            id_list += df.loc[df['% identity'] >= percent_identity, 'subject acc.ver'].drop_duplicates().tolist()

        # No minimum percent identity so take all records
        else:
            id_list = df['subject acc.ver'].drop_duplicates().tolist()

        # If we're not concatenating lets write this out to the
        if not concatenate:
            # Remove any existing previous file at this location
            # remove_file(output_dir + outpath)

            # Write out the sequences to the output path
            get_ids(id_list, output_dir + outpath)
            id_list = []

    if concatenate:
        # Get rid of any duplicates in id_list (i.e. hits that appeared in more than one input file)
        id_list = list(set(id_list))
        get_ids(id_list, output_dir + output_id + "_concatenated_output.fasta")


def load_ids(*args, percent_identity=""):
    id_list = []

    for file in args:
        # df = pd.read_csv(file, delimiter='\t', header=None)
        df = pd.read_csv(file, delimiter='\t', comment='#', names=["query acc.ver", "subject acc.ver", "% identity",
                                                                   "alignment length", "mismatches", "gap opens",
                                                                   "q. start", "q. end", "s. start", "s. end",
                                                                   "evalue", "bit score"])

        # If user has set a minimum percent identity then filter by that
        if percent_identity:
            id_list += df.loc[df['% identity'] >= percent_identity, 'subject acc.ver'].drop_duplicates().tolist()

        # No minimum percent identity so take all records
        else:
            id_list = df['subject acc.ver'].drop_duplicates().tolist()

    # Remove any duplicates in id_list (i.e. hits that appeared in more than one input file)
    id_list = list(set(id_list))
    handle = get_ids(id_list)
    return handle


def get_ids(id_list, filepath=""):
    for i in range(0, len(id_list), 500):

        final = i + 500 if (i + 500 < len(id_list)) else len(id_list) + 1
        handle = Entrez.efetch(db="protein", rettype="fasta", retmode="text", id=id_list[i:final])
        count = 0
        if filepath:
            os.chdir(filepath.rsplit('/', 1)[0] + "/")
            with open(filepath, "a") as outfile:
                for seq_record in SeqIO.parse(handle, "fasta"):
                    count += 1
                    outfile.write(">" + seq_record.description + "\n")
                    outfile.write(str(seq_record.seq + "\n"))
                handle.close()

        else:
            return SeqIO.parse(handle, "fasta")


def save_header_terms(header_terms, filepath):
    header_string = " ".join(header_terms)
    print(header_string)
    print(type(header_string))
    with open(filepath, "w") as text_file:
        text_file.write(header_string)


def load_header_terms(filepath):
    header_terms = []
    with open(filepath, "r") as text_file:
        for item in text_file.read().split():
            header_terms.append(item)
    return header_terms


def add_header_terms(sender, header_terms):
    for item in sender.value.split():
        header_terms.append(item)


def build_taxonomy_dict(seq_ids, seq_type="protein"):
    """
    Take a list of sequence ids and return a dictionary mapping those ids to their taxonomic ID

    :param seq_ids: List of sequence ids
    :param seq_type: Type of sequence / which database to query
    :return: Dictionary mapping seq ID to taxonomic ID
    """

    taxonomy_dict = {}

    for seq_id in seq_ids:
        print (seq_id)
        handle = Entrez.elink(dbfrom=seq_type, db="taxonomy", id=seq_id)
        records = Entrez.read(handle)
        if len(records[0]["LinkSetDb"]) > 0:
            taxonomy_dict[seq_id] = records[0]["LinkSetDb"][0]["Link"][0]['Id']
        else:
            taxonomy_dict[seq_id] = 0000

    return taxonomy_dict


def random_string(length=10):
    """
    Create a random string
    :param length: Length of string
    :return: Random string
    """
    valid_letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    return ''.join((random.choice(valid_letters) for i in range(length)))


def remove_file(*args):
    """
    Remove files in the list from the directory

    :param args: Files to remove
    :return:
    """
    for arg in args:
        print(arg)
        if os.path.exists(arg):
            print("exists")
            os.remove(arg)


def save_python_object(python_object, filepath):
    """
    Save a python object into the given filepath
    :param python_object: The object to save
    :param filepath: The filepath to write to
    :return:
    """
    with open (filepath, 'wb') as pickle_file:
        pickle.dump(python_object, pickle_file)


def open_python_object(filepath):
    """
    Open a python obbject from the given filepath
    :param filepath: The filepath to open from
    :return: The object
    """
    with open(filepath, 'rb') as pickle_file:
        python_object = pickle.load(pickle_file)

    return python_object

def get_mean(records):
    """
    Take a list or dictionary of values and return the mean
    :param records:
    :return:
    """
    if type(records) == dict:
        print ('dict')
        records_array = numpy.array([value for value in records.values()])
        print (records_array)
        print (records_array.mean())


def get_standard_deviation(records):
    """ Take a list or dictionary of values and return the standard deviation"""
    if type(records) == dict:
        records_array = numpy.array([value for value in records.values()])
        print (records_array)
        print (records_array.std())
