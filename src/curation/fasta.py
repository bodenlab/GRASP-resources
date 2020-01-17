from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import plotly
import utilities
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from plotly.graph_objs import *
# plotly.tools.set_credentials_file(username='gabefoley', api_key='xS8qT0kIbIKDWt0BalOd')
import re

init_notebook_mode(connected=True)

def print_record_overview(records):
    """
    Print the total number of sequences and the average length of sequences
    :param records: The records to print information about
    """
    total_len = 0
    for record in records.values():
        total_len += len(record)
    print("There are %s sequences and the average length of sequence is %d" % (len(records), total_len / len(records)))


def build_species_count(*exclude_on, records, length=0):
    """
    Builds up a dictionary with species as keys and lists of record ids that are from that species as values
    :param exclude_on:
    :param records:
    :param length:
    :return:
    """
    species_count = defaultdict(list)
    for record in records.values():
        if len(str(record.seq)) > length and not any(exclude in record.description for exclude in exclude_on):
            if "[" in record.description:
                name = record.description.split(">")[0].split("[")[1].split("]")[0]
            else:
                name = "_".join(record.description.split("|")[1].split("_")[0:2])
            species_count[name].append(record.id)
    return species_count


def subset_records(*header_terms, records, length=0, mode="exclude", ignore_case=True):
    """
    Return a subset of records, by defining terms that should or shouldn't be in the header and a minimum length

    :param header_terms: The terms to look for in the header
    :param records: The records to make a subset from
    :param length: The minimum length of sequences to keep
    :param mode: Whether to exclude or include on the presence of terms in the header
    :param ignore_case: Whether we should ignore case of the header terms we're checking
    :return:
    """
    new_records = {}

    for record in records.values():

        if mode == "exclude":
            if ignore_case and len(str(record.seq)) > length:
                header_lower = [i.lower() for i in header_terms]
                record_lower = record.description
                record_lower = record_lower.lower()
                if not any(term in record_lower for term in header_lower):
                    new_records[record.name] = record

            elif len(str(record.seq)) > length and not any(term in record.description for term in header_terms):
                    new_records[record.name] = record

        elif mode == "include":
            if ignore_case and len(str(record.seq)) > length:
                header_lower = [i.lower() for i in header_terms]
                record_lower = record.description
                record_lower = record_lower.lower()
                if any(term in record_lower for term in header_lower):
                    new_records[record.name] = record
            elif len(str(record.seq)) > length and any(term in record.description for term in header_terms):
                new_records[record.name] = record

    return new_records


def add_species_from_header(records):
    """ Add species names that exist within the header of the seq ID i.e. like [Amazona aestiva] to a dictionary of
    seqRecords
    :param records: The dictionary to update with species names
    :return:
    """

    for record in records:
        if "[" in records[record].description:
            species_name = records[record].description.split(">")[0].split("[")[1].split("]")[0]
            records[record].annotations["Species"] = species_name
    return records


def add_species_from_last_position(records):
    """ Add species names that exist within the last positiong of the seq ID to a dictionary of seqRecords
    :param records: The dictionary to update with species names
    :return:
    """

    for record in records:
        if "_" in records[record].description:
            species_name = records[record].description.split("_")[-1]
            records[record].annotations["Species"] = species_name
    return records


def get_records_without_species(records):
    """
    Get the entries in a seq records file that don't have a species annotated with them
    :param records: The records to check
    :return: The list of seq IDs missing species information
    """
    missing = []

    for record in records.values():
        if "Species" not in record.annotations:
            missing.append(record.id)

    return missing


def report_records_without_species(records):
    """
    Print out statements about whether entries in a seq record file have annotated species
    :param records:
    :return:
    """
    missing = get_records_without_species(records)
    if missing:
        print("The following records do not have species annotated", [x for x in missing])
    else:
        print("All of the records have a species annotated")


def check_cd_hit_output(records, cdhit_output):
    """
    Reports on the species information that occurs within clusters in CD-hit output
    :param records:
    :param cdhit_output:
    :return:
    """

    if get_records_without_species(records):
        records = add_species_from_header(records)
        if get_records_without_species(records):
            print("Some of your entries are missing their species annotation")
            return

    clusters = defaultdict(list)
    current = None

    # Read in and organise the CD-HIT output

    for line in cdhit_output:
        if "Cluster" in line:
            current = line

        if "Cluster" not in line:

            clusters[current].append(re.sub(r"\s+", " ", line))

    # Only interested in clusters with more than one sequence
    for cluster in clusters.values():
        if len(cluster) > 1:
            seq_list = []
            for seq in cluster:
                trim = re.search('>(.*)\.\.\.', seq)
                seq_id = trim.group(1).strip()
                if "*" in seq:
                    seq_info = "Retained " + seq_id + " " + records[seq_id].annotations["Species"]
                    seq_list.insert(0, seq_info)  # Add the retained sequence to the start of the seq list
                else:
                    seq_info = "Removed " + seq_id + " " + records[seq_id].annotations["Species"]
                    seq_list.append(seq_info)
            print("")
            for seq in seq_list:
                print(seq)


def subset_records_with_regex(*header_terms, records, length=0, mode="exclude"):
    search_term = ""
    for term in header_terms:
        search_term += "\s" + term + "\d*-like|" + term + "\d*|"

    regex = re.compile(search_term[0:len(search_term) - 1])

    new_records = {}
    for record in records.values():
        if mode == "exclude":
            if len(str(record.seq)) > length and not regex.search(record.description):
                new_records[record.name] = record
        elif mode == "include":
            if len(str(record.seq)) > length and regex.search(record.description):
                new_records[record.name] = record

    return new_records


def exclude_character(records, character, mode="exclude"):
    new_records = {}
    for record in records.values():
        if mode == "exclude":
            if character not in record.seq:
                new_records[record.name] = record
        elif mode == "include":
            if character in record.seq:
                new_records[record.name] = record
    return new_records


def build_species_count(*header_terms, records, length=0, mode="exclude"):
    species_count = defaultdict(list)
    for record in records.values():
        if mode == "exclude":
            if len(str(record.seq)) > length and not any(term in record.description for term in header_terms):
                if "[" in record.description:
                    name = record.description.split("[")[1].split("]")[0]
            else:
                print(record.description)
            species_count[name].append(record.id)
    return species_count


def count_ids(records, min_length=0):
    count = 0
    for record in records.values():
        if len(record) >= min_length:
            count += len(record)
    return count


def get_species_names(records, min_length=0, counts=False):
    species_names = []
    for k,v in records.items():
        if len(v) >= min_length:
            species_names.append(k) if not counts else species_names.append(k + " " + str(len(v)))
    return species_names


def plot_record_number(records, plot_type, min_length=0):
    plot_records = {}
    # If we need to restrict the records to plot based on a minimum number
    if min_length > 0:

        for k, v in records.items():
            if len(v) >= min_length:
                plot_records[k] = v
    else:
        plot_records = records

    if plot_type == "Bar":

        data = [Bar(x=list(plot_records.keys()),
                    y=[len(x) for x in plot_records.values()])]

    elif plot_type == "Pie":
        trace = go.Pie(labels=plot_records.keys(), values=plot_records.values())

        py.iplot([trace], filename='basic_pie_chart')

    else:
        print("Invalid plot type selected. Choose 'Bar' or 'Pie'")

    return data


def map_dict_to_records(records, full_dict={}, unique=False):
    # If full_dict hasn't been provided, take records as the full dictionary
    if not full_dict:
        full_dict = records
    out_records = []
    for record in records.values():
        out_records.append(full_dict[record.id])
    return out_records


def map_species_dict_to_records(records, full_dict, unique=False):
    out_records = []
    for id_list in records.values():
        if unique:
            out_records.append(full_dict[id_list[0]])
        else:
            for seq_id in id_list:
                out_records.append(full_dict[seq_id])
    return out_records


def map_list_to_records(records, full_dict, unique=False):
    out_records = []
    for seq_id in records:
        if unique and seq_id in full_dict:
            out_records.append(full_dict[seq_id])
        elif seq_id in full_dict:
            out_records.append(full_dict[seq_id])
    return out_records


def write_fasta(records, filename, full_details=None):

    if full_details:
        for seq in records:
            seq.id = seq.name



    with open(filename, "w") as handle:
            SeqIO.write(records, handle, "fasta")


def compare_fasta(records1, records2):
    compare_list = []
    for x in records1.keys():
        if x not in records2.keys():
            compare_list.append(x)
    return compare_list


def replace_words(changes, file):
    for k, v in changes:
        file.replace(k, v)
    return file


def get_different_ids(records1, records2):
    ids1 = set(records1.keys())
    ids2 = set(records2.keys())
    return ids1 - ids2


def get_different_records(records1, records2):
    ids = get_different_ids(records1, records2)
    different_dict = {k: records1[k] for k in ids}
    return different_dict


def subset_on_motif(records, motif, retain_motif_seqs = True):
    """
    Subset a set of records based on the presence of a motif.

    :param motif: The motif to check for
    :param retain_motif_seqs: Whether we should retain sequences with the motif or exclude them
    :return: The subseted records based on the presence of the motif
    """
    subset_motifs = {}

    pattern=re.compile(motif)

    for record in records.values():
        if (len(pattern.findall(str(record.seq))) > 0 and retain_motif_seqs) or (len(pattern.findall(str(record.seq))) == 0 and not retain_motif_seqs):
            subset_motifs[record.id] = record

    return subset_motifs


def correct_fastml_nodes(in_path, out_path):
    fastml_ancestors = utilities.load_sequences(in_path)

    for x in fastml_ancestors:
        if len(x) < 5 and x.startswith("N"):
            oldname = fastml_ancestors[x].name
            node = oldname.split("N")[1]
            newname = "N" + str(int(node) - 1)
            fastml_ancestors[x].id = newname
            fastml_ancestors[x].name = newname
            fastml_ancestors[x].description = newname

    records = [x for x in fastml_ancestors.values()]

    write_fasta(records=records, filename=out_path)

def add_gaps_to_fastml_joint(joint_path, marginal_path, out_path):

    joint = utilities.load_sequences(joint_path)
    marginal = utilities.load_sequences(marginal_path)

    for x in marginal:
        if len(x) < 5 and x.startswith("N"):
            updated_joint = ""

            for pos in zip(marginal[x].seq, joint[x].seq):
                if pos[0] == "-":
                    updated_joint += "-"
                else:
                    updated_joint += pos[1]

            joint[x].seq = Seq(updated_joint)

    records = [x for x in joint.values()]
    write_fasta(records=records, filename=out_path)

def write_out_fasta_without_info_after_spaces(records, out_path):
    for x in records:
        x.description = ""
        x.info = ""
    utilities.write_alignment(records, out_path)
