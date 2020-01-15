from ExonRecord import ExonRecord
from Bio import SeqIO, Entrez
from collections import defaultdict
import fasta
import re
import math
from urllib.request import urlopen
import utilities
from bs4 import BeautifulSoup as Soup
from urllib.error import HTTPError
import csv
import numpy

Entrez.email = "gabriel.foley@uqconnect.edu.au"


def map_exons(records, skipped_records_path):
    """
    Take a set of protein sequences and map them back to their exon coordinates from their genomic location
    :param records: The set of protein sequences
    :return: The exon coordinates
    """
    genomic_records = {}
    skipped_records = []

    for record in records:
        gene_id = genome_id = genome_record = exon_num = None
        mrna = False
        search_id = record.id
        print('Search id is ')
        print(search_id)

        # If it isn't an NCBI sequence, lets try and map it to the UniProt database
        if record.annotations["Database"] != "NCBI":
            search_id = map_from_uniprot(skipped_records, search_id, record)

            # If we can't map to a UniProt record, lets just check it isn't actually an NCBI ID
            if not search_id:
                search_id = record.id

        try:


            if not check_in_alternative_databases(genomic_records, skipped_records, search_id, record):


                protein_record = get_protein_record(search_id)

                coded_by = protein_record.find("gbseq_source-db").getText().split(" ")[-1]

                exon_location = get_exon_location_from_protein_record(protein_record, search_id)


                # Sanity check that our database source and the coded_by field on the protein record match (they should)

                exon_string = strip_text(exon_location, ["join(", "complement("]) # Remove leading info
                if exon_string[0:len(coded_by)] != coded_by:
                    # print (exon_string[0:len(coded_by)], coded_by)
                    raise RuntimeError("The record this protein is coded by doesn't match the database source")

                genomic_record = map_to_genomic_record(genomic_records, skipped_records, coded_by, record)

                if genomic_record:
                    # print ("Got a genomic record")
                    mrna_check = genomic_record.find('gbseq_moltype', text="mRNA")
                    if mrna_check:
                        # print ("This was an mRNA record")
                        gene_id = get_gene_id(search_id)

                        if gene_id:
                            # print ("Gene ID")

                            gene_record = get_gene_record(gene_id)

                            # print ('here')

                            exon_num = get_exon_number_from_gene_record(gene_record)
                            genome_id = get_genome_id(gene_record)

                        if genome_id:
                            genome_record = get_genome_record(gene_record, genome_id)

                        if genome_record:
                            exon_location = get_exon_location_from_genomic_record(genome_record, search_id, gene_id)

                        else:
                            skipped_records.append({record.id: "mRNA record"})
                            mrna = True

                # print ('state is')
                # print (exon_location, exon_num, mrna)
                if exon_location and not mrna:
                    # print ('here')
                    exon_record = build_exon_record(exon_location, search_id)
                    # print ('exon recortd is')
                    # print (exon_record)

                    if exon_record:
                        genomic_records[record.id] = exon_record

                    # We could get an exon number but not the locations
                    # elif exon_num:
                    #     # print ('got a num')
                    #     exon_record = get_dummy_record(search_id, exon_num)
                    #     # print (exon_record)
                    #     genomic_records[record.id] = exon_record

                    else:
                        print("Couldn't find an exon location in the genomic record")
                        skipped_records.append({record.id: "No exon location"})

        except:
            skipped_records.append({record.id: "Couldn't map record"})
    if skipped_records_path:
        write_skipped_records(skipped_records, skipped_records_path)

    return genomic_records


def get_protein_record(search_id):
    handle = Entrez.efetch(db="protein", id=search_id, rettype="gb", retmode="xml")
    protein_record = Soup(handle, "lxml")
    return protein_record


def get_gene_record(gene_id):
    handle = Entrez.efetch(db="gene", id=gene_id, rettype="gb", retmode="xml")
    gene_record = Soup(handle, "lxml")
    return gene_record

def get_genome_record(gene_record, genome_id):

    seq_from = gene_record.find("seq-interval_from")
    seq_to = gene_record.find("seq-interval_to")


    handle = Entrez.efetch(db="nuccore", id=genome_id, rettype="gb", retmode="xml",
                           seq_start=seq_from,
                           seq_stop=seq_to)
    genome_record = Soup(handle, "lxml")
    return genome_record

def get_gene_id(search_id):


    handle = Entrez.elink(dbfrom="protein", db='gene', id=search_id, rettype='xml')
    mapping = Entrez.read(handle, validate=False)


    # Retrieve the gene ID
    if mapping:
        for term in mapping:
            if term['LinkSetDb']:
                if term['LinkSetDb'][0]['Link'][0]['Id']:
                    gene_id = term['LinkSetDb'][0]['Link'][0]['Id']
            else:
                gene_id = None
                break
    print (gene_id)
    return gene_id


def get_genome_id(gene_record):
    genome_id_list = gene_record.find('gene-commentary_accession')
    # print (genome_id_list.getText())
    genome_id = genome_id_list.getText()
    # print (gene_record.prettify())
    # nc = nw = nt = None
    # for check_id in genome_id_list:
    #     print (check_id)
    #     id_text = check_id.getText()
    #     if id_text[0:2] == 'NC':
    #         nc = id_text
    #     elif id_text[0:2] == 'NW':
    #         nw = id_text
    #     elif id_text[0:2] == 'NT':
    #         nt = id_text
    #
    # if nc:
    #     genome_id = nc
    # elif nw:
    #     genome_id = nw
    # elif nt:
    #     genome_id = nt
    #
    # print (genome_id)

    return genome_id

def get_exon_location_from_protein_record(protein_record, search_id):
    coded_by = protein_record.find("gbqualifier_name", text="coded_by")

    for sibling in coded_by.next_siblings:
        if sibling != "\n":
            exon_location = sibling.text

    return exon_location

def get_exon_number_from_gene_record(gene_record):
    # print (gene_record.prettify())
    # print ('checking')
    exon_num = None
    feature = gene_record.find('gene-commentary_label', text="Exon count")
    # print (feature)
    num = feature.find_next('gene-commentary_text').getText()
    if (num):
        exon_num = num
    print ('exon num is', exon_num)
    return exon_num

def get_exon_location_from_genomic_record(genomic_record, search_id, gene_id):

    # print ('in here')
    # print (genomic_record.prettify())


    exon_location = None
    gb_features = genomic_record.find_all('gbfeature_key', text="CDS")

    for feature in gb_features:
        qualifier = feature.find_next('gbqualifier_value')
        # if qualifier == search_id:
        for qualifier in feature.find_all_next('gbqualifier_value'):
            if qualifier.getText() == search_id:
                for location in feature.find_next('gbfeature_location'):
                    exon_location = (location)
                    return exon_location

            gene_check = qualifier.getText().split(":")
            if gene_check[0] == "GeneID" and gene_check[1] == gene_id:
                # print ('found it')
                for location in feature.find_next('gbfeature_location'):
                    exon_location = (location)
                    return exon_location

    return exon_location



def strip_text(string_to_clean, text_list):
    for text in text_list:
        string_to_clean = string_to_clean.replace(text, "")

    return string_to_clean



def map_to_genomic_record(genomic_records, skipped_records, coded_by, record):


    handle = Entrez.efetch(db="nuccore", id=coded_by, rettype="gb", retmode="xml")

    soup = Soup(handle, "lxml")
    return soup


def map_from_uniprot(skipped_records, search_id, record):
    print ("Mapping from UniProt")

    try:
        handle = urlopen("https://www.uniprot.org/uniprot/" + search_id + ".xml")
        soup = Soup(handle, "lxml")
        print ('got it')
        search_id = soup.find('property', type="protein sequence ID")

        if search_id:
            search_id = search_id["value"]
            return search_id

        else:
            print("Couldn't map the record to an NCBI or Ensembl Fungi Genome automatically")
            return False

    except HTTPError as error:
        print("Couldn't access a Uniprot record for " + search_id)
        return False

def check_in_alternative_databases(genomic_records, skipped_records, search_id, record):
    # print ("Checking in alternative databases")

    found = False

    # Check if this is an ID that maps to the Ensembel Fungi Genome database
    if (search_id[:4] == "FOXG"):
        found = True
        handle = urlopen("http://www.ensemblgenomes.org/id/" + search_id)
        soup = Soup(handle, "lxml")

        transcript_info = soup.find('div', attrs={'class': 'lhs'}, text="About this transcript")

        for sibling in transcript_info.next_siblings:
            if sibling != "\n":
                exon_num = int(sibling.text.split("This transcript has ")[1].split(" ")[0])


        exon_record = get_dummy_record(search_id, exon_num)


        genomic_records[record.id] = exon_record

    # For Steph's CYP2D sequences that had ENSEMBL Genome records
    elif (search_id[:3] == "ENS"):
        found = True
        print("Couldn't find a genome ID")
        print("SKIPPING")
        skipped_records.append({record.id: "Couldn't map record"})

    return found

def get_dummy_record(search_id, exon_num):
    # print ('lets do')
    exons = generate_dummy_exons(int(exon_num))

    strand = "plus"

    exon_record = ExonRecord(protein_id=search_id, gene_id="", exons=exons,
                             strand=strand, calc_introns=False)

    # print (exon_record)

    return exon_record

def build_exon_record(exon_location, search_id):
    # print ('lets build')
    # print (exon_location)
    if 'join' in exon_location:
        exons = (exon_location.split('join(')[1].split(','))
        # print (exons)

        # If the exon has the record ID preceeding it, remove it
        for count, exon in enumerate(exons):
            if ":" in exon:
                exons[count] = re.split("[.]\d:", exon)[1]

    elif ":" in exon_location:
        # print ('lets')
        exons = [re.split("[.]\d:", exon_location)[1]]
        # print ('do it')

    else:
        return None

    if "complement" in exon_location:
        strand = "minus"
        for num, exon in enumerate(exons):
            exons[num] = re.sub(r'complement', r'', exon)
    else:
        strand = "plus"

    # print (exons)

    if exons:

        exon_record = ExonRecord(protein_id=search_id,
                                       exons=exons,
                                       strand=strand, calc_introns=True)
    else:
        return None

    print ("Exon count is %s" % (exon_record.exon_count))

    return exon_record



def get_feature_counts(records):
    """
    Takes an exon records dictionary and returns a dictionary with exon count as the key and a list of the speices that
    have that number exons as the value

    :param records: The dictionary containing all of the exon records
    :return: Dictionary with the exon counts as the key
    """

    exon_counts = defaultdict(list)



    for protein_id, genomic_record in records.items():
        exon_counts[genomic_record.exon_count].append(protein_id)
    return exon_counts


def map_exon_count_to_tree(exons, tree):
    pass


def append_tag_to_seqs_given_exon_count(*numbers, exon_counts, full_record, outpath):
    outfile = []

    for count in numbers:
        if count in exon_counts:
            for seq_id in exon_counts[count]:
                if seq_id in full_record:
                    full_record[seq_id].id += "*TAG " + str(count) + "*"

    for record in full_record.values():
        outfile.append(record)

    fasta.write_fasta(outfile, outpath)


def write_out_seqs_given_exon_count(*numbers, exon_counts, full_record, outpath):
    exon_records = []
    for count in numbers:
        exon_records += exon_counts[count]

    outfile = fasta.map_list_to_records(exon_records, full_record)

    fasta.write_fasta(outfile, outpath)


def save_genomic_records(records, filepath, skipped_records_path=None):
    # If records is a dictionary convert it into a list of the records

    if type(records) == dict:
        records = [record for record in records.values()]

    genomic_records = map_exons(records, skipped_records_path)

    utilities.save_python_object(genomic_records, filepath)



def open_genomic_records(filepath):
    genomic_records = utilities.open_python_object(filepath)
    return genomic_records


def get_exon_counts(records, genomic_records, filter_records=None, exclude=True):
    exon_counts = {}
    if not filter_records:
        filter_records = []

    if type(records) == dict:
        records = [record for record in records.values()]

    for record in records:

        if record.id in genomic_records:
            if (exclude and record.id not in filter_records) or (not exclude and record.id in filter_records):
                exon_counts[record.id] = len(genomic_records[record.id].exons)
    return exon_counts

def generate_dummy_exons(exon_num):
    exons = []

    for i in range(0, exon_num):
        exons.append("id.1:1..1")

    for count, exon in enumerate(exons):
        exons[count] = re.split("[.]\d:", exon)[1]
    return exons

def get_exon_array(exon_dict):
    """
    Convert an exon dictionary into an numpy exon array
    :param exon_dict:
    :return:
    """

    exon_array = numpy.array([record.exon_count for record in exon_dict.values()])
    return exon_array


def write_exon_counts_to_csv(records, filepath):
    print ("writing to ", filepath)
    with open(filepath, 'w+') as csvfile:
        csvfile.write("Sequence ID, Exon count\n")
        for seq_id, record in records.items():

            csvfile.write("%s,%s\n" % (str(seq_id), str(record.exon_count)))
    csvfile.close()

def write_skipped_records(records, filepath):
    with open(filepath, 'w+') as csvfile:
        csvfile.write("Sequence ID, Reason for exclusion\n")

        for record in records:

            for name, reason in record.items():
                csvfile.write("%s,%s\n" % (name, reason))



def write_exon_totals_to_csv(records, filepath):
    with open(filepath, 'w+') as csvfile:
        csvfile.write("Record ID, Number of sequences, Exon count (mean), Standard deviation of exon count\n")
        for record_id, record in records.items():
            csvfile.write("%s,%s,%s,%s\n" % (record_id, len(record), record.mean(), record.std()))


def map_exon_boundaries_to_alignment(records, genomic_records, filter_records=None, exclude=True):

    if not filter_records:
        filter_records = []

    # If records is a dictionary convert it into a list of the records
    if type(records) == dict:
        records = [record for record in records.values()]

    # exons = {'XP_019684690.2' : [178, 167, 161, 635, 489], 'XP_009883824.1' : [178, 167, 161, 705, 82],
    # 'XP_014065111.1' : [181, 167, 161, 334, 318, 441] }
    # exons = {'XP_019684690.2' : [3, 3, 2, 1], 'XP_009883824.1' : [2, 2, 1], 'XP_014065111.1' : [4] }
    #

    # genomic_records = map_exons(records)

    cols = ["\033[1;35;4m", "\033[1;34;4m", "\033[1;32;4m", "\033[1;33;4m", "\033[1;37;4m", "\033[1;31;4m",
            "\033[1;39;4m", "\033[1;29;4m", "\033[1;21;4m", "\033[1;22;4m", "\033[1;23;4m", "\033[1;24;4m", "\033[1;25;4m",]
    longest_header = 0
    for record in records:
        # print (record.id)
        if len(record.id) > longest_header:
            longest_header = len(record.id)

        buildseq = ""
        newseq = ""

        counter = 0

        # print (genomic_records[record].strand)
        # print (genomic_records[record].exon_lengths)
        # print (genomic_records[record].exons)




        if record.id in genomic_records:
            if (exclude and record.id not in filter_records) or (not exclude and record.id in filter_records):
                # print ('here')
                #

                # print (genomic_records[record].exon_lengths)
                # print (genomic_records[record].strand)
                for exon in genomic_records[record.id].exon_lengths:
                    # print (exon)
                    # print (exon / 3)
                    counter += (exon / 3)
                # print (counter)

                for count, exon in enumerate(genomic_records[record.id].exon_lengths):
                    # print (count, exon)
                    # print (count,int(math.ceil(exon/3)))
                    buildseq += str(count + 1) * int(math.ceil((exon / 3)))
                # print(record.id)
                # print (buildseq)
                # print (str(record.seq).replace("-", ""))
                # print (record.seq)

                try:
                    ind = 0
                    for pos in record.seq:
                        if pos == "-":
                            newseq += "-"
                        elif pos != "-":
                            # print (ind)
                            newseq += buildseq[ind]
                            ind += 1

                #     print(record.id)
                #     filter_records.append(record.id)
                #     print(filter_records)


                    # print (newseq)

                    # print (records[record].seq)
                    seq_with_exons = str(record.seq)
                    # print (newseq)
                    # print (len(genomic_records[record].exon_lengths))
                    # print (genomic_records[record.id].exon_lengths)
                    for exon in range(0, len(genomic_records[record.id].exon_lengths)):
                        # print (exon)
                        # print ('ind', ind)
                        # print ('exon', exon + 1)
                        # print (newseq)

                        ind = (newseq.index(str(exon + 1)))

                        # print (len(seq_with_exons))
                        seq_with_exons = seq_with_exons[:ind + exon] + "*" + seq_with_exons[ind + exon:]
                        # seq_with_exons =  str(exon) + "*" + seq_with_exons[:ind + 2 * exon]

                    for exon in range(0, len(genomic_records[record.id].exon_lengths)):
                        # print (seq_with_exons)
                        seq_with_exons = seq_with_exons.replace("*", cols[exon], 1)
                        # print ('hereeeee')
                        # print (seq_with_exons)
                    print(cols[-1] + '{message: <{fill}}'.format(message=record.id, fill=longest_header),
                              seq_with_exons)

                except IndexError:
                    pass


                # except:
                #     print (cols[-1] + '{message: <{fill}}'.format(message=record.id, fill=longest_header), record.seq)

    print (filter_records)