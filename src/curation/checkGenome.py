from Bio import Entrez
from lxml import etree as et
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib.units import cm

Entrez.email = "gabriel.foley@uqconnect.edu.au"


def check_genomic_location(records, min_length=0, visualise=None, file_path=None):
    for species, ids in records.items():
        if len(ids) >= min_length:
            records_query = ""
            record_dict = {}

            for record in ids:
                    # records_query += record + " OR "

                search = Entrez.esearch(term=record, db="gene", retmode="fasta", rettype="acc")
                result = Entrez.read(search)

                print (result)

                if len(result['IdList']) > 0:
                    gene_ids = result['IdList'][0]

                    record_dict[record] = gene_ids

                print ('record dict')
                print (record_dict)

                if len(record_dict) > 1:

                    if not visualise:

                        for protein_id, gene_id in record_dict.items():
                            handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
                            record = handle.read()
                            xml_parsed = et.fromstring(record)

                            start = xml_parsed.xpath(
                                '/Entrezgene-Set/Entrezgene/Entrezgene_track-info/Gene-track/Gene-track_geneid[contains('
                                '., \'' + gene_id +
                                "')]/../../../Entrezgene_locus/Gene-commentary/Gene-commentary_seqs/Seq-loc/Seq-loc_int"
                                "/Seq-interval/Seq-interval_from/text()")
                            finish = xml_parsed.xpath(
                                "/Entrezgene-Set/Entrezgene/Entrezgene_track-info/Gene-track/Gene-track_geneid[contains("
                                "., '" + gene_id +
                                "')]/../../../Entrezgene_locus/Gene-commentary/Gene-commentary_seqs/Seq-loc/Seq-loc_int"
                                "/Seq-interval/Seq-interval_to/text()")
                            chromosome = xml_parsed.xpath(
                                "/Entrezgene-Set/Entrezgene/Entrezgene_track-info/Gene-track/Gene-track_geneid[contains("
                                "., '" + gene_id +
                                "')]/../../../Entrezgene_source/BioSource/BioSource_subtype/SubSource/SubSource_name"
                                "/text()")

                            if file_path is None:
                                print()
                                print("Protein id is %s" % protein_id)
                                print("Gene id is %s " % gene_id)
                                print("Gene region starts at %s" % (start[0]))
                                print("Gene region ends at %s" % (finish[0]))
                                print("Gene region is %s nucleotides long " % (int(finish[0]) - int(start[0])))
                                if chromosome[0] == "Un":
                                    print("Chromosome is unassigned")
                                else:
                                    print("On chromosome %s" % (chromosome[0]))

                                print("-------------------------------------------------------------------")

                            else:
                                with open(file_path + " " + species, "a+") as text_file:
                                    text_file.write("Protein id is %s \n" % protein_id)
                                    text_file.write('Gene'
                                                    ' id is {} \n'.format(gene_id))
                                    text_file.write("Gene region starts at %s \n" % (start[0]))
                                    text_file.write("Gene region ends at %s \n" % (finish[0]))
                                    text_file.write(
                                        "Gene region is %s nucleotides long \n" % (int(finish[0]) - int(start[0])))
                                    print('Chromosome is ')
                                    print(chromosome)
                                    if chromosome:
                                        print("Has value")
                                    if chromosome:
                                        if chromosome[0] == "Un":
                                            text_file.write("Chromosome is unassigned \n")
                                        else:
                                            text_file.write("On chromosome %s \n" % (chromosome[0]))

                                    text_file.write("-------------------------------------------------------------------\n")

                    elif visualise == "circular" or "linear":
                        print ('triggered this')
                        draw_genome(species, record_dict, visualise)

                    else:
                        print("-visualise flag should be either 'circular' or 'linear'")


def draw_genome(species, record_dict, visualise):
    print("drawing genome")
    locations = []

    gdd = GenomeDiagram.Diagram('Test Diagram')
    gdt_features = gdd.new_track(1, greytrack=False)
    gds_features = gdt_features.new_set()

    for protein_id, gene_id in record_dict.items():
        handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
        record = handle.read()
        xml_parsed = et.fromstring(record)

        start = xml_parsed.xpath(
            '/Entrezgene-Set/Entrezgene/Entrezgene_track-info/Gene-track/Gene-track_geneid[contains(., \'' + gene_id +
            "')]/../../../Entrezgene_locus/Gene-commentary/Gene-commentary_seqs/Seq-loc/Seq-loc_int/"
            "Seq-interval/Seq-interval_from/text()")
        finish = xml_parsed.xpath(
            "/Entrezgene-Set/Entrezgene/Entrezgene_track-info/Gene-track/Gene-track_geneid[contains(., '" + gene_id +
            "')]/../../../Entrezgene_locus/Gene-commentary/Gene-commentary_seqs/Seq-loc/Seq-loc_int/"
            "Seq-interval/Seq-interval_to/text()")

        if start and finish:
            locations.append(int(start[0]))
            locations.append(int(finish[0]))
            feature = SeqFeature(FeatureLocation(int(start[0]), int(finish[0])), strand=+1)
            gds_features.add_feature(feature, name=protein_id + "(" + gene_id + ")", label=True)

    gdd.draw(format=visualise, pagesize=(15 * cm, 20 * cm), fragments=1,
             start=0, end=max(locations))
    gdd.write(species + ".pdf", "pdf")

def generate_multiple_hit_data(species_names, species_counts, full_record, file_path):
    id_dict = {}
    #     for name in species_names:
    #         seqs = map_species_to_records(species_counts[name], full_record)
    #         write_fasta(seqs, file_path + name + " sequences")
    #         alignmentFile = alignment.alignWithMAFFT(file_path + name + " sequences")
    #         alignment.writeAlignment(alignmentFile, file_path + name + ".aln", "fasta")


    #     checkGenome.check_genomic_location(species_counts, min_length=2, file_path=file_path +" gene locations ")
    check_genomic_location(species_counts, min_length=2, visualise="linear")
