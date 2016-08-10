import re
import xml.etree.ElementTree as ET
from urllib2 import urlopen
import requests
import sys


def get_pseudogene_ids(gene):
    ensg = ""
    request_url = urlopen(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=%s[Gene Name]+AND+"
        "Homo+sapiens[Organism]" % gene
    )
    tree = ET.parse(request_url)
    root = tree.getroot()
    pseudogene_ids = []
    pseudo_seqs = []
    pseudo_names = []

    for id_list in root.findall('IdList'):
        gene_id = id_list.find('Id').text
        gene_url = urlopen(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=related_functional_gene_%s[group]"
            % gene_id
        )
        pseudo_tree = ET.parse(gene_url)
        for parent in pseudo_tree.getiterator('IdList'):
            for child in parent:
                id_no = child.text
                pseudogene_ids.append(id_no)
        ensg_url = urlopen(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=%s&rettype=xml&retmode=text" % gene_id
        )
        ensg_tree = ET.parse(ensg_url)
        ensg_root = ensg_tree.getroot()
        for item in ensg_root.findall(
                'Entrezgene/Entrezgene_gene/Gene-ref/Gene-ref_db/Dbtag/Dbtag_tag/Object-id/Object-id_str'
        ):
            if re.match("ENSG(.*)", item.text):
                ensg = item.text

    no_of_pseudogenes = len(pseudogene_ids)
    gene_seq = get_sequence(ensg)
    for pseudogene in pseudogene_ids:
        name = get_pseudogene_name(pseudogene)
        pseudo_names.append(name)
        start, end, chrom = get_gi(pseudogene)
        pseudo_seq = get_pseudo_seq(chrom, start, end)
        pseudo_seqs.append(pseudo_seq)

    pseudogene_dict = dict(zip(pseudogene_ids, zip(pseudo_names, pseudo_seqs)))

    return no_of_pseudogenes, pseudogene_dict, gene_seq


def get_sequence(ensg):
    server = "https://rest.ensembl.org"
    ext = "/sequence/id/%s?" % ensg
    r = requests.get(server+ext, headers={"Content-Type": "text/plain"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    gene_seq = r.text

    return gene_seq


def get_pseudogene_name(pseudo_id):
    pseudogene_name = None
    request_url = urlopen("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=%s" % pseudo_id)
    tree = ET.parse(request_url)
    root = tree.getroot()
    for item in root.findall('DocumentSummarySet/DocumentSummary/Name'):
        pseudogene_name = item.text
    return pseudogene_name


def get_gi(pseudogene_id):
    seq_from = None
    seq_to = None
    chrom = None
    request_url = urlopen(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=%s&rettype=xml&retmode=text"
        % pseudogene_id
    )
    tree = ET.parse(request_url)
    root = tree.getroot()
    for item in root.find(
            'Entrezgene/Entrezgene_locus/Gene-commentary/Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval'
    ):
        if item.tag == "Seq-interval_from":
            seq_from = item.text
        elif item.tag == "Seq-interval_to":
            seq_to = item.text
    if seq_from < seq_to:
        start = seq_from
        end = seq_to
    else:
        start = seq_to
        end = seq_from
    for item in root.find(
        'Entrezgene/Entrezgene_source/BioSource/BioSource_subtype/SubSource'
    ):
        if item.tag == "SubSource_name":
            chrom = item.text
    return start, end, chrom


def get_pseudo_seq(chrom, start, end):
    pseudo_seq = None
    server = "https://rest.ensembl.org"
    ext = "/sequence/region/human/%s:%s..%s?" % (chrom, start, end)
    r = requests.get(server+ext, headers={"Content-Type": "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()
    for key, value in decoded.iteritems():
        if key == "seq":
            pseudo_seq = value

    return pseudo_seq


def create_fasta(pseudogene_dict, gene_seq, gene):
    with open("%s.txt" % gene, "w+") as fasta_file:
        for key, value in pseudogene_dict.iteritems():
            fasta_file.write(">%s_%s\n%s\n" % (key, value[0], value[1]))
        fasta_file.write(">%s_gene\n%s\n" % (gene, gene_seq))

no_of_pseudogenes, pseudogene_dict, gene_seq = get_pseudogene_ids("SDHD")
create_fasta(pseudogene_dict, gene_seq, "SDHD")
