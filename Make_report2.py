import sys
import re
import string
import Pseudo_lib
import pickle
from BCBio import GFF
import csv


def Get_relative_site(base_point, site, strand):
    if int(strand) == 1:
        if site >= base_point:
            relative_site = site - base_point + 1
        else:
            relative_site = site - base_point
    else:
        if site > base_point:
            relative_site = base_point - site
        else:
            relative_site = base_point - site + 1

    return relative_site


tag = sys.argv[1]
# tag = "S02_S01"
project_root = "/share/home/xuyuxing/Work/Pseudouracil/python/1.tab/"
report_file = project_root + tag + ".all_transcript.report"
genome_file = "/share/home/xuyuxing/Work/Pseudouracil/Ath/TAIR10_Chr.all.fasta"

in_file = "/share/home/xuyuxing/Work/Pseudouracil/Ath/TAIR10_GFF3_genes_transposons.gff"
in_handle = open(in_file)
gene_dict = {}
for rec in GFF.parse(in_handle):
    for gene in rec.features:
        gene_dict[gene.id] = gene, rec.id
in_handle.close()

transcript_dict = {}
for gene_id in gene_dict:
    gene, rec = gene_dict[gene_id]
    if not gene.type == "exon":
        for transcript in gene.sub_features:
            transcript_dict[transcript.id] = {}
            exon_num = 0
            for exon in transcript.sub_features:
                if not exon.type == "exon":
                    continue
                exon_num = exon_num + 1
                transcript_dict[transcript.id] = transcript, rec
            if exon_num == 0:
                del transcript_dict[transcript.id]

gene_type = {}
for i in transcript_dict:
    parent = transcript_dict[i][0].qualifiers['Parent'][0]
    type = gene_dict[parent][0].type
    Note = gene_dict[parent][0].qualifiers['Note'][0]
    gene_type[i] = Note

rRNA_dict = {}
rRNA_dict['chr2_rRNA'] = {}
rRNA_dict['chr2_rRNA']['position'] = range(3619, 9608 + 1)
rRNA_dict['chr2_rRNA']['ref'] = 'Chr2'
rRNA_dict['chr2_rRNA']['strand'] = 1
rRNA_dict['chr3_rRNA'] = {}
rRNA_dict['chr3_rRNA']['position'] = range(14197590, 14203579 + 1)
rRNA_dict['chr3_rRNA']['ref'] = 'Chr3'
rRNA_dict['chr3_rRNA']['strand'] = 1

TEMP = open(project_root + "transcript_dict.pyb", 'rb')
transcript_dict_1 = pickle.load(TEMP)
TEMP.close()

csv_reader = csv.reader(open('/share/home/xuyuxing/Work/Pseudouracil/Ath/Ath_anno.csv'))
Ath_anno = {}
for row in csv_reader:
    Ath_anno[row[0]] = row[1:]

##################################################################################################################

F1 = open(report_file)
dict_report1 = {}
for each_line in F1:
    each_line = re.sub('\n', '', each_line)
    info = string.split(each_line, "\t")
    (ID, contig, transcript, strand, up_site, site_raw, down_site, seq_interval, best_sup, best_rank, peak_p, peak_n,
     CMC_p_start, CMC_p_pass, CMC_n_start, CMC_n_pass) = info
    dict_report1[int(ID)] = (
    contig, transcript, strand, up_site, site_raw, down_site, seq_interval, best_sup, best_rank, peak_p, peak_n,
    CMC_p_start, CMC_p_pass, CMC_n_start, CMC_n_pass)
F1.close()

OUT = open(project_root + tag + ".all_transcript.report2", "w")

printer = "contig" + "\t" + "trans" + "\t" + "gene_type" + "\t" + "strand" + "\t" + "site" + "\t" + "relative_site" + "\t" + "peak_p" + "\t" + "peak_n" + "\t" + "CMC_p_start" + "\t" + "CMC_p_pass" + "\t" + "CMC_n_start" + "\t" + "CMC_n_pass" + "\t" + "seq_interval" + "\t" + 'five_prime_UTR' + "\t" + 'initiation_codon' + "\t" + 'CDS' + "\t" + 'termination_codon' + "\t" + 'three_prime_UTR' + "\t" + "gene_name" + "\t" + "description" + "\t" + "GO"

print >> OUT, printer

for i in dict_report1:
    (contig, transcript, strand, up_site, site_raw, down_site, seq_interval, best_sup, best_rank, peak_p, peak_n,
     CMC_p_start, CMC_p_pass, CMC_n_start, CMC_n_pass) = dict_report1[i]
    trans, site = transcript, int(best_sup)

    if not best_rank == "1":
        continue

    if site - 10 < 1:
        site__10 = 1
    else:
        site__10 = site - 10

    site_10 = site + 10
    pre_site_info = contig + ":" + str(site__10) + "-" + str(site - 1)
    sub_site_info = contig + ":" + str(site + 1) + "-" + str(site_10)
    site_info = contig + ":" + str(site) + "-" + str(site)
    pre_site_string = Pseudo_lib.sub_fasta(genome_file, pre_site_info)
    sub_site_string = Pseudo_lib.sub_fasta(genome_file, sub_site_info)
    site_string = Pseudo_lib.sub_fasta(genome_file, site_info)

    if int(strand) == -1:
        site_string = Pseudo_lib.reverse_complement(site_string)
        temp_string = Pseudo_lib.reverse_complement(pre_site_string)
        pre_site_string = Pseudo_lib.reverse_complement(sub_site_string)
        sub_site_string = temp_string

    seq_interval = pre_site_string + '"' + site_string + '"' + sub_site_string

    if trans == 'chr2_rRNA' or trans == 'chr3_rRNA':
        base_point = rRNA_dict[trans]['position'][0]
        relative_site = Get_relative_site(base_point, site, strand)
        printer = str(contig) + "\t" + str(trans) + "\trRNA\t" + str(strand) + "\t" + str(site) + "\t" + str(
            relative_site) + "\t" + str(peak_p) + "\t" + str(peak_n) + "\t" + str(CMC_p_start) + "\t" + str(
            CMC_p_pass) + "\t" + str(CMC_n_start) + "\t" + str(CMC_n_pass) + "\t" + str(
            seq_interval) + "\t" + "\t" + "\t" + "\t" + "\t"
    else:
        if gene_type[trans] == 'protein_coding_gene':
            base_point = transcript_dict_1[trans]["initiation_codon"][0]
            relative_site = Get_relative_site(base_point, site, strand)

            info_dict = {}
            types = ['five_prime_UTR', 'initiation_codon', 'termination_codon', 'three_prime_UTR']
            for j in types:
                if j == "initiation_codon":
                    if not transcript_dict_1[trans]['initiation_string'] == "ATG":
                        info_dict['initiation_codon'] = 0
                    else:
                        if site in transcript_dict_1[trans]['initiation_codon']:
                            info_dict['initiation_codon'] = 1
                        else:
                            info_dict['initiation_codon'] = 0
                elif j == "termination_codon":
                    if not transcript_dict_1[trans]['termination_string'] == "ATG":
                        info_dict['termination_codon'] = 0
                    else:
                        if site in transcript_dict_1[trans]['termination_codon']:
                            info_dict['termination_codon'] = 1
                        else:
                            info_dict['termination_codon'] = 0
                else:
                    if site in transcript_dict_1[trans][j]:
                        info_dict[j] = 1
                    else:
                        info_dict[j] = 0

            flag_cds = 0
            for j in transcript_dict_1[trans]['CDS']:
                if site <= max(j) and site >= min(j):
                    flag_cds = 1
            info_dict['CDS'] = flag_cds

            types = ['five_prime_UTR', 'initiation_codon', 'CDS', 'termination_codon', 'three_prime_UTR']

            printer = str(contig) + "\t" + str(trans) + "\tprotein_coding_gene\t" + str(strand) + "\t" + str(
                site) + "\t" + str(relative_site) + "\t" + str(peak_p) + "\t" + str(peak_n) + "\t" + str(
                CMC_p_start) + "\t" + str(CMC_p_pass) + "\t" + str(CMC_n_start) + "\t" + str(CMC_n_pass) + "\t" + str(
                seq_interval)

            for j in types:
                printer = printer + "\t" + str(info_dict[j])

            printer = printer.rstrip("\t")

        else:
            if int(strand) == 1:
                base_point = (transcript_dict[trans][0].location.start) + 1
            else:
                base_point = (transcript_dict[trans][0].location.end)

            relative_site = Get_relative_site(base_point, site, strand)

            printer = str(contig) + "\t" + str(trans) + "\t" + gene_type[trans] + "\t" + str(strand) + "\t" + str(
                site) + "\t" + str(relative_site) + "\t" + str(peak_p) + "\t" + str(peak_n) + "\t" + str(
                CMC_p_start) + "\t" + str(CMC_p_pass) + "\t" + str(CMC_n_start) + "\t" + str(CMC_n_pass) + "\t" + str(
                seq_interval + "\t" + "\t" + "\t" + "\t" + "\t")

    match = re.match(r'(AT.......)\..', trans)
    if match:
        Ath_id = match.group(1)
    if Ath_anno.has_key(Ath_id):
        anno = Ath_anno[Ath_id]
        for j in anno:
            printer = printer + "\t" + j
    print >> OUT, printer
OUT.close()
