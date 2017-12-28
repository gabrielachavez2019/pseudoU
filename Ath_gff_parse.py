from Pseudo_lib import *
import pickle
from BCBio import GFF 

Unknown_UTR_length = 200
windows = 150
in_file = "/share/home/xuyuxing/Work/Pseudouracil/Ath/TAIR10_GFF3_genes_transposons.gff"
genome_file = "/share/home/xuyuxing/Work/Pseudouracil/Ath/TAIR10_Chr.all.fasta"
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
                exon_num = exon_num+1
                transcript_dict[transcript.id] = transcript, rec
            if exon_num == 0:
                del transcript_dict[transcript.id]

num = 0
for transcript_id in transcript_dict:
    print num
    num = num+1
    transcript, rec = transcript_dict[transcript_id]

    if transcript.type == "mRNA":
        if gene_dict[transcript.qualifiers['Parent'][0]][0].qualifiers['Note'][0] == "protein_coding_gene":
            UTR_num = 0
            for UTR in transcript.sub_features:
                if UTR.type == "five_prime_UTR":
                    start = int(UTR.location.start)+1
                    end = int(UTR.location.end)
                    UTR_num = UTR_num+1
                    five_prime_UTR = range(start, end+1)
                elif UTR.type == "three_prime_UTR":
                    start = int(UTR.location.start)+1
                    end = int(UTR.location.end)
                    UTR_num = UTR_num+1
                    three_prime_UTR = range(start, end+1)

            cds_range = []
            for CDS in transcript.sub_features:
                if CDS.type == "CDS":
                    start = int(CDS.location.start)+1
                    end = int(CDS.location.end)
                    cds_range.append([start, end])

            cds_node = []
            for i in cds_range:
                cds_node.extend(i)
            cds_node = sorted(cds_node)

            if transcript.strand == 1:
                initiation_codon = range(min(cds_node), min(cds_node)+3)
                termination_codon = range(max(cds_node)-2, max(cds_node)+1)
            else:
                termination_codon = sorted(range(min(cds_node), min(cds_node)+3), reverse=True)
                initiation_codon = sorted(range(max(cds_node)-2, max(cds_node)+1), reverse=True)

            termination_string = ""
            for i in termination_codon:
                site_info = rec+":"+str(i)+"-"+str(i)
                now = sub_fasta(genome_file, site_info)
                if transcript.strand == -1:
                    now = reverse_complement(now)
                termination_string = termination_string+now

            initiation_string = ""
            for i in initiation_codon:
                site_info = rec+":"+str(i)+"-"+str(i)
                now = sub_fasta(genome_file, site_info)
                if transcript.strand == -1:
                    now = reverse_complement(now)
                initiation_string = initiation_string+now

            transcript_pos = []
            for exon in transcript.sub_features:
                if not exon.type == "exon":
                    continue
                start = exon.location.start.position+1
                end = exon.location.end.position
                transcript_pos.extend(range(start, end+1))
            transcript_pos = list(set(transcript_pos))
            transcript_pos = sorted(transcript_pos)
            if UTR_num == 0:
                if transcript.strand == 1:
                    five_prime_UTR = range(transcript_pos[0]-Unknown_UTR_length, transcript_pos[0])
                    three_prime_UTR = range(transcript_pos[-1]+1, transcript_pos[-1]+Unknown_UTR_length+1)
                elif transcript.strand == -1:
                    three_prime_UTR = range(transcript_pos[0]-Unknown_UTR_length, transcript_pos[0])
                    five_prime_UTR = range(transcript_pos[-1]+1, transcript_pos[-1]+Unknown_UTR_length+1)
                transcript_pos.extend(range(transcript_pos[0]-Unknown_UTR_length, transcript_pos[0]))
                transcript_pos.extend(range(transcript_pos[-1]+1, transcript_pos[-1]+Unknown_UTR_length+1))
                transcript_pos = list(set(transcript_pos))
            transcript_dict[transcript_id] = {}
            transcript_dict[transcript_id]['position'] = transcript_pos
            transcript_dict[transcript_id]['strand'] = transcript.strand
            transcript_dict[transcript_id]['ref'] = rec
            transcript_dict[transcript_id]['CDS'] = cds_range
            transcript_dict[transcript_id]['initiation_codon'] = initiation_codon
            transcript_dict[transcript_id]['initiation_string'] = initiation_string
            transcript_dict[transcript_id]['termination_codon'] = termination_codon
            transcript_dict[transcript_id]['termination_string'] = termination_string
            transcript_dict[transcript_id]['five_prime_UTR'] = five_prime_UTR
            transcript_dict[transcript_id]['three_prime_UTR'] = three_prime_UTR

OUT = open("transcript_dict.pyb", 'wb')
pickle.dump(transcript_dict, OUT)
OUT.close()