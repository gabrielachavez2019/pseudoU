import sys
import re
import string

def site_out_parse(project_root, tag, peak):
    #tag = "S02_S01"
    #peak = 1
    tag_list = string.split(tag, "_")
    CMC_p = tag_list[0]
    CMC_p_file = project_root + CMC_p + ".tab"
    CMC_n = tag_list[1]
    CMC_n_file = project_root + CMC_n + ".tab"

    tab_out_file = project_root + tag + ".all_transcript.tab"

    F1 = open(tab_out_file)
    peak_record = {}
    for each_line in F1:
        each_line = re.sub('\n', '', each_line)
        info = string.split(each_line, "\t")
        peak_record[(info[1], info[2])] = (info[0], float(info[3]), float(info[4]))
    F1.close()

    CMC_dict = {}
    CMC_dict[CMC_p] = {}
    CMC_dict[CMC_n] = {}
    
    F1 = open(CMC_p_file)
    for each_line in F1:
        each_line = re.sub('\n', '', each_line)
        chr,site,start_1,pass_1,start__1,pass__1 = string.split(each_line, "\t")
        CMC_dict[CMC_p][(chr, site)] = {}
        CMC_dict[CMC_p][(chr, site)][1] = [start_1, pass_1]
        CMC_dict[CMC_p][(chr, site)][-1] = [start__1, pass__1]
    F1.close()

    F1 = open(CMC_n_file)
    for each_line in F1:
        each_line = re.sub('\n', '', each_line)
        chr,site,start_1,pass_1,start__1,pass__1 = string.split(each_line, "\t")
        CMC_dict[CMC_n][(chr,site)]=[start_1,pass_1,start__1,pass__1]
        CMC_dict[CMC_n][(chr, site)][1] = [start_1, pass_1]
        CMC_dict[CMC_n][(chr, site)][-1] = [start__1, pass__1]
    F1.close()

    site_out_file = project_root + tag + ".all_transcript.site.out"
    record_now = {}
    F1 = open(site_out_file)
    for each_line in F1:
        each_line = re.sub('\n', '', each_line)
        info = string.split(each_line, "\t")
        ID = int(info[0])
        contig = info[1]
        transcript = info[2]
        strand = int(info[3])
        up_site = int(info[4])
        site_raw = int(info[5])
        down_site = int(info[6])
        seq_interval = info[7]
        best_sup = int(info[8])
        best_rank = int(info[9])
        search_peak = (info[2], info[5])
        search_chr = (info[1], info[5])
        if float(peak_record[search_peak][1]) >= peak:
            CMC_p_start, CMC_p_pass = CMC_dict[CMC_p][search_chr][strand][0], CMC_dict[CMC_p][search_chr][strand][1]
            if CMC_dict[CMC_n].has_key(search_chr):
                CMC_n_start, CMC_n_pass = CMC_dict[CMC_n][search_chr][strand][0], CMC_dict[CMC_n][search_chr][strand][1]
            else:
                CMC_n_start, CMC_n_pass = 0, 0
            record_now[ID] = (ID, contig, transcript, strand, up_site, site_raw, down_site, seq_interval, best_sup, best_rank, peak_record[search_peak][1], peak_record[search_peak][2], CMC_p_start, CMC_p_pass, CMC_n_start, CMC_n_pass)
    F1.close()
    return record_now

tag = sys.argv[1]
project_root = "/share/home/xuyuxing/Work/Pseudouracil/python/1.tab/" 
record_get = site_out_parse(project_root, tag, 1.0)

OUT = open(project_root+tag+".all_transcript.report", "w")
for ID in record_get:
    (ID, contig, transcript, strand, up_site, site_raw,down_site,seq_interval, best_sup, best_rank, peak_p, peak_n, CMC_p_start, CMC_p_pass, CMC_n_start, CMC_n_pass) = record_get[ID]
    printer = str(ID)+"\t"+str(contig)+"\t"+str(transcript)+"\t"+ str(strand)+"\t"+str(up_site)+"\t"+str(site_raw)+"\t"+str(down_site)+"\t"+ str(seq_interval)+"\t"+str(best_sup)+"\t"+str(best_rank)+"\t"+str(peak_p)+"\t"+str(peak_n)+"\t"+str(CMC_p_start)+"\t"+str(CMC_p_pass)+"\t"+str(CMC_n_start)+"\t"+str(CMC_n_pass)
    print >>OUT, printer
OUT.close()