import subprocess
import string
import re

def sub_fasta(genome_file,site_info):
    """
    site_info = "Chr1:1-2000"
    file_string = "../Ath/TAIR10_Chr.all.fasta"
    """
    match = re.match(r'(\S+):(\d+)-(\d+)',site_info)
    chr = match.group(1)
    start = int(match.group(2))
    end = int(match.group(3))

    chr_length={}
    F1 = open(genome_file+".fai")
    for each_line in F1:
        each_line = re.sub('\n', '', each_line)
        info = string.split(each_line, "\t")
        chr_length[info[0]]=int(info[1])
    F1.close()

    chr_length = chr_length[chr]

    if end > chr_length:
        end = chr_length

    site_info = chr + ":" + str(start) + "-" + str(end)

    cmd_string = "samtools faidx "+genome_file+" "+site_info
    samtools_output = subprocess.check_output(cmd_string,shell=True)
    samtools_output = string.split(samtools_output,"\n")
    output_string=""
    for i in samtools_output:
        if not re.match(r'^>',i):
            i = re.sub('\n','',i)
            output_string=output_string+i
    return output_string

def reverse_complement(seq):
    """
    #>>> seq = "TCGGinsGCCC"
    #>>> print "Reverse Complement:"
    #>>> print(reverse_complement(seq))
    #GGGCinsCCGA
    """
    alt_map = {'ins':'0'}
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    for k,v in alt_map.iteritems():
        seq = seq.replace(k,v)
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.iteritems():
        bases = bases.replace(v,k)
    return bases
