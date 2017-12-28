import string
import re

def parse_report2_file(file_name):
    F1 = open(file_name)
    dict_output = {}
    num=0
    for each_line in F1:
        num = num + 1
        each_line = re.sub('\n', '', each_line)
        info = string.split(each_line, "\t")
        if num == 1:
            title = info
            continue
        dict_output[(info[1],info[4])] = info
    return dict_output,title

def find_more_than(tag_list,more_than_num):
    import collections
    set_all_element = dict(collections.Counter(tag_list))
    output=[]
    for i in set_all_element:
        if set_all_element[i] >= more_than_num:
            output.append(i)
    return output

work_dir = "/share/home/xuyuxing/Work/Pseudouracil/python/1.tab/"

tag = "S28_S27"

suffix = ".all_transcript.report2"

dict_temp,title = parse_report2_file(work_dir+tag+suffix)

intersection_list={"WT":("S02_S01","S06_S05","S10_S26","S25_S09"),
"sT":("S04_S03","S08_S07","S12_S11","S28_S27"),
"Wm":("S18_S17","S22_S30","S29_S21","S32_S31"),
"sm":("S14_S13","S16_S15","S20_S19","S24_S23")}

dict_stat={}
for cond in intersection_list:
    dict_now = {}
    cond_list=[]
    for rep in intersection_list[cond]:
        dict_output, title = parse_report2_file(work_dir+rep+suffix)
        dict_now[rep]=dict_output
        for i in dict_output:
            cond_list.append(i)

    title_print = title[:6]+title[12:]
    for rep in intersection_list[cond]:
        for i in title[6:12]:
            title_print.append(rep+"_"+i)

    title_print_string = ""
    for i in title_print:
        title_print_string = title_print_string + i + "\t"
    title_print_string = title_print_string.strip("\t")

    OUT = open(work_dir+cond+".xls","w")
    print >> OUT, title_print_string

    candidate = find_more_than(cond_list, 3)
    for i in candidate:
        rep_printer = []
        for rep in intersection_list[cond]:
            if dict_now[rep].has_key(i):
                rep_list = dict_now[rep][i][6:12]
                common_list = dict_now[rep][i][:6]+dict_now[rep][i][12:]
            else:
                rep_list = [0 for x in range(0, 6)]
            rep_printer = rep_printer + rep_list
        #rep_printer = common_list + rep_printer

        printer = ""
        for i in common_list:
            printer = printer + str(i) + "\t"
        if len(common_list) < 15:
            printer = printer + "\t\t\t"
        for i in rep_printer:
            printer = printer + str(i) + "\t"

        printer = printer.strip("\t")

        print >> OUT, printer

    OUT.close()