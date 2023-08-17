import os
import argparse



parse = argparse.ArgumentParser(description="extract pair-mapped reads. Note: the pair-end sequence data need to be separeted by "_fwd" and "_rev".")
parse.add_argument("--input", help="the input list that contains the file to be processed, see example.txt", required=True)
parse.add_argument("--path", help="the path of file to be processed", required=True)
parse.add_argument("--out_o", help="output file of paired-mapped", required=True)
parse.add_argument("--out_s", help="output file of single-mapped", required=True)
args = parse.parse_args()


def proc_map_result(fn):
    info={}
    check={}
    with open(fn,'r') as fi:
         for line in fi:
             head=line.split()[1]
             read=line.split()[0]
             idt=float(line.split('\t')[2])
             cov=int(line.split('\t')[3])
             if idt<95 or cov<45:continue
             if read in check:
                check[read].append(head)
                info[read][head]=line
             else:
                check[read]=[head]
                info[read]={}
                info[read][head]=line
    return info,check



def get_info_overlap(listset,fwd_info,rev_info,head):
    data=[]
    for ele in listset:
        data.append(fwd_info[head][ele])
        data.append(rev_info[head][ele])
    return data

def get_info_u(listset,u_info,head):
    data=[]
    for ele in listset:
        data.append(u_info[head][ele])
    return data


def pair_match(fwd_info,fwd_check,rev_info,rev_check):
    count={}
    count_u={}
    for head in fwd_check:
        if head in rev_check:
           overlap=list(set(fwd_check[head])&set(rev_check[head]))
           fwd_u=list(set(fwd_check[head])-set(rev_check[head]))
           rev_u=list(set(rev_check[head])-set(fwd_check[head]))
           if len(overlap)>0:
              count[head]=get_info_overlap(overlap,fwd_info,rev_info,head)
              count_u[head]=get_info_u(fwd_u,fwd_info,head)+get_info_u(rev_u,rev_info,head)
    return count,count_u


obj=[line.split()[0] for line in open(args.input,'r')]


for ele in obj:
    map_cont_fwd,check_fwd=proc_map_result(args.path+'%s_fwd.txt'%ele)
    print('%s_fwd.txt'%ele)
    map_cont_rev,check_rev=proc_map_result(args.path+'%s_rev.txt'%ele)
    print('%s_rev.txt'%ele)
    map_pair,map_u=pair_match(map_cont_fwd,check_fwd,map_cont_rev,check_rev)
    fo=open(args.out_o+'%s.txt'%ele,'w')
    foo=open(args.out_s+'%s.txt'%ele,'w')
    for head in map_pair:
        for cont in map_pair[head]:
            fo.write(cont)
    fo.close()
    for head in map_u:
        for cont in map_u[head]:
            foo.write(cont)
    foo.close()
