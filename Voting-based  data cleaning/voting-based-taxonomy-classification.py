import numpy as np
import argparse


parse = argparse.ArgumentParser(description="voting-based data cleaning and taxonomy classification")
parse.add_argument("--in", help="input file", required=True)
parse.add_argument("--lineage", help="the reformated NCBI lineages", required=True)
parse.add_argument("--length", help="the file that contains the length of contigs or genomes. This file contains two columns seperated by Tab. the first column is the contig or genome name, the second column is the length (base pair) of the corresponding contigs or genomes.", required=True)
parse.add_argument("--out", help="output file", required=True)
args = parse.parse_args()



tax={}
with open(args.lineage,'r') as fi:
     for line in fi:
         try:
            tax[line.split('\t')[0]]=line.strip('\n').split('\t')[1]
         except IndexError:
            print(line)


leninfo={}
with open(args.length,'r') as fi:
     for line in fi:
         leninfo[line.split()[0]]=int(line.strip('\n').split('\t')[1])




info={}
info_1={}
with open(args.in,'r') as fi:
     for line in fi:
         evalue=float(line.split('\t')[10])
         if evalue>0.00001:continue
         if '___' in line.split('\t')[0]:
            tmp='___'.join(line.split('\t')[0].split('___')[:-1])
            _hh=line.split('\t')[0].split('___')[-1]
         else:
            tmp=line.split('\t')[0]
            _hh=line.split('\t')[0]
         if tmp.split('.')[-1]=='fasta':
            hh=tmp.split('.fasta')[0]
         elif tmp.split('.')[-1]=='fna':
            hh=tmp.split('.fna')[0]
         elif tmp.split('.')[-1]=='fa':
            hh=tmp.split('.fa')[0]
         else:
            hh=tmp
         qs=int(line.split('\t')[6])
         qe=int(line.split('\t')[7])
         ql=int(line.split('\t')[12])
         rag=set(range(qs,(qe+1)))
         kid=line.strip('\n').split('\t')[-3].split(';')
         if hh in info:
            for ele in kid:
                if ele not in tax:
                   continue
                else:
                   _tax=tax[ele]
                fam=','.join(_tax.split(',')[:5])
                if fam in info[hh]:
                   if _hh in info[hh][fam]:
                      info[hh][fam][_hh]=info[hh][fam][_hh]|rag
                      info_1[hh][fam.split(',')[-1]][_hh]=info_1[hh][fam.split(',')[-1]][_hh]|rag
                   else:
                      info[hh][fam][_hh]=rag
                      info_1[hh][fam.split(',')[-1]][_hh]=rag
                else:
                   info[hh][fam]={}
                   info_1[hh][fam.split(',')[-1]]={}
                   info[hh][fam][_hh]=rag
                   info_1[hh][fam.split(',')[-1]][_hh]=rag
         else:
            info[hh]={}
            info_1[hh]={}
            for ele in kid:
                if ele not in tax:
                   continue
                else:
                   _tax=tax[ele]
                fam=','.join(_tax.split(',')[:5])
                if fam in info[hh]:
                   if _hh in info[hh][fam]:
                      info[hh][fam][_hh]=info[hh][fam][_hh]|rag
                      info_1[hh][fam.split(',')[-1]][_hh]=info_1[hh][fam.split(',')[-1]][_hh]|rag
                   else:
                      info[hh][fam][_hh]=rag
                      info_1[hh][fam.split(',')[-1]][_hh]=rag
                else:
                   info[hh][fam]={}
                   info_1[hh][fam.split(',')[-1]]={}
                   info[hh][fam][_hh]=rag
                   info_1[hh][fam.split(',')[-1]][_hh]=rag



def rank_link(ranklist,high,low,data_high):
    for ele in ranklist:
        a=ele.split(',')[high]
        b=ele.split(',')[low]
        if a in data_high:
           data_high[a].append(b)
        else:
           data_high[a]=[b]
    return data_high




data_high={}
for hh in info:
    tmp=[]
    for fam in info[hh]:
        tmp.append(fam)
    tmp_high={}
    tmp_high=rank_link(tmp,0,1,tmp_high)
    tmp_high=rank_link(tmp,1,2,tmp_high)
    tmp_high=rank_link(tmp,2,3,tmp_high)
    data_high[hh]=rank_link(tmp,3,4,tmp_high)





def assign(ip,hh,data_high,info_1):
    for hrank in ip:
        if hrank not in info_1[hh]:
            info_1[hh][hrank]={}
        for lrank in data_high[hh][hrank]:
            for _hh in info_1[hh][lrank]:
                if _hh in info_1[hh][hrank]:
                   info_1[hh][hrank][_hh]=info_1[hh][hrank][_hh]|info_1[hh][lrank][_hh]
                else:
                   info_1[hh][hrank][_hh]=info_1[hh][lrank][_hh]
    return info_1


for hh in data_high:
    ip_o=[]
    ip_c=[]
    ip_p=[]
    ip_k=[]
    ip_f=[]
    for hrank in data_high[hh]:
        if 'f__' in hrank:
           ip_f.append(hrank)
        if 'o__' in hrank:
           ip_o.append(hrank)
        if 'c__' in hrank:
           ip_c.append(hrank)
        if 'p__' in hrank:
           ip_p.append(hrank)
        if 'k__' in hrank:
           ip_k.append(hrank)
    info_1=assign(ip_f,hh,data_high,info_1)
    info_1=assign(ip_o,hh,data_high,info_1)
    info_1=assign(ip_c,hh,data_high,info_1)
    info_1=assign(ip_p,hh,data_high,info_1)
    info_1=assign(ip_k,hh,data_high,info_1)




def cal_percent(dic,hh):
    tmp=0
    for i in dic:
        tmp+=len(list(dic[i]))
    if '.a___' in hh:
       _hh=hh.split('.a___')[1]
    else:
       _hh=hh
    return tmp/leninfo[_hh]*100


def sort(sortlist,info_2,hh):
    rank,score=sorted(sortlist,key=(lambda x:x[1]),reverse=True)[0]
    info_2[hh][rank]=score
    return info_2,rank


def dis(sortlist):
    kindom=[]
    tmp=sorted(sortlist,key=(lambda x:x[1]),reverse=True)
    for ele in tmp:
        kindom.append('%s:%f'%(ele[0],ele[1]))
    return kindom

def rfilter(ip,data_high,hh,record):
    _ip=[]
    for ele in ip:
        if ele[0] not in data_high[hh][record]:continue
        _ip.append(ele)
    return _ip



info_2={}
for hh in info_1:
    ip_o=[]
    ip_c=[]
    ip_p=[]
    ip_k=[]
    ip_f=[]
    info_2[hh]={}
    for rank in info_1[hh]:
        if 'f__' in rank:
           ip_f.append([rank,cal_percent(info_1[hh][rank],hh)])
        if 'o__' in rank:
           ip_o.append([rank,cal_percent(info_1[hh][rank],hh)])
        if 'c__' in rank:
           ip_c.append([rank,cal_percent(info_1[hh][rank],hh)])
        if 'p__' in rank:
           ip_p.append([rank,cal_percent(info_1[hh][rank],hh)])
        if 'k__' in rank:
           ip_k.append([rank,cal_percent(info_1[hh][rank],hh)])
    info_2,record=sort(ip_k,info_2,hh)
    kindom=dis(ip_k)
    ip_p=rfilter(ip_p,data_high,hh,record)
    info_2,record=sort(ip_p,info_2,hh)
    ip_c=rfilter(ip_c,data_high,hh,record)
    info_2,record=sort(ip_c,info_2,hh)
    ip_o=rfilter(ip_o,data_high,hh,record)
    info_2,record=sort(ip_o,info_2,hh)
    ip_f=rfilter(ip_f,data_high,hh,record)
    info_2,record=sort(ip_f,info_2,hh)






fo=open(args.out,'w')
for hh in info_2:
    tmp=[]
    for kin in info_2[hh]:
        tmp.append('%s:%f'%(kin,info_2[hh][kin]))
    fo.write('%s\t%s\t%s\n'%(hh,','.join(tmp),','.join(kindom)))
fo.close()

