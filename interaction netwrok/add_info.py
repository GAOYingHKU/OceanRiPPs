import argparse



parse = argparse.ArgumentParser(description="add information")
parse.add_argument("--input", help="input host-phage relationship, i.e. host_phage_pair_summary.00.txt", required=True)
parse.add_argument("--pair", help="the contigs and their belonging genomes, two columns, first is contigs, second is genomes", required=True)
parse.add_argument("--express", help="RiPPs exist in metatranscriptomic samples, two columns, first is RiPPs, second is contigs. Note: the format of RiPPs header need to be consistent throughout the whole process ", required=True)
parse.add_argument("--inform", help="i.e. supplementary table 1", required=True)
parse.add_argument("--vpc", help="information of viral protein families, two columns, first is the viral protein, second is viral protein families, seperated by ','. i.e. vpc.txt", required=True)
parse.add_argument("--node", help="WGCNA node information, i.e. example.txt", required=True)
parse.add_argument("--rep", help="RiPPs representation from CDHit results, two columns, first is RiPP family number, second is represntation, seperated by '\t'", required=True)
parse.add_argument("--out", help="output file", required=True)
args = parse.parse_args()




info={}
with open(args.pair,'r') as fi:
     for line in fi:
         info[line.split()[0]]=line.strip('\n').split('\t')[1]


data={}
with open(args.express,'r') as fi:
     for line in fi:
         ripp=line.split('\t')[0]
         cot=line.strip('\n').split('\t')[1]
         if cot in data:
            data[cot].append(ripp)
         else:
            data[cot]=[ripp]



df=pd.read_csv(args.inform)
RF,des,ptm=df['RiPP family'],df['description'],df['PTM']



info={}
ptm={}
for i,j,k in zip(RF,des,ptm):
    info[j]=i
    if k=='ptm':
       ptm[i]='high'
    else:
       ptm[i]='low'



query={}
with open(args.vpc,'r') as fi:
     for line in fi:
         a=line.split(',')[0]
         b='_'.join(a.split('_')[:-1])
         if 'PC_' in line:
            c=line.strip('\n').split(',')[1]
         else:
            c=a
         if b in query:
            query[b].append(c)
         else:
            query[b]=[c]



node={}
with open(args.node,'r') as fi:
     for line in fi:
         if 'nodeName' in line:continue
         color=line.strip('\n').split('\t')[-1]
         hh=line.split('\t')[0]
         node[hh]=color



rippclu={}
with open(args.rep,'r') as fi:
     for line in fi:
         clu=line.split('\t')[0]
         rep=line.strip('\n').split('\t')[1]
         rippclu[clu]=rep



score={}
with open(args.score,'r') as fi:
     for line in fi:
         a=line.split('\t')[0]
         b=line.split('\t')[1]
         c=line.split('\t')[2]
         d='%s|%s'%(b,a)
         score[d]=c

fo=open(args.out,'w')
with open(args.input,'r') as fi:
     for line in fi:
         if 'host' in line:continue
         a=line.split('\t')[0]
         if a not in data:continue
         for i in data[a]:
             for j in query[a]:
                 fo.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(info[i],ptm[info[i]],i,'\t'.join(line.split('\t')[:2]),j,'\t'.join(line.strip('\n').split('\t')[2:5]),score['%s|%s'%(info[i],j)]))


fo.close()

